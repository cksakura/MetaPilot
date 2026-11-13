package uk.ac.quadram.metapilot.core.opensearch;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;

/**
 * Rebuilds {@link OpenSearchPsm}s from {@code validated_psms.tsv} so post-analysis
 * can run without re-searching or re-validating. The Full Sequence and peptide
 * monoisotopic mass already written by the original run are kept: they are the
 * FlashLFQ join keys. Shift columns, when missing from that table, are filled
 * from {@code scored_psms.tsv} spectrum winners.
 */
public final class OpenSearchValidatedPsmReader {

	private OpenSearchValidatedPsmReader() {
	}

	public static final class Result {
		public final List<OpenSearchPsm> psms;
		public final Set<OpenSearchPsm> identificationOnly;

		Result(List<OpenSearchPsm> psms, Set<OpenSearchPsm> identificationOnly) {
			this.psms = psms;
			this.identificationOnly = identificationOnly;
		}
	}

	public static Result read(File validatedPsms) throws IOException {
		if (validatedPsms == null || !validatedPsms.isFile() || validatedPsms.length() == 0) {
			throw new IOException("Missing or empty validated_psms.tsv: " + validatedPsms);
		}
		List<OpenSearchPsm> psms = new ArrayList<>();
		Set<OpenSearchPsm> identificationOnly = Collections.newSetFromMap(new IdentityHashMap<>());
		try (BufferedReader reader = Files.newBufferedReader(validatedPsms.toPath(), StandardCharsets.UTF_8)) {
			String headerLine = reader.readLine();
			if (headerLine == null) {
				throw new IOException("Empty validated_psms.tsv: " + validatedPsms);
			}
			Map<String, Integer> col = index(headerLine);
			require(col, "File Name", "Scan number", "Base Sequence");
			String line;
			while ((line = reader.readLine()) != null) {
				if (line.isEmpty()) {
					continue;
				}
				String[] row = line.split("\t", -1);
				OpenSearchPsm psm = fromValidatedRow(col, row);
				psms.add(psm);
				if ("IDENTIFICATION_ONLY".equalsIgnoreCase(cell(col, row, "validation_tier"))) {
					identificationOnly.add(psm);
				}
			}
		}
		return new Result(psms, identificationOnly);
	}

	/**
	 * Copies fragment-supported localisation from {@code scored_psms.tsv} onto
	 * validated PSMs that do not already have a shift state. Join is
	 * file + scan + peptide on spectrum winners.
	 */
	public static int overlayShifts(List<OpenSearchPsm> psms, File scoredPsms) throws IOException {
		if (psms == null || psms.isEmpty() || scoredPsms == null || !scoredPsms.isFile()) {
			return 0;
		}
		boolean needsOverlay = false;
		for (OpenSearchPsm psm : psms) {
			if (psm.getShiftLocalisation() == null) {
				needsOverlay = true;
				break;
			}
		}
		if (!needsOverlay) {
			return 0;
		}
		Map<String, ShiftRecord> shifts = readScoredShifts(scoredPsms);
		int n = 0;
		for (OpenSearchPsm psm : psms) {
			if (psm.getShiftLocalisation() != null) {
				continue;
			}
			ShiftRecord rec = shifts.get(key(psm.getFileName(), psm.getScan(), psm.getPeptide()));
			if (rec == null) {
				continue;
			}
			apply(psm, rec);
			n++;
		}
		return n;
	}

	private static Map<String, ShiftRecord> readScoredShifts(File scoredPsms) throws IOException {
		Map<String, ShiftRecord> shifts = new HashMap<>();
		try (BufferedReader reader = Files.newBufferedReader(scoredPsms.toPath(), StandardCharsets.UTF_8)) {
			String headerLine = reader.readLine();
			if (headerLine == null) {
				return shifts;
			}
			Map<String, Integer> col = index(headerLine);
			if (!col.containsKey("shift_state") || !col.containsKey("is_spectrum_winner")) {
				return shifts;
			}
			String line;
			while ((line = reader.readLine()) != null) {
				if (line.isEmpty()) {
					continue;
				}
				String[] row = line.split("\t", -1);
				if (!Boolean.parseBoolean(cell(col, row, "is_spectrum_winner"))) {
					continue;
				}
				String state = cell(col, row, "shift_state");
				if (state.isEmpty()) {
					continue;
				}
				shifts.put(key(cell(col, row, "file"), parseInt(cell(col, row, "scan"), -1), cell(col, row, "peptide")),
						new ShiftRecord(state, parseInt(cell(col, row, "shift_site"), 0),
								parseDouble(cell(col, row, "shift_site_probability")),
								parseDouble(cell(col, row, "shift_support")),
								cell(col, row, "shift_explanation")));
			}
		}
		return shifts;
	}

	private static OpenSearchPsm fromValidatedRow(Map<String, Integer> col, String[] row) {
		OpenSearchPsm psm = new OpenSearchPsm();
		psm.setFileName(cell(col, row, "File Name"));
		psm.setScan(parseInt(cell(col, row, "Scan number"), 0));
		psm.setPeptide(cell(col, row, "Base Sequence"));
		String full = cell(col, row, "Full Sequence");
		if (!full.isEmpty()) {
			psm.setResolvedModifiedSequence(full);
		}
		double mono = parseDouble(cell(col, row, "Peptide Monoisotopic Mass"));
		if (!Double.isNaN(mono)) {
			psm.setResolvedMonoisotopicMass(mono);
		}
		double rtMin = parseDouble(cell(col, row, "Scan Retention Time"));
		if (!Double.isNaN(rtMin)) {
			psm.setRetentionTime(rtMin * 60.0);
		}
		psm.setCharge(parseInt(cell(col, row, "Precursor Charge"), 0));
		String proteins = cell(col, row, "Protein Accession");
		if (!proteins.isEmpty()) {
			psm.setAllProteins(proteins.split(";"));
		}
		psm.setClassScore(parseDouble(cell(col, row, "Score")));
		psm.setQValue(parseDouble(cell(col, row, "PEP")));
		double observed = parseDouble(cell(col, row, "Peptide Mass"));
		if (!Double.isNaN(observed)) {
			psm.setPrecursorNeutralMass(observed);
		}
		psm.setMassDiff(parseDouble(cell(col, row, "Mass difference")));
		psm.setNumMissedCleavages(parseInt(cell(col, row, "Missed cleavages"), 0));
		psm.setHitRank(parseInt(cell(col, row, "hit_rank"), 1));
		psm.setEngine(parseEngine(cell(col, row, "engine")));
		psm.setEngineSupport(cell(col, row, "engine_support"));
		psm.setLocalizedSequence(emptyToNull(cell(col, row, "localized_sequence")));
		psm.setProtein(emptyToNull(cell(col, row, "delegate_protein")));
		psm.setTarget(Boolean.parseBoolean(cell(col, row, "is_target")));
		psm.setNormalizedMassDiff(parseDouble(cell(col, row, "normalized_mass_diff")));
		psm.setObservedPeakId(parseInt(cell(col, row, "observed_peak_id"), -1));
		psm.setPeakId(parseInt(cell(col, row, "peak_id"), -1));
		psm.setPeakIsotopeOffset(parseInt(cell(col, row, "peak_isotope_offset"), 0));
		psm.setHyperscore(parseDouble(cell(col, row, "hyperscore")));
		psm.setExpect(parseDouble(cell(col, row, "expect")));
		psm.setClassLogit(parseDouble(cell(col, row, "class_logit")));
		if (!Double.isNaN(parseDouble(cell(col, row, "q_value")))) {
			psm.setQValue(parseDouble(cell(col, row, "q_value")));
		}
		psm.setIdentificationGroup(emptyToNull(cell(col, row, "id_group")));
		psm.setClosedStage("closed".equalsIgnoreCase(cell(col, row, "stage")));
		String state = cell(col, row, "shift_state");
		if (!state.isEmpty()) {
			apply(psm, new ShiftRecord(state, parseInt(cell(col, row, "shift_site"), 0),
					parseDouble(cell(col, row, "shift_site_probability")), Double.NaN,
					cell(col, row, "shift_explanation")));
		}
		return psm;
	}

	private static void apply(OpenSearchPsm psm, ShiftRecord rec) {
		ShiftLocaliser.State state = parseState(rec.state);
		if (state == null) {
			return;
		}
		int site = rec.site1Based > 0 ? rec.site1Based - 1 : -1;
		List<String> explanations = rec.explanation.isEmpty() ? List.of() : List.of(rec.explanation.split("\\|"));
		psm.setShiftLocalisation(new ShiftLocaliser.Result(rec.support, site, rec.probability, state, explanations,
				Double.NaN));
		if (!rec.explanation.isEmpty()) {
			psm.setResolvedShiftExplanation(rec.explanation);
			psm.setResolvedShiftStateAndClass(state.name().toLowerCase(Locale.ROOT),
					psm.getShiftLocalisation().explanationClassText());
		}
	}

	private static ShiftLocaliser.State parseState(String raw) {
		if (raw == null || raw.isBlank()) {
			return null;
		}
		try {
			return ShiftLocaliser.State.valueOf(raw.trim().toUpperCase(Locale.ROOT));
		} catch (IllegalArgumentException e) {
			return null;
		}
	}

	private static OpenSearchPsm.SearchEngine parseEngine(String raw) {
		if (raw == null || raw.isBlank()) {
			return OpenSearchPsm.SearchEngine.MSFRAGGER;
		}
		try {
			return OpenSearchPsm.SearchEngine.valueOf(raw.trim().toUpperCase(Locale.ROOT));
		} catch (IllegalArgumentException e) {
			return OpenSearchPsm.SearchEngine.MSFRAGGER;
		}
	}

	private static Map<String, Integer> index(String header) {
		String[] names = header.split("\t", -1);
		Map<String, Integer> col = new HashMap<>();
		for (int i = 0; i < names.length; i++) {
			col.put(names[i], i);
		}
		return col;
	}

	private static void require(Map<String, Integer> col, String... names) throws IOException {
		for (String name : names) {
			if (!col.containsKey(name)) {
				throw new IOException("validated_psms.tsv is missing column: " + name);
			}
		}
	}

	private static String cell(Map<String, Integer> col, String[] row, String name) {
		Integer i = col.get(name);
		if (i == null || i < 0 || i >= row.length) {
			return "";
		}
		return row[i] == null ? "" : row[i];
	}

	private static int parseInt(String raw, int fallback) {
		if (raw == null || raw.isBlank()) {
			return fallback;
		}
		try {
			return Integer.parseInt(raw.trim());
		} catch (NumberFormatException e) {
			return fallback;
		}
	}

	private static double parseDouble(String raw) {
		if (raw == null || raw.isBlank()) {
			return Double.NaN;
		}
		try {
			return Double.parseDouble(raw.trim());
		} catch (NumberFormatException e) {
			return Double.NaN;
		}
	}

	private static String emptyToNull(String raw) {
		return raw == null || raw.isEmpty() ? null : raw;
	}

	private static String key(String file, int scan, String peptide) {
		return file + '\0' + scan + '\0' + peptide;
	}

	private static final class ShiftRecord {
		final String state;
		final int site1Based;
		final double probability;
		final double support;
		final String explanation;

		ShiftRecord(String state, int site1Based, double probability, double support, String explanation) {
			this.state = state;
			this.site1Based = site1Based;
			this.probability = probability;
			this.support = support;
			this.explanation = explanation == null ? "" : explanation;
		}
	}
}
