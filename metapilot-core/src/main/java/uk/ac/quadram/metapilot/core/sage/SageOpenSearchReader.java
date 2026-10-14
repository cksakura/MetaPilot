package uk.ac.quadram.metapilot.core.sage;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.function.Consumer;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import uk.ac.quadram.metapilot.core.opensearch.OpenSearchPsm;

/**
 * Reads Sage's {@code results.sage.tsv} into {@link OpenSearchPsm}, so Sage
 * candidates flow through the same validation pipeline as MSFragger's.
 *
 * <h2>Columns are resolved by name, never by position</h2>
 * Sage's column set has changed between releases and the order is not
 * guaranteed. The existing {@code SagePsmReader} indexes some columns by
 * position ({@code cs[15]} for the target/decoy label), which silently
 * misreads a file whose layout differs. Everything here is looked up by header
 * name and absent columns become NaN rather than throwing, so a version that
 * drops or renames a column degrades one feature instead of failing the run.
 *
 * <h2>What Sage cannot supply</h2>
 * There is no counterpart to MSFragger's {@code localize_delta_mass} output, so
 * the delta-mass localization block stays unset. That is why Sage PSMs are
 * scored by their own classifier over
 * {@code OpenSearchFeatures.SHARED_FEATURE_NAMES}: mean-imputing those features
 * would give every Sage row an identical value and hand the model a perfect
 * marker of engine identity.
 *
 * <h2>Target/decoy</h2>
 * Decided from the protein accession against the shared decoy tag, not from
 * Sage's own {@code label} column, so both engines classify a hit the same way.
 * Sage's own {@code spectrum_q} is read but deliberately not used to filter:
 * FDR in this workflow is computed after the classifier and the spectrum-level
 * collapse, and pre-filtering by an engine's internal q-value would remove the
 * decoys that calculation needs.
 */
public final class SageOpenSearchReader {

	private static final Logger LOGGER = LogManager.getLogger(SageOpenSearchReader.class);

	private final String decoyTag;
	private final Set<String> unrecognisedColumns = new LinkedHashSet<>();
	private final Set<String> missingExpected = new LinkedHashSet<>();

	public SageOpenSearchReader(String decoyTag) {
		this.decoyTag = decoyTag == null || decoyTag.isBlank() ? "rev_" : decoyTag;
	}

	/** Columns this reader maps; anything else is reported once as unknown. */
	private static final Set<String> KNOWN = Set.of("psm_id", "filename", "scannr", "peptide", "proteins",
			"num_proteins", "rank", "label", "expmass", "calcmass", "charge", "peptide_len", "missed_cleavages",
			"semi_enzymatic", "isotope_error", "precursor_ppm", "fragment_ppm", "average_ppm", "hyperscore",
			"delta_next", "delta_best", "delta_mass", "rt", "aligned_rt", "predicted_rt", "delta_rt_model",
			"matched_peaks", "longest_b", "longest_y", "longest_y_pct", "matched_intensity_pct", "scored_candidates",
			"poisson", "sage_discriminant_score", "posterior_error", "spectrum_q", "peptide_q", "protein_q",
			"ms1_intensity", "ms2_intensity", "ion_mobility", "predicted_mobility", "delta_mobility");

	public List<OpenSearchPsm> read(File resultsTsv, String sampleBaseName, Consumer<String> progress)
			throws IOException {
		List<OpenSearchPsm> psms = new ArrayList<>();
		if (resultsTsv == null || !resultsTsv.exists() || resultsTsv.length() == 0) {
			throw new IOException("Sage results file missing or empty: " + resultsTsv);
		}
		try (BufferedReader reader = new BufferedReader(
				new InputStreamReader(new FileInputStream(resultsTsv), StandardCharsets.UTF_8))) {
			String header = reader.readLine();
			if (header == null) {
				throw new IOException("Sage results file has no header: " + resultsTsv);
			}
			Map<String, Integer> col = new HashMap<>();
			String[] names = header.split("\t", -1);
			for (int i = 0; i < names.length; i++) {
				String name = names[i].trim().toLowerCase(Locale.ROOT);
				col.put(name, i);
				if (!KNOWN.contains(name)) {
					unrecognisedColumns.add(names[i].trim());
				}
			}
			// Not delta_mass: Sage has no such column (checked against 0.14.x output).
			// The open-search delta mass is derived below from expmass - calcmass.
			for (String expected : new String[] { "peptide", "proteins", "hyperscore", "expmass", "calcmass",
					"charge" }) {
				if (!col.containsKey(expected)) {
					missingExpected.add(expected);
				}
			}
			if (!missingExpected.isEmpty()) {
				throw new IOException("Sage results file is missing required columns " + missingExpected
						+ " — cannot map it onto the shared PSM model: " + resultsTsv);
			}

			String line;
			int n = 0;
			while ((line = reader.readLine()) != null) {
				if (line.isBlank()) {
					continue;
				}
				String[] c = line.split("\t", -1);
				OpenSearchPsm psm = toPsm(c, col, sampleBaseName);
				if (psm != null) {
					psms.add(psm);
				}
				if (++n % 200000 == 0 && progress != null) {
					progress.accept(resultsTsv.getName() + ": " + n + " Sage rows read ...");
				}
			}
		}
		LOGGER.info("Read {} Sage PSM candidates from {}", psms.size(), resultsTsv.getName());
		return psms;
	}

	private OpenSearchPsm toPsm(String[] c, Map<String, Integer> col, String sampleBaseName) {
		OpenSearchPsm psm = new OpenSearchPsm();
		psm.setEngine(OpenSearchPsm.SearchEngine.SAGE);

		// Sage's filename is the mzML path; the workflow keys spectra on the
		// sample base name, so use the caller's rather than re-deriving it.
		psm.setFileName(sampleBaseName);
		psm.setScan(parseScan(str(c, col, "scannr")));

		String peptide = str(c, col, "peptide");
		psm.setPeptide(stripMods(peptide));
		// Sage encodes modifications inline, e.g. PEP[+15.9949]TIDE. Keeping the
		// decorated form as the localized sequence preserves the site information
		// Sage does have, without pretending it is MSFragger's localization score.
		psm.setLocalizedSequence(peptide);

		String proteins = str(c, col, "proteins");
		String[] all = proteins == null || proteins.isBlank() ? new String[0] : proteins.split(";");
		for (int i = 0; i < all.length; i++) {
			all[i] = all[i].trim();
		}
		psm.setAllProteins(all);
		String delegate = all.length > 0 ? all[0] : null;
		psm.setProtein(delegate);
		// From the accession, not Sage's label column, so both engines agree.
		psm.setTarget(delegate == null
				|| !delegate.toLowerCase(Locale.ROOT).startsWith(decoyTag.toLowerCase(Locale.ROOT)));

		psm.setHitRank((int) Math.max(1, num(c, col, "rank", 1)));
		psm.setCharge((int) num(c, col, "charge", 0));
		double expMass = num(c, col, "expmass", Double.NaN);
		double calcMass = num(c, col, "calcmass", Double.NaN);
		psm.setPrecursorNeutralMass(expMass);
		psm.setCalcNeutralPepMass(calcMass);
		// Sage writes no delta-mass column, so the open-search delta is derived the
		// same way MSFragger's massdiff is defined: observed minus theoretical.
		// A delta_mass column is still preferred if a future Sage adds one.
		double deltaMass = num(c, col, "delta_mass", Double.NaN);
		if (Double.isNaN(deltaMass)) {
			deltaMass = Double.isNaN(expMass) || Double.isNaN(calcMass) ? Double.NaN : expMass - calcMass;
		}
		psm.setMassDiff(deltaMass);
		// OpenSearchPsm.retentionTime is in SECONDS (MSFragger's retention_time_sec,
		// and MagOpenSearchHeadlessRunner.retentionTimeMinutes divides by 60 for
		// FlashLFQ). Sage's rt is in MINUTES - measured on 12mix YS-194-1-1: FragPipe
		// "Retention" / Sage rt = 60.0 (median over 63,544 shared scans). Stored
		// unconverted, Sage-only PSMs reached FlashLFQ at rt/60 minutes.
		double rtMinutes = num(c, col, "rt", Double.NaN);
		psm.setRetentionTime(Double.isNaN(rtMinutes) ? Double.NaN : rtMinutes * 60.0);
		psm.setNumMissedCleavages((int) num(c, col, "missed_cleavages", -1));

		// Sage reports semi_enzymatic as a boolean; num_tol_term is the count of
		// enzymatic termini, so map true -> 1, false -> 2.
		double semi = num(c, col, "semi_enzymatic", Double.NaN);
		if (!Double.isNaN(semi)) {
			psm.setNumTolTerm(semi > 0 ? 1 : 2);
		} else {
			psm.setNumTolTerm(-1);
		}

		psm.setNumMatchedIons((int) num(c, col, "matched_peaks", -1));
		// Sage has no total-theoretical-ion count. Left at -1 rather than guessed
		// from peptide length, which would make matched_ion_fraction mean something
		// different for each engine.
		psm.setTotNumIons(-1);

		psm.setHyperscore(num(c, col, "hyperscore", Double.NaN));
		// delta_next is (this hyperscore - next candidate's), so the next
		// candidate's score is hyperscore - delta_next. Reconstructing it rather
		// than storing delta_next directly keeps nextscore meaning the same thing
		// for both engines, and the derived hyperscore_minus_nextscore feature
		// then computes identically.
		double deltaNext = num(c, col, "delta_next", Double.NaN);
		double hyper = psm.getHyperscore();
		psm.setNextscore(Double.isNaN(deltaNext) || Double.isNaN(hyper) ? Double.NaN : hyper - deltaNext);

		// Sage's poisson is a log-scale score of match probability, playing the
		// same role as MSFragger's expect. Converted so neg_log10_expect is
		// comparable rather than on two different scales.
		double poisson = num(c, col, "poisson", Double.NaN);
		psm.setExpect(Double.isNaN(poisson) ? Double.NaN : Math.pow(10.0, poisson));

		Map<String, String> extras = new HashMap<>();
		putIfPresent(extras, c, col, "sage_discriminant_score");
		putIfPresent(extras, c, col, "spectrum_q");
		putIfPresent(extras, c, col, "peptide_q");
		putIfPresent(extras, c, col, "posterior_error");
		putIfPresent(extras, c, col, "longest_y_pct");
		putIfPresent(extras, c, col, "matched_intensity_pct");
		putIfPresent(extras, c, col, "scored_candidates");
		putIfPresent(extras, c, col, "isotope_error");
		psm.setSearchScores(extras);

		return psm;
	}

	/** Sage writes scan identifiers like {@code controllerType=0 ... scan=1234}. */
	static int parseScan(String scannr) {
		if (scannr == null) {
			return -1;
		}
		int idx = scannr.lastIndexOf("scan=");
		String digits = idx >= 0 ? scannr.substring(idx + 5) : scannr;
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < digits.length(); i++) {
			char ch = digits.charAt(i);
			if (Character.isDigit(ch)) {
				sb.append(ch);
			} else if (sb.length() > 0) {
				break;
			}
		}
		try {
			return sb.length() == 0 ? -1 : Integer.parseInt(sb.toString());
		} catch (NumberFormatException e) {
			return -1;
		}
	}

	/** Strips Sage's inline modification annotations, leaving the bare backbone. */
	static String stripMods(String peptide) {
		if (peptide == null) {
			return null;
		}
		StringBuilder sb = new StringBuilder(peptide.length());
		int depth = 0;
		for (int i = 0; i < peptide.length(); i++) {
			char ch = peptide.charAt(i);
			if (ch == '[' || ch == '(') {
				depth++;
			} else if (ch == ']' || ch == ')') {
				if (depth > 0) {
					depth--;
				}
			} else if (depth == 0 && Character.isLetter(ch)) {
				sb.append(Character.toUpperCase(ch));
			}
		}
		return sb.toString();
	}

	private static void putIfPresent(Map<String, String> into, String[] c, Map<String, Integer> col, String name) {
		String v = str(c, col, name);
		if (v != null && !v.isBlank()) {
			into.put(name, v);
		}
	}

	private static String str(String[] c, Map<String, Integer> col, String name) {
		Integer i = col.get(name);
		return i == null || i >= c.length ? null : c[i].trim();
	}

	private static double num(String[] c, Map<String, Integer> col, String name, double fallback) {
		String v = str(c, col, name);
		if (v == null || v.isBlank()) {
			return fallback;
		}
		if ("true".equalsIgnoreCase(v)) {
			return 1;
		}
		if ("false".equalsIgnoreCase(v)) {
			return 0;
		}
		try {
			return Double.parseDouble(v);
		} catch (NumberFormatException e) {
			return fallback;
		}
	}

	/** Reports columns this reader does not map, through the progress sink. */
	public void reportUnrecognisedColumns(Consumer<String> progress) {
		if (unrecognisedColumns.isEmpty()) {
			return;
		}
		String message = "Sage results contained columns this reader does not map: "
				+ String.join(", ", unrecognisedColumns)
				+ ". They are ignored; add them in SageOpenSearchReader if any should become features.";
		LOGGER.info(message);
		if (progress != null) {
			progress.accept(message);
		}
	}
}
