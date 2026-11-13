package uk.ac.quadram.metapilot.core.runner;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.IdentityHashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import uk.ac.quadram.metapilot.core.opensearch.OpenSearchPsm;

/**
 * The side quantification in {@code sequence_variants/}: substitution PSMs that are not
 * identifications, quantified by their own FlashLFQ run so that nothing in the main result moves.
 * <ul>
 * <li><b>Input.</b> {@code substitution_psms.tsv} holds the target PSMs with a substitution
 * verdict that are not in {@code validated_psms.tsv} ({@code quant_role=variant}), as the base
 * peptide plus its delta, never the guessed substitution. Beside them, the validated unmodified
 * target PSMs of the same base peptides ({@code quant_role=reference}), so each variant can be read
 * against its reference in one run, on one intensity scale.</li>
 * <li><b>Run.</b> No match-between-runs (these PSMs are not identifications, so they should not
 * transfer): FlashLFQ 2.0 cannot switch it off, so its {@code MBR} values are set to 0 and marked
 * {@code MBR_removed} in the peptide table and left out of the shared-feature check. Normalisation
 * off: FlashLFQ would otherwise fit per-file factors on the few hundred
 * variant and reference peptides of this file, which are unstable and on a different scale from
 * the main run's.</li>
 * <li><b>Output.</b> FlashLFQ's peptide table with {@code quant_role}, {@code variant_ids} (rows
 * of {@code sequence_variants.tsv}), {@code reference_sequence}, one {@code Ratio_<sample>} per
 * sample (variant over reference intensity, empty unless both are above zero), and
 * {@code shares_main_feature}: samples where the variant's MS1 peak is also a peak of a different
 * peptide in the main run (same file and charge, m/z within {@link #FEATURE_PPM} ppm, apex inside
 * the main peak). Two separate FlashLFQ runs cannot resolve a shared feature between them, so it is
 * flagged instead of hidden.</li>
 * </ul>
 */
final class SubstitutionQuantification {

	private static final Logger LOGGER = LogManager.getLogger(SubstitutionQuantification.class);

	/** Protein column. These PSMs are not identifications, so it must not be a searched protein. */
	static final String PROTEIN = "UNASSIGNED";

	static final String ROLE_VARIANT = "variant";
	static final String ROLE_REFERENCE = "reference";

	/** FlashLFQ's generic identification columns, then this table's own. */
	static final String HEADER = "File Name\tScan number\tBase Sequence\tFull Sequence\t"
			+ "Peptide Monoisotopic Mass\tScan Retention Time\tPrecursor Charge\tProtein Accession\tScore\tPEP\t"
			+ "Peptide Mass\tMass difference\tMissed cleavages\tquant_role\tvariant_id";

	private static final int FLASHLFQ_COLUMNS = 13;

	/** FlashLFQ's detection type for a match-between-runs transfer, and what it becomes here. */
	static final String MBR = "MBR";
	static final String MBR_REMOVED = "MBR_removed";

	/** FlashLFQ's own default peak-finding tolerance, so a shared feature is one it would also join. */
	static final double FEATURE_PPM = 10.0;

	/** Rows written to {@code substitution_psms.tsv}. */
	static final class Counts {
		final int variants;
		final int references;

		Counts(int variants, int references) {
			this.variants = variants;
			this.references = references;
		}
	}

	private SubstitutionQuantification() {
	}

	// --- input ------------------------------------------------------------

	/**
	 * Writes the side input from the PSMs in memory. {@code variantIds} maps a PSM to its
	 * {@code sequence_variants.tsv} id; null or missing means none.
	 */
	static Counts writePsms(List<OpenSearchPsm> candidates, List<OpenSearchPsm> validated,
			Map<OpenSearchPsm, String> variantIds, File out) throws IOException {
		Set<OpenSearchPsm> alreadyQuantified = Collections.newSetFromMap(new IdentityHashMap<>());
		if (validated != null) {
			alreadyQuantified.addAll(validated);
		}
		out.getParentFile().mkdirs();
		int variants = 0;
		int references = 0;
		Set<String> bases = new HashSet<>();
		try (PrintWriter writer = new PrintWriter(out, "UTF-8")) {
			writer.println(HEADER);
			if (candidates != null) {
				for (OpenSearchPsm psm : candidates) {
					if (psm == null || !psm.isTarget() || alreadyQuantified.contains(psm)) {
						continue;
					}
					String verdict = psm.getSubstitutionVerdict();
					if (verdict == null || verdict.isBlank()) {
						continue;
					}
					String id = variantIds == null ? null : variantIds.get(psm);
					writer.println(row(psm, ROLE_VARIANT, id == null ? "" : id));
					bases.add(psm.getPeptide());
					variants++;
				}
			}
			if (validated != null) {
				for (OpenSearchPsm psm : validated) {
					if (psm != null && psm.isTarget() && bases.contains(psm.getPeptide())
							&& psm.getPeptide().equals(OpenSearchPostAnalysisTsvExporter.canonicalModifiedSequenceOf(psm))) {
						writer.println(row(psm, ROLE_REFERENCE, ""));
						references++;
					}
				}
			}
		}
		return new Counts(variants, references);
	}

	private static String row(OpenSearchPsm psm, String role, String variantId) {
		return psm.getFileName() + "\t" + psm.getScan() + "\t" + psm.getPeptide() + "\t"
				+ OpenSearchPostAnalysisTsvExporter.canonicalModifiedSequenceOf(psm) + "\t"
				+ psm.getFlashLfqMonoisotopicMass() + "\t" + psm.getRetentionTime() / 60.0 + "\t" + psm.getCharge()
				+ "\t" + PROTEIN + "\t" + OpenSearchPostAnalysisTsvExporter.score(psm) + "\t" + psm.getQValue() + "\t"
				+ psm.getPrecursorNeutralMass() + "\t" + psm.getMassDiff() + "\t" + psm.getNumMissedCleavages() + "\t"
				+ role + "\t" + variantId;
	}

	/**
	 * Writes the side input from the tables on disk, for a resumed run whose validation did not run
	 * in this process. A variant row is a target row of {@code unexplained} that {@code forms} flags
	 * with a verdict and {@code validated} does not hold; a reference row is a target row of
	 * {@code validated} whose Full Sequence is its Base Sequence, for a base peptide with a variant row.
	 */
	static Counts writePsmsFromReports(File forms, File unexplained, File validated, File out) throws IOException {
		Map<String, String> flagged = new HashMap<>();
		List<String[]> formRows = rows(forms);
		if (!formRows.isEmpty()) {
			String[] t = formRows.get(0);
			int f = column(t, "file"), s = column(t, "scan"), r = column(t, "hit_rank"),
					v = column(t, "substitution_verdict"), id = column(t, "variant_id");
			for (String[] c : formRows.subList(1, formRows.size())) {
				if (!cell(c, v).isBlank()) {
					flagged.put(cell(c, f) + "	" + cell(c, s) + "	" + cell(c, r), cell(c, id));
				}
			}
		}
		// validated_psms.tsv can hold close to a million rows: streamed twice, never held.
		Set<String> quantified = new HashSet<>();
		stream(validated, (t, c) -> quantified.add(key(t, c)));

		out.getParentFile().mkdirs();
		int[] counts = new int[2];
		Set<String> bases = new HashSet<>();
		String[] names = HEADER.split("	");
		try (PrintWriter writer = new PrintWriter(out, "UTF-8")) {
			writer.println(HEADER);
			if (!flagged.isEmpty()) {
				stream(unexplained, (t, c) -> {
					String k = key(t, c);
					if ("true".equalsIgnoreCase(cell(c, column(t, "is_target"))) && flagged.containsKey(k)
							&& !quantified.contains(k)) {
						writer.println(copy(t, c, names, ROLE_VARIANT, flagged.get(k)));
						bases.add(cell(c, column(t, "Base Sequence")));
						counts[0]++;
					}
				});
			}
			if (!bases.isEmpty()) {
				stream(validated, (t, c) -> {
					String base = cell(c, column(t, "Base Sequence"));
					if ("true".equalsIgnoreCase(cell(c, column(t, "is_target"))) && bases.contains(base)
							&& base.equals(cell(c, column(t, "Full Sequence")))) {
						writer.println(copy(t, c, names, ROLE_REFERENCE, ""));
						counts[1]++;
					}
				});
			}
		}
		return new Counts(counts[0], counts[1]);
	}

	private interface RowSink {
		void accept(String[] title, String[] cells);
	}

	/** Calls {@code sink} for every data row of a tab-separated file; nothing when it is missing. */
	private static void stream(File file, RowSink sink) throws IOException {
		if (file == null || !file.isFile()) {
			return;
		}
		try (BufferedReader r = reader(file)) {
			String header = r.readLine();
			if (header == null) {
				return;
			}
			String[] title = header.split("	", -1);
			String line;
			while ((line = r.readLine()) != null) {
				if (!line.isBlank()) {
					sink.accept(title, line.split("	", -1));
				}
			}
		}
	}

	private static String key(String[] title, String[] c) {
		return cell(c, column(title, "File Name")) + "\t" + cell(c, column(title, "Scan number")) + "\t"
				+ cell(c, column(title, "hit_rank"));
	}

	private static String copy(String[] title, String[] c, String[] names, String role, String variantId) {
		StringBuilder row = new StringBuilder();
		for (int i = 0; i < FLASHLFQ_COLUMNS; i++) {
			if (i > 0) {
				row.append('\t');
			}
			row.append("Protein Accession".equals(names[i]) ? PROTEIN : cell(c, column(title, names[i])));
		}
		return row.append('\t').append(role).append('\t').append(variantId == null ? "" : variantId).toString();
	}

	/** True when {@code psmFile} is in this format (an older run wrote it without roles). */
	static boolean hasRoles(File psmFile) throws IOException {
		if (psmFile == null || !psmFile.isFile()) {
			return false;
		}
		try (BufferedReader r = reader(psmFile)) {
			String header = r.readLine();
			return header != null && Arrays.asList(header.split("\t", -1)).contains("quant_role");
		}
	}

	/** Variant rows in {@code psmFile}. */
	static int variantRows(File psmFile) throws IOException {
		List<String[]> all = rows(psmFile);
		if (all.isEmpty()) {
			return 0;
		}
		int role = column(all.get(0), "quant_role");
		int n = 0;
		for (String[] c : all.subList(1, all.size())) {
			if (ROLE_VARIANT.equals(cell(c, role))) {
				n++;
			}
		}
		return n;
	}

	// --- output -----------------------------------------------------------

	/** One MS1 peak of FlashLFQ's peak table. */
	private static final class Peak {
		final String sequence;
		final double mz;
		final double start;
		final double end;

		Peak(String sequence, double mz, double start, double end) {
			this.sequence = sequence;
			this.mz = mz;
			this.start = start;
			this.end = end;
		}
	}

	/**
	 * Rewrites {@code peptides} (FlashLFQ's peptide table from {@code psmFile}) with the role,
	 * variant ids, reference sequence, per-sample ratios and shared main-run features appended.
	 * {@code sidePeaks} and {@code mainPeaks} may be missing; {@code shares_main_feature} is then
	 * empty. Returns the number of variant peptides that share a main-run feature in any sample.
	 */
	static int annotate(File psmFile, File peptides, File sidePeaks, File mainPeaks) throws IOException {
		Map<String, String> roleOf = new HashMap<>();
		Map<String, TreeSet<String>> idsOf = new HashMap<>();
		Map<String, String> baseOf = new HashMap<>();
		List<String[]> psms = rows(psmFile);
		if (!psms.isEmpty()) {
			String[] t = psms.get(0);
			int full = column(t, "Full Sequence"), base = column(t, "Base Sequence"), role = column(t, "quant_role"),
					id = column(t, "variant_id");
			for (String[] c : psms.subList(1, psms.size())) {
				String seq = cell(c, full);
				// A sequence quantified as a reference is a reference, whatever else names it.
				if (!ROLE_REFERENCE.equals(roleOf.get(seq))) {
					roleOf.put(seq, cell(c, role));
				}
				baseOf.put(seq, cell(c, base));
				if (!cell(c, id).isBlank()) {
					idsOf.computeIfAbsent(seq, k -> new TreeSet<>()).add(cell(c, id));
				}
			}
		}
		Map<String, TreeSet<String>> shared = sharedFeatures(sidePeaks, mainPeaks, roleOf);

		List<String[]> table = rows(peptides);
		if (table.isEmpty()) {
			return 0;
		}
		String[] title = table.get(0);
		int seqCol = column(title, "Sequence");
		List<Integer> intensityCols = new ArrayList<>();
		for (int i = 0; i < title.length; i++) {
			if (title[i].startsWith("Intensity_")) {
				intensityCols.add(i);
			}
		}
		// FlashLFQ cannot be told not to match between runs (see FlashLfqTask.addSwitches); its MBR
		// values are dropped here. MBR fills only files without an MS/MS detection of the peptide, so
		// the MS/MS values are the same as in a run without it.
		int removed = 0;
		for (String[] c : table.subList(1, table.size())) {
			for (int i : intensityCols) {
				int type = column(title, "Detection Type_" + title[i].substring("Intensity_".length()));
				if (type >= 0 && type < c.length && i < c.length && MBR.equals(c[type])) {
					c[i] = "0";
					c[type] = MBR_REMOVED;
					removed++;
				}
			}
		}
		if (removed > 0) {
			LOGGER.info("Substitution quantification: {} match-between-runs value(s) dropped", removed);
		}
		Map<String, String[]> referenceRow = new HashMap<>();
		for (String[] c : table.subList(1, table.size())) {
			if (ROLE_REFERENCE.equals(roleOf.get(cell(c, seqCol)))) {
				referenceRow.put(cell(c, seqCol), c);
			}
		}

		int sharing = 0;
		File tmp = new File(peptides.getParentFile(), peptides.getName() + ".tmp");
		try (PrintWriter w = new PrintWriter(tmp, "UTF-8")) {
			StringBuilder head = new StringBuilder(String.join("\t", title));
			head.append("\tquant_role\tvariant_ids\treference_sequence");
			for (int i : intensityCols) {
				head.append("\tRatio_").append(title[i].substring("Intensity_".length()));
			}
			head.append("\tshares_main_feature");
			w.println(head);
			for (String[] c : table.subList(1, table.size())) {
				String seq = cell(c, seqCol);
				String role = roleOf.getOrDefault(seq, "");
				boolean variant = ROLE_VARIANT.equals(role);
				String base = baseOf.getOrDefault(seq, "");
				String[] ref = variant ? referenceRow.get(base) : null;
				StringBuilder line = new StringBuilder();
				for (int i = 0; i < title.length; i++) {
					line.append(i > 0 ? "\t" : "").append(cell(c, i));
				}
				line.append('\t').append(role);
				line.append('\t').append(String.join(";", idsOf.getOrDefault(seq, new TreeSet<>())));
				line.append('\t').append(variant && ref != null ? base : "");
				for (int i : intensityCols) {
					line.append('\t');
					if (ref != null) {
						double v = number(cell(c, i));
						double r = number(cell(ref, i));
						if (v > 0 && r > 0) {
							line.append(String.format(Locale.ROOT, "%.6g", v / r));
						}
					}
				}
				TreeSet<String> s = shared.get(seq);
				line.append('\t').append(s == null ? "" : String.join(";", s));
				if (s != null && variant) {
					sharing++;
				}
				w.println(line);
			}
		}
		Files.move(tmp.toPath(), peptides.toPath(), java.nio.file.StandardCopyOption.REPLACE_EXISTING);
		return sharing;
	}

	/**
	 * For each variant sequence, the {@code sample:main sequence} pairs where its peak is also a
	 * main-run peak of a different sequence.
	 */
	private static Map<String, TreeSet<String>> sharedFeatures(File sidePeaks, File mainPeaks,
			Map<String, String> roleOf) throws IOException {
		Map<String, TreeSet<String>> out = new HashMap<>();
		if (sidePeaks == null || !sidePeaks.isFile() || mainPeaks == null || !mainPeaks.isFile()) {
			return out;
		}
		// Only the side peaks are held as rows; the main table (often a million rows) is streamed once
		// and indexed by file and charge.
		Map<String, List<double[]>> wanted = new HashMap<>();
		Map<String, List<String>> wantedSeq = new HashMap<>();
		List<String[]> side = rows(sidePeaks);
		if (side.isEmpty()) {
			return out;
		}
		String[] st = side.get(0);
		int sf = column(st, "File Name"), ss = column(st, "Full Sequence"), smz = column(st, "Peak MZ"),
				sz = column(st, "Peak Charge"), sapex = column(st, "Peak RT Apex"), si = column(st, "Peak intensity"),
				stype = column(st, "Peak Detection Type");
		for (String[] c : side.subList(1, side.size())) {
			String seq = cell(c, ss);
			double mz = number(cell(c, smz));
			if (!ROLE_VARIANT.equals(roleOf.get(seq)) || MBR.equals(cell(c, stype)) || !(number(cell(c, si)) > 0)
					|| !(mz > 0)) {
				continue;
			}
			String k = cell(c, sf) + "\t" + cell(c, sz);
			wanted.computeIfAbsent(k, x -> new ArrayList<>()).add(new double[] { mz, number(cell(c, sapex)) });
			wantedSeq.computeIfAbsent(k, x -> new ArrayList<>()).add(seq);
		}
		if (wanted.isEmpty()) {
			return out;
		}
		Map<String, List<Peak>> main = new LinkedHashMap<>();
		try (BufferedReader r = reader(mainPeaks)) {
			String header = r.readLine();
			if (header == null) {
				return out;
			}
			String[] mt = header.split("\t", -1);
			int mf = column(mt, "File Name"), ms = column(mt, "Full Sequence"), mmz = column(mt, "Peak MZ"),
					mzc = column(mt, "Peak Charge"), mstart = column(mt, "Peak RT Start"),
					mend = column(mt, "Peak RT End");
			String line;
			while ((line = r.readLine()) != null) {
				String[] c = line.split("\t", -1);
				String k = cell(c, mf) + "\t" + cell(c, mzc);
				if (!wanted.containsKey(k)) {
					continue;
				}
				double mz = number(cell(c, mmz));
				if (mz > 0) {
					main.computeIfAbsent(k, x -> new ArrayList<>())
							.add(new Peak(cell(c, ms), mz, number(cell(c, mstart)), number(cell(c, mend))));
				}
			}
		}
		for (Map.Entry<String, List<double[]>> e : wanted.entrySet()) {
			List<Peak> peaks = main.get(e.getKey());
			if (peaks == null) {
				continue;
			}
			peaks.sort((a, b) -> Double.compare(a.mz, b.mz));
			double[] mzs = peaks.stream().mapToDouble(p -> p.mz).toArray();
			String sample = e.getKey().substring(0, e.getKey().indexOf('\t'));
			List<String> seqs = wantedSeq.get(e.getKey());
			for (int i = 0; i < e.getValue().size(); i++) {
				double mz = e.getValue().get(i)[0];
				double apex = e.getValue().get(i)[1];
				double tol = mz * FEATURE_PPM * 1e-6;
				int from = Arrays.binarySearch(mzs, mz - tol);
				from = from < 0 ? -from - 1 : from;
				for (int j = from; j < mzs.length && mzs[j] <= mz + tol; j++) {
					Peak p = peaks.get(j);
					if (!p.sequence.equals(seqs.get(i)) && apex >= p.start && apex <= p.end) {
						out.computeIfAbsent(seqs.get(i), x -> new TreeSet<>()).add(sample + ":" + p.sequence);
					}
				}
			}
		}
		return out;
	}

	// --- spectrum folder --------------------------------------------------

	/**
	 * {@code ExperimentalDesign.tsv} in the spectrum folder before a side FlashLFQ run.
	 * {@code FlashLfqTask} rewrites it for the run; {@link #restore} puts it back, or removes it
	 * when there was none. Nothing else in the folder is touched: another job may be writing there,
	 * so a file that appeared during the run is reported, never deleted.
	 */
	static final class SpectraFolderGuard {
		static final String DESIGN = "ExperimentalDesign.tsv";

		private final File dir;
		private final Set<String> names;
		private final byte[] design;
		private final long designModified;

		private SpectraFolderGuard(File dir, Set<String> names, byte[] design, long designModified) {
			this.dir = dir;
			this.names = names;
			this.design = design;
			this.designModified = designModified;
		}

		static SpectraFolderGuard take(File dir) throws IOException {
			Set<String> names = new HashSet<>();
			File[] children = dir.listFiles();
			if (children != null) {
				for (File child : children) {
					names.add(child.getName());
				}
			}
			File design = new File(dir, DESIGN);
			byte[] bytes = design.isFile() ? Files.readAllBytes(design.toPath()) : null;
			return new SpectraFolderGuard(dir, names, bytes, design.isFile() ? design.lastModified() : -1L);
		}

		/** Restores the design file; returns the other files that appeared, left in place. */
		List<String> restore() {
			List<String> appeared = new ArrayList<>();
			try {
				File designFile = new File(dir, DESIGN);
				if (design != null) {
					Files.write(designFile.toPath(), design);
					if (designModified >= 0) {
						designFile.setLastModified(designModified);
					}
				} else {
					Files.deleteIfExists(designFile.toPath());
				}
				File[] children = dir.listFiles();
				if (children != null) {
					for (File child : children) {
						if (!names.contains(child.getName()) && !DESIGN.equals(child.getName())) {
							appeared.add(child.getName());
						}
					}
				}
				if (!appeared.isEmpty()) {
					LOGGER.warn("Files appeared in {} during substitution quantification and were left in place: {}",
							dir.getAbsolutePath(), appeared);
				}
			} catch (IOException e) {
				LOGGER.warn("Could not restore the spectrum folder after substitution quantification: "
						+ dir.getAbsolutePath(), e);
			}
			Collections.sort(appeared);
			return appeared;
		}
	}

	// --- table helpers ----------------------------------------------------

	/** Every line split on tabs, header first; empty when the file is missing or empty. */
	private static List<String[]> rows(File file) throws IOException {
		List<String[]> out = new ArrayList<>();
		if (file == null || !file.isFile()) {
			return out;
		}
		try (BufferedReader r = reader(file)) {
			String line;
			while ((line = r.readLine()) != null) {
				if (line.startsWith("#") || (line.isBlank() && !out.isEmpty())) {
					continue;
				}
				out.add(line.split("\t", -1));
			}
		}
		return out;
	}

	private static BufferedReader reader(File file) throws IOException {
		return new BufferedReader(new InputStreamReader(new FileInputStream(file), StandardCharsets.UTF_8));
	}

	private static int column(String[] title, String name) {
		for (int i = 0; i < title.length; i++) {
			if (name.equals(title[i])) {
				return i;
			}
		}
		return -1;
	}

	private static String cell(String[] cells, int index) {
		return index >= 0 && index < cells.length ? cells[index] : "";
	}

	private static double number(String s) {
		if (s == null || s.isBlank()) {
			return Double.NaN;
		}
		try {
			return Double.parseDouble(s.trim());
		} catch (NumberFormatException e) {
			return Double.NaN;
		}
	}
}
