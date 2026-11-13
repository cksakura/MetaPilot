package uk.ac.quadram.metapilot.core.dia;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import uk.ac.quadram.metapilot.core.diann.DiaNNHeadlessRunner;
import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;
import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNParquetReader;
import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNPrecursor;
import uk.ac.quadram.metapilot.dbSearch.diann.DiannParameter;
import uk.ac.quadram.metapilot.dbSearch.diann.DiannPrMatrixReader;
import uk.ac.quadram.metapilot.quant.flashLFQ.FlashLfqQuanPeptide;

/**
 * Host handling shared by all DIA workflows (Marker and Seed, every mode, CLI and GUI).
 * <p>
 * Strategy: the host library takes no part in the early searches (marker/seed search, genome selection, Deep
 * rounds). It is searched together with the microbial library only in the final combined search, as further
 * {@code .parquet} libraries of the same DIA-NN run (several Standalone libraries are allowed), so host and microbial peptides compete for the spectra under one
 * FDR. The final search records the host library in {@code combined/host_library.txt}.
 * <p>
 * Export rule (as in DDA, where the host takes shared peptides): every identified precursor whose stripped sequence
 * is in the host library is a host precursor, also when the same peptide is in the microbial library. Host
 * precursors are removed before genome/protein scoring and written to {@code host_peptides.tsv}; per-run counts and
 * intensities of host and microbial peptides go to {@code host_summary.tsv}.
 */
public final class DiaHostLibrarySupport {

	private static final Logger LOGGER = LogManager.getLogger(DiaHostLibrarySupport.class);

	/** Written next to combined.parquet by the final search; one absolute host .parquet path per line. */
	public static final String HOST_RECORD = "host_library.txt";
	public static final String HOST_PEPTIDES_TSV = "host_peptides.tsv";
	public static final String HOST_SUMMARY_TSV = "host_summary.tsv";

	private DiaHostLibrarySupport() {
	}

	// ------------------------------------------------------------------ search side

	/**
	 * Returns {@code library} as a {@code .parquet} library: a {@code .parquet} is used as it is; a {@code .speclib}
	 * is converted to {@code <workDir>/<baseName>.lib.parquet}; a protein FASTA is first predicted into
	 * {@code <workDir>/<baseName>.speclib} with {@code libraryBuild} (normal protein rules), then converted.
	 * Existing outputs are reused.
	 */
	public static File ensureParquetLibrary(File library, File workDir, String baseName, DiannParameter libraryBuild,
			File diannExe, ProgressCallback cb, File condaPrefix) throws Exception {
		if (library == null || !library.exists()) {
			throw new IOException("Library not found: " + (library == null ? "null" : library.getAbsolutePath()));
		}
		String name = library.getName().toLowerCase(Locale.ROOT);
		if (name.endsWith(".parquet")) {
			return library;
		}
		if (!workDir.isDirectory() && !workDir.mkdirs()) {
			throw new IOException("Could not create " + workDir.getAbsolutePath());
		}
		File speclib = library;
		if (name.endsWith(".fasta") || name.endsWith(".fa")) {
			speclib = new File(workDir, baseName + ".speclib");
			if (!speclib.exists()) {
				cb.onMessage("Host library: predicting " + speclib.getName() + " from " + library.getAbsolutePath());
				int exit = DiaNNHeadlessRunner.run(diannExe,
						libraryBuild.generateFasta(library.getAbsolutePath(), speclib.getAbsolutePath(), false), workDir,
						cb, condaPrefix);
				if (exit != 0) LOGGER.warn("DIA-NN library prediction exit code {}", exit);
				File predicted = new File(workDir, baseName + ".predicted.speclib");
				if (!speclib.exists() && predicted.exists()) Files.move(predicted.toPath(), speclib.toPath());
			}
			if (!speclib.exists()) throw new IOException("DIA-NN did not produce " + speclib.getAbsolutePath());
		} else if (!name.endsWith(".speclib")) {
			throw new IOException("Unsupported library type (expected .parquet, .speclib, .fasta or .fa): "
					+ library.getAbsolutePath());
		}
		File parquet = new File(workDir, baseName + ".lib.parquet");
		if (!parquet.exists() || parquet.length() == 0) {
			cb.onMessage("Converting " + speclib.getName() + " to " + parquet.getName()
					+ " for the combined search with the host library");
			int exit = DiaNNHeadlessRunner.run(diannExe,
					libraryBuild.generateSpeclib2Parquet(speclib.getAbsolutePath(), parquet.getAbsolutePath()), workDir,
					cb, condaPrefix);
			if (exit != 0) LOGGER.warn("DIA-NN speclib to parquet exit code {}", exit);
		}
		if (!parquet.exists() || parquet.length() == 0) {
			throw new IOException("Failed to convert " + speclib.getAbsolutePath() + " to " + parquet.getAbsolutePath());
		}
		return parquet;
	}

	/**
	 * Libraries of the final combined search. Without a host: {@code { microbialLibrary }} unchanged. With host
	 * libraries: {@code { microbial .parquet, host .parquet, ... }}, converting the microbial library next to itself
	 * ({@code <name>.lib.parquet}) when needed. Writes (or, without a host, removes) {@code combinedDir/host_library.txt}.
	 */
	public static String[] finalSearchLibraries(File microbialLibrary, File hostParquet, File combinedDir,
			DiannParameter libraryBuild, File diannExe, ProgressCallback cb, File condaPrefix) throws Exception {
		return finalSearchLibraries(microbialLibrary, hostParquet == null ? List.of() : List.of(hostParquet),
				combinedDir, libraryBuild, diannExe, cb, condaPrefix);
	}

	public static String[] finalSearchLibraries(File microbialLibrary, List<File> hostParquets, File combinedDir,
			DiannParameter libraryBuild, File diannExe, ProgressCallback cb, File condaPrefix) throws Exception {
		return finalSearchLibraries(List.of(microbialLibrary), hostParquets, combinedDir, libraryBuild, diannExe, cb,
				condaPrefix);
	}

	/** Several microbial libraries (GUI Seed Direct may select more than one); each is converted as above. */
	public static String[] finalSearchLibraries(List<File> microbialLibraries, List<File> hostParquets,
			File combinedDir, DiannParameter libraryBuild, File diannExe, ProgressCallback cb, File condaPrefix)
			throws Exception {
		if (cb == null) cb = new ProgressCallback.NoOp();
		List<String> libs = new ArrayList<>();
		if (hostParquets == null || hostParquets.isEmpty()) {
			writeHostRecord(combinedDir, List.of());
			for (File micro : microbialLibraries) libs.add(micro.getAbsolutePath());
			return libs.toArray(new String[0]);
		}
		for (File host : hostParquets) {
			if (host == null || !host.isFile() || !host.getName().toLowerCase(Locale.ROOT).endsWith(".parquet")) {
				throw new IOException("Host library for the combined search must be an existing .parquet: "
						+ (host == null ? "null" : host.getAbsolutePath()));
			}
		}
		for (File micro : microbialLibraries) {
			libs.add(ensureParquetLibrary(micro, micro.getAbsoluteFile().getParentFile(), stripExtension(micro.getName()),
					libraryBuild, diannExe, cb, condaPrefix).getAbsolutePath());
		}
		writeHostRecord(combinedDir, hostParquets);
		for (File host : hostParquets) libs.add(host.getAbsolutePath());
		cb.onMessage("Combined search libraries: " + String.join(" + ", libs) + " (host: " + hostParquets.size() + ")");
		return libs.toArray(new String[0]);
	}

	public static void writeHostRecord(File combinedDir, List<File> hostParquets) throws IOException {
		File record = new File(combinedDir, HOST_RECORD);
		if (hostParquets == null || hostParquets.isEmpty()) {
			Files.deleteIfExists(record.toPath());
			return;
		}
		if (!combinedDir.isDirectory() && !combinedDir.mkdirs()) {
			throw new IOException("Could not create " + combinedDir.getAbsolutePath());
		}
		List<String> lines = new ArrayList<>();
		for (File host : hostParquets) lines.add(host.getAbsolutePath());
		Files.write(record.toPath(), lines, StandardCharsets.UTF_8);
	}

	/** Host .parquet libraries recorded by the final search in {@code combinedDir}; empty when none was searched. */
	public static List<File> readHostRecord(File combinedDir) throws IOException {
		File record = combinedDir == null ? null : new File(combinedDir, HOST_RECORD);
		List<File> hosts = new ArrayList<>();
		if (record == null || !record.isFile()) return hosts;
		for (String line : Files.readAllLines(record.toPath(), StandardCharsets.UTF_8)) {
			if (line.isBlank()) continue;
			File host = new File(line.trim());
			if (!host.isFile()) {
				throw new IOException("Host library recorded in " + record.getAbsolutePath() + " not found: "
						+ host.getAbsolutePath());
			}
			hosts.add(host);
		}
		return hosts;
	}

	/**
	 * For a library in the MAG {@code libraries/<name>/} folder: its {@code <name>.parquet} when the library resolved
	 * to {@code <name>.speclib} or {@code .tsv} (the library builder writes both), else the library itself. Only for
	 * library folders: in a result folder {@code combined.parquet} is the search output, not a library.
	 */
	public static File parquetInLibraryFolder(File library) {
		if (library == null || library.getName().toLowerCase(Locale.ROOT).endsWith(".parquet")) return library;
		File sibling = new File(library.getAbsoluteFile().getParentFile(), stripExtension(library.getName()) + ".parquet");
		return sibling.isFile() ? sibling : library;
	}

	static String stripExtension(String name) {
		String n = name;
		for (String ext : new String[] { ".predicted.speclib", ".speclib", ".lib.parquet", ".parquet", ".fasta", ".fa",
				".tsv" }) {
			if (n.toLowerCase(Locale.ROOT).endsWith(ext)) return n.substring(0, n.length() - ext.length());
		}
		return n;
	}

	// ------------------------------------------------------------------ export side

	/** Result of {@link #split}. */
	public static final class Split {
		public final DiaNNPrecursor[] microbial;
		public final DiaNNPrecursor[] host;

		Split(DiaNNPrecursor[] microbial, DiaNNPrecursor[] host) {
			this.microbial = microbial;
			this.host = host;
		}
	}

	/** Splits precursors by the host rule: stripped sequence in {@code hostPeptides} = host (host takes shared peptides). */
	public static Split split(DiaNNPrecursor[] precursors, Set<String> hostPeptides) {
		if (precursors == null) precursors = new DiaNNPrecursor[0];
		if (hostPeptides == null || hostPeptides.isEmpty()) return new Split(precursors, new DiaNNPrecursor[0]);
		List<DiaNNPrecursor> microbial = new ArrayList<>(precursors.length);
		List<DiaNNPrecursor> host = new ArrayList<>();
		for (DiaNNPrecursor p : precursors) {
			if (p == null) continue;
			String seq = p.getSeqString();
			if (seq != null && hostPeptides.contains(seq)) host.add(p);
			else microbial.add(p);
		}
		return new Split(microbial.toArray(new DiaNNPrecursor[0]), host.toArray(new DiaNNPrecursor[0]));
	}

	/**
	 * Export step for every DIA workflow: when {@code combinedDir/host_library.txt} exists, removes the host
	 * precursors (see class comment) and writes {@code host_peptides.tsv} and {@code host_summary.tsv} to
	 * {@code resultDir}. Without a host record the precursors are returned unchanged and stale host files are removed.
	 */
	public static DiaNNPrecursor[] removeHostPrecursors(DiaNNPrecursor[] precursors, File combinedDir, File resultDir,
			ProgressCallback cb) throws IOException {
		if (cb == null) cb = new ProgressCallback.NoOp();
		List<File> hostParquets = readHostRecord(combinedDir);
		File hostPeptidesTsv = new File(resultDir, HOST_PEPTIDES_TSV);
		File hostSummaryTsv = new File(resultDir, HOST_SUMMARY_TSV);
		if (hostParquets.isEmpty()) {
			Files.deleteIfExists(hostPeptidesTsv.toPath());
			Files.deleteIfExists(hostSummaryTsv.toPath());
			return precursors;
		}
		Set<String> hostPeptides = new HashSet<>();
		for (File hostParquet : hostParquets) {
			cb.onMessage("DIA post-analysis: reading host peptides from " + hostParquet.getAbsolutePath());
			hostPeptides.addAll(DiaNNParquetReader.getPeptideSet(hostParquet.getAbsolutePath()));
		}
		Split split = split(precursors, hostPeptides);
		cb.onMessage("DIA post-analysis: " + split.host.length + " host precursors removed (host takes shared peptides), "
				+ split.microbial.length + " microbial precursors kept");
		File prMatrix = new File(combinedDir, "combined.pr_matrix.tsv");
		writeHostReports(split, prMatrix, hostPeptidesTsv, hostSummaryTsv);
		return split.microbial;
	}

	static void writeHostReports(Split split, File prMatrix, File hostPeptidesTsv, File hostSummaryTsv)
			throws IOException {
		String[] runs = new String[0];
		FlashLfqQuanPeptide[] hostPeps = new FlashLfqQuanPeptide[0];
		FlashLfqQuanPeptide[] microPeps = new FlashLfqQuanPeptide[0];
		if (prMatrix.isFile()) {
			DiannPrMatrixReader hostReader = new DiannPrMatrixReader(prMatrix);
			hostPeps = nonNull(hostReader.getMetaPeptides(split.host));
			runs = hostReader.getIntensityTitle();
			DiannPrMatrixReader microReader = new DiannPrMatrixReader(prMatrix);
			microPeps = nonNull(microReader.getMetaPeptides(split.microbial));
			if (runs == null) runs = microReader.getIntensityTitle();
		}
		if (runs == null) runs = new String[0];

		// host peptides: one row per stripped sequence, intensities summed over its precursors
		Map<String, double[]> hostIntensity = new LinkedHashMap<>();
		Map<String, Set<String>> hostProteins = new LinkedHashMap<>();
		for (DiaNNPrecursor p : split.host) {
			String seq = p.getSeqString();
			Set<String> pros = hostProteins.computeIfAbsent(seq, k -> new TreeSet<>());
			if (p.getProGroups() != null) pros.addAll(Arrays.asList(p.getProGroups()));
		}
		for (FlashLfqQuanPeptide pep : hostPeps) {
			double[] in = pep.getIntensity();
			if (in == null) continue;
			double[] sum = hostIntensity.computeIfAbsent(pep.getSequence(), k -> new double[in.length]);
			for (int i = 0; i < Math.min(sum.length, in.length); i++) sum[i] += in[i];
		}
		try (PrintWriter w = new PrintWriter(hostPeptidesTsv, StandardCharsets.UTF_8)) {
			StringBuilder head = new StringBuilder("Sequence\tProteins");
			for (String r : runs) head.append("\tIntensity ").append(r);
			w.println(head);
			for (String seq : new TreeSet<>(hostProteins.keySet())) {
				StringBuilder sb = new StringBuilder(seq).append('\t').append(String.join(";", hostProteins.get(seq)));
				double[] in = hostIntensity.get(seq);
				for (int i = 0; i < runs.length; i++) sb.append('\t').append(in != null && i < in.length ? in[i] : 0.0);
				w.println(sb);
			}
		}

		int n = runs.length;
		int[] hostCount = new int[n], microCount = new int[n];
		double[] hostSum = new double[n], microSum = new double[n];
		Set<String>[] hostSeqs = newSets(n), microSeqs = newSets(n);
		accumulate(hostPeps, hostSeqs, hostSum);
		accumulate(microPeps, microSeqs, microSum);
		for (int i = 0; i < n; i++) {
			hostCount[i] = hostSeqs[i].size();
			microCount[i] = microSeqs[i].size();
		}
		try (PrintWriter w = new PrintWriter(hostSummaryTsv, StandardCharsets.UTF_8)) {
			w.println("Run\tHost peptides\tMicrobial peptides\tHost peptide fraction\tHost intensity\t"
					+ "Microbial intensity\tHost intensity fraction");
			for (int i = 0; i < n; i++) {
				int pt = hostCount[i] + microCount[i];
				double it = hostSum[i] + microSum[i];
				w.println(runs[i] + "\t" + hostCount[i] + "\t" + microCount[i] + "\t"
						+ (pt > 0 ? String.format(Locale.ROOT, "%.4f", (double) hostCount[i] / pt) : "0") + "\t"
						+ hostSum[i] + "\t" + microSum[i] + "\t"
						+ (it > 0 ? String.format(Locale.ROOT, "%.4f", hostSum[i] / it) : "0"));
			}
			w.println("# host precursors (all runs):\t" + split.host.length);
			w.println("# host peptides (all runs):\t" + hostProteins.size());
		}
	}

	private static void accumulate(FlashLfqQuanPeptide[] peps, Set<String>[] seqs, double[] sums) {
		for (FlashLfqQuanPeptide pep : peps) {
			double[] in = pep.getIntensity();
			if (in == null) continue;
			for (int i = 0; i < Math.min(in.length, sums.length); i++) {
				if (in[i] > 0) {
					sums[i] += in[i];
					seqs[i].add(pep.getSequence());
				}
			}
		}
	}

	@SuppressWarnings("unchecked")
	private static Set<String>[] newSets(int n) {
		Set<String>[] s = new Set[n];
		for (int i = 0; i < n; i++) s[i] = new HashSet<>();
		return s;
	}

	private static FlashLfqQuanPeptide[] nonNull(FlashLfqQuanPeptide[] peps) {
		if (peps == null) return new FlashLfqQuanPeptide[0];
		return Arrays.stream(peps).filter(p -> p != null && p.getSequence() != null).toArray(FlashLfqQuanPeptide[]::new);
	}
}
