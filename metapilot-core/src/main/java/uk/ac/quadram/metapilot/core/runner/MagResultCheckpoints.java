package uk.ac.quadram.metapilot.core.runner;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import uk.ac.quadram.metapilot.core.opensearch.OpenSearchLayout;
import uk.ac.quadram.metapilot.task.io.qpx.QpxPsmParquetWriter;

/**
 * Make-style checkpoints for the shared MAG post-analysis under a result
 * folder: the final tables, taxonomy and functional annotation, and the HTML
 * reports built from them. An output is current when it exists, is non-empty
 * and is no older than what it is derived from, as in
 * {@link OpenSearchLayout#isSatisfied(File, File...)}; a missing input is
 * ignored.
 * <p>
 * Parameters are not part of the check: a run with changed post-analysis or
 * report settings but unchanged inputs needs the stale outputs deleted.
 */
public final class MagResultCheckpoints {

	public static final String FINAL_SUMMARY = "final_summary.tsv";
	public static final String FINAL_PEPTIDES = "final_peptides.tsv";
	public static final String RAZOR_PROTEINS = "razor_proteins.tsv";
	public static final String ALL_PROTEINS = "all_proteins.tsv";
	public static final String TAXONOMY_DIR = "taxonomy_analysis";
	public static final String FUNCTION_DIR = "functional_annotation";
	public static final String METADATA = "metainfo.tsv";
	/** Which inputs the post-analysis outputs under a result folder were produced from. */
	public static final String POST_ANALYSIS_SOURCE = "post_analysis_source.tsv";

	private MagResultCheckpoints() {
	}

	/**
	 * True when {@link #POST_ANALYSIS_SOURCE} names exactly {@code inputs}.
	 * Outputs newer than the inputs can still come from other ones: an export
	 * interrupted after some tables were rewritten, or a folder produced before
	 * this record existed. Paths are recorded relative to {@code resultRoot},
	 * so the same folder reached through another drive letter still matches.
	 */
	public static boolean sameSource(File resultRoot, File... inputs) {
		File record = new File(resultRoot, POST_ANALYSIS_SOURCE);
		if (!record.isFile()) {
			return false;
		}
		try {
			List<String> recorded = Files.readAllLines(record.toPath(), StandardCharsets.UTF_8);
			return describeSource(resultRoot, inputs).equals(recorded);
		} catch (IOException e) {
			return false;
		}
	}

	/** Records {@code inputs} as the source of the post-analysis outputs under {@code resultRoot}. */
	public static void recordSource(File resultRoot, File... inputs) throws IOException {
		Files.write(new File(resultRoot, POST_ANALYSIS_SOURCE).toPath(), describeSource(resultRoot, inputs),
				StandardCharsets.UTF_8);
	}

	/** Forgets the recorded source, so an interrupted regeneration is not mistaken for a finished one. */
	public static void clearSource(File resultRoot) throws IOException {
		Files.deleteIfExists(new File(resultRoot, POST_ANALYSIS_SOURCE).toPath());
	}

	static List<String> describeSource(File resultRoot, File... inputs) {
		Path root = resultRoot.toPath().toAbsolutePath().normalize();
		List<String> lines = new ArrayList<>();
		for (File input : inputs) {
			if (input == null) {
				continue;
			}
			Path path = input.toPath().toAbsolutePath().normalize();
			String name = path.startsWith(root) ? root.relativize(path).toString() : path.toString();
			lines.add("input\t" + name.replace('\\', '/'));
		}
		return lines;
	}

	/**
	 * The first post-analysis output that is missing or older than what it is
	 * derived from, or null when the tables, taxonomy and (if expected)
	 * functional annotation are all current.
	 *
	 * @param inputs the files the post-analysis reads, e.g. the PSM table and
	 *               the quantified peptides
	 */
	public static File stalePostAnalysisOutput(File resultRoot, boolean functionExpected, File... inputs) {
		File peptides = new File(resultRoot, FINAL_PEPTIDES);
		File razor = new File(resultRoot, RAZOR_PROTEINS);
		for (File table : new File[] { new File(resultRoot, FINAL_SUMMARY), peptides, razor,
				new File(resultRoot, ALL_PROTEINS) }) {
			if (!OpenSearchLayout.isSatisfied(table, inputs)) {
				return table;
			}
		}
		File taxonomy = new File(resultRoot, TAXONOMY_DIR);
		for (File out : new File[] { new File(taxonomy, "genome_final_selected.tsv"),
				new File(taxonomy, "Taxa_report.tsv") }) {
			if (!OpenSearchLayout.isSatisfied(out, peptides, razor)) {
				return out;
			}
		}
		if (functionExpected) {
			File annotated = new File(new File(resultRoot, FUNCTION_DIR), "Annotated_proteins.tsv");
			if (!OpenSearchLayout.isSatisfied(annotated, razor)) {
				return annotated;
			}
		}
		return null;
	}

	/**
	 * True when the QPX PSM parquet beside the final tables exists and is no
	 * older than the PSM table it is written from.
	 */
	public static boolean qpxPsmsCurrent(File resultRoot, File psmTable) {
		return OpenSearchLayout.isSatisfied(new File(resultRoot, QpxPsmParquetWriter.DEFAULT_FILE_NAME), psmTable);
	}

	/**
	 * The first HTML report under {@code report/reports} that is missing or
	 * older than the table it shows or {@code metadata}, or null when every
	 * report with a source present is current. Mirrors the reports
	 * {@code MetaPilotMagTask.generateFinalReports} writes.
	 */
	public static File staleReport(File resultRoot, File metadata) {
		File reports = new File(new File(resultRoot, "report"), "reports");
		File taxonomy = new File(resultRoot, TAXONOMY_DIR);
		File function = new File(resultRoot, FUNCTION_DIR);
		File[][] sourceAndReport = {
				{ new File(resultRoot, FINAL_SUMMARY), new File(reports, "report_ID_summary.html") },
				{ new File(resultRoot, FINAL_PEPTIDES), new File(reports, "report_peptides_summary.html") },
				{ new File(resultRoot, RAZOR_PROTEINS), new File(reports, "report_proteinGroups_summary.html") },
				{ new File(taxonomy, "Taxa_report.tsv"), new File(reports, "report_taxonomy_summary.html") },
				{ new File(function, "Annotated_proteins.tsv"), new File(reports, "report_function_summary.html") } };
		for (File[] pair : sourceAndReport) {
			if (OpenSearchLayout.isSatisfied(pair[0]) && !OpenSearchLayout.isSatisfied(pair[1], pair[0], metadata)) {
				return pair[1];
			}
		}
		File[][] copies = { { new File(taxonomy, "Genome.html"), new File(reports, "Genome.html") },
				{ new File(function, "functions.html"), new File(reports, "functions.html") } };
		for (File[] pair : copies) {
			if (pair[0].isFile() && !OpenSearchLayout.isSatisfied(pair[1], pair[0])) {
				return pair[1];
			}
		}
		File manifest = new File(reports, "manifest.json");
		return OpenSearchLayout.isSatisfied(manifest) ? null : manifest;
	}

	/** True when both files exist and hold the same bytes. */
	public static boolean sameContent(File a, File b) throws IOException {
		if (!a.isFile() || !b.isFile() || a.length() != b.length()) {
			return false;
		}
		return Arrays.equals(digest(a), digest(b));
	}

	static byte[] digest(File file) throws IOException {
		MessageDigest md;
		try {
			md = MessageDigest.getInstance("SHA-256");
		} catch (NoSuchAlgorithmException e) {
			throw new IllegalStateException(e);
		}
		byte[] buffer = new byte[1 << 20];
		try (InputStream in = new BufferedInputStream(Files.newInputStream(file.toPath()), 1 << 20)) {
			int n;
			while ((n = in.read(buffer)) > 0) {
				md.update(buffer, 0, n);
			}
		}
		return md.digest();
	}

	/**
	 * A file a step is about to rewrite. When the rewrite produces the same
	 * bytes, {@link #restoreIfUnchanged()} puts the previous timestamp back so
	 * the steps derived from the file stay current and are not redone.
	 * <p>
	 * Only a file that was current against the step's inputs is remembered;
	 * restoring the timestamp of a stale one would leave it looking stale and
	 * redo the step on every run.
	 */
	public static final class UnchangedRewrite {
		private final File file;
		private final long lastModified;
		private final byte[] digest;

		private UnchangedRewrite(File file, long lastModified, byte[] digest) {
			this.file = file;
			this.lastModified = lastModified;
			this.digest = digest;
		}

		/** Remembers {@code file} if it is current against {@code inputs}. */
		public static UnchangedRewrite before(File file, File... inputs) {
			if (!OpenSearchLayout.isSatisfied(file, inputs)) {
				return new UnchangedRewrite(file, 0L, null);
			}
			try {
				return new UnchangedRewrite(file, file.lastModified(), digest(file));
			} catch (IOException e) {
				return new UnchangedRewrite(file, 0L, null);
			}
		}

		/** True when the previous version was remembered. */
		public boolean remembered() {
			return digest != null;
		}

		/**
		 * True when the file was remembered, still holds the same bytes and its
		 * previous timestamp was put back.
		 */
		public boolean restoreIfUnchanged() {
			if (digest == null || !file.isFile()) {
				return false;
			}
			try {
				return Arrays.equals(digest, digest(file)) && file.setLastModified(lastModified);
			} catch (IOException e) {
				return false;
			}
		}
	}
}
