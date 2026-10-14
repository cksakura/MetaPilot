package uk.ac.quadram.metapilot.core.runner;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class MagResultCheckpointsTest {

	private static final long T0 = 1_700_000_000_000L;

	@TempDir
	Path tmp;

	@Test
	void postAnalysisCurrentOnlyWhenEveryOutputIsNoOlderThanItsInputs() throws IOException {
		File root = tmp.toFile();
		File validated = write(new File(root, "validated_psms.tsv"), "v", T0);
		File quantified = write(new File(root, "QuantifiedPeptides.tsv"), "q", T0 + 1000);
		writePostAnalysis(root, T0 + 2000, true);

		assertNull(MagResultCheckpoints.stalePostAnalysisOutput(root, true, validated, quantified));

		quantified.setLastModified(T0 + 5000);
		assertEquals(MagResultCheckpoints.FINAL_SUMMARY,
				MagResultCheckpoints.stalePostAnalysisOutput(root, true, validated, quantified).getName());

		quantified.setLastModified(T0 + 1000);
		File taxa = new File(new File(root, MagResultCheckpoints.TAXONOMY_DIR), "Taxa_report.tsv");
		taxa.delete();
		assertEquals(taxa, MagResultCheckpoints.stalePostAnalysisOutput(root, true, validated, quantified));
	}

	@Test
	void functionalAnnotationIsOnlyRequiredWhenExpected() throws IOException {
		File root = tmp.toFile();
		File validated = write(new File(root, "validated_psms.tsv"), "v", T0);
		writePostAnalysis(root, T0 + 2000, false);

		assertNull(MagResultCheckpoints.stalePostAnalysisOutput(root, false, validated, null));
		assertEquals("Annotated_proteins.tsv",
				MagResultCheckpoints.stalePostAnalysisOutput(root, true, validated, null).getName());
	}

	@Test
	void reportStaleWhenOlderThanItsTableOrTheMetadata() throws IOException {
		File root = tmp.toFile();
		writePostAnalysis(root, T0, true);
		File metadata = write(new File(root, MagResultCheckpoints.METADATA), "m", T0);
		File reports = new File(new File(root, "report"), "reports");
		for (String name : new String[] { "report_ID_summary.html", "report_peptides_summary.html",
				"report_proteinGroups_summary.html", "report_taxonomy_summary.html", "report_function_summary.html",
				"manifest.json" }) {
			write(new File(reports, name), "r", T0 + 1000);
		}

		assertNull(MagResultCheckpoints.staleReport(root, metadata));

		new File(root, MagResultCheckpoints.FINAL_PEPTIDES).setLastModified(T0 + 2000);
		assertEquals("report_peptides_summary.html", MagResultCheckpoints.staleReport(root, metadata).getName());

		new File(root, MagResultCheckpoints.FINAL_PEPTIDES).setLastModified(T0);
		metadata.setLastModified(T0 + 2000);
		assertEquals("report_ID_summary.html", MagResultCheckpoints.staleReport(root, metadata).getName());

		metadata.setLastModified(T0);
		write(new File(new File(root, MagResultCheckpoints.TAXONOMY_DIR), "Genome.html"), "g", T0 + 2000);
		assertEquals("Genome.html", MagResultCheckpoints.staleReport(root, metadata).getName());
		write(new File(reports, "Genome.html"), "g", T0 + 2000);
		assertNull(MagResultCheckpoints.staleReport(root, metadata));

		new File(reports, "manifest.json").delete();
		assertEquals("manifest.json", MagResultCheckpoints.staleReport(root, metadata).getName());
	}

	@Test
	void unchangedRewriteKeepsTheTimestampAndAChangedOneDoesNot() throws IOException {
		File root = tmp.toFile();
		File input = write(new File(root, "input.pepXML"), "x", T0);
		File out = write(new File(root, "validated_psms.tsv"), "same", T0 + 1000);

		MagResultCheckpoints.UnchangedRewrite rewrite = MagResultCheckpoints.UnchangedRewrite.before(out, input);
		assertTrue(rewrite.remembered());
		write(out, "same", T0 + 9000);
		assertTrue(rewrite.restoreIfUnchanged());
		assertEquals(T0 + 1000, out.lastModified());

		rewrite = MagResultCheckpoints.UnchangedRewrite.before(out, input);
		write(out, "different", T0 + 9000);
		assertFalse(rewrite.restoreIfUnchanged());
		assertEquals(T0 + 9000, out.lastModified());
	}

	@Test
	void staleFileIsNotRemembered() throws IOException {
		File root = tmp.toFile();
		File out = write(new File(root, "validated_psms.tsv"), "same", T0);
		File input = write(new File(root, "input.pepXML"), "x", T0 + 1000);

		MagResultCheckpoints.UnchangedRewrite rewrite = MagResultCheckpoints.UnchangedRewrite.before(out, input);
		assertFalse(rewrite.remembered());
		write(out, "same", T0 + 9000);
		assertFalse(rewrite.restoreIfUnchanged());
		assertEquals(T0 + 9000, out.lastModified());
	}

	@Test
	void sameSourceOnlyForTheRecordedEngineInputs() throws IOException {
		File root = tmp.toFile();
		File fraggerValidated = write(new File(root, "combined_msfragger_only/validated_psms.tsv"), "v", T0);
		File fraggerQuantified = write(new File(root, "combined_msfragger_only/QuantifiedPeptides.tsv"), "q", T0);
		File sageValidated = write(new File(root, "combined_msfragger_sage/validated_psms.tsv"), "v", T0);
		File sageQuantified = write(new File(root, "combined_msfragger_sage/QuantifiedPeptides.tsv"), "q", T0);

		assertFalse(MagResultCheckpoints.sameSource(root, sageValidated, sageQuantified));

		MagResultCheckpoints.recordSource(root, sageValidated, sageQuantified);
		assertTrue(MagResultCheckpoints.sameSource(root, sageValidated, sageQuantified));
		assertFalse(MagResultCheckpoints.sameSource(root, fraggerValidated, fraggerQuantified));
		assertFalse(MagResultCheckpoints.sameSource(root, sageValidated, null));

		MagResultCheckpoints.clearSource(root);
		assertFalse(MagResultCheckpoints.sameSource(root, sageValidated, sageQuantified));
	}

	@Test
	void recordedSourceIsRelativeToTheResultFolder() throws IOException {
		File first = new File(tmp.toFile(), "Y");
		File second = new File(tmp.toFile(), "Z");
		File validated = write(new File(first, "combined_msfragger_sage/validated_psms.tsv"), "v", T0);
		MagResultCheckpoints.recordSource(first, validated, null);

		File moved = write(new File(second, "combined_msfragger_sage/validated_psms.tsv"), "v", T0);
		Files.copy(new File(first, MagResultCheckpoints.POST_ANALYSIS_SOURCE).toPath(),
				new File(second, MagResultCheckpoints.POST_ANALYSIS_SOURCE).toPath());
		assertTrue(MagResultCheckpoints.sameSource(second, moved, null));
	}

	@Test
	void qpxPsmsCurrentOnlyWhenPresentAndNoOlderThanThePsmTable() throws IOException {
		File root = tmp.toFile();
		File psms = write(new File(root, "combined/psms.tsv"), "p", T0 + 1000);
		assertFalse(MagResultCheckpoints.qpxPsmsCurrent(root, psms));

		File parquet = write(new File(root, "metapilot.psm.parquet"), "q", T0);
		assertFalse(MagResultCheckpoints.qpxPsmsCurrent(root, psms));

		parquet.setLastModified(T0 + 2000);
		assertTrue(MagResultCheckpoints.qpxPsmsCurrent(root, psms));
	}

	@Test
	void sameContentComparesBytes() throws IOException {
		File root = tmp.toFile();
		File a = write(new File(root, "a.tsv"), "abc", T0);
		File b = write(new File(root, "b.tsv"), "abc", T0 + 1000);
		File c = write(new File(root, "c.tsv"), "abd", T0);
		assertTrue(MagResultCheckpoints.sameContent(a, b));
		assertFalse(MagResultCheckpoints.sameContent(a, c));
		assertFalse(MagResultCheckpoints.sameContent(a, new File(root, "missing.tsv")));
	}

	private static void writePostAnalysis(File root, long time, boolean function) throws IOException {
		for (String name : new String[] { MagResultCheckpoints.FINAL_SUMMARY, MagResultCheckpoints.FINAL_PEPTIDES,
				MagResultCheckpoints.RAZOR_PROTEINS, MagResultCheckpoints.ALL_PROTEINS }) {
			write(new File(root, name), name, time);
		}
		File taxonomy = new File(root, MagResultCheckpoints.TAXONOMY_DIR);
		write(new File(taxonomy, "genome_final_selected.tsv"), "s", time);
		write(new File(taxonomy, "Taxa_report.tsv"), "t", time);
		if (function) {
			write(new File(new File(root, MagResultCheckpoints.FUNCTION_DIR), "Annotated_proteins.tsv"), "a", time);
		}
	}

	private static File write(File file, String content, long time) throws IOException {
		file.getParentFile().mkdirs();
		Files.write(file.toPath(), content.getBytes(StandardCharsets.UTF_8));
		file.setLastModified(time);
		return file;
	}
}
