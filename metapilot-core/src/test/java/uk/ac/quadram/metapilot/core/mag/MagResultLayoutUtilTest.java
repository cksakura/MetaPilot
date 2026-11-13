package uk.ac.quadram.metapilot.core.mag;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class MagResultLayoutUtilTest {

	@TempDir
	Path tmp;

	@Test
	void nestedLayoutLeavesSpectraSummaryAtTheEngineRoot() throws Exception {
		File root = tmp.resolve("fragpipe_MAG_search").toFile();
		File summary = MagResultLayoutUtil.spectraSummaryDir(root);
		assertTrue(summary.mkdirs());
		Files.writeString(new File(summary, "ms2Count.txt").toPath(), "FileName\tMS2 count\nS1\t10\n");
		File sample = new File(root, "S1");
		assertTrue(sample.mkdirs());

		File searchResults = MagResultLayoutUtil.resolveSearchResultsRoot(root);

		assertTrue(summary.isDirectory(), "spectra_summary stays at the engine root");
		assertTrue(new File(summary, "ms2Count.txt").isFile());
		assertFalse(new File(new File(searchResults, MagResultLayoutUtil.SPECTRA_SUMMARY_DIR_NAME), "ms2Count.txt")
				.exists());
		assertTrue(new File(searchResults, "S1").isDirectory());
		assertFalse(sample.exists());
	}

	@Test
	void nestedLayoutLeavesExcludedEvidenceAtTheEngineRoot() throws Exception {
		File root = tmp.resolve("pfind_MAG_search").toFile();
		File evidence = new File(root, MagResultLayoutUtil.EXCLUDED_EVIDENCE_DIR_NAME);
		assertTrue(evidence.mkdirs());
		Files.writeString(new File(evidence, "peptides.tsv").toPath(),
				"Base sequence\tModified sequence\tOriginal proteins\tReason\n");
		File sample = new File(root, "S1");
		assertTrue(sample.mkdirs());

		File searchResults = MagResultLayoutUtil.resolveSearchResultsRoot(root);

		assertTrue(new File(evidence, "peptides.tsv").isFile());
		assertFalse(new File(searchResults, MagResultLayoutUtil.EXCLUDED_EVIDENCE_DIR_NAME).exists());
		assertTrue(new File(searchResults, "S1").isDirectory());
	}

	@Test
	void restoreMovesAPreviouslyNestedCopyBack() throws Exception {
		File root = tmp.resolve("opensearch_MAG_search").toFile();
		File nested = new File(new File(root, MagResultLayoutUtil.SEARCH_RESULTS_DIR_NAME),
				MagResultLayoutUtil.SPECTRA_SUMMARY_DIR_NAME);
		assertTrue(nested.mkdirs());
		Files.writeString(new File(nested, "ms2Count.txt").toPath(), "FileName\tMS2 count\nS1\t10\n");

		File restored = MagResultLayoutUtil.spectraCountDir(root);

		assertEquals(MagResultLayoutUtil.spectraDir(root), restored);
		assertTrue(new File(restored, "ms2Count.txt").isFile());
		assertFalse(new File(nested, "ms2Count.txt").exists());
		assertEquals(new File(restored, "ms2Count.txt"), MagResultLayoutUtil.findMs2Count(root));
	}

	@Test
	void findMs2CountReadsACopyStillUnderSearchResults() throws Exception {
		File root = tmp.resolve("engine").toFile();
		File nested = new File(new File(root, MagResultLayoutUtil.SEARCH_RESULTS_DIR_NAME),
				MagResultLayoutUtil.SPECTRA_SUMMARY_DIR_NAME);
		assertTrue(nested.mkdirs());
		File counts = new File(nested, "ms2Count.txt");
		Files.writeString(counts.toPath(), "FileName\tMS2 count\nS1\t10\n");

		assertEquals(counts, MagResultLayoutUtil.findMs2Count(root));
	}
}
