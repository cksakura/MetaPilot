package uk.ac.quadram.metapilot.core.opensearch;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import uk.ac.quadram.metapilot.task.dda.mag.par.DDAMagOpenSearchParameters;

class OpenSearchLayoutTest {

	private static final File ROOT = OpenSearchLayout.engineSearchRoot(new File("results"));
	private static final String ONLY = DDAMagOpenSearchParameters.ENGINES_MSFRAGGER;
	private static final String SAGE = DDAMagOpenSearchParameters.ENGINES_MSFRAGGER_AND_SAGE;

	@Test
	void deliverablesSitInTheEngineRootLikeOtherMagWorkflows() {
		assertEquals(ROOT, OpenSearchLayout.modeRoot(ROOT, ONLY));
		assertEquals(ROOT, OpenSearchLayout.modeRoot(ROOT, SAGE));
		assertEquals(ROOT, OpenSearchLayout.modeRoot(ROOT, null));
		assertEquals(Arrays.asList("msfragger_only", "msfragger_sage"), OpenSearchLayout.modeFolderNames());

		assertEquals(new File(new File(ROOT, "combined"), "validated_psms.tsv"),
				OpenSearchLayout.validatedPsms(ROOT, ONLY));
		assertEquals(OpenSearchLayout.validatedPsms(ROOT, ONLY), OpenSearchLayout.validatedPsms(ROOT, SAGE));
		assertEquals(new File(ROOT, "final_peptides.tsv"), OpenSearchLayout.finalPeptides(ROOT));
	}

	@Test
	void sequenceVariantTablesLeaveCombined(@TempDir Path tmp) throws Exception {
		File root = OpenSearchLayout.engineSearchRoot(tmp.toFile());
		File combined = new File(root, "combined");
		assertTrue(combined.mkdirs());
		Files.writeString(new File(combined, "validated_psms.tsv").toPath(), "scan\n");
		Files.writeString(new File(combined, "sequence_variants.tsv").toPath(), "site\n");
		Files.writeString(new File(combined, "variant_candidates.tsv").toPath(), "psm\n");

		File table = OpenSearchLayout.sequenceVariants(root, ONLY);
		File candidates = OpenSearchLayout.variantCandidates(root, ONLY);
		File searchResults = OpenSearchLayout.searchResultsRoot(root);

		assertEquals(new File(new File(root, "sequence_variants"), "sequence_variants.tsv"), table);
		assertTrue(table.isFile());
		assertTrue(candidates.isFile());
		assertFalse(new File(combined, "sequence_variants.tsv").exists());
		assertFalse(new File(combined, "variant_candidates.tsv").exists());
		assertTrue(new File(combined, "validated_psms.tsv").isFile());
		assertFalse(new File(searchResults, "sequence_variants").exists());
	}

	@Test
	void leftoverSageFolderIsNotTreatedAsASample(@TempDir Path tmp) throws Exception {
		File root = OpenSearchLayout.engineSearchRoot(tmp.toFile());
		File sageMode = new File(root, OpenSearchLayout.LEGACY_MSFRAGGER_SAGE_DIR);
		assertTrue(sageMode.mkdirs());
		Files.writeString(new File(sageMode, "keep.txt").toPath(), "sage\n");
		File legacySample = new File(root, "S1");
		assertTrue(legacySample.mkdirs());

		File searchResults = OpenSearchLayout.searchResultsRoot(root);

		assertTrue(new File(sageMode, "keep.txt").isFile());
		assertFalse(legacySample.exists());
		assertTrue(new File(searchResults, "S1").isDirectory());
	}

	@Test
	void nestedLayoutLeavesSpectraSummaryBesideCombined(@TempDir Path tmp) throws Exception {
		File root = OpenSearchLayout.engineSearchRoot(tmp.toFile());
		File summary = new File(root, uk.ac.quadram.metapilot.core.mag.MagResultLayoutUtil.SPECTRA_SUMMARY_DIR_NAME);
		assertTrue(summary.mkdirs());
		java.nio.file.Files.writeString(new File(summary, "ms2Count.txt").toPath(),
				"FileName\tMS2 count\nS1\t10\n");
		File legacySample = new File(root, "S1");
		assertTrue(legacySample.mkdirs());

		File searchResults = OpenSearchLayout.searchResultsRoot(root);

		assertTrue(summary.isDirectory());
		assertTrue(new File(summary, "ms2Count.txt").isFile());
		assertFalse(legacySample.exists());
		assertTrue(new File(searchResults, "S1").isDirectory());
	}

	@Test
	void liftsMsfraggerOnlyIntoTheEngineRoot(@TempDir Path tmp) throws Exception {
		File root = OpenSearchLayout.engineSearchRoot(tmp.toFile());
		File legacy = new File(root, OpenSearchLayout.LEGACY_MSFRAGGER_ONLY_DIR);
		File combined = new File(legacy, "combined");
		assertTrue(combined.mkdirs());
		File validated = new File(combined, "validated_psms.tsv");
		Files.writeString(validated.toPath(), "file\tscan\n");
		File peptides = new File(legacy, "final_peptides.tsv");
		Files.writeString(peptides.toPath(), "Base sequence\n");

		OpenSearchLayout.migrateLegacyModeRoot(root);

		assertTrue(new File(new File(root, "combined"), "validated_psms.tsv").isFile());
		assertTrue(new File(root, "final_peptides.tsv").isFile());
		assertFalse(legacy.exists());
		assertEquals(root, OpenSearchLayout.modeRoot(root, ONLY));
	}

	@Test
	void doesNotLiftWhenTheNewCombinedAlreadyExists(@TempDir Path tmp) throws Exception {
		File root = OpenSearchLayout.engineSearchRoot(tmp.toFile());
		File newCombined = new File(root, "combined");
		assertTrue(newCombined.mkdirs());
		Files.writeString(new File(newCombined, "validated_psms.tsv").toPath(), "new\n");
		File legacy = new File(new File(root, OpenSearchLayout.LEGACY_MSFRAGGER_ONLY_DIR), "combined");
		assertTrue(legacy.mkdirs());
		Files.writeString(new File(legacy, "validated_psms.tsv").toPath(), "old\n");

		OpenSearchLayout.migrateLegacyModeRoot(root);

		assertEquals("new\n", Files.readString(new File(newCombined, "validated_psms.tsv").toPath()));
		assertTrue(legacy.isDirectory());
	}
}
