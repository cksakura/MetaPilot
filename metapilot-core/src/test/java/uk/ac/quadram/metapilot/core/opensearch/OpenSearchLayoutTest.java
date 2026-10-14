package uk.ac.quadram.metapilot.core.opensearch;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
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
	void eachEngineConfigurationHasItsOwnFolderBesideTheSharedSearches() {
		assertEquals(new File(ROOT, "msfragger_only"), OpenSearchLayout.modeRoot(ROOT, ONLY));
		assertEquals(new File(ROOT, "msfragger_sage"), OpenSearchLayout.modeRoot(ROOT, SAGE));
		assertEquals(OpenSearchLayout.modeRoot(ROOT, ONLY), OpenSearchLayout.modeRoot(ROOT, null));
		assertEquals(OpenSearchLayout.modeRoot(ROOT, ONLY), OpenSearchLayout.modeRoot(ROOT, " "));
		assertEquals(Arrays.asList("msfragger_only", "msfragger_sage"), OpenSearchLayout.modeFolderNames());

		assertEquals(new File(ROOT, "search_results"), OpenSearchLayout.searchResultsRoot(ROOT));
		assertEquals(new File(new File(new File(ROOT, "search_results"), "s1"), "sample_db"),
				OpenSearchLayout.sampleDbDir(ROOT, "s1"));
	}

	@Test
	void legacyMigrationLeavesThisWorkflowsFoldersInPlace(@TempDir Path tmp) {
		File root = OpenSearchLayout.engineSearchRoot(tmp.toFile());
		File onlyMode = OpenSearchLayout.modeRoot(root, ONLY);
		File sageMode = OpenSearchLayout.modeRoot(root, SAGE);
		File legacySample = new File(root, "S1");
		for (File dir : new File[] { onlyMode, sageMode, legacySample }) {
			assertTrue(dir.mkdirs());
		}

		File searchResults = OpenSearchLayout.searchResultsRoot(root);

		assertTrue(onlyMode.isDirectory());
		assertTrue(sageMode.isDirectory());
		assertFalse(legacySample.exists());
		assertTrue(new File(searchResults, "S1").isDirectory());
	}

	@Test
	void pooledOutputsLiveInTheModeFolderAndDoNotCollide() {
		File onlyCombined = new File(OpenSearchLayout.modeRoot(ROOT, ONLY), "combined");
		assertEquals(onlyCombined, OpenSearchLayout.combinedDir(ROOT, ONLY));
		assertEquals(OpenSearchLayout.combinedDir(OpenSearchLayout.modeRoot(ROOT, SAGE)),
				OpenSearchLayout.combinedDir(ROOT, SAGE));
		assertEquals(new File(onlyCombined, "validated_psms.tsv"), OpenSearchLayout.validatedPsms(ROOT, ONLY));
		assertEquals(new File(onlyCombined, "variant_species_summary.tsv"),
				OpenSearchLayout.variantSpeciesSummary(ROOT, ONLY));
		assertNotEquals(OpenSearchLayout.validatedPsms(ROOT, ONLY), OpenSearchLayout.validatedPsms(ROOT, SAGE));
		assertNotEquals(OpenSearchLayout.quantifiedPeptides(ROOT, ONLY),
				OpenSearchLayout.quantifiedPeptides(ROOT, SAGE));

		File sageMode = OpenSearchLayout.modeRoot(ROOT, SAGE);
		assertEquals(new File(sageMode, "final_peptides.tsv"), OpenSearchLayout.finalPeptides(sageMode));
		assertEquals(new File(sageMode, "razor_proteins.tsv"), OpenSearchLayout.razorProteins(sageMode));
	}
}
