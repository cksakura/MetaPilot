package uk.ac.quadram.metapilot.task.dda.mag.par;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.nio.file.Files;
import java.nio.file.Path;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import uk.ac.quadram.metapilot.core.opensearch.OpenSearchLayout;
import uk.ac.quadram.metapilot.task.mag.GenomeSelectionCut;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;

class DDAMagOpenSearchGenomeSelectionParameterTest {

	@TempDir
	Path tempDir;

	@Test
	void oldFileWithoutTheNewKeysKeepsFocus() throws Exception {
		Path file = tempDir.resolve("old.json");
		Files.writeString(file, "{\"additional_parameters\":{\"genome_selection_strategy\":\"Focus\"}}");
		DDAMagOpenSearchParameters params = new DDAMagOpenSearchParameterIO().loadParameters(file.toFile());
		assertEquals(GenomeSelectionCut.STRATEGY_FOCUS, params.getGenomeSelectionStrategy());
		assertEquals(GenomeSelectionCut.DEFAULT_PSM_COVERAGE, params.getGenomeSelectionPsmCoverage(), 1e-12);
		assertEquals(GenomeSelectionCut.DEFAULT_STEP_GUARD, params.getGenomeSelectionStepGuard(), 1e-12);
		assertEquals(GenomeSelectionCut.DEFAULT_MAX_GENOMES, params.getGenomeSelectionMaxGenomes());
	}

	@Test
	void savedFileRestoresTheCoverageCutAndTheRunnerSeesIt() throws Exception {
		DDAMagOpenSearchParameters params = new DDAMagOpenSearchParameters();
		params.setResult(tempDir.toString());
		params.setEnzyme("Trypsin/P");
		params.setMissCleavages(1);
		params.setThreadCount(6);
		params.setMsfraggerMaxHeap("16g");
		params.setGenomeSelectorMode(MagDbItem.MARKER4);
		params.setGenomeSelectionStrategy(GenomeSelectionCut.STRATEGY_PSM_COVERAGE);
		params.setGenomeSelectionPsmCoverage(0.99);
		params.setGenomeSelectionStepGuard(0);
		params.setGenomeSelectionMaxGenomes(500);
		params.setEntrapmentEnabled(true);
		DDAMagOpenSearchParameterIO io = new DDAMagOpenSearchParameterIO();
		assertTrue(io.saveAnalysisParameters(tempDir.toString(), params));

		DDAMagOpenSearchParameters loaded = io
				.loadParameters(tempDir.resolve(OpenSearchLayout.PARAMETER_FILE).toFile());
		assertEquals("Trypsin/P", loaded.getEnzyme());
		assertEquals(1, loaded.getMissCleavages());
		assertEquals(6, loaded.getThreadCount());
		assertEquals("16g", loaded.getMsfraggerMaxHeap());
		assertEquals(MagDbItem.MARKER4, loaded.getGenomeSelectorMode());
		assertEquals(GenomeSelectionCut.STRATEGY_PSM_COVERAGE, loaded.getGenomeSelectionStrategy());
		assertEquals(0.99, loaded.getGenomeSelectionPsmCoverage(), 1e-12);
		assertEquals(0.0, loaded.getGenomeSelectionStepGuard(), 1e-12);
		assertEquals(500, loaded.getGenomeSelectionMaxGenomes());
		assertTrue(loaded.isEntrapmentEnabled());
		assertFalse(GenomeSelectionCut.isPsmCoverage(new DDAMagOpenSearchParameters().getGenomeSelectionStrategy()));
	}

	@Test
	void outOfRangeCoverageBecomesTheDefault() throws Exception {
		Path file = tempDir.resolve("bad.json");
		Files.writeString(file, "{\"additional_parameters\":{\"genome_selection_psm_coverage\":1.5}}");
		DDAMagOpenSearchParameters loaded = new DDAMagOpenSearchParameterIO().loadParameters(file.toFile());
		assertEquals(GenomeSelectionCut.DEFAULT_PSM_COVERAGE, loaded.getGenomeSelectionPsmCoverage(), 1e-12);
	}
}
