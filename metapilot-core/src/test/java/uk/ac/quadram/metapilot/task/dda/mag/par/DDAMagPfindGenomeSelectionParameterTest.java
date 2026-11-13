package uk.ac.quadram.metapilot.task.dda.mag.par;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.nio.file.Files;
import java.nio.file.Path;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import uk.ac.quadram.metapilot.task.mag.GenomeSelectionCut;

class DDAMagPfindGenomeSelectionParameterTest {

	@TempDir
	Path tempDir;

	@Test
	void oldFileWithoutTheNewKeysKeepsFocus() throws Exception {
		Path file = tempDir.resolve("old.json");
		Files.writeString(file, "{\"additional_parameters\":{\"genome_selection_strategy\":\"Focus\"}}");
		DDAMagPfindParameters params = new DDAMagPfindParameterIO().loadParameters(file.toFile());
		assertEquals(GenomeSelectionCut.STRATEGY_FOCUS, params.getGenomeSelectionStrategy());
		assertEquals(GenomeSelectionCut.DEFAULT_PSM_COVERAGE, params.getGenomeSelectionPsmCoverage(), 1e-12);
		assertEquals(GenomeSelectionCut.DEFAULT_STEP_GUARD, params.getGenomeSelectionStepGuard(), 1e-12);
		assertEquals(GenomeSelectionCut.DEFAULT_MAX_GENOMES, params.getGenomeSelectionMaxGenomes());
	}

	@Test
	void savedFileRestoresTheCoverageCut(@TempDir Path dir) throws Exception {
		DDAMagPfindParameters params = new DDAMagPfindParameters();
		params.setGenomeSelectionStrategy(GenomeSelectionCut.STRATEGY_PSM_COVERAGE);
		params.setGenomeSelectionPsmCoverage(0.99);
		params.setGenomeSelectionStepGuard(0);
		params.setGenomeSelectionMaxGenomes(500);
		DDAMagPfindParameterIO io = new DDAMagPfindParameterIO();
		assertTrue(io.saveAnalysisParameters(dir.toString(), params));
		DDAMagPfindParameters loaded = io.loadParameters(dir.resolve("mag_pfind_analysis_parameters.json").toFile());
		assertEquals(GenomeSelectionCut.STRATEGY_PSM_COVERAGE, loaded.getGenomeSelectionStrategy());
		assertEquals(0.99, loaded.getGenomeSelectionPsmCoverage(), 1e-12);
		assertEquals(0.0, loaded.getGenomeSelectionStepGuard(), 1e-12);
		assertEquals(500, loaded.getGenomeSelectionMaxGenomes());
	}

	@Test
	void outOfRangeCoverageBecomesTheDefault() throws Exception {
		Path file = tempDir.resolve("bad.json");
		Files.writeString(file, "{\"additional_parameters\":{\"genome_selection_psm_coverage\":1.5,"
				+ "\"genome_selection_step_guard\":0.01,\"genome_selection_max_genomes\":0,"
				+ "\"genome_selection_strategy\":\"not-a-mode\"}}");
		DDAMagPfindParameters loaded = new DDAMagPfindParameterIO().loadParameters(file.toFile());
		assertEquals(GenomeSelectionCut.STRATEGY_FOCUS, loaded.getGenomeSelectionStrategy());
		assertEquals(GenomeSelectionCut.DEFAULT_PSM_COVERAGE, loaded.getGenomeSelectionPsmCoverage(), 1e-12);
		assertEquals(GenomeSelectionCut.DEFAULT_STEP_GUARD, loaded.getGenomeSelectionStepGuard(), 1e-12);
		assertEquals(GenomeSelectionCut.DEFAULT_MAX_GENOMES, loaded.getGenomeSelectionMaxGenomes());
	}
}
