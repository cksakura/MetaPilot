package uk.ac.quadram.metapilot.task.mag;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.List;

import org.json.JSONObject;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class GenomeSelectionCutTest {

	@Test
	void seedsAlreadyAboveTheThresholdAreKept() {
		double[] coverage = { 0.92, 0.93, 0.94 };
		int[] increment = { 10, 4, 3 };
		GenomeSelectionCut.Decision decision = GenomeSelectionCut.decide(coverage, increment, 1000, 3, 0.90, 0.001, 300);
		assertEquals(3, decision.selectedCount);
		assertEquals(GenomeSelectionCut.STOP_SEEDS_ABOVE_THRESHOLD, decision.stopReason);
		assertEquals(0.94, decision.coverageReached, 1e-9);
	}

	@Test
	void guardIsTheLastIndexAboveEitherLimit() {
		// Increments are not monotonic. The guard is the earlier of the two last indexes,
		// not the first genome that falls below a limit.
		double[] coverage = { 0.10, 0.20, 0.30, 0.40, 0.50 };
		int[] increment = { 5, 0, 4, 0, 2 };
		GenomeSelectionCut.Decision decision = GenomeSelectionCut.decide(coverage, increment, 1000, 1, 0.99, 0.003, 300);
		assertEquals(3, decision.selectedCount);
		assertEquals(GenomeSelectionCut.STOP_STEP_GUARD, decision.stopReason);
	}

	@Test
	void coverageBeatsTheGuardAtTheSameIndex() {
		double[] coverage = { 0.10, 0.20, 0.96 };
		int[] increment = { 5, 0, 4 };
		GenomeSelectionCut.Decision decision = GenomeSelectionCut.decide(coverage, increment, 1000, 1, 0.95, 0.003, 300);
		assertEquals(3, decision.selectedCount);
		assertEquals(GenomeSelectionCut.STOP_COVERAGE, decision.stopReason);
	}

	@Test
	void zeroGuardDoesNotStopEarly() {
		double[] coverage = { 0.10, 0.20 };
		int[] increment = { 1, 0 };
		GenomeSelectionCut.Decision decision = GenomeSelectionCut.decide(coverage, increment, 1000, 1, 0.99, 0.0, 300);
		assertEquals(2, decision.selectedCount);
		assertEquals(GenomeSelectionCut.STOP_END_OF_RANKING, decision.stopReason);
	}

	@Test
	void genomeCapStopsWhenTheTargetIsStillAhead() {
		double[] coverage = { 0.10, 0.20, 0.30, 0.40 };
		int[] increment = { 5, 5, 5, 5 };
		GenomeSelectionCut.Decision decision = GenomeSelectionCut.decide(coverage, increment, 1000, 1, 0.99, 0.0, 3);
		assertEquals(3, decision.selectedCount);
		assertEquals(GenomeSelectionCut.STOP_GENOME_CAP, decision.stopReason);
	}

	@Test
	void outOfRangeValuesAreRejected() {
		assertThrows(IllegalArgumentException.class, () -> GenomeSelectionCut.requireCoverage(0));
		assertThrows(IllegalArgumentException.class, () -> GenomeSelectionCut.requireCoverage(1.1));
		assertThrows(IllegalArgumentException.class, () -> GenomeSelectionCut.requireStepGuard(-0.001));
		assertThrows(IllegalArgumentException.class, () -> GenomeSelectionCut.requireStepGuard(0.01));
		assertThrows(IllegalArgumentException.class, () -> GenomeSelectionCut.requireStepGuard(0.0021));
		GenomeSelectionCut.requireStepGuard(0.002);
		assertEquals(10, GenomeSelectionCut.stepGuardPercentSteps(0.001));
		assertEquals(0.001, GenomeSelectionCut.stepGuardFromPercentSteps(10), 1e-12);
		assertEquals("0.1%", GenomeSelectionCut.formatStepGuardPercent(10));
		assertEquals("0.01%", GenomeSelectionCut.formatStepGuardPercent(1));
		assertEquals("0.2%", GenomeSelectionCut.formatStepGuardPercent(20));
		assertEquals(10, GenomeSelectionCut.parseStepGuardPercent("0.1%"));
		assertEquals(1, GenomeSelectionCut.parseStepGuardPercent("0.01"));
		assertThrows(IllegalArgumentException.class, () -> GenomeSelectionCut.requireMaxGenomes(0));
		assertThrows(IllegalArgumentException.class, () -> GenomeSelectionCut.requireStrategy("MetaLab"));
		assertEquals(GenomeSelectionCut.STRATEGY_PSM_COVERAGE, GenomeSelectionCut.requireStrategy("psm coverage"));
	}

	@Test
	void missedTargetMessages() {
		assertEquals(
				"Reached 96.5% of the 99% target: stopped by the step guard. Lower or disable it to continue.",
				GenomeSelectionCut.missedTargetMessage(0.99, 0.965, GenomeSelectionCut.STOP_STEP_GUARD, 300));
		assertEquals("Stopped at Max genomes (300) with 98.6% coverage.",
				GenomeSelectionCut.missedTargetMessage(0.99, 0.986, GenomeSelectionCut.STOP_GENOME_CAP, 300));
		assertEquals("", GenomeSelectionCut.missedTargetMessage(0.95, 0.951, GenomeSelectionCut.STOP_COVERAGE, 300));
	}

	@Test
	void curveSkipsGenomesListedAfterTheCap(@TempDir File dir) throws Exception {
		File ranked = new File(dir, "genome_ranked.tsv");
		Files.writeString(ranked.toPath(), String.join("\n",
				"Rank\tGenome\tReason\tGreedyGain\tPsmCoverage\tCutStep\tInRanking",
				"0\tseed\tPRIVILEGE\t0\t0.80\t1\t1",
				"1\tg2\tGREEDY_COVERAGE\t4\t0.95\t1\t1",
				"2\tafter\t\t0\t0.99\t0\t0"), StandardCharsets.UTF_8);
		List<GenomeSelectionSummary.CurvePoint> curve = GenomeSelectionSummary.readRankingCurve(ranked);
		assertEquals(2, curve.size());
		assertEquals(1, curve.get(0).genomes);
		assertEquals(2, curve.get(1).genomes);
		assertTrue(curve.get(1).selected);
		assertEquals(0.95, curve.get(1).coverage, 1e-9);
	}

	@Test
	void oldStepGuardOfOnePercentBecomesTheDefault() {
		JSONObject additional = new JSONObject();
		additional.put(GenomeSelectionCut.KEY_STEP_GUARD, 0.01);
		GenomeSelectionCut.LoadedCoverage loaded = GenomeSelectionCut.readAdditionalKeys(additional, 0.95, 0.001, 300);
		assertEquals(0.001, loaded.stepGuard, 1e-12);
		assertTrue(loaded.replacementWarning.contains("0.01"));
		assertTrue(loaded.replacementWarning.contains("0.001"));
	}
}
