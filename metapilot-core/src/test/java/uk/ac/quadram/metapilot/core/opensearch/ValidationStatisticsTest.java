package uk.ac.quadram.metapilot.core.opensearch;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

/**
 * Tests for the statistical machinery behind PSM validation: monotonised
 * q-values, the peak screen's demotion of rejected peaks, and spectrum-level
 * fold assignment.
 * <p>
 * The binomial peak gate and its Benjamini-Hochberg correction were removed on
 * 2026-09-13, with the tests that exercised them. A p-value depends on how many
 * PSMs a peak holds, so the same peak quality drew opposite verdicts at
 * different sizes - on 12mix it accepted peak 1 at 23.5% decoys (n=1,355) and
 * rejected peak 36 at 23.8% (n=21) - and the problem grows with the experiment.
 */
class ValidationStatisticsTest {

	// ---------------------------------------------------------------- q-values

	@Test
	@DisplayName("q-values are monotone non-increasing with score")
	void qValuesAreMonotone() {
		Random rng = new Random(7);
		int n = 2000;
		double[] scores = new double[n];
		boolean[] isTarget = new boolean[n];
		for (int i = 0; i < n; i++) {
			isTarget[i] = rng.nextBoolean();
			// Targets score better on average; decoys are the null.
			scores[i] = rng.nextGaussian() + (isTarget[i] ? 1.5 : 0.0);
		}

		double[] q = TargetDecoyFdr.qValues(scores, isTarget);

		Integer[] order = new Integer[n];
		for (int i = 0; i < n; i++) {
			order[i] = i;
		}
		java.util.Arrays.sort(order, (a, b) -> Double.compare(scores[b], scores[a]));

		double previous = -1;
		for (int rank = n - 1; rank >= 0; rank--) {
			double current = q[order[rank]];
			if (previous >= 0) {
				assertTrue(current <= previous + 1e-12,
						"q must not increase as score increases (rank " + rank + ")");
			}
			previous = current;
		}
	}

	@Test
	@DisplayName("An early decoy spike no longer truncates the list, unlike break-on-first-violation")
	void monotonisationSurvivesAnEarlyDecoy() {
		// Hand-built ranking with one unlucky decoy at rank 1. The raw running
		// ratio there is 1/1 = 1.0, so the original implementation - which walked
		// down and broke at the first rank exceeding the threshold - stopped dead
		// after a single PSM. Deeper ranks are far cleaner (6 targets, 1 decoy),
		// so the monotone q-value must reflect what is actually achievable.
		double[] scores = { 10, 9, 8, 7, 6, 5, 4, 3 };
		boolean[] isTarget = { true, false, true, true, true, true, true, false };

		double[] q = TargetDecoyFdr.qValues(scores, isTarget);

		// Raw d/t at rank 1 is 1.0; monotonised it becomes 1/6 = 0.1667.
		assertEquals(0.0, q[0], 1e-12, "the leading target should be q=0");
		assertEquals(1.0 / 6.0, q[1], 1e-9,
				"the early decoy's q must be pulled down to the best value achievable deeper in the list");

		// The concrete consequence at a 20% local FDR: monotone + deepest-passing
		// accepts seven PSMs, whereas break-on-first-violation accepted only one.
		int monotoneAccepted = 0;
		for (double value : q) {
			if (value <= 0.2) {
				monotoneAccepted++;
			}
		}
		assertEquals(7, monotoneAccepted, "monotone q at 20% local FDR should accept ranks 0..6");

		int breakOnFirstViolation = 0;
		int targets = 0;
		int decoys = 0;
		for (int rank = 0; rank < scores.length; rank++) {
			if (isTarget[rank]) {
				targets++;
			} else {
				decoys++;
			}
			if (targets > 0 && (double) decoys / targets <= 0.2) {
				breakOnFirstViolation++;
			} else {
				break;
			}
		}
		assertEquals(1, breakOnFirstViolation, "sanity: the old behaviour truncates at the first breach");
	}

	@Test
	@DisplayName("A 20%-decoy known-shift peak is rejected, and its PSMs are demoted rather than deleted")
	void decoyRichPeakIsRejected() {
		List<OpenSearchPsm> psms = new ArrayList<>();
		for (int i = 0; i < 4; i++) {
			OpenSearchPsm target = new OpenSearchPsm();
			target.setPeakId(3);
			target.setPeakCentroid(15.9949);
			target.setTarget(true);
			target.setClassLogit(10 - i);
			psms.add(target);
		}
		OpenSearchPsm decoy = new OpenSearchPsm();
		decoy.setPeakId(3);
		decoy.setPeakCentroid(15.9949);
		decoy.setTarget(false);
		decoy.setClassLogit(5);
		psms.add(decoy);

		ModificationValidator.Result result = ModificationValidator.validate(psms,
				new ModificationValidator.Config());

		assertEquals(1, result.verdicts.size());
		assertFalse(result.verdicts.get(0).accepted,
				"1 decoy in 5 has a Wilson upper bound of 0.62, above the 0.5 null background");
		assertTrue(result.peakValidated.isEmpty(), "a rejected peak makes no modification claim");

		// Rejecting the peak rejects the mass shift, not the spectra. The PSMs fall
		// through to the rejected-peak identification walk, where the four clean targets sit
		// above the single decoy and clear the walk. Deleting them here would let
		// peak membership decide identification - the failure this test pins down.
		assertEquals(4, result.unexplainedIdentifications.size(),
				"rejected-peak PSMs are demoted to the identification-only pool, not discarded");
		for (OpenSearchPsm psm : result.unexplainedIdentifications) {
			assertEquals(3, psm.getPeakId(), "demoted PSMs keep their peak assignment for audit");
			assertTrue(psm.isTarget());
		}
		assertEquals(4, result.validated.size(), "the reporting union carries the demoted identifications");
	}

	// ------------------------------------------------------- fold assignment

	@Test
	@DisplayName("Folds are keyed on the spectrum, so competing candidates never straddle a split")
	void foldsAreAssignedPerSpectrum() {
		List<OpenSearchPsm> psms = new ArrayList<>();
		for (int scan = 1; scan <= 300; scan++) {
			for (int rank = 1; rank <= 4; rank++) {
				OpenSearchPsm psm = new OpenSearchPsm();
				psm.setFileName("runA");
				psm.setScan(scan);
				psm.setHitRank(rank);
				psms.add(psm);
			}
		}

		int[] folds = SemiSupervisedPsmClassifier.assignFolds(psms, 3);

		Map<String, Integer> seen = new HashMap<>();
		for (int i = 0; i < psms.size(); i++) {
			String key = psms.get(i).getFileName() + ":" + psms.get(i).getScan();
			Integer existing = seen.get(key);
			if (existing == null) {
				seen.put(key, folds[i]);
			} else {
				// This is the leak the split exists to prevent: rank-1 in train
				// and rank-2 in test would be near-duplicate feature vectors.
				assertEquals(existing.intValue(), folds[i],
						"all candidates of spectrum " + key + " must share a fold");
			}
		}

		Set<Integer> used = new HashSet<>();
		for (int f : folds) {
			used.add(f);
		}
		assertEquals(3, used.size(), "all folds should be populated");
	}

	@Test
	@DisplayName("Fold assignment is reproducible across runs")
	void foldAssignmentIsDeterministic() {
		List<OpenSearchPsm> psms = new ArrayList<>();
		for (int scan = 1; scan <= 100; scan++) {
			OpenSearchPsm psm = new OpenSearchPsm();
			psm.setFileName("runB");
			psm.setScan(scan);
			psms.add(psm);
		}
		int[] first = SemiSupervisedPsmClassifier.assignFolds(psms, 4);
		int[] second = SemiSupervisedPsmClassifier.assignFolds(psms, 4);
		org.junit.jupiter.api.Assertions.assertArrayEquals(first, second,
				"fold assignment must not depend on run order or RNG state");
	}

	// ------------------------------------------------------------------ logit

	@Test
	@DisplayName("Logit is finite at the probability boundaries")
	void logitIsFiniteAtBoundaries() {
		assertTrue(Double.isFinite(SemiSupervisedPsmClassifier.toLogit(0.0)));
		assertTrue(Double.isFinite(SemiSupervisedPsmClassifier.toLogit(1.0)));
		assertEquals(0.0, SemiSupervisedPsmClassifier.toLogit(0.5), 1e-12);
		assertTrue(SemiSupervisedPsmClassifier.toLogit(0.9) > SemiSupervisedPsmClassifier.toLogit(0.1));
	}

	@Test
	@DisplayName("Unassigned mass is identification-only, never a pooled modification hypothesis")
	void unassignedMassIsKeptInASeparatePopulation() {
		List<OpenSearchPsm> psms = new ArrayList<>();
		for (int i = 0; i < 10; i++) {
			OpenSearchPsm psm = new OpenSearchPsm();
			psm.setPeakId(7);
			psm.setPeakCentroid(0.0);
			psm.setTarget(true);
			psm.setClassLogit(10 - i);
			psm.setQValue(0.005);
			psms.add(psm);
		}
		OpenSearchPsm retainedUnexplained = new OpenSearchPsm();
		retainedUnexplained.setPeakId(-1);
		retainedUnexplained.setTarget(true);
		retainedUnexplained.setClassLogit(3.0);
		psms.add(retainedUnexplained);

		OpenSearchPsm unexplainedDecoy = new OpenSearchPsm();
		unexplainedDecoy.setPeakId(-1);
		unexplainedDecoy.setTarget(false);
		unexplainedDecoy.setClassLogit(2.0);
		psms.add(unexplainedDecoy);

		OpenSearchPsm rejectedUnexplained = new OpenSearchPsm();
		rejectedUnexplained.setPeakId(-1);
		rejectedUnexplained.setTarget(true);
		rejectedUnexplained.setClassLogit(1.0);
		psms.add(rejectedUnexplained);

		ModificationValidator.Config cfg = new ModificationValidator.Config();
		cfg.unexplainedFdr = 0.4;
		ModificationValidator.Result result = ModificationValidator.validate(psms, cfg);

		assertEquals(1, result.verdicts.size(), "only fitted peaks receive modification verdicts");
		assertEquals(7, result.verdicts.get(0).peakId);
		assertEquals(10, result.peakValidated.size());
		assertEquals(1, result.unexplainedIdentifications.size());
		assertTrue(result.unexplainedIdentifications.contains(retainedUnexplained));
		assertFalse(result.unexplainedIdentifications.contains(rejectedUnexplained));
		assertEquals(11, result.validated.size(), "reporting union preserves confident peptide coverage");
	}
}
