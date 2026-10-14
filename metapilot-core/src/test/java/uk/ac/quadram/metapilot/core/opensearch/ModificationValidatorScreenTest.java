package uk.ac.quadram.metapilot.core.opensearch;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

/**
 * Tests for the per-peak screen in {@link ModificationValidator#filterPeaks},
 * the separate identification-only walks, and the allow-list exemption in
 * {@link DeltaMassEvidenceGate}.
 */
class ModificationValidatorScreenTest {

	private static void addPeak(List<OpenSearchPsm> psms, int peakId, double centroid, int targets, int decoys) {
		for (int i = 0; i < targets + decoys; i++) {
			OpenSearchPsm psm = new OpenSearchPsm();
			psm.setPeakId(peakId);
			psm.setPeakCentroid(centroid);
			psm.setHitRank(1);
			psm.setTarget(i < targets);
			psm.setClassLogit(100 - i);
			psms.add(psm);
		}
	}

	private static Map<Integer, ModificationValidator.PeakVerdict> byId(ModificationValidator.Gate gate) {
		Map<Integer, ModificationValidator.PeakVerdict> out = new HashMap<>();
		for (ModificationValidator.PeakVerdict v : gate.verdicts) {
			out.put(v.peakId, v);
		}
		return out;
	}

	@Test
	@DisplayName("The zero peak is accepted even when small clean peaks would fill a cumulative budget first")
	void zeroPeakIsAlwaysAccepted() {
		List<OpenSearchPsm> psms = new ArrayList<>();
		// 3% decoys: a cumulative 1.5% walk over 30 clean 20/0 peaks plus this one
		// reaches 30/1,600 = 1.9% and stops before admitting it - the PXD050781 failure.
		addPeak(psms, 1, -0.0004, 970, 30);
		for (int i = 0; i < 30; i++) {
			addPeak(psms, 10 + i, 100.0 + i, 20, 0);
		}

		ModificationValidator.Gate gate = ModificationValidator.filterPeaks(psms, new ModificationValidator.Config(),
				null, null);

		Map<Integer, ModificationValidator.PeakVerdict> verdicts = byId(gate);
		assertTrue(verdicts.get(1).accepted);
		assertEquals(ModificationValidator.MassClass.ZERO, verdicts.get(1).massClass);
		assertEquals(Set.of(1), gate.acceptedPeakIds, "unexplained masses make no modification claim");
	}

	@Test
	@DisplayName("A known shift is accepted when target-enriched against the other non-zero peaks, not otherwise")
	void knownShiftJudgedAgainstBackground() {
		List<OpenSearchPsm> psms = new ArrayList<>();
		addPeak(psms, 1, 0.0, 1000, 20);
		addPeak(psms, 2, 15.9949, 200, 4);
		addPeak(psms, 3, 7.5, 160, 40);
		addPeak(psms, 4, 21.9819, 30, 10);

		Map<Integer, ModificationValidator.PeakVerdict> verdicts = byId(
				ModificationValidator.filterPeaks(psms, new ModificationValidator.Config(), null, null));

		ModificationValidator.PeakVerdict oxidation = verdicts.get(2);
		assertEquals(ModificationValidator.MassClass.KNOWN, oxidation.massClass);
		assertEquals(50.0 / 240.0, oxidation.backgroundDecoyFraction, 1e-12,
				"background is every other non-zero peak, leave-one-out");
		assertTrue(oxidation.decoyUpperBound < oxidation.backgroundDecoyFraction);
		assertTrue(oxidation.accepted);

		ModificationValidator.PeakVerdict sodium = verdicts.get(4);
		assertEquals(ModificationValidator.MassClass.KNOWN, sodium.massClass);
		assertFalse(sodium.accepted, "25% decoys is not enriched against an 11% background");
	}

	@Test
	@DisplayName("A clean, fragment-supported unexplained peak is still not accepted")
	void unexplainedMassIsNeverAccepted() {
		List<OpenSearchPsm> psms = new ArrayList<>();
		addPeak(psms, 1, 0.0, 1000, 20);
		addPeak(psms, 2, 2.1, 500, 0);

		ModificationValidator.PeakVerdict verdict = byId(
				ModificationValidator.filterPeaks(psms, new ModificationValidator.Config(), Set.of(2), null)).get(2);

		assertEquals(ModificationValidator.MassClass.UNEXPLAINED, verdict.massClass);
		assertTrue(verdict.fragmentSupported, "fragment support is recorded for audit");
		assertFalse(verdict.accepted);
	}

	@Test
	@DisplayName("A small zero-decoy known peak is judged by its upper bound, not its raw fraction")
	void smallKnownPeakNeedsEnoughEvidence() {
		assertEquals(3.8415 / (15 + 3.8415), ModificationValidator.wilsonUpperBound(0, 15), 1e-4);

		List<OpenSearchPsm> psms = new ArrayList<>();
		addPeak(psms, 1, 0.0, 1000, 20);
		addPeak(psms, 2, 0.9840, 3, 0);
		addPeak(psms, 3, -17.0265, 30, 0);
		addPeak(psms, 4, 7.5, 160, 40);

		Map<Integer, ModificationValidator.PeakVerdict> verdicts = byId(
				ModificationValidator.filterPeaks(psms, new ModificationValidator.Config(), null, null));

		assertFalse(verdicts.get(2).accepted, "0/3 decoys bounds at 0.56, above the 17% background");
		assertTrue(verdicts.get(3).accepted, "0/30 decoys bounds at 0.11, below the 20% background");
	}

	@Test
	@DisplayName("Rejected-peak and no-peak identifications are walked separately")
	void identificationPopulationsAreSeparate() {
		List<OpenSearchPsm> psms = new ArrayList<>();
		List<OpenSearchPsm> rejectedPeak = new ArrayList<>();
		for (int i = 0; i < 10; i++) {
			OpenSearchPsm psm = new OpenSearchPsm();
			psm.setPeakId(5);
			psm.setPeakCentroid(2.16);
			psm.setTarget(true);
			psm.setClassLogit(20 - i);
			psms.add(psm);
			rejectedPeak.add(psm);
		}
		OpenSearchPsm topDecoy = new OpenSearchPsm();
		topDecoy.setPeakId(-1);
		topDecoy.setTarget(false);
		topDecoy.setClassLogit(30);
		psms.add(topDecoy);
		for (int i = 0; i < 5; i++) {
			OpenSearchPsm psm = new OpenSearchPsm();
			psm.setPeakId(-1);
			psm.setTarget(true);
			psm.setClassLogit(9 - i);
			psms.add(psm);
		}

		ModificationValidator.Result result = ModificationValidator.validate(psms, new ModificationValidator.Config());

		// Pooled, the decoy at the top would give every target q >= 1/15 and none
		// would pass at 1%. Separately, the rejected peak's walk has no decoy.
		assertEquals(10, result.unexplainedIdentifications.size());
		assertTrue(result.unexplainedIdentifications.containsAll(rejectedPeak));
		for (OpenSearchPsm psm : rejectedPeak) {
			assertEquals(0.0, psm.getQValue(), 1e-12);
		}
	}

	@Test
	@DisplayName("The evidence gate never folds an allow-list mass, and reports supported peaks")
	void evidenceGateExemptsKnownShifts() {
		List<OpenSearchPsm> psms = new ArrayList<>();
		addGainPeak(psms, 0, 0.0, 10, 0.0);
		addGainPeak(psms, 2, 0.9840, 60, 0.0);
		addGainPeak(psms, 3, 4.0, 60, 0.0);
		addGainPeak(psms, 4, 100.0, 60, 5.0);

		DeltaMassEvidenceGate.Result result = DeltaMassEvidenceGate.apply(psms, null,
				new DeltaMassEvidenceGate.Config(), null);

		assertEquals(Set.of(3), result.rejectedPeakIds, "deamidation sits in MSFragger's excluded range");
		assertEquals(Set.of(4), result.supportedPeakIds);
		int deamidationKept = 0;
		int foldedOntoZero = 0;
		for (OpenSearchPsm psm : psms) {
			if (psm.getObservedPeakId() == 2 && psm.getPeakId() == 2) {
				deamidationKept++;
			}
			if (psm.getObservedPeakId() == 3 && psm.getPeakId() == 0) {
				foldedOntoZero++;
			}
		}
		assertEquals(60, deamidationKept);
		assertEquals(60, foldedOntoZero);
	}

	private static void addGainPeak(List<OpenSearchPsm> psms, int peakId, double centroid, int n, double gain) {
		for (int i = 0; i < n; i++) {
			OpenSearchPsm psm = new OpenSearchPsm();
			psm.setPeakId(peakId);
			psm.setObservedPeakId(peakId);
			psm.setPeakCentroid(centroid);
			psm.setTarget(true);
			psm.setScoreWithoutDeltaMass(10.0);
			psm.setScoreWithDeltaMass(10.0 + gain);
			psms.add(psm);
		}
	}
}
