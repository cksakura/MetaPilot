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
 * the separate identification-only walks, and the engine exclude-range exemption of
 * {@link DeltaMassEvidenceGate}.
 */
class ModificationValidatorScreenTest {

	private static final String LONG_PEPTIDE = "ACDEFGHIKLMNPQR";
	private static final String SHORT_PEPTIDE = "ACDEFGHK";

	private static void addPeak(List<OpenSearchPsm> psms, int peakId, double centroid, int targets, int decoys) {
		addPeak(psms, peakId, centroid, targets, decoys, LONG_PEPTIDE);
	}

	private static void addPeak(List<OpenSearchPsm> psms, int peakId, double centroid, int targets, int decoys,
			String peptide) {
		for (int i = 0; i < targets + decoys; i++) {
			OpenSearchPsm psm = new OpenSearchPsm();
			psm.setPeakId(peakId);
			psm.setPeakCentroid(centroid);
			psm.setHitRank(1);
			psm.setTarget(i < targets);
			psm.setClassLogit(100 - i);
			psm.setPeptide(peptide);
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
	@DisplayName("Any mass with the evidence is established, whether or not it has a name")
	void evidenceDecidesNotTheMass() {
		List<OpenSearchPsm> psms = new ArrayList<>();
		addPeak(psms, 1, 0.0, 1000, 20);
		addPeak(psms, 2, 42.0106, 200, 4);
		addPeak(psms, 3, 123.4567, 300, 5);
		addPeak(psms, 4, 7.5, 160, 40);

		Map<Integer, ModificationValidator.PeakVerdict> verdicts = byId(
				ModificationValidator.filterPeaks(psms, new ModificationValidator.Config(), null, null));

		assertEquals(ModificationValidator.MassClass.ESTABLISHED, verdicts.get(2).massClass);
		assertTrue(verdicts.get(2).accepted);
		assertTrue(verdicts.get(3).accepted, "a mass no list ever named");
		assertEquals(204, verdicts.get(2).longPeptidePsms);
		assertEquals(4.0 / 204, verdicts.get(2).longPeptideDecoyFraction, 1e-12);
		assertEquals(ModificationValidator.MassClass.UNEXPLAINED, verdicts.get(4).massClass);
		assertFalse(verdicts.get(4).accepted, "20% long-peptide decoys");
	}

	@Test
	@DisplayName("A peak of short peptides fails however large it is: the residue-mass case")
	void shortPeptidePeakFails() {
		List<OpenSearchPsm> psms = new ArrayList<>();
		addPeak(psms, 1, 0.0, 1000, 20);
		addPeak(psms, 2, 15.9949, 300, 3);
		// -Q as on faecal: nearly all short, and the few long ones clean.
		addPeak(psms, 3, -128.0586, 1500, 1400, SHORT_PEPTIDE);
		addPeak(psms, 3, -128.0586, 60, 1);

		ModificationValidator.PeakVerdict verdict = byId(
				ModificationValidator.filterPeaks(psms, new ModificationValidator.Config(), null, null)).get(3);

		assertEquals(61, verdict.longPeptidePsms);
		assertTrue(verdict.longPeptideDecoyFraction <= 0.05, "clean among long peptides");
		assertTrue(verdict.longPeptideShareRatio < 0.75, "but skewed to short ones");
		assertEquals(ModificationValidator.MassClass.UNEXPLAINED, verdict.massClass);
		assertTrue(verdict.reason.startsWith("skewed to short peptides"), verdict.reason);
	}

	@Test
	@DisplayName("A peak with too few long peptides is not judged, however clean")
	void tooFewLongPeptides() {
		List<OpenSearchPsm> psms = new ArrayList<>();
		addPeak(psms, 1, 0.0, 1000, 20);
		addPeak(psms, 2, 2.1, 49, 0);
		addPeak(psms, 3, 7.5, 50, 0);

		Map<Integer, ModificationValidator.PeakVerdict> verdicts = byId(
				ModificationValidator.filterPeaks(psms, new ModificationValidator.Config(), Set.of(2), null));

		assertFalse(verdicts.get(2).accepted);
		assertTrue(verdicts.get(2).fragmentSupported, "fragment support is recorded for audit");
		assertTrue(verdicts.get(3).accepted, "the minimum is inclusive");
	}

	@Test
	@DisplayName("A substitution-shaped peak that passes is withdrawn, with the evidence kept in its verdict")
	void substitutionPeakIsWithdrawn() {
		List<OpenSearchPsm> psms = new ArrayList<>();
		addPeak(psms, 1, 0.0, 1000, 20);
		addPeak(psms, 2, 14.0157, 500, 5);

		ModificationValidator.Gate gate = ModificationValidator.filterPeaks(psms, new ModificationValidator.Config(),
				null, null);
		assertTrue(gate.acceptedPeakIds.contains(2));
		ModificationValidator.Gate withdrawn = gate.withoutPeaks(Set.of(2));

		assertFalse(withdrawn.acceptedPeakIds.contains(2));
		ModificationValidator.PeakVerdict verdict = byId(withdrawn).get(2);
		assertFalse(verdict.accepted);
		assertEquals(505, verdict.longPeptidePsms);
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
	@DisplayName("The evidence gate never folds a peak inside the engine's exclude range, and reports supported peaks")
	void evidenceGateExemptsOnlyTheExcludeRange() {
		List<OpenSearchPsm> psms = new ArrayList<>();
		addGainPeak(psms, 0, 0.0, 10, 0.0);
		addGainPeak(psms, 2, 0.9840, 60, 0.0);
		addGainPeak(psms, 3, 4.0, 60, 0.0);
		addGainPeak(psms, 4, 100.0, 60, 5.0);
		addGainPeak(psms, 5, -1.2, 60, 0.0);
		addGainPeak(psms, 6, 15.9949, 60, 0.0);

		DeltaMassEvidenceGate.Result result = DeltaMassEvidenceGate.apply(psms, null,
				new DeltaMassEvidenceGate.Config(), null);

		// 0.984 and -1.2 sit inside (-1.5, 3.5), where MSFragger scores no shifted fragment, so no
		// gain is no evidence. A named mass outside it, here oxidation, gets no exemption.
		assertEquals(Set.of(3, 6), result.rejectedPeakIds);
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
