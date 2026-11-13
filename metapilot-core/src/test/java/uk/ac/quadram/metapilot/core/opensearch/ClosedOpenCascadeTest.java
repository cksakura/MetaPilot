package uk.ac.quadram.metapilot.core.opensearch;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.EnumMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.dom4j.DocumentHelper;
import org.dom4j.Element;
import org.junit.jupiter.api.Test;

import uk.ac.quadram.metapilot.core.opensearch.ClosedOpenCascade.IdClass;
import uk.ac.quadram.metapilot.core.opensearch.ModificationValidator.MassClass;
import uk.ac.quadram.metapilot.core.opensearch.ModificationValidator.PeakVerdict;

class ClosedOpenCascadeTest {

	private static int nextScan = 1;

	private static OpenSearchPsm psm(String peptide, boolean target, double score, int peakId) {
		OpenSearchPsm psm = new OpenSearchPsm();
		psm.setFileName("run1");
		psm.setScan(nextScan++);
		psm.setPeptide(peptide);
		psm.setTarget(target);
		psm.setClassScore(score);
		psm.setPeakId(peakId);
		return psm;
	}

	private static PeakVerdict verdict(int peakId, double centroid, MassClass massClass, boolean accepted) {
		return new PeakVerdict(peakId, centroid, massClass, 0, 0, 0.0, 0.0, Double.NaN, true, accepted, "test",
				0);
	}

	/**
	 * Synthetic closed search, two candidates per spectrum, ranked by hyperscore as the engine ranks
	 * them. The predicted-spectrum entropy tells the right peptide; the hyperscore does not always: in
	 * a fifth of the spectra a near-identical wrong target outscores the right peptide.
	 */
	@Test
	void closedRerankingPromotesTheRank2CandidateThePredictedSpectrumSupports() throws Exception {
		java.util.Random random = new java.util.Random(7);
		List<OpenSearchPsm> closed = new ArrayList<>();
		Set<OpenSearchPsm> right = java.util.Collections.newSetFromMap(new java.util.IdentityHashMap<>());
		int spectra = 3000;
		for (int s = 0; s < spectra; s++) {
			double u = random.nextDouble();
			boolean hasTrue = u < 0.8;
			boolean swap = u >= 0.6 && u < 0.8;
			// {hyperscore, right, target, entropy}
			double[][] c = new double[2][];
			double rightScore = 30 + 4 * random.nextGaussian();
			for (int k = 0; k < 2; k++) {
				boolean isRight = hasTrue && k == 0;
				double hyperscore = isRight ? rightScore
						// A wrong candidate that outscores the right one is a near-identical target peptide.
						: swap ? rightScore + 0.5 + Math.abs(random.nextGaussian()) : 21 + 4 * random.nextGaussian();
				boolean target = isRight || swap || random.nextBoolean();
				double entropy = isRight ? 0.6 + 0.4 * random.nextDouble() : 0.45 * random.nextDouble();
				c[k] = new double[] { hyperscore, isRight ? 1 : 0, target ? 1 : 0, entropy };
			}
			java.util.Arrays.sort(c, (a, b) -> Double.compare(b[0], a[0]));
			double nextscore = 17 + 2 * random.nextGaussian();
			for (int rank = 1; rank <= 2; rank++) {
				double[] x = c[rank - 1];
				OpenSearchPsm psm = closedPsm(s, rank, x[2] == 1, x[0], nextscore, x[3], x[3] - c[2 - rank][3], random);
				closed.add(psm);
				if (x[1] == 1) {
					right.add(psm);
				}
			}
		}
		int rightAtRankOne = 0;
		for (OpenSearchPsm psm : right) {
			rightAtRankOne += psm.getHitRank() == 1 ? 1 : 0;
		}
		SemiSupervisedPsmClassifier.Config cfg = new SemiSupervisedPsmClassifier.Config();
		ClosedOpenCascade.ClosedStage plain = ClosedOpenCascade.validateClosed(closed, cfg, 0.01, null);
		for (OpenSearchPsm psm : plain.winners) {
			assertEquals(1, psm.getHitRank());
		}
		ClosedOpenCascade.ClosedStage reranked = ClosedOpenCascade.validateClosed(closed, cfg, 0.01, null, true,
				PredictedSpectrumFeatures.PSM_FEATURES);
		assertEquals(spectra, reranked.winners.size());
		int rightChosen = 0;
		int rightRankOneKept = 0;
		for (OpenSearchPsm psm : reranked.winners) {
			if (right.contains(psm)) {
				rightChosen++;
				rightRankOneKept += psm.getHitRank() == 1 ? 1 : 0;
			}
		}
		assertTrue(right.size() - rightAtRankOne > 0.15 * right.size(), "the test needs rank-2 right answers");
		assertTrue(rightChosen > 0.97 * right.size(), rightChosen + " of " + right.size());
		assertTrue(rightRankOneKept > 0.99 * rightAtRankOne, rightRankOneKept + " of " + rightAtRankOne);
		assertTrue(reranked.acceptedTargets > plain.acceptedTargets, reranked.acceptedTargets + " vs " + plain.acceptedTargets);
	}

	private static OpenSearchPsm closedPsm(int scan, int rank, boolean target, double hyperscore, double nextscore,
			double entropy, double entropyMargin, java.util.Random random) {
		StringBuilder peptide = new StringBuilder();
		String residues = "ACDEFGHKLMNPQRSTVWY";
		for (int k = 0; k < 9 + random.nextInt(6); k++) {
			peptide.append(residues.charAt(random.nextInt(residues.length())));
		}
		peptide.append('K');
		OpenSearchPsm psm = new OpenSearchPsm();
		psm.setFileName("run1");
		psm.setScan(scan);
		psm.setHitRank(rank);
		psm.setClosedStage(true);
		psm.setPeptide(peptide.toString());
		psm.setProtein((target ? "P" : "rev_P") + scan + "_" + rank);
		psm.setTarget(target);
		psm.setCharge(2);
		psm.setPrecursorNeutralMass(1200 + 400 * random.nextDouble());
		psm.setMassDiff(0.001 * random.nextGaussian());
		psm.setNumMatchedIons((int) Math.max(1, Math.round(hyperscore / 3 + random.nextGaussian())));
		psm.setTotNumIons(18 + random.nextInt(8));
		psm.setHyperscore(hyperscore);
		psm.setNextscore(nextscore);
		psm.setExpect(Math.pow(10, -(hyperscore - nextscore) / 3 + 0.3 * random.nextGaussian()));
		Map<String, Double> external = new java.util.HashMap<>();
		external.put(PredictedSpectrumFeatures.ENTROPY, entropy);
		external.put(PredictedSpectrumFeatures.MATCHED_FRACTION, Math.min(1, entropy + 0.2 * random.nextDouble()));
		external.put(PredictedSpectrumFeatures.ENTROPY_MARGIN, entropyMargin);
		external.put(PredictedSpectrumFeatures.RT_ERROR, Math.abs(random.nextGaussian()) * (entropy > 0.5 ? 1 : 3));
		psm.setExternalFeatures(external);
		return psm;
	}

	@Test
	void peptideFilterTightensEveryClassWhenTheUnionExceedsTheFdr() {
		Map<IdClass, List<OpenSearchPsm>> pools = new EnumMap<>(IdClass.class);
		for (IdClass c : IdClass.values()) {
			List<OpenSearchPsm> pool = new ArrayList<>();
			for (int i = 0; i < 100; i++) {
				pool.add(psm("PEPTK" + i, true, 200 - i, 0));
			}
			pool.add(psm("DECOY" + c.name(), false, 50, 0));
			pools.put(c, pool);
		}

		ClosedOpenCascade.PeptideFilter filter = ClosedOpenCascade.filterPeptides(pools, pools, 0.01);

		for (IdClass c : IdClass.values()) {
			assertArrayEquals(new int[] { 100, 1 }, filter.peptidesBefore.get(c), c.label);
			assertArrayEquals(new int[] { 100, 0 }, filter.peptidesAfter.get(c), c.label);
			assertArrayEquals(new int[] { 100, 0 }, filter.psmsKept.get(c), c.label);
		}
		assertEquals(100, filter.unionTargetsBefore);
		assertEquals(3, filter.unionDecoysBefore);
		assertEquals(0, filter.unionDecoysAfter);
		assertEquals(0.9, filter.factor, 1e-12);
		assertEquals(300, filter.kept.size());
		assertTrue(filter.kept.stream().allMatch(OpenSearchPsm::isTarget));
	}

	@Test
	void peptideFilterLeavesTheThresholdAloneWhenTheUnionPasses() {
		List<OpenSearchPsm> pool = new ArrayList<>();
		for (int i = 0; i < 100; i++) {
			pool.add(psm("PEPTK" + i, true, 200 - i, 0));
		}
		pool.add(psm("DECOYK", false, 50, 0));
		Map<IdClass, List<OpenSearchPsm>> pools = new EnumMap<>(IdClass.class);
		pools.put(IdClass.CLOSED, pool);

		ClosedOpenCascade.PeptideFilter filter = ClosedOpenCascade.filterPeptides(pools, pools, 0.01);

		assertEquals(1.0, filter.factor, 0.0);
		assertEquals(101, filter.kept.size());
		assertArrayEquals(new int[] { 0, 0 }, filter.psmsKept.get(IdClass.OPEN_ACCEPTED_SHIFT));
	}

	@Test
	void peptideUnitsMergeIsoleucineAndLeucine() {
		List<OpenSearchPsm> pool = List.of(psm("PEPTIK", true, 10, 0), psm("PEPTLK", true, 5, 0));
		Map<String, Double> q = ClosedOpenCascade.peptideQValues(pool);
		assertEquals(1, q.size());
	}

	@Test
	void openStageIdentifiesOnlyAcceptedZeroAndKnownPeaks() {
		List<OpenSearchPsm> remaining = new ArrayList<>();
		for (int i = 0; i < 20; i++) {
			remaining.add(psm("ZERO" + i, true, 100 - i, 0));
		}
		for (int i = 0; i < 10; i++) {
			remaining.add(psm("KNOWN" + i, true, 100 - i, 1));
		}
		for (int i = 0; i < 10; i++) {
			remaining.add(psm("SHIFT" + i, true, 100 - i, 2));
		}
		for (int i = 0; i < 5; i++) {
			remaining.add(psm("NOPEAK" + i, true, 100 - i, -1));
		}
		ModificationValidator.Gate gate = new ModificationValidator.Gate(Set.of(0, 1),
				List.of(verdict(0, 0.0, MassClass.ZERO, true), verdict(1, 15.9949, MassClass.ESTABLISHED, true),
						verdict(2, 7.3, MassClass.UNEXPLAINED, false)),
				Set.of());

		ModificationValidator.Result result = ModificationValidator.validateRemaining(remaining, null, gate);

		assertEquals(30, result.validated.size());
		assertTrue(result.validated.stream().allMatch(p -> p.getPeakId() == 0 || p.getPeakId() == 1));
		assertEquals(15, result.unexplainedIdentifications.size());
		assertTrue(result.unexplainedIdentifications.stream().noneMatch(result.validated::contains));
		assertEquals(45, result.diagnostic.size());
		PeakVerdict shifted = result.verdicts.stream().filter(v -> v.peakId == 2).findFirst().orElseThrow();
		assertEquals(0, shifted.psmsRetained);
		assertFalse(shifted.accepted);
	}

	@Test
	void gateRejectsWithoutFoldingWhenAsked() {
		List<OpenSearchPsm> psms = new ArrayList<>();
		for (int i = 0; i < 10; i++) {
			OpenSearchPsm p = psm("ZERO" + i, true, 10, 0);
			p.setPeakCentroid(0.0);
			psms.add(p);
		}
		List<OpenSearchPsm> shifted = new ArrayList<>();
		for (int i = 0; i < 60; i++) {
			OpenSearchPsm p = psm("SHIFT" + i, true, 10, 1);
			p.setPeakCentroid(4.718);
			p.setScoreWithDeltaMass(5.0);
			p.setScoreWithoutDeltaMass(5.0);
			shifted.add(p);
		}
		psms.addAll(shifted);
		DeltaMassEvidenceGate.Config cfg = new DeltaMassEvidenceGate.Config();
		cfg.foldRejected = false;

		DeltaMassEvidenceGate.Result result = DeltaMassEvidenceGate.apply(psms, null, cfg, null);

		assertTrue(result.rejectedPeakIds.contains(1));
		assertEquals(0, result.psmsReassigned);
		for (OpenSearchPsm p : shifted) {
			assertEquals(1, p.getPeakId());
			assertEquals(4.718, p.getPeakCentroid(), 0.0);
		}
	}

	@Test
	void closedPsmsCarryTheirVariableModificationsAndIsotopeOffset() {
		OpenSearchPsm p = psm("PEPMK", true, 10, -1);
		p.setClosedStage(true);
		p.setClosedModifications(new double[] { 42.0106, 0, 0, 0, 15.9949, 0 });
		p.setMassDiff(1.0040);
		p.setCalcNeutralPepMass(700.0);

		ClosedOpenCascade.prepareClosed(List.of(p));

		assertEquals(1, p.getPeakIsotopeOffset());
		assertEquals(1.0040 - ClosedOpenCascade.ISOTOPE_SPACING, p.getNormalizedMassDiff(), 1e-9);
		assertEquals(2.0, p.getExternalFeatures().get(ClosedOpenCascade.VARIABLE_MODIFICATIONS_FEATURE), 0.0);
		assertEquals(2, p.getClosedModificationCount());
		assertEquals(42.0106 + 15.9949, p.getClosedModificationMass(), 1e-9);
		assertEquals("[+42.0106]-PEPM[+15.9949]K", p.getClosedModifiedSequence());
		assertEquals(p.getClosedModifiedSequence(), p.getModificationFamilyKey());
		assertEquals(-1, p.getPeakId());
	}

	@Test
	void readerTakesVariableModificationsAndIgnoresFixedCysteine() {
		Element hit = DocumentHelper.createElement("search_hit");
		Element info = hit.addElement("modification_info");
		info.addAttribute("mod_nterm_mass", "43.0184");
		info.addElement("mod_aminoacid_mass").addAttribute("position", "2").addAttribute("mass", "160.0307");
		info.addElement("mod_aminoacid_mass").addAttribute("position", "4").addAttribute("mass", "147.0354")
				.addAttribute("variable", "15.9949");
		info.addElement("mod_aminoacid_mass").addAttribute("position", "5").addAttribute("mass", "147.0354");

		double[] mods = OpenSearchPepxmlReader.readVariableModifications(hit, "ACDMMK");

		assertNotNull(mods);
		assertEquals(7, mods.length);
		assertEquals(42.0106, mods[0], 1e-3);
		assertEquals(0.0, mods[2], 0.0);
		assertEquals(15.9949, mods[4], 1e-4);
		assertEquals(15.9949, mods[5], 1e-3);
	}

	@Test
	void readerReturnsNullForAnUnmodifiedPeptide() {
		Element hit = DocumentHelper.createElement("search_hit");
		hit.addElement("modification_info").addElement("mod_aminoacid_mass").addAttribute("position", "1")
				.addAttribute("mass", "160.0307");
		assertNull(OpenSearchPepxmlReader.readVariableModifications(hit, "CPEPK"));
	}
}
