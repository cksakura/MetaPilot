package uk.ac.quadram.metapilot.core.opensearch;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.EnumMap;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;

import org.junit.jupiter.api.Test;

class JointIdentificationTest {

	private static final String AA = "ACDEFGHKLMNPQRSTVWY";

	@Test
	void sameSpectrumAndPeptideFromBothSearchesIsOneExactCandidate() {
		OpenSearchPsm open = psm("a", 1, "PEPTIDEK", true, 1, 30, false, 0.0, 0);
		OpenSearchPsm closed = psm("a", 1, "PEPTIDEK", true, 1, 31, true, 0.0, -1);
		OpenSearchPsm other = psm("a", 1, "LLLLLLLK", true, 2, 12, false, 57.02, -1);
		List<JointCandidate> merged = JointCandidate.merge(Arrays.asList(open, other), Collections.singletonList(closed),
				new int[1]);
		assertEquals(2, merged.size());
		JointCandidate both = merged.get(0);
		assertTrue(both.foundByClosed() && both.foundByOpen());
		assertTrue(both.representative() == closed, "the closed row is reported: its interpretation is exact");
		ModificationValidator.Gate gate = gate();
		Map<Integer, ModificationValidator.MassClass> classes = ModificationValidator.massClasses(gate);
		assertEquals(JointCandidate.Group.EXACT, JointIdentification.groupOf(both, gate, classes));
		assertEquals(JointCandidate.Group.OPEN_UNEXPLAINED, JointIdentification.groupOf(merged.get(1), gate, classes));
	}

	@Test
	void groupsAreWalkedSeparately() {
		// Group EXACT: clean; group OPEN_UNMODIFIED: half decoys. Pooled, the open targets would ride on
		// the exact group's decoy count; walked alone they cannot pass.
		List<JointCandidate> candidates = new ArrayList<>();
		List<Double> scores = new ArrayList<>();
		for (int i = 0; i < 200; i++) {
			candidates.add(candidate("e" + i, true, JointCandidate.Group.EXACT));
			scores.add(10.0 + i);
		}
		for (int i = 0; i < 20; i++) {
			candidates.add(candidate("o" + i, i % 2 == 0, JointCandidate.Group.OPEN_UNMODIFIED));
			scores.add(100.0 + i);
		}
		double[] s = scores.stream().mapToDouble(Double::doubleValue).toArray();
		boolean[] winner = new boolean[s.length];
		Arrays.fill(winner, true);
		double[] q = JointIdentification.groupQValues(candidates, s, winner);
		for (int i = 0; i < 200; i++) {
			// (0 decoys + 1) / 200 targets: the conservative estimate, never 0.
			assertEquals(1.0 / 200, q[i], 1e-12);
		}
		for (int i = 200; i < 220; i++) {
			assertTrue(q[i] > 0.5, "an open group that is half decoys must not pass on the exact group's decoys");
		}
	}

	@Test
	void conservativeQValuesAddOneDecoyAndWalkTiesAsOneBlock() {
		// 40 targets, no decoy: q = 1/40, not 0.
		double[] s = new double[40];
		boolean[] t = new boolean[40];
		for (int i = 0; i < 40; i++) {
			s[i] = 100 - i;
			t[i] = true;
		}
		for (double v : TargetDecoyFdr.qValuesConservative(s, t)) {
			assertEquals(1.0 / 40, v, 1e-12);
		}
		// A target and a decoy at one score get one q-value, whichever row comes first.
		double[] q1 = TargetDecoyFdr.qValuesConservative(new double[] { 5, 5, 4 }, new boolean[] { true, false, true });
		double[] q2 = TargetDecoyFdr.qValuesConservative(new double[] { 5, 5, 4 }, new boolean[] { false, true, true });
		assertEquals(q1[0], q1[1], 1e-12);
		assertEquals(q2[0], q2[1], 1e-12);
		assertEquals(q1[0], q2[0], 1e-12);
	}

	@Test
	void logitIsRecoveredFromBothClassProbabilities() {
		// What weka's Logistic returns at a linear predictor of 50: the target probability has
		// rounded to 1.0, the decoy probability has not.
		double z = 50;
		double[] distribution = { 1.0 / (1.0 + Math.exp(-z)), 1.0 / (1.0 + Math.exp(z)) };
		assertEquals(1.0, distribution[0], 0.0);
		assertEquals(z, SemiSupervisedPsmClassifier.logitOf(distribution), 1e-9);
		assertEquals(-z, SemiSupervisedPsmClassifier.logitOf(new double[] { distribution[1], distribution[0] }), 1e-9);
		assertEquals(SemiSupervisedPsmClassifier.LOGIT_LIMIT,
				SemiSupervisedPsmClassifier.logitOf(new double[] { 1.0, 0.0 }), 0.0);
		assertEquals(0.0, SemiSupervisedPsmClassifier.logitOf(new double[] { 0.5, 0.5 }), 1e-12);
		assertEquals(0.5, SemiSupervisedPsmClassifier.sigmoid(0.0), 1e-12);
		assertEquals(1.0, SemiSupervisedPsmClassifier.sigmoid(800), 0.0);
	}

	@Test
	void mergeCountsLabelConflicts() {
		OpenSearchPsm open = psm("a", 1, "PEPTIDEK", true, 1, 30, false, 0.0, 0);
		OpenSearchPsm closed = psm("a", 1, "PEPTIDEK", true, 1, 31, true, 0.0, -1);
		closed.setTarget(false);
		int[] dropped = new int[2];
		List<JointCandidate> merged = JointCandidate.merge(Collections.singletonList(open),
				Collections.singletonList(closed), dropped);
		assertEquals(1, merged.size());
		assertTrue(merged.get(0).isTarget(), "the first row's label is kept");
		assertEquals(1, dropped[1]);
	}

	@Test
	void peptideBelongsToItsHighestPriorityGroup() {
		List<JointCandidate> candidates = new ArrayList<>();
		candidates.add(candidate("PEPTLDEK", true, JointCandidate.Group.OPEN_UNMODIFIED));
		candidates.add(candidate("PEPTIDEK", true, JointCandidate.Group.EXACT));
		candidates.add(candidate("NEWPEPK", true, JointCandidate.Group.OPEN_ACCEPTED_SHIFT));
		double[] scores = { 5, 1, 3 };
		boolean[] accepted = { true, true, true };
		JointIdentification.PeptideWalk walk = JointIdentification.peptideWalk(candidates, scores, accepted, 0.01);
		assertEquals(JointCandidate.Group.EXACT, walk.unitGroup.get("T:PEPTLDEK"), "I and L merge; exact wins");
		assertEquals(JointCandidate.Group.OPEN_ACCEPTED_SHIFT, walk.unitGroup.get("T:NEWPEPK"));
	}

	@Test
	void winnerIsTheBestCandidateOfEachSpectrum() {
		boolean[] w = JointIdentification.winners(new int[] { 0, 0, 1, 1, 1 }, new double[] { 1, 2, 5, Double.NaN, 4 });
		assertFalse(w[0]);
		assertTrue(w[1]);
		assertTrue(w[2]);
		assertFalse(w[3]);
		assertFalse(w[4]);
	}

	@Test
	void fragmentMassesAndSimilarity() {
		double[] none = new double[8];
		// PEPTIDE: b2 = P + E + proton, y1 = E + water + proton.
		assertEquals(227.1026, PredictedSpectrumFeatures.fragmentMz("PEPTIDE", none, 2, false, 1), 1e-3);
		assertEquals(148.0604, PredictedSpectrumFeatures.fragmentMz("PEPTIDE", none, 1, true, 1), 1e-3);
		double[] ox = new double[8];
		ox[3] = 15.9949; // residue 3 (P) carries the shift: b2 unchanged, b3 shifted
		assertEquals(227.1026, PredictedSpectrumFeatures.fragmentMz("PEPTIDE", ox, 2, false, 1), 1e-3);
		assertEquals(PredictedSpectrumFeatures.fragmentMz("PEPTIDE", none, 3, false, 1) + 15.9949,
				PredictedSpectrumFeatures.fragmentMz("PEPTIDE", ox, 3, false, 1), 1e-6);
		double[] p = { 1, 2, 3 };
		assertEquals(1.0, PredictedSpectrumFeatures.entropySimilarity(p, p.clone(), 3), 1e-9);
		assertEquals(0.0, PredictedSpectrumFeatures.entropySimilarity(new double[] { 1, 0 }, new double[] { 0, 1 }, 2),
				1e-9);
		double[] y = { 1, 3, 2, 4 };
		PredictedSpectrumFeatures.isotonic(y);
		assertTrue(y[0] <= y[1] && y[1] <= y[2] && y[2] <= y[3]);
		assertEquals(2.5, PredictedSpectrumFeatures.interpolate(new double[] { 0, 1 }, new double[] { 0, 1 }, 2.5), 1e-9);
	}

	@Test
	void predictedSpectrumScorePlacesTheShiftAndPrefersTheRightSite() {
		String peptide = "PEPTMDEK";
		// Predicted b2..b7 and y1..y7 at charge 1, equal intensities.
		List<Integer> ions = new ArrayList<>();
		List<Boolean> ys = new ArrayList<>();
		for (int ion = 1; ion < peptide.length(); ion++) {
			ions.add(ion);
			ys.add(false);
			ions.add(ion);
			ys.add(true);
		}
		int m = ions.size();
		int[] ion = new int[m];
		boolean[] y = new boolean[m];
		int[] z = new int[m];
		float[] inten = new float[m];
		for (int k = 0; k < m; k++) {
			ion[k] = ions.get(k);
			y[k] = ys.get(k);
			z[k] = 1;
			inten[k] = 100f;
		}
		uk.ac.quadram.metapilot.core.variant.DiannFragmentPredictor.Prediction p =
				new uk.ac.quadram.metapilot.core.variant.DiannFragmentPredictor.Prediction(ion, y, z, inten, 30f);
		// Observed: every fragment with +15.9949 on M (residue 5).
		double[] mods = new double[peptide.length() + 1];
		mods[5] = 15.9949;
		double[] mz = new double[m];
		for (int k = 0; k < m; k++) {
			mz[k] = PredictedSpectrumFeatures.fragmentMz(peptide, mods, ion[k], y[k], 1);
		}
		Arrays.sort(mz);
		float[] obs = new float[m];
		Arrays.fill(obs, 50f);
		uk.ac.quadram.metapilot.core.variant.MzmlScanReader.Peaks peaks =
				new uk.ac.quadram.metapilot.core.variant.MzmlScanReader.Peaks(mz, obs);
		// Unlocalized open candidate: the best placement is found.
		OpenSearchPsm open = psm("f", 1, peptide, true, 1, 20, false, 15.9949, 1);
		JointCandidate c = JointCandidate.merge(Collections.singletonList(open), null, new int[1]).get(0);
		double[] s = PredictedSpectrumFeatures.score(c, p, peaks, 20.0);
		assertEquals(1.0, s[0], 1e-9);
		assertEquals(1.0, s[1], 1e-9);
		// Localized on the wrong residue: worse.
		open.setLocalizedSequence("pEPTMDEK");
		double[] wrong = PredictedSpectrumFeatures.score(c, p, peaks, 20.0);
		assertTrue(wrong[0] < s[0]);
		// Closed candidate with the oxidation as a variable modification: exact.
		OpenSearchPsm closed = psm("f", 2, peptide, true, 1, 20, true, 15.9949, -1);
		closed.setClosedModifications(mods);
		JointCandidate cc = JointCandidate.merge(null, Collections.singletonList(closed), new int[1]).get(0);
		assertEquals(1.0, PredictedSpectrumFeatures.score(cc, p, peaks, 20.0)[0], 1e-9);
	}

	/**
	 * Synthetic run with known truth: the exact group's and the open groups' measured false discovery
	 * proportions must stay near their decoy estimates. Exact-mass peptides are found by the closed
	 * search in 70% of identifiable spectra; the rest only by the open one.
	 */
	@Test
	void groupFdrHoldsOnASyntheticRun() throws Exception {
		Random r = new Random(11);
		String[] present = peptides(r, 800, 'K');
		String[] background = peptides(r, 50_000, 'R');
		List<OpenSearchPsm> open = new ArrayList<>();
		List<OpenSearchPsm> closed = new ArrayList<>();
		Map<String, String> truth = new HashMap<>();
		for (int s = 0; s < 6000; s++) {
			String file = "run" + (s % 3);
			double u = r.nextDouble();
			String t = present[r.nextInt(present.length)];
			int type = u < 0.55 ? 0 : u < 0.62 ? 1 : 3;
			if (type <= 1) {
				truth.put(file + ":" + s, t);
			}
			if (type == 0 && r.nextDouble() < 0.7) {
				closed.add(psm(file, s, t, true, 1, 28 + 6 * r.nextGaussian(), true, 0, -1));
			} else if (r.nextDouble() < 0.8) {
				closed.add(psm(file, s, background[r.nextInt(background.length)], r.nextBoolean(), 1,
						13 + 4 * r.nextGaussian(), true, 0, -1));
			}
			int rank = 1;
			int first = open.size();
			if (type == 0) {
				open.add(psm(file, s, t, true, rank++, 25 + 6 * r.nextGaussian(), false, 0.0, 0));
			} else if (type == 1) {
				open.add(psm(file, s, t, true, rank++, 24 + 6 * r.nextGaussian(), false, 15.9949, 1));
			}
			for (int k = 0; k < 3; k++) {
				double dm = r.nextDouble() < 0.3 ? 0.0 : (r.nextDouble() - 0.5) * 400;
				open.add(psm(file, s, background[r.nextInt(background.length)], r.nextBoolean(), rank++,
						13 + 4 * r.nextGaussian(), false, dm, dm == 0.0 ? 0 : -1));
			}
		}
		ModificationValidator.Gate gate = gate();
		for (OpenSearchPsm p : open) {
			p.setInValidatedPeak(p.getPeakId() >= 0 && gate.acceptedPeakIds.contains(p.getPeakId()));
		}
		JointIdentification.Config cfg = new JointIdentification.Config();
		cfg.classifier.maxIterations = 5;
		JointIdentification.Result result = JointIdentification.run(open, closed, gate, cfg, null, null, null);
		assertNotNull(result.validation);
		EnumMap<JointCandidate.Group, int[]> fdp = new EnumMap<>(JointCandidate.Group.class);
		for (int i = 0; i < result.candidates.size(); i++) {
			JointCandidate c = result.candidates.get(i);
			if (!result.kept[i] || !c.isTarget()) {
				continue;
			}
			boolean right = c.getPeptide().equals(truth.get(c.spectrumKey()));
			fdp.computeIfAbsent(c.getGroup(), g -> new int[2])[right ? 0 : 1]++;
		}
		for (Map.Entry<JointCandidate.Group, int[]> e : fdp.entrySet()) {
			int total = e.getValue()[0] + e.getValue()[1];
			assertTrue(e.getValue()[1] <= Math.max(3, 0.03 * total),
					e.getKey() + ": " + e.getValue()[1] + " wrong of " + total + " kept");
		}
		assertTrue(fdp.containsKey(JointCandidate.Group.OPEN_UNMODIFIED),
				"open-only unmodified identifications exist when the closed search misses spectra");
		// Every row of a candidate carries the joint q-value: a winner's group q, 1.0 otherwise.
		for (int i = 0; i < result.candidates.size(); i++) {
			JointCandidate c = result.candidates.get(i);
			double expected = result.winner[i] ? result.psmQ[i] : 1.0;
			for (OpenSearchPsm p : new OpenSearchPsm[] { c.closed(), c.fragger(), c.sage() }) {
				if (p != null) {
					assertEquals(expected, p.getQValue(), 0.0);
				}
			}
		}
	}

	/**
	 * Confident candidates must not tie. Recovering the logit from the rounded target probability
	 * put half the winners of this run at one value (+27.63).
	 */
	@Test
	void confidentCandidatesKeepDistinctScores() throws Exception {
		Random r = new Random(11);
		String[] present = peptides(r, 800, 'K');
		String[] background = peptides(r, 50_000, 'R');
		List<OpenSearchPsm> open = new ArrayList<>();
		List<OpenSearchPsm> closed = new ArrayList<>();
		for (int s = 0; s < 3000; s++) {
			String file = "run" + (s % 3);
			String t = present[r.nextInt(present.length)];
			boolean identifiable = r.nextDouble() < 0.55;
			if (identifiable) {
				closed.add(psm(file, s, t, true, 1, 28 + 6 * r.nextGaussian(), true, 0, -1));
				open.add(psm(file, s, t, true, 1, 26 + 6 * r.nextGaussian(), false, 0.0, 0));
			} else {
				closed.add(psm(file, s, background[r.nextInt(background.length)], r.nextBoolean(), 1,
						13 + 4 * r.nextGaussian(), true, 0, -1));
			}
			for (int k = 0; k < 2; k++) {
				open.add(psm(file, s, background[r.nextInt(background.length)], r.nextBoolean(), 2 + k,
						13 + 4 * r.nextGaussian(), false, 0.0, 0));
			}
		}
		ModificationValidator.Gate gate = gate();
		for (OpenSearchPsm p : open) {
			p.setInValidatedPeak(true);
		}
		JointIdentification.Config cfg = new JointIdentification.Config();
		cfg.classifier.maxIterations = 3;
		JointIdentification.Result result = JointIdentification.run(open, closed, gate, cfg, null, null, null);
		int accepted = 0;
		java.util.Set<Double> distinct = new HashSet<>();
		for (int i = 0; i < result.candidates.size(); i++) {
			if (result.accepted[i]) {
				accepted++;
				distinct.add(result.logits[i]);
			}
		}
		assertTrue(accepted > 1000, "accepted " + accepted);
		assertTrue(distinct.size() > 0.99 * accepted, distinct.size() + " distinct scores among " + accepted);
	}

	// ------------------------------------------------------------ fixtures

	// ------------------------------------------------------------ identification v2

	@Test
	void precursorQueriesSeparateReassignedPrecursors() {
		OpenSearchPsm closed = psm("a", 1, "PEPTIDEK", true, 1, 31, true, 0.0, -1);
		OpenSearchPsm sameMass = psm("a", 1, "PEPTIDEK", true, 1, 30, false, 0.0, 0);
		OpenSearchPsm isotope = psm("a", 1, "LLLLLLLK", true, 2, 20, false, 0.0, 0);
		isotope.setPrecursorNeutralMass(closed.getPrecursorNeutralMass() + ClosedOpenCascade.ISOTOPE_SPACING);
		OpenSearchPsm otherCharge = psm("a", 1, "PEPTIDEK", true, 3, 18, false, 0.0, 0);
		otherCharge.setCharge(3);
		List<JointCandidate> merged = JointCandidate.merge(Arrays.asList(sameMass, isotope, otherCharge),
				Collections.singletonList(closed), new int[2], true);
		assertEquals(3, merged.size(), "same precursor merges with the closed row; the reassigned one does not");
		JointCandidate both = merged.get(0);
		assertTrue(both.foundByClosed() && both.foundByOpen() && !both.isReassigned());
		assertFalse(merged.get(1).isReassigned(), "an isotope offset is the same precursor");
		JointCandidate reassigned = merged.get(2);
		assertTrue(reassigned.isReassigned() && !reassigned.foundByClosed());
		assertEquals("a:1:R", reassigned.queryKey());
		assertEquals("a:1", reassigned.spectrumKey(), "the fold key stays the scan");
		// v1 merge is unchanged: all three rows of PEPTIDEK and the other peptide make two candidates.
		assertEquals(2, JointCandidate.merge(Arrays.asList(sameMass, isotope, otherCharge),
				Collections.singletonList(closed), new int[2]).size());
		OpenSearchPsm second = psm("a", 1, "QQQQQQQK", true, 1, 15, false, 0.0, 0);
		second.setCharge(4);
		org.junit.jupiter.api.Assertions.assertThrows(IllegalStateException.class,
				() -> JointCandidate.merge(Arrays.asList(otherCharge, second), Collections.singletonList(closed),
						new int[2], true), "a second reassigned precursor in one scan is refused");
	}

	@Test
	void weightedConservativeQValuesCountEachDecoyTwice() {
		double[] s = { 10, 9, 8, 7, 6, 5 };
		boolean[] t = { true, true, false, true, true, true };
		double[] plain = TargetDecoyFdr.qValuesConservative(s, t);
		double[] reset = TargetDecoyFdr.qValuesConservative(s, t, 2.0);
		assertEquals(2.0 / 5, plain[5], 1e-12);
		assertEquals(Math.min(1.0, 2 * 2.0 / 5), reset[5], 1e-12);
	}

	@Test
	void decoyWinsATiedQueryAndResetRolesAreStable() {
		JointCandidate target = candidate("PEPTIDEK", true, JointCandidate.Group.EXACT);
		JointCandidate decoy = candidate("KEDITPEP", false, JointCandidate.Group.EXACT);
		boolean[] w = JointIdentification.winners(new int[] { 0, 0 }, new double[] { 3.0, 3.0 },
				Arrays.asList(target, decoy));
		assertFalse(w[0]);
		assertTrue(w[1]);
		w = JointIdentification.winners(new int[] { 0, 0 }, new double[] { 3.0, 3.0 }, Arrays.asList(decoy, target));
		assertTrue(w[0], "whatever the row order");
		List<JointCandidate> many = new ArrayList<>();
		Random r = new Random(3);
		for (String p : peptides(r, 400, 'K')) {
			many.add(candidate(p, r.nextBoolean(), JointCandidate.Group.EXACT));
		}
		boolean[] roles = JointIdentification.resetRoles(many);
		int training = 0;
		int decoys = 0;
		for (int i = 0; i < roles.length; i++) {
			assertFalse(roles[i] && many.get(i).isTarget(), "only decoys train");
			training += roles[i] ? 1 : 0;
			decoys += many.get(i).isTarget() ? 0 : 1;
		}
		assertTrue(training > decoys / 4 && training < 3 * decoys / 4, training + " of " + decoys);
		assertTrue(Arrays.equals(roles, JointIdentification.resetRoles(many)), "deterministic");
	}

	@Test
	void precursorClassesAreSeparateFdrGroups() {
		List<JointCandidate> cs = new ArrayList<>();
		double[] scores = new double[40];
		boolean[] winner = new boolean[40];
		for (int i = 0; i < 40; i++) {
			JointCandidate c = candidate("PEP" + i + "K", i % 10 != 0, JointCandidate.Group.OPEN_UNMODIFIED);
			c.reassigned = i >= 20;
			cs.add(c);
			scores[i] = 100 - i;
			winner[i] = true;
		}
		double[] split = JointIdentification.groupQValues(cs, scores, winner, true, 1.0);
		double[] pooled = JointIdentification.groupQValues(cs, scores, winner, false, 1.0);
		// Each group ends at 2 decoys and 18 targets: (2 + 1) / 18; pooled, (4 + 1) / 36.
		assertEquals(3.0 / 18, split[20], 1e-12);
		assertEquals(3.0 / 18, split[0], 1e-12);
		assertEquals(5.0 / 36, pooled[39], 1e-12);
		assertTrue(pooled[20] < split[20]);
	}

	@Test
	void treeModelIsMonotoneAndScoresBatchesLikeRows() throws Exception {
		java.util.ArrayList<weka.core.Attribute> atts = new java.util.ArrayList<>();
		atts.add(new weka.core.Attribute("open_hyperscore"));
		atts.add(new weka.core.Attribute("noise"));
		atts.add(new weka.core.Attribute("label", Arrays.asList("target", "decoy")));
		weka.core.Instances data = new weka.core.Instances("t", atts, 4000);
		data.setClassIndex(2);
		Random r = new Random(5);
		for (int i = 0; i < 4000; i++) {
			boolean target = r.nextBoolean();
			double hs = (target ? 25 : 15) + 5 * r.nextGaussian();
			data.add(new weka.core.DenseInstance(1.0, new double[] { hs, r.nextGaussian(), target ? 0 : 1 }));
		}
		XGBoostTreeClassifier model = new XGBoostTreeClassifier(java.util.Set.of("open_hyperscore"));
		model.rounds = 50;
		model.buildClassifier(data);
		float[] grid = new float[2 * 60];
		for (int k = 0; k < 60; k++) {
			grid[2 * k] = 5 + k * 0.5f;
			grid[2 * k + 1] = 0f;
		}
		double[] logits = model.targetLogits(grid, 60);
		for (int k = 1; k < 60; k++) {
			assertTrue(logits[k] >= logits[k - 1] - 1e-9, "monotone in the constrained score at " + k);
		}
		assertTrue(logits[59] > logits[0] + 1, "the score separates targets from decoys");
		weka.core.Instance one = data.instance(0);
		double p = model.distributionForInstance(one)[0];
		double batch = model.targetLogits(new float[] { (float) one.value(0), (float) one.value(1) }, 1)[0];
		assertEquals(batch, Math.log(p / (1 - p)), 1e-4);
	}

	/** The v2 joint identification (precursor queries, RESET, trees) on the synthetic run of {@link #groupFdrHoldsOnASyntheticRun}. */
	@Test
	void identificationV2ControlsGroupFdrOnASyntheticRun() throws Exception {
		Random r = new Random(11);
		String[] present = peptides(r, 800, 'K');
		String[] background = peptides(r, 50_000, 'R');
		List<OpenSearchPsm> open = new ArrayList<>();
		List<OpenSearchPsm> closed = new ArrayList<>();
		Map<String, String> truth = new HashMap<>();
		for (int s = 0; s < 6000; s++) {
			String file = "run" + (s % 3);
			double u = r.nextDouble();
			String t = present[r.nextInt(present.length)];
			int type = u < 0.55 ? 0 : u < 0.62 ? 1 : 3;
			if (type <= 1) {
				truth.put(file + ":" + s, t);
			}
			if (type == 0 && r.nextDouble() < 0.7) {
				closed.add(psm(file, s, t, true, 1, 28 + 6 * r.nextGaussian(), true, 0, -1));
			} else if (r.nextDouble() < 0.8) {
				closed.add(psm(file, s, background[r.nextInt(background.length)], r.nextBoolean(), 1,
						13 + 4 * r.nextGaussian(), true, 0, -1));
			}
			int rank = 1;
			int first = open.size();
			if (type == 0) {
				open.add(psm(file, s, t, true, rank++, 25 + 6 * r.nextGaussian(), false, 0.0, 0));
			} else if (type == 1) {
				open.add(psm(file, s, t, true, rank++, 24 + 6 * r.nextGaussian(), false, 15.9949, 1));
			}
			for (int k = 0; k < 3; k++) {
				double dm = r.nextDouble() < 0.3 ? 0.0 : (r.nextDouble() - 0.5) * 400;
				OpenSearchPsm b = psm(file, s, background[r.nextInt(background.length)], r.nextBoolean(), rank++,
						13 + 4 * r.nextGaussian(), false, dm, dm == 0.0 ? 0 : -1);
				if (k == 2 && s % 7 == 0) {
					b.setCharge(3); // a co-isolated precursor Crystal-C reassigned the scan to
				}
				open.add(b);
			}
			// The open rows of a scan share its observed precursor.
			double observed = open.get(first).getPrecursorNeutralMass();
			for (int k = first; k < open.size(); k++) {
				open.get(k).setPrecursorNeutralMass(observed);
			}
			if (!closed.isEmpty() && closed.get(closed.size() - 1).getScan() == s) {
				closed.get(closed.size() - 1).setPrecursorNeutralMass(observed);
			}
		}
		ModificationValidator.Gate gate = gate();
		for (OpenSearchPsm p : open) {
			p.setInValidatedPeak(p.getPeakId() >= 0 && gate.acceptedPeakIds.contains(p.getPeakId()));
		}
		JointIdentification.Config cfg = new JointIdentification.Config();
		cfg.classifier.maxIterations = 4;
		cfg.precursorQueries = true;
		cfg.reset = true;
		cfg.classifier.unclampedLogits = true;
		cfg.classifier.modelFactory = new SemiSupervisedPsmClassifier.ModelFactory() {
			@Override
			public weka.classifiers.Classifier create() {
				XGBoostTreeClassifier m = new XGBoostTreeClassifier(
						java.util.Set.of("closed_hyperscore", "open_hyperscore", "closed_neg_log10_expect",
								"open_neg_log10_expect", "hyperscore_margin"));
				m.rounds = 60;
				return m;
			}

			@Override
			public String describe() {
				return "test trees";
			}
		};
		JointIdentification.Result result = JointIdentification.run(open, closed, gate, cfg, null, null, null);
		assertNotNull(result.trainingDecoy);
		assertEquals(2.0, result.decoyWeight, 0.0);
		int reassignedQueries = 0;
		EnumMap<JointCandidate.Group, int[]> fdp = new EnumMap<>(JointCandidate.Group.class);
		for (int i = 0; i < result.candidates.size(); i++) {
			JointCandidate c = result.candidates.get(i);
			assertFalse(result.winner[i] && result.trainingDecoy[i], "a training decoy is never reported");
			reassignedQueries += result.winner[i] && c.isReassigned() ? 1 : 0;
			if (!result.kept[i] || !c.isTarget() || c.isReassigned()) {
				continue;
			}
			boolean right = c.getPeptide().equals(truth.get(c.spectrumKey()));
			fdp.computeIfAbsent(c.getGroup(), g -> new int[2])[right ? 0 : 1]++;
		}
		assertTrue(reassignedQueries > 0, "reassigned precursors are their own queries");
		for (Map.Entry<JointCandidate.Group, int[]> e : fdp.entrySet()) {
			int total = e.getValue()[0] + e.getValue()[1];
			assertTrue(e.getValue()[1] <= Math.max(3, 0.03 * total),
					e.getKey() + ": " + e.getValue()[1] + " wrong of " + total + " kept");
		}
		assertTrue(fdp.containsKey(JointCandidate.Group.EXACT) && fdp.get(JointCandidate.Group.EXACT)[0] > 1000,
				"the exact group identifies most present peptides");
	}

	private static ModificationValidator.Gate gate() {
		List<ModificationValidator.PeakVerdict> verdicts = new ArrayList<>();
		verdicts.add(new ModificationValidator.PeakVerdict(0, 0.0, ModificationValidator.MassClass.ZERO, 0, 0, 0,
				0, Double.NaN, false, true, "", 0));
		verdicts.add(new ModificationValidator.PeakVerdict(1, 15.9949, ModificationValidator.MassClass.ESTABLISHED,
				0, 0, 0, 0, 0, false, true, "", 0));
		return new ModificationValidator.Gate(new HashSet<>(Arrays.asList(0, 1)), verdicts, Collections.emptySet());
	}

	private static JointCandidate candidate(String peptide, boolean target, JointCandidate.Group group) {
		JointCandidate c = new JointCandidate("f", peptide.hashCode(), peptide, target);
		c.group = group;
		return c;
	}

	private static String[] peptides(Random r, int n, char last) {
		String[] out = new String[n];
		for (int i = 0; i < n; i++) {
			StringBuilder sb = new StringBuilder();
			int len = 8 + r.nextInt(10);
			for (int k = 0; k < len - 1; k++) {
				sb.append(AA.charAt(r.nextInt(AA.length())));
			}
			out[i] = sb.append(last).toString();
		}
		return out;
	}

	private static OpenSearchPsm psm(String file, int scan, String peptide, boolean target, int rank, double hyper,
			boolean closedStage, double dm, int peak) {
		OpenSearchPsm p = new OpenSearchPsm();
		p.setFileName(file);
		p.setScan(scan);
		p.setPeptide(target ? peptide : new StringBuilder(peptide).reverse().toString());
		p.setTarget(target);
		p.setHitRank(rank);
		p.setHyperscore(hyper);
		p.setNextscore(hyper - 2);
		p.setExpect(Math.pow(10, -(hyper - (closedStage ? 8 : 12)) / 4.0));
		p.setCharge(2);
		p.setNumMatchedIons((int) Math.max(1, hyper / 2));
		p.setTotNumIons(30);
		p.setCalcNeutralPepMass(1000 + peptide.length() * 110);
		p.setPrecursorNeutralMass(1000 + peptide.length() * 110 + dm);
		p.setMassDiff(dm);
		p.setNormalizedMassDiff(dm);
		p.setPeakId(peak);
		p.setClosedStage(closedStage);
		p.setNumTolTerm(2);
		p.setNumMissedCleavages(0);
		p.setProtein((target ? "" : "rev_") + "MGYG000000001_" + Math.floorMod(peptide.hashCode(), 1000));
		return p;
	}
}
