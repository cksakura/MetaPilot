package uk.ac.quadram.metapilot.core.unified;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

import org.junit.jupiter.api.Test;

import uk.ac.quadram.metapilot.core.opensearch.TargetDecoyFdr;
import uk.ac.quadram.metapilot.core.variant.MzmlScanReader;

class UnifiedIdentificationTest {

	@Test
	void groupWalkWithOneGroupIsTargetDecoyCompetitionWithPlusOne() {
		Random r = new Random(3);
		int n = 500;
		double[] s = new double[n];
		boolean[] t = new boolean[n];
		int[] g = new int[n];
		for (int i = 0; i < n; i++) {
			t[i] = r.nextDouble() < 0.7;
			s[i] = (t[i] ? 1.5 : 0) + r.nextGaussian();
		}
		double[] walk = GroupWalk.qValues(s, t, g, 1, 40, 1.0, 1.0);
		double[] tdc = TargetDecoyFdr.qValuesConservative(s, t);
		for (int i = 0; i < n; i++) {
			if (t[i]) {
				assertEquals(tdc[i], walk[i], 1e-12, "row " + i);
			}
		}
	}

	/**
	 * Two groups, one clean and one mostly null. A false target is as likely as a decoy, so the
	 * share of reported targets that are false (FDP), averaged over repetitions, must stay at alpha.
	 */
	@Test
	void groupWalkControlsTheFdrAcrossGroups() {
		Random r = new Random(7);
		double alpha = 0.05;
		double fdpSum = 0;
		int reps = 200;
		for (int rep = 0; rep < reps; rep++) {
			List<double[]> items = new ArrayList<>(); // {score, target, group, false}
			for (int i = 0; i < 1500; i++) {
				boolean group1 = i % 3 == 0;
				boolean correct = r.nextDouble() < (group1 ? 0.1 : 0.6);
				if (correct) {
					items.add(new double[] { 2.5 + r.nextGaussian(), 1, group1 ? 1 : 0, 0 });
				} else {
					boolean target = r.nextBoolean();
					items.add(new double[] { r.nextGaussian(), target ? 1 : 0, group1 ? 1 : 0, 1 });
				}
			}
			double[] s = new double[items.size()];
			boolean[] t = new boolean[items.size()];
			int[] g = new int[items.size()];
			for (int i = 0; i < s.length; i++) {
				s[i] = items.get(i)[0];
				t[i] = items.get(i)[1] == 1;
				g[i] = (int) items.get(i)[2];
			}
			double[] q = GroupWalk.qValues(s, t, g, 2, 40, 1.0, 1.0);
			int reported = 0;
			int falseReported = 0;
			for (int i = 0; i < s.length; i++) {
				if (t[i] && q[i] <= alpha) {
					reported++;
					if (items.get(i)[3] == 1) {
						falseReported++;
					}
				}
			}
			fdpSum += reported == 0 ? 0 : falseReported / (double) reported;
		}
		double meanFdp = fdpSum / reps;
		assertTrue(meanFdp <= alpha + 0.01, "mean FDP " + meanFdp);
	}

	@Test
	void listwiseModelLearnsWhichCandidateIsRightAndWhenNoneIs() {
		Random r = new Random(1);
		int queries = 4000;
		List<float[]> rows = new ArrayList<>();
		int[] qStart = new int[queries + 1];
		int[] positive = new int[queries];
		for (int q = 0; q < queries; q++) {
			qStart[q] = rows.size();
			boolean identifiable = q % 2 == 0;
			int right = identifiable ? r.nextInt(3) : -1;
			for (int c = 0; c < 3; c++) {
				float quality = (float) (c == right ? 2.0 + r.nextGaussian() * 0.5 : r.nextGaussian() * 0.5);
				rows.add(new float[] { quality, (float) r.nextGaussian() });
				if (c == right) {
					positive[q] = rows.size() - 1;
				}
			}
			if (!identifiable) {
				positive[q] = -1;
			}
		}
		qStart[queries] = rows.size();
		float[][] x = rows.toArray(new float[0][]);
		ListwiseGbdt.Params p = new ListwiseGbdt.Params();
		p.trees = 60;
		ListwiseGbdt.Model model = ListwiseGbdt.train(x, qStart, positive, new int[] { 1, 0 }, p);
		int correct = 0;
		double noneUtility = 0;
		int none = 0;
		for (int q = 0; q < queries; q++) {
			int best = qStart[q];
			for (int row = qStart[q]; row < qStart[q + 1]; row++) {
				if (model.utility(x[row]) > model.utility(x[best])) {
					best = row;
				}
			}
			if (positive[q] >= 0 && best == positive[q]) {
				correct++;
			}
			if (positive[q] < 0) {
				noneUtility += model.utility(x[best]);
				none++;
			}
		}
		assertTrue(correct > 0.95 * queries / 2, "picked the right candidate in " + correct);
		assertTrue(noneUtility / none < 0, "the best candidate of an unidentifiable spectrum sits below 'none'");
		// Monotone in feature 0.
		double prev = Double.NEGATIVE_INFINITY;
		for (float v = -2; v <= 4; v += 0.25f) {
			double u = model.utility(new float[] { v, 0 });
			assertTrue(u >= prev - 1e-12, "monotone at " + v);
			prev = u;
		}
	}

	@Test
	void annotatorPlacesTheDeltaOnTheResidueThatExplainsTheSpectrum() {
		String peptide = "PEPTMDEKAR";
		double delta = 15.9949;
		int site = 4; // M
		SpectrumAnnotator.Ladder ladder = new SpectrumAnnotator.Ladder(peptide, 2);
		double[] mz = new double[ladder.size()];
		for (int k = 0; k < ladder.size(); k++) {
			mz[k] = ladder.mz(k, ladder.contains(k, site) ? delta : 0.0);
		}
		Arrays.sort(mz);
		float[] in = new float[mz.length];
		Arrays.fill(in, 100f);
		SpectrumAnnotator.Spectrum s = SpectrumAnnotator.prepare(new MzmlScanReader.Peaks(mz, in));
		SpectrumAnnotator annotator = new SpectrumAnnotator(20);
		SpectrumAnnotator.Annotation a = annotator.annotate(s, peptide, 2, null, delta);
		assertEquals(site, a.bestSite);
		assertEquals(1.0, a.eIntTic, 1e-6);
		assertTrue(a.uIntTic < a.eIntTic);
		assertTrue(a.shiftedIonsGained > 0);
		// The closed search's variable modification at the same residue explains it exactly too.
		double[] mods = new double[peptide.length() + 1];
		mods[site + 1] = delta;
		assertEquals(1.0, annotator.annotate(s, peptide, 2, mods, 0).eIntTic, 1e-6);
	}

	@Test
	void isoleucineLeucineVariantsAreNeighbours() {
		double[] a = SpectrumAnnotator.theoreticalIons("PEPTIDEK", 2, null, 0, -1);
		double[] b = SpectrumAnnotator.theoreticalIons("PEPTLDEK", 2, null, 0, -1);
		assertEquals(1.0, SpectrumAnnotator.sharedIonFraction(a, b, 20), 1e-12);
		double[] c = SpectrumAnnotator.theoreticalIons("QWERTYAK", 2, null, 0, -1);
		assertTrue(SpectrumAnnotator.sharedIonFraction(a, c, 20) < 0.3);
	}

	@Test
	void batchesWithIndistinguishableScoresMergeAndSmallOnesJoinTheMostSimilar() {
		Random r = new Random(9);
		List<double[]> batches = new ArrayList<>();
		batches.add(gaussian(r, 5000, 3.0)); // 0: exact, always its own group
		batches.add(gaussian(r, 2000, 0.0)); // 1
		batches.add(gaussian(r, 1500, 0.0)); // 2: same distribution as 1 -> merged
		batches.add(gaussian(r, 1500, 2.0)); // 3: different -> own group
		batches.add(gaussian(r, 30, 2.0)); // 4: too small -> joins 3, its most similar group
		int[] g = GroupWalk.mergeBatches(batches, 0, 0.01, 80);
		assertEquals(0, g[0]);
		assertEquals(g[1], g[2]);
		assertTrue(g[3] != g[1] && g[3] != 0);
		assertEquals(g[3], g[4]);
	}

	private static double[] gaussian(Random r, int n, double mean) {
		double[] v = new double[n];
		for (int i = 0; i < n; i++) {
			v[i] = mean + r.nextGaussian();
		}
		return v;
	}

	@Test
	void binomialTail() {
		assertEquals(1.0, DeltaMassBatches.binomialUpperTail(0, 10, 0.3), 1e-12);
		assertEquals(Math.pow(0.5, 10), DeltaMassBatches.binomialUpperTail(10, 10, 0.5), 1e-9);
		assertEquals(0.5, DeltaMassBatches.binomialUpperTail(1, 1, 0.5), 1e-9);
		assertTrue(DeltaMassBatches.binomialUpperTail(60, 100, 0.2) < 1e-10);
	}

	@Test
	void openCandidateMergesIntoAClosedOneOnlyWhenBothExplainTheSameMass() {
		double[] oxM2 = { 0, 0, 15.9949, 0, 0, 0 };
		Hypothesis closedOx = hyp(0, "PEMTLDK", true, 1, oxM2, 0.0);
		Hypothesis closedPlain = hyp(0, "PEMTLDK", true, 2, null, 1.00335);
		Hypothesis openOx = hyp(0, "PEMTLDK", true, -1, null, 15.9952);
		Hypothesis openIso = hyp(0, "PEMTLDK", true, -1, null, 1.0040);
		Hypothesis openOther = hyp(0, "PEMTLDK", true, -1, null, 42.0106);
		Hypothesis openElsewhere = hyp(1, "PEMTLDK", true, -1, null, 15.9949);
		List<Hypothesis> hyps = new ArrayList<>(
				Arrays.asList(closedOx, closedPlain, openOx, openIso, openOther, openElsewhere));
		assertEquals(2, UnifiedIdentification.mergeEquivalent(hyps));
		assertEquals(Arrays.asList(closedOx, closedPlain, openOther, openElsewhere), hyps);
		assertTrue(closedOx.foundByOpen && closedPlain.foundByOpen);
		assertEquals(15.9952, closedOx.openMassDiff, 1e-9);
		assertEquals(1.0040, closedPlain.openMassDiff, 1e-9);
		// A distinct peptidoform (acetylation) stays a candidate of its own, not "found by both".
		assertTrue(!openOther.foundByClosed);
	}

	private static Hypothesis hyp(int spectrum, String peptide, boolean target, int closedRank, double[] mods,
			double massDiff) {
		Hypothesis h = new Hypothesis(spectrum, peptide, target, 2, 800, 600, 2, 0, new String[] { "P1" });
		if (closedRank > 0) {
			h.foundByClosed = true;
			h.closedRank = closedRank;
			h.closedMods = mods;
			h.closedMassDiff = massDiff;
		} else {
			h.foundByOpen = true;
			h.openRank = 1;
			h.openMassDiff = massDiff;
		}
		return h;
	}
}
