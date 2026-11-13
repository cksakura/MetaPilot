package uk.ac.quadram.metapilot.core.opensearch;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.Random;
import java.util.Arrays;
import java.util.List;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

import uk.ac.quadram.metapilot.core.aminoacid.AminoAcidProperty;

/**
 * Ground-truth tests for {@link DeltaMassPeakFitter}: synthetic delta-mass
 * distributions with planted peaks of known centroid, width and population, so
 * the fit can be checked against values it cannot have derived circularly.
 */
class DeltaMassPeakFitterTest {

	/** Systematic calibration offset planted in the data (Da). */
	private static final double CALIBRATION_OFFSET = 0.0012;
	private static final double OXIDATION = 15.9949;
	private static final double PHOSPHO = 79.96633;

	private static final long SEED = 42L;

	/**
	 * Zero peak, its +1 isotope artefact, two real modifications, and a uniform
	 * background standing in for random/decoy matches.
	 */
	private static double[] syntheticDeltaMasses() {
		Random rng = new Random(SEED);
		int zeroN = 40000, isoN = 4000, oxN = 6000, phN = 1500, bgN = 30000;
		double[] out = new double[zeroN + isoN + oxN + phN + bgN];
		int k = 0;
		for (int i = 0; i < zeroN; i++) {
			out[k++] = CALIBRATION_OFFSET + rng.nextGaussian() * 0.0030;
		}
		for (int i = 0; i < isoN; i++) {
			out[k++] = CALIBRATION_OFFSET + AminoAcidProperty.isotope_deltaM + rng.nextGaussian() * 0.0032;
		}
		for (int i = 0; i < oxN; i++) {
			out[k++] = CALIBRATION_OFFSET + OXIDATION + rng.nextGaussian() * 0.0040;
		}
		for (int i = 0; i < phN; i++) {
			out[k++] = CALIBRATION_OFFSET + PHOSPHO + rng.nextGaussian() * 0.0045;
		}
		for (int i = 0; i < bgN; i++) {
			out[k++] = -150.0 + rng.nextDouble() * 650.0;
		}
		return out;
	}

	private static DeltaMassPeakFitter.Peak nearest(DeltaMassPeakFitter.Result result, double expected) {
		DeltaMassPeakFitter.Peak best = null;
		double bestDistance = Double.MAX_VALUE;
		for (DeltaMassPeakFitter.Peak peak : result.peaks) {
			double distance = Math.abs(peak.centroid - expected);
			if (distance < bestDistance) {
				bestDistance = distance;
				best = peak;
			}
		}
		return bestDistance < 0.05 ? best : null;
	}

	private static DeltaMassPeakFitter.Config gaussianConfig() {
		DeltaMassPeakFitter.Config cfg = new DeltaMassPeakFitter.Config();
		cfg.useKde = false;
		return cfg;
	}

	@Test
	@DisplayName("Savitzky-Golay coefficients match the textbook 5-point quadratic filter")
	void savitzkyGolayMatchesTextbook() {
		double[] c = DeltaMassPeakFitter.savitzkyGolayCoefficients(2, 2);
		double[] expected = { -3.0 / 35, 12.0 / 35, 17.0 / 35, 12.0 / 35, -3.0 / 35 };
		assertEquals(expected.length, c.length);
		for (int i = 0; i < expected.length; i++) {
			assertEquals(expected[i], c[i], 1e-12, "coefficient " + i);
		}
		double sum = 0;
		for (double v : c) {
			sum += v;
		}
		// A smoothing filter must preserve a constant signal.
		assertEquals(1.0, sum, 1e-12, "coefficients must sum to 1");
	}

	@Test
	@DisplayName("Zero peak recovers the systematic calibration offset")
	void zeroPeakRecoversCalibrationOffset() {
		DeltaMassPeakFitter.Result result = DeltaMassPeakFitter.fit(syntheticDeltaMasses(),
				gaussianConfig());

		assertNotEquals(-1, result.zeroPeakId, "a zero peak should be found");
		DeltaMassPeakFitter.Peak zero = result.peaks.get(result.zeroPeakId);

		// This is why stage 1 fits rather than using a window centred on 0.0:
		// a fixed window would clip one side of an offset peak.
		assertEquals(CALIBRATION_OFFSET, zero.centroid, 0.0005, "calibration offset should be recovered");
		assertEquals(CALIBRATION_OFFSET, result.calibrationOffset, 0.0005);
		assertEquals(0.0030, zero.sigma, 0.0008, "fitted sigma should match the planted width");
		assertTrue(zero.r2 >= 0.9, "zero peak R2 should be high, was " + zero.r2);
	}

	@Test
	@DisplayName("Stage-1 zero-family discovery uses KDE by default")
	void stageOneZeroFamilyUsesKde() {
		Random rng = new Random(SEED);
		double isotope = AminoAcidProperty.isotope_deltaM;
		double[] masses = new double[1230];
		boolean[] targets = new boolean[masses.length];
		int k = 0;
		for (int i = 0; i < 900; i++) {
			masses[k] = CALIBRATION_OFFSET + rng.nextGaussian() * 0.003;
			targets[k++] = true;
		}
		for (int i = 0; i < 250; i++) {
			masses[k] = CALIBRATION_OFFSET + isotope + rng.nextGaussian() * 0.003;
			targets[k++] = true;
		}
		for (int i = 0; i < 80; i++) {
			masses[k] = -0.2 + rng.nextDouble() * 0.4;
			targets[k++] = false;
		}

		DeltaMassPeakFitter.Result result = DeltaMassPeakFitter.fitZeroPeakOnly(masses, targets,
				new DeltaMassPeakFitter.Config(), null);

		assertEquals(0, result.zeroPeakId);
		assertEquals("KDE", result.peaks.get(0).discoveryMethod);
		assertTrue(Double.isNaN(result.peaks.get(0).r2), "KDE discovery must not depend on Gaussian R2");
		assertEquals(CALIBRATION_OFFSET, result.calibrationOffset, 0.001);
	}

	@Test
	@DisplayName("Stage-1 KDE normalizes isotope PSMs without creating modification peaks")
	void stageOneKdeNormalizesOnlyTheZeroIsotopeFamily() {
		Random rng = new Random(SEED);
		double isotope = AminoAcidProperty.isotope_deltaM;
		double[] masses = new double[1200];
		boolean[] targets = new boolean[masses.length];
		int k = 0;
		for (int i = 0; i < 800; i++) {
			masses[k] = CALIBRATION_OFFSET + rng.nextGaussian() * 0.0025;
			targets[k++] = true;
		}
		int firstIsotope = k;
		for (int i = 0; i < 250; i++) {
			masses[k] = CALIBRATION_OFFSET + isotope + rng.nextGaussian() * 0.0025;
			targets[k++] = true;
		}
		int firstDeamidation = k;
		for (int i = 0; i < 150; i++) {
			masses[k] = 0.984016 + rng.nextGaussian() * 0.0025;
			targets[k++] = true;
		}

		DeltaMassPeakFitter.Result result = DeltaMassPeakFitter.fitZeroPeakOnly(masses, targets,
				new DeltaMassPeakFitter.Config(), null);

		assertEquals(1, result.peaks.size(), "isotope observations must not become separate modification peaks");
		assertEquals(0, result.familyAssignment[firstIsotope]);
		assertEquals(1, result.isotopeOffset[firstIsotope]);
		assertEquals(-1, result.observedAssignment[firstIsotope]);
		assertEquals(CALIBRATION_OFFSET, result.normalizedMassDiff[firstIsotope], 0.01);
		assertEquals(-1, result.familyAssignment[firstDeamidation],
				"deamidation near +1 Da must remain outside the zero isotope family");
	}

	@Test
	@DisplayName("Stage-1 folding stops at +2 by default and follows round1MaxIsotopeOffset")
	void stageOneIsotopeFoldingIsBoundedSeparatelyFromStageTwo() {
		Random rng = new Random(SEED);
		double isotope = AminoAcidProperty.isotope_deltaM;
		double[] masses = new double[1250];
		boolean[] targets = new boolean[masses.length];
		int k = 0;
		for (int i = 0; i < 900; i++) {
			masses[k] = CALIBRATION_OFFSET + rng.nextGaussian() * 0.0025;
			targets[k++] = true;
		}
		int firstSecond = k;
		for (int i = 0; i < 200; i++) {
			masses[k] = CALIBRATION_OFFSET + 2 * isotope + rng.nextGaussian() * 0.0025;
			targets[k++] = true;
		}
		int firstThird = k;
		for (int i = 0; i < 150; i++) {
			masses[k] = CALIBRATION_OFFSET + 3 * isotope + rng.nextGaussian() * 0.0025;
			targets[k++] = true;
		}

		DeltaMassPeakFitter.Config defaults = new DeltaMassPeakFitter.Config();
		assertEquals(2, defaults.round1MaxIsotopeOffset, "stage 1 must not inherit the stage-2 bound of 3");

		DeltaMassPeakFitter.Result bounded = DeltaMassPeakFitter.fitZeroPeakOnly(masses, targets, defaults, null);
		assertEquals(0, bounded.familyAssignment[firstSecond], "+2 is measured and must still be folded");
		assertEquals(2, bounded.isotopeOffset[firstSecond]);
		assertEquals(-1, bounded.familyAssignment[firstThird],
				"+3 is unmeasured and must be left out of the stage-1 set");
		assertEquals(0, bounded.isotopeOffset[firstThird], "an unfolded PSM carries no isotope offset");

		// The bound is a knob, not a hard-coded 2: widening it readmits +3.
		DeltaMassPeakFitter.Config widened = new DeltaMassPeakFitter.Config();
		widened.round1MaxIsotopeOffset = 3;
		DeltaMassPeakFitter.Result wide = DeltaMassPeakFitter.fitZeroPeakOnly(masses, targets, widened, null);
		assertEquals(0, wide.familyAssignment[firstThird]);
		assertEquals(3, wide.isotopeOffset[firstThird]);

		// ...and it can never exceed what stage 2 would be willing to link.
		DeltaMassPeakFitter.Config clamped = new DeltaMassPeakFitter.Config();
		clamped.maxIsotopeOffset = 1;
		clamped.round1MaxIsotopeOffset = 3;
		DeltaMassPeakFitter.Result narrow = DeltaMassPeakFitter.fitZeroPeakOnly(masses, targets, clamped, null);
		assertEquals(-1, narrow.familyAssignment[firstSecond],
				"round1MaxIsotopeOffset must be clamped to maxIsotopeOffset");
	}

	@Test
	@DisplayName("+1 isotope peak is labelled as an isotope, not reported as a modification")
	void isotopePeakIsLabelled() {
		DeltaMassPeakFitter.Result result = DeltaMassPeakFitter.fit(syntheticDeltaMasses(),
				new DeltaMassPeakFitter.Config());

		DeltaMassPeakFitter.Peak isotope = nearest(result,
				CALIBRATION_OFFSET + AminoAcidProperty.isotope_deltaM);
		assertNotNull(isotope, "the +1 isotope peak should be fitted");
		assertEquals(1, isotope.isotopeOffset, "should be linked as a +1 isotope");
		assertEquals(result.zeroPeakId, isotope.monoPeakId, "should be linked to the zero peak");
		assertTrue(isotope.isIsotope());

		// The observed satellite mode remains auditable after normalization.
		assertNotEquals(result.zeroPeakId, isotope.id);
		double[] masses = syntheticDeltaMasses();
		int isotopeMember = -1;
		for (int i = 0; i < masses.length; i++) {
			if (result.observedAssignment[i] == isotope.id) {
				isotopeMember = i;
				break;
			}
		}
		assertTrue(isotopeMember >= 0);
		assertEquals(result.zeroPeakId, result.familyAssignment[isotopeMember]);
		assertEquals(1, result.isotopeOffset[isotopeMember]);
		assertEquals(CALIBRATION_OFFSET, result.normalizedMassDiff[isotopeMember], 0.015);
	}

	@Test
	@DisplayName("Zero peak is the canonical root even when a -1 Da artefact exists")
	void zeroPeakCannotBecomeAnIsotopeOfNegativeArtefact() {
		double isotope = AminoAcidProperty.isotope_deltaM;
		List<DeltaMassPeakFitter.Peak> peaks = Arrays.asList(
				new DeltaMassPeakFitter.Peak(0, -isotope, 0.002, 20, 0.99, 20, 0, -1),
				new DeltaMassPeakFitter.Peak(1, 0.0001, 0.002, 100, 0.99, 100, 0, -1),
				new DeltaMassPeakFitter.Peak(2, isotope, 0.002, 40, 0.99, 40, 0, -1),
				new DeltaMassPeakFitter.Peak(3, 2 * isotope, 0.002, 10, 0.99, 10, 0, -1));

		List<DeltaMassPeakFitter.Peak> linked = DeltaMassPeakFitter.linkIsotopes(peaks,
				new DeltaMassPeakFitter.Config(), 1);

		assertEquals(0, linked.get(1).isotopeOffset);
		assertEquals(-1, linked.get(1).monoPeakId);
		assertFalse(linked.get(1).isIsotope(), "zero must remain the family root");
		assertEquals(-1, linked.get(0).isotopeOffset, "negative artefact is related to, not the root of, zero");
		assertEquals(1, linked.get(0).monoPeakId);
		assertEquals(1, linked.get(2).isotopeOffset);
		assertEquals(1, linked.get(2).monoPeakId);
		assertEquals(2, linked.get(3).isotopeOffset);
		assertEquals(1, linked.get(3).monoPeakId, "family members point directly to the canonical root");
	}

	@Test
	@DisplayName("Isotope families never grow through transitive satellite chains")
	void isotopeFamiliesAreDirectAndBounded() {
		double isotope = AminoAcidProperty.isotope_deltaM;
		List<DeltaMassPeakFitter.Peak> peaks = Arrays.asList(
				new DeltaMassPeakFitter.Peak(0, 39.0, 0.002, 100, 0.99, 100, 0, -1),
				new DeltaMassPeakFitter.Peak(1, 39.0 + 2 * isotope, 0.002, 50, 0.99, 50, 0, -1),
				new DeltaMassPeakFitter.Peak(2, 39.0 + 4 * isotope, 0.002, 20, 0.99, 20, 0, -1));

		List<DeltaMassPeakFitter.Peak> linked = DeltaMassPeakFitter.linkIsotopes(peaks,
				new DeltaMassPeakFitter.Config(), -1);

		assertEquals(0, linked.get(0).isotopeOffset);
		assertEquals(0, linked.get(1).isotopeOffset,
				"a lone +2 spacing without a +1 mode is insufficient family evidence");
		assertEquals(-1, linked.get(1).monoPeakId);
		assertEquals(0, linked.get(2).isotopeOffset,
				"an offset beyond the configured range must not be reached through a satellite");
		assertEquals(-1, linked.get(2).monoPeakId);
	}

	@Test
	@DisplayName("Real modification peaks are found and are not mistaken for isotopes")
	void modificationPeaksAreFound() {
		DeltaMassPeakFitter.Result result = DeltaMassPeakFitter.fit(syntheticDeltaMasses(),
				new DeltaMassPeakFitter.Config());

		DeltaMassPeakFitter.Peak oxidation = nearest(result, CALIBRATION_OFFSET + OXIDATION);
		assertNotNull(oxidation, "oxidation peak should be fitted");
		assertEquals(CALIBRATION_OFFSET + OXIDATION, oxidation.centroid, 0.001);
		assertEquals(0, oxidation.isotopeOffset, "oxidation is not an isotope artefact");

		DeltaMassPeakFitter.Peak phospho = nearest(result, CALIBRATION_OFFSET + PHOSPHO);
		assertNotNull(phospho, "phospho peak should be fitted");
		assertEquals(CALIBRATION_OFFSET + PHOSPHO, phospho.centroid, 0.001);
		assertEquals(0, phospho.isotopeOffset);
	}

	@Test
	@DisplayName("Membership at +/-3 sigma captures the expected Gaussian coverage")
	void membershipCoverageMatchesThreeSigma() {
		double[] deltaMasses = syntheticDeltaMasses();
		DeltaMassPeakFitter.Result result = DeltaMassPeakFitter.fit(deltaMasses, gaussianConfig());

		int assignedToZero = 0;
		for (int observed : result.observedAssignment) {
			if (observed == result.zeroPeakId) {
				assignedToZero++;
			}
		}
		// 40000 planted; +/-3 sigma covers 99.73% of a Gaussian.
		double coverage = assignedToZero / 40000.0;
		assertTrue(coverage > 0.985 && coverage < 1.02,
				"zero-peak membership should be near the 99.7% three-sigma coverage, was " + coverage);
	}

	@Test
	@DisplayName("Nothing is discarded: every input is scored, out-of-peak PSMs get -1 rather than being dropped")
	void nothingIsDiscarded() {
		double[] deltaMasses = syntheticDeltaMasses();
		DeltaMassPeakFitter.Result result = DeltaMassPeakFitter.fit(deltaMasses, new DeltaMassPeakFitter.Config());

		// The pipeline's governing principle: peak membership is a classifier
		// feature, never a filter. Every input must survive with a feature value.
		assertEquals(deltaMasses.length, result.assignment.length);
		assertEquals(deltaMasses.length, result.distanceToNearestCentroid.length);
		assertEquals(deltaMasses.length, result.distanceInSigmas.length);
		assertEquals(deltaMasses.length, result.localDensity.length);

		int outsideAnyPeak = 0;
		for (int a : result.assignment) {
			if (a < 0) {
				outsideAnyPeak++;
			}
		}
		// The uniform background is mostly outside every peak; those PSMs are
		// retained with assignment == -1 and still reach the classifier.
		assertTrue(outsideAnyPeak > 1000,
				"background PSMs should be retained outside peaks, not dropped; got " + outsideAnyPeak);
	}

	@Test
	@DisplayName("R-squared filter rejects background segments rather than fitting noise")
	void backgroundDoesNotProduceSpuriousPeaks() {
		DeltaMassPeakFitter.Result result = DeltaMassPeakFitter.fit(syntheticDeltaMasses(),
				gaussianConfig());

		// Four peaks were planted. A handful of extra low-quality fits is
		// tolerable, but the R2 gate should keep this close to four rather than
		// fitting every fluctuation in the uniform background.
		assertTrue(result.peaks.size() >= 4, "all four planted peaks should be found, got " + result.peaks.size());
		assertTrue(result.peaks.size() <= 8,
				"R2 filter should suppress background segments; got " + result.peaks.size() + " peaks");
		for (DeltaMassPeakFitter.Peak peak : result.peaks) {
			assertTrue(peak.r2 >= 0.9, "kept peak below R2 threshold: " + peak);
			assertTrue(peak.sigma > 0, "kept peak with non-positive sigma: " + peak);
		}
	}

	@Test
	@DisplayName("An empty input yields no peaks and no zero peak rather than throwing")
	void emptyInputIsHandled() {
		DeltaMassPeakFitter.Result result = DeltaMassPeakFitter.fit(new double[0], new DeltaMassPeakFitter.Config());
		assertEquals(0, result.peaks.size());
		assertEquals(-1, result.zeroPeakId);
		assertTrue(Double.isNaN(result.calibrationOffset));
		assertEquals(0, result.assignment.length);
	}

	@Test
	@DisplayName("A narrow Gaussian is recovered on top of a structured local background")
	void gaussianWithBackgroundIsRecovered() {
		DeltaMassPeakFitter.Config cfg = gaussianConfig();
		cfg.windowLower = -0.05;
		cfg.windowUpper = 0.05;
		cfg.binWidth = 0.001;
		cfg.maxPeaksToFit = 20;

		List<Double> values = new java.util.ArrayList<>();
		for (int bin = 0; bin < 100; bin++) {
			double x = cfg.windowLower + (bin + 0.5) * cfg.binWidth;
			int count = (int) Math.round(20 + 80 * Math.exp(-(x * x) / (2 * 0.003 * 0.003)));
			for (int i = 0; i < count; i++) {
				values.add(x);
			}
		}
		double[] masses = new double[values.size()];
		for (int i = 0; i < masses.length; i++) {
			masses[i] = values.get(i);
		}

		DeltaMassPeakFitter.Peak peak = nearest(DeltaMassPeakFitter.fit(masses, cfg), 0.0);
		assertNotNull(peak, "the narrow excess must not be rejected because the continuum is non-zero");
		assertEquals(0.0, peak.centroid, 0.0006);
		assertEquals(0.003, peak.sigma, 0.0008);
		assertEquals(20.0, peak.background, 2.0);
		assertTrue(peak.r2 >= 0.9);

		boolean[] targets = new boolean[masses.length];
		double[] weights = new double[masses.length];
		Arrays.fill(targets, true);
		Arrays.fill(weights, 0.2);
		cfg.minPeakCount = 1;
		cfg.minPeakArea = 3;
		DeltaMassPeakFitter.Peak weighted = nearest(DeltaMassPeakFitter.fit(masses, targets, weights, cfg, null), 0.0);
		assertNotNull(weighted);
		assertEquals(weighted.psmCount * 0.2, weighted.discoverySupport, 1e-8,
				"reported discovery support must reflect spectrum-balancing weights");
	}

	@Test
	@DisplayName("Isotope links must be plausible relative to the fitted peak widths")
	void isotopeLinkUsesFittedWidthNotOnlyAbsoluteTolerance() {
		double isotope = AminoAcidProperty.isotope_deltaM;
		List<DeltaMassPeakFitter.Peak> peaks = Arrays.asList(
				new DeltaMassPeakFitter.Peak(0, 10.0, 0.0005, 100, 0.99, 100, 0, -1),
				new DeltaMassPeakFitter.Peak(1, 10.0 + isotope + 0.008, 0.0005, 50, 0.99, 50, 0, -1));

		List<DeltaMassPeakFitter.Peak> linked = DeltaMassPeakFitter.linkIsotopes(peaks,
				new DeltaMassPeakFitter.Config(), -1);

		assertEquals(0, linked.get(1).isotopeOffset,
				"an 8 mDa residual is inside the old absolute window but far outside the fitted widths");
		assertEquals(-1, linked.get(1).monoPeakId);
	}

	@Test
	@DisplayName("Deamidation remains a distinct root rather than a zero-family isotope")
	void deamidationIsNotCollapsedIntoZeroIsotopeFamily() {
		List<DeltaMassPeakFitter.Peak> peaks = Arrays.asList(
				new DeltaMassPeakFitter.Peak(0, 0.0, 0.001, 100, 0.99, 100, 0, -1),
				new DeltaMassPeakFitter.Peak(1, 0.984016, 0.001, 80, 0.99, 80, 0, -1),
				new DeltaMassPeakFitter.Peak(2, AminoAcidProperty.isotope_deltaM, 0.001, 50, 0.99, 50, 0, -1));

		List<DeltaMassPeakFitter.Peak> linked = DeltaMassPeakFitter.linkIsotopes(peaks,
				new DeltaMassPeakFitter.Config(), 0);

		assertEquals(0, linked.get(1).isotopeOffset,
				"0.984016 differs from +1 isotope by about 19.3 mDa");
		assertEquals(-1, linked.get(1).monoPeakId);
		assertEquals(1, linked.get(2).isotopeOffset);
		assertEquals(0, linked.get(2).monoPeakId);
	}

	@Test
	@DisplayName("PSM context preserves the observed mode while exposing the canonical family")
	void peakContextCarriesObservedAndNormalizedMasses() {
		DeltaMassPeakFitter.Peak root = new DeltaMassPeakFitter.Peak(0, 0.001, 0.001, 100, 0.99, 50, 0, -1);
		DeltaMassPeakFitter.Peak satellite = new DeltaMassPeakFitter.Peak(1,
				0.001 + AminoAcidProperty.isotope_deltaM, 0.001, 20, 0.99, 10, 1, 0);
		DeltaMassPeakFitter.Result result = new DeltaMassPeakFitter.Result(Arrays.asList(root, satellite), 0,
				0.001, new int[] { 1 }, new int[] { 0 }, new double[] { 0.001 }, new int[] { 1 },
				new double[] { 0.0 }, new double[] { 0.0 }, new double[] { 10.0 });
		OpenSearchPsm psm = new OpenSearchPsm();
		psm.setMassDiff(0.001 + AminoAcidProperty.isotope_deltaM);
		psm.setCalcNeutralPepMass(1000.0);

		OpenSearchFeatures.applyPeakContext(Arrays.asList(psm), result);

		assertEquals(1, psm.getObservedPeakId());
		assertEquals(0, psm.getPeakId());
		assertEquals(1, psm.getPeakIsotopeOffset());
		assertEquals(0.001, psm.getNormalizedMassDiff(), 1e-12);
		assertEquals(root.centroid, psm.getPeakCentroid(), 1e-12);
	}

	@Test
	@DisplayName("Sub-bin-width fits are rejected as unidentifiable")
	void subBinWidthFitIsRejected() {
		DeltaMassPeakFitter.Config cfg = gaussianConfig();
		cfg.windowLower = -0.05;
		cfg.windowUpper = 0.05;
		cfg.binWidth = 0.001;
		cfg.maxPeaksToFit = 20;

		List<Double> values = new java.util.ArrayList<>();
		for (int bin = 0; bin < 100; bin++) {
			double x = cfg.windowLower + (bin + 0.5) * cfg.binWidth;
			int count = 20 + (bin == 50 ? 100 : 0);
			for (int i = 0; i < count; i++) {
				values.add(x);
			}
		}
		double[] masses = values.stream().mapToDouble(Double::doubleValue).toArray();

		assertTrue(DeltaMassPeakFitter.fit(masses, cfg).peaks.isEmpty(),
				"a one-bin spike has no measurable width and must not become a modification hypothesis");
	}

	@Test
	@DisplayName("A fitted bump below the local-background SNR threshold is rejected")
	void lowSnrFitIsRejected() {
		DeltaMassPeakFitter.Config cfg = gaussianConfig();
		cfg.windowLower = -0.05;
		cfg.windowUpper = 0.05;
		cfg.binWidth = 0.001;
		cfg.maxPeaksToFit = 20;

		List<Double> values = new java.util.ArrayList<>();
		for (int bin = 0; bin < 100; bin++) {
			double x = cfg.windowLower + (bin + 0.5) * cfg.binWidth;
			int count = (int) Math.round(100 + 20 * Math.exp(-(x * x) / (2 * 0.003 * 0.003)));
			for (int i = 0; i < count; i++) {
				values.add(x);
			}
		}
		double[] masses = values.stream().mapToDouble(Double::doubleValue).toArray();

		assertTrue(DeltaMassPeakFitter.fit(masses, cfg).peaks.isEmpty(),
				"a small fluctuation above a large continuum must fail the Poisson SNR gate");
	}

	@Test
	@DisplayName("Resolution-aware consolidation keeps the stronger unresolved fit")
	void unresolvedFitsAreConsolidated() {
		DeltaMassPeakFitter.Config cfg = new DeltaMassPeakFitter.Config();
		List<DeltaMassPeakFitter.Peak> peaks = Arrays.asList(
				new DeltaMassPeakFitter.Peak(0, 10.000, 0.001, 100, 5, 0.99, 100, 0, -1),
				new DeltaMassPeakFitter.Peak(1, 10.002, 0.001, 40, 5, 0.99, 40, 0, -1),
				new DeltaMassPeakFitter.Peak(2, 10.020, 0.001, 30, 5, 0.99, 30, 0, -1));

		List<DeltaMassPeakFitter.Peak> kept = DeltaMassPeakFitter.consolidateUnresolvedPeaks(peaks, cfg);

		assertEquals(2, kept.size());
		assertEquals(10.000, kept.get(0).centroid, 1e-12,
				"the higher-area fit should represent the unresolved pair");
		assertEquals(10.020, kept.get(1).centroid, 1e-12,
				"a peak beyond the resolution limit should remain distinct");
	}

	@Test
	@DisplayName("Membership checks every eligible peak instead of leaving a false gap")
	void assignmentFallsBackToAnotherEligibleMembershipWindow() {
		List<DeltaMassPeakFitter.Peak> peaks = Arrays.asList(
				new DeltaMassPeakFitter.Peak(0, 1.002714, 0.001616, 800, 64, 0.98, 1000, 1, 2),
				new DeltaMassPeakFitter.Peak(1, 1.020338, 0.000504, 30, 24, 0.98, 100, 0, -1));
		DeltaMassPeakFitter.PeakIndex index = new DeltaMassPeakFitter.PeakIndex(peaks);
		int[] assignment = { -1 };
		double[] distance = { Double.NaN };
		double[] sigmas = { Double.NaN };

		index.assign(1.029663, new DeltaMassPeakFitter.Config(), 0, assignment, distance, sigmas);

		assertEquals(1, assignment[0],
				"the sigma-nearest broad peak is outside membership, but the narrow peak contains the PSM");
		assertEquals(1.029663 - 1.020338, distance[0], 1e-12);
	}

	/**
	 * Zero peak; carbamyl with its real +1 and +2 isotope satellites; acetyl 0.9953 Da below it, as on
	 * PXD013427; and a uniform background.
	 */
	private static double[] acetylCarbamylDeltaMasses() {
		Random rng = new Random(SEED);
		double step = AminoAcidProperty.isotope_deltaM;
		int zeroN = 40000, carbamylN = 27000, plusOneN = 4000, plusTwoN = 1200, acetylN = 19000, bgN = 30000;
		double[] out = new double[zeroN + carbamylN + plusOneN + plusTwoN + acetylN + bgN];
		int k = 0;
		for (int i = 0; i < zeroN; i++) {
			out[k++] = rng.nextGaussian() * 0.0015;
		}
		for (int i = 0; i < carbamylN; i++) {
			out[k++] = 43.005814 + rng.nextGaussian() * 0.0015;
		}
		for (int i = 0; i < plusOneN; i++) {
			out[k++] = 43.005814 + step + rng.nextGaussian() * 0.0015;
		}
		for (int i = 0; i < plusTwoN; i++) {
			out[k++] = 43.005814 + 2 * step + rng.nextGaussian() * 0.0015;
		}
		for (int i = 0; i < acetylN; i++) {
			out[k++] = 42.010565 + rng.nextGaussian() * 0.0015;
		}
		for (int i = 0; i < bgN; i++) {
			out[k++] = -150.0 + rng.nextDouble() * 650.0;
		}
		return out;
	}

	@Test
	@DisplayName("Acetyl is not linked as an isotope of carbamyl; without the guard it was")
	void acetylIsNotAnIsotopeOfCarbamyl() {
		DeltaMassPeakFitter.Result result = DeltaMassPeakFitter.fit(acetylCarbamylDeltaMasses(),
				new DeltaMassPeakFitter.Config());
		DeltaMassPeakFitter.Peak acetyl = nearest(result, 42.010565);
		DeltaMassPeakFitter.Peak carbamyl = nearest(result, 43.005814);
		assertNotNull(acetyl);
		assertNotNull(carbamyl);
		assertEquals(0, acetyl.isotopeOffset, "acetyl is its own modification");
		assertEquals(0, carbamyl.isotopeOffset, "carbamyl is its own modification");
		assertEquals(1, nearest(result, 43.005814 + AminoAcidProperty.isotope_deltaM).isotopeOffset,
				"carbamyl keeps its real +1 satellite");

		DeltaMassPeakFitter.Config unguarded = new DeltaMassPeakFitter.Config();
		unguarded.isotopeNamedMatchDa = 0;
		DeltaMassPeakFitter.Result before = DeltaMassPeakFitter.fit(acetylCarbamylDeltaMasses(), unguarded);
		assertEquals(-1, nearest(before, 42.010565).isotopeOffset, "the width rule alone links them");
	}

	@Test
	@DisplayName("The named-modification guard leaves real isotope spacings alone")
	void namedGuardKeepsRealIsotopes() {
		DeltaMassPeakFitter.Config cfg = new DeltaMassPeakFitter.Config();
		double step = AminoAcidProperty.isotope_deltaM;
		// Acetyl as the -1 of carbamyl: 7.6 mDa off the isotope step, 0.04 mDa from acetyl.
		assertTrue(DeltaMassPeakFitter.isNamedRatherThanIsotope(42.0105, Math.abs(42.0105 - (43.0058 - step)), 1,
				cfg));
		// A real +1 of the zero peak, and a +1 of oxidation landing within 1 mDa of the spacing.
		assertFalse(DeltaMassPeakFitter.isNamedRatherThanIsotope(step + 0.0003, 0.0003, 1, cfg));
		assertFalse(DeltaMassPeakFitter.isNamedRatherThanIsotope(15.9949 + step, 0.0008, 1, cfg));
		// A mass named nowhere is never blocked, however far off.
		assertFalse(DeltaMassPeakFitter.isNamedRatherThanIsotope(123.4567, 0.009, 1, cfg));
	}
}
