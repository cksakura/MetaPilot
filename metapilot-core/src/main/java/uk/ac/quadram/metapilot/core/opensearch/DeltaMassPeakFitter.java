package uk.ac.quadram.metapilot.core.opensearch;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.analysis.ParametricUnivariateFunction;
import org.apache.commons.math3.fitting.GaussianCurveFitter;
import org.apache.commons.math3.fitting.SimpleCurveFitter;
import org.apache.commons.math3.fitting.WeightedObservedPoints;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import uk.ac.quadram.metapilot.core.aminoacid.AminoAcidProperty;

/**
 * Discovers modes in the delta-mass distribution of an open search.
 * <p>
 * Used by <em>both</em> stages of the standalone open-search pipeline (see
 * {@code docs/developer/open-search-standalone-design.md}):
 * <ul>
 * <li><b>Stage 1 (genome discovery)</b> consumes only the zero family. Fitting
 * rather than applying a fixed tolerance window matters here: the fitted
 * centroid absorbs any systematic mass calibration offset, which a window
 * hard-centred on 0.0 would silently clip one side of. Exact isotope modes
 * remain observable but their PSMs are normalized onto the canonical zero
 * root before genome mapping.</li>
 * <li><b>Stage 2 (PTM validation)</b> consumes the full peak set. Peak
 * membership, distance to the nearest centroid and local density become
 * classifier <em>features</em>; nothing is discarded here. A PSM outside every
 * peak gets {@code assignment[i] == -1}, which is a feature value, not a
 * rejection.</li>
 * </ul>
 * <p>
 * Deliberately decoupled from any PSM type: it takes a {@code double[]} of
 * delta masses and returns per-index results, so it is unit-testable in
 * isolation and reusable if the engine changes.
 * <p>
 * The original MetaLab 2.0-style Gaussian-plus-background implementation is
 * retained as a compatibility path. The default uses multi-bandwidth binned
 * KDE, target-excess support, cross-file recurrence and empirical valleys,
 * followed by conservative isotope-family inference.
 */
public final class DeltaMassPeakFitter {

	private static final Logger LOGGER = LogManager.getLogger(DeltaMassPeakFitter.class);

	/** Serialises progress reporting from the parallel fitting stream. */
	private static final Object PROGRESS_LOCK = new Object();

	/** Gaussian signal on a locally constant continuum: [background, height, mean, sigma]. */
	private static final ParametricUnivariateFunction GAUSSIAN_WITH_BASELINE = new ParametricUnivariateFunction() {
		@Override
		public double value(double x, double... p) {
			double sigma = p[3];
			double z = (x - p[2]) / sigma;
			return p[0] + p[1] * Math.exp(-0.5 * z * z);
		}

		@Override
		public double[] gradient(double x, double... p) {
			double sigma = p[3];
			double delta = x - p[2];
			double exponential = Math.exp(-(delta * delta) / (2.0 * sigma * sigma));
			return new double[] { 1.0, exponential, p[1] * exponential * delta / (sigma * sigma),
					p[1] * exponential * delta * delta / (sigma * sigma * sigma) };
		}
	};

	private DeltaMassPeakFitter() {
	}

	/** Tunable parameters. Defaults are the MetaLab 2.0 values where one exists. */
	public static final class Config {
		/** Use persistent multi-bandwidth KDE instead of Gaussian curve acceptance. */
		public boolean useKde = true;
		/** Smallest KDE bandwidth (Da); wider scales are integer multiples of this. */
		public double kdeBandwidthDa = 0.001;
		/** Slash-separated KDE scale multipliers, for example {@code 1/2/4}. */
		public String kdeScaleMultipliers = "1/2/4";
		/** Minimum number of bandwidths in which a candidate mode must recur. */
		public int kdeMinPersistentScales = 2;
		/** Maximum distance from a KDE apex at which an empirical valley is sought. */
		public double kdeMaxMembershipDa = 0.050;
		/** Minimum independent files in which a mode must recur when file labels exist. */
		public int kdeMinReproducibleFiles = 2;
		/** Minimum spectrum-equivalent target support in a contributing file. */
		public double kdeMinFileSupport = 3.0;
		/** Hard ceiling for KDE modes after evidence filters; 0 disables the ceiling. */
		/**
		 * Upper bound on retained modes. Raised from 500 on 2026-09-11: the 12mix
		 * run returned exactly 500 peaks, with the weakest still carrying 161
		 * target PSMs at SNR 6 - far above every evidence gate - so the peak set
		 * was being decided by this number rather than by the data, and nothing in
		 * the task log said so. It is a runaway guard, not a tuning knob; the
		 * finder now reports when it binds.
		 */
		public int kdeMaxModes = 2000;
		/** Lower bound of the delta-mass axis (Da). Should match the search window. */
		public double windowLower = -150.0;
		/** Upper bound of the delta-mass axis (Da). Should match the search window. */
		public double windowUpper = 500.0;
		/** Histogram bin width (Da). */
		public double binWidth = 0.001;
		/** Savitzky-Golay half-width, in bins. Smoothing is for detection only. */
		public int smoothHalfWidth = 4;
		/** Savitzky-Golay polynomial order. */
		public int smoothPolyOrder = 2;
		/** Minimum raw count at a local maximum before a segment is considered. */
		public int minPeakCount = 3;

		/**
		 * Minimum spectrum-equivalent discovery support across a segment before a
		 * Gaussian is fitted to it.
		 * <p>
		 * This is the load-bearing filter, and its default is derived from the
		 * test the peak will later face rather than chosen by taste. With 1:1
		 * target:decoy, a peak of {@code n} PSMs has best-case binomial
		 * p = 0.5^n — that is, the p-value it gets if <em>every single PSM in it
		 * is a target</em>. At n=3 that is 0.125 and at n=7 it is 0.0078, so a
		 * peak below roughly this size cannot reach significance however clean it
		 * is, and cannot be reported.
		 * <p>
		 * Fitting such segments is therefore not merely wasted work, it is
		 * harmful: every one raises the multiple-testing burden m that the
		 * Benjamini-Hochberg step applies to the real peaks. At m=24,491 (an
		 * observed run) a genuine peak needs >=22 all-target PSMs to stay
		 * callable; at m=200 it needs 15. Fitting noise suppresses real
		 * modifications.
		 * <p>
		 * {@link #minPeakCount} gates the apex bin only. Over a 650,000-bin
		 * histogram holding ~1M PSMs — averaging ~1.5 per bin — random
		 * fluctuation clears an apex threshold of 3 constantly, which is how a
		 * run reached 24,491 "peaks". Area is the honest measure of a peak.
		 */
		public int minPeakArea = 15;

		/**
		 * Hard cap on Gaussian fits per run; segments are ranked by area and the
		 * largest kept. Bounds the worst case no matter how pathological the
		 * histogram, so the step can never again run for an unknown number of
		 * hours. 0 or less means no cap.
		 */
		public int maxPeaksToFit = 500;

		/**
		 * Iteration ceiling for a single Gaussian fit.
		 * <p>
		 * A segment that is genuinely peak-shaped converges quickly; one that is
		 * not runs to the optimiser's default limit and is discarded anyway. An
		 * observed run averaged ~147 ms per fit, which is the signature of
		 * Levenberg-Marquardt grinding to that limit on non-peaks.
		 */
		public int maxFitIterations = 200;
		/** Minimum R² for a fitted Gaussian to be kept. */
		public double minR2 = 0.9;
		/**
		 * Minimum number of raw histogram bins in a fitted segment. A Gaussian plus
		 * background has four parameters; seven bins leave at least three residual
		 * degrees of freedom and prevent four-point fits from achieving a meaningless
		 * R² of 1.0.
		 */
		public int minFitBins = 7;
		/**
		 * Minimum fitted sigma in histogram-bin units. Widths below half a bin are not
		 * identifiable from a binned histogram and should not become modification
		 * hypotheses.
		 */
		public double minSigmaBins = 0.5;
		/**
		 * Minimum peak height relative to Poisson fluctuation in the fitted local
		 * background: {@code amplitude / sqrt(background + 1)}.
		 */
		public double minPeakSnr = 5.0;
		/** Absolute resolution floor, in histogram bins, for distinct fitted peaks. */
		public double minPeakSeparationBins = 2.0;
		/**
		 * Width-scaled resolution floor for distinct fitted peaks, in combined fitted
		 * sigmas. Unresolved alternatives are reduced to the stronger fit.
		 */
		public double minPeakSeparationSigmas = 3.0;

		/**
		 * Reject a mode whose width exceeds this multiple of the median mode
		 * width; 0 disables it. A mode much broader than the instrument's own
		 * precursor accuracy is an unresolved blend rather than a peak. Off by
		 * default so the resolution gate above can be measured on its own first.
		 */
		public double kdeMaxSigmaMultiple = 0.0;

		/**
		 * Minimum topographic prominence, as a fraction of a mode's own height; 0
		 * disables it. A local maximum sitting on a broad elevated region is a
		 * ripple, not a peak, and the absolute support test cannot tell them apart
		 * because a large background passes it at every wiggle. Off by default: the
		 * finder reports the prominence distribution of the modes it retained, so
		 * the threshold is chosen from the data rather than guessed.
		 */
		public double kdeMinProminenceRatio = 0.0;

		/**
		 * Per-row mask restricting which candidates build the detection histogram,
		 * applied to targets AND decoys alike - unlike {@code useForPeakDetection},
		 * which gates targets only. Null means every row counts.
		 * <p>
		 * The runner sets it to rank-1 membership: before rescoring, the engine's own
		 * ranking is the only ordering available, so admitting ranks 2-5 here would
		 * contradict the evidence the step relies on. Measured on 12mix, 57.4% of
		 * rank-1 candidates fall inside a fitted peak against 16.9% of ranks 2-5, and
		 * under per-spectrum weighting those lower ranks carried roughly 80% of each
		 * spectrum's vote into diffuse background.
		 * <p>
		 * Kept here rather than added to the {@code fit} signature so the existing
		 * overloads keep their semantics exactly: a caller passing
		 * {@code useForPeakDetection} meaning "targets" behaves as it always did.
		 * Overloading that parameter instead would have excluded every decoy from the
		 * density - which DeltaMassKdePeakFinderTest catches.
		 */
		public boolean[] detectionRows = null;

		/** Half-width of the window searched for the zero peak (Da). */
		public double zeroWindow = 0.2;
		/** Membership half-width, in fitted sigmas. */
		public double membershipSigmas = 3.0;

		/**
		 * Absolute lower bound on the membership half-width (Da). A PSM belongs to
		 * a peak when it is within {@code max(membershipSigmas * sigma, this)}.
		 *
		 * <h2>Why a floor is needed</h2>
		 * Measured delta-mass peaks are not Gaussian. They have a sharp core a few
		 * mDa wide and heavy, long tails: on real data the zero peak's core is
		 * ~3 mDa FWHM while appreciable signal continues past 14 mDa. A Gaussian
		 * fitted to that shape recovers the <em>core</em> faithfully — measured
		 * sigma 0.9-1.3 mDa, which is correct — and then a pure
		 * {@code 3 * sigma} rule cuts membership at +/-3 mDa and discards the
		 * tails. Those tails are the same peptides measured with slightly worse
		 * mass accuracy, so they are real members being lost.
		 * <p>
		 * The fix is not to widen the fit, which would mis-describe the core, but
		 * to recognise that the relevant scale for membership is instrument mass
		 * accuracy rather than fit geometry.
		 * <p>
		 * Empirically on a three-sample run: a 10 mDa floor reabsorbed 16,160
		 * previously unassigned PSMs at a target:decoy ratio of 204:1, against
		 * 299:1 inside fitted peaks and 100:1 for the unassigned pool as a whole
		 * — i.e. the reabsorbed PSMs look like peak members, not background.
		 * Beyond ~15 mDa the ratio decays toward background, which is where the
		 * window stops collecting tails and starts collecting the broad
		 * inter-peak continuum.
		 */
		public double membershipFloorDa = 0.010;
		/** Tolerance for matching a peak to an isotope of another peak (Da). */
		public double isotopeTolerance = 0.01;
		/** Minimum isotope-link tolerance after scaling by fitted peak widths (Da). */
		public double isotopeToleranceFloorDa = 0.001;
		/** Number of combined fitted sigmas allowed for an isotope link. */
		public double isotopeSigmaMultiplier = 4.0;
		/** Half-width of the window used for the local-density feature (Da). */
		public double densityWindow = 0.05;
		/** Highest isotope offset considered (+1, +2, ...). */
		public int maxIsotopeOffset = 3;
		/**
		 * Highest isotope offset folded onto the zero family in stage 1.
		 * <p>
		 * Deliberately tighter than {@link #maxIsotopeOffset}. Stage 2 builds
		 * isotope families across the whole delta-mass axis and has independent
		 * evidence for each link (consecutive satellites, monotone amplitude,
		 * spacing residual), so a generous bound there costs nothing. Stage 1 has
		 * no such evidence: a PSM is folded onto zero on mass proximity alone, so
		 * every extra offset is an extra chance to absorb an unrelated shift.
		 * <p>
		 * Measured on the 12mix marker stage against
		 * {@code 12mix_genomes_mapping.tsv}: the monoisotopic window is 93.9%
		 * correct, +1 is 94.1% and +2 is 91.8% — already declining. +3 was never
		 * measured; it was inherited from {@code maxIsotopeOffset} rather than
		 * chosen, and contributes ~2% of the stage-1 set. Raise this only with a
		 * correct-ratio measurement for the offsets being admitted.
		 */
		public int round1MaxIsotopeOffset = 2;
		/** Minimum satellites required to collapse a non-zero modification family. */
		public int minIsotopeFamilySatellites = 2;
		/** Maximum permitted KDE-amplitude increase between successive satellites. */
		public double isotopeMaxAbundanceIncrease = 1.25;

		/**
		 * Effective stage-1 isotope bound. Clamped to {@code maxIsotopeOffset} so
		 * stage 1 can never fold an offset stage 2 would refuse to link, and to
		 * zero so a negative setting disables folding rather than inverting the
		 * loop.
		 */
		int round1IsotopeBound() {
			return Math.max(0, Math.min(maxIsotopeOffset, round1MaxIsotopeOffset));
		}

		double[] kdeBandwidthsDa() {
			String[] values = kdeScaleMultipliers == null ? new String[0] : kdeScaleMultipliers.split("/");
			List<Double> bandwidths = new ArrayList<>();
			for (String value : values) {
				try {
					double multiplier = Double.parseDouble(value.trim());
					if (multiplier > 0.0 && Double.isFinite(multiplier)) {
						bandwidths.add(kdeBandwidthDa * multiplier);
					}
				} catch (NumberFormatException ignored) {
					// Invalid entries are ignored; the base bandwidth below is mandatory.
				}
			}
			if (bandwidths.isEmpty()) {
				bandwidths.add(kdeBandwidthDa);
			}
			bandwidths.sort(null);
			double[] out = new double[bandwidths.size()];
			for (int i = 0; i < out.length; i++) {
				out[i] = bandwidths.get(i);
			}
			return out;
		}
	}

	/** One fitted delta-mass peak. */
	public static final class Peak {
		public final int id;
		/** Fitted Gaussian mean (Da). For the zero peak this is the calibration offset. */
		public final double centroid;
		/** Fitted Gaussian standard deviation (Da). */
		public final double sigma;
		/** Fitted Gaussian height. */
		public final double amplitude;
		/** Constant local background fitted underneath the Gaussian. */
		public final double background;
		/** Goodness of fit against the raw binned counts. */
		public final double r2;
		/** Number of input delta masses assigned to this peak. */
		public final int psmCount;
		/** Spectrum-equivalent target support used by the discovery histogram. */
		public final double discoverySupport;
		/**
		 * 0 when this peak is monoisotopic (or unlinked), k &gt; 0 when it sits
		 * k isotope spacings above {@link #monoPeakId}. A negative value marks a
		 * peak below a protected family root (most importantly the zero peak); it
		 * is isotope-related selection error, not a second monoisotopic root.
		 */
		public final int isotopeOffset;
		/** Id of the monoisotopic peak this is an isotope of, or -1. */
		public final int monoPeakId;
		/** Empirical lower KDE-basin boundary; NaN for legacy Gaussian fits. */
		public final double lowerBoundary;
		/** Empirical upper KDE-basin boundary; NaN for legacy Gaussian fits. */
		public final double upperBoundary;
		/** Number of KDE bandwidths supporting this mode. */
		public final int persistentScales;
		/** Number of independent files supporting this mode. */
		public final int reproducibleFiles;
		/** Discovery model name, currently KDE or GAUSSIAN. */
		public final String discoveryMethod;

		Peak(int id, double centroid, double sigma, double amplitude, double r2, int psmCount, int isotopeOffset,
				int monoPeakId) {
			this(id, centroid, sigma, amplitude, 0.0, r2, psmCount, psmCount, isotopeOffset, monoPeakId);
		}

		Peak(int id, double centroid, double sigma, double amplitude, double background, double r2, int psmCount,
				int isotopeOffset, int monoPeakId) {
			this(id, centroid, sigma, amplitude, background, r2, psmCount, psmCount, isotopeOffset, monoPeakId);
		}

		Peak(int id, double centroid, double sigma, double amplitude, double background, double r2, int psmCount,
				double discoverySupport, int isotopeOffset, int monoPeakId) {
			this(id, centroid, sigma, amplitude, background, r2, psmCount, discoverySupport, isotopeOffset,
					monoPeakId, Double.NaN, Double.NaN, 0, 0, "GAUSSIAN");
		}

		Peak(int id, double centroid, double sigma, double amplitude, double background, double r2, int psmCount,
				double discoverySupport, int isotopeOffset, int monoPeakId, double lowerBoundary,
				double upperBoundary, int persistentScales, int reproducibleFiles, String discoveryMethod) {
			this.id = id;
			this.centroid = centroid;
			this.sigma = sigma;
			this.amplitude = amplitude;
			this.background = background;
			this.r2 = r2;
			this.psmCount = psmCount;
			this.discoverySupport = discoverySupport;
			this.isotopeOffset = isotopeOffset;
			this.monoPeakId = monoPeakId;
			this.lowerBoundary = lowerBoundary;
			this.upperBoundary = upperBoundary;
			this.persistentScales = persistentScales;
			this.reproducibleFiles = reproducibleFiles;
			this.discoveryMethod = discoveryMethod;
		}

		public int familyPeakId() {
			return isotopeOffset == 0 || monoPeakId < 0 ? id : monoPeakId;
		}

		public boolean isIsotope() {
			return isotopeOffset != 0;
		}

		/** Peak height measured in local-background Poisson standard deviations. */
		public double signalToNoise() {
			return amplitude / Math.sqrt(background + 1.0);
		}

		/** Peak height relative to the fitted local background. */
		public double signalToBackground() {
			return background > 0.0 ? amplitude / background : Double.POSITIVE_INFINITY;
		}

		@Override
		public String toString() {
			return String.format("Peak[%d] centroid=%.5f sigma=%.5f amp=%.1f bg=%.1f r2=%.3f n=%d support=%.1f%s",
					id, centroid, sigma, amplitude, background, r2, psmCount, discoverySupport,
					isotopeOffset != 0 ? (" isotope" + (isotopeOffset > 0 ? "+" : "") + isotopeOffset + " of "
							+ monoPeakId) : "");
		}
	}

	/** Fit output, indexed parallel to the input delta-mass array. */
	public static final class Result {
		public final List<Peak> peaks;
		/** Id of the peak containing zero, or -1 if none was fitted. */
		public final int zeroPeakId;
		/**
		 * Systematic mass offset (Da), i.e. the zero peak's centroid. NaN when
		 * no zero peak was fitted.
		 */
		public final double calibrationOffset;
		/** Per input index: assigned peak id, or -1 when inside no peak. */
		public final int[] assignment;
		/** Per input index: observed KDE mode before isotope-family collapse. */
		public final int[] observedAssignment;
		/** Per input index: canonical root mode after isotope-family collapse. */
		public final int[] familyAssignment;
		/** Raw mass difference normalized onto the canonical isotope-family root. */
		public final double[] normalizedMassDiff;
		/** Per-PSM isotope offset used during family normalization. */
		public final int[] isotopeOffset;
		/** Per input index: signed distance to the nearest peak centroid (Da). */
		public final double[] distanceToNearestCentroid;
		/** Per input index: |distance| to the nearest centroid, in that peak's sigmas. */
		public final double[] distanceInSigmas;
		/**
		 * Per input index: spectrum-weighted mass of delta masses within
		 * {@code densityWindow}.
		 * <p>
		 * Weighted with the same {@code detectionWeights} the peak fit uses, so a
		 * spectrum contributes 1 however many candidates it reported and however
		 * many engines reported it. It was a raw count until 2026-09-12, which made
		 * it scale with the size of the candidate pool: pooling a second engine
		 * doubled every value (40,269 -> 81,684 mean on 12mix) while the fitted
		 * peaks, being weighted, did not move. The logistic regression then
		 * saturated for the low-density population and scored all 4,383,711
		 * unexplained candidates at exactly 0, destroying every distinction inside
		 * that group.
		 */
		public final double[] localDensity;

		Result(List<Peak> peaks, int zeroPeakId, double calibrationOffset, int[] observedAssignment,
				int[] familyAssignment, double[] normalizedMassDiff, int[] isotopeOffset,
				double[] distanceToNearestCentroid, double[] distanceInSigmas, double[] localDensity) {
			this.peaks = peaks;
			this.zeroPeakId = zeroPeakId;
			this.calibrationOffset = calibrationOffset;
			this.assignment = observedAssignment;
			this.observedAssignment = observedAssignment;
			this.familyAssignment = familyAssignment;
			this.normalizedMassDiff = normalizedMassDiff;
			this.isotopeOffset = isotopeOffset;
			this.distanceToNearestCentroid = distanceToNearestCentroid;
			this.distanceInSigmas = distanceInSigmas;
			this.localDensity = localDensity;
		}

		public Peak peakOf(int index) {
			int id = observedAssignment[index];
			return id < 0 ? null : peaks.get(id);
		}

		public Peak familyPeakOf(int index) {
			int id = familyAssignment[index];
			return id < 0 ? null : peaks.get(id);
		}

		/**
		 * Indices belonging to the zero peak - stage 1's "unmodified" set.
		 * Observed isotope modes remain separate in {@link #observedAssignment},
		 * while this family view includes members normalized onto the zero root.
		 */
		public int[] zeroPeakIndices() {
			if (zeroPeakId < 0) {
				return new int[0];
			}
			int n = 0;
			for (int a : familyAssignment) {
				if (a == zeroPeakId) {
					n++;
				}
			}
			int[] out = new int[n];
			int k = 0;
			for (int i = 0; i < familyAssignment.length; i++) {
				if (familyAssignment[i] == zeroPeakId) {
					out[k++] = i;
				}
			}
			return out;
		}
	}

	public static Result fit(double[] deltaMasses, Config cfg) {
		return fit(deltaMasses, null, cfg, null);
	}

	public static Result fit(double[] deltaMasses, boolean[] useForPeakDetection, Config cfg) {
		return fit(deltaMasses, useForPeakDetection, cfg, null);
	}

	/**
	 * Fits peaks from a weighted discovery histogram while still assigning every
	 * input PSM. The weights affect only peak discovery and curve fitting; they do
	 * not alter downstream target/decoy counts or per-PSM membership.
	 */
	public static Result fit(double[] deltaMasses, boolean[] useForPeakDetection, double[] detectionWeights,
			Config cfg, java.util.function.Consumer<String> progress) {
		return fitInternal(deltaMasses, useForPeakDetection, detectionWeights, null, null, cfg, progress);
	}

	/**
	 * Full KDE discovery entry point. Target labels create a target-excess density;
	 * file groups provide independent reproducibility evidence for each mode.
	 */
	public static Result fit(double[] deltaMasses, boolean[] useForPeakDetection, double[] detectionWeights,
			boolean[] isTarget, String[] fileGroups, Config cfg,
			java.util.function.Consumer<String> progress) {
		return fitInternal(deltaMasses, useForPeakDetection, detectionWeights, isTarget, fileGroups, cfg, progress);
	}

	/**
	 * Discovers <b>only</b> the zero family, over
	 * {@code +/-cfg.zeroWindow} Da after isotope folding.
	 * <p>
	 * This is what stage 1 needs and all it needs. Stage 1 selects genomes; it
	 * does not identify modifications, so the delta-mass landscape above and
	 * below zero is irrelevant to it — the only question is where the
	 * unmodified peak starts and stops.
	 * <p>
	 * The default path uses the same non-parametric, target-excess KDE principle
	 * as stage 2, but folds the finite set of expected isotope offsets onto the
	 * zero window first. It therefore examines only a few constrained mass
	 * regions rather than discovering modes over the full [-150, 500] Da axis.
	 * The legacy Gaussian implementation remains available when
	 * {@code cfg.useKde == false} for controlled comparisons.
	 * <p>
	 * The returned {@link Result} carries exactly one peak, so
	 * {@code zeroPeakId} is 0 whenever discovery succeeds. Isotope-offset PSMs
	 * are normalized onto that root and retain their offset as audit metadata;
	 * they are not emitted as separate modification peaks.
	 */
	public static Result fitZeroPeakOnly(double[] deltaMasses, boolean[] useForPeakDetection, Config cfg,
			java.util.function.Consumer<String> progress) {
		if (deltaMasses == null) {
			throw new IllegalArgumentException("deltaMasses is null");
		}
		if (cfg == null) {
			cfg = new Config();
		}
		if (useForPeakDetection != null && useForPeakDetection.length != deltaMasses.length) {
			throw new IllegalArgumentException("useForPeakDetection length does not match deltaMasses");
		}
		if (cfg.useKde) {
			return discoverZeroFamilyKde(deltaMasses, useForPeakDetection, cfg, progress);
		}

		Config windowed = copyWithWindow(cfg, -cfg.zeroWindow, cfg.zeroWindow);
		int binCount = (int) Math.ceil((windowed.windowUpper - windowed.windowLower) / windowed.binWidth);

		double[] histogram = new double[binCount];
		int inWindow = 0;
		for (int i = 0; i < deltaMasses.length; i++) {
			if (useForPeakDetection != null && !useForPeakDetection[i]) {
				continue;
			}
			int bin = binOf(deltaMasses[i], windowed);
			if (bin >= 0 && bin < binCount) {
				histogram[bin]++;
				inWindow++;
			}
		}

		report(progress, String.format(
				"Zero-peak fit: %d of %d delta masses fall within +/-%.3f Da, binned into %d bins of %.4f Da",
				inWindow, deltaMasses.length, cfg.zeroWindow, binCount, windowed.binWidth));

		double[] smoothed = smooth(histogram, windowed.smoothHalfWidth, windowed.smoothPolyOrder);

		// The zero peak is the tallest thing in this window by a wide margin, so
		// take the global maximum rather than enumerating local maxima.
		int apex = -1;
		for (int bin = 0; bin < binCount; bin++) {
			if (apex < 0 || smoothed[bin] > smoothed[apex]) {
				apex = bin;
			}
		}
		if (apex < 0 || histogram[apex] < windowed.minPeakCount) {
			report(progress, "Zero-peak fit: no apex above minPeakCount=" + windowed.minPeakCount
					+ " within the zero window");
			return emptyResult(deltaMasses.length);
		}

		// Walk out to where the smoothed curve stops descending: the start and
		// stop of the curve, which is the whole of what has to be determined.
		int left = apex;
		while (left > 0 && smoothed[left - 1] < smoothed[left] && smoothed[left - 1] > 0) {
			left--;
		}
		int right = apex;
		while (right < binCount - 1 && smoothed[right + 1] < smoothed[right] && smoothed[right + 1] > 0) {
			right++;
		}

		Peak fittedPeak = fitSegment(histogram, left, right, 0, windowed);
		if (fittedPeak == null) {
			report(progress, String.format(
					"Zero-peak fit: Gaussian fit failed or fell below R2 >= %.2f over bins [%.4f, %.4f] Da",
					windowed.minR2, binCentre(left, windowed), binCentre(right, windowed)));
			return emptyResult(deltaMasses.length);
		}

		// Same floored half-width as the full fit: the zero peak has the same
		// sharp-core / heavy-tail shape as any other, so a pure 3-sigma window
		// would drop well-measured unmodified peptides out of stage 1's
		// genome-discovery set.
		double halfWidth = Math.max(windowed.membershipSigmas * fittedPeak.sigma, windowed.membershipFloorDa);
		double lo = fittedPeak.centroid - halfWidth;
		double hi = fittedPeak.centroid + halfWidth;

		int[] assignment = new int[deltaMasses.length];
		double[] distance = new double[deltaMasses.length];
		double[] sigmas = new double[deltaMasses.length];
		Arrays.fill(assignment, -1);
		Arrays.fill(distance, Double.NaN);
		Arrays.fill(sigmas, Double.NaN);

		// One interval test per PSM. No scan over a peak list, because there is
		// only one peak.
		int assigned = 0;
		for (int i = 0; i < deltaMasses.length; i++) {
			double d = deltaMasses[i] - fittedPeak.centroid;
			distance[i] = d;
			sigmas[i] = Math.abs(d) / fittedPeak.sigma;
			if (deltaMasses[i] >= lo && deltaMasses[i] <= hi) {
				assignment[i] = 0;
				assigned++;
			}
		}

		int[] familyAssignment = Arrays.copyOf(assignment, assignment.length);
		double[] normalizedMassDiff = Arrays.copyOf(deltaMasses, deltaMasses.length);
		int[] isotopeOffset = new int[deltaMasses.length];
		int bound = windowed.round1IsotopeBound();
		for (int i = 0; i < deltaMasses.length; i++) {
			if (familyAssignment[i] == 0) {
				continue;
			}
			for (int offset = -bound; offset <= bound; offset++) {
				if (offset == 0) {
					continue;
				}
				double normalized = deltaMasses[i] - offset * AminoAcidProperty.isotope_deltaM;
				if (normalized >= lo && normalized <= hi) {
					familyAssignment[i] = 0;
					normalizedMassDiff[i] = normalized;
					isotopeOffset[i] = offset;
					distance[i] = normalized - fittedPeak.centroid;
					sigmas[i] = Math.abs(distance[i]) / fittedPeak.sigma;
					assigned++;
					break;
				}
			}
		}

		List<Peak> peaks = new ArrayList<>();
		peaks.add(new Peak(0, fittedPeak.centroid, fittedPeak.sigma, fittedPeak.amplitude, fittedPeak.background,
				fittedPeak.r2, assigned, 0, -1));

		report(progress, String.format(
				"Zero-peak fit: centroid %.5f Da (calibration offset), sigma %.5f Da, R2 %.3f, "
						+ "curve spans [%.5f, %.5f] Da, %d PSMs within +/-%.1f sigma",
				fittedPeak.centroid, fittedPeak.sigma, fittedPeak.r2, lo, hi, assigned, windowed.membershipSigmas));

		// localDensity is a stage-2 classifier feature; stage 1 does not use it,
		// and computing it would mean walking the full axis again.
		return new Result(peaks, 0, fittedPeak.centroid, assignment, familyAssignment, normalizedMassDiff,
				isotopeOffset, distance, sigmas,
				new double[deltaMasses.length]);
	}

	/**
	 * Constrained KDE discovery for stage 1. Each PSM that lies near one of the
	 * finite expected isotope offsets is folded onto the monoisotopic coordinate
	 * before density estimation. This pools evidence for the same chemical state
	 * without treating isotope-selection artefacts as separate modifications.
	 */
	private static Result discoverZeroFamilyKde(double[] deltaMasses, boolean[] isTarget, Config cfg,
			java.util.function.Consumer<String> progress) {
		Config windowed = copyWithWindow(cfg, -cfg.zeroWindow, cfg.zeroWindow);
		double[] foldedMasses = Arrays.copyOf(deltaMasses, deltaMasses.length);
		int[] candidateOffsets = new int[deltaMasses.length];
		int foldedCount = 0;
		int isotopeBound = windowed.round1IsotopeBound();
		for (int i = 0; i < deltaMasses.length; i++) {
			double mass = deltaMasses[i];
			double bestResidual = mass;
			int bestOffset = 0;
			for (int offset = -isotopeBound; offset <= isotopeBound; offset++) {
				double residual = mass - offset * AminoAcidProperty.isotope_deltaM;
				if (Math.abs(residual) < Math.abs(bestResidual)) {
					bestResidual = residual;
					bestOffset = offset;
				}
			}
			if (Math.abs(bestResidual) <= windowed.zeroWindow) {
				foldedMasses[i] = bestResidual;
				candidateOffsets[i] = bestOffset;
				if (bestOffset != 0) {
					foldedCount++;
				}
			}
		}

		// In stage 1 the supplied mask is exactly the target/decoy label. Targets
		// contribute signal and decoys estimate local background, matching the
		// target-excess KDE used by full stage-2 discovery.
		List<Peak> modes = DeltaMassKdePeakFinder.discover(foldedMasses, null, null, isTarget, null, windowed);
		int discoveredZeroId = findZeroPeak(modes, windowed);
		if (discoveredZeroId < 0) {
			report(progress, String.format(
					"Zero-family KDE: no persistent target-excess mode found within +/-%.3f Da",
					windowed.zeroWindow));
			return emptyResult(deltaMasses.length);
		}

		Peak discovered = modes.get(discoveredZeroId);
		double lo = Double.isFinite(discovered.lowerBoundary) ? discovered.lowerBoundary
				: discovered.centroid - windowed.membershipFloorDa;
		double hi = Double.isFinite(discovered.upperBoundary) ? discovered.upperBoundary
				: discovered.centroid + windowed.membershipFloorDa;

		int[] observedAssignment = new int[deltaMasses.length];
		int[] familyAssignment = new int[deltaMasses.length];
		double[] normalizedMassDiff = Arrays.copyOf(deltaMasses, deltaMasses.length);
		int[] isotopeOffset = new int[deltaMasses.length];
		double[] distance = new double[deltaMasses.length];
		double[] sigmas = new double[deltaMasses.length];
		Arrays.fill(observedAssignment, -1);
		Arrays.fill(familyAssignment, -1);
		Arrays.fill(distance, Double.NaN);
		Arrays.fill(sigmas, Double.NaN);

		int assigned = 0;
		int isotopeMembers = 0;
		for (int i = 0; i < deltaMasses.length; i++) {
			double normalized = foldedMasses[i];
			if (normalized < lo || normalized > hi) {
				continue;
			}
			int offset = candidateOffsets[i];
			familyAssignment[i] = 0;
			normalizedMassDiff[i] = normalized;
			isotopeOffset[i] = offset;
			distance[i] = normalized - discovered.centroid;
			sigmas[i] = Math.abs(distance[i]) / discovered.sigma;
			// Only the monoisotopic observations belong directly to the root mode.
			// Satellites remain auditable through isotopeOffset, but are not peaks.
			if (offset == 0) {
				observedAssignment[i] = 0;
			} else {
				isotopeMembers++;
			}
			assigned++;
		}

		Peak zeroPeak = new Peak(0, discovered.centroid, discovered.sigma, discovered.amplitude,
				discovered.background, Double.NaN, assigned, discovered.discoverySupport, 0, -1, lo, hi,
				discovered.persistentScales, discovered.reproducibleFiles, "KDE");
		report(progress, String.format(
				"Zero-family KDE: centroid %.5f Da (calibration offset), empirical basin [%.5f, %.5f] Da, "
						+ "%d PSMs assigned (%d isotope-normalized from %d candidates)",
				zeroPeak.centroid, lo, hi, assigned, isotopeMembers, foldedCount));

		return new Result(List.of(zeroPeak), 0, zeroPeak.centroid, observedAssignment, familyAssignment,
				normalizedMassDiff, isotopeOffset, distance, sigmas, new double[deltaMasses.length]);
	}

	/**
	 * Counts delta masses lying within {@code isotopeTolerance} of +1, +2, ...
	 * {@code maxIsotopeOffset} carbon-13 spacings.
	 * <p>
	 * Replaces the old "did we fit any isotope peaks?" check, which required
	 * fitting the whole axis. A run with {@code isotope_error != 0} folds these
	 * precursors into the zero peak, so a count of zero here is the signal that
	 * the exclusion stage 1 relies on is not actually happening.
	 */
	public static int countNearIsotopeOffsets(double[] deltaMasses, Config cfg) {
		if (cfg == null) {
			cfg = new Config();
		}
		int count = 0;
		for (double delta : deltaMasses) {
			for (int k = 1; k <= cfg.maxIsotopeOffset; k++) {
				if (Math.abs(delta - k * AminoAcidProperty.isotope_deltaM) <= cfg.isotopeTolerance) {
					count++;
					break;
				}
			}
		}
		return count;
	}

	private static Config copyWithWindow(Config source, double lower, double upper) {
		Config copy = new Config();
		copy.useKde = source.useKde;
		copy.kdeBandwidthDa = source.kdeBandwidthDa;
		copy.kdeScaleMultipliers = source.kdeScaleMultipliers;
		copy.kdeMinPersistentScales = source.kdeMinPersistentScales;
		copy.kdeMaxMembershipDa = source.kdeMaxMembershipDa;
		copy.kdeMinReproducibleFiles = source.kdeMinReproducibleFiles;
		copy.kdeMinFileSupport = source.kdeMinFileSupport;
		copy.detectionRows = source.detectionRows;
		copy.kdeMaxModes = source.kdeMaxModes;
		copy.kdeMaxSigmaMultiple = source.kdeMaxSigmaMultiple;
		copy.kdeMinProminenceRatio = source.kdeMinProminenceRatio;
		copy.windowLower = lower;
		copy.windowUpper = upper;
		copy.binWidth = source.binWidth;
		copy.smoothHalfWidth = source.smoothHalfWidth;
		copy.smoothPolyOrder = source.smoothPolyOrder;
		copy.minPeakCount = source.minPeakCount;
		copy.minPeakArea = source.minPeakArea;
		copy.maxPeaksToFit = source.maxPeaksToFit;
		// Must be carried: fitZeroPeakOnly calls fitSegment with this copy, and a
		// zero iteration ceiling would fail the one fit stage 1 depends on.
		copy.maxFitIterations = source.maxFitIterations;
		copy.minR2 = source.minR2;
		copy.minFitBins = source.minFitBins;
		copy.minSigmaBins = source.minSigmaBins;
		copy.minPeakSnr = source.minPeakSnr;
		copy.minPeakSeparationBins = source.minPeakSeparationBins;
		copy.minPeakSeparationSigmas = source.minPeakSeparationSigmas;
		copy.zeroWindow = source.zeroWindow;
		copy.membershipSigmas = source.membershipSigmas;
		copy.membershipFloorDa = source.membershipFloorDa;
		copy.isotopeTolerance = source.isotopeTolerance;
		copy.isotopeToleranceFloorDa = source.isotopeToleranceFloorDa;
		copy.isotopeSigmaMultiplier = source.isotopeSigmaMultiplier;
		copy.densityWindow = source.densityWindow;
		copy.maxIsotopeOffset = source.maxIsotopeOffset;
		// Must be carried: both stage-1 fold loops read it off this copy.
		copy.round1MaxIsotopeOffset = source.round1MaxIsotopeOffset;
		copy.minIsotopeFamilySatellites = source.minIsotopeFamilySatellites;
		copy.isotopeMaxAbundanceIncrease = source.isotopeMaxAbundanceIncrease;
		return copy;
	}

	private static Result emptyResult(int n) {
		int[] assignment = new int[n];
		double[] distance = new double[n];
		double[] sigmas = new double[n];
		Arrays.fill(assignment, -1);
		Arrays.fill(distance, Double.NaN);
		Arrays.fill(sigmas, Double.NaN);
		double[] normalized = new double[n];
		Arrays.fill(normalized, Double.NaN);
		return new Result(new ArrayList<>(), -1, Double.NaN, assignment, Arrays.copyOf(assignment, n), normalized,
				new int[n], distance, sigmas, new double[n]);
	}

	private static void report(java.util.function.Consumer<String> progress, String message) {
		LOGGER.info(message);
		if (progress != null) {
			progress.accept(message);
		}
	}

	/**
	 * @param deltaMasses          delta mass per PSM (Da)
	 * @param useForPeakDetection  optional mask selecting which entries build the
	 *                             histogram peaks are detected from; {@code null}
	 *                             uses all. Passing targets-only produces cleaner
	 *                             peaks (decoys contribute a broadly uniform
	 *                             background) at the cost of making downstream
	 *                             peak-level target-enrichment statistics
	 *                             optimistic - see the design doc's calibration
	 *                             caveat. Assignment always covers every entry.
	 */
	public static Result fit(double[] deltaMasses, boolean[] useForPeakDetection, Config cfg,
			java.util.function.Consumer<String> progress) {
		return fitInternal(deltaMasses, useForPeakDetection, null, null, null, cfg, progress);
	}

	private static Result fitInternal(double[] deltaMasses, boolean[] useForPeakDetection, double[] detectionWeights,
			boolean[] isTarget, String[] fileGroups, Config cfg,
			java.util.function.Consumer<String> progress) {
		if (deltaMasses == null) {
			throw new IllegalArgumentException("deltaMasses is null");
		}
		if (cfg == null) {
			cfg = new Config();
		}
		if (useForPeakDetection != null && useForPeakDetection.length != deltaMasses.length) {
			throw new IllegalArgumentException("useForPeakDetection length does not match deltaMasses");
		}
		if (detectionWeights != null && detectionWeights.length != deltaMasses.length) {
			throw new IllegalArgumentException("detectionWeights length does not match deltaMasses");
		}
		if (isTarget != null && isTarget.length != deltaMasses.length) {
			throw new IllegalArgumentException("isTarget length does not match deltaMasses");
		}
		if (fileGroups != null && fileGroups.length != deltaMasses.length) {
			throw new IllegalArgumentException("fileGroups length does not match deltaMasses");
		}
		if (cfg.useKde) {
			return fitKdeInternal(deltaMasses, useForPeakDetection, detectionWeights, isTarget, fileGroups, cfg,
					progress);
		}

		int binCount = (int) Math.ceil((cfg.windowUpper - cfg.windowLower) / cfg.binWidth);
		if (binCount <= 0) {
			throw new IllegalArgumentException("Invalid delta-mass window: " + cfg.windowLower + ".." + cfg.windowUpper);
		}

		double[] histogram = new double[binCount];
		double effectiveWeight = 0.0;
		for (int i = 0; i < deltaMasses.length; i++) {
			if (useForPeakDetection != null && !useForPeakDetection[i]) {
				continue;
			}
			if (cfg.detectionRows != null && !cfg.detectionRows[i]) {
				continue;
			}
			int bin = binOf(deltaMasses[i], cfg);
			if (bin >= 0 && bin < binCount) {
				double weight = detectionWeights == null ? 1.0 : detectionWeights[i];
				if (!Double.isFinite(weight) || weight < 0) {
					throw new IllegalArgumentException("detectionWeights must be finite and non-negative");
				}
				histogram[bin] += weight;
				effectiveWeight += weight;
			}
		}

		report(progress, String.format(
				"Delta-mass fit: %d candidates contributed %.1f spectrum-equivalent votes into %d bins of %.4f Da",
				deltaMasses.length, effectiveWeight, binCount, cfg.binWidth));

		double[] smoothed = smooth(histogram, cfg.smoothHalfWidth, cfg.smoothPolyOrder);
		List<int[]> segments = detectSegments(histogram, smoothed, cfg);
		report(progress, "Delta-mass fit: " + segments.size() + " candidate peak segments detected");

		// Discard segments too small to be reportable BEFORE fitting them. See
		// Config.minPeakArea: a peak below this size cannot pass the binomial
		// gate even if every PSM in it is a target, so fitting it buys nothing
		// and inflates the multiple-testing burden on the peaks that matter.
		List<int[]> candidates = new ArrayList<>();
		for (int[] segment : segments) {
			if (segmentArea(histogram, segment[0], segment[1]) >= cfg.minPeakArea) {
				candidates.add(segment);
			}
		}
		int belowArea = segments.size() - candidates.size();

		// Largest first, so that if the cap bites it removes the weakest
		// candidates rather than an arbitrary tail of the mass axis.
		candidates.sort((a, b) -> Double.compare(segmentArea(histogram, b[0], b[1]),
				segmentArea(histogram, a[0], a[1])));
		int overCap = 0;
		if (cfg.maxPeaksToFit > 0 && candidates.size() > cfg.maxPeaksToFit) {
			overCap = candidates.size() - cfg.maxPeaksToFit;
			candidates = new ArrayList<>(candidates.subList(0, cfg.maxPeaksToFit));
		}

		report(progress, String.format(
				"Delta-mass fit: fitting %d segments (%d dropped below minPeakArea=%d, %d beyond maxPeaksToFit=%d) ...",
				candidates.size(), belowArea, cfg.minPeakArea, overCap, cfg.maxPeaksToFit));

		// Fits are independent, so run them in parallel. Ids are assigned after
		// the fact because parallel completion order is not deterministic and
		// the peak list is sorted by centroid immediately below anyway.
		final Config fitCfg = cfg;
		final double[] histogramRef = histogram;
		final java.util.concurrent.atomic.AtomicInteger done = new java.util.concurrent.atomic.AtomicInteger();
		final int totalToFit = candidates.size();
		final int reportEvery = Math.max(1, totalToFit / 10);
		List<Peak> fitted = candidates.parallelStream().map(segment -> {
			Peak peak = fitSegment(histogramRef, segment[0], segment[1], 0, fitCfg);
			int completed = done.incrementAndGet();
			if (completed % reportEvery == 0) {
				// Without this the step is silent for its whole duration and
				// cannot be told apart from a hang. Serialised because the
				// progress sink is a UI/log callback with no thread-safety
				// contract, and this fires ~10 times per run so contention is nil.
				synchronized (PROGRESS_LOCK) {
					report(progress, "Delta-mass fit: " + completed + "/" + totalToFit + " segments fitted");
				}
			}
			return peak;
		}).filter(java.util.Objects::nonNull).collect(java.util.stream.Collectors.toList());

		report(progress, "Delta-mass fit: " + fitted.size() + " identifiable peaks kept (R2 >= " + cfg.minR2
				+ ", sigma >= " + cfg.minSigmaBins + " bins, SNR >= " + cfg.minPeakSnr + ") from "
				+ candidates.size() + " fitted segments");

		fitted.sort((a, b) -> Double.compare(a.centroid, b.centroid));
		int beforeConsolidation = fitted.size();
		fitted = consolidateUnresolvedPeaks(fitted, cfg);
		if (fitted.size() < beforeConsolidation) {
			report(progress, "Delta-mass fit: consolidated " + (beforeConsolidation - fitted.size())
					+ " unresolved alternative fits");
		}
		fitted = reindex(fitted);
		int zeroPeakId = findZeroPeak(fitted, cfg);
		fitted = linkIsotopes(fitted, cfg, zeroPeakId);

		int[] assignment = new int[deltaMasses.length];
		double[] distance = new double[deltaMasses.length];
		double[] sigmas = new double[deltaMasses.length];
		Arrays.fill(assignment, -1);
		Arrays.fill(distance, Double.NaN);
		Arrays.fill(sigmas, Double.NaN);

		report(progress, "Delta-mass fit: assigning " + deltaMasses.length + " PSMs across " + fitted.size()
				+ " peaks ...");
		PeakIndex peakIndex = new PeakIndex(fitted);
		for (int i = 0; i < deltaMasses.length; i++) {
			peakIndex.assign(deltaMasses[i], cfg, i, assignment, distance, sigmas);
		}

		double[] density = localDensity(deltaMasses, detectionWeights, cfg, binCount);

		List<Peak> withCounts = withPsmCounts(fitted, assignment, useForPeakDetection, detectionWeights);

		LOGGER.info("Delta-mass fit: " + withCounts.size() + " peaks over [" + cfg.windowLower + ", " + cfg.windowUpper
				+ "] Da; zero peak "
				+ (zeroPeakId < 0 ? "NOT FOUND" : String.format("centroid=%.5f Da sigma=%.5f Da",
						withCounts.get(zeroPeakId).centroid, withCounts.get(zeroPeakId).sigma)));
		if (zeroPeakId < 0) {
			LOGGER.warn("No zero peak fitted within +/-" + cfg.zeroWindow
					+ " Da. Stage 1 genome discovery has no unmodified set to work from; check the search actually ran "
					+ "open and that calibration is not grossly off.");
		}

		FamilyContext family = normalizeFamilies(deltaMasses, withCounts, assignment, distance, sigmas);
		return new Result(withCounts, zeroPeakId,
				zeroPeakId < 0 ? Double.NaN : withCounts.get(zeroPeakId).centroid, assignment,
				family.assignment, family.normalizedMassDiff, family.isotopeOffset, family.distance,
				family.distanceInSigmas, density);
	}

	private static Result fitKdeInternal(double[] deltaMasses, boolean[] useForPeakDetection,
			double[] detectionWeights, boolean[] isTarget, String[] fileGroups, Config cfg,
			java.util.function.Consumer<String> progress) {
		int binCount = (int) Math.ceil((cfg.windowUpper - cfg.windowLower) / cfg.binWidth);
		if (binCount <= 0) {
			throw new IllegalArgumentException("Invalid delta-mass window: " + cfg.windowLower + ".." + cfg.windowUpper);
		}
		report(progress, "Delta-mass KDE: discovering persistent modes at bandwidths "
				+ Arrays.toString(cfg.kdeBandwidthsDa()) + " Da");
		List<Peak> peaks = DeltaMassKdePeakFinder.discover(deltaMasses, useForPeakDetection, detectionWeights,
				isTarget, fileGroups, cfg, progress);
		peaks = reindex(peaks);
		int zeroPeakId = findZeroPeak(peaks, cfg);
		peaks = linkIsotopes(peaks, cfg, zeroPeakId);

		int[] assignment = new int[deltaMasses.length];
		double[] distance = new double[deltaMasses.length];
		double[] sigmas = new double[deltaMasses.length];
		Arrays.fill(assignment, -1);
		Arrays.fill(distance, Double.NaN);
		Arrays.fill(sigmas, Double.NaN);
		PeakIndex index = new PeakIndex(peaks);
		for (int i = 0; i < deltaMasses.length; i++) {
			index.assign(deltaMasses[i], cfg, i, assignment, distance, sigmas);
		}

		List<Peak> withCounts = withPsmCounts(peaks, assignment, useForPeakDetection, detectionWeights);
		FamilyContext family = normalizeFamilies(deltaMasses, withCounts, assignment, distance, sigmas);
		double[] density = localDensity(deltaMasses, detectionWeights, cfg, binCount);
		long isotopeMembers = withCounts.stream().filter(Peak::isIsotope).count();
		report(progress, "Delta-mass KDE: retained " + withCounts.size() + " persistent modes; "
				+ isotopeMembers + " satellite modes normalized into canonical families");
		return new Result(withCounts, zeroPeakId,
				zeroPeakId < 0 ? Double.NaN : withCounts.get(zeroPeakId).centroid, assignment,
				family.assignment, family.normalizedMassDiff, family.isotopeOffset, family.distance,
				family.distanceInSigmas, density);
	}

	private static FamilyContext normalizeFamilies(double[] deltaMasses, List<Peak> peaks, int[] observed,
			double[] observedDistance, double[] observedSigmas) {
		int[] family = Arrays.copyOf(observed, observed.length);
		double[] normalized = Arrays.copyOf(deltaMasses, deltaMasses.length);
		int[] offsets = new int[observed.length];
		double[] distance = Arrays.copyOf(observedDistance, observedDistance.length);
		double[] sigmas = Arrays.copyOf(observedSigmas, observedSigmas.length);
		for (int i = 0; i < observed.length; i++) {
			if (observed[i] < 0) {
				continue;
			}
			Peak mode = peaks.get(observed[i]);
			int familyId = mode.familyPeakId();
			Peak root = peaks.get(familyId);
			family[i] = familyId;
			offsets[i] = mode.isotopeOffset;
			normalized[i] = deltaMasses[i] - mode.isotopeOffset * AminoAcidProperty.isotope_deltaM;
			distance[i] = normalized[i] - root.centroid;
			sigmas[i] = Math.abs(distance[i]) / Math.max(root.sigma, 1e-12);
		}
		return new FamilyContext(family, normalized, offsets, distance, sigmas);
	}

	private static final class FamilyContext {
		final int[] assignment;
		final double[] normalizedMassDiff;
		final int[] isotopeOffset;
		final double[] distance;
		final double[] distanceInSigmas;

		FamilyContext(int[] assignment, double[] normalizedMassDiff, int[] isotopeOffset, double[] distance,
				double[] distanceInSigmas) {
			this.assignment = assignment;
			this.normalizedMassDiff = normalizedMassDiff;
			this.isotopeOffset = isotopeOffset;
			this.distance = distance;
			this.distanceInSigmas = distanceInSigmas;
		}
	}

	private static int binOf(double deltaMass, Config cfg) {
		return (int) Math.floor((deltaMass - cfg.windowLower) / cfg.binWidth);
	}

	private static double binCentre(int bin, Config cfg) {
		return cfg.windowLower + (bin + 0.5) * cfg.binWidth;
	}

	/**
	 * Savitzky-Golay smoothing. Preferred over a moving average because it
	 * preserves peak height and width, which a box filter flattens - and the
	 * fitted amplitude and sigma are themselves classifier features. Coefficients
	 * are derived here rather than pulled from a table or an extra dependency:
	 * the smoothing weights are the first row of (AᵀA)⁻¹Aᵀ for the Vandermonde
	 * matrix A of the window offsets.
	 */
	static double[] smooth(double[] histogram, int halfWidth, int polyOrder) {
		if (halfWidth <= 0 || polyOrder < 0 || polyOrder >= 2 * halfWidth + 1) {
			double[] copy = new double[histogram.length];
			for (int i = 0; i < histogram.length; i++) {
				copy[i] = histogram[i];
			}
			return copy;
		}
		double[] coefficients = savitzkyGolayCoefficients(halfWidth, polyOrder);
		double[] out = new double[histogram.length];
		for (int i = 0; i < histogram.length; i++) {
			double acc = 0;
			for (int k = -halfWidth; k <= halfWidth; k++) {
				int j = i + k;
				if (j < 0 || j >= histogram.length) {
					continue; // zero-padded at the edges
				}
				acc += coefficients[k + halfWidth] * histogram[j];
			}
			out[i] = acc;
		}
		return out;
	}

	static double[] savitzkyGolayCoefficients(int halfWidth, int polyOrder) {
		int n = 2 * halfWidth + 1;
		double[][] a = new double[n][polyOrder + 1];
		for (int i = 0; i < n; i++) {
			int x = i - halfWidth;
			double v = 1.0;
			for (int j = 0; j <= polyOrder; j++) {
				a[i][j] = v;
				v *= x;
			}
		}
		RealMatrix matrixA = new Array2DRowRealMatrix(a, false);
		RealMatrix ata = matrixA.transpose().multiply(matrixA);
		RealMatrix pseudo = new LUDecomposition(ata).getSolver().getInverse().multiply(matrixA.transpose());
		double[] coefficients = new double[n];
		for (int i = 0; i < n; i++) {
			coefficients[i] = pseudo.getEntry(0, i);
		}
		return coefficients;
	}

	/**
	 * Finds candidate peak segments as {@code [startBin, endBin]} pairs: each
	 * local maximum of the smoothed histogram that clears {@code minPeakCount}
	 * in the raw data, expanded outwards until the smoothed signal stops
	 * descending (a valley) or leaves the axis.
	 */
	private static List<int[]> detectSegments(double[] histogram, double[] smoothed, Config cfg) {
		List<int[]> segments = new ArrayList<>();
		int i = 1;
		while (i < smoothed.length - 1) {
			if (!(smoothed[i] > smoothed[i - 1] && smoothed[i] >= smoothed[i + 1])) {
				i++;
				continue;
			}
			if (histogram[i] < cfg.minPeakCount) {
				i++;
				continue;
			}
			int left = i;
			while (left > 0 && smoothed[left - 1] < smoothed[left] && smoothed[left - 1] > 0) {
				left--;
			}
			int right = i;
			while (right < smoothed.length - 1 && smoothed[right + 1] < smoothed[right] && smoothed[right + 1] > 0) {
				right++;
			}
			if (right - left >= 2) {
				segments.add(new int[] { left, right });
			}
			i = Math.max(right, i + 1);
		}
		return segments;
	}

	/**
	 * Fits a Gaussian signal plus a constant local background to the
	 * <em>raw</em> support across a segment. Smoothed data is used only to find the
	 * segment; fitting the smoothed curve would bias sigma upwards. The legacy
	 * zero-background Gaussian is retained as a numerical fallback.
	 */
	/** Total spectrum-equivalent support across a segment. */
	private static double segmentArea(double[] histogram, int startBin, int endBin) {
		double area = 0;
		for (int bin = startBin; bin <= endBin; bin++) {
			area += histogram[bin];
		}
		return area;
	}

	private static Peak fitSegment(double[] histogram, int startBin, int endBin, int id, Config cfg) {
		WeightedObservedPoints observed = new WeightedObservedPoints();
		int points = endBin - startBin + 1;
		if (points < Math.max(4, cfg.minFitBins)) {
			return null;
		}
		for (int bin = startBin; bin <= endBin; bin++) {
			observed.add(binCentre(bin, cfg), histogram[bin]);
		}

		double[] params;
		try {
			double left = histogram[startBin];
			double right = histogram[endBin];
			double background = Math.max(0.0, (left + right) / 2.0);
			int apex = startBin;
			for (int bin = startBin + 1; bin <= endBin; bin++) {
				if (histogram[bin] > histogram[apex]) {
					apex = bin;
				}
			}
			double amplitude = Math.max(cfg.binWidth, histogram[apex] - background);
			double sigma = Math.max(cfg.binWidth, (endBin - startBin + 1) * cfg.binWidth / 6.0);
			params = SimpleCurveFitter.create(GAUSSIAN_WITH_BASELINE,
					new double[] { background, amplitude, binCentre(apex, cfg), sigma })
					.withMaxIterations(cfg.maxFitIterations).fit(observed.toList());
		} catch (RuntimeException e) {
			// Commons Math throws a family of unchecked convergence/dimension
			// exceptions here; a segment that will not fit is simply not a peak.
			LOGGER.debug("Gaussian fit failed for segment [" + startBin + "," + endBin + "]: " + e.getMessage());
			return fitPureGaussian(histogram, startBin, endBin, id, cfg);
		}

		double background = params[0];
		double amplitude = params[1];
		double centroid = params[2];
		double sigma = Math.abs(params[3]);
		if (!(sigma > 0) || !(amplitude > 0) || background < 0 || !Double.isFinite(centroid)
				|| !Double.isFinite(amplitude) || !Double.isFinite(background)) {
			return fitPureGaussian(histogram, startBin, endBin, id, cfg);
		}
		double segmentLo = binCentre(startBin, cfg) - cfg.binWidth;
		double segmentHi = binCentre(endBin, cfg) + cfg.binWidth;
		if (centroid < segmentLo || centroid > segmentHi || sigma > segmentHi - segmentLo) {
			return fitPureGaussian(histogram, startBin, endBin, id, cfg);
		}

		double r2 = rSquared(histogram, startBin, endBin, background, amplitude, centroid, sigma, cfg);
		if (r2 < cfg.minR2) {
			return fitPureGaussian(histogram, startBin, endBin, id, cfg);
		}
		if (!passesShapeGate(sigma, amplitude, background, cfg)) {
			return null;
		}
		return new Peak(id, centroid, sigma, amplitude, background, r2, 0, 0, -1);
	}

	private static Peak fitPureGaussian(double[] histogram, int startBin, int endBin, int id, Config cfg) {
		WeightedObservedPoints observed = new WeightedObservedPoints();
		int points = 0;
		for (int bin = startBin; bin <= endBin; bin++) {
			if (histogram[bin] > 0) {
				observed.add(binCentre(bin, cfg), histogram[bin]);
				points++;
			}
		}
		if (points < 3) {
			return null;
		}
		try {
			double[] params = GaussianCurveFitter.create().withMaxIterations(cfg.maxFitIterations)
					.fit(observed.toList());
			double amplitude = params[0];
			double centroid = params[1];
			double sigma = Math.abs(params[2]);
			if (!(sigma > 0) || !(amplitude > 0) || !Double.isFinite(centroid)) {
				return null;
			}
			double r2 = rSquared(histogram, startBin, endBin, 0.0, amplitude, centroid, sigma, cfg);
			return r2 >= cfg.minR2 && passesShapeGate(sigma, amplitude, 0.0, cfg)
					? new Peak(id, centroid, sigma, amplitude, 0.0, r2, 0, 0, -1)
					: null;
		} catch (RuntimeException e) {
			return null;
		}
	}

	private static boolean passesShapeGate(double sigma, double amplitude, double background, Config cfg) {
		if (sigma < cfg.minSigmaBins * cfg.binWidth) {
			return false;
		}
		double snr = amplitude / Math.sqrt(background + 1.0);
		return Double.isFinite(snr) && snr >= cfg.minPeakSnr;
	}

	private static double rSquared(double[] histogram, int startBin, int endBin, double background, double amplitude,
			double centroid, double sigma, Config cfg) {
		double mean = 0;
		int n = 0;
		for (int bin = startBin; bin <= endBin; bin++) {
			mean += histogram[bin];
			n++;
		}
		if (n == 0) {
			return 0;
		}
		mean /= n;

		double ssTot = 0;
		double ssRes = 0;
		for (int bin = startBin; bin <= endBin; bin++) {
			double x = binCentre(bin, cfg);
			double predicted = background
					+ amplitude * Math.exp(-((x - centroid) * (x - centroid)) / (2 * sigma * sigma));
			double actual = histogram[bin];
			ssTot += (actual - mean) * (actual - mean);
			ssRes += (actual - predicted) * (actual - predicted);
		}
		if (ssTot <= 0) {
			return 0;
		}
		return 1.0 - (ssRes / ssTot);
	}

	private static List<Peak> reindex(List<Peak> peaks) {
		List<Peak> out = new ArrayList<>(peaks.size());
		for (int i = 0; i < peaks.size(); i++) {
			Peak p = peaks.get(i);
			out.add(new Peak(i, p.centroid, p.sigma, p.amplitude, p.background, p.r2, p.psmCount,
					p.discoverySupport, p.isotopeOffset, p.monoPeakId, p.lowerBoundary, p.upperBoundary,
					p.persistentScales, p.reproducibleFiles, p.discoveryMethod));
		}
		return out;
	}

	/**
	 * Removes alternative fits that the histogram cannot resolve as separate
	 * peaks. This deliberately uses fitted widths, not the wider membership
	 * floor: two peaks may be resolvable even when their instrument-accuracy
	 * assignment windows overlap.
	 */
	static List<Peak> consolidateUnresolvedPeaks(List<Peak> peaks, Config cfg) {
		if (peaks.size() < 2) {
			return new ArrayList<>(peaks);
		}
		List<Peak> strongestFirst = new ArrayList<>(peaks);
		strongestFirst.sort((a, b) -> {
			double scoreA = a.amplitude * a.sigma * Math.max(0.0, a.r2);
			double scoreB = b.amplitude * b.sigma * Math.max(0.0, b.r2);
			int byScore = Double.compare(scoreB, scoreA);
			return byScore != 0 ? byScore : Double.compare(a.centroid, b.centroid);
		});

		List<Peak> kept = new ArrayList<>();
		for (Peak candidate : strongestFirst) {
			boolean unresolved = false;
			for (Peak accepted : kept) {
				double absoluteFloor = cfg.minPeakSeparationBins * cfg.binWidth;
				double widthFloor = cfg.minPeakSeparationSigmas
						* Math.hypot(candidate.sigma, accepted.sigma);
				double required = Math.max(absoluteFloor, widthFloor);
				if (Math.abs(candidate.centroid - accepted.centroid) < required) {
					unresolved = true;
					break;
				}
			}
			if (!unresolved) {
				kept.add(candidate);
			}
		}
		kept.sort((a, b) -> Double.compare(a.centroid, b.centroid));
		return kept;
	}

	/**
	 * Marks peaks that sit one or more isotope spacings above another peak.
	 * <p>
	 * A precursor triggered on the +1 (or +2) isotopic peak rather than the
	 * monoisotopic peak yields a measured mass ~{@code isotope_deltaM} higher
	 * than the true peptide mass, so it appears as a delta-mass peak at
	 * {@code centroid(mono) + k * isotope_deltaM}. These are selection
	 * artefacts, not independent modifications. Their observed modes remain
	 * auditable while family assignments and normalized masses point to one root.
	 */
	/**
	 * @param peaks sorted by centroid and reindexed, so {@code peak.id} is the
	 *              list position — which is what makes the binary-searched
	 *              lookup below able to return an id directly.
	 */
	static List<Peak> linkIsotopes(List<Peak> peaks, Config cfg, int zeroPeakId) {
		if (peaks.isEmpty()) {
			return peaks;
		}

		/*
		 * Families are deliberately stars, never transitive graphs. On a dense open
		 * search axis, unrelated modifications frequently differ by approximately
		 * one nominal neutron. Connecting every such pair and then taking connected
		 * components percolates across the axis and produces impossible offsets such
		 * as +20. A canonical mode may therefore recruit only direct satellites in
		 * [-maxIsotopeOffset,+maxIsotopeOffset], and a recruited satellite can never
		 * recruit another mode.
		 *
		 * The calibrated zero mode is considered first and is always its own root.
		 * Remaining roots are considered by decreasing target discovery support, so
		 * the dominant mode becomes canonical while weaker isotope-selection modes
		 * are evidence for it. Exact spacing is necessary but still not proof; final
		 * family acceptance remains target/decoy and localization based.
		 */
		PeakIndex index = new PeakIndex(peaks);
		int[] offsets = new int[peaks.size()];
		int[] roots = new int[peaks.size()];
		Arrays.fill(roots, -1);
		boolean[] assigned = new boolean[peaks.size()];
		List<Integer> rootOrder = new ArrayList<>(peaks.size());
		if (zeroPeakId >= 0) {
			rootOrder.add(zeroPeakId);
		}
		peaks.stream().map(peak -> peak.id).filter(id -> id != zeroPeakId)
				.sorted((a, b) -> Double.compare(isotopeRootStrength(peaks.get(b)), isotopeRootStrength(peaks.get(a))))
				.forEach(rootOrder::add);

		for (int rootId : rootOrder) {
			if (assigned[rootId]) {
				continue;
			}
			Peak root = peaks.get(rootId);
			assigned[rootId] = true;
			roots[rootId] = rootId;
			List<IsotopeEdge> satellites = new ArrayList<>();
			collectIsotopeDirection(root, -1, peaks, index, assigned, cfg, satellites);
			collectIsotopeDirection(root, 1, peaks, index, assigned, cfg, satellites);
			if (rootId == zeroPeakId || satellites.size() >= cfg.minIsotopeFamilySatellites) {
				for (IsotopeEdge satellite : satellites) {
					assigned[satellite.to] = true;
					offsets[satellite.to] = satellite.offset;
					roots[satellite.to] = rootId;
				}
			}
		}

		List<Peak> out = new ArrayList<>(peaks.size());
		for (Peak peak : peaks) {
			int offset = offsets[peak.id];
			int monoPeakId = offset == 0 ? -1 : roots[peak.id];
			out.add(new Peak(peak.id, peak.centroid, peak.sigma, peak.amplitude, peak.background, peak.r2,
					peak.psmCount, peak.discoverySupport, offset, monoPeakId, peak.lowerBoundary,
					peak.upperBoundary, peak.persistentScales, peak.reproducibleFiles, peak.discoveryMethod));
		}
		return out;
	}

	private static double isotopeRootStrength(Peak peak) {
		if (peak.amplitude > 0.0 && Double.isFinite(peak.amplitude)) {
			return peak.amplitude;
		}
		if (peak.discoverySupport > 0.0) {
			return peak.discoverySupport;
		}
		if (peak.psmCount > 0) {
			return peak.psmCount;
		}
		return peak.amplitude;
	}

	private static void collectIsotopeDirection(Peak root, int direction, List<Peak> peaks, PeakIndex index,
			boolean[] assigned, Config cfg, List<IsotopeEdge> out) {
		double previousStrength = isotopeRootStrength(root);
		for (int magnitude = 1; magnitude <= cfg.maxIsotopeOffset; magnitude++) {
			int offset = direction * magnitude;
			double expected = root.centroid + offset * AminoAcidProperty.isotope_deltaM;
			int match = index.findNearestWithin(expected, cfg.isotopeTolerance, root.id);
			if (match < 0 || assigned[match]
					|| !isPlausibleIsotopeLink(root, peaks.get(match), expected, cfg)) {
				break;
			}
			double strength = isotopeRootStrength(peaks.get(match));
			if (strength > previousStrength * cfg.isotopeMaxAbundanceIncrease) {
				break;
			}
			out.add(new IsotopeEdge(match, offset));
			previousStrength = strength;
		}
	}

	private static boolean isPlausibleIsotopeLink(Peak first, Peak second, double expected, Config cfg) {
		double combinedSigma = Math.hypot(first.sigma, second.sigma);
		double widthScaled = cfg.isotopeSigmaMultiplier * combinedSigma;
		double allowed = Math.min(cfg.isotopeTolerance, Math.max(cfg.isotopeToleranceFloorDa, widthScaled));
		return Math.abs(second.centroid - expected) <= allowed;
	}

	private static final class IsotopeEdge {
		final int to;
		final int offset;

		IsotopeEdge(int to, int offset) {
			this.to = to;
			this.offset = offset;
		}
	}

	/**
	 * The zero peak is the non-isotope peak whose centroid lies closest to 0
	 * within {@code zeroWindow}. Its centroid is the systematic calibration
	 * offset, which is why stage 1 selects on the fitted peak rather than on a
	 * window centred on exactly 0.0.
	 */
	private static int findZeroPeak(List<Peak> peaks, Config cfg) {
		int best = -1;
		double bestDistance = Double.MAX_VALUE;
		for (Peak peak : peaks) {
			if (peak.isIsotope()) {
				continue;
			}
			double distance = Math.abs(peak.centroid);
			if (distance <= cfg.zeroWindow && distance < bestDistance) {
				bestDistance = distance;
				best = peak.id;
			}
		}
		return best;
	}

	/**
	 * Centroid-sorted view of the fitted peaks, supporting nearest-peak and
	 * isotope lookups in logarithmic rather than linear time.
	 * <p>
	 * Both lookups were linear scans over the peak list. An open search over
	 * [-150, 500] Da resolves peaks in the thousands, so nearest-peak assignment
	 * cost peaks x PSMs — of order 10^10 for a million candidates against ten
	 * thousand peaks, which is hours of wall clock and the dominant cost of the
	 * whole fitting step. The list is already sorted by centroid before either
	 * runs, so the ordering needed for binary search was there to be used.
	 *
	 * <h2>Why nearest-in-Da is not enough</h2>
	 * Proximity is measured in <em>sigmas</em>, and peaks have different widths,
	 * so the closest peak in daltons need not be the closest in sigmas: a broad
	 * peak further away can still be nearer once scaled. Checking only the
	 * immediate neighbours would therefore change results, not just speed.
	 * <p>
	 * The search instead walks outwards from the insertion point and stops using
	 * a bound that cannot skip a better peak: once the dalton distance to the
	 * next candidate, divided by the <em>widest</em> sigma in the whole set,
	 * already exceeds the best sigma distance found so far, no remaining peak on
	 * that side can improve on it. Results are identical to the linear scan.
	 */
	static final class PeakIndex {
		private final double[] centroids;
		private final double[] sigmas;
		private final double[] lowerBoundaries;
		private final double[] upperBoundaries;
		private final double maxSigma;
		private final double maxMembershipHalfWidth;
		private final int size;

		PeakIndex(List<Peak> sortedByCentroid) {
			this.size = sortedByCentroid.size();
			this.centroids = new double[size];
			this.sigmas = new double[size];
			this.lowerBoundaries = new double[size];
			this.upperBoundaries = new double[size];
			double widest = 0;
			double widestMembership = 0;
			for (int i = 0; i < size; i++) {
				Peak peak = sortedByCentroid.get(i);
				centroids[i] = peak.centroid;
				sigmas[i] = peak.sigma;
				lowerBoundaries[i] = peak.lowerBoundary;
				upperBoundaries[i] = peak.upperBoundary;
				widest = Math.max(widest, peak.sigma);
				if (Double.isFinite(peak.lowerBoundary) && Double.isFinite(peak.upperBoundary)) {
					widestMembership = Math.max(widestMembership,
							Math.max(peak.centroid - peak.lowerBoundary, peak.upperBoundary - peak.centroid));
				}
			}
			this.maxSigma = widest;
			this.maxMembershipHalfWidth = widestMembership;
		}

		/** Index of the first centroid >= value, or {@code size} if none. */
		private int lowerBound(double value) {
			int lo = 0;
			int hi = size;
			while (lo < hi) {
				int mid = (lo + hi) >>> 1;
				if (centroids[mid] < value) {
					lo = mid + 1;
				} else {
					hi = mid;
				}
			}
			return lo;
		}

		void assign(double deltaMass, Config cfg, int index, int[] assignment, double[] distance,
				double[] sigmasOut) {
			if (size == 0) {
				return;
			}
			int start = lowerBound(deltaMass);

			int bestId = -1;
			double bestSigmas = Double.MAX_VALUE;
			double bestSigned = Double.NaN;

			// Walk outwards from the insertion point in both directions.
			for (int right = start; right < size; right++) {
				double signed = deltaMass - centroids[right];
				if (Math.abs(signed) / maxSigma >= bestSigmas) {
					break; // no peak further right can beat the incumbent
				}
				double inSigmas = Math.abs(signed) / sigmas[right];
				if (inSigmas < bestSigmas) {
					bestSigmas = inSigmas;
					bestSigned = signed;
					bestId = right;
				}
			}
			for (int left = start - 1; left >= 0; left--) {
				double signed = deltaMass - centroids[left];
				if (Math.abs(signed) / maxSigma >= bestSigmas) {
					break; // nor any further left
				}
				double inSigmas = Math.abs(signed) / sigmas[left];
				if (inSigmas < bestSigmas) {
					bestSigmas = inSigmas;
					bestSigned = signed;
					bestId = left;
				}
			}

			if (bestId < 0) {
				return;
			}
			/*
			 * Choose among peaks whose actual membership intervals contain the PSM.
			 * The former implementation chose the globally smallest sigma-distance
			 * first and then tested only that peak's interval. A broad peak could be
			 * sigma-nearer but outside its interval while a neighbouring narrow
			 * peak was inside the 10 mDa floor, leaving a false assignment gap.
			 *
			 * We still rank eligible peaks by sigma distance, so a broad peak cannot
			 * steal a PSM from a narrower peak's core merely because their floored
			 * membership intervals overlap.
			 */
			double maxHalfWidth = Math.max(maxMembershipHalfWidth,
					Math.max(cfg.membershipFloorDa, cfg.membershipSigmas * maxSigma));
			int eligibleStart = lowerBound(deltaMass - maxHalfWidth);
			int bestEligibleId = -1;
			double bestEligibleSigmas = Double.MAX_VALUE;
			double bestEligibleSigned = Double.NaN;
			for (int candidate = eligibleStart; candidate < size
					&& centroids[candidate] <= deltaMass + maxHalfWidth; candidate++) {
				double signed = deltaMass - centroids[candidate];
				boolean empirical = Double.isFinite(lowerBoundaries[candidate])
						&& Double.isFinite(upperBoundaries[candidate]);
				boolean contained = empirical
						? deltaMass >= lowerBoundaries[candidate] && deltaMass <= upperBoundaries[candidate]
						: Math.abs(signed) <= Math.max(cfg.membershipSigmas * sigmas[candidate],
								cfg.membershipFloorDa);
				if (!contained) {
					continue;
				}
				double inSigmas = Math.abs(signed) / sigmas[candidate];
				if (inSigmas < bestEligibleSigmas) {
					bestEligibleSigmas = inSigmas;
					bestEligibleSigned = signed;
					bestEligibleId = candidate;
				}
			}

			if (bestEligibleId >= 0) {
				assignment[index] = bestEligibleId;
				distance[index] = bestEligibleSigned;
				sigmasOut[index] = bestEligibleSigmas;
			} else {
				// Membership is a feature, not a filter. Preserve distance to the
				// nearest standardized peak for the classifier even when unassigned.
				distance[index] = bestSigned;
				sigmasOut[index] = bestSigmas;
			}
		}

		/** Closest centroid within tolerance, excluding self; -1 if none. */
		int findNearestWithin(double target, double tolerance, int selfIndex) {
			int best = -1;
			double bestDistance = Double.MAX_VALUE;
			for (int i = lowerBound(target - tolerance); i < size; i++) {
				if (centroids[i] > target + tolerance) {
					break;
				}
				if (i != selfIndex) {
					double candidateDistance = Math.abs(centroids[i] - target);
					if (candidateDistance < bestDistance) {
						bestDistance = candidateDistance;
						best = i;
					}
				}
			}
			return best;
		}
	}

	/**
	 * Local delta-mass density, weighted exactly as the fitting histogram is.
	 * <p>
	 * {@code detectionWeights} may be null, in which case every candidate counts
	 * 1 - the old behaviour, kept only for callers that have no weights.
	 */
	private static double[] localDensity(double[] deltaMasses, double[] detectionWeights, Config cfg,
			int binCount) {
		double[] histogram = new double[binCount];
		for (int i = 0; i < deltaMasses.length; i++) {
			int bin = binOf(deltaMasses[i], cfg);
			if (bin >= 0 && bin < binCount) {
				double weight = detectionWeights == null ? 1.0 : detectionWeights[i];
				if (Double.isFinite(weight) && weight > 0) {
					histogram[bin] += weight;
				}
			}
		}
		int halfBins = Math.max(1, (int) Math.round(cfg.densityWindow / cfg.binWidth));
		double[] prefix = new double[binCount + 1];
		for (int i = 0; i < binCount; i++) {
			prefix[i + 1] = prefix[i] + histogram[i];
		}
		double[] density = new double[deltaMasses.length];
		for (int i = 0; i < deltaMasses.length; i++) {
			int bin = binOf(deltaMasses[i], cfg);
			if (bin < 0 || bin >= binCount) {
				density[i] = 0;
				continue;
			}
			int lo = Math.max(0, bin - halfBins);
			int hi = Math.min(binCount, bin + halfBins + 1);
			density[i] = prefix[hi] - prefix[lo];
		}
		return density;
	}

	private static List<Peak> withPsmCounts(List<Peak> peaks, int[] assignment, boolean[] useForPeakDetection,
			double[] detectionWeights) {
		int[] counts = new int[peaks.size()];
		double[] support = new double[peaks.size()];
		for (int i = 0; i < assignment.length; i++) {
			int a = assignment[i];
			if (a >= 0) {
				counts[a]++;
				if (useForPeakDetection == null || useForPeakDetection[i]) {
					support[a] += detectionWeights == null ? 1.0 : detectionWeights[i];
				}
			}
		}
		List<Peak> out = new ArrayList<>(peaks.size());
		for (Peak p : peaks) {
			out.add(new Peak(p.id, p.centroid, p.sigma, p.amplitude, p.background, p.r2, counts[p.id],
					support[p.id], p.isotopeOffset, p.monoPeakId, p.lowerBoundary, p.upperBoundary,
					p.persistentScales, p.reproducibleFiles, p.discoveryMethod));
		}
		return out;
	}
}
