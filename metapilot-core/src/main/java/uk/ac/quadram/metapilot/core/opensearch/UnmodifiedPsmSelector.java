package uk.ac.quadram.metapilot.core.opensearch;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * Stage 1 of the open-search pipeline: selects the unmodified PSMs that feed
 * genome discovery.
 *
 * <h2>Why stage 1 runs an open search at all</h2>
 * Only unmodified PSMs are used to map genomes, so a closed search looks
 * sufficient — but it is not, and the difference is not merely one of yield. In
 * a <b>closed</b> search a spectrum from a modified peptide either fails to
 * match or, worse, is force-matched to whatever wrong peptide happens to fall
 * inside the narrow precursor window; that wrong match can pull an incorrect
 * genome into the refined database. In an <b>open</b> search the same spectrum
 * matches its true peptide carrying a mass shift, so it is diverted out of the
 * zero peak rather than contaminating it. The surviving zero-peak set is
 * therefore <em>purer</em> for genome mapping, not just smaller.
 * <p>
 * Measured, the gain did not justify the cost, and the default marker search is
 * now closed ({@code round1_search_mode}). On the 12-species mix the closed and
 * open marker searches selected the same spiked genomes (11.3 vs 11.4 per run)
 * at about a fifth of the run time; the closed lists carried more non-spiked
 * genomes (32 vs 20 distinct genomes across 9 runs), which the genome test
 * downstream has to remove. Open mode remains available.
 *
 * <h2>Why the fitted peak, not a fixed tolerance</h2>
 * The zero peak's fitted centroid absorbs any systematic mass calibration
 * offset. A window hard-centred on 0.0 would silently clip one side of the peak
 * on any run with calibration drift; a fitted peak tracks it.
 *
 * <h2>Why isotope modes are normalized</h2>
 * A precursor triggered on the +1 or +2 isotopic peak yields a delta mass near
 * +1.003 or +2.006 Da. Those are precursor-selection artefacts, not separate
 * modifications. Their raw mass and observed offset stay auditable, while the
 * PSM is normalized onto the calibrated zero-family root before genome mapping.
 * <p>
 * The offset is read from precursor minus calculated mass, so it stays visible
 * whether the search found the spectrum through a wide window
 * ({@code isotope_error = 0}, open) or through MSFragger's isotope correction
 * ({@code isotope_error = 0/1/2}, closed). {@link #warnIfIsotopePeaksMissing}
 * reports a run in which none are seen.
 *
 * <h2>Selection is a hard cut here, deliberately</h2>
 * Elsewhere in this pipeline every condition is a classifier feature rather
 * than a filter. Stage 1 is the exception, and consciously so: its job is
 * genome mapping, not PTM discovery, and running the full validation machinery
 * here would be circular — a well-trained classifier wants the refined database
 * that this step exists to produce. The no-hard-filters rule governs stage 2,
 * where validation decisions are actually made.
 */
public final class UnmodifiedPsmSelector {

	private static final Logger LOGGER = LogManager.getLogger(UnmodifiedPsmSelector.class);

	public static final class Result {
		/** PSMs inside the fitted zero family; the input to genome mapping. */
		public final List<OpenSearchPsm> unmodified;
		/**
		 * Every zero-family member, before the q-value cut. Exposed so a second
		 * stage-1 pass can make its assignment exclusive against this family rather
		 * than re-deriving membership from the centroid and re-folding isotopes,
		 * which would not reproduce the fitter's decisions exactly. See
		 * {@link ExpectedPeakSelector}.
		 */
		public final List<OpenSearchPsm> zeroFamily;
		/**
		 * Peptides of the accepted TARGET PSMs. The corroboration criterion in
		 * {@link ExpectedPeakSelector} tests a shifted peak's peptides against this
		 * set, so it is built here where the accepted set is already known.
		 */
		public final Set<String> acceptedTargetPeptides;
		/** Fitted centroid of the zero peak = the systematic calibration offset (Da). */
		public final double calibrationOffset;
		/** Fitted width of the zero peak (Da). */
		public final double zeroPeakSigma;
		/** Empirical lower boundary of the zero-family KDE basin (Da). */
		public final double zeroPeakLowerBoundary;
		/** Empirical upper boundary of the zero-family KDE basin (Da). */
		public final double zeroPeakUpperBoundary;
		/** Discovery model used for the zero family. */
		public final String discoveryMethod;
		/**
		 * PSMs outside the zero peak. Not broken down by cause: distinguishing
		 * "isotope" from "modified" would require fitting the rest of the
		 * delta-mass axis, which stage 1 deliberately does not do.
		 */
		public final int excluded;
		/**
		 * PSMs sitting near a +1/+2 isotope offset. Informational; exact members may
		 * be normalized into the zero family. Also detects an
		 * {@code isotope_error != 0} misconfiguration.
		 */
		public final int nearIsotopeOffsets;
		/** Target PSMs in the zero family before the q-value cut. */
		public final int targetsBeforeFdr;
		/** Decoy PSMs in the zero family before the q-value cut. */
		public final int decoysBeforeFdr;
		/** Target PSMs in the zero peak surviving the stage-1 FDR cut. */
		public final int targetsAfterFdr;
		/** Decoy PSMs surviving the stage-1 FDR cut. */
		public final int decoysAfterFdr;
		/** Achieved target-decoy estimate, decoys/targets, at the accepted boundary. */
		public final double estimatedFdr;
		/** q-value cut applied to the zero-peak set. */
		public final double fdrThreshold;

		Result(List<OpenSearchPsm> unmodified, double calibrationOffset, double zeroPeakSigma,
				double zeroPeakLowerBoundary, double zeroPeakUpperBoundary, String discoveryMethod, int excluded,
				int nearIsotopeOffsets, int targetsBeforeFdr, int decoysBeforeFdr, int targetsAfterFdr,
				int decoysAfterFdr, double estimatedFdr, double fdrThreshold) {
			this(unmodified, calibrationOffset, zeroPeakSigma, zeroPeakLowerBoundary, zeroPeakUpperBoundary,
					discoveryMethod, excluded, nearIsotopeOffsets, targetsBeforeFdr, decoysBeforeFdr, targetsAfterFdr,
					decoysAfterFdr, estimatedFdr, fdrThreshold, new ArrayList<>(), Collections.emptySet());
		}

		Result(List<OpenSearchPsm> unmodified, double calibrationOffset, double zeroPeakSigma,
				double zeroPeakLowerBoundary, double zeroPeakUpperBoundary, String discoveryMethod, int excluded,
				int nearIsotopeOffsets, int targetsBeforeFdr, int decoysBeforeFdr, int targetsAfterFdr,
				int decoysAfterFdr, double estimatedFdr, double fdrThreshold, List<OpenSearchPsm> zeroFamily,
				Set<String> acceptedTargetPeptides) {
			this.zeroFamily = zeroFamily;
			this.acceptedTargetPeptides = acceptedTargetPeptides;
			this.unmodified = unmodified;
			this.calibrationOffset = calibrationOffset;
			this.zeroPeakSigma = zeroPeakSigma;
			this.zeroPeakLowerBoundary = zeroPeakLowerBoundary;
			this.zeroPeakUpperBoundary = zeroPeakUpperBoundary;
			this.discoveryMethod = discoveryMethod;
			this.excluded = excluded;
			this.nearIsotopeOffsets = nearIsotopeOffsets;
			this.targetsBeforeFdr = targetsBeforeFdr;
			this.decoysBeforeFdr = decoysBeforeFdr;
			this.targetsAfterFdr = targetsAfterFdr;
			this.decoysAfterFdr = decoysAfterFdr;
			this.estimatedFdr = estimatedFdr;
			this.fdrThreshold = fdrThreshold;
		}
	}

	private UnmodifiedPsmSelector() {
	}

	/** Default stage-1 q-value cut on the zero-peak set. */
	public static final double DEFAULT_FDR = 0.01;

	public static Result select(List<OpenSearchPsm> psms, DeltaMassPeakFitter.Config fitConfig) {
		return select(psms, fitConfig, DEFAULT_FDR, null);
	}

	public static Result select(List<OpenSearchPsm> psms, DeltaMassPeakFitter.Config fitConfig,
			java.util.function.Consumer<String> progress) {
		return select(psms, fitConfig, DEFAULT_FDR, progress);
	}

	public static Result select(List<OpenSearchPsm> psms, DeltaMassPeakFitter.Config fitConfig, double fdrThreshold,
			java.util.function.Consumer<String> progress) {
		if (fitConfig == null) {
			fitConfig = new DeltaMassPeakFitter.Config();
		}

		double[] deltaMasses = new double[psms.size()];
		boolean[] targetsOnly = new boolean[psms.size()];
		for (int i = 0; i < psms.size(); i++) {
			deltaMasses[i] = psms.get(i).getMassDiff();
			targetsOnly[i] = psms.get(i).isTarget();
		}

		// ONLY the zero family is discovered by constrained KDE. Stage 1 selects genomes, so the rest of
		// the delta-mass axis is of no interest to it, and fitting the axis in
		// full is unnecessary. Expected isotope offsets are folded onto zero;
		// all that has to be determined is this family's empirical KDE basin.
		//
		// Peaks are detected from targets only: decoys contribute a broadly
		// uniform background that blurs peak edges. Assignment still covers every
		// PSM, so decoys keep their peak membership for the genome-level FDR that
		// consumes this output.
		DeltaMassPeakFitter.Result fit = DeltaMassPeakFitter.fitZeroPeakOnly(deltaMasses, targetsOnly, fitConfig,
				progress);

		if (fit.zeroPeakId < 0) {
			LOGGER.error("Stage 1: no zero peak could be fitted within +/-" + fitConfig.zeroWindow
					+ " Da, so no unmodified PSM set can be derived for genome discovery. "
					+ "Check the search ran open, and that mass calibration is not grossly off.");
			return new Result(new ArrayList<>(), Double.NaN, Double.NaN, Double.NaN, Double.NaN, "NONE",
					psms.size(), 0, 0, 0, 0, 0, Double.NaN, fdrThreshold);
		}

		DeltaMassPeakFitter.Peak zeroPeak = fit.peaks.get(fit.zeroPeakId);
		double legacyHalfWidth = Math.max(fitConfig.membershipSigmas * zeroPeak.sigma,
				fitConfig.membershipFloorDa);
		double zeroLower = Double.isFinite(zeroPeak.lowerBoundary) ? zeroPeak.lowerBoundary
				: zeroPeak.centroid - legacyHalfWidth;
		double zeroUpper = Double.isFinite(zeroPeak.upperBoundary) ? zeroPeak.upperBoundary
				: zeroPeak.centroid + legacyHalfWidth;

		List<OpenSearchPsm> inZeroPeak = new ArrayList<>();
		for (int i = 0; i < psms.size(); i++) {
			if (fit.familyAssignment[i] == fit.zeroPeakId) {
				OpenSearchPsm psm = psms.get(i);
				psm.setObservedPeakId(fit.observedAssignment[i]);
				psm.setPeakId(fit.familyAssignment[i]);
				psm.setNormalizedMassDiff(fit.normalizedMassDiff[i]);
				psm.setPeakIsotopeOffset(fit.isotopeOffset[i]);
				inZeroPeak.add(psm);
			}
		}

		int nearIsotope = DeltaMassPeakFitter.countNearIsotopeOffsets(deltaMasses, fitConfig);
		warnIfIsotopePeaksMissing(nearIsotope);

		// Target-decoy FDR on the zero-peak set, then hand the surviving target
		// peptides to genome selection. Decoys are kept in the list because the
		// genome-level FDR downstream needs them; only the q-value cut is applied.
		double[] scores = new double[inZeroPeak.size()];
		boolean[] isTarget = new boolean[inZeroPeak.size()];
		int targetsBeforeFdr = 0;
		for (int i = 0; i < inZeroPeak.size(); i++) {
			OpenSearchPsm psm = inZeroPeak.get(i);
			double expect = psm.getExpect();
			scores[i] = Double.isNaN(expect) ? psm.getHyperscore() : -Math.log10(Math.max(expect, 1e-300));
			isTarget[i] = psm.isTarget();
			if (isTarget[i]) {
				targetsBeforeFdr++;
			}
		}
		int decoysBeforeFdr = inZeroPeak.size() - targetsBeforeFdr;
		double[] q = TargetDecoyFdr.qValues(scores, isTarget);

		List<OpenSearchPsm> unmodified = new ArrayList<>();
		Set<String> acceptedTargetPeptides = new HashSet<>();
		int targetsAfterFdr = 0;
		int decoysAfterFdr = 0;
		for (int i = 0; i < inZeroPeak.size(); i++) {
			if (q[i] <= fdrThreshold) {
				inZeroPeak.get(i).setQValue(q[i]);
				unmodified.add(inZeroPeak.get(i));
				if (isTarget[i]) {
					targetsAfterFdr++;
					String peptide = inZeroPeak.get(i).getPeptide();
					if (peptide != null) {
						acceptedTargetPeptides.add(peptide);
					}
				} else {
					decoysAfterFdr++;
				}
			}
		}
		double estimatedFdr = targetsAfterFdr == 0 ? Double.NaN
				: (double) decoysAfterFdr / (double) targetsAfterFdr;

		LOGGER.info(String.format(
				"Stage 1 unmodified selection: %d/%d PSMs inside the zero family (%s centroid %.5f Da = calibration "
						+ "offset, empirical basin [%.5f, %.5f] Da, local width %.5f Da); before FDR: %d targets/%d decoys; "
						+ "after q<=%.3f: %d targets/%d decoys (estimated FDR=%.6f, estimator=decoys/targets). "
						+ "%d PSMs lie near a +1/+2 isotope offset; exact zero-family satellites are normalized into the root.",
				inZeroPeak.size(), psms.size(), zeroPeak.discoveryMethod, zeroPeak.centroid,
				zeroLower, zeroUpper, zeroPeak.sigma,
				targetsBeforeFdr, decoysBeforeFdr, fdrThreshold, targetsAfterFdr, decoysAfterFdr, estimatedFdr,
				nearIsotope));

		return new Result(unmodified, zeroPeak.centroid, zeroPeak.sigma, zeroLower,
				zeroUpper, zeroPeak.discoveryMethod, psms.size() - inZeroPeak.size(), nearIsotope,
				targetsBeforeFdr, decoysBeforeFdr, targetsAfterFdr, decoysAfterFdr, estimatedFdr, fdrThreshold,
				inZeroPeak, acceptedTargetPeptides);
	}

	/**
	 * Both marker-search modes should show +1/+2 offsets: the open search through
	 * its Da-wide window, the closed one through {@code isotope_error = 0/1/2}.
	 * Their complete absence means a closed window without isotope correction,
	 * which loses every isotope-triggered spectrum.
	 */
	static void warnIfIsotopePeaksMissing(int nearIsotopeCount) {
		if (nearIsotopeCount > 0) {
			return;
		}
		LOGGER.warn("Stage 1: no PSMs were found near the +1/+2 isotope offsets. A closed marker search needs "
				+ "isotope_error = 0/1/2 and an open one a Da-wide window; check the precursor settings in the "
				+ "marker search parameters.");
	}
}
