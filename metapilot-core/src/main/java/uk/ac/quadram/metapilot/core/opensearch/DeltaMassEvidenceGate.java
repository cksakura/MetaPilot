package uk.ac.quadram.metapilot.core.opensearch;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.function.Consumer;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;


/**
 * Rejects fitted peaks whose members show no fragment-level evidence of
 * carrying a mass shift, and folds their PSMs back onto the zero family.
 *
 * <h2>The test</h2>
 * Because the searches set {@code localize_delta_mass = 1}, MSFragger scores
 * every spectrum twice - once with the delta mass placed on fragment ions and
 * once without - and the difference reaches this workflow as
 * {@link OpenSearchPsm#getDeltaMassScoreGain()}. If the delta is a real
 * modification, shifted fragments exist and the gain is positive. If the
 * precursor mass was simply wrong, the peptide is unmodified, there are no
 * shifted fragments, and the gain is zero or negative.
 *
 * <h2>Why it is not a tuned threshold</h2>
 * Measured on 12mix (2026-09-14), the fraction of each peak's PSMs with gain
 * &lt;= 0, for peaks with at least 200 members, zero family excluded, sorted
 * ascending:
 * <pre>
 * 3 6 6 8 8 9 10 15 18 23 28 37 47 58 65 66 | 86 87 88 91 92 ... 99 100
 * </pre>
 * There is a twenty-point gap with nothing in it. Gates at 0.70, 0.80 and 0.85
 * reject exactly the same 42 peaks and 40,014 PSMs and keep exactly the same 16
 * peaks and 37,637 PSMs. The population is bimodal; the threshold only has to
 * land in the gap.
 * <p>
 * The rejected side is every peak in the +3.9..+4.3 and +4.8..+5.4 Da bands -
 * 94-96% no-gain, median gain negative, indistinguishable from the unmodified
 * population. The kept side is chemically coherent: pyroGlu (-17.0265, 8.3%
 * no-gain), Ca-2H (+37.9472, 5.9%), oxidation (+15.9948, 18.2%), Na-H
 * (+21.9799, 37.2%), carbamidomethyl (+57.0214, 58.4%).
 *
 * <h2>Why folding onto zero rather than demoting</h2>
 * {@link ModificationValidator} demotes a rejected peak's PSMs to the
 * identification-only pool, which is right for a peak that failed an FDR test -
 * its members might be anything. It is the wrong answer here, because the
 * fragment evidence positively identifies what these PSMs are: the unmodified
 * peptide, carrying a precursor mass that Crystal-C did not correct. Demoting
 * them would file a known unmodified peptide in the "no accepted mass
 * explanation" tier and keep it out of the quantified unmodified form. So this
 * gate reassigns rather than demotes: {@code normalizedMassDiff} goes to zero,
 * the peak id becomes the zero family's, and the reason is recorded on the PSM.
 * <p>
 * {@code massDiff} - the raw, observed difference - is never touched. It stays
 * as the record of what the instrument actually reported, exactly as it does
 * after isotope-family folding.
 *
 * <h2>Why peaks inside the exclude range are never folded</h2>
 * MSFragger does not search shifted fragment ions inside
 * {@code delta_mass_exclude_ranges} (-1.5 to +3.5 Da in the bundled open-search
 * parameters), so a peak there reports no gain whether or not the shift is
 * real. On PXD050781 this gate folded the deamidation peak (+0.9836 Da, 96.5%
 * no-gain of 254 scored) onto zero. A peak inside that range is therefore left to
 * {@link ModificationValidator}; its evidence is still measured and reported.
 * <p>
 * This replaced an exemption by a hand-picked list of seven masses (2026-11-08).
 * A catalogue such as UNIMOD cannot take its place: the +3.8..+4.2 and +4.8..+5.2
 * Da artefact bands are a comb of fitted peaks about 10 mDa apart, and even
 * without isotopic labels 5 of 73 such peaks on 12mix and 5 of 52 on faecal lie
 * within 8 mDa of a UNIMOD mass (Val-&gt;Cys, Formylasparagine, Trp-&gt;Kynurenin),
 * so a catalogue exemption would keep artefacts this gate exists to remove.
 * <p>
 * In the joint and cascade identifications ({@link Config#foldRejected} off) the
 * gate folds nothing and its rejections are read by no later step; only the
 * per-peak support it records reaches {@code modification_verdicts.tsv}.
 */
public final class DeltaMassEvidenceGate {

	private static final Logger LOGGER = LogManager.getLogger(DeltaMassEvidenceGate.class);

	/** Recorded on a PSM this gate moved; surfaces in {@code validated_psms.tsv}. */
	public static final String REASSIGNMENT_TAG = "PRECURSOR_ARTIFACT";

	/** {@code delta_mass_exclude_ranges = (-1.5,3.5)} in the bundled MSFragger open-search parameters. */
	public static final double DELTA_MASS_EXCLUDE_LOWER_DA = -1.5;
	public static final double DELTA_MASS_EXCLUDE_UPPER_DA = 3.5;

	public static final class Config {
		/**
		 * Reject a peak when at least this share of its members gain nothing from
		 * placing the delta mass on fragment ions. Anywhere in 0.70-0.85 gives the
		 * same answer on 12mix; see the class docs.
		 */
		public double maxNoGainFraction = 0.80;
		/**
		 * Below this many scored members the fraction is too noisy to act on, and
		 * the peak is left to the normal FDR machinery.
		 */
		public int minScoredPsms = 50;
		/** A gain at or below this counts as "no gain"; guards float noise. */
		public double gainEpsilon = 0.001;
		/** Peaks this close to zero are the zero family and are never tested. */
		public double zeroWindowDa = 0.05;
		/**
		 * Lower bound, Da, of the search engine's {@code delta_mass_exclude_ranges}: a peak whose
		 * centroid lies inside it cannot show shifted-fragment gain and is never rejected or folded.
		 * Must equal the value in the bundled {@code msfragger/open_fragger*.params}, which a test checks.
		 */
		public double exemptRangeLowerDa = DELTA_MASS_EXCLUDE_LOWER_DA;
		/** Upper bound, Da, of that range. */
		public double exemptRangeUpperDa = DELTA_MASS_EXCLUDE_UPPER_DA;
		/** Off restores the previous behaviour exactly. */
		public boolean enabled = true;
		/**
		 * Move a rejected peak's PSMs onto the zero family. Off with the closed stage: the peak is
		 * still rejected as a modification, its PSMs keep their delta and are unexplained shifts,
		 * and a peptide that is really unmodified is found by the closed search at its mass.
		 */
		public boolean foldRejected = true;
	}

	public static final class PeakEvidence {
		public final int peakId;
		public final double centroid;
		public final int scored;
		public final double noGainFraction;
		public final boolean exempt;
		public final boolean rejected;

		PeakEvidence(int peakId, double centroid, int scored, double noGainFraction, boolean exempt,
				boolean rejected) {
			this.peakId = peakId;
			this.centroid = centroid;
			this.scored = scored;
			this.noGainFraction = noGainFraction;
			this.exempt = exempt;
			this.rejected = rejected;
		}
	}

	public static final class Result {
		public final List<PeakEvidence> evidence;
		public final Set<Integer> rejectedPeakIds;
		/**
		 * Peaks with at least {@link Config#minScoredPsms} scored members, fewer
		 * than {@link Config#maxNoGainFraction} of which gain nothing.
		 */
		public final Set<Integer> supportedPeakIds;
		public final int psmsReassigned;

		Result(List<PeakEvidence> evidence, Set<Integer> rejectedPeakIds, Set<Integer> supportedPeakIds,
				int psmsReassigned) {
			this.evidence = evidence;
			this.rejectedPeakIds = rejectedPeakIds;
			this.supportedPeakIds = supportedPeakIds;
			this.psmsReassigned = psmsReassigned;
		}
	}

	private DeltaMassEvidenceGate() {
	}

	/**
	 * Tests every non-zero peak and folds the artifact ones onto the zero family.
	 * <p>
	 * Must run after {@link OpenSearchFeatures#applyPeakContext} - it needs peak
	 * assignments - and before the classifier, so that the reassigned PSMs carry
	 * a zero delta into the peak-context features and the modification verdict.
	 *
	 * @param psms      every candidate, mutated in place for the reassigned ones
	 * @param fit       the fitted peaks, for centroids and the zero family's id
	 * @param cfg       thresholds
	 * @param progress  optional sink
	 */
	public static Result apply(List<OpenSearchPsm> psms, DeltaMassPeakFitter.Result fit, Config cfg,
			Consumer<String> progress) {
		if (cfg == null) {
			cfg = new Config();
		}
		if (!cfg.enabled) {
			return new Result(new ArrayList<>(), new HashSet<>(), new HashSet<>(), 0);
		}
		Map<Integer, int[]> counts = new HashMap<>();     // peakId -> {scored, noGain}
		Map<Integer, Double> centroid = new HashMap<>();
		int zeroPeakId = -1;
		double bestZero = Double.MAX_VALUE;

		for (OpenSearchPsm psm : psms) {
			int pid = psm.getPeakId();
			if (pid < 0) {
				continue;
			}
			double c = psm.getPeakCentroid();
			centroid.putIfAbsent(pid, c);
			if (Math.abs(c) < cfg.zeroWindowDa && Math.abs(c) < bestZero) {
				bestZero = Math.abs(c);
				zeroPeakId = pid;
			}
			if (Math.abs(c) < cfg.zeroWindowDa) {
				continue;
			}
			double gain = psm.getDeltaMassScoreGain();
			if (Double.isNaN(gain)) {
				continue;
			}
			int[] a = counts.computeIfAbsent(pid, k -> new int[2]);
			a[0]++;
			if (gain <= cfg.gainEpsilon) {
				a[1]++;
			}
		}

		List<PeakEvidence> evidence = new ArrayList<>();
		Set<Integer> rejected = new HashSet<>();
		Set<Integer> supported = new HashSet<>();
		int exemptedFromRejection = 0;
		for (Map.Entry<Integer, int[]> e : counts.entrySet()) {
			int scored = e.getValue()[0];
			double frac = scored == 0 ? 0.0 : (double) e.getValue()[1] / scored;
			double c = centroid.getOrDefault(e.getKey(), Double.NaN);
			boolean exempt = c >= cfg.exemptRangeLowerDa && c <= cfg.exemptRangeUpperDa;
			boolean enoughScored = scored >= cfg.minScoredPsms;
			boolean noSupport = enoughScored && frac >= cfg.maxNoGainFraction;
			boolean rej = noSupport && !exempt;
			evidence.add(new PeakEvidence(e.getKey(), c, scored, frac, exempt, rej));
			if (rej) {
				rejected.add(e.getKey());
			}
			if (noSupport && exempt) {
				exemptedFromRejection++;
			}
			if (enoughScored && frac < cfg.maxNoGainFraction) {
				supported.add(e.getKey());
			}
		}
		evidence.sort((a, b) -> Double.compare(b.noGainFraction, a.noGainFraction));

		if (!cfg.foldRejected) {
			int kept = 0;
			for (OpenSearchPsm psm : psms) {
				if (rejected.contains(psm.getPeakId())) {
					kept++;
				}
			}
			if (progress != null) {
				progress.accept(String.format(Locale.ROOT,
						"Delta-mass evidence gate: %d of %d non-zero peaks show no shifted-fragment support "
								+ "(>=%.0f%% of members gain nothing) and are rejected; their %,d PSMs keep "
								+ "their delta as unexplained shifts; %d peak(s) inside the engine exclude range without support kept",
						rejected.size(), counts.size(), 100 * cfg.maxNoGainFraction, kept, exemptedFromRejection));
			}
			return new Result(evidence, rejected, supported, 0);
		}

		if (zeroPeakId < 0) {
			String warning = "Delta-mass evidence gate found no zero-family peak within +/-" + cfg.zeroWindowDa
					+ " Da, so there is nowhere to fold precursor artifacts onto. Leaving every peak as fitted; "
					+ rejected.size() + " peak(s) would otherwise have been rejected.";
			LOGGER.warn(warning);
			if (progress != null) {
				progress.accept("WARNING: " + warning);
			}
			return new Result(evidence, new HashSet<>(), supported, 0);
		}

		int moved = 0;
		double zeroCentroid = centroid.getOrDefault(zeroPeakId, 0.0);
		for (OpenSearchPsm psm : psms) {
			if (!rejected.contains(psm.getPeakId())) {
				continue;
			}
			psm.setPeakId(zeroPeakId);
			psm.setPeakCentroid(zeroCentroid);
			psm.setPeakIsotopeOffset(0);
			// The observed peak is left alone deliberately: it records which fitted
			// mode this PSM actually landed in, which is what makes the
			// reassignment auditable after the fact.
			psm.setNormalizedMassDiff(0.0);
			// The peak-context features describe distance to the peak this PSM was
			// in. It is now in the zero family with a zero delta, so leaving the old
			// distances would hand the classifier a feature that contradicts the
			// assignment it is scoring.
			psm.setDistanceToPeakCentroid(0.0);
			psm.setDistanceInSigmas(0.0);
			psm.setDeltaMassReassignment(REASSIGNMENT_TAG);
			moved++;
		}

		if (progress != null) {
			progress.accept(String.format(Locale.ROOT,
					"Delta-mass evidence gate: %d of %d non-zero peaks show no shifted-fragment support "
							+ "(>=%.0f%% of members gain nothing); %,d PSMs folded onto the zero family as "
							+ "precursor artifacts; %d peak(s) inside the engine exclude range without support kept as fitted",
					rejected.size(), counts.size(), 100 * cfg.maxNoGainFraction, moved, exemptedFromRejection));
			int shown = 0;
			for (PeakEvidence pe : evidence) {
				if (!pe.rejected || shown++ >= 8) {
					continue;
				}
				progress.accept(String.format(Locale.ROOT, "    rejected peak %d at %+.4f Da: %.1f%% no-gain of %d",
						pe.peakId, pe.centroid, 100 * pe.noGainFraction, pe.scored));
			}
		}
		LOGGER.info("Delta-mass evidence gate rejected {} peaks and reassigned {} PSMs to the zero family; "
				+ "{} peaks inside the engine exclude range exempted", rejected.size(), moved, exemptedFromRejection);
		return new Result(evidence, rejected, supported, moved);
	}
}
