package uk.ac.quadram.metapilot.core.opensearch;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.EnumMap;
import java.util.HashMap;
import java.util.HashSet;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.function.Consumer;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import uk.ac.quadram.metapilot.core.variant.MzmlScanReader;

/**
 * Joint closed+open identification: one candidate pool per spectrum from both searches, one
 * semi-supervised score over both searches' evidence, one winner per spectrum, and target-decoy
 * FDR controlled separately within each evidence group, at PSM and at peptide level.
 *
 * <h2>Why not a cascade</h2>
 * A cascade hands each spectrum to the closed search first: a strong open candidate can never
 * overturn a weak closed answer, a confident closed answer never informs the open ranking, and the
 * two stages score on two scales that later steps have to reconcile. Here every candidate of a
 * spectrum, whichever search proposed it, competes on one score, and the score sees what the other
 * search found (presence, the best exact-mass rival, the hyperscore margin) and, with DIA-NN, how
 * well each candidate's predicted spectrum and retention time match.
 *
 * <h2>FDR by group</h2>
 * A winner's {@link JointCandidate.Group} is decided by the evidence it carries before any score
 * exists, identically for targets and decoys, so a target-decoy walk inside each group estimates
 * that group's own false discovery rate. Pooling the groups would let the large, clean exact-mass
 * group carry open-only identifications at a much higher error rate under one overall 1%: on the
 * 12mix cascade, open-only peptides ran at 4.66% decoy-estimated FDR inside a class that read 1%.
 * Each group at 1% keeps the union at 1% or below.
 * <p>
 * Peptides likewise: each peptide belongs to the highest-priority group in which it has an accepted
 * PSM (exact before open unmodified before open accepted shift), and is walked in that group alone,
 * scored by its best accepted PSM there. A peptide confirmed at exact mass anywhere is an exact
 * peptide; a peptide only the open search identifies has to pass the open group's own walk - which
 * is what makes the reported open-only peptides confident.
 * <p>
 * The {@link JointCandidate.Group#OPEN_UNEXPLAINED} group is walked and reported but its PSMs are
 * not identifications unless {@link Config#acceptUnexplained} is set: decoys do not model a
 * relative's peptide absorbing a sequence difference as a mass shift.
 */
public final class JointIdentification {

	private static final Logger LOGGER = LogManager.getLogger(JointIdentification.class);

	/** 13C - 12C. */
	static final double ISOTOPE_SPACING = 1.0033548;

	public static final class Config {
		/** PSM- and peptide-level FDR, per group. */
		public double fdr = 0.01;
		/** Treat the unexplained-shift group as identifications. Off: reported only. */
		public boolean acceptUnexplained = false;
		/** Feature names never given to the model. */
		public Set<String> withheldFeatures = new HashSet<>();
		public SemiSupervisedPsmClassifier.Config classifier = new SemiSupervisedPsmClassifier.Config();
		/**
		 * Identification v2, C1: one query per (spectrum, precursor); closed and open rows merge only
		 * for the same precursor, and FDR groups split by precursor class (same, reassigned).
		 */
		public boolean precursorQueries = false;
		/**
		 * Identification v2, C2: RESET. Decoy peptides are split once into training and estimating
		 * decoys; a query won by a training decoy reports nothing; every estimate counts each estimating
		 * decoy {@link #decoyWeight()} times. A tied query goes to the decoy.
		 */
		public boolean reset = false;
		/**
		 * Shift localisation (identification v2): unexplained-shift spectrum winners whose shift the
		 * fragments place on one residue, where a substitution or UNIMOD entry explains it, move to
		 * {@link JointCandidate.Group#OPEN_LOCALISED}, an FDR group of its own. Null: off.
		 */
		public ShiftLocaliser shiftLocaliser;
		/** Spectra by sample for {@link #shiftLocaliser}. */
		public Map<String, File> localisationSpectra = new java.util.LinkedHashMap<>();
		/** {@link PeakFragmentTest}: largest chance share at which an established peak's shift is on the fragments. */
		public double peakTestMaxChanceShare = PeakFragmentTest.DEFAULT_MAX_CHANCE_SHARE;
		/** {@link PeakFragmentTest}: unmodified zero-peak winners scored at each peak's mass. */
		public int peakTestNullSample = PeakFragmentTest.DEFAULT_NULL_SAMPLE;

		public double decoyWeight() {
			return reset ? 2.0 : 1.0;
		}
	}

	public static final class Result {
		public final List<JointCandidate> candidates;
		public final JointFeatures features;
		public final double[] scores;
		public final double[] logits;
		public final boolean[] winner;
		/** Group-wise PSM q-value of each winner; NaN for other rows. */
		public final double[] psmQ;
		/** Group-wise peptide q-value of each accepted winner's peptide; NaN otherwise. */
		public final double[] peptideQ;
		/** Winner in an identification group at {@code psmQ <= fdr}. */
		public final boolean[] accepted;
		/** Accepted, and its peptide passed its group's peptide walk: the identifications. */
		public final boolean[] kept;
		public final ModificationValidator.Result validation;
		public final List<String> summary;
		/** RESET training decoys (never reported); null without RESET. */
		public final boolean[] trainingDecoy;
		/** Every decoy-based estimate counts each (estimating) decoy this many times. */
		public final double decoyWeight;
		/** {@link PeakFragmentTest} of every established peak with localised winners; empty when localisation is off. */
		public List<PeakFragmentTest.PeakResult> peakFragmentTests = List.of();

		Result(List<JointCandidate> candidates, JointFeatures features, double[] scores, double[] logits,
				boolean[] winner, double[] psmQ, double[] peptideQ, boolean[] accepted, boolean[] kept,
				ModificationValidator.Result validation, List<String> summary, boolean[] trainingDecoy,
				double decoyWeight) {
			this.trainingDecoy = trainingDecoy;
			this.decoyWeight = decoyWeight;
			this.candidates = candidates;
			this.features = features;
			this.scores = scores;
			this.logits = logits;
			this.winner = winner;
			this.psmQ = psmQ;
			this.peptideQ = peptideQ;
			this.accepted = accepted;
			this.kept = kept;
			this.validation = validation;
			this.summary = summary;
		}
	}

	private JointIdentification() {
	}

	/**
	 * @param open      every open candidate (any rank, either engine), peak context and the peak
	 *                  filter already applied ({@code isInValidatedPeak} set from {@code gate})
	 * @param closed    the closed search's candidates (any rank)
	 * @param gate      the open search's peak decisions
	 * @param ms2Inputs DIA-NN and spectra for the predicted-spectrum features, or null for none
	 */
	public static Result run(List<OpenSearchPsm> open, List<OpenSearchPsm> closed, ModificationValidator.Gate gate,
			Config cfg, PredictedSpectrumFeatures.Inputs ms2Inputs, Consumer<String> progress, Consumer<String> warn)
			throws Exception {
		if (cfg == null) {
			cfg = new Config();
		}
		Consumer<String> info = message -> report(progress, message);
		Consumer<String> warning = warn == null ? info : warn;

		prepareClosed(closed);
		int[] dropped = new int[2];
		List<JointCandidate> candidates = JointCandidate.merge(open, closed, dropped, cfg.precursorQueries);
		if (dropped[0] > 0) {
			warning.accept("Joint identification: " + dropped[0] + " rows had no peptide or file name and were not merged");
		}
		if (dropped[1] > 0) {
			warning.accept("Joint identification: " + dropped[1] + " rows disagree with another search on whether the "
					+ "same peptide is a target or a decoy; the first row's label was kept");
		}
		Map<Integer, ModificationValidator.MassClass> massClasses = ModificationValidator.massClasses(gate);
		Set<Integer> zeroPeaks = new HashSet<>();
		for (Map.Entry<Integer, ModificationValidator.MassClass> e : massClasses.entrySet()) {
			if (e.getValue() == ModificationValidator.MassClass.ZERO && gate.acceptedPeakIds.contains(e.getKey())) {
				zeroPeaks.add(e.getKey());
			}
		}
		EnumMap<JointCandidate.Group, Integer> groupSizes = new EnumMap<>(JointCandidate.Group.class);
		int both = 0;
		for (JointCandidate c : candidates) {
			c.group = groupOf(c, gate, massClasses);
			groupSizes.merge(c.group, 1, Integer::sum);
			if (c.closed != null && c.foundByOpen()) {
				both++;
			}
		}
		info.accept(String.format(Locale.ROOT,
				"Joint identification: %,d open + %,d closed rows -> %,d candidates (%,d proposed by both searches); "
						+ "candidates by group %s",
				open.size(), closed.size(), candidates.size(), both, groupSizes));

		JointFeatures.SpectrumScores ms2 = ms2Inputs == null ? null
				: PredictedSpectrumFeatures.compute(candidates, ms2Inputs, info, warning);
		if (cfg.withheldFeatures != null) {
			Set<String> unknown = new java.util.TreeSet<>(cfg.withheldFeatures);
			unknown.removeAll(Arrays.asList(JointFeatures.ALL_FEATURES));
			if (!unknown.isEmpty()) {
				warning.accept("Joint identification: joint_withheld_features names no such feature: " + unknown
						+ "; known features: " + String.join(",", JointFeatures.ALL_FEATURES));
			}
		}
		JointFeatures features = JointFeatures.of(candidates, ms2, zeroPeaks, cfg.withheldFeatures);
		boolean[] trainingDecoy = cfg.reset ? resetRoles(candidates) : null;
		if (trainingDecoy != null) {
			int t = 0;
			for (boolean b : trainingDecoy) {
				t += b ? 1 : 0;
			}
			info.accept(String.format(Locale.ROOT, "Joint identification: RESET - %,d training-decoy candidates, "
					+ "decoy weight %.0f", t, cfg.decoyWeight()));
		}
		JointPsmClassifier.Result scored = JointPsmClassifier.classify(features, candidates, cfg.classifier, progress,
				trainingDecoy);

		int n = candidates.size();
		for (int i = 0; i < n; i++) {
			JointCandidate c = candidates.get(i);
			String support = c.support();
			for (OpenSearchPsm psm : new OpenSearchPsm[] { c.closed, c.fragger, c.sage }) {
				if (psm != null) {
					psm.setClassScore(scored.scores[i]);
					psm.setClassLogit(scored.logits[i]);
					psm.setEngineSupport(support);
					psm.setIdentificationGroup(c.group.name());
				}
			}
		}

		boolean[] winner = cfg.reset ? winners(features.spectrumIds(), scored.logits, candidates)
				: winners(features.spectrumIds(), scored.logits);
		if (trainingDecoy != null) {
			// A query won by a training decoy reports nothing.
			for (int i = 0; i < winner.length; i++) {
				winner[i] &= !trainingDecoy[i];
			}
		}
		Localisation localisation = null;
		if (cfg.shiftLocaliser != null) {
			localisation = localise(candidates, winner, cfg, PeakContext.of(open, gate, massClasses), info, warning);
		}
		// The peak test first, on every peak's winners as the peak was accepted: it counts the winners
		// accepted with the whole peak as one group, exactly as before the per-PSM rule existed, so that
		// rule cannot change which peaks pass. Then the per-PSM rule (Localisation#decideAcceptedShifts)
		// takes winners out of the accepted-shift group, and only then come the q-values that count.
		List<PeakFragmentTest.PeakResult> peakTests = List.of();
		if (localisation != null) {
			double[] peakQ = groupQValues(candidates, scored.logits, winner, cfg.precursorQueries, cfg.decoyWeight());
			boolean[] peakAccepted = new boolean[n];
			for (int i = 0; i < n; i++) {
				peakAccepted[i] = winner[i] && peakQ[i] <= cfg.fdr && identifies(candidates.get(i).group, cfg);
			}
			peakTests = localisation.testPeaks(peakAccepted, info);
			localisation.decideAcceptedShifts(info);
			localisation.apply();
		}
		double[] psmQ = groupQValues(candidates, scored.logits, winner, cfg.precursorQueries, cfg.decoyWeight());
		boolean[] accepted = new boolean[n];
		for (int i = 0; i < n; i++) {
			accepted[i] = winner[i] && psmQ[i] <= cfg.fdr && identifies(candidates.get(i).group, cfg);
		}
		double[] peptideQ = new double[n];
		Arrays.fill(peptideQ, Double.NaN);
		PeptideWalk walk = peptideWalk(candidates, scored.logits, accepted, cfg.fdr, cfg.decoyWeight());
		boolean[] kept = new boolean[n];
		for (int i = 0; i < n; i++) {
			if (!accepted[i]) {
				continue;
			}
			Double q = walk.unitQ.get(unitKey(candidates.get(i)));
			peptideQ[i] = q == null ? Double.NaN : q;
			kept[i] = q != null && q <= cfg.fdr;
		}

		ModificationValidator.Result validation = toValidation(candidates, winner, psmQ, accepted, kept, gate, cfg);
		List<String> summary = summarise(candidates, winner, psmQ, accepted, kept, walk, cfg, both, features, scored);
		for (String line : summary) {
			LOGGER.info("Joint identification: " + line.replace('\t', ' '));
		}
		Result result = new Result(candidates, features, scored.scores, scored.logits, winner, psmQ, peptideQ,
				accepted, kept, validation, summary, trainingDecoy, cfg.decoyWeight());
		result.peakFragmentTests = peakTests;
		return result;
	}

	/**
	 * Isotope offset ({@code 0..2}, the closed search's {@code isotope_error}) and the mass error left
	 * after it. A closed row has no delta-mass peak.
	 */
	static void prepareClosed(List<OpenSearchPsm> closed) {
		for (OpenSearchPsm psm : closed) {
			double massDiff = psm.getMassDiff();
			int k = Double.isFinite(massDiff) ? (int) Math.round(massDiff / ISOTOPE_SPACING) : 0;
			k = Math.max(0, Math.min(2, k));
			psm.setPeakIsotopeOffset(k);
			psm.setNormalizedMassDiff(Double.isFinite(massDiff) ? massDiff - k * ISOTOPE_SPACING : Double.NaN);
			psm.setPeakId(-1);
		}
	}

	static JointCandidate.Group groupOf(JointCandidate c, ModificationValidator.Gate gate,
			Map<Integer, ModificationValidator.MassClass> massClasses) {
		if (c.closed != null) {
			return JointCandidate.Group.EXACT;
		}
		OpenSearchPsm open = c.openPrimary();
		int peak = open.getPeakId();
		if (peak >= 0 && gate.acceptedPeakIds.contains(peak)) {
			ModificationValidator.MassClass mc = massClasses.get(peak);
			if (mc == ModificationValidator.MassClass.ZERO) {
				return JointCandidate.Group.OPEN_UNMODIFIED;
			}
			if (mc == ModificationValidator.MassClass.ESTABLISHED) {
				return JointCandidate.Group.OPEN_ACCEPTED_SHIFT;
			}
		}
		return JointCandidate.Group.OPEN_UNEXPLAINED;
	}

	/**
	 * The fitted peaks as the localisation step needs them: each peak family's centroid, the accepted
	 * zero and established families, and every observed peak's mass range, family and isotope offset, so
	 * a closed PSM's mass difference lands in the peak an open PSM of that mass is in.
	 */
	static final class PeakContext {
		final Map<Integer, Double> centroids = new HashMap<>();
		final Set<Integer> zeroFamilies = new HashSet<>();
		final Set<Integer> established = new HashSet<>();
		private double[] lo = new double[0];
		private double[] hi = new double[0];
		private double[] mean = new double[0];
		private int[] family = new int[0];
		private int[] isotope = new int[0];

		/** Where a mass falls: a peak family, and the isotope offset of the observed peak within it. */
		static final class Place {
			final int family;
			final int isotopeOffset;

			Place(int family, int isotopeOffset) {
				this.family = family;
				this.isotopeOffset = isotopeOffset;
			}
		}

		static PeakContext of(List<OpenSearchPsm> open, ModificationValidator.Gate gate,
				Map<Integer, ModificationValidator.MassClass> massClasses) {
			PeakContext ctx = new PeakContext();
			for (ModificationValidator.PeakVerdict v : gate.verdicts) {
				ctx.centroids.put(v.peakId, v.centroid);
				if (!gate.acceptedPeakIds.contains(v.peakId)) {
					continue;
				}
				ModificationValidator.MassClass mc = massClasses.get(v.peakId);
				if (mc == ModificationValidator.MassClass.ZERO) {
					ctx.zeroFamilies.add(v.peakId);
				} else if (mc == ModificationValidator.MassClass.ESTABLISHED) {
					ctx.established.add(v.peakId);
				}
			}
			// observed peak -> {min, max, sum, n}, and votes for its (family, isotope offset)
			Map<Integer, double[]> range = new java.util.TreeMap<>();
			Map<Integer, Map<List<Integer>, Integer>> votes = new HashMap<>();
			for (OpenSearchPsm psm : open) {
				int obs = psm.getObservedPeakId();
				double d = psm.getMassDiff();
				if (obs < 0 || !Double.isFinite(d)) {
					continue;
				}
				double[] r = range.computeIfAbsent(obs,
						k -> new double[] { Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, 0, 0 });
				r[0] = Math.min(r[0], d);
				r[1] = Math.max(r[1], d);
				r[2] += d;
				r[3]++;
				votes.computeIfAbsent(obs, k -> new HashMap<>())
						.merge(List.of(psm.getPeakId(), psm.getPeakIsotopeOffset()), 1, Integer::sum);
			}
			int n = range.size();
			ctx.lo = new double[n];
			ctx.hi = new double[n];
			ctx.mean = new double[n];
			ctx.family = new int[n];
			ctx.isotope = new int[n];
			int i = 0;
			for (Map.Entry<Integer, double[]> e : range.entrySet()) {
				double[] r = e.getValue();
				ctx.lo[i] = r[0];
				ctx.hi[i] = r[1];
				ctx.mean[i] = r[2] / r[3];
				List<Integer> best = null;
				int bestCount = -1;
				for (Map.Entry<List<Integer>, Integer> v : votes.get(e.getKey()).entrySet()) {
					List<Integer> k = v.getKey();
					if (v.getValue() > bestCount || v.getValue() == bestCount
							&& (k.get(0) < best.get(0) || k.get(0).equals(best.get(0)) && k.get(1) < best.get(1))) {
						best = k;
						bestCount = v.getValue();
					}
				}
				ctx.family[i] = best.get(0);
				ctx.isotope[i] = best.get(1);
				i++;
			}
			return ctx;
		}

		/** The observed peak whose mass range holds {@code massDiff} (nearest mean on overlap); null outside all. */
		Place locate(double massDiff) {
			int best = -1;
			for (int i = 0; i < lo.length; i++) {
				if (massDiff >= lo[i] && massDiff <= hi[i]
						&& (best < 0 || Math.abs(mean[i] - massDiff) < Math.abs(mean[best] - massDiff))) {
					best = i;
				}
			}
			return best < 0 || family[best] < 0 ? null : new Place(family[best], isotope[best]);
		}
	}

	/** A winner whose shift is localised: the PSM scored, at which mass, and the peak family it belongs to. */
	private static final class ShiftItem {
		/** Row of {@link #candidate} in the candidate list. */
		final int index;
		final JointCandidate candidate;
		final OpenSearchPsm psm;
		final double delta;
		/** Peak family the peak test judges; -1 for unexplained shifts, which have none. */
		final int family;

		ShiftItem(int index, JointCandidate candidate, OpenSearchPsm psm, double delta, int family) {
			this.index = index;
			this.candidate = candidate;
			this.psm = psm;
			this.delta = delta;
			this.family = family;
		}
	}

	/** An unmodified zero-peak winner of the per-peak test's null pool. */
	private static final class NullPsm {
		final int index;
		final String file;
		final int scan;
		final String peptide;
		final int charge;
		MzmlScanReader.Peaks peaks;

		NullPsm(int index, String file, int scan, String peptide, int charge) {
			this.index = index;
			this.file = file;
			this.scan = scan;
			this.peptide = peptide;
			this.charge = charge;
		}
	}

	/** As {@link #localiseWinners(List, boolean[], Config, PeakContext, Consumer, Consumer)}, without peaks. */
	static List<PeakFragmentTest.PeakResult> localiseWinners(List<JointCandidate> candidates, boolean[] winner,
			Config cfg, Consumer<String> info, Consumer<String> warning) {
		return localiseWinners(candidates, winner, cfg, null, info, warning);
	}

	/**
	 * {@link #localise}, then the peak test counting every winner as accepted, then the labels: in one
	 * step, for callers without q-values. {@link #run} tests on the accepted winners instead.
	 */
	static List<PeakFragmentTest.PeakResult> localiseWinners(List<JointCandidate> candidates, boolean[] winner,
			Config cfg, PeakContext peaks, Consumer<String> info, Consumer<String> warning) {
		Localisation l = localise(candidates, winner, cfg, peaks, info, warning);
		List<PeakFragmentTest.PeakResult> tests = peaks == null ? List.of() : l.testPeaks(null, info);
		l.decideAcceptedShifts(info);
		l.apply();
		return tests;
	}

	/**
	 * Localises every unexplained-shift spectrum winner, target or decoy alike, and moves the
	 * localised-explained ones to {@link JointCandidate.Group#OPEN_LOCALISED}, or to
	 * {@link JointCandidate.Group#OPEN_SUBSTITUTION_SITE} where a substitution could produce the
	 * shift at the site. Runs after the winners
	 * are chosen: the model and the competition are unchanged, only the FDR group and the reported
	 * explanation are.
	 * <p>
	 * Accepted-shift winners are checked too ({@link ShiftLocaliser#checkAcceptedShift}, design doc
	 * section 8), on the folded shift. One explained on a residue a substitution fits moves to
	 * {@code OPEN_SUBSTITUTION_SITE} like an off-peak winner; every other stays
	 * {@code OPEN_ACCEPTED_SHIFT}, and the result decides only what is exported: the unmodified peptide
	 * without fragment support, the modification on an explained residue, or the shift on its own site
	 * when the fragments contradict every explained residue.
	 * <p>
	 * One rule for every shift, whichever search found it. An exact-mass winner that the closed search
	 * accepted with an isotope offset and no variable modification is checked the same way when its
	 * mass difference lies in a fitted peak other than an isotope peak of the zero family, accepted as an
	 * identification group or not (deamidation is withdrawn as substitution-shaped); its group stays
	 * {@code EXACT}. On 12mix, 18,849 such PSMs sat at 0.984 Da (deamidation, 19 mDa from
	 * the +1 isotope step, inside a 20 ppm precursor tolerance) and their fragments carried +0.984 on
	 * N/Q in 79% of them against 3% for unmodified peptides.
	 * <p>
	 * Nothing is labelled yet: {@link Localisation#testPeaks} judges each peak as a whole
	 * ({@link PeakFragmentTest}) once q-values exist, and {@link Localisation#apply} then labels the PSMs.
	 * Groups are set here, as before, because the q-value walks need them.
	 */
	static Localisation localise(List<JointCandidate> candidates, boolean[] winner, Config cfg, PeakContext peaks,
			Consumer<String> info, Consumer<String> warning) {
		Map<String, List<ShiftItem>> byFile = new java.util.TreeMap<>();
		List<NullPsm> nullPool = new ArrayList<>();
		int closedIsotope = 0;
		int closedInIsotopePeak = 0;
		int closedChecked = 0;
		for (int i = 0; i < candidates.size(); i++) {
			if (!winner[i]) {
				continue;
			}
			JointCandidate c = candidates.get(i);
			ShiftItem item = null;
			if (c.group == JointCandidate.Group.OPEN_UNEXPLAINED) {
				OpenSearchPsm open = c.openPrimary();
				item = new ShiftItem(i, c, open, open.getMassDiff(), -1);
			} else if (c.group == JointCandidate.Group.OPEN_ACCEPTED_SHIFT) {
				OpenSearchPsm open = c.openPrimary();
				item = new ShiftItem(i, c, open, open.getNormalizedMassDiff(), open.getPeakId());
			} else if (peaks != null && c.group == JointCandidate.Group.EXACT && c.closed != null
					&& c.closed.getPeakIsotopeOffset() >= 1 && c.closed.getClosedModificationCount() == 0) {
				closedIsotope++;
				PeakContext.Place place = peaks.locate(c.closed.getMassDiff());
				if (place != null && peaks.zeroFamilies.contains(place.family)) {
					closedInIsotopePeak++;
				} else if (place != null && peaks.centroids.containsKey(place.family)) {
					// Any other fitted peak, accepted as an identification group or not: whether the PSM
					// is identified is settled (it is an exact-mass winner); only its label is open.
					closedChecked++;
					item = new ShiftItem(i, c, c.closed,
							c.closed.getMassDiff() - place.isotopeOffset * ISOTOPE_SPACING, place.family);
				}
			}
			if (item != null) {
				byFile.computeIfAbsent(c.getFileName(), k -> new ArrayList<>()).add(item);
			}
			if (peaks != null) {
				NullPsm z = nullCandidate(i, c);
				if (z != null) {
					nullPool.add(z);
				}
			}
		}
		// The null pool: a fixed, label-blind subset of unmodified zero-peak winners, NULL_POOL_FACTOR times
		// the null sample, so that enough of them are accepted once q-values exist (Localisation#testPeaks).
		nullPool.sort((a, b) -> {
			int h = Long.compare(SemiSupervisedPsmClassifier.stableHash(a.file + ":" + a.scan),
					SemiSupervisedPsmClassifier.stableHash(b.file + ":" + b.scan));
			return h != 0 ? h : a.file.equals(b.file) ? Integer.compare(a.scan, b.scan) : a.file.compareTo(b.file);
		});
		List<NullPsm> nullSample = nullPool.subList(0,
				(int) Math.min((long) cfg.peakTestNullSample * NULL_POOL_FACTOR, nullPool.size()));
		Map<String, List<NullPsm>> nullByFile = new HashMap<>();
		for (NullPsm z : nullSample) {
			nullByFile.computeIfAbsent(z.file, k -> new ArrayList<>()).add(z);
		}
		Set<String> files = new java.util.TreeSet<>(byFile.keySet());
		files.addAll(nullByFile.keySet());

		EnumMap<ShiftLocaliser.State, Integer> states = new EnumMap<>(ShiftLocaliser.State.class);
		EnumMap<ShiftLocaliser.State, Integer> acceptedShift = new EnumMap<>(ShiftLocaliser.State.class);
		int noSpectra = 0;
		for (String file : files) {
			List<ShiftItem> items = byFile.getOrDefault(file, List.of());
			List<NullPsm> nulls = nullByFile.getOrDefault(file, List.of());
			File mzml = cfg.localisationSpectra.get(file);
			if (mzml == null || !mzml.isFile()) {
				noSpectra += items.size();
				continue;
			}
			Set<Integer> scans = new HashSet<>();
			for (ShiftItem it : items) {
				scans.add(it.candidate.getScan());
			}
			for (NullPsm z : nulls) {
				scans.add(z.scan);
			}
			MzmlScanReader.Result spectra;
			try {
				spectra = MzmlScanReader.read(mzml, scans);
			} catch (IOException ex) {
				warning.accept("Joint identification: shift localisation could not read " + mzml + ": "
						+ ex.getMessage() + "; its unexplained-shift winners stay unexplained");
				noSpectra += items.size();
				continue;
			}
			for (NullPsm z : nulls) {
				z.peaks = spectra.scans.get(z.scan);
			}
			for (ShiftItem it : items) {
				JointCandidate c = it.candidate;
				MzmlScanReader.Peaks scan = spectra.scans.get(c.getScan());
				if (it.family >= 0) {
					c.localisation = cfg.shiftLocaliser.checkAcceptedShift(c.getPeptide(), it.psm.getCharge(), it.delta,
							scan);
					acceptedShift.merge(c.localisation.state, 1, Integer::sum);
				} else {
					c.localisation = cfg.shiftLocaliser.localise(c.getPeptide(), it.psm.getCharge(), it.delta, scan);
					states.merge(c.localisation.state, 1, Integer::sum);
					if (c.localisation.state == ShiftLocaliser.State.LOCALISED_EXPLAINED) {
						c.group = PeakSiteShare.isSubstitutionSite(c.getPeptide().charAt(c.localisation.site),
								it.delta, SITE_CONFIG) ? JointCandidate.Group.OPEN_SUBSTITUTION_SITE
										: JointCandidate.Group.OPEN_LOCALISED;
					}
				}
			}
		}
		if (noSpectra > 0) {
			warning.accept("Joint identification: no readable mzML for " + noSpectra
					+ " unexplained- or accepted-shift winners; they keep their group's export");
		}

		info.accept("Joint identification: shift localisation (" + cfg.shiftLocaliser + ") of unexplained-shift "
				+ "winners " + states + "; check of accepted-shift and isotope-offset winners " + acceptedShift);
		if (peaks != null) {
			info.accept(String.format(Locale.ROOT, "Joint identification: %,d exact-mass winners with an isotope offset: "
					+ "%,d in an isotope peak of the zero family, %,d in another fitted peak and checked like an "
					+ "open shift", closedIsotope, closedInIsotopePeak, closedChecked));
		}
		return new Localisation(byFile, new ArrayList<>(nullSample), peaks, cfg);
	}

	/** What {@link #localise} found: scored winners and the null pool, not yet tested or labelled. */
	static final class Localisation {
		private final Map<String, List<ShiftItem>> byFile;
		private final List<NullPsm> nullSample;
		private final PeakContext peaks;
		private final Config cfg;

		Localisation(Map<String, List<ShiftItem>> byFile, List<NullPsm> nullSample, PeakContext peaks, Config cfg) {
			this.byFile = byFile;
			this.nullSample = nullSample;
			this.peaks = peaks;
			this.cfg = cfg;
		}

		/**
		 * {@link PeakFragmentTest} on every peak family with localised winners. {@code r1} counts the
		 * winners {@code accepted} marks (all winners when null); {@code r0} the first
		 * {@link Config#peakTestNullSample} null-pool PSMs it marks, so both rates come from one population.
		 * Below the FDR cut the zero-peak winners include spectra of modified peptides (12mix: the jar's
		 * null at +0.984 was 10% from every winner against 1.07% from 1,500 accepted ones, scored by the
		 * same {@link ShiftLocaliser#supportAt}). Every winner of a peak that fails,
		 * accepted or not, becomes {@link ShiftLocaliser.State#UNSUPPORTED}.
		 */
		List<PeakFragmentTest.PeakResult> testPeaks(boolean[] accepted, Consumer<String> info) {
			if (peaks == null) {
				return List.of();
			}
			return JointIdentification.testPeaks(byFile, nullSample, peaks, cfg, accepted, info);
		}

		/**
		 * The per-PSM rule for accepted-shift winners, after the peak test has had its say: one rule for
		 * what stays in the main result, whatever the peak as a whole looks like.
		 * <ul>
		 * <li>Explained on a residue a substitution fits: {@code OPEN_SUBSTITUTION_SITE}, as off a peak.
		 * {@link PeakSiteShare} withdraws a peak at a 0.6 share of such sites; below it they stayed
		 * identified (round 5: 508 faecal and 4,336 12mix PSMs, +40.031 above all).</li>
		 * <li>Fragments place the shift where no explanation fits: {@code OPEN_UNEXPLAINED} (1,627 and
		 * 4,685 PSMs).</li>
		 * <li>Fragments do not place the shift (this PSM's, or its whole peak's by the peak test): the
		 * modification UNIMOD allows for the peak's mass on this peptide ({@link ShiftLocaliser#assign}),
		 * exported as that modified peptide; without one, {@code OPEN_UNEXPLAINED}. These PSMs were
		 * exported as the unmodified peptide and quantified at a mass their precursor did not have
		 * (6,219 and 51,251 PSMs; 425 and 8,628 have a UNIMOD site).</li>
		 * </ul>
		 * Label-blind, and before the q-values that count, so each group's walk sees its final members.
		 * Exact-mass winners and winners without a readable spectrum are left as they are.
		 */
		void decideAcceptedShifts(Consumer<String> info) {
			int substitutionSite = 0;
			int contradicted = 0;
			int assigned = 0;
			int unaccounted = 0;
			for (List<ShiftItem> items : byFile.values()) {
				for (ShiftItem it : items) {
					JointCandidate c = it.candidate;
					ShiftLocaliser.Result loc = c.localisation;
					if (c.group != JointCandidate.Group.OPEN_ACCEPTED_SHIFT || loc == null || loc.site < 0) {
						continue;
					}
					if (loc.state == ShiftLocaliser.State.LOCALISED_EXPLAINED) {
						if (PeakSiteShare.isSubstitutionSite(c.getPeptide().charAt(loc.site), it.delta, SITE_CONFIG)) {
							c.group = JointCandidate.Group.OPEN_SUBSTITUTION_SITE;
							substitutionSite++;
						}
					} else if (loc.state == ShiftLocaliser.State.LOCALISED_UNEXPLAINED) {
						c.group = JointCandidate.Group.OPEN_UNEXPLAINED;
						contradicted++;
					} else if (loc.state == ShiftLocaliser.State.UNSUPPORTED) {
						Double centroid = peaks == null ? null : peaks.centroids.get(it.family);
						double mass = centroid != null && Double.isFinite(centroid) ? centroid : it.delta;
						ShiftLocaliser.Result byUnimod = cfg.shiftLocaliser.assign(c.getPeptide(), mass, loc.support);
						if (byUnimod != null) {
							c.localisation = byUnimod;
							assigned++;
						} else {
							c.group = JointCandidate.Group.OPEN_UNEXPLAINED;
							unaccounted++;
						}
					}
				}
			}
			info.accept("Joint identification: accepted-shift winners by PSM: " + substitutionSite
					+ " accepted-shift winners on a substitution site moved to " + JointCandidate.Group.OPEN_SUBSTITUTION_SITE
					+ "; " + contradicted + " with the shift on a residue nothing explains and " + unaccounted
					+ " unplaced with no UNIMOD site in the peptide moved to " + JointCandidate.Group.OPEN_UNEXPLAINED
					+ "; " + assigned + " unplaced given the UNIMOD modification of their peak's mass");
		}

		/** Gives each localised winner's PSMs its group and localisation. */
		void apply() {
			for (List<ShiftItem> items : byFile.values()) {
				for (ShiftItem it : items) {
					JointCandidate c = it.candidate;
					if (c.localisation == null) {
						continue;
					}
					for (OpenSearchPsm psm : new OpenSearchPsm[] { c.closed, c.fragger, c.sage }) {
						if (psm != null) {
							psm.setIdentificationGroup(c.group.name());
							psm.setShiftLocalisation(c.localisation);
						}
					}
				}
			}
		}
	}

	/** The winner as an unmodified zero-peak PSM for the null pool, or null when it is not one. */
	private static NullPsm nullCandidate(int index, JointCandidate c) {
		if (c.group == JointCandidate.Group.EXACT && c.closed != null) {
			OpenSearchPsm z = c.closed;
			if (z.getClosedModificationCount() == 0 && Math.abs(z.getMassDiff()) < NULL_ZERO_WINDOW_DA) {
				return new NullPsm(index, c.getFileName(), c.getScan(), c.getPeptide(), z.getCharge());
			}
			return null;
		}
		if (c.group == JointCandidate.Group.OPEN_UNMODIFIED) {
			OpenSearchPsm z = c.openPrimary();
			if (z != null && Math.abs(z.getMassDiff()) < NULL_ZERO_WINDOW_DA) {
				return new NullPsm(index, c.getFileName(), c.getScan(), c.getPeptide(), z.getCharge());
			}
		}
		return null;
	}

	/** Mass difference within which a zero-peak winner counts as unmodified for the null sample (Da). */
	static final double NULL_ZERO_WINDOW_DA = 0.02;

	/**
	 * Null pool size as a multiple of the null sample: spectra are read for this many zero-peak winners
	 * before q-values exist, and the first {@link Config#peakTestNullSample} accepted ones are the null.
	 */
	static final int NULL_POOL_FACTOR = 3;

	/** See {@link Localisation#testPeaks}. */
	private static List<PeakFragmentTest.PeakResult> testPeaks(Map<String, List<ShiftItem>> byFile,
			List<NullPsm> nullPool, PeakContext peaks, Config cfg, boolean[] accepted, Consumer<String> info) {
		// The null: the first peakTestNullSample pool PSMs (stable-hash order) with spectra, accepted as r1's are.
		List<NullPsm> nullSample = new ArrayList<>();
		for (NullPsm z : nullPool) {
			if (nullSample.size() >= cfg.peakTestNullSample) {
				break;
			}
			if (z.peaks != null && (accepted == null || accepted[z.index])) {
				nullSample.add(z);
			}
		}
		Map<Integer, List<ShiftItem>> byFamily = new java.util.TreeMap<>();
		for (List<ShiftItem> items : byFile.values()) {
			for (ShiftItem it : items) {
				if (it.family >= 0 && it.candidate.localisation != null) {
					byFamily.computeIfAbsent(it.family, k -> new ArrayList<>()).add(it);
				}
			}
		}
		double min = cfg.shiftLocaliser.minSupport();
		List<PeakFragmentTest.PeakResult> out = new ArrayList<>();
		int failed = 0;
		int withdrawn = 0;
		for (Map.Entry<Integer, List<ShiftItem>> e : byFamily.entrySet()) {
			int family = e.getKey();
			double centroid = peaks.centroids.getOrDefault(family, Double.NaN);
			int scored = 0;
			int supported = 0;
			for (ShiftItem it : e.getValue()) {
				ShiftLocaliser.Result r = it.candidate.localisation;
				if ((accepted == null || accepted[it.index]) && r.site >= 0 && Double.isFinite(r.support)) {
					scored++;
					supported += r.support >= min ? 1 : 0;
				}
			}
			int nullScored = 0;
			int nullSupported = 0;
			if (Double.isFinite(centroid)) {
				for (NullPsm z : nullSample) {
					double s = cfg.shiftLocaliser.supportAt(z.peptide, z.charge, centroid, z.peaks);
					if (Double.isFinite(s)) {
						nullScored++;
						nullSupported += s >= min ? 1 : 0;
					}
				}
			}
			PeakFragmentTest.PeakResult r = PeakFragmentTest.decide(family, centroid, scored, supported, nullScored,
					nullSupported, cfg.peakTestMaxChanceShare);
			out.add(r);
			if (!r.passed) {
				failed++;
				for (ShiftItem it : e.getValue()) {
					ShiftLocaliser.Result before = it.candidate.localisation;
					it.candidate.localisation = ShiftLocaliser.asUnsupported(before);
					withdrawn += before.state != ShiftLocaliser.State.UNSUPPORTED ? 1 : 0;
				}
			}
		}
		info.accept(String.format(Locale.ROOT, "Joint identification: per-peak fragment test (bounded chance share <= "
				+ "%.2f, %,d %s unmodified PSMs as the null, %s): %d of %d peaks not on the fragments; %,d localisations "
				+ "there withdrawn", cfg.peakTestMaxChanceShare, nullSample.size(),
				accepted == null ? "winning" : "accepted", accepted == null ? "every winner" : "accepted winners", failed, out.size(), withdrawn));
		return out;
	}

	private static final PeakSiteShare.Config SITE_CONFIG = new PeakSiteShare.Config();

	/**
	 * Whether a group's accepted PSMs are identifications. {@code OPEN_SUBSTITUTION_SITE} never is,
	 * even with fragment support: the fragments confirm the shift sits on a residue a substitution
	 * could change, not that the rest of the peptide is right. On 12mix (2026-11-08) these PSMs were
	 * 12-19% outside the true genomes against 2-4% for exact-mass and {@code OPEN_LOCALISED} PSMs, and
	 * their one-residue variant was almost never a true peptide (5 of 378): wrong peptides, not
	 * variants. They are reported in {@code unexplained_psms.tsv} and read by the variant outputs.
	 */
	static boolean identifies(JointCandidate.Group group, Config cfg) {
		return group != JointCandidate.Group.OPEN_SUBSTITUTION_SITE
				&& (group != JointCandidate.Group.OPEN_UNEXPLAINED || cfg.acceptUnexplained);
	}

	/** Highest score per spectrum; ties to the earlier row, so reruns agree. */
	/** RESET role of a decoy peptide (I = L): training when its stable hash is even, as the unified method. */
	static boolean[] resetRoles(List<JointCandidate> candidates) {
		boolean[] training = new boolean[candidates.size()];
		for (int i = 0; i < training.length; i++) {
			JointCandidate c = candidates.get(i);
			training[i] = !c.isTarget() && (SemiSupervisedPsmClassifier.stableHash(
					"decoy:" + c.getPeptide().toUpperCase(Locale.ROOT).replace('I', 'L')) & 1L) == 0L;
		}
		return training;
	}

	/** As {@link #winners(int[], double[])}; a decoy wins a tied query, whatever the row order. */
	static boolean[] winners(int[] spectrumId, double[] scores, List<JointCandidate> candidates) {
		int spectra = 0;
		for (int id : spectrumId) {
			spectra = Math.max(spectra, id + 1);
		}
		int[] best = new int[spectra];
		Arrays.fill(best, -1);
		for (int i = 0; i < scores.length; i++) {
			int s = spectrumId[i];
			int b = best[s];
			if (b < 0 || rank(scores[i]) > rank(scores[b]) || (rank(scores[i]) == rank(scores[b])
					&& !candidates.get(i).isTarget() && candidates.get(b).isTarget())) {
				best[s] = i;
			}
		}
		boolean[] out = new boolean[scores.length];
		for (int b : best) {
			if (b >= 0) {
				out[b] = true;
			}
		}
		return out;
	}

	static boolean[] winners(int[] spectrumId, double[] scores) {
		int spectra = 0;
		for (int id : spectrumId) {
			spectra = Math.max(spectra, id + 1);
		}
		int[] best = new int[spectra];
		Arrays.fill(best, -1);
		for (int i = 0; i < scores.length; i++) {
			int s = spectrumId[i];
			if (best[s] < 0 || rank(scores[i]) > rank(scores[best[s]])) {
				best[s] = i;
			}
		}
		boolean[] out = new boolean[scores.length];
		for (int b : best) {
			if (b >= 0) {
				out[b] = true;
			}
		}
		return out;
	}

	private static double rank(double score) {
		return Double.isNaN(score) ? Double.NEGATIVE_INFINITY : score;
	}

	/**
	 * Target-decoy q-value of each winner, walked within its own group: {@code (decoys + 1) /
	 * targets}, tied scores as one block ({@link TargetDecoyFdr#qValuesConservative}), because the
	 * open groups can be small.
	 */
	static double[] groupQValues(List<JointCandidate> candidates, double[] scores, boolean[] winner) {
		return groupQValues(candidates, scores, winner, false, 1.0);
	}

	/**
	 * @param byPrecursorClass FDR groups split by precursor class (same, reassigned): identification v2
	 * @param decoyWeight      each decoy counts this many times (2 under RESET)
	 */
	static double[] groupQValues(List<JointCandidate> candidates, double[] scores, boolean[] winner,
			boolean byPrecursorClass, double decoyWeight) {
		double[] q = new double[scores.length];
		Arrays.fill(q, Double.NaN);
		Map<String, List<Integer>> byGroup = new java.util.TreeMap<>();
		for (int i = 0; i < scores.length; i++) {
			if (winner[i]) {
				JointCandidate c = candidates.get(i);
				String key = c.group.name() + (byPrecursorClass && c.isReassigned() ? ":reassigned" : "");
				byGroup.computeIfAbsent(key, g -> new ArrayList<>()).add(i);
			}
		}
		for (List<Integer> rows : byGroup.values()) {
			double[] s = new double[rows.size()];
			boolean[] t = new boolean[rows.size()];
			for (int k = 0; k < rows.size(); k++) {
				s[k] = scores[rows.get(k)];
				t[k] = candidates.get(rows.get(k)).isTarget();
			}
			double[] gq = TargetDecoyFdr.qValuesConservative(s, t, decoyWeight);
			for (int k = 0; k < rows.size(); k++) {
				q[rows.get(k)] = gq[k];
			}
		}
		return q;
	}

	/** Peptide unit: target or decoy, base sequence with I and L merged. */
	static String unitKey(JointCandidate c) {
		return (c.isTarget() ? "T:" : "D:") + c.getPeptide().toUpperCase(Locale.ROOT).replace('I', 'L');
	}

	static final class PeptideWalk {
		final Map<String, JointCandidate.Group> unitGroup = new HashMap<>();
		final Map<String, Double> unitQ = new HashMap<>();
		/** Per group {targets, decoys} among walked units, and at q <= fdr. */
		final EnumMap<JointCandidate.Group, int[]> walked = new EnumMap<>(JointCandidate.Group.class);
		final EnumMap<JointCandidate.Group, int[]> passed = new EnumMap<>(JointCandidate.Group.class);
	}

	/**
	 * Each peptide unit in the highest-priority group where it has an accepted PSM, scored by its
	 * best accepted PSM in that group, then a target-decoy walk per group.
	 */
	static PeptideWalk peptideWalk(List<JointCandidate> candidates, double[] scores, boolean[] accepted, double fdr) {
		return peptideWalk(candidates, scores, accepted, fdr, 1.0);
	}

	static PeptideWalk peptideWalk(List<JointCandidate> candidates, double[] scores, boolean[] accepted, double fdr,
			double decoyWeight) {
		PeptideWalk walk = new PeptideWalk();
		for (int i = 0; i < scores.length; i++) {
			if (!accepted[i]) {
				continue;
			}
			String key = unitKey(candidates.get(i));
			JointCandidate.Group g = candidates.get(i).group;
			JointCandidate.Group held = walk.unitGroup.get(key);
			if (held == null || g.ordinal() < held.ordinal()) {
				walk.unitGroup.put(key, g);
			}
		}
		Map<String, Double> best = new HashMap<>();
		for (int i = 0; i < scores.length; i++) {
			if (!accepted[i]) {
				continue;
			}
			String key = unitKey(candidates.get(i));
			if (walk.unitGroup.get(key) == candidates.get(i).group) {
				best.merge(key, scores[i], Math::max);
			}
		}
		EnumMap<JointCandidate.Group, List<String>> byGroup = new EnumMap<>(JointCandidate.Group.class);
		for (Map.Entry<String, JointCandidate.Group> e : walk.unitGroup.entrySet()) {
			byGroup.computeIfAbsent(e.getValue(), g -> new ArrayList<>()).add(e.getKey());
		}
		for (Map.Entry<JointCandidate.Group, List<String>> e : byGroup.entrySet()) {
			List<String> keys = e.getValue();
			Collections.sort(keys);
			double[] s = new double[keys.size()];
			boolean[] t = new boolean[keys.size()];
			int[] w = new int[2];
			for (int k = 0; k < keys.size(); k++) {
				s[k] = best.get(keys.get(k));
				t[k] = keys.get(k).startsWith("T:");
				w[t[k] ? 0 : 1]++;
			}
			double[] q = TargetDecoyFdr.qValuesConservative(s, t, decoyWeight);
			int[] p = new int[2];
			for (int k = 0; k < keys.size(); k++) {
				walk.unitQ.put(keys.get(k), q[k]);
				if (q[k] <= fdr) {
					p[t[k] ? 0 : 1]++;
				}
			}
			walk.walked.put(e.getKey(), w);
			walk.passed.put(e.getKey(), p);
		}
		return walk;
	}

	/**
	 * The workflow's validation result. Identifications are the kept winners' reported rows (the
	 * closed row for an exact candidate); the open groups' kept rows are the peak-validated set; the
	 * unexplained group's winners at {@code psmQ <= fdr} are the reported, unidentified shifts; and
	 * the variant outputs read the open rows of both.
	 */
	static ModificationValidator.Result toValidation(List<JointCandidate> candidates, boolean[] winner,
			double[] psmQ, boolean[] accepted, boolean[] kept, ModificationValidator.Gate gate, Config cfg) {
		List<OpenSearchPsm> validated = new ArrayList<>();
		List<OpenSearchPsm> peakValidated = new ArrayList<>();
		List<OpenSearchPsm> unexplained = new ArrayList<>();
		Map<Integer, Integer> retained = new HashMap<>();
		// Every row of every candidate gets this method's q-value: a winner's group q on all of its
		// searches' rows, 1.0 (not accepted) on the rest. A row left alone would keep whatever an
		// earlier step - the comparison arms walk the same objects - last wrote there.
		for (int i = 0; i < candidates.size(); i++) {
			JointCandidate c = candidates.get(i);
			double q = winner[i] ? psmQ[i] : 1.0;
			for (OpenSearchPsm psm : new OpenSearchPsm[] { c.closed, c.fragger, c.sage }) {
				if (psm != null) {
					psm.setQValue(q);
				}
			}
		}
		for (int i = 0; i < candidates.size(); i++) {
			if (!winner[i]) {
				continue;
			}
			JointCandidate c = candidates.get(i);
			OpenSearchPsm row = c.representative();
			if (kept[i]) {
				validated.add(row);
				if (c.group != JointCandidate.Group.EXACT) {
					peakValidated.add(row);
					if (row.getPeakId() >= 0) {
						retained.merge(row.getPeakId(), 1, Integer::sum);
					}
				}
			} else if ((c.group == JointCandidate.Group.OPEN_UNEXPLAINED
					|| c.group == JointCandidate.Group.OPEN_SUBSTITUTION_SITE) && !accepted[i] && psmQ[i] <= cfg.fdr) {
				unexplained.add(row);
			}
		}
		List<ModificationValidator.PeakVerdict> verdicts = new ArrayList<>(gate.verdicts.size());
		for (ModificationValidator.PeakVerdict v : gate.verdicts) {
			verdicts.add(v.withRetained(retained.getOrDefault(v.peakId, 0)));
		}
		List<OpenSearchPsm> diagnostic = new ArrayList<>(peakValidated.size() + unexplained.size());
		diagnostic.addAll(peakValidated);
		diagnostic.addAll(unexplained);
		return new ModificationValidator.Result(validated, peakValidated, unexplained, verdicts, diagnostic);
	}

	// ------------------------------------------------------------ reporting

	private static List<String> summarise(List<JointCandidate> candidates, boolean[] winner, double[] psmQ,
			boolean[] accepted, boolean[] kept, PeptideWalk walk, Config cfg, int both, JointFeatures features,
			JointPsmClassifier.Result scored) {
		List<String> lines = new ArrayList<>();
		lines.add("candidates\t" + candidates.size());
		lines.add("candidates_proposed_by_both_searches\t" + both);
		lines.add("fdr_per_group\t" + cfg.fdr);
		lines.add("unexplained_shift_identifications\t" + cfg.acceptUnexplained);
		lines.add("features\t" + String.join(",", features.featureNames()));
		lines.add("features_withheld\t" + String.join(",", features.withheldColumns()));
		lines.add("classifier_iterations\t" + scored.iterationsUsed);
		EnumMap<JointCandidate.Group, int[]> winners = new EnumMap<>(JointCandidate.Group.class);
		EnumMap<JointCandidate.Group, int[]> acc = new EnumMap<>(JointCandidate.Group.class);
		EnumMap<JointCandidate.Group, int[]> keep = new EnumMap<>(JointCandidate.Group.class);
		int overruledClosed = 0;
		int closedRankOneSpectra = 0;
		for (int i = 0; i < candidates.size(); i++) {
			if (!winner[i]) {
				continue;
			}
			JointCandidate c = candidates.get(i);
			int side = c.isTarget() ? 0 : 1;
			winners.computeIfAbsent(c.group, g -> new int[2])[side]++;
			if (accepted[i]) {
				acc.computeIfAbsent(c.group, g -> new int[2])[side]++;
			}
			if (kept[i]) {
				keep.computeIfAbsent(c.group, g -> new int[2])[side]++;
			}
		}
		// How often the winner is not the closed search's own rank-1 candidate, on spectra that have one.
		Map<Integer, Integer> closedRankOne = new HashMap<>();
		int[] sid = features.spectrumIds();
		for (int i = 0; i < candidates.size(); i++) {
			JointCandidate c = candidates.get(i);
			if (c.closed != null && c.closed.getHitRank() == 1) {
				closedRankOne.put(sid[i], i);
			}
		}
		for (int i = 0; i < candidates.size(); i++) {
			if (winner[i] && closedRankOne.containsKey(sid[i])) {
				closedRankOneSpectra++;
				if (closedRankOne.get(sid[i]) != i) {
					overruledClosed++;
				}
			}
		}
		lines.add("spectra_with_closed_rank1\t" + closedRankOneSpectra);
		lines.add("spectra_where_winner_is_not_closed_rank1\t" + overruledClosed);
		int[] unionPsm = new int[2];
		for (JointCandidate.Group g : JointCandidate.Group.values()) {
			int[] w = winners.getOrDefault(g, new int[2]);
			int[] a = acc.getOrDefault(g, new int[2]);
			int[] k = keep.getOrDefault(g, new int[2]);
			int[] pw = walk.walked.getOrDefault(g, new int[2]);
			int[] pp = walk.passed.getOrDefault(g, new int[2]);
			lines.add(g.name().toLowerCase(Locale.ROOT) + "_spectrum_winners\t[" + w[0] + ", " + w[1] + "]");
			lines.add(g.name().toLowerCase(Locale.ROOT) + "_psms_at_fdr\t[" + a[0] + ", " + a[1] + "]");
			lines.add(g.name().toLowerCase(Locale.ROOT) + "_peptides_walked\t[" + pw[0] + ", " + pw[1] + "]");
			lines.add(g.name().toLowerCase(Locale.ROOT) + "_peptides_at_fdr\t[" + pp[0] + ", " + pp[1] + "]");
			lines.add(g.name().toLowerCase(Locale.ROOT) + "_psms_kept\t[" + k[0] + ", " + k[1] + "]");
			unionPsm[0] += k[0];
			unionPsm[1] += k[1];
		}
		int[] unionPep = new int[2];
		for (int[] p : walk.passed.values()) {
			unionPep[0] += p[0];
			unionPep[1] += p[1];
		}
		lines.add("identifications_psms\t[" + unionPsm[0] + ", " + unionPsm[1] + "]");
		lines.add(String.format(Locale.ROOT, "identifications_psm_decoy_ratio\t%.5f",
				unionPsm[1] / (double) Math.max(1, unionPsm[0])));
		lines.add("identifications_peptides\t[" + unionPep[0] + ", " + unionPep[1] + "]");
		lines.add(String.format(Locale.ROOT, "identifications_peptide_decoy_ratio\t%.5f",
				unionPep[1] / (double) Math.max(1, unionPep[0])));
		int[] openOnly = new int[2];
		for (JointCandidate.Group g : new JointCandidate.Group[] { JointCandidate.Group.OPEN_UNMODIFIED,
				JointCandidate.Group.OPEN_ACCEPTED_SHIFT, JointCandidate.Group.OPEN_LOCALISED,
				JointCandidate.Group.OPEN_UNEXPLAINED }) {
			int[] p = walk.passed.getOrDefault(g, new int[2]);
			openOnly[0] += p[0];
			openOnly[1] += p[1];
		}
		lines.add("open_only_peptides\t[" + openOnly[0] + ", " + openOnly[1] + "]");
		lines.addAll(scored.modelDescriptions);
		return lines;
	}

	/**
	 * Per feature, how well it separates target from decoy spectrum winners: ROC AUC per group and
	 * overall, the share of winners it is present for. A diagnostic for choosing features - decoys
	 * are known negatives, so this uses no information the FDR walks do not.
	 */
	public static void writeFeatureReport(Result result, File out) throws IOException {
		JointFeatures f = result.features;
		String[] names = f.featureNames();
		JointCandidate.Group[] groups = JointCandidate.Group.values();
		List<Integer> rows = new ArrayList<>();
		for (int i = 0; i < result.candidates.size(); i++) {
			if (result.winner[i]) {
				rows.add(i);
			}
		}
		int m = rows.size();
		boolean[] targets = new boolean[m];
		int[] group = new int[m];
		for (int r = 0; r < m; r++) {
			JointCandidate c = result.candidates.get(rows.get(r));
			targets[r] = c.isTarget();
			group[r] = c.group.ordinal();
		}
		out.getParentFile().mkdirs();
		try (PrintWriter w = new PrintWriter(out, "UTF-8")) {
			StringBuilder header = new StringBuilder("feature");
			for (JointCandidate.Group g : groups) {
				header.append('\t').append(g.name().toLowerCase(Locale.ROOT)).append("_auc");
			}
			header.append("\tall_auc\tpresent_fraction");
			w.println("# ROC AUC of each active feature between target and decoy spectrum winners (0.5 = no "
					+ "separation; below 0.5 = higher values mark decoys). Coefficients are in the summary.");
			w.println(header);
			// One column at a time: a winners x features matrix would hold hundreds of megabytes.
			double[] values = new double[m];
			for (int k = 0; k < names.length; k++) {
				int present = 0;
				for (int r = 0; r < m; r++) {
					values[r] = f.values(rows.get(r))[k];
					if (!Double.isNaN(values[r])) {
						present++;
					}
				}
				StringBuilder row = new StringBuilder(names[k]);
				for (JointCandidate.Group g : groups) {
					int size = 0;
					for (int r = 0; r < m; r++) {
						if (group[r] == g.ordinal()) {
							size++;
						}
					}
					double[] v = new double[size];
					boolean[] t = new boolean[size];
					int j = 0;
					for (int r = 0; r < m; r++) {
						if (group[r] == g.ordinal()) {
							v[j] = values[r];
							t[j++] = targets[r];
						}
					}
					row.append('\t').append(format(auc(v, t)));
				}
				row.append('\t').append(format(auc(values, targets)));
				row.append('\t').append(format(present / (double) Math.max(1, m)));
				w.println(row);
			}
		}
	}

	/** Mann-Whitney AUC of targets over decoys on non-missing values; NaN without both classes. */
	static double auc(double[] values, boolean[] targets) {
		List<Integer> idx = new ArrayList<>();
		for (int i = 0; i < values.length; i++) {
			if (!Double.isNaN(values[i])) {
				idx.add(i);
			}
		}
		idx.sort((a, b) -> Double.compare(values[a], values[b]));
		double rankSumTargets = 0;
		long nt = 0;
		long nd = 0;
		int i = 0;
		while (i < idx.size()) {
			int j = i;
			while (j + 1 < idx.size() && values[idx.get(j + 1)] == values[idx.get(i)]) {
				j++;
			}
			double midRank = (i + j) / 2.0 + 1.0;
			for (int k = i; k <= j; k++) {
				if (targets[idx.get(k)]) {
					rankSumTargets += midRank;
					nt++;
				} else {
					nd++;
				}
			}
			i = j + 1;
		}
		if (nt == 0 || nd == 0) {
			return Double.NaN;
		}
		return (rankSumTargets - nt * (nt + 1) / 2.0) / ((double) nt * nd);
	}

	private static String format(double v) {
		return Double.isNaN(v) ? "" : String.format(Locale.ROOT, "%.4f", v);
	}

	/** Rows of every search, each once, for the audit table. */
	public static List<OpenSearchPsm> allRows(List<JointCandidate> candidates) {
		Set<OpenSearchPsm> seen = Collections.newSetFromMap(new IdentityHashMap<>());
		List<OpenSearchPsm> rows = new ArrayList<>();
		for (JointCandidate c : candidates) {
			for (OpenSearchPsm psm : new OpenSearchPsm[] { c.closed, c.fragger, c.sage }) {
				if (psm != null && seen.add(psm)) {
					rows.add(psm);
				}
			}
		}
		return rows;
	}

	private static void report(Consumer<String> progress, String message) {
		LOGGER.info(message);
		if (progress != null) {
			progress.accept(message);
		}
	}
}
