package uk.ac.quadram.metapilot.core.opensearch;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import weka.core.Attribute;
import weka.core.DenseInstance;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.Utils;

/**
 * Feature view over {@link JointCandidate}: one row per (spectrum, peptide), with one block per
 * search, presence indicators, the spectrum's competition and, when DIA-NN ran, predicted-spectrum
 * and retention-time agreement.
 *
 * <h2>Blocks</h2>
 * <ul>
 * <li><b>Shared</b> - properties of the candidate whichever search found it.</li>
 * <li><b>Presence</b> - which searches proposed it. Agreement between the closed and the open
 * search is evidence: the same peptide survives both a narrow and a wide precursor window.</li>
 * <li><b>Closed</b>, <b>open MSFragger</b>, <b>open Sage</b> - each search's own scores, missing
 * when that search did not propose the candidate. Never merged into one column: an open
 * expect value is computed over a far larger candidate pool than a closed one.</li>
 * <li><b>Open peak context</b> - the delta-mass hypothesis of an open candidate.</li>
 * <li><b>Spectrum competition</b> - what else explains the spectrum: whether the closed search
 * found any exact-mass peptide for it, how strong the best <em>other</em> closed candidate is,
 * and the hyperscore margin over the best other candidate. This is what lets an exact-mass
 * alternative demote an open near-miss, and a weak closed answer lose to a strong open one.</li>
 * <li><b>Predicted spectrum</b> - spectral-entropy similarity to DIA-NN's predicted fragments,
 * the fraction of predicted intensity observed, the similarity margin over the spectrum's other
 * candidates, and the calibrated retention-time error.</li>
 * </ul>
 * Every value is a property of the candidate or of its spectrum, computed identically for
 * targets and decoys; the target/decoy label is the class and never an input.
 *
 * <h2>Missing values and column selection</h2>
 * An absent block is {@link Utils#missingValue()}, never imputed here; weka's {@code Logistic}
 * replaces missing values with training means, which with the presence indicators beside them is
 * a per-block offset. A column with no value on any row, or one value on every row, carries no
 * information and is withheld, so a run without Sage or without DIA-NN trains on what it has.
 *
 * <p><b>Not thread safe.</b> One scratch buffer is reused across {@link #instance} calls.
 */
public final class JointFeatures {

	/**
	 * No same-peptide or same-protein count. Those aggregate over targets and decoys separately, and a
	 * false target lands on a peptide or protein that true identifications also hit, a decoy never
	 * does: the count marks target identity rather than correctness. On a synthetic run with known
	 * truth (refined databases are mostly the genomes present, so random target matches hit present
	 * proteins) the two counts took the exact group's measured false discovery proportion to 7.2%
	 * at a decoy estimate of 0.3%; without them it was 0.38%.
	 */
	static final String[] SHARED = { "charge", "precursor_mass", "peptide_length", "num_tol_term",
			"num_missed_cleavages" };
	static final String[] PRESENCE = { "found_by_closed", "found_by_msfragger_open", "found_by_sage_open" };
	static final String[] CLOSED = { "closed_hit_rank", "closed_num_matched_ions", "closed_matched_ion_fraction",
			"closed_hyperscore", "closed_hyperscore_minus_nextscore", "closed_neg_log10_expect",
			"closed_abs_mass_error_ppm", "closed_isotope_offset", "closed_variable_modifications" };
	static final String[] OPEN = { "open_hit_rank", "open_num_matched_ions", "open_matched_ion_fraction",
			"open_hyperscore", "open_hyperscore_minus_nextscore", "open_neg_log10_expect",
			"open_delta_mass_score_gain", "open_localization_delta_score" };
	static final String[] PEAK = { "open_abs_normalized_mass_diff", "open_zero_family", "open_accepted_shift",
			"open_in_validated_peak", "open_distance_in_sigmas", "open_log_local_density",
			"open_log_same_peak_count", "open_peak_isotope_offset" };
	static final String[] SAGE = { "sage_hit_rank", "sage_matched_ion_fraction", "sage_hyperscore",
			"sage_hyperscore_minus_nextscore", "sage_neg_log10_expect" };
	static final String[] SPECTRUM = { "spectrum_has_closed", "closed_rival_neg_log10_expect",
			"hyperscore_margin" };
	static final String[] MS2 = { "ms2_entropy", "ms2_matched_fraction", "ms2_entropy_margin", "rt_error_norm" };

	/** Every column this class can compute, in fill order. */
	public static final String[] ALL_FEATURES = concat(SHARED, PRESENCE, CLOSED, OPEN, PEAK, SAGE, SPECTRUM, MS2);

	public static final String CLASS_ATTRIBUTE = "label";
	public static final String TARGET_LABEL = "target";
	public static final String DECOY_LABEL = "decoy";

	/** Predicted-spectrum values per candidate, NaN where not computed. */
	public static final class SpectrumScores {
		public final float[] entropy;
		public final float[] matchedFraction;
		public final float[] rtErrorNorm;

		public SpectrumScores(int n) {
			entropy = nanArray(n);
			matchedFraction = nanArray(n);
			rtErrorNorm = nanArray(n);
		}

		private static float[] nanArray(int n) {
			float[] a = new float[n];
			Arrays.fill(a, Float.NaN);
			return a;
		}
	}

	private final List<JointCandidate> candidates;
	private final float[] spectrumHasClosed;
	private final float[] closedRival;
	private final float[] hyperscoreMargin;
	private final SpectrumScores ms2;
	private final float[] ms2Margin;
	private final int[] spectrumId;
	private final java.util.Set<Integer> zeroPeakIds;
	private int[] activeColumns;
	private String[] activeNames;
	private final double[] scratch = new double[ALL_FEATURES.length];

	private JointFeatures(List<JointCandidate> candidates, SpectrumScores ms2, java.util.Set<Integer> zeroPeakIds) {
		this.candidates = candidates;
		this.zeroPeakIds = zeroPeakIds == null ? java.util.Collections.emptySet() : zeroPeakIds;
		int n = candidates.size();
		this.ms2 = ms2;
		this.spectrumId = assignSpectrumIds(candidates);
		this.spectrumHasClosed = new float[n];
		this.closedRival = new float[n];
		this.hyperscoreMargin = new float[n];
		this.ms2Margin = new float[n];
		computeContext();
	}

	/**
	 * Builds the view and chooses its columns: every computable column that has at least two
	 * distinct values, minus {@code withheld}.
	 *
	 * @param ms2         predicted-spectrum values, or null when DIA-NN did not run
	 * @param zeroPeakIds peak ids of the accepted zero family (unmodified, isotope errors folded in)
	 * @param withheld    feature names never to use (may be null)
	 */
	public static JointFeatures of(List<JointCandidate> candidates, SpectrumScores ms2,
			java.util.Set<Integer> zeroPeakIds, java.util.Set<String> withheld) {
		JointFeatures f = new JointFeatures(candidates, ms2, zeroPeakIds);
		f.chooseColumns(withheld);
		return f;
	}

	public int size() {
		return candidates.size();
	}

	public boolean isTarget(int row) {
		return candidates.get(row).isTarget();
	}

	public int spectrumId(int row) {
		return spectrumId[row];
	}

	public int[] spectrumIds() {
		return spectrumId;
	}

	public String[] featureNames() {
		return activeNames.clone();
	}

	/** Columns the view could compute but withheld: no information on this data, or withheld by name. */
	public List<String> withheldColumns() {
		List<String> out = new ArrayList<>();
		java.util.Set<String> active = new java.util.HashSet<>(Arrays.asList(activeNames));
		for (String name : ALL_FEATURES) {
			if (!active.contains(name)) {
				out.add(name);
			}
		}
		return out;
	}

	public Instances emptyDataset(String relationName, int capacity) {
		ArrayList<Attribute> attributes = new ArrayList<>(activeNames.length + 1);
		for (String name : activeNames) {
			attributes.add(new Attribute(name));
		}
		attributes.add(new Attribute(CLASS_ATTRIBUTE, Arrays.asList(TARGET_LABEL, DECOY_LABEL)));
		Instances instances = new Instances(relationName, attributes, Math.max(0, capacity));
		instances.setClassIndex(attributes.size() - 1);
		return instances;
	}

	public Instance instance(int row, Instances dataset) {
		fill(row, scratch);
		double[] out = new double[activeColumns.length + 1];
		for (int k = 0; k < activeColumns.length; k++) {
			out[k] = scratch[activeColumns[k]];
		}
		out[activeColumns.length] = candidates.get(row).isTarget() ? 0.0 : 1.0;
		DenseInstance instance = new DenseInstance(1.0, out);
		instance.setDataset(dataset);
		return instance;
	}

	public Instances subset(int[] rows, String relationName) {
		return subset(rows, relationName, null);
	}

	/**
	 * As {@link #subset(int[], String)}, with each row's class from {@code asTarget} instead of the
	 * candidate's label when given: RESET trains estimating decoys as pseudo-targets.
	 */
	public Instances subset(int[] rows, String relationName, boolean[] asTarget) {
		Instances dataset = emptyDataset(relationName, rows.length);
		for (int row : rows) {
			Instance instance = instance(row, dataset);
			if (asTarget != null) {
				instance.setClassValue(asTarget[row] ? 0.0 : 1.0);
			}
			dataset.add(instance);
		}
		return dataset;
	}

	/** Active feature values of one row, NaN for missing; for the feature report. */
	public double[] values(int row) {
		fill(row, scratch);
		double[] out = new double[activeColumns.length];
		for (int k = 0; k < activeColumns.length; k++) {
			double v = scratch[activeColumns[k]];
			out[k] = Utils.isMissingValue(v) ? Double.NaN : v;
		}
		return out;
	}

	// ------------------------------------------------------------ column choice

	private void chooseColumns(java.util.Set<String> withheld) {
		int width = ALL_FEATURES.length;
		double[] first = new double[width];
		boolean[] seen = new boolean[width];
		boolean[] varies = new boolean[width];
		for (int row = 0; row < candidates.size(); row++) {
			fill(row, scratch);
			for (int k = 0; k < width; k++) {
				double v = scratch[k];
				if (Utils.isMissingValue(v) || varies[k]) {
					continue;
				}
				if (!seen[k]) {
					seen[k] = true;
					first[k] = v;
				} else if (v != first[k]) {
					varies[k] = true;
				}
			}
		}
		List<Integer> cols = new ArrayList<>();
		List<String> names = new ArrayList<>();
		for (int k = 0; k < width; k++) {
			if (!varies[k] || (withheld != null && withheld.contains(ALL_FEATURES[k]))) {
				continue;
			}
			cols.add(k);
			names.add(ALL_FEATURES[k]);
		}
		activeColumns = new int[cols.size()];
		for (int i = 0; i < activeColumns.length; i++) {
			activeColumns[i] = cols.get(i);
		}
		activeNames = names.toArray(new String[0]);
	}

	// ------------------------------------------------------------ context

	/**
	 * Per-spectrum competition, from the best and second-best value per spectrum: "the best other
	 * candidate" is the second best for the spectrum's top candidate and the best for every other.
	 */
	private void computeContext() {
		int n = candidates.size();
		int spectra = 0;
		for (int id : spectrumId) {
			spectra = Math.max(spectra, id + 1);
		}
		Top2 closedExpect = new Top2(spectra);
		Top2 hyper = new Top2(spectra);
		Top2 entropy = new Top2(spectra);
		boolean[] hasClosed = new boolean[spectra];
		for (int i = 0; i < n; i++) {
			JointCandidate c = candidates.get(i);
			int s = spectrumId[i];
			if (c.closed != null) {
				hasClosed[s] = true;
				closedExpect.offer(s, i, negLog10(c.closed.getExpect()));
			}
			hyper.offer(s, i, msfraggerHyperscore(c));
			if (ms2 != null) {
				entropy.offer(s, i, ms2.entropy[i]);
			}
		}
		for (int i = 0; i < n; i++) {
			int s = spectrumId[i];
			spectrumHasClosed[i] = hasClosed[s] ? 1f : 0f;
			closedRival[i] = (float) closedExpect.bestOther(s, i);
			double own = msfraggerHyperscore(candidates.get(i));
			double other = hyper.bestOther(s, i);
			hyperscoreMargin[i] = Double.isNaN(own) || Double.isNaN(other) ? Float.NaN : (float) (own - other);
			if (ms2 != null) {
				double e = ms2.entropy[i];
				double o = entropy.bestOther(s, i);
				ms2Margin[i] = Double.isNaN(e) || Double.isNaN(o) ? Float.NaN : (float) (e - o);
			} else {
				ms2Margin[i] = Float.NaN;
			}
		}
	}

	/** Best and second-best value per spectrum, NaN ignored. */
	private static final class Top2 {
		final double[] best;
		final double[] second;
		final int[] bestIndex;

		Top2(int spectra) {
			best = new double[spectra];
			second = new double[spectra];
			bestIndex = new int[spectra];
			Arrays.fill(best, Double.NaN);
			Arrays.fill(second, Double.NaN);
			Arrays.fill(bestIndex, -1);
		}

		void offer(int s, int index, double value) {
			if (Double.isNaN(value)) {
				return;
			}
			if (Double.isNaN(best[s]) || value > best[s]) {
				second[s] = best[s];
				best[s] = value;
				bestIndex[s] = index;
			} else if (Double.isNaN(second[s]) || value > second[s]) {
				second[s] = value;
			}
		}

		double bestOther(int s, int index) {
			return bestIndex[s] == index ? second[s] : best[s];
		}
	}

	/** Best MSFragger hyperscore of a candidate over the closed and open searches (one scoring function). */
	static double msfraggerHyperscore(JointCandidate c) {
		double best = Double.NaN;
		for (OpenSearchPsm psm : new OpenSearchPsm[] { c.closed, c.fragger }) {
			if (psm != null && !Double.isNaN(psm.getHyperscore())) {
				best = Double.isNaN(best) ? psm.getHyperscore() : Math.max(best, psm.getHyperscore());
			}
		}
		return best;
	}

	// ------------------------------------------------------------ fill

	private void fill(int row, double[] out) {
		JointCandidate c = candidates.get(row);
		OpenSearchPsm shared = c.representative();
		OpenSearchPsm closed = c.closed;
		OpenSearchPsm fragger = c.fragger;
		OpenSearchPsm sage = c.sage;
		OpenSearchPsm peak = c.openPrimary();
		int i = 0;

		out[i++] = shared.getCharge();
		out[i++] = orMissing(shared.getPrecursorNeutralMass());
		out[i++] = shared.getPeptideLength();
		out[i++] = shared.getNumTolTerm() < 0 ? Utils.missingValue() : shared.getNumTolTerm();
		out[i++] = shared.getNumMissedCleavages() < 0 ? Utils.missingValue() : shared.getNumMissedCleavages();

		out[i++] = closed != null ? 1.0 : 0.0;
		out[i++] = fragger != null ? 1.0 : 0.0;
		out[i++] = sage != null ? 1.0 : 0.0;

		if (closed != null) {
			out[i++] = closed.getHitRank();
			out[i++] = closed.getNumMatchedIons() < 0 ? Utils.missingValue() : closed.getNumMatchedIons();
			out[i++] = orMissing(closed.getMatchedIonFraction());
			out[i++] = orMissing(closed.getHyperscore());
			out[i++] = orMissing(subtract(closed.getHyperscore(), closed.getNextscore()));
			out[i++] = orMissing(negLog10(closed.getExpect()));
			out[i++] = orMissing(Math.abs(closed.getNormalizedMassDiffPpm()));
			out[i++] = closed.getPeakIsotopeOffset();
			out[i++] = closed.getClosedModificationCount();
		} else {
			i = missing(out, i, CLOSED.length);
		}

		if (fragger != null) {
			out[i++] = fragger.getHitRank();
			out[i++] = fragger.getNumMatchedIons() < 0 ? Utils.missingValue() : fragger.getNumMatchedIons();
			out[i++] = orMissing(fragger.getMatchedIonFraction());
			out[i++] = orMissing(fragger.getHyperscore());
			out[i++] = orMissing(subtract(fragger.getHyperscore(), fragger.getNextscore()));
			out[i++] = orMissing(negLog10(fragger.getExpect()));
			out[i++] = orMissing(fragger.getDeltaMassScoreGain());
			out[i++] = orMissing(fragger.getLocalizationDeltaScore());
		} else {
			i = missing(out, i, OPEN.length);
		}

		if (peak != null) {
			out[i++] = orMissing(Math.abs(peak.getNormalizedMassDiff()));
			out[i++] = onZeroFamily(peak) ? 1.0 : 0.0;
			out[i++] = peak.isInValidatedPeak() && !onZeroFamily(peak) ? 1.0 : 0.0;
			out[i++] = peak.isInValidatedPeak() ? 1.0 : 0.0;
			out[i++] = orMissing(peak.getDistanceInSigmas());
			out[i++] = Math.log1p(Math.max(0.0, peak.getLocalDensity()));
			out[i++] = Math.log1p(Math.max(0, peak.getSamePeakCount()));
			out[i++] = peak.getPeakIsotopeOffset();
		} else {
			i = missing(out, i, PEAK.length);
		}

		if (sage != null) {
			out[i++] = sage.getHitRank();
			out[i++] = orMissing(sage.getMatchedIonFraction());
			out[i++] = orMissing(sage.getHyperscore());
			out[i++] = orMissing(subtract(sage.getHyperscore(), sage.getNextscore()));
			out[i++] = orMissing(negLog10(sage.getExpect()));
		} else {
			i = missing(out, i, SAGE.length);
		}

		out[i++] = spectrumHasClosed[row];
		out[i++] = orMissing(closedRival[row]);
		out[i++] = orMissing(hyperscoreMargin[row]);

		if (ms2 != null) {
			out[i++] = orMissing(ms2.entropy[row]);
			out[i++] = orMissing(ms2.matchedFraction[row]);
			out[i++] = orMissing(ms2Margin[row]);
			out[i] = orMissing(ms2.rtErrorNorm[row]);
		} else {
			missing(out, i, MS2.length);
		}
	}

	/** Whether an open row sits on an accepted zero-family peak. */
	private boolean onZeroFamily(OpenSearchPsm psm) {
		return psm.isInValidatedPeak() && zeroPeakIds.contains(psm.getPeakId());
	}

	private static int missing(double[] out, int i, int width) {
		for (int k = 0; k < width; k++) {
			out[i + k] = Utils.missingValue();
		}
		return i + width;
	}

	static int[] assignSpectrumIds(List<JointCandidate> candidates) {
		Map<String, Integer> ids = new HashMap<>();
		int[] out = new int[candidates.size()];
		for (int i = 0; i < out.length; i++) {
			// One id per query: the spectrum, or (identification v2) its reassigned precursor.
			String key = candidates.get(i).queryKey();
			Integer id = ids.get(key);
			if (id == null) {
				id = ids.size();
				ids.put(key, id);
			}
			out[i] = id;
		}
		return out;
	}

	private static double orMissing(double value) {
		return Double.isNaN(value) || Double.isInfinite(value) ? Utils.missingValue() : value;
	}

	private static double subtract(double a, double b) {
		return Double.isNaN(a) || Double.isNaN(b) ? Double.NaN : a - b;
	}

	static double negLog10(double expect) {
		if (Double.isNaN(expect)) {
			return Double.NaN;
		}
		return -Math.log10(Math.max(expect, 1e-300));
	}

	private static String[] concat(String[]... blocks) {
		List<String> all = new ArrayList<>();
		for (String[] block : blocks) {
			all.addAll(Arrays.asList(block));
		}
		return all.toArray(new String[0]);
	}
}
