package uk.ac.quadram.metapilot.core.opensearch;

import java.util.ArrayList;
import java.util.List;

import weka.core.Attribute;
import weka.core.DenseInstance;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.Utils;

/**
 * Feature view over {@link CombinedPsmCandidate}: one row per (spectrum,
 * peptide), with a shared block, one block per engine, and explicit presence
 * indicators where an engine did not propose the candidate.
 *
 * <h2>Nothing is stored</h2>
 * Feature values are computed when an instance is requested and discarded with
 * it. An earlier version kept the whole matrix as float, which was the wrong
 * trade: measured on real object shapes it cost 199 bytes per row - 1.99 GB at
 * 10.0M candidates - on top of a {@code List<OpenSearchPsm>} that is already
 * 462 bytes per row, or 5.62 GB at 12.17M rows. Adding two gigabytes to save
 * arithmetic is backwards when the heap is the binding constraint, and it was
 * what pushed the 2026-09-13 combined run into {@code OutOfMemoryError}.
 * <p>
 * The cost is recomputing {@link #fill} once per instance built. That is about
 * forty arithmetic operations against a
 * {@link weka.classifiers.Classifier#distributionForInstance} call that
 * dominates it by orders of magnitude.
 *
 * <h2>Missing values</h2>
 * An absent engine block is written as {@link Utils#missingValue()}, never
 * imputed: a mean-filled block would assert a score the engine never reported.
 * Weka's {@code Logistic} replaces missing values internally, using means taken
 * from the training set, which is the right place for that decision.
 *
 * <h2>Withholding a feature</h2>
 * Every feature is always computed; {@code activeColumns} decides which reach
 * the model. One fill path means a withheld feature cannot silently become an
 * uncomputed one. See
 * {@link OpenSearchFeatures#withoutLabelDefiningFeatures}.
 *
 * <p><b>Not thread safe.</b> One scratch buffer is reused across
 * {@link #instance} calls.
 */
public final class CombinedPsmFeatures {

	/** Shared: properties of the candidate, identical whichever engine found it. */
	private static final String[] SHARED = {
			"charge", "precursor_mass", "peptide_length", "num_tol_term", "num_missed_cleavages",
			"mass_diff", "mass_diff_ppm", "normalized_mass_diff", "normalized_mass_diff_ppm",
			"in_peak", "in_validated_peak", "distance_to_peak_centroid", "distance_in_sigmas",
			"peak_sigma", "peak_amplitude", "peak_isotope_offset", "log_local_density",
			"log_same_peptide_count", "log_same_protein_count", "log_same_peak_count" };

	/**
	 * Engine agreement. These two indicators encode the whole agreement state -
	 * both set means both engines proposed this exact peptide for this exact
	 * spectrum - and they are the reason this class exists.
	 */
	private static final String[] PRESENCE = { "found_by_msfragger", "found_by_sage" };

	/** MSFragger's block, including the localize_delta_mass columns Sage has no counterpart for. */
	private static final String[] MSFRAGGER = {
			"msfragger_hit_rank", "msfragger_num_matched_ions", "msfragger_tot_num_ions",
			"msfragger_matched_ion_fraction", "msfragger_hyperscore", "msfragger_nextscore",
			"msfragger_hyperscore_minus_nextscore", "msfragger_neg_log10_expect",
			"msfragger_score_with_delta_mass", "msfragger_score_without_delta_mass",
			"msfragger_delta_mass_score_gain", "msfragger_second_best_score_with_delta_mass",
			"msfragger_localization_delta_score", "msfragger_ions_with_delta_mass",
			"msfragger_ions_without_delta_mass", "msfragger_localization_site_count" };

	/** Sage's block. Deliberately not merged with MSFragger's: a hyperscore and a Sage score are not the same quantity and must never share a column. */
	private static final String[] SAGE = {
			"sage_hit_rank", "sage_num_matched_ions", "sage_tot_num_ions", "sage_matched_ion_fraction",
			"sage_hyperscore", "sage_nextscore", "sage_hyperscore_minus_nextscore", "sage_neg_log10_expect" };

	/** Cross-engine, missing unless both engines proposed the candidate. */
	private static final String[] CROSS = { "rank_difference", "best_rank" };

	/**
	 * DIA-NN predicted-spectrum and retention-time agreement ({@link PredictedSpectrumFeatures#attach}),
	 * read from the primary row's external features; missing when not computed. The classifier uses
	 * this block only when some candidate carries it.
	 */
	static final String[] PREDICTED = PredictedSpectrumFeatures.PSM_FEATURES;

	public static final String[] FEATURE_NAMES = concat(SHARED, PRESENCE, MSFRAGGER, SAGE, CROSS, PREDICTED);

	public static final String CLASS_ATTRIBUTE = "label";
	public static final String TARGET_LABEL = "target";
	public static final String DECOY_LABEL = "decoy";

	private final List<CombinedPsmCandidate> candidates;
	private final int[] activeColumns;
	private final String[] activeNames;
	private final double[] scratch = new double[FEATURE_NAMES.length];

	private CombinedPsmFeatures(List<CombinedPsmCandidate> candidates, int[] activeColumns, String[] activeNames) {
		this.candidates = candidates;
		this.activeColumns = activeColumns;
		this.activeNames = activeNames;
	}

	/** View over every feature. */
	public static CombinedPsmFeatures of(List<CombinedPsmCandidate> candidates) {
		return of(candidates, FEATURE_NAMES);
	}

	/** View exposing only {@code activeFeatureNames}, in that order. */
	public static CombinedPsmFeatures of(List<CombinedPsmCandidate> candidates, String[] activeFeatureNames) {
		int[] active = new int[activeFeatureNames.length];
		for (int k = 0; k < activeFeatureNames.length; k++) {
			int found = -1;
			for (int j = 0; j < FEATURE_NAMES.length; j++) {
				if (FEATURE_NAMES[j].equals(activeFeatureNames[k])) {
					found = j;
					break;
				}
			}
			if (found < 0) {
				throw new IllegalArgumentException("Not a combined feature: " + activeFeatureNames[k]);
			}
			active[k] = found;
		}
		return new CombinedPsmFeatures(candidates, active, activeFeatureNames.clone());
	}

	public int size() {
		return candidates.size();
	}

	public boolean isTarget(int row) {
		return candidates.get(row).isTarget();
	}

	/** The feature names this view exposes. */
	public String[] featureNames() {
		return activeNames.clone();
	}

	public ArrayList<Attribute> attributes() {
		ArrayList<Attribute> attributes = new ArrayList<>(activeNames.length + 1);
		for (String name : activeNames) {
			attributes.add(new Attribute(name));
		}
		attributes.add(new Attribute(CLASS_ATTRIBUTE, java.util.Arrays.asList(TARGET_LABEL, DECOY_LABEL)));
		return attributes;
	}

	/** An empty {@link Instances} with this view's feature set, ready to be filled. */
	public Instances emptyDataset(String relationName, int capacity) {
		ArrayList<Attribute> attributes = attributes();
		Instances instances = new Instances(relationName, attributes, Math.max(0, capacity));
		instances.setClassIndex(attributes.size() - 1);
		return instances;
	}

	/** Computes one row's instance. Nothing is retained; the caller owns the result. */
	public Instance instance(int row, Instances dataset) {
		fill(candidates.get(row), scratch);
		double[] out = new double[activeColumns.length + 1];
		for (int k = 0; k < activeColumns.length; k++) {
			out[k] = scratch[activeColumns[k]];
		}
		out[activeColumns.length] = candidates.get(row).isTarget() ? 0.0 : 1.0;
		DenseInstance instance = new DenseInstance(1.0, out);
		instance.setDataset(dataset);
		return instance;
	}

	/** Materialises the given rows as a dataset - the training set, and nothing larger. */
	public Instances subset(int[] rowIndices, String relationName) {
		Instances dataset = emptyDataset(relationName, rowIndices.length);
		for (int row : rowIndices) {
			dataset.add(instance(row, dataset));
		}
		return dataset;
	}

	// ------------------------------------------------------------- feature fill

	private static void fill(CombinedPsmCandidate candidate, double[] out) {
		OpenSearchPsm shared = candidate.primary();
		OpenSearchPsm fragger = candidate.msfragger();
		OpenSearchPsm sage = candidate.sage();
		int i = 0;

		out[i++] = shared.getCharge();
		out[i++] = shared.getPrecursorNeutralMass();
		out[i++] = shared.getPeptideLength();
		out[i++] = shared.getNumTolTerm();
		out[i++] = shared.getNumMissedCleavages();

		out[i++] = orMissing(shared.getMassDiff());
		out[i++] = orMissing(shared.getMassDiffPpm());
		out[i++] = orMissing(shared.getNormalizedMassDiff());
		out[i++] = orMissing(shared.getNormalizedMassDiffPpm());

		out[i++] = shared.isInPeak() ? 1.0 : 0.0;
		out[i++] = shared.isInValidatedPeak() ? 1.0 : 0.0;

		out[i++] = orMissing(shared.getDistanceToPeakCentroid());
		out[i++] = orMissing(shared.getDistanceInSigmas());
		out[i++] = orMissing(shared.getPeakSigma());
		out[i++] = orMissing(shared.getPeakAmplitude());
		out[i++] = shared.getPeakIsotopeOffset();
		out[i++] = Math.log1p(Math.max(0.0, shared.getLocalDensity()));

		out[i++] = Math.log1p(Math.max(0, shared.getSamePeptideCount()));
		out[i++] = Math.log1p(Math.max(0, shared.getSameProteinCount()));
		out[i++] = Math.log1p(Math.max(0, shared.getSamePeakCount()));

		out[i++] = fragger != null ? 1.0 : 0.0;
		out[i++] = sage != null ? 1.0 : 0.0;

		i = fillEngine(fragger, out, i, true);
		i = fillEngine(sage, out, i, false);

		if (fragger != null && sage != null) {
			out[i++] = Math.abs(fragger.getHitRank() - sage.getHitRank());
			out[i++] = Math.min(fragger.getHitRank(), sage.getHitRank());
		} else {
			// Not "no difference" - undefined. A zero here would tell the model the
			// two engines agreed perfectly on rank for a candidate only one of them
			// ever saw.
			out[i++] = Utils.missingValue();
			out[i++] = fragger != null ? fragger.getHitRank() : sage.getHitRank();
		}

		for (String name : PREDICTED) {
			Double v = shared.getExternalFeature(name);
			out[i++] = v == null ? Utils.missingValue() : orMissing(v);
		}
	}

	/**
	 * Writes one engine's block, or a run of missing values where that engine did
	 * not propose the candidate. The MSFragger block is eight columns longer: the
	 * localize_delta_mass outputs have no Sage counterpart, and those columns are
	 * simply absent from Sage's block rather than shared and left empty.
	 */
	private static int fillEngine(OpenSearchPsm psm, double[] out, int i, boolean msfragger) {
		int width = msfragger ? MSFRAGGER.length : SAGE.length;
		if (psm == null) {
			for (int k = 0; k < width; k++) {
				out[i + k] = Utils.missingValue();
			}
			return i + width;
		}
		out[i++] = psm.getHitRank();
		out[i++] = psm.getNumMatchedIons();
		out[i++] = psm.getTotNumIons();
		out[i++] = orMissing(psm.getMatchedIonFraction());
		out[i++] = orMissing(psm.getHyperscore());
		out[i++] = orMissing(psm.getNextscore());
		out[i++] = orMissing(subtract(psm.getHyperscore(), psm.getNextscore()));
		out[i++] = orMissing(negLog10(psm.getExpect()));
		if (!msfragger) {
			return i;
		}
		out[i++] = orMissing(psm.getScoreWithDeltaMass());
		out[i++] = orMissing(psm.getScoreWithoutDeltaMass());
		out[i++] = orMissing(psm.getDeltaMassScoreGain());
		out[i++] = orMissing(psm.getSecondBestScoreWithDeltaMass());
		out[i++] = orMissing(psm.getLocalizationDeltaScore());
		out[i++] = orMissing(psm.getIonsWithDeltaMass());
		out[i++] = orMissing(psm.getIonsWithoutDeltaMass());
		int siteCount = psm.getLocalizationSiteCount();
		out[i++] = siteCount < 0 ? Utils.missingValue() : siteCount;
		return i;
	}

	private static double orMissing(double value) {
		return Double.isNaN(value) || Double.isInfinite(value) ? Utils.missingValue() : value;
	}

	private static double subtract(double a, double b) {
		return Double.isNaN(a) || Double.isNaN(b) ? Double.NaN : a - b;
	}

	private static double negLog10(double expect) {
		if (Double.isNaN(expect)) {
			return Double.NaN;
		}
		return -Math.log10(Math.max(expect, 1e-300));
	}

	private static String[] concat(String[]... blocks) {
		List<String> all = new ArrayList<>();
		for (String[] block : blocks) {
			for (String name : block) {
				all.add(name);
			}
		}
		return all.toArray(new String[0]);
	}
}
