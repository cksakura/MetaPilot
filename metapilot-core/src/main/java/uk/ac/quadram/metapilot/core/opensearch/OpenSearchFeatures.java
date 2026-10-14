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
 * Turns {@link OpenSearchPsm} candidates into WEKA instances for the validation
 * classifier.
 * <p>
 * Governing principle: <b>every condition is a feature, never a filter</b>.
 * Delta-mass peak membership, isotope status and localization ambiguity all
 * appear here as feature values — including for PSMs that fall outside every
 * fitted peak, which are scored rather than discarded.
 * <p>
 * The target/decoy label is the class, never a feature.
 */
public final class OpenSearchFeatures {

	/** Feature names, in the order they are written into each instance. */
	public static final String[] FEATURE_NAMES = {
			// identification quality
			"charge", "precursor_mass", "peptide_length", "num_tol_term", "num_missed_cleavages", "num_matched_ions",
			"tot_num_ions", "matched_ion_fraction", "hyperscore", "nextscore", "hyperscore_minus_nextscore",
			"neg_log10_expect", "hit_rank",
			// delta-mass evidence
			"mass_diff", "mass_diff_ppm", "normalized_mass_diff", "normalized_mass_diff_ppm",
			"score_with_delta_mass", "score_without_delta_mass",
			"delta_mass_score_gain", "second_best_score_with_delta_mass", "localization_delta_score",
			"ions_with_delta_mass", "ions_without_delta_mass", "localization_site_count",
			// delta-mass peak context
			"in_peak", "in_validated_peak", "distance_to_peak_centroid", "distance_in_sigmas", "peak_sigma", "peak_amplitude", "peak_r2",
			"peak_isotope_offset", "log_local_density",
			// aggregate counts (log-scaled; computed once from raw data)
			"log_same_peptide_count", "log_same_protein_count", "log_same_peak_count" };

	/**
	 * Features MSFragger alone can supply, from {@code localize_delta_mass}.
	 * <p>
	 * Sage has no counterpart for any of these, so they are excluded when
	 * building a Sage feature table. They are not merely absent — they are
	 * absent <em>for every Sage row</em>, which is why a pooled classifier
	 * cannot simply impute them: the imputed value would be a perfect marker of
	 * engine identity, and the two engines differ in target:decoy behaviour, so
	 * the model could use "which engine" as a stand-in for "how good".
	 */
	public static final String[] MSFRAGGER_ONLY_FEATURES = { "score_with_delta_mass", "score_without_delta_mass",
			"delta_mass_score_gain", "second_best_score_with_delta_mass", "localization_delta_score",
			"ions_with_delta_mass", "ions_without_delta_mass", "localization_site_count" };

	/** Feature names usable by any engine; the MSFragger-only block removed. */
	public static final String[] SHARED_FEATURE_NAMES = sharedFeatureNames();

	private static String[] sharedFeatureNames() {
		List<String> only = Arrays.asList(MSFRAGGER_ONLY_FEATURES);
		List<String> shared = new ArrayList<>(FEATURE_NAMES.length);
		for (String name : FEATURE_NAMES) {
			if (!only.contains(name)) {
				shared.add(name);
			}
		}
		return shared.toArray(new String[0]);
	}

	/**
	 * The feature whose value defines which targets may seed confident positives
	 * when {@code restrictPositivesToValidatedPeaks} is on.
	 */
	public static final String LABEL_DEFINING_FEATURE = "in_validated_peak";

	/**
	 * Removes any feature that also decides the training label.
	 * <p>
	 * A quantity used to select the positive class must not also be an input: the
	 * model can then separate its training set by reading that one column, learns
	 * nothing else, and returns a flat score for everything on the other side of
	 * it. On the 2026-09-13 MSFragger-only run that produced within-group AUCs of
	 * 0.4986 outside every peak and 0.5035 inside rejected peaks - coin flips -
	 * with the two score ranges not overlapping at all, and the
	 * identification-only walk retaining 11 PSMs of 403,393.
	 * <p>
	 * Called wherever the classifier's feature set is chosen, so the leak cannot
	 * reappear by someone switching the filter on.
	 */
	public static String[] withoutLabelDefiningFeatures(String[] featureNames, boolean restrictPositives) {
		if (!restrictPositives) {
			return featureNames;
		}
		List<String> kept = new ArrayList<>(featureNames.length);
		for (String name : featureNames) {
			if (!LABEL_DEFINING_FEATURE.equals(name)) {
				kept.add(name);
			}
		}
		return kept.toArray(new String[0]);
	}

	/** Feature set appropriate to an engine. */
	public static String[] featureNamesFor(OpenSearchPsm.SearchEngine engine) {
		return engine == OpenSearchPsm.SearchEngine.SAGE ? SHARED_FEATURE_NAMES : FEATURE_NAMES;
	}

	public static final String CLASS_ATTRIBUTE = "label";
	public static final String TARGET_LABEL = "target";
	public static final String DECOY_LABEL = "decoy";

	private OpenSearchFeatures() {
	}

	public static ArrayList<Attribute> attributes() {
		return attributes(FEATURE_NAMES);
	}

	/** Attribute list for a given feature set; the class attribute is appended. */
	public static ArrayList<Attribute> attributes(String[] featureNames) {
		ArrayList<Attribute> attributes = new ArrayList<>(featureNames.length + 1);
		for (String name : featureNames) {
			attributes.add(new Attribute(name));
		}
		attributes.add(new Attribute(CLASS_ATTRIBUTE, Arrays.asList(TARGET_LABEL, DECOY_LABEL)));
		return attributes;
	}

	/**
	 * Fills the aggregate count fields from the raw PSM list.
	 * <p>
	 * Computed <b>once, from the data</b> — deliberately not recomputed from
	 * classifier scores inside the training loop as the original implementation
	 * did. Recomputing them per round lets PSMs reinforce one another's features
	 * across iterations, inflating confidence without any new evidence.
	 */
	public static void computeAggregateCounts(List<OpenSearchPsm> psms) {
		Map<String, Integer> peptideCounts = new HashMap<>();
		Map<String, Integer> proteinCounts = new HashMap<>();
		Map<Integer, Integer> peakCounts = new HashMap<>();

		for (OpenSearchPsm psm : psms) {
			if (psm.getPeptide() != null) {
				peptideCounts.merge(psm.getPeptide(), 1, Integer::sum);
			}
			if (psm.getProtein() != null) {
				proteinCounts.merge(psm.getProtein(), 1, Integer::sum);
			}
			peakCounts.merge(psm.getPeakId(), 1, Integer::sum);
		}

		for (OpenSearchPsm psm : psms) {
			psm.setSamePeptideCount(psm.getPeptide() == null ? 0 : peptideCounts.getOrDefault(psm.getPeptide(), 0));
			psm.setSameProteinCount(psm.getProtein() == null ? 0 : proteinCounts.getOrDefault(psm.getProtein(), 0));
			psm.setSamePeakCount(peakCounts.getOrDefault(psm.getPeakId(), 0));
		}
	}

	/** Copies fitted peak context onto each PSM, by parallel index. */
	public static void applyPeakContext(List<OpenSearchPsm> psms, DeltaMassPeakFitter.Result fit) {
		if (psms.size() != fit.observedAssignment.length) {
			throw new IllegalArgumentException("PSM count does not match fit result length");
		}
		for (int i = 0; i < psms.size(); i++) {
			OpenSearchPsm psm = psms.get(i);
			psm.setObservedPeakId(fit.observedAssignment[i]);
			psm.setPeakId(fit.familyAssignment[i]);
			psm.setNormalizedMassDiff(fit.normalizedMassDiff[i]);
			psm.setDistanceToPeakCentroid(fit.distanceToNearestCentroid[i]);
			psm.setDistanceInSigmas(fit.distanceInSigmas[i]);
			psm.setLocalDensity(fit.localDensity[i]);
			DeltaMassPeakFitter.Peak peak = fit.familyPeakOf(i);
			if (peak != null) {
				psm.setPeakCentroid(peak.centroid);
				psm.setPeakIsotopeOffset(fit.isotopeOffset[i]);
				psm.setPeakSigma(peak.sigma);
				psm.setPeakAmplitude(peak.amplitude);
				psm.setPeakR2(peak.r2);
			} else {
				// Outside every peak: still scored, with peak features absent.
				psm.setPeakCentroid(Double.NaN);
				psm.setPeakIsotopeOffset(0);
				psm.setPeakSigma(Double.NaN);
				psm.setPeakAmplitude(Double.NaN);
				psm.setPeakR2(Double.NaN);
			}
		}
	}

	public static Instances buildInstances(List<OpenSearchPsm> psms, String relationName) {
		return buildInstances(psms, relationName, FEATURE_NAMES);
	}

	/**
	 * Builds a table over an explicit feature set, so each engine is scored on
	 * features it can actually supply rather than on mean-imputed placeholders.
	 */
	public static Instances buildInstances(List<OpenSearchPsm> psms, String relationName, String[] featureNames) {
		ArrayList<Attribute> attributes = attributes(featureNames);
		Instances instances = new Instances(relationName, attributes, psms.size());
		instances.setClassIndex(attributes.size() - 1);
		for (OpenSearchPsm psm : psms) {
			instances.add(toInstance(psm, instances, featureNames));
		}
		return instances;
	}

	public static Instance toInstance(OpenSearchPsm psm, Instances dataset) {
		return toInstance(psm, dataset, FEATURE_NAMES);
	}

	/**
	 * Computes the full feature vector, then projects it onto {@code featureNames}
	 * by name. Projecting rather than branching keeps one definition of how each
	 * feature is derived, so a Sage row and an MSFragger row compute any shared
	 * feature identically.
	 */
	public static Instance toInstance(OpenSearchPsm psm, Instances dataset, String[] featureNames) {
		double[] full = fullFeatureVector(psm);
		// Identity, not length. A caller-supplied array that happens to be the
		// same length as FEATURE_NAMES is not the same array, and taking the fast
		// path for it would silently feed the classifier the wrong columns under
		// the right names.
		if (featureNames == FEATURE_NAMES) {
			double[] values = Arrays.copyOf(full, FEATURE_NAMES.length + 1);
			values[FEATURE_NAMES.length] = psm.isTarget() ? 0.0 : 1.0;
			DenseInstance instance = new DenseInstance(1.0, values);
			instance.setDataset(dataset);
			return instance;
		}
		Map<String, Double> byName = new HashMap<>();
		for (int k = 0; k < FEATURE_NAMES.length; k++) {
			byName.put(FEATURE_NAMES[k], full[k]);
		}
		double[] values = new double[featureNames.length + 1];
		for (int k = 0; k < featureNames.length; k++) {
			Double v = byName.get(featureNames[k]);
			if (v == null) {
				// Not one this class computes: MSBooster's predicted-spectrum and
				// retention-time columns, and stage-1's derived ratios, arrive
				// through the PSM's external channel instead.
				v = psm.getExternalFeature(featureNames[k]);
			}
			values[k] = v == null ? Utils.missingValue() : v;
		}
		values[featureNames.length] = psm.isTarget() ? 0.0 : 1.0;
		DenseInstance instance = new DenseInstance(1.0, values);
		instance.setDataset(dataset);
		return instance;
	}

	private static double[] fullFeatureVector(OpenSearchPsm psm) {
		double[] values = new double[FEATURE_NAMES.length + 1];
		int i = 0;

		values[i++] = psm.getCharge();
		values[i++] = psm.getPrecursorNeutralMass();
		values[i++] = psm.getPeptideLength();
		values[i++] = psm.getNumTolTerm();
		values[i++] = psm.getNumMissedCleavages();
		values[i++] = psm.getNumMatchedIons();
		values[i++] = psm.getTotNumIons();
		values[i++] = orMissing(psm.getMatchedIonFraction());
		values[i++] = orMissing(psm.getHyperscore());
		values[i++] = orMissing(psm.getNextscore());
		values[i++] = orMissing(subtract(psm.getHyperscore(), psm.getNextscore()));
		values[i++] = orMissing(negLog10(psm.getExpect()));
		values[i++] = psm.getHitRank();

		values[i++] = orMissing(psm.getMassDiff());
		values[i++] = orMissing(psm.getMassDiffPpm());
		values[i++] = orMissing(psm.getNormalizedMassDiff());
		values[i++] = orMissing(psm.getNormalizedMassDiffPpm());
		values[i++] = orMissing(psm.getScoreWithDeltaMass());
		values[i++] = orMissing(psm.getScoreWithoutDeltaMass());
		values[i++] = orMissing(psm.getDeltaMassScoreGain());
		values[i++] = orMissing(psm.getSecondBestScoreWithDeltaMass());
		values[i++] = orMissing(psm.getLocalizationDeltaScore());
		values[i++] = orMissing(psm.getIonsWithDeltaMass());
		values[i++] = orMissing(psm.getIonsWithoutDeltaMass());
		int siteCount = psm.getLocalizationSiteCount();
		values[i++] = siteCount < 0 ? Utils.missingValue() : siteCount;

		values[i++] = psm.isInPeak() ? 1.0 : 0.0;
		// Membership of a peak the held-out files confirmed, as opposed to one the
		// discovery step merely proposed. Both are kept: "in a peak that failed the
		// gate" is itself informative, and collapsing the two would hide it.
		values[i++] = psm.isInValidatedPeak() ? 1.0 : 0.0;
		values[i++] = orMissing(psm.getDistanceToPeakCentroid());
		values[i++] = orMissing(psm.getDistanceInSigmas());
		values[i++] = orMissing(psm.getPeakSigma());
		values[i++] = orMissing(psm.getPeakAmplitude());
		values[i++] = orMissing(psm.getPeakR2());
		values[i++] = psm.getPeakIsotopeOffset();
		// Log-scaled, like the log_same_* counts below: a raw density spans several
		// orders of magnitude (3,879 to 289,614 mean on 12mix between the
		// unexplained and in-peak populations), which is what pushed the logistic
		// regression into saturation once a second engine doubled it.
		values[i++] = Math.log1p(Math.max(0.0, psm.getLocalDensity()));

		values[i++] = log1p(psm.getSamePeptideCount());
		values[i++] = log1p(psm.getSameProteinCount());
		values[i++] = log1p(psm.getSamePeakCount());

		values[i] = psm.isTarget() ? 0.0 : 1.0; // index of TARGET_LABEL / DECOY_LABEL
		return values;
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
		// Guard the zero/denormal tail: MSFragger reports 0.0 for very strong
		// matches, and -log10(0) would be infinite.
		double floored = Math.max(expect, 1e-300);
		return -Math.log10(floored);
	}

	private static double log1p(int count) {
		return Math.log1p(Math.max(0, count));
	}
}
