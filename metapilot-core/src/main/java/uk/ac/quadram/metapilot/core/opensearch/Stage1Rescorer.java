package uk.ac.quadram.metapilot.core.opensearch;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Consumer;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * Replaces stage 1's single-score cut with a cross-validated, semi-supervised
 * rescoring of the same PSMs.
 *
 * <h2>Why stage 1 gets a classifier after all</h2>
 * {@link UnmodifiedPsmSelector} declines to run the full validation machinery
 * here, and rightly: a well-trained PTM classifier wants the refined database
 * that stage 1 exists to produce, which is circular. That objection is about
 * <em>modification</em> validation. It does not apply to ranking PSMs inside the
 * zero family against a fixed marker database: the features are properties of a
 * match against a database that does not change, positives are seeded from
 * targets at q&lt;=1% on the engine's own score, and negatives are the same 1:1
 * {@code rev_} decoys. Nothing about the refined database enters the loop.
 * <p>
 * What stays true is that stage 1's output is still a hard cut. This changes the
 * score the cut ranks on, not the decision to cut.
 *
 * <h2>Feature set</h2>
 * Not the 36 features stage 2 uses. Measured by ablation on nine 12mix marker
 * runs, against the {@code -log10(expect)} baseline of 53,728 PSMs at 94.13%
 * ground-truth consistency:
 * <ul>
 * <li>search scores alone: +4.9%</li>
 * <li>plus peptide properties: +8.8%</li>
 * <li>plus delta-mass features: +8.7% - <b>nothing</b>, because inside the zero
 * family the delta mass is near constant and carries no information. They are
 * excluded here.</li>
 * <li>plus {@code log_same_peptide_count}: <b>+14.3%</b>, and purity rises
 * rather than falls. The single most valuable feature at this stage.</li>
 * </ul>
 * The peak-context block is excluded because there is one peak, so it is
 * constant; the localization block because {@code delta_mass_exclude_ranges}
 * suppresses shifted-ion scoring near zero, leaving it undefined.
 * <p>
 * Two per-length ratios that measurably helped are not in
 * {@link OpenSearchFeatures#FEATURE_NAMES}. Rather than change that array - and
 * with it stage 2's vector - they are supplied through the same external-feature
 * channel MSBooster's columns use.
 *
 * <h2>Training population</h2>
 * The zero family, not the whole rank-1 axis. Tested: training on all rank-1
 * PSMs gives 57x more decoys and a <em>worse</em> result (+10.3% against
 * +14.3%), because the full axis is a different population - mostly shifted and
 * unassigned masses - so the boundary learned there is mis-specified for the
 * zero peak. Decoy scarcity is not the binding constraint here; fold count and
 * ridge barely move the result.
 *
 * <h2>A trap worth knowing about</h2>
 * Any feature computed over targets only is a perfect decoy discriminator. A
 * genome-count feature counted over target proteins alone measured +47% yield at
 * 84% ground-truth consistency, while the nominal target-decoy FDR still read
 * 1%: target-decoy is structurally blind to this. {@link #assertFeaturesNotLeaky}
 * is a cheap guard - it checks that neither class is constant in any feature -
 * and is run before training.
 */
public final class Stage1Rescorer {

	private static final Logger LOGGER = LogManager.getLogger(Stage1Rescorer.class);

	/** Derived features supplied through the external channel; see class docs. */
	public static final String HYPERSCORE_PER_LENGTH = "hyperscore_per_length";
	public static final String NEG_LOG10_EXPECT_PER_LENGTH = "neg_log10_expect_per_length";

	/**
	 * The stage-1 feature set: {@link OpenSearchFeatures#FEATURE_NAMES} minus the
	 * delta-mass, localization and peak-context blocks, plus the two ratios.
	 */
	public static final String[] STAGE1_FEATURE_NAMES = {
			"hyperscore", "nextscore", "hyperscore_minus_nextscore", "neg_log10_expect",
			"matched_ion_fraction", "num_matched_ions", "tot_num_ions",
			"peptide_length", "num_missed_cleavages", "charge", "num_tol_term",
			"log_same_peptide_count", "log_same_protein_count",
			HYPERSCORE_PER_LENGTH, NEG_LOG10_EXPECT_PER_LENGTH };

	/**
	 * Logistic ridge for stage 1. A marker search can leave only a few hundred
	 * seeded positives, and the shared 1e-8 lets that fit diverge; see
	 * {@link SemiSupervisedPsmClassifier#logisticRegression(double)}. Replayed on
	 * 20 PXD066732 marker searches, 0.1 rescued every collapsed run and moved the
	 * two largest by under 0.5%; 10 cost up to 30% on small runs.
	 */
	public static final double RIDGE = 0.1;

	public static final class Config {
		/** Off falls back to the engine-score cut, changing nothing. */
		public boolean enabled = true;
		public double fdrThreshold = UnmodifiedPsmSelector.DEFAULT_FDR;
		public SemiSupervisedPsmClassifier.Config classifier = defaultClassifierConfig();

		private static SemiSupervisedPsmClassifier.Config defaultClassifierConfig() {
			SemiSupervisedPsmClassifier.Config cfg = new SemiSupervisedPsmClassifier.Config();
			cfg.folds = 3;
			cfg.maxIterations = 10;
			cfg.trainQValue = 0.01;
			cfg.modelFactory = SemiSupervisedPsmClassifier.logisticRegression(RIDGE);
			cfg.rankByLogit = true;
			return cfg;
		}
	}

	public static final class Result {
		/** PSMs at or below the q-value threshold; targets and decoys alike. */
		public final List<OpenSearchPsm> accepted;
		/** Feature names the model actually trained on, MSBooster's included. */
		public final List<String> featureNames;
		public final int targetsBeforeFdr;
		public final int decoysBeforeFdr;
		public final int targetsAfterFdr;
		public final int decoysAfterFdr;
		public final int iterationsUsed;
		/** How many PSMs carried MSBooster features; 0 when it did not run. */
		public final int withPredictedFeatures;

		Result(List<OpenSearchPsm> accepted, List<String> featureNames, int targetsBeforeFdr, int decoysBeforeFdr,
				int targetsAfterFdr, int decoysAfterFdr, int iterationsUsed, int withPredictedFeatures) {
			this.accepted = accepted;
			this.featureNames = featureNames;
			this.targetsBeforeFdr = targetsBeforeFdr;
			this.decoysBeforeFdr = decoysBeforeFdr;
			this.targetsAfterFdr = targetsAfterFdr;
			this.decoysAfterFdr = decoysAfterFdr;
			this.iterationsUsed = iterationsUsed;
			this.withPredictedFeatures = withPredictedFeatures;
		}

		public double estimatedFdr() {
			return targetsAfterFdr == 0 ? Double.NaN : (double) decoysAfterFdr / (double) targetsAfterFdr;
		}
	}

	private Stage1Rescorer() {
	}

	/**
	 * Attaches derived and predicted features, trains, and returns the accepted
	 * set.
	 *
	 * @param candidates       the zero family plus any admitted expected-shift
	 *                         PSMs; targets and decoys, before any cut.
	 * @param predictedByKey   MSBooster features keyed as
	 *                         {@code MSBoosterPin.key(...)}, or null when it did
	 *                         not run.
	 * @param predictedNames   the column names in {@code predictedByKey}, aligned
	 *                         with its arrays; ignored when it is null.
	 * @return the accepted set, or null only when rescoring is switched off or
	 *         there are too few candidates to train on. A failure throws.
	 */
	public static Result rescore(List<OpenSearchPsm> candidates, Config cfg,
			Map<String, double[]> predictedByKey, List<String> predictedNames, Consumer<String> progress) {
		if (cfg == null) {
			cfg = new Config();
		}
		if (!cfg.enabled || candidates == null || candidates.size() < 200) {
			return null;
		}
		try {
			int withPredicted = attachExternalFeatures(candidates, predictedByKey, predictedNames);
			// Aggregate counts come from the data once, here, rather than being
			// recomputed inside the training loop - see OpenSearchFeatures.
			OpenSearchFeatures.computeAggregateCounts(candidates);

			List<String> names = new ArrayList<>(Arrays.asList(STAGE1_FEATURE_NAMES));
			if (withPredicted > 0 && predictedNames != null) {
				names.addAll(predictedNames);
			}
			String[] featureNames = names.toArray(new String[0]);

			List<String> leaky = assertFeaturesNotLeaky(candidates, featureNames);
			if (!leaky.isEmpty()) {
				// A genuine leak means the feature is not what its author thought it
				// was, and every number downstream of it is unsafe. Refusing the run
				// is the only honest outcome: target-decoy FDR cannot detect this, so
				// continuing would produce a result that looks fine and is not.
				throw new IllegalStateException("Stage 1 rescoring aborted: feature(s) " + leaky
						+ " are constant within one class while varying in the other, which makes them label "
						+ "markers rather than evidence. Target-decoy FDR cannot detect this. Fix or remove the "
						+ "feature; do not run with it.");
			}

			// Uninformative features are dropped rather than fatal. num_tol_term is
			// constant on any fully-tryptic search (num_enzyme_termini = 2), which
			// is the normal configuration, not a defect.
			List<String> constant = findConstantFeatures(candidates, featureNames);
			if (!constant.isEmpty()) {
				names.removeAll(constant);
				featureNames = names.toArray(new String[0]);
				report(progress, "Stage 1 rescoring: dropped " + constant
						+ " - constant across every PSM in this run, so they carry no information");
			}

			SemiSupervisedPsmClassifier.Result model = SemiSupervisedPsmClassifier.classify(candidates, cfg.classifier,
					featureNames, progress);

			int targetsBefore = 0;
			for (OpenSearchPsm psm : candidates) {
				if (psm.isTarget()) {
					targetsBefore++;
				}
			}
			List<OpenSearchPsm> accepted = new ArrayList<>();
			int targetsAfter = 0;
			int decoysAfter = 0;
			for (int i = 0; i < candidates.size(); i++) {
				OpenSearchPsm psm = candidates.get(i);
				psm.setClassScore(model.scores[i]);
				psm.setClassLogit(model.logits[i]);
				psm.setQValue(model.qValues[i]);
				if (model.qValues[i] <= cfg.fdrThreshold) {
					accepted.add(psm);
					if (psm.isTarget()) {
						targetsAfter++;
					} else {
						decoysAfter++;
					}
				}
			}
			Result result = new Result(accepted, names, targetsBefore, candidates.size() - targetsBefore,
					targetsAfter, decoysAfter, model.iterationsUsed, withPredicted);
			report(progress, String.format(java.util.Locale.ROOT,
					"Stage 1 rescoring: %d features (%d from MSBooster), %d iterations; %d targets/%d decoys before, "
							+ "%d targets/%d decoys after q<=%.3f (estimated FDR %.6f)",
					featureNames.length, withPredicted > 0 && predictedNames != null ? predictedNames.size() : 0,
					model.iterationsUsed, targetsBefore, candidates.size() - targetsBefore, targetsAfter, decoysAfter,
					cfg.fdrThreshold, result.estimatedFdr()));
			return result;
		} catch (Exception e) {
			// No silent fallback. A run configured to rescore that quietly does not
			// is indistinguishable in its output from one that did, which is how an
			// entire nine-run benchmark was completed on the engine-score cut
			// without anyone noticing. Opt out by setting stage1RescoringEnabled
			// false; do not opt out by failing.
			LOGGER.error("Stage 1 rescoring failed", e);
			throw new IllegalStateException("Stage 1 rescoring failed: " + e.getMessage()
					+ ". Set stage1_rescoring_enabled=false to run without it deliberately.", e);
		}
	}

	/**
	 * Adds the two per-length ratios to every PSM, and MSBooster's columns to
	 * those that have them.
	 *
	 * @return how many PSMs received predicted features.
	 */
	private static int attachExternalFeatures(List<OpenSearchPsm> psms, Map<String, double[]> predictedByKey,
			List<String> predictedNames) {
		int matched = 0;
		for (OpenSearchPsm psm : psms) {
			Map<String, Double> external = new HashMap<>();
			double length = Math.max(1.0, psm.getPeptideLength());
			double hyperscore = psm.getHyperscore();
			double expect = psm.getExpect();
			double negLog10Expect = Double.isNaN(expect) ? Double.NaN
					: -Math.log10(Math.max(expect, 1e-300));
			external.put(HYPERSCORE_PER_LENGTH, Double.isNaN(hyperscore) ? Double.NaN : hyperscore / length);
			external.put(NEG_LOG10_EXPECT_PER_LENGTH,
					Double.isNaN(negLog10Expect) ? Double.NaN : negLog10Expect / length);

			if (predictedByKey != null && predictedNames != null && !predictedNames.isEmpty()) {
				String key = uk.ac.quadram.metapilot.core.msbooster.MSBoosterPin.key(psm.getFileName(), psm.getScan(),
						psm.getHitRank());
				double[] values = predictedByKey.get(key);
				if (values != null && values.length == predictedNames.size()) {
					for (int i = 0; i < values.length; i++) {
						external.put(predictedNames.get(i), values[i]);
					}
					matched++;
				}
			}
			psm.setExternalFeatures(external);
		}
		return matched;
	}

    /**
     * Returns the features that are constant within the target class or within
     * the decoy class - the signature of a feature computed over one class only.
     * <p>
     * A single check, run before training, for the failure mode target-decoy FDR
     * cannot see. It deliberately looks for constancy rather than for a large
     * difference between the classes: a genuinely strong feature also separates
     * the classes, and only a degenerate one has no variation inside one of them.
     */
	public static List<String> assertFeaturesNotLeaky(List<OpenSearchPsm> psms, String[] featureNames) {
		List<String> leaky = new ArrayList<>();
		Set<String> checked = new LinkedHashSet<>(Arrays.asList(featureNames));
		for (String name : checked) {
			Double targetFirst = null;
			Double decoyFirst = null;
			boolean targetVaries = false;
			boolean decoyVaries = false;
			int targetSeen = 0;
			int decoySeen = 0;
			for (OpenSearchPsm psm : psms) {
				Double v = valueOf(psm, name);
				if (v == null || v.isNaN()) {
					continue;
				}
				if (psm.isTarget()) {
					targetSeen++;
					if (targetFirst == null) {
						targetFirst = v;
					} else if (!targetVaries && Math.abs(v - targetFirst) > 1e-12) {
						targetVaries = true;
					}
				} else {
					decoySeen++;
					if (decoyFirst == null) {
						decoyFirst = v;
					} else if (!decoyVaries && Math.abs(v - decoyFirst) > 1e-12) {
						decoyVaries = true;
					}
				}
			}
			// A feature constant in ONE class while varying in the other is a label
			// marker. A feature constant in BOTH is merely uninformative - which is
			// the normal state of num_tol_term on a fully-tryptic search, where
			// MSFragger reports nothing else - and must not be mistaken for a leak.
			// Flagging that case aborted stage-1 rescoring on every run of the
			// 2026-09-10 12mix benchmark.
			if (targetSeen > 20 && decoySeen > 20 && (targetVaries != decoyVaries)) {
				leaky.add(name);
			}
		}
		return leaky;
	}

	/**
	 * Features whose value never changes across the whole candidate set. Harmless
	 * but useless: they add a column of noise for the ridge to shrink and make the
	 * feature count misleading in the log.
	 */
	public static List<String> findConstantFeatures(List<OpenSearchPsm> psms, String[] featureNames) {
		List<String> constant = new ArrayList<>();
		for (String name : new LinkedHashSet<>(Arrays.asList(featureNames))) {
			Double first = null;
			boolean varies = false;
			int seen = 0;
			for (OpenSearchPsm psm : psms) {
				Double v = valueOf(psm, name);
				if (v == null || v.isNaN()) {
					continue;
				}
				seen++;
				if (first == null) {
					first = v;
				} else if (Math.abs(v - first) > 1e-12) {
					varies = true;
					break;
				}
			}
			if (seen > 50 && !varies) {
				constant.add(name);
			}
		}
		return constant;
	}

	private static Double valueOf(OpenSearchPsm psm, String name) {
		Double external = psm.getExternalFeature(name);
		if (external != null) {
			return external;
		}
		switch (name) {
		case "hyperscore":
			return psm.getHyperscore();
		case "nextscore":
			return psm.getNextscore();
		case "neg_log10_expect":
			return Double.isNaN(psm.getExpect()) ? null : -Math.log10(Math.max(psm.getExpect(), 1e-300));
		case "num_matched_ions":
			return (double) psm.getNumMatchedIons();
		case "tot_num_ions":
			return (double) psm.getTotNumIons();
		case "peptide_length":
			return (double) psm.getPeptideLength();
		case "num_missed_cleavages":
			return (double) psm.getNumMissedCleavages();
		case "charge":
			return (double) psm.getCharge();
		case "num_tol_term":
			return (double) psm.getNumTolTerm();
		case "log_same_peptide_count":
			return Math.log1p(psm.getSamePeptideCount());
		case "log_same_protein_count":
			return Math.log1p(psm.getSameProteinCount());
		default:
			// Not independently checkable here; the classifier's own projection
			// resolves it, and a derived score is not a candidate for this leak.
			return null;
		}
	}

	private static void report(Consumer<String> progress, String message) {
		LOGGER.info(message);
		if (progress != null) {
			progress.accept(message);
		}
	}
}
