package uk.ac.quadram.metapilot.core.opensearch;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Consumer;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import weka.classifiers.Classifier;
import weka.core.Instance;
import weka.core.Instances;

/**
 * Semi-supervised rescoring over {@link CombinedPsmCandidate} - one instance per
 * (spectrum, peptide), carrying both engines' evidence - replacing
 * {@link PerEngineClassifier}'s train-per-engine-then-calibrate arrangement.
 *
 * <h2>What this fixes</h2>
 * Scoring per engine requires mapping two score scales onto one before the
 * spectrum-level collapse can compare them. {@code PerEngineClassifier} does
 * that with {@code -log10(engine-local q)}, which is comparable across engines
 * only if each engine's q-value is accurate. Measured on the 2026-09-10
 * combined run it is not: at the same calibrated score, Sage rows are about
 * 1.8x more likely to be decoys.
 * <pre>
 *   band   MSFragger targets/decoys      Sage targets/decoys
 *   1-3       58,715 / 3,202  5.17%       37,548 / 3,927  9.47%
 *   3-6      349,200 /   196  0.06%       57,472 /   155  0.27%
 * </pre>
 * Sage's classifier sees 25 features against MSFragger's 33 - it has no
 * counterpart for the {@code localize_delta_mass} block - so the q-value it
 * reports is better than the q-value it delivers, and the pooled walk then
 * trusts that number. One model over one feature space removes both the
 * engine-local q-value and the scale that had to be reconciled.
 *
 * <h2>Memory</h2>
 * The candidate matrix is held by {@link CombinedPsmFeatures} as float and
 * instances are built on demand, so only the training subset - spectrum winners
 * in the training partition - exists as Weka objects at any time. The full
 * table is never materialised as {@link Instances}.
 * <p>
 * One deliberate deviation from {@link SemiSupervisedPsmClassifier}: during the
 * iteration loop only the training partition is scored, because
 * {@code spectrumWinners} discards held-out rows and the convergence count is
 * computed from the same winners. Held-out rows are scored once, by the final
 * model. This is exactly equivalent and saves one fold's worth of scoring per
 * iteration.
 */
public final class CombinedPsmClassifier {

	private static final Logger LOGGER = LogManager.getLogger(CombinedPsmClassifier.class);

	public static final class Result {
		/** Probability of being a correct identification, per candidate index. */
		public final double[] scores;
		/** Logit of {@link #scores}. */
		public final double[] logits;
		/** Monotonised q-values over the cross-validated scores. */
		public final double[] qValues;
		public final int acceptedTargetsAtOnePercent;
		public final int iterationsUsed;
		public final int candidateCount;
		public final int bothEngineCount;

		Result(double[] scores, double[] logits, double[] qValues, int acceptedTargetsAtOnePercent,
				int iterationsUsed, int candidateCount, int bothEngineCount) {
			this.scores = scores;
			this.logits = logits;
			this.qValues = qValues;
			this.acceptedTargetsAtOnePercent = acceptedTargetsAtOnePercent;
			this.iterationsUsed = iterationsUsed;
			this.candidateCount = candidateCount;
			this.bothEngineCount = bothEngineCount;
		}
	}

	private CombinedPsmClassifier() {
	}

	public static Result classify(List<CombinedPsmCandidate> candidates, SemiSupervisedPsmClassifier.Config cfg,
			Consumer<String> progress) throws Exception {
		if (cfg == null) {
			cfg = new SemiSupervisedPsmClassifier.Config();
		}
		int n = candidates.size();
		double[] scores = new double[n];
		double[] logits = new double[n];
		java.util.Arrays.fill(scores, Double.NaN);
		java.util.Arrays.fill(logits, Double.NaN);
		if (n == 0) {
			return new Result(scores, logits, new double[0], 0, 0, 0, 0);
		}

		boolean[] isTarget = new boolean[n];
		boolean[] positiveEligible = new boolean[n];
		int both = 0;
		int eligible = 0;
		for (int i = 0; i < n; i++) {
			CombinedPsmCandidate candidate = candidates.get(i);
			isTarget[i] = candidate.isTarget();
			if (candidate.foundByBoth()) {
				both++;
			}
			// The filter is defined over a PSM (restrictPositivesToValidatedPeaks
			// reads in_validated_peak), which is a property of the candidate rather
			// than of the search that found it, so the shared row answers for both.
			positiveEligible[i] = cfg.positiveFilter == null || cfg.positiveFilter.test(candidate.primary());
			if (positiveEligible[i]) {
				eligible++;
			}
		}

		report(progress, String.format(java.util.Locale.ROOT,
				"Combined classifier: %,d (spectrum, peptide) candidates; %,d (%.1f%%) proposed by both engines",
				n, both, 100.0 * both / n));
		if (cfg.positiveFilter != null) {
			report(progress, String.format(java.util.Locale.ROOT,
					"Combined classifier: %,d of %,d candidates (%.1f%%) eligible to seed confident positives",
					eligible, n, 100.0 * eligible / n));
		}

		// A feature that also decides the training label is withheld from the model.
		String[] featureNames = OpenSearchFeatures.withoutLabelDefiningFeatures(CombinedPsmFeatures.FEATURE_NAMES,
				cfg.positiveFilter != null);
		if (featureNames.length != CombinedPsmFeatures.FEATURE_NAMES.length) {
			report(progress, "Combined classifier: withholding " + OpenSearchFeatures.LABEL_DEFINING_FEATURE
					+ " - it selects the training positives, so using it as an input as well collapses the model"
					+ " onto that one column");
		}
		report(progress, "Combined classifier: scoring on " + featureNames.length
				+ " features (computed per instance, not stored) ...");
		CombinedPsmFeatures table = CombinedPsmFeatures.of(candidates, featureNames);
		Instances header = table.emptyDataset("combined_psm_candidates", 0);

		int[] fold = assignFolds(candidates, cfg.folds);
		int[] spectrumId = assignSpectrumIds(candidates);
		int[] stratum = new int[n];
		for (int i = 0; i < n; i++) {
			stratum[i] = SemiSupervisedPsmClassifier.stratumOf(candidates.get(i).primary());
		}
		double[] initial = initialRanking(candidates);

		int iterationsUsed = 0;
		for (int f = 0; f < cfg.folds; f++) {
			report(progress, "Combined classifier: fold " + (f + 1) + "/" + cfg.folds + " - training ...");
			FoldOutcome outcome = trainFold(table, header, initial, isTarget, fold, spectrumId, stratum, f, cfg,
					positiveEligible, progress);
			iterationsUsed = Math.max(iterationsUsed, outcome.iterations);
			for (int i = 0; i < n; i++) {
				if (fold[i] == f) {
					scores[i] = outcome.scores[i];
					logits[i] = SemiSupervisedPsmClassifier.toLogit(outcome.scores[i]);
				}
			}
		}

		// Per-CANDIDATE diagnostic only. The reportable count comes from
		// SpectrumLevelCollapse + TargetDecoyFdr on one PSM per spectrum.
		double[] q = TargetDecoyFdr.qValues(scores, isTarget);
		int accepted = TargetDecoyFdr.countAcceptedTargets(scores, isTarget, 0.01);

		LOGGER.info("Combined PSM classifier: " + cfg.modelFactory.describe() + ", " + cfg.folds
				+ "-fold cross-validated (split by spectrum), " + iterationsUsed + " iterations; " + accepted
				+ " target candidates at 1% candidate-level FDR out of " + n
				+ " combined candidates (diagnostic only)");

		return new Result(scores, logits, q, accepted, iterationsUsed, n, both);
	}

	private static final class FoldOutcome {
		final double[] scores;
		final int iterations;

		FoldOutcome(double[] scores, int iterations) {
			this.scores = scores;
			this.iterations = iterations;
		}
	}

	private static FoldOutcome trainFold(CombinedPsmFeatures table, Instances header, double[] initialRanking,
			boolean[] isTarget, int[] fold, int[] spectrumId, int[] stratum, int heldOutFold,
			SemiSupervisedPsmClassifier.Config cfg, boolean[] positiveEligible, Consumer<String> progress)
			throws Exception {
		int n = table.size();
		double[] working = initialRanking.clone();

		int previousAccepted = -1;
		int iterations = 0;
		Classifier model = null;

		for (int iteration = 0; iteration < cfg.maxIterations; iteration++) {
			iterations = iteration + 1;

			boolean[] isWinner = SemiSupervisedPsmClassifier.spectrumWinners(working, spectrumId, fold, heldOutFold);
			// Same rules as the single-engine path, from the same helpers, so the
			// two cannot drift apart.
			boolean[] seeded = SemiSupervisedPsmClassifier.seedPositives(working, isTarget, isWinner, stratum, cfg);

			List<Integer> trainingRows = new ArrayList<>();
			int[] posByStratum = new int[SemiSupervisedPsmClassifier.STRATUM_NAMES.length];
			int[] negByStratum = new int[SemiSupervisedPsmClassifier.STRATUM_NAMES.length];
			for (int i = 0; i < n; i++) {
				if (isTarget[i]) {
					if (!seeded[i] || !positiveEligible[i]) {
						continue;
					}
					trainingRows.add(i);
					posByStratum[stratum[i]]++;
				} else {
					boolean useAsNegative = cfg.negativeSource == SemiSupervisedPsmClassifier.Config.NegativeSource.ALL_DECOYS
							? fold[i] != heldOutFold
							: isWinner[i];
					if (!useAsNegative) {
						continue;
					}
					trainingRows.add(i);
					negByStratum[stratum[i]]++;
				}
			}
			int positives = 0;
			int negatives = 0;
			for (int s = 0; s < SemiSupervisedPsmClassifier.STRATUM_NAMES.length; s++) {
				positives += posByStratum[s];
				negatives += negByStratum[s];
			}
			if (iteration == 0) {
				report(progress, "Training set [combined fold " + (heldOutFold + 1) + "]: "
						+ SemiSupervisedPsmClassifier.describeComposition(posByStratum, negByStratum));
			}

			if (positives < cfg.minPositives || negatives < cfg.minPositives) {
				LOGGER.warn("Combined fold " + heldOutFold + " iteration " + iterations + ": only " + positives
						+ " positives / " + negatives + " negatives; stopping iteration and using the last usable"
						+ " model. Check that the search produced enough confident identifications.");
				break;
			}

			// Materialised here and nowhere else. Capped because weka's Logistic
			// holds about three copies of whatever it is given while fitting.
			int[] trainingIndices = SemiSupervisedPsmClassifier.capTrainingRows(trainingRows, isTarget,
					cfg.maxTrainingInstances, heldOutFold, "combined fold " + (heldOutFold + 1),
					iteration == 0 ? progress : null);
			Instances trainingSet = table.subset(trainingIndices, "combined_training_fold");
			Classifier fitted = cfg.modelFactory.create();
			fitted.buildClassifier(trainingSet);
			model = fitted;
			trainingSet = null;

			// Only the training partition is scored here - see the class javadoc.
			for (int i = 0; i < n; i++) {
				if (fold[i] == heldOutFold) {
					continue;
				}
				working[i] = probabilityOfTarget(model, table.instance(i, header));
			}

			boolean[] refreshedWinners = SemiSupervisedPsmClassifier.spectrumWinners(working, spectrumId, fold,
					heldOutFold);
			int accepted = TargetDecoyFdr.countAcceptedTargets(
					SemiSupervisedPsmClassifier.subset(working, refreshedWinners),
					SemiSupervisedPsmClassifier.subset(isTarget, refreshedWinners), cfg.trainQValue);
			report(progress, String.format(java.util.Locale.ROOT,
					"Combined classifier: fold %d iteration %d - %d positives / %d decoys, %d target spectra"
							+ " accepted at q<=%.3f",
					heldOutFold + 1, iterations, positives, negatives, accepted, cfg.trainQValue));
			if (previousAccepted > 0) {
				double change = Math.abs(accepted - previousAccepted) / (double) previousAccepted;
				if (change < cfg.convergenceTolerance) {
					previousAccepted = accepted;
					break;
				}
			}
			previousAccepted = accepted;
		}

		if (model == null) {
			LOGGER.warn("Combined fold " + heldOutFold
					+ ": no model could be trained; falling back to the engine ranking");
			return new FoldOutcome(SemiSupervisedPsmClassifier.normalise(initialRanking), iterations);
		}

		SemiSupervisedPsmClassifier.reportCoefficients(model, table.featureNames(),
				"combined fold " + (heldOutFold + 1), progress);

		double[] finalScores = new double[n];
		for (int i = 0; i < n; i++) {
			finalScores[i] = probabilityOfTarget(model, table.instance(i, header));
		}
		return new FoldOutcome(finalScores, iterations);
	}

	private static double probabilityOfTarget(Classifier model, Instance instance) throws Exception {
		double[] distribution = model.distributionForInstance(instance);
		// Class attribute is declared {target, decoy}, so index 0 is target.
		return distribution.length > 0 ? distribution[0] : Double.NaN;
	}

	/**
	 * Bootstrap ranking before any model exists: the best evidence either engine
	 * offers for this candidate. Taking the maximum rather than one engine's
	 * value keeps a candidate only Sage proposed from starting at the bottom for
	 * a reason that has nothing to do with its quality.
	 */
	private static double[] initialRanking(List<CombinedPsmCandidate> candidates) {
		double[] initial = new double[candidates.size()];
		for (int i = 0; i < initial.length; i++) {
			CombinedPsmCandidate candidate = candidates.get(i);
			double best = Double.NEGATIVE_INFINITY;
			best = Math.max(best, rankOf(candidate.msfragger()));
			best = Math.max(best, rankOf(candidate.sage()));
			initial[i] = best == Double.NEGATIVE_INFINITY ? 0.0 : best;
		}
		return initial;
	}

	private static double rankOf(OpenSearchPsm psm) {
		if (psm == null) {
			return Double.NEGATIVE_INFINITY;
		}
		double expect = psm.getExpect();
		if (!Double.isNaN(expect)) {
			return -Math.log10(Math.max(expect, 1e-300));
		}
		double hyperscore = psm.getHyperscore();
		return Double.isNaN(hyperscore) ? Double.NEGATIVE_INFINITY : hyperscore;
	}

	/** Folds keyed on the spectrum, so a spectrum's candidates never straddle the split. */
	private static int[] assignFolds(List<CombinedPsmCandidate> candidates, int folds) {
		Map<String, Integer> spectrumFold = new HashMap<>();
		int[] assignment = new int[candidates.size()];
		for (int i = 0; i < candidates.size(); i++) {
			String key = candidates.get(i).spectrumKey();
			Integer f = spectrumFold.get(key);
			if (f == null) {
				f = (int) Math.floorMod(SemiSupervisedPsmClassifier.stableHash(key), folds);
				spectrumFold.put(key, f);
			}
			assignment[i] = f;
		}
		return assignment;
	}

	/** Dense per-spectrum ids, in order of first appearance. */
	private static int[] assignSpectrumIds(List<CombinedPsmCandidate> candidates) {
		Map<String, Integer> ids = new HashMap<>();
		int[] assignment = new int[candidates.size()];
		for (int i = 0; i < candidates.size(); i++) {
			String key = candidates.get(i).spectrumKey();
			Integer id = ids.get(key);
			if (id == null) {
				id = ids.size();
				ids.put(key, id);
			}
			assignment[i] = id;
		}
		return assignment;
	}

	private static int[] toArray(List<Integer> values) {
		int[] out = new int[values.size()];
		for (int i = 0; i < out.length; i++) {
			out[i] = values.get(i);
		}
		return out;
	}

	private static void report(Consumer<String> progress, String message) {
		LOGGER.info(message);
		if (progress != null) {
			progress.accept(message);
		}
	}
}
