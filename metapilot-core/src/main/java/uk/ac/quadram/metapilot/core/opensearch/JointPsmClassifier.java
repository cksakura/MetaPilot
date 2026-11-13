package uk.ac.quadram.metapilot.core.opensearch;

import java.util.ArrayList;
import java.util.List;
import java.util.function.Consumer;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import weka.classifiers.Classifier;
import weka.core.Instance;
import weka.core.Instances;

/**
 * Semi-supervised, cross-validated scoring of {@link JointCandidate}s: one model over the closed,
 * open and predicted-spectrum evidence, so every spectrum's candidates from both searches are ranked
 * on one scale.
 * <p>
 * The loop is {@link CombinedPsmClassifier}'s (Percolator's scheme): positives re-seeded each
 * iteration from target spectrum winners at {@code trainQValue} in the training partition, every
 * decoy in the training partition a negative, folds keyed on the spectrum, held-out rows scored
 * only by a model that never saw them. The helpers are {@link SemiSupervisedPsmClassifier}'s, so
 * the three classifiers cannot drift apart.
 * <p>
 * Positives are not restricted to accepted delta-mass peaks here: the FDR groups already keep
 * unexplained shifts apart, and restricting positives to a set the features describe would let the
 * model separate the training set by reading that set's indicator.
 */
public final class JointPsmClassifier {

	private static final Logger LOGGER = LogManager.getLogger(JointPsmClassifier.class);

	public static final class Result {
		public final double[] scores;
		public final double[] logits;
		public final int iterationsUsed;
		/** Per fold, the coefficient line {@link SemiSupervisedPsmClassifier#reportCoefficients} logged. */
		public final List<String> modelDescriptions;

		Result(double[] scores, double[] logits, int iterationsUsed, List<String> modelDescriptions) {
			this.scores = scores;
			this.logits = logits;
			this.iterationsUsed = iterationsUsed;
			this.modelDescriptions = modelDescriptions;
		}
	}

	private JointPsmClassifier() {
	}

	public static Result classify(JointFeatures table, List<JointCandidate> candidates,
			SemiSupervisedPsmClassifier.Config cfg, Consumer<String> progress) throws Exception {
		return classify(table, candidates, cfg, progress, null);
	}

	/**
	 * @param trainingDecoy RESET ({@code docs/developer/identification-v2-design.md}, C2): per candidate,
	 *                      whether it is a training decoy. Training decoys are the only negatives; the
	 *                      other decoys (estimating decoys) are trained as targets (pseudo-targets), and
	 *                      positives are re-seeded among pseudo-targets against training decoys. Null:
	 *                      every decoy trains, as before.
	 */
	public static Result classify(JointFeatures table, List<JointCandidate> candidates,
			SemiSupervisedPsmClassifier.Config cfg, Consumer<String> progress, boolean[] trainingDecoy)
			throws Exception {
		if (cfg == null) {
			cfg = new SemiSupervisedPsmClassifier.Config();
		}
		int n = table.size();
		double[] scores = new double[n];
		double[] logits = new double[n];
		java.util.Arrays.fill(scores, Double.NaN);
		java.util.Arrays.fill(logits, Double.NaN);
		List<String> descriptions = new ArrayList<>();
		if (n == 0) {
			return new Result(scores, logits, 0, descriptions);
		}
		boolean[] isTarget = new boolean[n];
		int[] fold = new int[n];
		int[] stratum = new int[n];
		int folds = Math.max(2, cfg.folds);
		for (int i = 0; i < n; i++) {
			JointCandidate c = candidates.get(i);
			// Pseudo-target under RESET: targets and estimating decoys train as one class.
			isTarget[i] = c.isTarget() || (trainingDecoy != null && !trainingDecoy[i]);
			fold[i] = (int) Math.floorMod(SemiSupervisedPsmClassifier.stableHash(c.spectrumKey()), folds);
		}
		int[] spectrumId = table.spectrumIds();
		double[] initial = initialRanking(candidates);
		Instances header = table.emptyDataset("joint_psm_candidates", 0);
		report(progress, "Joint classifier: " + cfg.modelFactory.describe() + " on " + table.featureNames().length
				+ " features: " + String.join(", ", table.featureNames()));
		List<String> withheld = table.withheldColumns();
		if (!withheld.isEmpty()) {
			report(progress, "Joint classifier: withheld (no information on this data, or excluded): "
					+ String.join(", ", withheld));
		}

		int iterationsUsed = 0;
		for (int f = 0; f < folds; f++) {
			report(progress, "Joint classifier: fold " + (f + 1) + "/" + folds + " - training ...");
			double[] held = trainFold(table, header, initial, isTarget, fold, spectrumId, stratum, f, cfg, progress,
					descriptions, trainingDecoy != null);
			iterationsUsed = Math.max(iterationsUsed, (int) held[n]);
			for (int i = 0; i < n; i++) {
				if (fold[i] == f) {
					logits[i] = held[i];
					scores[i] = SemiSupervisedPsmClassifier.sigmoid(held[i]);
				}
			}
		}
		LOGGER.info("Joint PSM classifier: " + cfg.modelFactory.describe() + ", " + folds + "-fold, "
				+ iterationsUsed + " iterations over " + n + " candidates");
		return new Result(scores, logits, iterationsUsed, descriptions);
	}

	/**
	 * Logits of the held-out rows from the fold's final model (other rows NaN), with the iteration
	 * count appended at index n.
	 */
	private static double[] trainFold(JointFeatures table, Instances header, double[] initialRanking,
			boolean[] isTarget, int[] fold, int[] spectrumId, int[] stratum, int heldOutFold,
			SemiSupervisedPsmClassifier.Config cfg, Consumer<String> progress, List<String> descriptions,
			boolean relabel) throws Exception {
		int n = table.size();
		double[] working = initialRanking.clone();
		int previousAccepted = -1;
		int iterations = 0;
		Classifier model = null;
		for (int iteration = 0; iteration < cfg.maxIterations; iteration++) {
			iterations = iteration + 1;
			boolean[] isWinner = SemiSupervisedPsmClassifier.spectrumWinners(working, spectrumId, fold, heldOutFold);
			boolean[] seeded = SemiSupervisedPsmClassifier.seedPositives(working, isTarget, isWinner, stratum, cfg);
			List<Integer> rows = new ArrayList<>();
			int positives = 0;
			int negatives = 0;
			for (int i = 0; i < n; i++) {
				if (isTarget[i]) {
					if (seeded[i]) {
						rows.add(i);
						positives++;
					}
				} else {
					boolean negative = cfg.negativeSource == SemiSupervisedPsmClassifier.Config.NegativeSource.ALL_DECOYS
							? fold[i] != heldOutFold
							: isWinner[i];
					if (negative) {
						rows.add(i);
						negatives++;
					}
				}
			}
			if (positives < cfg.minPositives || negatives < cfg.minPositives) {
				LOGGER.warn("Joint fold " + heldOutFold + " iteration " + iterations + ": only " + positives
						+ " positives / " + negatives + " negatives; keeping the last usable model");
				break;
			}
			int[] training = SemiSupervisedPsmClassifier.capTrainingRows(rows, isTarget, cfg.maxTrainingInstances,
					heldOutFold, "joint fold " + (heldOutFold + 1), iteration == 0 ? progress : null);
			// isTarget is the pseudo-target label under RESET; the table's own label is the candidate's.
			Instances trainingSet = table.subset(training, "joint_training_fold", relabel ? isTarget : null);
			Classifier fitted = cfg.modelFactory.create();
			fitted.buildClassifier(trainingSet);
			model = fitted;
			trainingSet = null;
			if (model instanceof XGBoostTreeClassifier) {
				scoreInBatches((XGBoostTreeClassifier) model, table, fold, heldOutFold, false, working);
			}
			for (int i = 0; i < n && !(model instanceof XGBoostTreeClassifier); i++) {
				if (fold[i] != heldOutFold) {
					working[i] = logitOfTarget(model, table.instance(i, header));
				}
			}
			boolean[] refreshed = SemiSupervisedPsmClassifier.spectrumWinners(working, spectrumId, fold, heldOutFold);
			int accepted = TargetDecoyFdr.countAcceptedTargets(SemiSupervisedPsmClassifier.subset(working, refreshed),
					SemiSupervisedPsmClassifier.subset(isTarget, refreshed), cfg.trainQValue);
			report(progress, String.format(java.util.Locale.ROOT,
					"Joint classifier: fold %d iteration %d - %,d positives / %,d decoys, %,d target spectra at q<=%.3f",
					heldOutFold + 1, iterations, positives, negatives, accepted, cfg.trainQValue));
			if (previousAccepted > 0
					&& Math.abs(accepted - previousAccepted) / (double) previousAccepted < cfg.convergenceTolerance) {
				break;
			}
			previousAccepted = accepted;
		}
		double[] out = new double[n + 1];
		java.util.Arrays.fill(out, Double.NaN);
		if (model == null) {
			LOGGER.warn("Joint fold " + heldOutFold + ": no model could be trained; using the search ranking");
			double[] fallback = SemiSupervisedPsmClassifier.normalise(initialRanking);
			for (int i = 0; i < n; i++) {
				out[i] = SemiSupervisedPsmClassifier.toLogit(fallback[i]);
			}
		} else {
			SemiSupervisedPsmClassifier.reportCoefficients(model, table.featureNames(),
					"joint fold " + (heldOutFold + 1), line -> {
						descriptions.add(line);
						if (progress != null) {
							progress.accept(line);
						}
					});
			// Only the held-out rows are this fold's to score.
			if (model instanceof XGBoostTreeClassifier) {
				scoreInBatches((XGBoostTreeClassifier) model, table, fold, heldOutFold, true, out);
			}
			for (int i = 0; i < n && !(model instanceof XGBoostTreeClassifier); i++) {
				if (fold[i] == heldOutFold) {
					out[i] = logitOfTarget(model, table.instance(i, header));
				}
			}
		}
		out[n] = iterations;
		return out;
	}

	/**
	 * Target log-odds from both class probabilities ({@link SemiSupervisedPsmClassifier#logitOf}):
	 * continuous where the target probability alone has rounded to 1.0, so confident candidates do
	 * not tie.
	 */
	private static double logitOfTarget(Classifier model, Instance instance) throws Exception {
		return SemiSupervisedPsmClassifier.logitOf(model.distributionForInstance(instance));
	}

	/** Rows scored per XGBoost call. */
	private static final int SCORING_BATCH = 200_000;

	/**
	 * Target logits from the tree model for the rows of the held-out fold ({@code heldOut}) or of the
	 * training folds, in batches: one prediction call per batch instead of one per row.
	 */
	private static void scoreInBatches(XGBoostTreeClassifier model, JointFeatures table, int[] fold, int heldOutFold,
			boolean heldOut, double[] out) throws Exception {
		int n = table.size();
		int nf = model.numFeatures();
		float[] buffer = new float[SCORING_BATCH * nf];
		int[] rows = new int[SCORING_BATCH];
		int k = 0;
		for (int i = 0; i <= n; i++) {
			if (i < n && (fold[i] == heldOutFold) == heldOut) {
				double[] v = table.values(i);
				for (int j = 0; j < nf; j++) {
					buffer[k * nf + j] = Double.isFinite(v[j]) ? (float) v[j] : Float.NaN;
				}
				rows[k++] = i;
			}
			if (k > 0 && (k == SCORING_BATCH || i == n)) {
				double[] logits = model.targetLogits(k == SCORING_BATCH ? buffer : java.util.Arrays.copyOf(buffer, k * nf), k);
				for (int r = 0; r < k; r++) {
					out[rows[r]] = logits[r];
				}
				k = 0;
			}
		}
	}

	/** Before any model: the best -log10(expect) any search gives the candidate. */
	static double[] initialRanking(List<JointCandidate> candidates) {
		double[] initial = new double[candidates.size()];
		for (int i = 0; i < initial.length; i++) {
			JointCandidate c = candidates.get(i);
			double best = Double.NEGATIVE_INFINITY;
			for (OpenSearchPsm psm : new OpenSearchPsm[] { c.closed, c.fragger, c.sage }) {
				if (psm == null) {
					continue;
				}
				double v = JointFeatures.negLog10(psm.getExpect());
				if (Double.isNaN(v)) {
					v = psm.getHyperscore();
				}
				if (!Double.isNaN(v)) {
					best = Math.max(best, v);
				}
			}
			initial[i] = best == Double.NEGATIVE_INFINITY ? 0.0 : best;
		}
		return initial;
	}

	private static void report(Consumer<String> progress, String message) {
		LOGGER.info(message);
		if (progress != null) {
			progress.accept(message);
		}
	}
}
