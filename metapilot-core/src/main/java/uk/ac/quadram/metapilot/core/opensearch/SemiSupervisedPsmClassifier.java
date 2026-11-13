package uk.ac.quadram.metapilot.core.opensearch;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import weka.classifiers.Classifier;
import weka.classifiers.functions.Logistic;
import weka.classifiers.meta.FilteredClassifier;
import weka.filters.unsupervised.attribute.Standardize;
import weka.core.Instance;
import weka.core.Instances;

/**
 * Semi-supervised, cross-validated discriminant scoring of open-search PSMs.
 * <p>
 * One classifier scores <b>every</b> PSM — modified and unmodified alike. The
 * original MetaLab implementation routed unmodified PSMs around the classifier
 * through a separate decoy-ratio scoring path; the statistical premise behind
 * that (targets concentrate at zero while decoys spread across the whole
 * delta-mass axis) is sound and is preserved, but it is expressed here as
 * classifier <em>features</em> — peak membership, distance to centroid, local
 * density — rather than as a bypass.
 *
 * <h2>Why not train on "all targets = positive"</h2>
 * Targets are a mixture of correct and incorrect identifications, so labelling
 * them all positive teaches the model to accept the incorrect ones. Instead,
 * confident positives are re-seeded each round from targets passing a q-value
 * threshold, with decoys as negatives — the Percolator approach.
 *
 * <h2>Four defects in {@code core.ml.SemiSuperClassifier} this avoids</h2>
 * <ol>
 * <li><b>Inverted bootstrap.</b> That implementation seeds initial positives
 * from <em>high</em>-scoring targets but then, in its iteration, re-seeds them
 * from targets with <em>low</em> classifier scores
 * ({@code scores[i] <= scoreThreshold1}, clamped to 0.05–0.1). Training on
 * low-scoring targets as positives degrades discrimination. Here positives are
 * always the confident end.</li>
 * <li><b>Non-determinism.</b> It assigns folds with {@code Math.random()}. Here
 * folds are a deterministic function of the spectrum key, so reruns agree.</li>
 * <li><b>Evaluation on training data.</b> It scores every PSM with models that
 * were trained on that same PSM. Here the split is cross-validated: a PSM is
 * only ever scored by a model that never saw it.</li>
 * <li><b>A non-ensemble ensemble.</b> It adds every decoy to every fold, so its
 * five models share almost identical negatives. Here each fold trains on a
 * genuinely disjoint held-out partition.</li>
 * </ol>
 *
 * <h2>Splitting by spectrum, not by PSM</h2>
 * With {@code output_report_topN > 1} a spectrum contributes several competing
 * candidates. Splitting per PSM would put a spectrum's rank-1 candidate in the
 * training set and its rank-2 candidate in the test set — near-duplicate
 * features, leaking the answer across the fold boundary. Folds are therefore
 * keyed on {@code fileName:scan}.
 */
public final class SemiSupervisedPsmClassifier {

	private static final Logger LOGGER = LogManager.getLogger(SemiSupervisedPsmClassifier.class);

	/** Builds a fresh, untrained classifier. Swap here to compare model families. */
	public interface ModelFactory {
		Classifier create();

		String describe();
	}

	/**
	 * Logistic regression: the field standard for PSM rescoring (Percolator uses
	 * a linear SVM). The features are already well engineered, probabilities come
	 * out calibrated, and correlated inputs are handled gracefully — unlike naive
	 * Bayes, whose independence assumption is badly violated here, since
	 * hyperscore, nextscore, expect and matched-ion counts all covary strongly.
	 */
	public static final ModelFactory LOGISTIC_REGRESSION = logisticRegression(1e-8);

	/**
	 * {@link #LOGISTIC_REGRESSION} with its own ridge.
	 * <p>
	 * 1e-8 is safe only while the positives overlap the decoys. Positives seeded
	 * by a q-value cut sit above a score threshold with about 1% decoys among
	 * them, so with a few hundred positives almost nothing overlaps and the fit
	 * is unbounded: on PXD066732 stage 1 (2026-11) the coefficients reached
	 * thousands, signs flipped between folds, and held-out decoys tied with
	 * targets at a probability of exactly 1.0.
	 */
	public static ModelFactory logisticRegression(double ridge) {
		return new ModelFactory() {
			@Override
			public Classifier create() {
				Logistic logistic = new Logistic();
				logistic.setRidge(ridge);
				// Standardise every attribute to zero mean and unit variance before
				// fitting. weka's Logistic does not do this itself, and with ridge at
				// 1e-8 there is effectively nothing bounding the coefficients, so one
				// wide-ranging feature can drive the linear predictor past the point
				// where the sigmoid underflows: on 2026-09-12 a MSFragger+Sage run
				// returned a probability of EXACTLY 0.0 for 74.0% of candidates
				// (4,391,627 rows, including every one of the 4,383,711 outside any
				// peak), collapsing target/decoy AUC in that group from 0.784 to
				// 0.5000 - a column of perfect ties. The feature scales that provoked
				// it are real: local_density alone spanned 0 to 289,614.
				//
				// This is a conditioning change, not a modelling one. Standardising
				// inputs to a linear model leaves the fitted decision function
				// unchanged in exact arithmetic; it only stops the optimiser and the
				// sigmoid being asked to work at scales where double precision runs
				// out.
				FilteredClassifier filtered = new FilteredClassifier();
				filtered.setClassifier(logistic);
				filtered.setFilter(new Standardize());
				return filtered;
			}

			@Override
			public String describe() {
				return "logistic regression (weka.classifiers.functions.Logistic, standardised inputs"
						+ (ridge == 1e-8 ? "" : ", ridge " + ridge) + ")";
			}
		};
	}

	public static final class Config {
		public int folds = 3;
		public int maxIterations = 10;
		/** q-value defining a confident positive during training. */
		public double trainQValue = 0.01;

		/**
		 * Optional extra condition a target must meet to be seeded as a confident
		 * positive. Negatives are unaffected - every decoy is still a negative.
		 * <p>
		 * Used to require that a positive sit inside a peak the held-out gate
		 * confirmed, so the model learns what a correct <em>and explainable</em>
		 * PSM looks like rather than what a high-scoring one looks like. Training
		 * still runs over every candidate, which is what keeps the peak-context
		 * features varying and leaves non-peak PSMs with a calibrated score.
		 */
		public java.util.function.Predicate<OpenSearchPsm> positiveFilter = null;
		/** Stop when the accepted-target count changes by less than this fraction. */
		public double convergenceTolerance = 0.005;
		/** Minimum confident positives needed to attempt a fit. */
		public int minPositives = 50;

		/**
		 * Seed confident positives within each peak-membership stratum rather than
		 * from one pooled ranking.
		 * <p>
		 * Pooled seeding produced a training set in which <b>no</b> positive lay
		 * outside a peak. Measured on the 2026-09-13 MSFragger-only run:
		 * <pre>
		 *   POSITIVE (target, q&lt;=0.01)   791,939   accepted_peak 98.4%  rejected_peak 1.6%  no_peak 0.0%
		 *   NEGATIVE (decoy winners)     180,916   accepted_peak 18.7%  rejected_peak 17.9% no_peak 63.4%
		 * </pre>
		 * Every no-peak training row was therefore a negative - 114,639 of them,
		 * not one positive - so {@code in_peak == 0} predicted "decoy" in the
		 * training data without exception. The loop is self-reinforcing: a no-peak
		 * target has to outrank in-peak targets to be seeded, which it cannot do
		 * under a model trained only on in-peak positives.
		 * <p>
		 * Stratifying computes the seeding q-value inside each stratum, so each
		 * contributes its own top slice at the same {@link #trainQValue}.
		 * <p>
		 * <b>Default off: measured, and it makes the problem worse.</b> On a
		 * synthetic set shaped like the real run it gave the no-peak stratum FEWER
		 * positives, not more, and its within-group AUC fell:
		 * <pre>
		 *   stratify  negatives        no-peak AUC   no_peak positives/negatives
		 *   false     DECOY_WINNERS       0.9564       +206 / -640
		 *   false     ALL_DECOYS          0.9636       +170 / -6,000
		 *   true      DECOY_WINNERS       0.5528         +0 / -639
		 *   true      ALL_DECOYS          0.7287        +12 / -6,000
		 * </pre>
		 * The reason is straightforward once seen: the no-peak stratum is roughly
		 * a third decoy, so computing the q-value inside it makes the FDR estimate
		 * much STRICTER than pooling, where the same targets ride on the low
		 * overall decoy count. A stratum with no confident members by target-decoy
		 * standards cannot be made to yield confident positives by thresholding it
		 * separately.
		 * <p>
		 * Kept as a switch because the negative result is worth preserving - the
		 * idea is an obvious one to have twice.
		 */
		public boolean stratifyPositiveSeeding = false;

		/** Where training negatives come from. */
		public enum NegativeSource {
			/** Decoys that won their spectrum - one per spectrum, the hardest decoys. */
			DECOY_WINNERS,
			/** Every decoy candidate, as Percolator and mokapot use. */
			ALL_DECOYS
		}

		/**
		 * Negatives from every decoy candidate rather than decoy winners only.
		 * <p>
		 * Winner-only negatives are a biased and much smaller sample: on the
		 * 2026-09-13 MSFragger-only run there were 1,924,256 decoy candidates
		 * against 180,916 decoy winners, and the resulting training set ran at
		 * 4.38 positives per negative - backwards for this problem. Using every
		 * decoy gives 0.41:1.
		 * <p>
		 * This does <b>not</b> touch the seeding q-value, which is still computed
		 * over spectrum winners: target-decoy FDR is only defined per spectrum, so
		 * the set used to <em>seed</em> positives and the set used to <em>train</em>
		 * on negatives are deliberately different.
		 */
		public NegativeSource negativeSource = NegativeSource.ALL_DECOYS;


		/**
		 * Upper bound on the training set handed to the model, per fold.
		 * <p>
		 * weka's {@code Logistic} runs {@code ReplaceMissingValues} over the
		 * training set with {@code Filter.useFilter}, which materialises a full
		 * second copy, and then builds its own {@code double[][]} on top - so a
		 * training set costs roughly three times its own size in heap while
		 * fitting. On the 2026-09-13 combined run that was where the
		 * {@code OutOfMemoryError} landed.
		 * <p>
		 * Beyond a few hundred thousand rows a logistic regression over ~48
		 * features gains almost nothing from more data, so the excess is paid
		 * entirely in memory. The subsample is stratified - the positive:negative
		 * ratio is preserved - and seeded, so a rerun trains on the same rows.
		 * Zero or negative disables the cap.
		 */
		public int maxTrainingInstances = 500_000;
		public ModelFactory modelFactory = LOGISTIC_REGRESSION;
		/**
		 * Logits from both class probabilities ({@link #logitOf}) instead of the logit of the
		 * clamped target probability ({@link #toLogit}, +-27.63). Needed where candidates of one
		 * spectrum compete on the logit: confident candidates otherwise tie at the clamp. Off by
		 * default, so the existing classifiers' scores are unchanged.
		 */
		public boolean unclampedLogits = false;
		/**
		 * Seed positives and compute q-values on the unclamped logit instead of the target
		 * probability. Probabilities round to exactly 1.0 once the linear predictor passes
		 * about 37, and a decoy inside that tie block makes the whole block one q-value.
		 * Implies unclamped logits. Off by default, so the existing classifiers' scores are
		 * unchanged.
		 */
		public boolean rankByLogit = false;
	}

	public static final class Result {
		/** Probability of being a correct identification, per input index. */
		public final double[] scores;
		/** Logit of {@link #scores} — unbounded, so better behaved for downstream tests. */
		public final double[] logits;
		/** Monotonised q-values from the cross-validated scores. */
		public final double[] qValues;
		public final int acceptedTargetsAtOnePercent;
		public final int iterationsUsed;

		Result(double[] scores, double[] logits, double[] qValues, int acceptedTargetsAtOnePercent,
				int iterationsUsed) {
			this.scores = scores;
			this.logits = logits;
			this.qValues = qValues;
			this.acceptedTargetsAtOnePercent = acceptedTargetsAtOnePercent;
			this.iterationsUsed = iterationsUsed;
		}
	}

	private SemiSupervisedPsmClassifier() {
	}

	public static Result classify(List<OpenSearchPsm> psms, Config cfg) throws Exception {
		return classify(psms, cfg, null);
	}

	/**
	 * @param progress optional sink reporting each fold and training iteration.
	 *                 Training is the longest uninterrupted computation in the
	 *                 workflow — several folds times several iterations of model
	 *                 fitting over every PSM — so without this it runs silently
	 *                 for minutes and cannot be distinguished from a hang.
	 */
	public static Result classify(List<OpenSearchPsm> psms, Config cfg, java.util.function.Consumer<String> progress)
			throws Exception {
		return classify(psms, cfg, OpenSearchFeatures.FEATURE_NAMES, progress);
	}

	/**
	 * @param featureNames the feature set to train on. Supplied explicitly so a
	 *                     mixed-engine pool can be scored per engine on features
	 *                     that engine can actually populate — see
	 *                     {@link PerEngineClassifier} for why imputing the
	 *                     missing block instead would be unsound.
	 */
	public static Result classify(List<OpenSearchPsm> psms, Config cfg, String[] featureNames,
			java.util.function.Consumer<String> progress) throws Exception {
		if (cfg == null) {
			cfg = new Config();
		}
		int n = psms.size();
		double[] scores = new double[n];
		double[] logits = new double[n];
		java.util.Arrays.fill(scores, Double.NaN);
		java.util.Arrays.fill(logits, Double.NaN);

		if (n == 0) {
			return new Result(scores, logits, new double[0], 0, 0);
		}

		boolean[] isTarget = new boolean[n];
		for (int i = 0; i < n; i++) {
			isTarget[i] = psms.get(i).isTarget();
		}

		report(progress, "Classifier: building feature table for " + n + " PSM candidates ...");
		Instances all = OpenSearchFeatures.buildInstances(psms, "open_search_psms", featureNames);
		int[] fold = assignFolds(psms, cfg.folds);
		// Used to collapse to one candidate per spectrum when seeding confident
		// positives; see trainFold. Folds are already split by the same key, so a
		// spectrum's candidates never straddle the train/held-out boundary.
		int[] spectrumId = assignSpectrumIds(psms);
		int[] stratum = new int[n];
		for (int i = 0; i < n; i++) {
			stratum[i] = stratumOf(psms.get(i));
		}

		// Bootstrap ranking: the single best raw discriminant available before any
		// model exists. -log10(expect) is what MSFragger itself ranks by.
		double[] initial = new double[n];
		for (int i = 0; i < n; i++) {
			double expect = psms.get(i).getExpect();
			initial[i] = Double.isNaN(expect) ? psms.get(i).getHyperscore()
					: -Math.log10(Math.max(expect, 1e-300));
		}

		// Which targets may be seeded as confident positives. Computed once, over
		// every candidate: the filter narrows what the model learns FROM, never
		// what it scores, so peak-context features keep varying and PSMs outside
		// every peak still get a calibrated score.
		boolean[] positiveEligible = new boolean[n];
		int eligible = 0;
		for (int i = 0; i < n; i++) {
			positiveEligible[i] = cfg.positiveFilter == null || cfg.positiveFilter.test(psms.get(i));
			if (positiveEligible[i]) {
				eligible++;
			}
		}
		if (cfg.positiveFilter != null) {
			report(progress, String.format(java.util.Locale.ROOT,
					"Classifier: %d of %d candidates (%.1f%%) are eligible to seed confident positives",
					eligible, n, 100.0 * eligible / n));
		}

		int iterationsUsed = 0;
		for (int f = 0; f < cfg.folds; f++) {
			report(progress, "Classifier: fold " + (f + 1) + "/" + cfg.folds + " - training ...");
			FoldOutcome outcome = trainFold(all, initial, isTarget, fold, spectrumId, stratum, f, cfg, progress, f,
					positiveEligible);
			iterationsUsed = Math.max(iterationsUsed, outcome.iterations);
			for (int i = 0; i < n; i++) {
				if (fold[i] == f) {
					scores[i] = outcome.scores[i];
					logits[i] = outcome.logits != null ? outcome.logits[i] : toLogit(outcome.scores[i]);
				}
			}
		}

		// NOTE: these are per-CANDIDATE figures, reported only as a diagnostic of
		// the scoring step. They are not the workflow's FDR and must not be quoted
		// as an identification count: with output_report_topN > 1 one spectrum
		// contributes several rows here. The reportable numbers come from
		// SpectrumLevelCollapse + TargetDecoyFdr on one PSM per spectrum.
		double[] ranking = cfg.rankByLogit ? logits : scores;
		double[] q = TargetDecoyFdr.qValues(ranking, isTarget);
		int accepted = TargetDecoyFdr.countAcceptedTargets(ranking, isTarget, 0.01);

		LOGGER.info("PSM classifier: " + cfg.modelFactory.describe() + ", " + cfg.folds
				+ "-fold cross-validated (split by spectrum), " + iterationsUsed + " iterations; " + accepted
				+ " target candidates at 1% candidate-level FDR out of " + n
				+ " candidates (diagnostic only - see the per-spectrum figure for the reportable count)");

		return new Result(scores, logits, q, accepted, iterationsUsed);
	}

	private static final class FoldOutcome {
		final double[] scores;
		/** Unclamped logits ({@link Config#unclampedLogits}), else null. */
		final double[] logits;
		final int iterations;

		FoldOutcome(double[] scores, int iterations) {
			this(scores, null, iterations);
		}

		FoldOutcome(double[] scores, double[] logits, int iterations) {
			this.scores = scores;
			this.logits = logits;
			this.iterations = iterations;
		}
	}

	/**
	 * Trains on every fold except {@code heldOutFold}, iterating the confident
	 * positive set, then scores all instances with the final model. Only the
	 * held-out scores are kept by the caller.
	 */
	private static FoldOutcome trainFold(Instances all, double[] initialRanking, boolean[] isTarget, int[] fold,
			int[] spectrumId, int[] stratum, int heldOutFold, Config cfg,
			java.util.function.Consumer<String> progress, int foldIndex, boolean[] positiveEligible) throws Exception {
		int n = all.numInstances();
		double[] working = initialRanking.clone();

		int previousAccepted = -1;
		int iterations = 0;
		Classifier model = null;

		for (int iteration = 0; iteration < cfg.maxIterations; iteration++) {
			iterations = iteration + 1;

			// Training examples are the spectrum winners in the TRAINING partition.
			// Collapsing here matters for the same reason it matters at reporting
			// time: the seeding q-value is a target-decoy FDR, which is only defined
			// per spectrum. Without the collapse a rank-3 target could be seeded as
			// a confident positive even though a better candidate exists for that
			// same spectrum, so the model would be taught the wrong peptide.
			// Winners are recomputed every iteration because `working` changes as
			// the model improves — a candidate the engine ranked third can take over
			// its spectrum once the classifier has learned.
			boolean[] isWinner = spectrumWinners(working, spectrumId, fold, heldOutFold);

			// Positives: confident target winners. Negatives: decoys, from either the
			// winners or the whole candidate pool. The two sets are chosen by
			// different rules on purpose - the seeding q-value is only defined per
			// spectrum, while nothing requires a decoy to have won one.
			boolean[] seeded = seedPositives(working, isTarget, isWinner, stratum, cfg);

			List<Integer> trainingRows = new ArrayList<>();
			int[] posByStratum = new int[STRATUM_NAMES.length];
			int[] negByStratum = new int[STRATUM_NAMES.length];
			for (int i = 0; i < n; i++) {
				if (isTarget[i]) {
					if (!seeded[i] || !positiveEligible[i]) {
						continue;
					}
					trainingRows.add(i);
					posByStratum[stratum[i]]++;
				} else {
					boolean useAsNegative = cfg.negativeSource == Config.NegativeSource.ALL_DECOYS
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
			for (int s = 0; s < STRATUM_NAMES.length; s++) {
				positives += posByStratum[s];
				negatives += negByStratum[s];
			}
			if (iteration == 0) {
				report(progress, "Training set [fold " + (foldIndex + 1) + "]: "
						+ describeComposition(posByStratum, negByStratum));
			}

			Instances trainingSet = new Instances(all, 0);
			for (int row : capTrainingRows(trainingRows, isTarget, cfg.maxTrainingInstances, heldOutFold,
					"fold " + (foldIndex + 1), iteration == 0 ? progress : null)) {
				trainingSet.add(all.instance(row));
			}

			if (positives < cfg.minPositives || negatives < cfg.minPositives) {
				LOGGER.warn("Fold " + heldOutFold + " iteration " + iterations + ": only " + positives
						+ " positives / " + negatives + " negatives; stopping iteration and using the last usable model."
						+ " Check that the search produced enough confident identifications.");
				break;
			}

			trainingSet.setClassIndex(all.classIndex());
			Classifier candidate = cfg.modelFactory.create();
			candidate.buildClassifier(trainingSet);
			model = candidate;

			for (int i = 0; i < n; i++) {
				working[i] = cfg.rankByLogit ? logitOf(model.distributionForInstance(all.instance(i)))
						: probabilityOfTarget(model, all.instance(i));
			}

			// Convergence is measured on the same quantity that seeds the next
			// iteration — accepted target SPECTRA, winners recomputed under the
			// updated scores — so the loop stops when the thing it optimises stops
			// moving. Counting candidates here instead would let the number drift
			// with reshuffling among a spectrum's own candidates.
			boolean[] refreshedWinners = spectrumWinners(working, spectrumId, fold, heldOutFold);
			int accepted = TargetDecoyFdr.countAcceptedTargets(subset(working, refreshedWinners),
					subset(isTarget, refreshedWinners), cfg.trainQValue);
			report(progress, String.format("Classifier: fold %d iteration %d - %d positives / %d decoys, "
					+ "%d target spectra accepted at q<=%.3f", foldIndex + 1, iterations, positives, negatives,
					accepted, cfg.trainQValue));
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
			// No usable model: fall back to the raw ranking rather than returning
			// nothing, so a weak run degrades instead of failing outright.
			LOGGER.warn("Fold " + heldOutFold + ": no model could be trained; falling back to -log10(expect) ranking");
			return new FoldOutcome(normalise(initialRanking), iterations);
		}

		reportCoefficients(model, featureNamesOf(all), "fold " + (heldOutFold + 1), progress);

		double[] finalScores = new double[n];
		double[] finalLogits = cfg.unclampedLogits || cfg.rankByLogit ? new double[n] : null;
		for (int i = 0; i < n; i++) {
			if (finalLogits == null) {
				finalScores[i] = probabilityOfTarget(model, all.instance(i));
			} else {
				double[] distribution = model.distributionForInstance(all.instance(i));
				finalScores[i] = distribution.length > 0 ? distribution[0] : Double.NaN;
				finalLogits[i] = logitOf(distribution);
			}
		}
		return new FoldOutcome(finalScores, finalLogits, iterations);
	}

	/** Attribute names of a training table, class attribute excluded. */
	static String[] featureNamesOf(Instances data) {
		String[] names = new String[Math.max(0, data.numAttributes() - 1)];
		for (int i = 0; i < names.length; i++) {
			names[i] = data.attribute(i).name();
		}
		return names;
	}

	private static double probabilityOfTarget(Classifier model, Instance instance) throws Exception {
		double[] distribution = model.distributionForInstance(instance);
		// Class attribute is declared {target, decoy}, so index 0 is target.
		return distribution.length > 0 ? distribution[0] : Double.NaN;
	}

	/**
	 * Folds keyed on the spectrum, so all candidates of one spectrum land
	 * together. Uses a stable string hash rather than {@code String.hashCode()}
	 * for reproducibility across JVMs.
	 */
	static int[] assignFolds(List<OpenSearchPsm> psms, int folds) {
		Map<String, Integer> spectrumFold = new HashMap<>();
		int[] assignment = new int[psms.size()];
		for (int i = 0; i < psms.size(); i++) {
			OpenSearchPsm psm = psms.get(i);
			String key = psm.getFileName() + ":" + psm.getScan();
			Integer f = spectrumFold.get(key);
			if (f == null) {
				f = (int) Math.floorMod(stableHash(key), folds);
				spectrumFold.put(key, f);
			}
			assignment[i] = f;
		}
		return assignment;
	}

	/** Peak-membership strata, in the order the logs report them. */
	static final String[] STRATUM_NAMES = { "accepted_peak", "rejected_peak", "no_peak" };

	/** Which stratum a PSM belongs to; the filter has already run at this point. */
	static int stratumOf(OpenSearchPsm psm) {
		if (psm.getPeakId() < 0) {
			return 2;
		}
		return psm.isInValidatedPeak() ? 0 : 1;
	}

	/**
	 * Marks the confident positives for this iteration.
	 * <p>
	 * Only spectrum winners are eligible: the seeding q-value is a target-decoy
	 * FDR, which is defined per spectrum, so a rank-3 target must not be seeded
	 * while a better candidate for the same spectrum exists. When stratified, the
	 * q-value is computed separately inside each stratum at the same threshold.
	 */
	static boolean[] seedPositives(double[] scores, boolean[] isTarget, boolean[] isWinner, int[] stratum,
			Config cfg) {
		int n = scores.length;
		boolean[] seeded = new boolean[n];
		int strata = cfg.stratifyPositiveSeeding ? STRATUM_NAMES.length : 1;
		for (int s = 0; s < strata; s++) {
			List<Integer> rows = new ArrayList<>();
			for (int i = 0; i < n; i++) {
				if (!isWinner[i]) {
					continue;
				}
				if (cfg.stratifyPositiveSeeding && stratum[i] != s) {
					continue;
				}
				rows.add(i);
			}
			if (rows.isEmpty()) {
				continue;
			}
			double[] sc = new double[rows.size()];
			boolean[] tg = new boolean[rows.size()];
			for (int k = 0; k < rows.size(); k++) {
				sc[k] = scores[rows.get(k)];
				tg[k] = isTarget[rows.get(k)];
			}
			boolean[] accepted = TargetDecoyFdr.acceptedTargets(sc, tg, cfg.trainQValue);
			for (int k = 0; k < rows.size(); k++) {
				if (tg[k] && accepted[k]) {
					seeded[rows.get(k)] = true;
				}
			}
		}
		return seeded;
	}

	/**
	 * Stratified, seeded subsample when the training rows exceed the cap; the
	 * positive:negative ratio is preserved and the same rows come back on a rerun.
	 */
	static int[] capTrainingRows(List<Integer> rows, boolean[] isTarget, int cap, int seed, String context,
			java.util.function.Consumer<String> progress) {
		if (cap <= 0 || rows.size() <= cap) {
			return toIntArray(rows);
		}
		List<Integer> positives = new ArrayList<>();
		List<Integer> negatives = new ArrayList<>();
		for (int row : rows) {
			(isTarget[row] ? positives : negatives).add(row);
		}
		double keep = cap / (double) rows.size();
		int keepPositives = Math.max(1, (int) Math.round(positives.size() * keep));
		int keepNegatives = Math.max(1, (int) Math.round(negatives.size() * keep));
		java.util.Random random = new java.util.Random(0x5EEDL + seed);
		java.util.Collections.shuffle(positives, random);
		java.util.Collections.shuffle(negatives, random);
		List<Integer> kept = new ArrayList<>(keepPositives + keepNegatives);
		kept.addAll(positives.subList(0, Math.min(keepPositives, positives.size())));
		kept.addAll(negatives.subList(0, Math.min(keepNegatives, negatives.size())));
		java.util.Collections.sort(kept);
		report(progress, String.format(java.util.Locale.ROOT,
				"Training set [%s] capped %,d -> %,d rows (%,d positives / %,d decoys, ratio preserved)",
				context, rows.size(), kept.size(), keepPositives, keepNegatives));
		return toIntArray(kept);
	}

	static int[] toIntArray(List<Integer> values) {
		int[] out = new int[values.size()];
		for (int i = 0; i < out.length; i++) {
			out[i] = values.get(i);
		}
		return out;
	}

	/** One line describing what the fold is about to train on. */
	static String describeComposition(int[] posByStratum, int[] negByStratum) {
		int pos = 0;
		int neg = 0;
		for (int s = 0; s < STRATUM_NAMES.length; s++) {
			pos += posByStratum[s];
			neg += negByStratum[s];
		}
		StringBuilder sb = new StringBuilder();
		sb.append(String.format(java.util.Locale.ROOT, "%,d positives / %,d decoys (%.2f:1)", pos, neg,
				neg == 0 ? Double.NaN : pos / (double) neg));
		for (int s = 0; s < STRATUM_NAMES.length; s++) {
			sb.append(String.format(java.util.Locale.ROOT, "; %s +%,d/-%,d", STRATUM_NAMES[s], posByStratum[s],
					negByStratum[s]));
		}
		return sb.toString();
	}

	/**
	 * Logs the fitted logistic coefficients, largest magnitude first.
	 * <p>
	 * Inputs are standardised before fitting, so these are directly comparable:
	 * a coefficient's size is how much that feature moves the linear predictor
	 * per standard deviation. Logged because the alternative is guessing. On the
	 * 2026-09-13 MSFragger-only run 84% of the PSMs outside every peak came back
	 * with a probability of exactly 0.0 - the linear predictor below -745 - and
	 * nothing in the output said which terms put it there, so two candidate
	 * explanations had to be tested against synthetic data and both were wrong.
	 * <p>
	 * Silent for model families that are not logistic regression.
	 */
	static void reportCoefficients(Classifier model, String[] featureNames, String context,
			java.util.function.Consumer<String> progress) {
		Logistic logistic = null;
		if (model instanceof Logistic) {
			logistic = (Logistic) model;
		} else if (model instanceof FilteredClassifier) {
			Classifier inner = ((FilteredClassifier) model).getClassifier();
			if (inner instanceof Logistic) {
				logistic = (Logistic) inner;
			}
		}
		if (logistic == null) {
			return;
		}
		double[][] par;
		try {
			par = logistic.coefficients();
		} catch (Exception e) {
			return;
		}
		if (par == null || par.length == 0) {
			return;
		}
		// Row 0 is the intercept; rows 1..n follow the filtered attribute order,
		// which Standardize leaves unchanged. Column 0 is the only one for a
		// two-class problem.
		Integer[] order = new Integer[Math.max(0, par.length - 1)];
		for (int i = 0; i < order.length; i++) {
			order[i] = i + 1;
		}
		java.util.Arrays.sort(order, (a, b) -> Double.compare(Math.abs(par[b][0]), Math.abs(par[a][0])));
		double total = 0.0;
		double largest = 0.0;
		for (int i = 1; i < par.length; i++) {
			total += Math.abs(par[i][0]);
			largest = Math.max(largest, Math.abs(par[i][0]));
		}
		StringBuilder top = new StringBuilder();
		int shown = Math.min(10, order.length);
		for (int k = 0; k < shown; k++) {
			int attribute = order[k];
			String name = attribute - 1 < featureNames.length ? featureNames[attribute - 1]
					: "attribute_" + attribute;
			if (top.length() > 0) {
				top.append(", ");
			}
			top.append(String.format(java.util.Locale.ROOT, "%s %+.3f", name, par[attribute][0]));
		}
		report(progress, String.format(java.util.Locale.ROOT,
				"Coefficients [%s]: intercept %+.3f, sum|coef| %.1f, largest |coef| %.3f; top %d: %s",
				context, par[0][0], total, largest, shown, top));
	}

	private static void report(java.util.function.Consumer<String> progress, String message) {
		LOGGER.info(message);
		if (progress != null) {
			progress.accept(message);
		}
	}

	/** FNV-1a: deterministic across JVM versions, unlike String.hashCode(). */
	static long stableHash(String value) {
		long hash = 0xcbf29ce484222325L;
		for (int i = 0; i < value.length(); i++) {
			hash ^= value.charAt(i);
			hash *= 0x100000001b3L;
		}
		return hash;
	}

	/**
	 * Dense per-PSM spectrum identifier, keyed the same way as
	 * {@link #assignFolds}, so all candidates of one spectrum share an id and sit
	 * in one fold.
	 */
	static int[] assignSpectrumIds(List<OpenSearchPsm> psms) {
		Map<String, Integer> ids = new HashMap<>();
		int[] out = new int[psms.size()];
		for (int i = 0; i < psms.size(); i++) {
			OpenSearchPsm psm = psms.get(i);
			String key = psm.getFileName() + ":" + psm.getScan();
			Integer id = ids.get(key);
			if (id == null) {
				id = ids.size();
				ids.put(key, id);
			}
			out[i] = id;
		}
		return out;
	}

	/**
	 * Marks the highest-scoring candidate of each spectrum within the training
	 * partition (everything except {@code heldOutFold}). Held-out PSMs are never
	 * marked, so the returned mask doubles as the training-partition selector.
	 */
	static boolean[] spectrumWinners(double[] scores, int[] spectrumId, int[] fold, int heldOutFold) {
		int maxId = -1;
		for (int i = 0; i < spectrumId.length; i++) {
			maxId = Math.max(maxId, spectrumId[i]);
		}
		int[] bestIndex = new int[maxId + 1];
		java.util.Arrays.fill(bestIndex, -1);

		for (int i = 0; i < scores.length; i++) {
			if (fold[i] == heldOutFold) {
				continue;
			}
			int id = spectrumId[i];
			int incumbent = bestIndex[id];
			if (incumbent < 0 || rank(scores[i]) > rank(scores[incumbent])) {
				bestIndex[id] = i;
			}
		}

		boolean[] winner = new boolean[scores.length];
		for (int id = 0; id <= maxId; id++) {
			if (bestIndex[id] >= 0) {
				winner[bestIndex[id]] = true;
			}
		}
		return winner;
	}

	/** NaN sorts below every real score rather than losing every comparison. */
	private static double rank(double score) {
		return Double.isNaN(score) ? Double.NEGATIVE_INFINITY : score;
	}

	static double[] subset(double[] values, boolean[] keep) {
		int count = 0;
		for (boolean k : keep) {
			if (k) {
				count++;
			}
		}
		double[] out = new double[count];
		int cursor = 0;
		for (int i = 0; i < values.length; i++) {
			if (keep[i]) {
				out[cursor++] = values[i];
			}
		}
		return out;
	}

	static boolean[] subset(boolean[] values, boolean[] keep) {
		int count = 0;
		for (boolean k : keep) {
			if (k) {
				count++;
			}
		}
		boolean[] out = new boolean[count];
		int cursor = 0;
		for (int i = 0; i < values.length; i++) {
			if (keep[i]) {
				out[cursor++] = values[i];
			}
		}
		return out;
	}


	static double[] normalise(double[] raw) {
		double min = Double.MAX_VALUE;
		double max = -Double.MAX_VALUE;
		for (double v : raw) {
			if (Double.isNaN(v)) {
				continue;
			}
			min = Math.min(min, v);
			max = Math.max(max, v);
		}
		double[] out = new double[raw.length];
		double span = max - min;
		for (int i = 0; i < raw.length; i++) {
			out[i] = Double.isNaN(raw[i]) || span <= 0 ? 0.5 : (raw[i] - min) / span;
		}
		return out;
	}

	/**
	 * Probability to logit. Downstream per-peak testing compares distributions,
	 * and a bounded [0,1] probability piles up at both ends; the logit is
	 * unbounded and far closer to symmetric for a linear model.
	 */
	static double toLogit(double probability) {
		if (Double.isNaN(probability)) {
			return Double.NaN;
		}
		double clamped = Math.min(1 - 1e-12, Math.max(1e-12, probability));
		return Math.log(clamped / (1 - clamped));
	}

	/** Below {@code exp(-LOGIT_LIMIT)} a probability is 0 in double precision. */
	static final double LOGIT_LIMIT = 745.0;

	/**
	 * The log-odds of a two-class distribution {target, decoy}, without the clamp of
	 * {@link #toLogit}.
	 * <p>
	 * weka's {@code Logistic} computes each class probability separately, so the
	 * decoy probability keeps full relative precision (2.4e-303 at a linear
	 * predictor of 697) while the target probability has long since rounded to
	 * exactly 1.0. {@code log(p_target) - log(p_decoy)} is therefore the linear
	 * predictor itself up to about +-700. Recovering the logit from one rounded
	 * probability instead tied 251,164 of 1,132,941 closed winners at +27.63 on the
	 * 12mix cascade run (2026-10-20), which a target-decoy walk then breaks by row
	 * order.
	 */
	static double logitOf(double[] distribution) {
		if (distribution == null || distribution.length < 2 || Double.isNaN(distribution[0])
				|| Double.isNaN(distribution[1])) {
			return distribution != null && distribution.length == 1 ? toLogit(distribution[0]) : Double.NaN;
		}
		double target = distribution[0];
		double decoy = distribution[1];
		if (target <= 0 && decoy <= 0) {
			return Double.NaN;
		}
		if (decoy <= 0) {
			return LOGIT_LIMIT;
		}
		if (target <= 0) {
			return -LOGIT_LIMIT;
		}
		return Math.max(-LOGIT_LIMIT, Math.min(LOGIT_LIMIT, Math.log(target) - Math.log(decoy)));
	}

	/** Inverse of the logit, for reporting a probability beside it. */
	static double sigmoid(double logit) {
		if (Double.isNaN(logit)) {
			return Double.NaN;
		}
		return logit >= 0 ? 1.0 / (1.0 + Math.exp(-logit)) : Math.exp(logit) / (1.0 + Math.exp(logit));
	}
}
