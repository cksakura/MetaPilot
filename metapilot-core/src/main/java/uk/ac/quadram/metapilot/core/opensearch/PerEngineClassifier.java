package uk.ac.quadram.metapilot.core.opensearch;

import java.util.ArrayList;
import java.util.EnumMap;
import java.util.List;
import java.util.Map;
import java.util.function.Consumer;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * Scores a mixed-engine PSM pool by training one classifier per engine and
 * calibrating the results onto a common scale.
 *
 * <h2>Why not one classifier over everything</h2>
 * Sage cannot supply MSFragger's eight delta-mass localization features. Weka's
 * logistic regression has no native handling for missing values, so a pooled
 * table would mean-impute them — giving <em>every</em> Sage row the same value
 * for those eight features. That is a perfect marker of engine identity, and
 * because the engines differ in target:decoy behaviour the model could learn
 * "which engine" as a proxy for "how good". The resulting scores would look
 * reasonable and be untrustworthy, and would defeat the very question the
 * two-engine option exists to answer.
 * <p>
 * Each engine is therefore scored on features it can actually supply:
 * MSFragger on the full set, Sage on
 * {@link OpenSearchFeatures#SHARED_FEATURE_NAMES}.
 *
 * <h2>Why the scores need calibrating before they are pooled</h2>
 * Two classifiers trained on different feature sets produce logits on different
 * scales. Comparing them directly during the spectrum-level collapse would let
 * the engine with the wider score range win contests it should lose — an
 * artefact of scale, not evidence.
 * <p>
 * Each PSM therefore also gets a calibrated score
 * ({@link OpenSearchPsm#getCalibratedScore()}) that means the same
 * thing for both engines: {@code -log10(q)} of that PSM's own engine-local
 * target-decoy q-value. A q-value is a probability statement about the
 * identification, computed against that engine's own decoys, so it is
 * comparable across engines by construction in a way that a raw logit is not.
 * Ties within an engine are broken by the underlying logit, preserving the
 * classifier's ordering inside each engine.
 */
public final class PerEngineClassifier {

	private static final Logger LOGGER = LogManager.getLogger(PerEngineClassifier.class);

	/** Score assigned when q is 0; corresponds to q = 1e-12. */
	private static final double MAX_CALIBRATED_SCORE = 12.0;

	private PerEngineClassifier() {
	}

	public static final class Result {
		/** Calibrated, cross-engine-comparable score per input index. */
		public final double[] scores;
		/** Engine-local q-value per input index. */
		public final double[] qValues;
		public final Map<OpenSearchPsm.SearchEngine, Integer> countByEngine;
		public final Map<OpenSearchPsm.SearchEngine, Integer> acceptedByEngine;

		Result(double[] scores, double[] qValues, Map<OpenSearchPsm.SearchEngine, Integer> countByEngine,
				Map<OpenSearchPsm.SearchEngine, Integer> acceptedByEngine) {
			this.scores = scores;
			this.qValues = qValues;
			this.countByEngine = countByEngine;
			this.acceptedByEngine = acceptedByEngine;
		}
	}

	/**
	 * Trains and scores per engine, writing each PSM's calibrated score, logit
	 * and probability onto it. The class logit stays the classifier's own:
	 * genome scoring weights peptides by it, and the calibrated score, at the
	 * cap for most accepted PSMs, would leave real genomes indistinguishable
	 * from the decoy background.
	 * <p>
	 * With a single engine present this is equivalent to the previous behaviour
	 * apart from the calibration step, so a MSFragger-only run is unaffected in
	 * its ordering.
	 */
	public static Result classify(List<OpenSearchPsm> psms, SemiSupervisedPsmClassifier.Config cfg,
			Consumer<String> progress) throws Exception {
		int n = psms.size();
		double[] scores = new double[n];
		double[] qValues = new double[n];
		java.util.Arrays.fill(scores, Double.NaN);
		java.util.Arrays.fill(qValues, 1.0);

		Map<OpenSearchPsm.SearchEngine, List<Integer>> byEngine = new EnumMap<>(OpenSearchPsm.SearchEngine.class);
		for (int i = 0; i < n; i++) {
			byEngine.computeIfAbsent(psms.get(i).getEngine(), k -> new ArrayList<>()).add(i);
		}

		Map<OpenSearchPsm.SearchEngine, Integer> counts = new EnumMap<>(OpenSearchPsm.SearchEngine.class);
		Map<OpenSearchPsm.SearchEngine, Integer> accepted = new EnumMap<>(OpenSearchPsm.SearchEngine.class);

		for (Map.Entry<OpenSearchPsm.SearchEngine, List<Integer>> entry : byEngine.entrySet()) {
			OpenSearchPsm.SearchEngine engine = entry.getKey();
			List<Integer> idx = entry.getValue();
			counts.put(engine, idx.size());

			List<OpenSearchPsm> subset = new ArrayList<>(idx.size());
			for (int i : idx) {
				subset.add(psms.get(i));
			}
			String[] features = OpenSearchFeatures.withoutLabelDefiningFeatures(
					OpenSearchFeatures.featureNamesFor(engine), cfg.positiveFilter != null);
			report(progress, String.format("Classifier [%s]: %,d candidates, %d features", engine, subset.size(),
					features.length));

			SemiSupervisedPsmClassifier.Result scored = SemiSupervisedPsmClassifier.classify(subset, cfg, features,
					progress);

			double[] logits = scored.logits;
			boolean[] isTarget = new boolean[subset.size()];
			for (int k = 0; k < subset.size(); k++) {
				isTarget[k] = subset.get(k).isTarget();
			}
			// Engine-local q-values: computed against this engine's own decoys,
			// which is what makes the calibrated score comparable across engines.
			double[] q = TargetDecoyFdr.qValues(logits, isTarget);

			int acc = 0;
			for (int k = 0; k < subset.size(); k++) {
				int i = idx.get(k);
				qValues[i] = q[k];
				scores[i] = calibrate(q[k], logits[k]);
				subset.get(k).setCalibratedScore(scores[i]);
				subset.get(k).setClassLogit(logits[k]);
				subset.get(k).setClassScore(scored.scores[k]);
				if (isTarget[k] && q[k] <= 0.01) {
					acc++;
				}
			}
			accepted.put(engine, acc);
			report(progress, String.format("Classifier [%s]: %,d target candidates at q<=0.01 (engine-local)", engine,
					acc));
		}

		if (byEngine.size() > 1) {
			LOGGER.info("Per-engine scoring complete; scores calibrated to -log10(engine-local q) so the "
					+ "spectrum-level collapse compares like with like across {} engines", byEngine.size());
		}
		return new Result(scores, qValues, counts, accepted);
	}

	/**
	 * Maps an engine-local q-value onto a common scale, with the classifier's own
	 * logit as a tiebreak so ordering within an engine is preserved where many
	 * PSMs share a q-value.
	 */
	private static double calibrate(double q, double logit) {
		double base = q <= 0 ? MAX_CALIBRATED_SCORE : Math.min(MAX_CALIBRATED_SCORE, -Math.log10(q));
		double tiebreak = Double.isNaN(logit) ? 0 : Math.tanh(logit / 50.0) * 1e-3;
		return base + tiebreak;
	}

	private static void report(Consumer<String> progress, String message) {
		LOGGER.info(message);
		if (progress != null) {
			progress.accept(message);
		}
	}
}
