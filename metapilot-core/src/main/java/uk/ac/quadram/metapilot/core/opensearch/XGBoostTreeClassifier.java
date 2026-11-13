package uk.ac.quadram.metapilot.core.opensearch;

import java.util.HashMap;
import java.util.Locale;
import java.util.Map;
import java.util.Set;

import ml.dmlc.xgboost4j.java.Booster;
import ml.dmlc.xgboost4j.java.DMatrix;
import ml.dmlc.xgboost4j.java.XGBoost;
import ml.dmlc.xgboost4j.java.XGBoostError;
import weka.classifiers.AbstractClassifier;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.Utils;

/**
 * Gradient-boosted trees (XGBoost, {@code hist}) behind weka's {@code Classifier} interface, so the
 * semi-supervised loops use it like the logistic model (identification v2, C3,
 * {@code docs/developer/identification-v2-design.md}).
 * <p>
 * Class value 0 is the target ({@link JointFeatures#TARGET_LABEL}); the booster learns
 * {@code P(target)}. Attributes named in {@code increasing} get a monotone increasing constraint, so a
 * better search score never lowers the target score. Missing values stay missing (XGBoost learns a
 * default direction for them). {@link #targetLogits} scores many rows at once; the weka per-instance
 * path is kept for small uses.
 */
public final class XGBoostTreeClassifier extends AbstractClassifier {

	private static final long serialVersionUID = 1L;

	/** Defaults as the audit harness: depth 4, 200 rounds, learning rate 0.1, L2 1. */
	public int maxDepth = 4;
	public int rounds = 200;
	public double learningRate = 0.1;
	public double l2 = 1.0;
	public int minChildWeight = 1;
	public int threads = 1;
	public long seed = 1;
	private final Set<String> increasing;
	private transient Booster booster;
	private int numFeatures;

	public XGBoostTreeClassifier(Set<String> increasing) {
		this.increasing = increasing;
	}

	@Override
	public void buildClassifier(Instances data) throws Exception {
		int classIndex = data.classIndex();
		numFeatures = data.numAttributes() - 1;
		int n = data.numInstances();
		float[] x = new float[n * numFeatures];
		float[] y = new float[n];
		for (int r = 0; r < n; r++) {
			Instance in = data.instance(r);
			fill(in, classIndex, x, r * numFeatures);
			y[r] = in.classValue() == 0.0 ? 1f : 0f;
		}
		StringBuilder constraints = new StringBuilder("(");
		for (int a = 0, k = 0; a < data.numAttributes(); a++) {
			if (a == classIndex) {
				continue;
			}
			constraints.append(k++ == 0 ? "" : ",").append(increasing.contains(data.attribute(a).name()) ? 1 : 0);
		}
		constraints.append(')');
		Map<String, Object> params = new HashMap<>();
		params.put("objective", "binary:logistic");
		params.put("tree_method", "hist");
		params.put("max_depth", maxDepth);
		params.put("eta", learningRate);
		params.put("lambda", l2);
		params.put("min_child_weight", minChildWeight);
		params.put("nthread", Math.max(1, threads));
		params.put("seed", seed);
		params.put("monotone_constraints", constraints.toString());
		params.put("verbosity", 0);
		DMatrix train = new DMatrix(x, n, numFeatures, Float.NaN);
		try {
			train.setLabel(y);
			booster = XGBoost.train(train, params, rounds, new HashMap<>(), null, null);
		} finally {
			train.dispose();
		}
	}

	private void fill(Instance in, int classIndex, float[] x, int offset) {
		for (int a = 0, k = 0; a < in.numAttributes(); a++) {
			if (a == classIndex) {
				continue;
			}
			double v = in.value(a);
			x[offset + k++] = Utils.isMissingValue(v) || !Double.isFinite(v) ? Float.NaN : (float) v;
		}
	}

	/** Target log-odds (the booster's margin) for rows given as a row-major matrix, NaN for missing. */
	public double[] targetLogits(float[] rows, int n) throws XGBoostError {
		DMatrix m = new DMatrix(rows, n, numFeatures, Float.NaN);
		try {
			float[][] margin = booster.predict(m, true);
			double[] out = new double[n];
			for (int r = 0; r < n; r++) {
				out[r] = margin[r][0];
			}
			return out;
		} finally {
			m.dispose();
		}
	}

	@Override
	public double[] distributionForInstance(Instance instance) throws Exception {
		float[] x = new float[numFeatures];
		fill(instance, instance.classIndex(), x, 0);
		double logit = targetLogits(x, 1)[0];
		double p = 1.0 / (1.0 + Math.exp(-logit));
		return new double[] { p, 1.0 - p };
	}

	public int numFeatures() {
		return numFeatures;
	}

	@Override
	public String toString() {
		return String.format(Locale.ROOT, "XGBoost hist, depth %d, %d rounds, eta %.2f, lambda %.1f, monotone on %s",
				maxDepth, rounds, learningRate, l2, increasing);
	}
}
