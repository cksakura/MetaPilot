package uk.ac.quadram.metapilot.core.unified;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.stream.IntStream;

/**
 * Gradient-boosted trees for a listwise choice with a "none" option: one query is one spectrum,
 * its rows are the spectrum's candidates, and a virtual "none" row has fixed utility 0. With
 * utilities {@code u_j}, {@code P(j) = exp(u_j) / (1 + sum_k exp(u_k))} and
 * {@code P(none) = 1 / (1 + sum_k exp(u_k))}. A query's label is either one candidate (it is
 * correct) or none (no candidate is). The fixed "none" utility anchors every spectrum on one scale.
 * <p>
 * Second-order boosting (gradient {@code p - y}, hessian {@code p (1 - p)}) on histogram-binned
 * features, depth-wise trees, L2 on leaf values, and XGBoost-style monotone constraints: a feature
 * marked +1 can only raise a candidate's utility as it increases. The output is the raw utility, a
 * continuous margin, so confident candidates do not tie.
 */
public final class ListwiseGbdt {

	public static final class Params {
		public int trees = 120;
		public int maxDepth = 5;
		public double learningRate = 0.1;
		public double lambda = 1.0;
		public double minChildHessian = 2.0;
		public double minGain = 1e-6;
		public double featureFraction = 0.8;
		public double queryFraction = 0.8;
		public int bins = 64;
		public long seed = 0x5EEDL;
	}

	/** One tree, flattened. Internal nodes: feature >= 0; leaves: feature = -1. */
	static final class Tree {
		final int[] feature;
		final double[] threshold;
		final int[] left;
		final int[] right;
		final double[] value;

		Tree(int[] feature, double[] threshold, int[] left, int[] right, double[] value) {
			this.feature = feature;
			this.threshold = threshold;
			this.left = left;
			this.right = right;
			this.value = value;
		}

		double predict(float[] x) {
			int node = 0;
			while (feature[node] >= 0) {
				float v = x[feature[node]];
				node = !(v > threshold[node]) ? left[node] : right[node];
			}
			return value[node];
		}
	}

	public static final class Model {
		final List<Tree> trees;

		Model(List<Tree> trees) {
			this.trees = trees;
		}

		public double utility(float[] x) {
			double u = 0;
			for (Tree t : trees) {
				u += t.predict(x);
			}
			return u;
		}

		public int size() {
			return trees.size();
		}
	}

	private ListwiseGbdt() {
	}

	/**
	 * @param x         training rows (candidate features)
	 * @param queryStart CSR offsets: query q owns rows {@code queryStart[q] .. queryStart[q+1]-1}
	 * @param positive  per query, the row index of the correct candidate, or -1 for "none"
	 * @param monotone  per feature: +1, -1 or 0
	 */
	public static Model train(float[][] x, int[] queryStart, int[] positive, int[] monotone, Params p) {
		int rows = x.length;
		int features = rows == 0 ? 0 : x[0].length;
		int queries = queryStart.length - 1;
		double[][] cuts = new double[features][];
		byte[][] binned = new byte[features][rows];
		IntStream.range(0, features).parallel().forEach(f -> {
			cuts[f] = cutPoints(x, f, p.bins);
			for (int r = 0; r < rows; r++) {
				binned[f][r] = (byte) binOf(cuts[f], x[r][f]);
			}
		});
		double[] u = new double[rows];
		double[] g = new double[rows];
		double[] h = new double[rows];
		Random random = new Random(p.seed);
		List<Tree> trees = new ArrayList<>(p.trees);
		for (int t = 0; t < p.trees; t++) {
			gradients(u, queryStart, positive, g, h);
			boolean[] useQuery = new boolean[queries];
			int used = 0;
			for (int q = 0; q < queries; q++) {
				useQuery[q] = random.nextDouble() < p.queryFraction;
				if (useQuery[q]) {
					used += queryStart[q + 1] - queryStart[q];
				}
			}
			int[] sample = new int[used];
			int n = 0;
			for (int q = 0; q < queries; q++) {
				if (useQuery[q]) {
					for (int r = queryStart[q]; r < queryStart[q + 1]; r++) {
						sample[n++] = r;
					}
				}
			}
			boolean[] useFeature = new boolean[features];
			for (int f = 0; f < features; f++) {
				useFeature[f] = random.nextDouble() < p.featureFraction;
			}
			Tree tree = grow(sample, binned, cuts, g, h, useFeature, monotone, p);
			trees.add(tree);
			for (int r = 0; r < rows; r++) {
				u[r] += tree.predict(x[r]);
			}
		}
		return new Model(trees);
	}

	/** Per-row gradient and hessian of the listwise loss with a fixed "none" utility 0. */
	static void gradients(double[] u, int[] queryStart, int[] positive, double[] g, double[] h) {
		int queries = queryStart.length - 1;
		IntStream.range(0, queries).parallel().forEach(q -> {
			int from = queryStart[q];
			int to = queryStart[q + 1];
			double max = 0; // the "none" utility
			for (int r = from; r < to; r++) {
				max = Math.max(max, u[r]);
			}
			double z = Math.exp(-max);
			for (int r = from; r < to; r++) {
				z += Math.exp(u[r] - max);
			}
			for (int r = from; r < to; r++) {
				double pr = Math.exp(u[r] - max) / z;
				double y = positive[q] == r ? 1.0 : 0.0;
				g[r] = pr - y;
				h[r] = Math.max(pr * (1 - pr), 1e-6);
			}
		});
	}

	private static Tree grow(int[] sample, byte[][] binned, double[][] cuts, double[] g, double[] h,
			boolean[] useFeature, int[] monotone, Params p) {
		List<int[]> nodeRows = new ArrayList<>();
		List<double[]> nodeBounds = new ArrayList<>();
		List<Integer> nodeDepth = new ArrayList<>();
		IntList feature = new IntList();
		DoubleList threshold = new DoubleList();
		IntList left = new IntList();
		IntList right = new IntList();
		DoubleList value = new DoubleList();
		nodeRows.add(sample);
		nodeBounds.add(new double[] { Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY });
		nodeDepth.add(0);
		feature.add(-1);
		threshold.add(0);
		left.add(-1);
		right.add(-1);
		value.add(0);
		for (int node = 0; node < nodeRows.size(); node++) {
			int[] rows = nodeRows.get(node);
			double[] bounds = nodeBounds.get(node);
			double gs = 0;
			double hs = 0;
			for (int r : rows) {
				gs += g[r];
				hs += h[r];
			}
			double w = clip(-gs / (hs + p.lambda), bounds);
			value.set(node, w * p.learningRate);
			if (nodeDepth.get(node) >= p.maxDepth || hs < 2 * p.minChildHessian) {
				continue;
			}
			final double gTotal = gs;
			final double hTotal = hs;
			Split[] best = new Split[binned.length];
			IntStream.range(0, binned.length).parallel().forEach(f -> {
				if (useFeature[f] && cuts[f].length > 0) {
					best[f] = bestSplit(f, rows, binned[f], cuts[f].length + 1, g, h, gTotal, hTotal, bounds,
							monotone[f], p);
				}
			});
			Split split = null;
			for (Split s : best) {
				if (s != null && (split == null || s.gain > split.gain)) {
					split = s;
				}
			}
			if (split == null || split.gain <= p.minGain) {
				continue;
			}
			int nl = 0;
			for (int r : rows) {
				if ((binned[split.feature][r] & 0xFF) <= split.bin) {
					nl++;
				}
			}
			int[] lr = new int[nl];
			int[] rr = new int[rows.length - nl];
			int a = 0;
			int b = 0;
			for (int r : rows) {
				if ((binned[split.feature][r] & 0xFF) <= split.bin) {
					lr[a++] = r;
				} else {
					rr[b++] = r;
				}
			}
			double[] lb = bounds.clone();
			double[] rb = bounds.clone();
			if (monotone[split.feature] != 0) {
				double mid = (split.leftWeight + split.rightWeight) / 2;
				if (monotone[split.feature] > 0) {
					lb[1] = Math.min(lb[1], mid);
					rb[0] = Math.max(rb[0], mid);
				} else {
					lb[0] = Math.max(lb[0], mid);
					rb[1] = Math.min(rb[1], mid);
				}
			}
			feature.set(node, split.feature);
			threshold.set(node, cuts[split.feature][split.bin]);
			for (int side = 0; side < 2; side++) {
				int child = nodeRows.size();
				nodeRows.add(side == 0 ? lr : rr);
				nodeBounds.add(side == 0 ? lb : rb);
				nodeDepth.add(nodeDepth.get(node) + 1);
				feature.add(-1);
				threshold.add(0);
				left.add(-1);
				right.add(-1);
				value.add(0);
				if (side == 0) {
					left.set(node, child);
				} else {
					right.set(node, child);
				}
			}
			nodeRows.set(node, null);
		}
		return new Tree(feature.toArray(), threshold.toArray(), left.toArray(), right.toArray(), value.toArray());
	}

	private static final class Split {
		int feature;
		int bin;
		double gain;
		double leftWeight;
		double rightWeight;
	}

	private static Split bestSplit(int f, int[] rows, byte[] bins, int nBins, double[] g, double[] h, double gTotal,
			double hTotal, double[] bounds, int mono, Params p) {
		double[] gh = new double[nBins];
		double[] hh = new double[nBins];
		for (int r : rows) {
			int b = bins[r] & 0xFF;
			gh[b] += g[r];
			hh[b] += h[r];
		}
		double parent = gTotal * gTotal / (hTotal + p.lambda);
		Split best = null;
		double gl = 0;
		double hl = 0;
		for (int b = 0; b < nBins - 1; b++) {
			gl += gh[b];
			hl += hh[b];
			double gr = gTotal - gl;
			double hr = hTotal - hl;
			if (hl < p.minChildHessian || hr < p.minChildHessian) {
				continue;
			}
			double wl = clip(-gl / (hl + p.lambda), bounds);
			double wr = clip(-gr / (hr + p.lambda), bounds);
			if ((mono > 0 && wl > wr) || (mono < 0 && wl < wr)) {
				continue;
			}
			// Gain of the (possibly clipped) weights: -(2 G w + (H + lambda) w^2) summed over children.
			double gain = -(2 * gl * wl + (hl + p.lambda) * wl * wl) - (2 * gr * wr + (hr + p.lambda) * wr * wr)
					- parent;
			if (best == null || gain > best.gain) {
				best = new Split();
				best.feature = f;
				best.bin = b;
				best.gain = gain;
				best.leftWeight = wl;
				best.rightWeight = wr;
			}
		}
		return best;
	}

	private static double clip(double w, double[] bounds) {
		return Math.max(bounds[0], Math.min(bounds[1], w));
	}

	/** Up to {@code bins - 1} distinct quantile cut points of feature f; value <= cut[b] is bin <= b. */
	static double[] cutPoints(float[][] x, int f, int bins) {
		int rows = x.length;
		int sample = Math.min(rows, 200_000);
		float[] v = new float[sample];
		long step = Math.max(1, rows / sample);
		int n = 0;
		for (long r = 0; r < rows && n < sample; r += step) {
			v[n++] = x[(int) r][f];
		}
		v = Arrays.copyOf(v, n);
		Arrays.sort(v);
		double[] cuts = new double[bins - 1];
		int c = 0;
		for (int b = 1; b < bins; b++) {
			double q = v[Math.min(n - 1, (int) ((long) b * n / bins))];
			if (c == 0 || q > cuts[c - 1]) {
				cuts[c++] = q;
			}
		}
		// The last cut must leave something on its right, or it is not a split.
		if (c > 0 && cuts[c - 1] >= v[n - 1]) {
			c--;
		}
		return Arrays.copyOf(cuts, c);
	}

	static int binOf(double[] cuts, float v) {
		int lo = 0;
		int hi = cuts.length;
		while (lo < hi) {
			int mid = (lo + hi) >>> 1;
			if (v <= cuts[mid]) {
				hi = mid;
			} else {
				lo = mid + 1;
			}
		}
		return lo;
	}

	private static final class IntList {
		int[] a = new int[16];
		int n;

		void add(int v) {
			if (n == a.length) {
				a = Arrays.copyOf(a, n * 2);
			}
			a[n++] = v;
		}

		void set(int i, int v) {
			a[i] = v;
		}

		int[] toArray() {
			return Arrays.copyOf(a, n);
		}
	}

	private static final class DoubleList {
		double[] a = new double[16];
		int n;

		void add(double v) {
			if (n == a.length) {
				a = Arrays.copyOf(a, n * 2);
			}
			a[n++] = v;
		}

		void set(int i, double v) {
			a[i] = v;
		}

		double[] toArray() {
			return Arrays.copyOf(a, n);
		}
	}
}
