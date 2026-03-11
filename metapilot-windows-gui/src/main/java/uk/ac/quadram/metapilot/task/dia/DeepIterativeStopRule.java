package uk.ac.quadram.metapilot.task.dia;

import java.util.Arrays;
import java.util.List;

/**
 * Shared stop-rule computation for deep iterative DIA workflows.
 * Supports legacy algorithm variants to keep behavior unchanged per workflow.
 */
public final class DeepIterativeStopRule {

	public enum AlgorithmProfile {
		DEEP_COMPREHENSIVE,
		DEEP_MARKER
	}

	public static final class StopRuleResult {
		public final int stopRound;
		public final int nRounds;
		public final int kChord;
		public final int kPiece;
		public final int kKnee;
		public final int k0;
		public final int tMin;

		public StopRuleResult(int stopRound, int nRounds, int kChord, int kPiece, int kKnee, int k0, int tMin) {
			this.stopRound = stopRound;
			this.nRounds = nRounds;
			this.kChord = kChord;
			this.kPiece = kPiece;
			this.kKnee = kKnee;
			this.k0 = k0;
			this.tMin = tMin;
		}
	}

	private DeepIterativeStopRule() {
	}

	public static StopRuleResult computeStopRoundFromIncrements(List<Integer> inc, double relativeGainThreshold,
			int consecutiveRoundsForPlateau, AlgorithmProfile profile) {
		if (inc == null) {
			return new StopRuleResult(0, 0, -1, -1, -1, -1, -1);
		}
		int n = inc.size();
		if (n == 0) {
			return new StopRuleResult(0, 0, -1, -1, -1, -1, -1);
		}
		if (n < 8) {
			return new StopRuleResult(0, n, -1, -1, -1, -1, -1);
		}

		double[] cum = new double[n];
		cum[0] = inc.get(0);
		for (int i = 1; i < n; i++) {
			cum[i] = cum[i - 1] + inc.get(i);
		}

		int kChord = (profile == AlgorithmProfile.DEEP_COMPREHENSIVE) ? computeKneeIndexComprehensive(cum)
				: computeKneeIndexMarker(cum);
		int kPiece = (profile == AlgorithmProfile.DEEP_COMPREHENSIVE) ? computePiecewiseComprehensive(cum)
				: computePiecewiseMarker(cum);
		int kKnee = (profile == AlgorithmProfile.DEEP_COMPREHENSIVE) ? computeKneedleComprehensive(cum)
				: computeKneedleMarker(cum);
		int k0 = (profile == AlgorithmProfile.DEEP_COMPREHENSIVE) ? medianOfThreeComprehensive(kChord, kPiece, kKnee)
				: medianOfThreeMarker(kChord, kPiece, kKnee);
		if (k0 < 1) {
			k0 = 1;
		}
		if (k0 > n) {
			k0 = n;
		}
		int tMin = Math.max(8, k0);

		for (int t = tMin; t <= n; t++) {
			boolean allBelow = true;
			for (int k = 0; k < consecutiveRoundsForPlateau && allBelow; k++) {
				int roundIndex = t - 1 - k;
				if (roundIndex < 0) {
					allBelow = false;
					break;
				}
				double cumulative = cum[roundIndex];
				if (cumulative <= 0) {
					allBelow = false;
					break;
				}
				int incVal = inc.get(roundIndex) != null ? inc.get(roundIndex).intValue() : 0;
				double relativeGain = incVal / cumulative;
				if (relativeGain >= relativeGainThreshold) {
					allBelow = false;
				}
			}
			if (allBelow) {
				return new StopRuleResult(t, n, kChord, kPiece, kKnee, k0, tMin);
			}
		}
		return new StopRuleResult(0, n, kChord, kPiece, kKnee, k0, tMin);
	}

	// ---- Deep Comprehensive legacy helpers ----

	private static int computeKneeIndexComprehensive(double[] y) {
		int n = y.length;
		if (n <= 2) {
			return n;
		}
		double x0 = 0.0;
		double y0 = y[0];
		double x1 = 1.0;
		double y1 = y[n - 1];
		double dx = x1 - x0;
		double dy = y1 - y0;
		double denom = Math.sqrt(dx * dx + dy * dy);
		if (denom == 0.0) {
			return 1;
		}
		double maxDist = -1.0;
		int idx = 1;
		for (int i = 1; i < n - 1; i++) {
			double x = (double) i / (double) (n - 1);
			double yy = y[i];
			double num = Math.abs(dy * x - dx * yy + x1 * y0 - y1 * x0);
			double dist = num / denom;
			if (dist > maxDist) {
				maxDist = dist;
				idx = i;
			}
		}
		return idx + 1;
	}

	private static int computePiecewiseComprehensive(double[] y) {
		int n = y.length;
		if (n <= 2) {
			return n;
		}
		double[] x = new double[n];
		for (int i = 0; i < n; i++) {
			x[i] = i;
		}
		double bestSse = Double.POSITIVE_INFINITY;
		int bestK = n;
		for (int k = 1; k < n - 1; k++) {
			double sse = 0.0;
			sse += segmentSSEComprehensive(x, y, 0, k);
			sse += segmentSSEComprehensive(x, y, k, n - 1);
			if (sse < bestSse) {
				bestSse = sse;
				bestK = k;
			}
		}
		return bestK + 1;
	}

	private static double segmentSSEComprehensive(double[] x, double[] y, int i0, int i1) {
		if (i1 <= i0) {
			return 0.0;
		}
		double sumX = 0.0;
		double sumY = 0.0;
		double sumXY = 0.0;
		double sumXX = 0.0;
		int n = i1 - i0 + 1;
		for (int i = i0; i <= i1; i++) {
			double xi = x[i];
			double yi = y[i];
			sumX += xi;
			sumY += yi;
			sumXY += xi * yi;
			sumXX += xi * xi;
		}
		double denom = (n * sumXX - sumX * sumX);
		double a = 0.0;
		double b = 0.0;
		if (denom != 0.0) {
			a = (n * sumXY - sumX * sumY) / denom;
			b = (sumY - a * sumX) / n;
		} else {
			b = sumY / n;
		}
		double sse = 0.0;
		for (int i = i0; i <= i1; i++) {
			double yi = y[i];
			double yhat = a * x[i] + b;
			double diff = yi - yhat;
			sse += diff * diff;
		}
		return sse;
	}

	private static int computeKneedleComprehensive(double[] y) {
		int n = y.length;
		if (n <= 2) {
			return n;
		}
		double minY = Double.POSITIVE_INFINITY;
		double maxY = Double.NEGATIVE_INFINITY;
		for (double v : y) {
			if (v < minY) {
				minY = v;
			}
			if (v > maxY) {
				maxY = v;
			}
		}
		double range = maxY - minY;
		if (range == 0.0) {
			return 1;
		}
		double[] yn = new double[n];
		for (int i = 0; i < n; i++) {
			yn[i] = (y[i] - minY) / range;
		}
		double[] diff = new double[n];
		for (int i = 0; i < n; i++) {
			double x = (double) i / (double) (n - 1);
			diff[i] = yn[i] - x;
		}
		double maxDiff = Double.NEGATIVE_INFINITY;
		int idx = 1;
		for (int i = 1; i < n - 1; i++) {
			if (diff[i] > maxDiff) {
				maxDiff = diff[i];
				idx = i;
			}
		}
		return idx + 1;
	}

	private static int medianOfThreeComprehensive(int a, int b, int c) {
		if ((a >= b && a <= c) || (a <= b && a >= c)) {
			return a;
		}
		if ((b >= a && b <= c) || (b <= a && b >= c)) {
			return b;
		}
		return c;
	}

	// ---- Deep Marker legacy helpers ----

	private static int computeKneeIndexMarker(double[] y) {
		if (y == null || y.length == 0) {
			return 0;
		}
		if (y.length <= 2) {
			return y.length;
		}

		int n = y.length;
		double x1 = 0.0;
		double y1 = y[0];
		double x2 = 1.0;
		double y2 = y[n - 1];
		double denom = Math.sqrt((y2 - y1) * (y2 - y1) + (x2 - x1) * (x2 - x1));
		if (denom == 0.0) {
			return n;
		}

		double maxDist = -1.0;
		int bestIndex = 1;
		for (int i = 0; i < n; i++) {
			double x = (n == 1) ? 0.0 : ((double) i / (double) (n - 1));
			double yi = y[i];
			double dist = Math.abs((y2 - y1) * x - (x2 - x1) * yi + x2 * y1 - y2 * x1) / denom;
			if (dist > maxDist) {
				maxDist = dist;
				bestIndex = i + 1;
			}
		}
		return bestIndex;
	}

	private static int computeKneedleMarker(double[] y) {
		if (y == null || y.length == 0) {
			return 0;
		}
		if (y.length == 1) {
			return 1;
		}
		double min = Arrays.stream(y).min().orElse(0.0);
		double max = Arrays.stream(y).max().orElse(0.0);
		double range = max - min;
		double bestDiff = -Double.MAX_VALUE;
		int bestIdx = 1;
		for (int i = 0; i < y.length; i++) {
			double x = (y.length == 1) ? 0.0 : ((double) i / (double) (y.length - 1));
			double yn = (range == 0.0) ? 0.0 : ((y[i] - min) / range);
			double diff = yn - (1.0 - x);
			if (diff > bestDiff) {
				bestDiff = diff;
				bestIdx = i + 1;
			}
		}
		return bestIdx;
	}

	private static int computePiecewiseMarker(double[] y) {
		if (y == null || y.length == 0) {
			return 0;
		}
		if (y.length <= 3) {
			return y.length;
		}
		double bestSse = Double.MAX_VALUE;
		int bestK = 2;
		for (int k = 2; k <= y.length - 1; k++) {
			double sse1 = linearFitSseMarker(y, 0, k - 1);
			double sse2 = linearFitSseMarker(y, k - 1, y.length - 1);
			double sse = sse1 + sse2;
			if (sse < bestSse) {
				bestSse = sse;
				bestK = k;
			}
		}
		return bestK;
	}

	private static double linearFitSseMarker(double[] y, int start, int end) {
		int n = end - start + 1;
		if (n <= 1) {
			return 0.0;
		}
		double sumX = 0.0;
		double sumY = 0.0;
		double sumXY = 0.0;
		double sumXX = 0.0;
		for (int i = start; i <= end; i++) {
			double x = i - start;
			sumX += x;
			sumY += y[i];
			sumXY += x * y[i];
			sumXX += x * x;
		}
		double denom = n * sumXX - sumX * sumX;
		double b = denom == 0.0 ? 0.0 : (n * sumXY - sumX * sumY) / denom;
		double a = (sumY - b * sumX) / n;
		double sse = 0.0;
		for (int i = start; i <= end; i++) {
			double x = i - start;
			double fit = a + b * x;
			double err = y[i] - fit;
			sse += err * err;
		}
		return sse;
	}

	private static int medianOfThreeMarker(int a, int b, int c) {
		int[] v = new int[] { a, b, c };
		Arrays.sort(v);
		return v[1];
	}
}
