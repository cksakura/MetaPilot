package uk.ac.quadram.metapilot.core.unified;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

/**
 * Group-walk (Freestone, Short, Noble and Keich, Bioinformatics 2022): FDR control by
 * target-decoy competition with a group structure, rigorous where separate walks per group are
 * not.
 * <p>
 * Each group is sorted by increasing score; a frontier per group marks the active set (everything
 * at or above it). After every advance the FDR of the active set is estimated as
 * {@code (decoyWeight * decoys + correction) / targets}; at level alpha the walk stops at the first
 * estimate {@code <= alpha} and reports the active targets. The frontier advances round-robin until
 * every group has passed {@code K} items, then in the group with the most decoys among its last
 * {@code K} passed items. The path does not depend on alpha, so a target's q-value is the smallest
 * estimate reached while it was still active.
 * <p>
 * {@code decoyWeight = 2} is the RESET estimate (only half the decoys, the estimating ones, take
 * part); 1 is plain target-decoy competition.
 */
public final class GroupWalk {

	private GroupWalk() {
	}

	/**
	 * @param scores  winning score per hypothesis (higher is better)
	 * @param target  true for a target win, false for a (estimating) decoy win
	 * @param group   group id per hypothesis, 0..groups-1
	 * @return q-value per hypothesis (targets meaningful; decoys get the value at their removal)
	 */
	public static double[] qValues(double[] scores, boolean[] target, int[] group, int groups, int k,
			double decoyWeight, double correction) {
		int n = scores.length;
		double[] q = new double[n];
		Arrays.fill(q, 1.0);
		if (n == 0) {
			return q;
		}
		List<List<Integer>> members = new ArrayList<>();
		for (int g = 0; g < groups; g++) {
			members.add(new ArrayList<>());
		}
		for (int i = 0; i < n; i++) {
			members.get(group[i]).add(i);
		}
		int[][] sorted = new int[groups][];
		for (int g = 0; g < groups; g++) {
			List<Integer> m = members.get(g);
			// Ascending score; ties in a pseudo-random but reproducible order. Ordering ties by index let
			// the caller's row order decide them - with rows sorted decoys first, every tied block passed
			// its decoys before its targets, and the decoy-count rule that moves the front then read a
			// clean window and left decoys behind: 576 peptides at 1% where 73,843 were reachable
			// (12mix, boosted-tree peptide scores, 53,776 distinct values over 314,771 units).
			m.sort(Comparator.<Integer>comparingDouble(i -> rank(scores[i])).thenComparingLong(GroupWalk::tieKey));
			sorted[g] = m.stream().mapToInt(Integer::intValue).toArray();
		}
		int targets = 0;
		int decoys = 0;
		for (int i = 0; i < n; i++) {
			if (target[i]) {
				targets++;
			} else {
				decoys++;
			}
		}
		int[] front = new int[groups];
		double runningMin = Double.POSITIVE_INFINITY;
		int remaining = n;
		int roundRobin = 0;
		while (remaining > 0) {
			double fdr = targets == 0 ? 1.0 : Math.min(1.0, (decoyWeight * decoys + correction) / targets);
			runningMin = Math.min(runningMin, fdr);
			int g = chooseGroup(sorted, front, target, k, roundRobin);
			if (g < 0) {
				break;
			}
			if (front[g] < k) {
				roundRobin = g + 1;
			}
			int item = sorted[g][front[g]++];
			q[item] = Math.min(1.0, runningMin);
			if (target[item]) {
				targets--;
			} else {
				decoys--;
			}
			remaining--;
		}
		return q;
	}

	private static int chooseGroup(int[][] sorted, int[] front, boolean[] target, int k, int roundRobin) {
		int groups = sorted.length;
		// Round-robin while any unfinished group has passed fewer than K items.
		for (int step = 0; step < groups; step++) {
			int g = (roundRobin + step) % groups;
			if (front[g] < sorted[g].length && front[g] < k) {
				return g;
			}
		}
		int best = -1;
		int bestDecoys = -1;
		for (int g = 0; g < groups; g++) {
			if (front[g] >= sorted[g].length) {
				continue;
			}
			int d = 0;
			for (int j = Math.max(0, front[g] - k); j < front[g]; j++) {
				if (!target[sorted[g][j]]) {
					d++;
				}
			}
			if (d > bestDecoys) {
				bestDecoys = d;
				best = g;
			}
		}
		return best;
	}

	private static double rank(double score) {
		return Double.isNaN(score) ? Double.NEGATIVE_INFINITY : score;
	}

	/** SplitMix64 of the index: a fixed pseudo-random order, independent of labels and row order. */
	static long tieKey(int index) {
		long z = index + 0x9E3779B97F4A7C15L;
		z = (z ^ (z >>> 30)) * 0xBF58476D1CE4E5B9L;
		z = (z ^ (z >>> 27)) * 0x94D049BB133111EBL;
		return z ^ (z >>> 31);
	}

	/**
	 * CONGA's group construction (Freestone, Noble and Keich, JPR 2024), from scores alone: a fixed
	 * first batch, then each other batch, largest first, joins the existing group whose score
	 * distribution a two-sample Kolmogorov-Smirnov test cannot tell from its own (p > {@code beta},
	 * the most similar one) or starts a new group; finally every group below {@code minSize} joins the
	 * group it is most similar to. Group-walk needs this: many small groups of noisy scores stall its
	 * front (on 12mix, 50 small peak batches left no peptide at q <= 1% although a pooled walk reached
	 * 74,050).
	 *
	 * @param batchScores scores per batch, targets and decoys pooled (no label is read)
	 * @param fixed       index of the batch that is always its own group (the exact one), or -1
	 * @return group index per batch
	 */
	public static int[] mergeBatches(List<double[]> batchScores, int fixed, double beta, int minSize) {
		int b = batchScores.size();
		int[] group = new int[b];
		Arrays.fill(group, -1);
		List<List<Integer>> groups = new ArrayList<>();
		List<double[]> pooled = new ArrayList<>();
		if (fixed >= 0) {
			group[fixed] = 0;
			groups.add(new ArrayList<>(List.of(fixed)));
			pooled.add(sortedCopy(batchScores.get(fixed)));
		}
		Integer[] order = new Integer[b];
		for (int i = 0; i < b; i++) {
			order[i] = i;
		}
		Arrays.sort(order, (x, y) -> Integer.compare(batchScores.get(y).length, batchScores.get(x).length));
		int firstOpen = fixed >= 0 ? 1 : 0;
		for (int i : order) {
			if (i == fixed) {
				continue;
			}
			double[] s = sortedCopy(batchScores.get(i));
			int best = -1;
			double bestP = beta;
			for (int g = firstOpen; g < groups.size(); g++) {
				double p = ksPValue(s, pooled.get(g));
				if (p > bestP) {
					bestP = p;
					best = g;
				}
			}
			if (best < 0) {
				group[i] = groups.size();
				groups.add(new ArrayList<>(List.of(i)));
				pooled.add(s);
			} else {
				group[i] = best;
				groups.get(best).add(i);
				pooled.set(best, merge(pooled.get(best), s));
			}
		}
		// Groups below the minimum size join their most similar open group.
		boolean changed = true;
		while (changed) {
			changed = false;
			int smallest = -1;
			for (int g = firstOpen; g < groups.size(); g++) {
				if (!groups.get(g).isEmpty() && pooled.get(g).length < minSize
						&& (smallest < 0 || pooled.get(g).length < pooled.get(smallest).length)) {
					smallest = g;
				}
			}
			if (smallest < 0) {
				break;
			}
			int target = -1;
			double bestP = -1;
			for (int g = firstOpen; g < groups.size(); g++) {
				if (g != smallest && !groups.get(g).isEmpty()) {
					double p = ksPValue(pooled.get(smallest), pooled.get(g));
					if (p > bestP) {
						bestP = p;
						target = g;
					}
				}
			}
			if (target < 0) {
				break; // the only open group; keep it
			}
			for (int i : groups.get(smallest)) {
				group[i] = target;
			}
			groups.get(target).addAll(groups.get(smallest));
			groups.get(smallest).clear();
			pooled.set(target, merge(pooled.get(target), pooled.get(smallest)));
			pooled.set(smallest, new double[0]);
			changed = true;
		}
		// Renumber densely.
		int[] dense = new int[groups.size()];
		Arrays.fill(dense, -1);
		int next = 0;
		for (int i = 0; i < b; i++) {
			if (dense[group[i]] < 0) {
				dense[group[i]] = next++;
			}
		}
		if (fixed >= 0 && dense[group[fixed]] != 0) {
			// Keep the fixed batch as group 0.
			int f = dense[group[fixed]];
			for (int g = 0; g < dense.length; g++) {
				if (dense[g] == 0) {
					dense[g] = f;
				} else if (dense[g] == f) {
					dense[g] = 0;
				}
			}
		}
		for (int i = 0; i < b; i++) {
			group[i] = dense[group[i]];
		}
		return group;
	}

	/** Asymptotic two-sample Kolmogorov-Smirnov p-value of two sorted samples. */
	static double ksPValue(double[] a, double[] b) {
		if (a.length == 0 || b.length == 0) {
			return 1.0;
		}
		int i = 0;
		int j = 0;
		double d = 0;
		while (i < a.length && j < b.length) {
			double v = Math.min(a[i], b[j]);
			while (i < a.length && a[i] <= v) {
				i++;
			}
			while (j < b.length && b[j] <= v) {
				j++;
			}
			d = Math.max(d, Math.abs(i / (double) a.length - j / (double) b.length));
		}
		double ne = a.length * (double) b.length / (a.length + b.length);
		double lambda = (Math.sqrt(ne) + 0.12 + 0.11 / Math.sqrt(ne)) * d;
		if (lambda < 0.27) {
			return 1.0; // the series is 1 to double precision here and converges too slowly to sum
		}
		double sum = 0;
		for (int k = 1; k <= 100; k++) {
			double term = 2 * (k % 2 == 1 ? 1 : -1) * Math.exp(-2.0 * k * k * lambda * lambda);
			sum += term;
			if (Math.abs(term) < 1e-12) {
				break;
			}
		}
		return Math.max(0, Math.min(1, sum));
	}

	private static double[] sortedCopy(double[] v) {
		double[] c = v.clone();
		Arrays.sort(c);
		return c;
	}

	private static double[] merge(double[] a, double[] b) {
		double[] out = new double[a.length + b.length];
		int i = 0;
		int j = 0;
		int k = 0;
		while (i < a.length || j < b.length) {
			out[k++] = j >= b.length || (i < a.length && a[i] <= b[j]) ? a[i++] : b[j++];
		}
		return out;
	}
}
