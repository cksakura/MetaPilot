package uk.ac.quadram.metapilot.core.opensearch;

import java.util.Arrays;
import java.util.Comparator;

/**
 * Target-decoy FDR and q-value estimation for a 1:1 target/decoy database.
 * <p>
 * Shared by the classifier's semi-supervised training loop (to pick confident
 * positives each round) and by per-peak validation (to cut each modification's
 * PSM list), so both use the same estimator.
 * <p>
 * <b>Assumes a 1:1 target:decoy ratio</b>, which {@link
 * uk.ac.quadram.metapilot.core.prodb.TargetDecoyFasta} enforces. The estimator
 * {@code FDR = decoys / targets} holds only at 1:1; against the FragPipe
 * combined stage's 1:3 database it would overstate FDR roughly threefold.
 * <p>
 * <b>q-values are monotonised.</b> The running ratio {@code d/t} fluctuates
 * rank to rank, so used directly — as the original implementation did, stopping
 * at the first rank to breach the threshold — a single unlucky decoy near the
 * top truncates a list that should have continued. The q-value at rank <i>i</i>
 * is instead the minimum running FDR over all ranks at or below <i>i</i>, which
 * is non-increasing with score, and the cut is taken at the <em>last</em> rank
 * satisfying the threshold.
 */
public final class TargetDecoyFdr {

	private TargetDecoyFdr() {
	}

	/**
	 * Monotonised q-values, returned in the caller's original index order.
	 * Higher score means better; ties are broken deterministically by index so
	 * repeated runs agree.
	 *
	 * @param scores   discriminant score per PSM (higher = more confident)
	 * @param isTarget target/decoy flag per PSM
	 */
	public static double[] qValues(double[] scores, boolean[] isTarget) {
		if (scores.length != isTarget.length) {
			throw new IllegalArgumentException("scores and isTarget differ in length");
		}
		int n = scores.length;
		double[] q = new double[n];
		if (n == 0) {
			return q;
		}

		Integer[] order = sortedByScoreDesc(scores);

		// Running FDR walking down the ranking.
		double[] running = new double[n];
		int targets = 0;
		int decoys = 0;
		for (int rank = 0; rank < n; rank++) {
			if (isTarget[order[rank]]) {
				targets++;
			} else {
				decoys++;
			}
			running[rank] = targets == 0 ? 1.0 : Math.min(1.0, (double) decoys / (double) targets);
		}

		// Monotonise from the bottom so q is non-increasing with score.
		double best = 1.0;
		for (int rank = n - 1; rank >= 0; rank--) {
			best = Math.min(best, running[rank]);
			q[order[rank]] = best;
		}
		return q;
	}

	/**
	 * Monotonised q-values with the {@code (decoys + 1) / targets} estimate and tied scores walked
	 * as one block, in the caller's original index order.
	 * <p>
	 * For small populations, such as the joint identification's per-group walks: with 0 decoys
	 * among 40 targets, {@link #qValues} reports q = 0 where this reports 1/40. The +1 makes the
	 * estimate conservative (He et al. 2015; Levitsky et al. 2017) and costs nothing on a large
	 * group. A block of tied scores is accepted or rejected together, so a target and a decoy at
	 * the same score never get different q-values from their row order.
	 */
	public static double[] qValuesConservative(double[] scores, boolean[] isTarget) {
		return qValuesConservative(scores, isTarget, 1.0);
	}

	/**
	 * As {@link #qValuesConservative(double[], boolean[])}, each decoy standing for
	 * {@code decoyWeight} decoys: {@code decoyWeight * (decoys + 1) / targets}. RESET passes only the
	 * estimating half of the decoys and weight 2.
	 */
	public static double[] qValuesConservative(double[] scores, boolean[] isTarget, double decoyWeight) {
		if (scores.length != isTarget.length) {
			throw new IllegalArgumentException("scores and isTarget differ in length");
		}
		int n = scores.length;
		double[] q = new double[n];
		if (n == 0) {
			return q;
		}
		Integer[] order = sortedByScoreDesc(scores);
		double[] running = new double[n];
		int targets = 0;
		int decoys = 0;
		int rank = 0;
		while (rank < n) {
			int end = rank;
			while (end + 1 < n && sameScore(scores[order[end + 1]], scores[order[rank]])) {
				end++;
			}
			for (int k = rank; k <= end; k++) {
				if (isTarget[order[k]]) {
					targets++;
				} else {
					decoys++;
				}
			}
			double fdr = targets == 0 ? 1.0 : Math.min(1.0, decoyWeight * (decoys + 1.0) / targets);
			for (int k = rank; k <= end; k++) {
				running[k] = fdr;
			}
			rank = end + 1;
		}
		double best = 1.0;
		for (int k = n - 1; k >= 0; k--) {
			best = Math.min(best, running[k]);
			q[order[k]] = best;
		}
		return q;
	}

	private static boolean sameScore(double a, double b) {
		return a == b || (Double.isNaN(a) && Double.isNaN(b));
	}

	/**
	 * The score at the last rank whose q-value still satisfies
	 * {@code qThreshold}. PSMs scoring at or above this are accepted.
	 *
	 * @return {@code Double.POSITIVE_INFINITY} when no rank qualifies, so that
	 *         comparing {@code score >= threshold} accepts nothing
	 */
	public static double scoreThresholdAt(double[] scores, boolean[] isTarget, double qThreshold) {
		double[] q = qValues(scores, isTarget);
		Integer[] order = sortedByScoreDesc(scores);
		double threshold = Double.POSITIVE_INFINITY;
		for (int rank = 0; rank < order.length; rank++) {
			if (q[order[rank]] <= qThreshold) {
				threshold = scores[order[rank]];
			}
		}
		return threshold;
	}

	/**
	 * Indices accepted at {@code qThreshold}, restricted to targets. This is the
	 * confident-positive set the semi-supervised training loop re-seeds from
	 * each round.
	 */
	public static boolean[] acceptedTargets(double[] scores, boolean[] isTarget, double qThreshold) {
		double[] q = qValues(scores, isTarget);
		boolean[] accepted = new boolean[scores.length];
		for (int i = 0; i < scores.length; i++) {
			accepted[i] = isTarget[i] && q[i] <= qThreshold;
		}
		return accepted;
	}

	/** Number of targets accepted at {@code qThreshold} — the usual yield metric. */
	public static int countAcceptedTargets(double[] scores, boolean[] isTarget, double qThreshold) {
		boolean[] accepted = acceptedTargets(scores, isTarget, qThreshold);
		int n = 0;
		for (boolean a : accepted) {
			if (a) {
				n++;
			}
		}
		return n;
	}

	private static Integer[] sortedByScoreDesc(double[] scores) {
		Integer[] order = new Integer[scores.length];
		for (int i = 0; i < scores.length; i++) {
			order[i] = i;
		}
		Arrays.sort(order, new Comparator<Integer>() {
			@Override
			public int compare(Integer a, Integer b) {
				double sa = scores[a];
				double sb = scores[b];
				// NaN sorts last so unscored PSMs never occupy confident ranks.
				boolean na = Double.isNaN(sa);
				boolean nb = Double.isNaN(sb);
				if (na && nb) {
					return Integer.compare(a, b);
				}
				if (na) {
					return 1;
				}
				if (nb) {
					return -1;
				}
				int c = Double.compare(sb, sa);
				return c != 0 ? c : Integer.compare(a, b);
			}
		});
		return order;
	}
}
