package uk.ac.quadram.metapilot.core.opensearch;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * Reduces the classifier's per-candidate scores to one PSM per spectrum, which
 * is the unit target-decoy FDR is defined over.
 *
 * <h2>Why this is required</h2>
 * With {@code output_report_topN > 1} a single spectrum contributes several
 * competing candidates. Computing FDR across all of them is wrong in three
 * ways: decoy counts stop being counts of independent false discoveries (one
 * spectrum can contribute several); the accepted-target count reports
 * <em>candidates</em> rather than identified spectra, inflating the headline
 * number; and nothing prevents two candidates from the same spectrum both
 * being accepted, which is meaningless as an identification.
 *
 * <h2>Why this is not a "hard filter"</h2>
 * The governing principle elsewhere in this pipeline is that the classifier,
 * not the engine, decides between candidates — so every candidate must reach
 * the classifier, and every candidate does. This collapse happens strictly
 * <em>after</em> scoring: the classifier has already been free to promote a
 * rank-3 candidate over the engine's rank-1 pick. What is being chosen here is
 * the unit of analysis for FDR, not which evidence the model may see. It is
 * the same step Percolator and MS²Rescore perform before reporting.
 */
public final class SpectrumLevelCollapse {

	private static final Logger LOGGER = LogManager.getLogger(SpectrumLevelCollapse.class);

	private SpectrumLevelCollapse() {
	}

	public static final class Result {
		/** One PSM per spectrum: the highest-scoring candidate. */
		public final List<OpenSearchPsm> best;
		public final int spectraCount;
		public final int candidatesDropped;
		/**
		 * How often the winner came from each engine rank (index 0 = rank 1).
		 * This is the empirical answer to "is topN too high?": if the tail is
		 * always zero, the extra candidates cost time and buy nothing.
		 */
		public final int[] winnerRankHistogram;

		Result(List<OpenSearchPsm> best, int spectraCount, int candidatesDropped, int[] winnerRankHistogram) {
			this.best = best;
			this.spectraCount = spectraCount;
			this.candidatesDropped = candidatesDropped;
			this.winnerRankHistogram = winnerRankHistogram;
		}

		/** Fraction of spectra where the classifier overruled the engine's rank-1 pick. */
		public double overruleFraction() {
			if (spectraCount == 0) {
				return 0;
			}
			int overruled = 0;
			for (int rank = 1; rank < winnerRankHistogram.length; rank++) {
				overruled += winnerRankHistogram[rank];
			}
			return (double) overruled / spectraCount;
		}

		public String describeRanks() {
			StringBuilder sb = new StringBuilder();
			for (int i = 0; i < winnerRankHistogram.length; i++) {
				if (winnerRankHistogram[i] > 0) {
					if (sb.length() > 0) {
						sb.append(", ");
					}
					sb.append("rank ").append(i + 1).append(": ").append(winnerRankHistogram[i]);
				}
			}
			return sb.toString();
		}
	}

	/**
	 * Keeps the best-scoring candidate per {@code fileName:scan}.
	 * <p>
	 * Ranking uses {@link OpenSearchPsm#getComparableScore()} (the per-engine
	 * calibrated score, else the classifier logit, else its probability) and then
	 * hyperscore, so this is still meaningful if called before scoring.
	 */
	public static Result bestPerSpectrum(List<OpenSearchPsm> psms) {
		Map<String, OpenSearchPsm> bestBySpectrum = new LinkedHashMap<>();

		for (OpenSearchPsm psm : psms) {
			String key = psm.getFileName() + ":" + psm.getScan();
			OpenSearchPsm incumbent = bestBySpectrum.get(key);
			if (incumbent == null || rankingScore(psm) > rankingScore(incumbent)) {
				bestBySpectrum.put(key, psm);
			}
		}

		List<OpenSearchPsm> best = new ArrayList<>(bestBySpectrum.values());

		int maxRank = 1;
		for (OpenSearchPsm psm : best) {
			maxRank = Math.max(maxRank, Math.max(1, psm.getHitRank()));
		}
		int[] histogram = new int[maxRank];
		for (OpenSearchPsm psm : best) {
			int rank = Math.max(1, psm.getHitRank());
			histogram[rank - 1]++;
		}

		Result result = new Result(best, best.size(), psms.size() - best.size(), histogram);
		LOGGER.info(String.format(
				"Spectrum-level collapse: %d candidates -> %d spectra (%d dropped). Winner engine rank: %s. "
						+ "The classifier overruled the engine's top hit for %.1f%% of spectra.",
				psms.size(), result.spectraCount, result.candidatesDropped, result.describeRanks(),
				result.overruleFraction() * 100));
		return result;
	}

	private static double rankingScore(OpenSearchPsm psm) {
		double score = psm.getComparableScore();
		if (!Double.isNaN(score)) {
			return score;
		}
		return Double.isNaN(psm.getHyperscore()) ? Double.NEGATIVE_INFINITY : psm.getHyperscore();
	}
}
