package uk.ac.quadram.metapilot.core.runner;

import uk.ac.quadram.metapilot.core.tools.FormatTool;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;

/**
 * Writes final MAG genome ranking artifacts for taxonomy/report labels.
 *
 * This ranking is intentionally separate from marker-based genome selection for
 * refined searches. It assumes confidence filtering has already been applied by
 * post-analysis and ranks the remaining genomes by final all-protein abundance.
 */
public final class MagFinalGenomeRankingService {

	public static final double DEFAULT_PRIVILEGE_UNIQUE_PEPTIDE_FRACTION = 0.01;

	private MagFinalGenomeRankingService() {
	}

	public static final class Request {
		public File taxonomyAnalysisDir;
		public Map<String, double[]> genomeIntensityMap;
		public Map<String, Integer> genomePsmCountMap;
		public Map<String, ? extends Set<String>> totalPepCountMap;
		public Map<String, ? extends Set<String>> razorPepCountMap;
		public Map<String, ? extends Set<String>> uniquePepCountMap;
		public double privilegeUniquePeptideFraction = DEFAULT_PRIVILEGE_UNIQUE_PEPTIDE_FRACTION;
	}

	public static final class Result {
		public File evidenceTsv;
		public File rankedTsv;
		public File selectedTsv;
		public int genomeCount;
		public int privilegeCount;
		public int greedyCount;
		public int focusKnee;
		public int comprehensiveKnee;
	}

	private static final class GenomeRow {
		String genome;
		double intensity;
		double relativeAbundance;
		int totalPeptides;
		int razorPeptides;
		int uniquePeptides;
		int psmCount;
		HashSet<String> coveragePeptides = new HashSet<>();
		HashSet<String> uniquePeptideSet = new HashSet<>();
	}

	public static Result write(Request request) throws IOException {
		validate(request);
		request.taxonomyAnalysisDir.mkdirs();

		List<GenomeRow> rows = buildRows(request);
		double totalIntensity = totalIntensity(rows);
		if (totalIntensity <= 0.0 && !rows.isEmpty()) {
			int genomesWithPeptides = 0;
			for (GenomeRow row : rows) {
				if (!row.coveragePeptides.isEmpty()) {
					genomesWithPeptides++;
				}
			}
			throw new IOException("Final genome intensity total is zero across " + rows.size() + " genomes ("
					+ genomesWithPeptides + " with assigned peptides); RelativeAbundance cannot be calculated. "
					+ "The genomes were identified but every peptide intensity reaching post-analysis was zero: "
					+ "check that the quantified peptide table read by post-analysis is the right file and that its "
					+ "intensity columns matched (DDA: FlashLFQ/isobaric peptide reader column match in the task log; "
					+ "DIA: parquet/PR-matrix intensity transfer).");
		}
		rows.sort(abundanceComparator());

		File evidenceTsv = new File(request.taxonomyAnalysisDir, "genome_final_evidence.tsv");
		File rankedTsv = new File(request.taxonomyAnalysisDir, "genome_final_ranked.tsv");
		File selectedTsv = new File(request.taxonomyAnalysisDir, "genome_final_selected.tsv");

		writeEvidence(evidenceTsv, rows);
		Ranking ranking = rank(rows, request.privilegeUniquePeptideFraction);
		writeRanked(rankedTsv, rows, ranking);
		writeSelected(selectedTsv, ranking);

		Result result = new Result();
		result.evidenceTsv = evidenceTsv;
		result.rankedTsv = rankedTsv;
		result.selectedTsv = selectedTsv;
		result.genomeCount = rows.size();
		result.privilegeCount = ranking.privileged.size();
		result.greedyCount = ranking.greedyOrder.size();
		result.focusKnee = ranking.focusKnee;
		result.comprehensiveKnee = ranking.comprehensiveKnee;
		return result;
	}

	private static void validate(Request request) throws IOException {
		if (request == null) {
			throw new IllegalArgumentException("request is null");
		}
		if (request.taxonomyAnalysisDir == null) {
			throw new IOException("taxonomyAnalysisDir is null");
		}
		if (request.genomeIntensityMap == null) {
			throw new IOException("genomeIntensityMap is null");
		}
	}

	private static List<GenomeRow> buildRows(Request request) {
		ArrayList<GenomeRow> rows = new ArrayList<>();
		double totalIntensity = 0.0;
		for (double[] values : request.genomeIntensityMap.values()) {
			totalIntensity += sum(values);
		}
		for (String genome : request.genomeIntensityMap.keySet()) {
			GenomeRow row = new GenomeRow();
			row.genome = genome;
			row.intensity = sum(request.genomeIntensityMap.get(genome));
			row.relativeAbundance = totalIntensity > 0.0 ? row.intensity / totalIntensity : 0.0;
			row.totalPeptides = sizeOf(request.totalPepCountMap, genome);
			row.razorPeptides = sizeOf(request.razorPepCountMap, genome);
			row.uniquePeptides = sizeOf(request.uniquePepCountMap, genome);
			row.psmCount = request.genomePsmCountMap != null ? request.genomePsmCountMap.getOrDefault(genome, 0) : 0;
			addAll(row.coveragePeptides, request.totalPepCountMap, genome);
			addAll(row.uniquePeptideSet, request.uniquePepCountMap, genome);
			if (row.coveragePeptides.isEmpty()) {
				addAll(row.coveragePeptides, request.razorPepCountMap, genome);
			}
			if (row.coveragePeptides.isEmpty()) {
				addAll(row.coveragePeptides, request.uniquePepCountMap, genome);
			}
			rows.add(row);
		}
		return rows;
	}

	/** In-memory final selection, identical to genome_final_selected.tsv. Zero-gain
	 * leftovers are evidence, not selected genomes. */
	public static java.util.Set<String> selectGenomes(Request request) {
		List<GenomeRow> rows = buildRows(request);
		rows.sort(abundanceComparator());
		Ranking ranking = rank(rows, request.privilegeUniquePeptideFraction);
		java.util.Set<String> selected = new java.util.LinkedHashSet<>(ranking.privileged);
		selected.addAll(ranking.greedyOrder);
		return java.util.Collections.unmodifiableSet(selected);
	}

	private static Ranking rank(List<GenomeRow> rows, double privilegeFraction) {
		Ranking ranking = new Ranking();
		HashMap<String, GenomeRow> byGenome = new HashMap<>();
		HashSet<String> totalUniquePeptides = new HashSet<>();
		for (GenomeRow row : rows) {
			byGenome.put(row.genome, row);
			totalUniquePeptides.addAll(row.uniquePeptideSet);
		}
		int minUnique = Math.max(1, (int) Math.ceil(totalUniquePeptides.size() * privilegeFraction));

		HashSet<String> covered = new HashSet<>();
		for (GenomeRow row : rows) {
			if (row.uniquePeptides >= minUnique) {
				ranking.privileged.add(row.genome);
				ranking.reason.put(row.genome, "PRIVILEGE_UNIQUE_PEPTIDE_1PCT");
				covered.addAll(row.coveragePeptides);
			}
		}

		LinkedHashSet<String> remaining = new LinkedHashSet<>();
		for (GenomeRow row : rows) {
			if (!ranking.privileged.contains(row.genome)) {
				remaining.add(row.genome);
			}
		}

		while (!remaining.isEmpty()) {
			String best = null;
			int bestGain = -1;
			for (String genome : remaining) {
				GenomeRow row = byGenome.get(genome);
				int gain = incrementalGain(covered, row.coveragePeptides);
				if (best == null || gain > bestGain || (gain == bestGain
						&& abundanceComparator().compare(row, byGenome.get(best)) < 0)) {
					best = genome;
					bestGain = gain;
				}
			}
			if (best == null || bestGain <= 0) {
				break;
			}
			ranking.greedyOrder.add(best);
			ranking.greedyGain.put(best, bestGain);
			covered.addAll(byGenome.get(best).coveragePeptides);
			remaining.remove(best);
		}

		ranking.leftovers.addAll(remaining);
		ranking.leftovers.sort((a, b) -> abundanceComparator().compare(byGenome.get(a), byGenome.get(b)));
		int knee = computeKnee(ranking.greedyOrder, ranking.greedyGain);
		ranking.focusKnee = knee;
		int tailLength = Math.max(0, ranking.greedyOrder.size() - knee);
		ranking.comprehensiveKnee = Math.min(ranking.greedyOrder.size(),
				knee + (int) Math.ceil(tailLength * 0.5));
		return ranking;
	}

	private static final class Ranking {
		LinkedHashSet<String> privileged = new LinkedHashSet<>();
		ArrayList<String> greedyOrder = new ArrayList<>();
		ArrayList<String> leftovers = new ArrayList<>();
		HashMap<String, String> reason = new HashMap<>();
		HashMap<String, Integer> greedyGain = new HashMap<>();
		int focusKnee;
		int comprehensiveKnee;
	}

	private static void writeEvidence(File out, List<GenomeRow> rows) throws IOException {
		try (PrintWriter writer = new PrintWriter(out)) {
			writer.println("Genome\tIntensity\tRelativeAbundance\tTotalPeptides\tRazorPeptides\tUniquePeptides\tTotalPSM");
			for (GenomeRow row : rows) {
				writer.printf(Locale.ROOT, "%s\t%s\t%.8f\t%d\t%d\t%d\t%d%n", row.genome,
						FormatTool.getDFE4().format(row.intensity), row.relativeAbundance, row.totalPeptides,
						row.razorPeptides, row.uniquePeptides, row.psmCount);
			}
		}
	}

	private static void writeRanked(File out, List<GenomeRow> rows, Ranking ranking) throws IOException {
		HashMap<String, GenomeRow> byGenome = new HashMap<>();
		for (GenomeRow row : rows) {
			byGenome.put(row.genome, row);
		}
		ArrayList<String> order = new ArrayList<>();
		order.addAll(ranking.privileged);
		order.addAll(ranking.greedyOrder);
		order.addAll(ranking.leftovers);

		HashSet<String> covered = new HashSet<>();
		try (PrintWriter writer = new PrintWriter(out)) {
			writer.println("Rank\tGenome\tIntensity\tRelativeAbundance\tTotalPeptides\tRazorPeptides\tUniquePeptides\tTotalPSM\tReason\tGreedyGain\tGreedyStep\tGreedyRatio\tGreedyMA\tKneeStep\tKneeKneedle\tKneeCurvature\tKneePiecewise\tKneeSlope\tTotalCoveredPeptides\tCoverage\tIncrement");
			int rank = 0;
			int greedyStep = 0;
			int totalPeptides = totalCoveragePeptides(rows);
			ArrayList<Double> ratios = new ArrayList<>();
			for (String genome : order) {
				GenomeRow row = byGenome.get(genome);
				if (row == null) {
					continue;
				}
				int before = covered.size();
				covered.addAll(row.coveragePeptides);
				int increment = covered.size() - before;
				int gain = ranking.greedyGain.getOrDefault(genome, 0);
				String reason = ranking.reason.getOrDefault(genome, gain > 0 ? "GREEDY" : "");
				String stepText = "";
				String ratioText = "";
				String maText = "";
				if (gain > 0) {
					greedyStep++;
					double ratio = totalPeptides > 0 ? (double) gain / (double) totalPeptides : 0.0;
					ratios.add(ratio);
					stepText = String.valueOf(greedyStep);
					ratioText = String.format(Locale.ROOT, "%.6f", ratio);
					maText = movingAverageText(ratios, 5);
				}
				double coverage = totalPeptides > 0 ? (double) covered.size() / (double) totalPeptides : 1.0;
				writer.printf(Locale.ROOT,
						"%d\t%s\t%s\t%.8f\t%d\t%d\t%d\t%d\t%s\t%d\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%.6f\t%d%n",
						rank++, row.genome, FormatTool.getDFE4().format(row.intensity), row.relativeAbundance,
						row.totalPeptides, row.razorPeptides, row.uniquePeptides, row.psmCount, reason, gain,
						stepText, ratioText, maText, ranking.focusKnee, ranking.focusKnee, ranking.focusKnee,
						ranking.comprehensiveKnee, ranking.focusKnee, covered.size(), coverage, increment);
			}
		}
	}

	private static void writeSelected(File out, Ranking ranking) throws IOException {
		try (PrintWriter writer = new PrintWriter(out)) {
			writer.println("Genome\tReason");
			for (String genome : ranking.privileged) {
				writer.println(genome + "\t" + ranking.reason.getOrDefault(genome, ""));
			}
			int step = 0;
			for (String genome : ranking.greedyOrder) {
				step++;
				String label = step <= ranking.focusKnee ? "core"
						: (step <= ranking.comprehensiveKnee ? "mid" : "tail");
				writer.println(genome + "\tGREEDY_" + label);
			}
			writer.println();
			writer.println("#PrivilegeGenomes\t" + ranking.privileged.size());
			writer.println("#GreedyGenomes\t" + ranking.greedyOrder.size());
			writer.println("#FocusKnee\t" + ranking.focusKnee);
			writer.println("#ComprehensiveKnee\t" + ranking.comprehensiveKnee);
		}
	}

	private static Comparator<GenomeRow> abundanceComparator() {
		return (a, b) -> {
			int cmp = Double.compare(b.intensity, a.intensity);
			if (cmp != 0) return cmp;
			cmp = Integer.compare(b.razorPeptides, a.razorPeptides);
			if (cmp != 0) return cmp;
			cmp = Integer.compare(b.uniquePeptides, a.uniquePeptides);
			if (cmp != 0) return cmp;
			cmp = Integer.compare(b.psmCount, a.psmCount);
			if (cmp != 0) return cmp;
			return a.genome.compareTo(b.genome);
		};
	}

	private static int computeKnee(List<String> greedyOrder, Map<String, Integer> gainMap) {
		int n = greedyOrder.size();
		if (n <= 0) {
			return 0;
		}
		if (n < 5) {
			return n;
		}
		double[] y = new double[n];
		int maxGain = 1;
		for (String genome : greedyOrder) {
			maxGain = Math.max(maxGain, gainMap.getOrDefault(genome, 0));
		}
		for (int i = 0; i < n; i++) {
			y[i] = (double) gainMap.getOrDefault(greedyOrder.get(i), 0) / (double) maxGain;
		}
		double bestDistance = -1.0;
		int bestStep = n;
		double x1 = 0.0, y1 = y[0], x2 = 1.0, y2 = y[n - 1];
		double denom = Math.sqrt(Math.pow(y2 - y1, 2.0) + Math.pow(x2 - x1, 2.0));
		for (int i = 1; i < n - 1; i++) {
			double x = (double) i / (double) (n - 1);
			double distance = Math.abs((y2 - y1) * x - (x2 - x1) * y[i] + x2 * y1 - y2 * x1) / denom;
			if (distance > bestDistance) {
				bestDistance = distance;
				bestStep = i + 1;
			}
		}
		return bestStep;
	}

	private static String movingAverageText(List<Double> values, int window) {
		if (values.size() < window) {
			return "";
		}
		double sum = 0.0;
		for (int i = values.size() - window; i < values.size(); i++) {
			sum += values.get(i);
		}
		return String.format(Locale.ROOT, "%.6f", sum / (double) window);
	}

	private static int totalCoveragePeptides(List<GenomeRow> rows) {
		HashSet<String> total = new HashSet<>();
		for (GenomeRow row : rows) {
			total.addAll(row.coveragePeptides);
		}
		return total.size();
	}

	private static double totalIntensity(List<GenomeRow> rows) {
		double total = 0.0;
		for (GenomeRow row : rows) {
			total += row.intensity;
		}
		return total;
	}

	private static int incrementalGain(Set<String> covered, Set<String> candidate) {
		int gain = 0;
		for (String peptide : candidate) {
			if (!covered.contains(peptide)) {
				gain++;
			}
		}
		return gain;
	}

	private static int sizeOf(Map<String, ? extends Set<String>> map, String key) {
		if (map == null || !map.containsKey(key) || map.get(key) == null) {
			return 0;
		}
		return map.get(key).size();
	}

	private static void addAll(Set<String> target, Map<String, ? extends Set<String>> map, String key) {
		if (map != null && map.get(key) != null) {
			target.addAll(map.get(key));
		}
	}

	private static double sum(double[] values) {
		double sum = 0.0;
		if (values != null) {
			for (double value : values) {
				sum += value;
			}
		}
		return sum;
	}
}
