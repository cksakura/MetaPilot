package uk.ac.quadram.metapilot.task.mag;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

/**
 * One row per raw file describing the genome-selection cut.
 */
public final class GenomeSelectionSummary {

	public static final String FILE_NAME = "genome_selection_summary.tsv";
	public static final String HEADER = "RawFile\tMarkerSet\tStrategy\tGenomesSelected\tPsmCoverageReached\tStopReason"
			+ "\tFocusGenomes\tFocusCoverage\tComprehensiveGenomes\tComprehensiveCoverage"
			+ "\tGenomesFor0.80\tGenomesFor0.90\tGenomesFor0.95\tGenomesFor0.99"
			+ "\tPsmCoverageThreshold\tStepGuard\tMaxGenomes";

	private GenomeSelectionSummary() {
	}

	public static final class Row {
		public String rawFile = "";
		public String markerSet = "";
		public String strategy = "";
		public int genomesSelected = -1;
		public double psmCoverageReached = Double.NaN;
		public String stopReason = "";
		public int focusGenomes = -1;
		public double focusCoverage = Double.NaN;
		public int comprehensiveGenomes = -1;
		public double comprehensiveCoverage = Double.NaN;
		public int genomesFor080 = -1;
		public int genomesFor090 = -1;
		public int genomesFor095 = -1;
		public int genomesFor099 = -1;
		public double psmCoverageThreshold = Double.NaN;
		public double stepGuard = Double.NaN;
		public int maxGenomes = -1;
	}

	public static final class CurvePoint {
		public final int genomes;
		public final double coverage;
		public final boolean selected;

		public CurvePoint(int genomes, double coverage, boolean selected) {
			this.genomes = genomes;
			this.coverage = coverage;
			this.selected = selected;
		}
	}

	public static void writeSample(File file, String rawFile, MarkerGenomeSelector.SelectionResult result)
			throws IOException {
		Row row = fromResult(rawFile, result);
		File parent = file.getParentFile();
		if (parent != null && !parent.exists()) {
			parent.mkdirs();
		}
		try (BufferedWriter writer = new BufferedWriter(new FileWriter(file))) {
			writer.write(HEADER);
			writer.newLine();
			writer.write(format(row));
			writer.newLine();
		}
	}

	/**
	 * Writes one combined table from each sample's marker-folder summary.
	 * Samples without a summary are omitted.
	 */
	public static void combine(File[] sampleDirs, String markerDirName, File destination) throws IOException {
		Map<String, String> lines = new LinkedHashMap<>();
		if (sampleDirs != null) {
			for (File sample : sampleDirs) {
				if (sample == null) {
					continue;
				}
				File summary = new File(new File(sample, markerDirName), FILE_NAME);
				if (!summary.isFile()) {
					continue;
				}
				for (Row row : read(summary)) {
					if (row.rawFile != null && !row.rawFile.isEmpty()) {
						lines.put(row.rawFile, format(row));
					}
				}
			}
		}
		File parent = destination.getParentFile();
		if (parent != null && !parent.exists()) {
			parent.mkdirs();
		}
		try (BufferedWriter writer = new BufferedWriter(new FileWriter(destination))) {
			writer.write(HEADER);
			writer.newLine();
			for (String line : lines.values()) {
				writer.write(line);
				writer.newLine();
			}
		}
	}

	public static List<Row> read(File file) throws IOException {
		List<Row> rows = new ArrayList<>();
		if (file == null || !file.isFile()) {
			return rows;
		}
		try (BufferedReader reader = new BufferedReader(new FileReader(file))) {
			String header = reader.readLine();
			if (header == null) {
				return rows;
			}
			String[] columns = header.split("\t", -1);
			Map<String, Integer> index = new LinkedHashMap<>();
			for (int i = 0; i < columns.length; i++) {
				index.put(columns[i].trim(), i);
			}
			String line;
			while ((line = reader.readLine()) != null) {
				line = line.trim();
				if (line.isEmpty() || line.startsWith("#")) {
					continue;
				}
				String[] cells = line.split("\t", -1);
				Row row = new Row();
				row.rawFile = cell(cells, index, "RawFile");
				row.markerSet = cell(cells, index, "MarkerSet");
				row.strategy = cell(cells, index, "Strategy");
				row.genomesSelected = cellInt(cells, index, "GenomesSelected");
				row.psmCoverageReached = cellDouble(cells, index, "PsmCoverageReached");
				row.stopReason = cell(cells, index, "StopReason");
				row.focusGenomes = cellInt(cells, index, "FocusGenomes");
				row.focusCoverage = cellDouble(cells, index, "FocusCoverage");
				row.comprehensiveGenomes = cellInt(cells, index, "ComprehensiveGenomes");
				row.comprehensiveCoverage = cellDouble(cells, index, "ComprehensiveCoverage");
				row.genomesFor080 = cellInt(cells, index, "GenomesFor0.80");
				row.genomesFor090 = cellInt(cells, index, "GenomesFor0.90");
				row.genomesFor095 = cellInt(cells, index, "GenomesFor0.95");
				row.genomesFor099 = cellInt(cells, index, "GenomesFor0.99");
				row.psmCoverageThreshold = cellDouble(cells, index, "PsmCoverageThreshold");
				row.stepGuard = cellDouble(cells, index, "StepGuard");
				row.maxGenomes = cellInt(cells, index, "MaxGenomes");
				rows.add(row);
			}
		}
		return rows;
	}

	/**
	 * PSM coverage against genomes in rank order. Rows after the greedy cap are omitted:
	 * their coverage values are not a ranking.
	 */
	public static List<CurvePoint> readRankingCurve(File rankedFile) throws IOException {
		List<CurvePoint> points = new ArrayList<>();
		if (rankedFile == null || !rankedFile.isFile()) {
			return points;
		}
		try (BufferedReader reader = new BufferedReader(new FileReader(rankedFile))) {
			String header = reader.readLine();
			if (header == null) {
				return points;
			}
			String[] columns = header.split("\t", -1);
			int coverageIdx = -1;
			int inRankingIdx = -1;
			int cutIdx = -1;
			int gainIdx = -1;
			int reasonIdx = -1;
			for (int i = 0; i < columns.length; i++) {
				String name = columns[i].trim();
				if ("PsmCoverage".equals(name)) {
					coverageIdx = i;
				} else if ("InRanking".equals(name)) {
					inRankingIdx = i;
				} else if ("CutStep".equals(name)) {
					cutIdx = i;
				} else if ("GreedyGain".equals(name)) {
					gainIdx = i;
				} else if ("Reason".equals(name)) {
					reasonIdx = i;
				}
			}
			if (coverageIdx < 0) {
				return points;
			}
			String line;
			int kept = 0;
			while ((line = reader.readLine()) != null) {
				if (line.isEmpty() || line.startsWith("#")) {
					continue;
				}
				String[] cells = line.split("\t", -1);
				if (coverageIdx >= cells.length) {
					continue;
				}
				String reason = reasonIdx >= 0 && reasonIdx < cells.length ? cells[reasonIdx].trim() : "";
				double gain = gainIdx >= 0 && gainIdx < cells.length ? parseDouble(cells[gainIdx]) : Double.NaN;
				boolean inRanking;
				if (inRankingIdx >= 0 && inRankingIdx < cells.length && !cells[inRankingIdx].trim().isEmpty()) {
					inRanking = !"0".equals(cells[inRankingIdx].trim());
				} else {
					inRanking = (reason != null && !reason.isEmpty()) || (!Double.isNaN(gain) && gain > 0);
				}
				if (!inRanking) {
					continue;
				}
				double coverage = parseDouble(cells[coverageIdx]);
				if (Double.isNaN(coverage)) {
					continue;
				}
				boolean selected;
				if (cutIdx >= 0 && cutIdx < cells.length && !cells[cutIdx].trim().isEmpty()) {
					selected = "1".equals(cells[cutIdx].trim());
				} else {
					selected = reason != null && !reason.isEmpty();
				}
				kept++;
				points.add(new CurvePoint(kept, coverage, selected));
			}
		}
		return points;
	}

	public static Row fromResult(String rawFile, MarkerGenomeSelector.SelectionResult result) {
		Row row = new Row();
		row.rawFile = rawFile == null ? "" : rawFile;
		if (result == null) {
			return row;
		}
		row.markerSet = result.markerSet == null ? "" : result.markerSet;
		row.strategy = result.selectionStrategy == null ? "" : result.selectionStrategy;
		row.genomesSelected = result.selected == null ? 0 : result.selected.size();
		row.psmCoverageReached = result.psmCoverageReached;
		row.stopReason = result.stopReason == null ? "" : result.stopReason;
		row.focusGenomes = result.focusGenomes;
		row.focusCoverage = result.focusCoverage;
		row.comprehensiveGenomes = result.comprehensiveGenomes;
		row.comprehensiveCoverage = result.comprehensiveCoverage;
		row.genomesFor080 = result.genomesFor080;
		row.genomesFor090 = result.genomesFor090;
		row.genomesFor095 = result.genomesFor095;
		row.genomesFor099 = result.genomesFor099;
		row.psmCoverageThreshold = parseDouble(result.psmCoverageThresholdText);
		row.stepGuard = parseDouble(result.stepGuardText);
		row.maxGenomes = parseInt(result.maxGenomesText);
		return row;
	}

	static String format(Row row) {
		return text(row.rawFile) + "\t" + text(row.markerSet) + "\t" + text(row.strategy) + "\t"
				+ integer(row.genomesSelected) + "\t" + decimal(row.psmCoverageReached) + "\t" + text(row.stopReason)
				+ "\t" + integer(row.focusGenomes) + "\t" + decimal(row.focusCoverage) + "\t"
				+ integer(row.comprehensiveGenomes) + "\t" + decimal(row.comprehensiveCoverage) + "\t"
				+ integer(row.genomesFor080) + "\t" + integer(row.genomesFor090) + "\t" + integer(row.genomesFor095)
				+ "\t" + integer(row.genomesFor099) + "\t" + decimal(row.psmCoverageThreshold) + "\t"
				+ decimal(row.stepGuard) + "\t" + integer(row.maxGenomes);
	}

	private static String text(String value) {
		return value == null ? "" : value;
	}

	private static String integer(int value) {
		return value < 0 ? "" : Integer.toString(value);
	}

	private static String decimal(double value) {
		return Double.isNaN(value) ? "" : String.format(Locale.ROOT, "%.6f", value);
	}

	private static String cell(String[] cells, Map<String, Integer> index, String name) {
		Integer at = index.get(name);
		if (at == null || at < 0 || at >= cells.length) {
			return "";
		}
		return cells[at].trim();
	}

	private static int cellInt(String[] cells, Map<String, Integer> index, String name) {
		return parseInt(cell(cells, index, name));
	}

	private static double cellDouble(String[] cells, Map<String, Integer> index, String name) {
		return parseDouble(cell(cells, index, name));
	}

	private static int parseInt(String text) {
		if (text == null || text.trim().isEmpty()) {
			return -1;
		}
		try {
			return Integer.parseInt(text.trim());
		} catch (NumberFormatException e) {
			return -1;
		}
	}

	private static double parseDouble(String text) {
		if (text == null || text.trim().isEmpty()) {
			return Double.NaN;
		}
		try {
			return Double.parseDouble(text.trim());
		} catch (NumberFormatException e) {
			return Double.NaN;
		}
	}
}
