package uk.ac.quadram.metapilot.task.dda.mag.gui;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import javafx.scene.Scene;
import javafx.scene.chart.LineChart;
import javafx.scene.chart.NumberAxis;
import javafx.scene.chart.XYChart;
import javafx.scene.control.ComboBox;
import javafx.scene.control.Label;
import javafx.scene.layout.Priority;
import javafx.scene.layout.VBox;
import javafx.stage.Window;
import javafx.util.StringConverter;
import uk.ac.quadram.metapilot.core.mag.MagResultLayoutUtil;
import uk.ac.quadram.metapilot.task.mag.GenomeSelectionCut;
import uk.ac.quadram.metapilot.task.mag.GenomeSelectionSummary;

/**
 * PSM coverage against genomes in rank order, for one raw file at a time.
 */
public final class GenomeSelectionCurveDialog {

	private GenomeSelectionCurveDialog() {
	}

	public static List<String> notes(File resultDirectory) {
		List<String> notes = new ArrayList<>();
		for (GenomeSelectionSummary.Row row : load(resultDirectory)) {
			if (!GenomeSelectionCut.isPsmCoverage(row.strategy)) {
				continue;
			}
			String missed = GenomeSelectionCut.missedTargetMessage(row.psmCoverageThreshold, row.psmCoverageReached,
					row.stopReason, row.maxGenomes);
			if (!missed.isEmpty()) {
				notes.add(row.rawFile + ": " + missed);
			}
		}
		return notes;
	}

	public static void show(Window owner, File resultDirectory) {
		List<GenomeSelectionSummary.Row> rows = load(resultDirectory);
		if (rows.isEmpty()) {
			return;
		}
		ComboBox<GenomeSelectionSummary.Row> files = new ComboBox<>();
		files.getItems().addAll(rows);
		files.setConverter(new StringConverter<>() {
			@Override
			public String toString(GenomeSelectionSummary.Row row) {
				return row == null ? "" : row.rawFile + " (" + row.markerSet + ")";
			}

			@Override
			public GenomeSelectionSummary.Row fromString(String text) {
				return null;
			}
		});
		files.setMaxWidth(Double.MAX_VALUE);

		Label summary = new Label();
		summary.setWrapText(true);
		NumberAxis xAxis = new NumberAxis();
		xAxis.setLabel("Genomes in rank order (log scale)");
		xAxis.setTickLabelFormatter(new StringConverter<>() {
			@Override
			public String toString(Number value) {
				if (value == null) {
					return "";
				}
				double genomes = Math.pow(10.0, value.doubleValue());
				if (genomes < 10) {
					return String.valueOf(Math.max(1, (int) Math.round(genomes)));
				}
				return String.valueOf((int) Math.round(genomes));
			}

			@Override
			public Number fromString(String text) {
				return Math.log10(Double.parseDouble(text));
			}
		});
		NumberAxis yAxis = new NumberAxis(0, 1, 0.1);
		yAxis.setLabel("PSM coverage");
		LineChart<Number, Number> chart = new LineChart<>(xAxis, yAxis);
		chart.setAnimated(false);
		chart.setCreateSymbols(true);
		chart.setLegendVisible(true);
		VBox.setVgrow(chart, Priority.ALWAYS);

		files.setOnAction(e -> draw(resultDirectory, files.getValue(), chart, summary));
		files.getSelectionModel().selectFirst();
		draw(resultDirectory, files.getValue(), chart, summary);

		VBox root = new VBox(8, files, summary, chart);
		root.setStyle("-fx-padding: 12;");
		javafx.stage.Stage stage = new javafx.stage.Stage();
		if (owner != null) {
			stage.initOwner(owner);
		}
		stage.setTitle("Genome selection");
		stage.setScene(new Scene(root, 760, 520));
		stage.show();
	}

	private static void draw(File resultDirectory, GenomeSelectionSummary.Row row, LineChart<Number, Number> chart,
			Label summary) {
		chart.getData().clear();
		if (row == null) {
			summary.setText("");
			return;
		}
		summary.setText(describe(row));
		List<GenomeSelectionSummary.CurvePoint> curve = List.of();
		try {
			curve = GenomeSelectionSummary.readRankingCurve(rankedFile(resultDirectory, row));
		} catch (IOException ignored) {
		}
		XYChart.Series<Number, Number> coverage = new XYChart.Series<>();
		coverage.setName("PSM coverage");
		for (GenomeSelectionSummary.CurvePoint point : curve) {
			coverage.getData().add(new XYChart.Data<>(Math.log10(point.genomes), point.coverage));
		}
		chart.getData().add(coverage);
		addPoint(chart, "Cut", row.genomesSelected, row.psmCoverageReached);
		addPoint(chart, "Focus", row.focusGenomes, row.focusCoverage);
		addPoint(chart, "Comprehensive", row.comprehensiveGenomes, row.comprehensiveCoverage);
		if (GenomeSelectionCut.isPsmCoverage(row.strategy) && !Double.isNaN(row.psmCoverageThreshold) && !curve.isEmpty()) {
			XYChart.Series<Number, Number> target = new XYChart.Series<>();
			target.setName("Target");
			double first = Math.log10(curve.get(0).genomes);
			double last = Math.log10(curve.get(curve.size() - 1).genomes);
			target.getData().add(new XYChart.Data<>(first, row.psmCoverageThreshold));
			target.getData().add(new XYChart.Data<>(last, row.psmCoverageThreshold));
			chart.getData().add(target);
		}
	}

	private static void addPoint(LineChart<Number, Number> chart, String name, int genomes, double coverage) {
		if (genomes < 1 || Double.isNaN(coverage)) {
			return;
		}
		XYChart.Series<Number, Number> series = new XYChart.Series<>();
		series.setName(name);
		series.getData().add(new XYChart.Data<>(Math.log10(genomes), coverage));
		chart.getData().add(series);
	}

	private static String describe(GenomeSelectionSummary.Row row) {
		String coverage = Double.isNaN(row.psmCoverageReached) ? "n/a"
				: String.format(java.util.Locale.ROOT, "%.1f%%", row.psmCoverageReached * 100.0);
		StringBuilder text = new StringBuilder();
		text.append(row.genomesSelected < 0 ? "n/a" : row.genomesSelected).append(" genomes, coverage ")
				.append(coverage).append(", ").append(row.stopReason == null ? "" : row.stopReason).append(".");
		if (GenomeSelectionCut.isPsmCoverage(row.strategy)) {
			String missed = GenomeSelectionCut.missedTargetMessage(row.psmCoverageThreshold, row.psmCoverageReached,
					row.stopReason, row.maxGenomes);
			if (!missed.isEmpty()) {
				text.append(" ").append(missed);
			}
		}
		return text.toString();
	}

	private static List<GenomeSelectionSummary.Row> load(File resultDirectory) {
		if (resultDirectory == null) {
			return List.of();
		}
		File summary = new File(new File(resultDirectory, "pfind_MAG_search"), GenomeSelectionSummary.FILE_NAME);
		try {
			return GenomeSelectionSummary.read(summary);
		} catch (IOException e) {
			return List.of();
		}
	}

	private static File rankedFile(File resultDirectory, GenomeSelectionSummary.Row row) {
		String prefix = "Marker4".equalsIgnoreCase(row.markerSet) ? "marker_genome" : "genome";
		File marker = new File(new File(new File(new File(resultDirectory, "pfind_MAG_search"),
				MagResultLayoutUtil.SEARCH_RESULTS_DIR_NAME), row.rawFile), "marker");
		File ranked = new File(marker, prefix + "_ranked.tsv");
		if (ranked.isFile()) {
			return ranked;
		}
		return new File(new File(new File(new File(resultDirectory, "pfind_MAG_search"), row.rawFile), "marker"),
				prefix + "_ranked.tsv");
	}
}
