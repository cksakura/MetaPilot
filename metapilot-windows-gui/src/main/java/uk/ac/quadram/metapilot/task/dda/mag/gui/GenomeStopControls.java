package uk.ac.quadram.metapilot.task.dda.mag.gui;

import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.control.Label;
import javafx.scene.control.RadioButton;
import javafx.scene.control.Separator;
import javafx.scene.control.Spinner;
import javafx.scene.control.SpinnerValueFactory;
import javafx.scene.control.TitledPane;
import javafx.scene.control.ToggleGroup;
import javafx.scene.control.Tooltip;
import javafx.util.StringConverter;
import javafx.scene.layout.ColumnConstraints;
import javafx.scene.layout.FlowPane;
import javafx.scene.layout.GridPane;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Region;
import javafx.scene.layout.VBox;
import uk.ac.quadram.metapilot.task.mag.GenomeSelectionCut;
import uk.ac.quadram.metapilot.task.mag.MarkerGenomeSelector;

/**
 * "When to stop adding genomes" for Marker2 and Marker4.
 * Construction does not notify the listener; the rest of the form may not exist yet.
 */
public final class GenomeStopControls extends GridPane {

	private final Runnable onChange;
	private boolean notify = true;

	private final RadioButton automaticRadio;
	private final RadioButton focusRadio;
	private final RadioButton comprehensiveRadio;
	private final RadioButton coverageRadio;
	private final Spinner<Double> targetCoverageSpinner;
	private final Spinner<Integer> stepGuardSpinner;
	private final Spinner<Integer> maxGenomesSpinner;
	private final Label warningLabel;
	private final TitledPane advanced;

	public GenomeStopControls(Runnable onChange) {
		this.onChange = onChange;
		setHgap(15);
		setVgap(6);
		setPadding(new Insets(4, 0, 0, 0));

		ColumnConstraints left = new ColumnConstraints();
		left.setPercentWidth(50);
		ColumnConstraints right = new ColumnConstraints();
		right.setPercentWidth(50);
		getColumnConstraints().addAll(left, right);

		Label title = new Label("When to stop adding genomes");
		title.setStyle("-fx-font-weight: bold;");

		ToggleGroup stopGroup = new ToggleGroup();
		automaticRadio = new RadioButton("Automatic");
		automaticRadio.setToggleGroup(stopGroup);
		automaticRadio.setSelected(true);
		automaticRadio.setWrapText(true);
		automaticRadio.setTooltip(new Tooltip(
				"Detects the tail of the curve, where adding genomes stops paying off. Focus favours depth on the main genomes; Comprehensive adds more low-abundance genomes."));

		ToggleGroup modeGroup = new ToggleGroup();
		focusRadio = new RadioButton("Focus");
		comprehensiveRadio = new RadioButton("Comprehensive");
		focusRadio.setToggleGroup(modeGroup);
		comprehensiveRadio.setToggleGroup(modeGroup);
		focusRadio.setSelected(true);
		focusRadio.setTooltip(automaticRadio.getTooltip());
		comprehensiveRadio.setTooltip(automaticRadio.getTooltip());
		FlowPane modeBox = new FlowPane(8, 4, new Label("Mode"), focusRadio, comprehensiveRadio);

		coverageRadio = new RadioButton("Stop at target PSM coverage");
		coverageRadio.setToggleGroup(stopGroup);
		coverageRadio.setWrapText(true);
		coverageRadio.setTooltip(new Tooltip(
				"Selects genomes until this share of the marker PSMs is explained. Higher values add more genomes, and larger databases can lose low-abundance proteins from the other genomes."));

		targetCoverageSpinner = new Spinner<>(new SpinnerValueFactory.DoubleSpinnerValueFactory(0.50, 1.00, 0.95, 0.01));
		targetCoverageSpinner.setEditable(true);
		targetCoverageSpinner.setPrefWidth(90);
		targetCoverageSpinner.setTooltip(coverageRadio.getTooltip());
		HBox targetBox = new HBox(8, new Label("Target coverage"), targetCoverageSpinner);
		targetBox.setAlignment(Pos.CENTER_LEFT);

		int defaultSteps = GenomeSelectionCut.stepGuardPercentSteps(GenomeSelectionCut.DEFAULT_STEP_GUARD);
		int maxSteps = GenomeSelectionCut.stepGuardPercentSteps(GenomeSelectionCut.MAX_STEP_GUARD);
		SpinnerValueFactory.IntegerSpinnerValueFactory stepGuardValues =
				new SpinnerValueFactory.IntegerSpinnerValueFactory(0, maxSteps, defaultSteps, 1);
		stepGuardValues.setConverter(stepGuardPercentConverter());
		stepGuardSpinner = new Spinner<>(stepGuardValues);
		stepGuardSpinner.setEditable(true);
		stepGuardSpinner.setPrefWidth(100);
		stepGuardSpinner.setTooltip(new Tooltip(
				"Stops when further genomes add almost nothing, even if the target is not reached. "
						+ "Shown as a percentage of marker PSMs: 0.1% is saved as 0.001. "
						+ "Range 0–0.2%, in steps of 0.01 percentage points. Set to 0% to reach the target exactly."));
		maxGenomesSpinner = new Spinner<>(new SpinnerValueFactory.IntegerSpinnerValueFactory(1, 100000, 300, 1));
		maxGenomesSpinner.setEditable(true);
		maxGenomesSpinner.setPrefWidth(100);
		maxGenomesSpinner.setTooltip(new Tooltip("Upper limit on genomes per raw file."));
		GridPane advancedFields = new GridPane();
		advancedFields.setHgap(8);
		advancedFields.setVgap(6);
		advancedFields.addRow(0, new Label("Step guard"), stepGuardSpinner);
		advancedFields.addRow(1, new Label("Max genomes"), maxGenomesSpinner);
		advanced = new TitledPane("Advanced", advancedFields);
		advanced.setExpanded(false);
		advanced.setMaxWidth(Region.USE_PREF_SIZE);
		advanced.visibleProperty().bind(coverageRadio.selectedProperty());
		advanced.managedProperty().bind(advanced.visibleProperty());

		warningLabel = new Label();
		warningLabel.setWrapText(true);
		warningLabel.setMaxWidth(280);
		warningLabel.setStyle("-fx-text-fill: #8a5a00;");
		warningLabel.visibleProperty().bind(warningLabel.textProperty().isNotEmpty());
		warningLabel.managedProperty().bind(warningLabel.visibleProperty());

		modeBox.setPadding(new Insets(0, 0, 0, 18));
		VBox coverageDetails = new VBox(6, targetBox, advanced, warningLabel);
		coverageDetails.setPadding(new Insets(0, 0, 0, 18));
		add(new Separator(), 0, 0, 2, 1);
		add(title, 0, 1, 2, 1);
		add(automaticRadio, 0, 2);
		add(coverageRadio, 1, 2);
		add(modeBox, 0, 3);
		add(coverageDetails, 1, 3);
		GridPane.setValignment(modeBox, javafx.geometry.VPos.TOP);
		GridPane.setValignment(coverageDetails, javafx.geometry.VPos.TOP);

		automaticRadio.setOnAction(e -> changed());
		coverageRadio.setOnAction(e -> changed());
		focusRadio.setOnAction(e -> changed());
		comprehensiveRadio.setOnAction(e -> changed());
		targetCoverageSpinner.valueProperty().addListener((obs, oldValue, value) -> {
			roundTarget(value);
			changed();
		});
		stepGuardSpinner.valueProperty().addListener((obs, oldValue, value) -> changed());
		maxGenomesSpinner.valueProperty().addListener((obs, oldValue, value) -> changed());
		syncEnabled();
	}

	public String strategy() {
		commitSpinners();
		if (coverageRadio.isSelected()) {
			return MarkerGenomeSelector.KNEE_MODE_PSM_COVERAGE;
		}
		if (comprehensiveRadio.isSelected()) {
			return MarkerGenomeSelector.KNEE_MODE_COMPREHENSIVE;
		}
		return MarkerGenomeSelector.KNEE_MODE_FOCUS;
	}

	public double psmCoverage() {
		commitSpinners();
		return targetCoverageSpinner.getValue() == null ? GenomeSelectionCut.DEFAULT_PSM_COVERAGE
				: targetCoverageSpinner.getValue();
	}

	public double stepGuard() {
		commitSpinners();
		Integer steps = stepGuardSpinner.getValue();
		return steps == null ? GenomeSelectionCut.DEFAULT_STEP_GUARD
				: GenomeSelectionCut.stepGuardFromPercentSteps(steps);
	}

	public int maxGenomes() {
		commitSpinners();
		return maxGenomesSpinner.getValue() == null ? GenomeSelectionCut.DEFAULT_MAX_GENOMES : maxGenomesSpinner.getValue();
	}

	public String summary() {
		if (coverageRadio.isSelected()) {
			return "At a PSM coverage (target " + psmCoverage() + ", step guard "
					+ GenomeSelectionCut.formatStepGuardPercent(
							GenomeSelectionCut.stepGuardPercentSteps(stepGuard()))
					+ ", max genomes " + maxGenomes() + ")";
		}
		if (comprehensiveRadio.isSelected()) {
			return "Automatically, Comprehensive";
		}
		return "Automatically, Focus";
	}

	public void apply(String strategy, double psmCoverage, double stepGuard, int maxGenomes) {
		notify = false;
		try {
			if (GenomeSelectionCut.isPsmCoverage(strategy)) {
				coverageRadio.setSelected(true);
			} else if (GenomeSelectionCut.isComprehensive(strategy)) {
				automaticRadio.setSelected(true);
				comprehensiveRadio.setSelected(true);
			} else {
				automaticRadio.setSelected(true);
				focusRadio.setSelected(true);
			}
			setDouble(targetCoverageSpinner, psmCoverage);
			if (stepGuardSpinner.getValueFactory() != null) {
				stepGuardSpinner.getValueFactory().setValue(GenomeSelectionCut.stepGuardPercentSteps(stepGuard));
			}
			if (maxGenomesSpinner.getValueFactory() != null) {
				maxGenomesSpinner.getValueFactory().setValue(Math.max(1, maxGenomes));
			}
			syncEnabled();
		} finally {
			notify = true;
		}
	}

	private void changed() {
		syncEnabled();
		if (notify && onChange != null) {
			onChange.run();
		}
	}

	private void syncEnabled() {
		boolean coverage = coverageRadio.isSelected();
		focusRadio.setDisable(coverage);
		comprehensiveRadio.setDisable(coverage);
		targetCoverageSpinner.setDisable(!coverage);
		stepGuardSpinner.setDisable(!coverage);
		maxGenomesSpinner.setDisable(!coverage);
		if (!coverage || targetCoverageSpinner.getValue() == null || stepGuardSpinner.getValue() == null) {
			warningLabel.setText("");
		} else {
			warningLabel.setText(GenomeSelectionCut.highTargetWarning(targetCoverageSpinner.getValue(),
					GenomeSelectionCut.stepGuardFromPercentSteps(stepGuardSpinner.getValue())));
		}
	}

	private void roundTarget(Double value) {
		if (value == null || targetCoverageSpinner.getValueFactory() == null) {
			return;
		}
		double rounded = Math.round(value * 100.0) / 100.0;
		if (Math.abs(rounded - value) > 1e-9) {
			notify = false;
			try {
				targetCoverageSpinner.getValueFactory().setValue(rounded);
			} finally {
				notify = true;
			}
		}
	}

	/**
	 * Typed spinner text is not the value until Enter or focus loss. Collecting
	 * for a run or a parameter file has to commit first, or the JSON keeps the
	 * previous number.
	 */
	private void commitSpinners() {
		commit(targetCoverageSpinner);
		commit(stepGuardSpinner);
		commit(maxGenomesSpinner);
	}

	private static void commit(Spinner<?> spinner) {
		if (spinner == null || !spinner.isEditable() || spinner.getEditor() == null
				|| spinner.getValueFactory() == null || spinner.getValueFactory().getConverter() == null) {
			return;
		}
		String text = spinner.getEditor().getText();
		if (text == null) {
			return;
		}
		try {
			Object parsed = spinner.getValueFactory().getConverter().fromString(text);
			@SuppressWarnings("unchecked")
			SpinnerValueFactory<Object> factory = (SpinnerValueFactory<Object>) spinner.getValueFactory();
			factory.setValue(parsed);
		} catch (RuntimeException e) {
			Object current = spinner.getValue();
			String shown = "";
			if (current != null && spinner.getValueFactory().getConverter() != null) {
				@SuppressWarnings("unchecked")
				StringConverter<Object> converter = (StringConverter<Object>) spinner.getValueFactory().getConverter();
				shown = converter.toString(current);
			}
			spinner.getEditor().setText(shown == null ? "" : shown);
		}
	}

	private static StringConverter<Integer> stepGuardPercentConverter() {
		return new StringConverter<>() {
			@Override
			public String toString(Integer steps) {
				return steps == null ? "" : GenomeSelectionCut.formatStepGuardPercent(steps);
			}

			@Override
			public Integer fromString(String text) {
				return GenomeSelectionCut.parseStepGuardPercent(text);
			}
		};
	}

	private static void setDouble(Spinner<Double> spinner, double value) {
		if (spinner.getValueFactory() == null) {
			return;
		}
		SpinnerValueFactory.DoubleSpinnerValueFactory factory = (SpinnerValueFactory.DoubleSpinnerValueFactory) spinner
				.getValueFactory();
		if (value < factory.getMin()) {
			factory.setMin(value);
		}
		if (value > factory.getMax()) {
			factory.setMax(value);
		}
		factory.setValue(value);
	}
}
