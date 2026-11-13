package uk.ac.quadram.metapilot.task.dia.gui;

import javafx.application.Platform;
import javafx.geometry.Insets;
import javafx.scene.control.*;
import uk.ac.quadram.metapilot.task.common.gui.ParameterSection;
import javafx.scene.control.Alert.AlertType;
import javafx.scene.layout.*;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.scene.text.Font;
import javafx.scene.text.FontWeight;
import javafx.scene.paint.Color;
import uk.ac.quadram.metapilot.task.common.gui.MetaPilotIOPanelFX;
import uk.ac.quadram.metapilot.task.common.gui.MagDatabaseSettingsPanel;
import uk.ac.quadram.metapilot.task.common.gui.MetadataConfigurationPanel;
import uk.ac.quadram.metapilot.task.common.gui.AbstractWorkflowPanel;
import uk.ac.quadram.metapilot.task.par.MetaData;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;
import uk.ac.quadram.metapilot.task.dia.par.DiaMarkerParameterIO;
import uk.ac.quadram.metapilot.task.dia.par.MetaParameterDiaMarker;
import uk.ac.quadram.metapilot.task.dia.DiaMarkerComprehensiveTask;
import uk.ac.quadram.metapilot.task.dia.DiaMarkerCoreTask;
import uk.ac.quadram.metapilot.task.dia.DiaMarkerDeepPerFileTask;
import uk.ac.quadram.metapilot.task.dia.DiaMarkerDeepTogetherTask;
import uk.ac.quadram.metapilot.task.dia.DiaRefinementMode;
import uk.ac.quadram.metapilot.core.runner.DiaMarkerDeepEntrapmentHeadlessRunner;
import uk.ac.quadram.metapilot.task.mag.MagDbItem.MagVersion;
import uk.ac.quadram.metapilot.task.dda.igc.gui.MetaPilotDdaIgcMaxQuantPanel;
import uk.ac.quadram.metapilot.resources.MetaResourceIO;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.File;
import java.util.Optional;

/**
 * DIA Marker workflow panel: marker-based genome selection (Core / Focus or Comprehensive)
 * using Marker4 and DiaMarkerCoreTask.
 */
public class MetaPilotDiaMarkerPanel extends AbstractWorkflowPanel<MetaParameterDiaMarker> {
	private static final Logger LOGGER = LogManager.getLogger(MetaPilotDiaMarkerPanel.class);

	// UI Components
	private MetaPilotIOPanelFX ioPanel;
	private MagDatabaseSettingsPanel magDatabasePanel;
	private MetadataConfigurationPanel metadataPanel;
	// Note: consolePanel is inherited from AbstractWorkflowPanel

	/** Marker workflows: Core/Focus, Comprehensive, Deep Together, Deep Per File. */
	private static final String WORKFLOW_FOCUS = DiaRefinementMode.CORE.getDisplayLabel();
	private static final String WORKFLOW_COMPREHENSIVE = DiaRefinementMode.COMPREHENSIVE.getDisplayLabel();
	private static final String WORKFLOW_DEEP = DiaRefinementMode.DEEP_TOGETHER.getDisplayLabel();
	private static final String WORKFLOW_DEEP_PER_FILE = DiaRefinementMode.DEEP_PER_FILE.getDisplayLabel();
	private static final String REP_CLUSTER_OFF = "Off";
	private static final String REP_CLUSTER_CONSERVATIVE = "Conservative";
	private static final String REP_CLUSTER_BALANCED = "Balanced";
	private static final String REP_CLUSTER_AGGRESSIVE = "Aggressive";
	private static final String REP_CLUSTER_CUSTOM = "Custom";
	// Experimental UI: keep false until representative clustering is ready for release.
	private static final boolean SHOW_REPRESENTATIVE_CLUSTERING_UI = true;

	// Parameter components
	private ComboBox<String> workflowModeComboBox;
	private ComboBox<Integer> threadsComboBox;
	private ComboBox<String> representativeClusteringComboBox;
	private TextField representativeJaccardField;
	private TextField representativeContainmentField;
	private CheckBox representativeDryRunCheckBox;
	private CheckBox entrapmentCheckBox;
	private VBox representativeClusteringBox;
	private VBox representativeCustomThresholdBox;
	private Label marker4StatusLabel;

	// Modification checkboxes (from DiannParameter lines 28-33)
	private CheckBox mExcisionCheckBox;
	private CheckBox carbamiCCheckBox;
	private CheckBox oxidationMCheckBox;
	private CheckBox acetyProNCheckBox;
	private CheckBox phosphoSTYCheckBox;
	private CheckBox kGGCheckBox;

	// Task management
	private DiaMarkerCoreTask currentTask;
	private MetaResourceParameters resourceParams;
	// currentAnalysisParams inherited from AbstractWorkflowPanel; set in startAnalysis for viewReport()


	public MetaPilotDiaMarkerPanel(MagDbItem[] magDbItems, MetaResourceParameters resourceParams) {
		super();
		parameterIO = new DiaMarkerParameterIO();
		this.resourceParams = resourceParams;
		initializeSpecificUI();
		setAvailableMagDatabases(magDbItems);
		loadLastUsedParameters();
		updateMarker4Status();
	}

	// Load available MAG databases
	public void setAvailableMagDatabases(MagDbItem[] magDbItems) {
		if (magDatabasePanel != null) {
			magDatabasePanel.setAvailableMagDatabases(magDbItems);
		}
	}

	private void initializeSpecificUI() {
		// Initialize reusable components
		// Note: consolePanel is already initialized by AbstractWorkflowPanel
		metadataPanel = new MetadataConfigurationPanel(true);
		magDatabasePanel = new MagDatabaseSettingsPanel();
		magDatabasePanel.setHostDatabaseExtensionFilter("DIA host database files",
				"*.fasta", "*.fa", "*.speclib", "*.parquet");

		// Set up callbacks
		metadataPanel.setParameterUpdateCallback(this::updateParameterView);
		magDatabasePanel.setParameterUpdateCallback(() -> {
			updateParameterView();
			updateMarker4Status();
		});
		// consolePanel callbacks are already set by AbstractWorkflowPanel

		// Create tab pane
		TabPane tabPane = new TabPane();
		tabPane.setTabClosingPolicy(TabPane.TabClosingPolicy.UNAVAILABLE);
		tabPane.setMinWidth(980);
		tabPane.setMinHeight(780);
		VBox.setVgrow(tabPane, Priority.ALWAYS);

		// Create tabs
		Tab ioTab = new Tab("Input/Output");
		ioTab.setContent(createIOTab());

		Tab paramTab = new Tab("Parameters");
		paramTab.setContent(createParametersTab());

		Tab consoleTab = new Tab("Console");
		consoleTab.setContent(consolePanel);

		tabPane.getTabs().addAll(ioTab, paramTab, consoleTab);

		// Add tab change listener
		tabPane.getSelectionModel().selectedItemProperty().addListener((obs, oldTab, newTab) -> {
			if (newTab == paramTab) {
				updateMetadataTableFromIOPanel();
				updateParameterView();
				// Update Marker4 status when switching to Parameters tab
					updateMarker4Status();
			} else if (newTab == consoleTab) {
				updateParameterView();
				// Auto-save without validation or alerts (validation/alert only when user clicks Run)
				try {
					MetaParameterDiaMarker currentParams = collectParametersFromUIForAutoSave();
					if (currentParams != null) {
						parameterIO.saveLastUsedParameters(currentParams);
					}
				} catch (Exception e) {
					LOGGER.warn("Failed to auto-save parameters: " + e.getMessage(), e);
				}
			}
		});

		getChildren().add(tabPane);
		updateParameterView();
	}

	@Override
	protected String getWorkflowName() {
		return "DIA Marker";
	}

	@Override
	protected boolean isHpcExportSupported() {
		return true;
	}

	@Override
	public VBox createIOTab() {
		VBox tabContent = new VBox(10);
		tabContent.setPadding(new Insets(10));
		ioPanel = new MetaPilotIOPanelFX();
		tabContent.getChildren().addAll(
			createSectionTitle("Input/Output Settings"),
			ioPanel,
			createSectionTitle("Database Settings"),
			magDatabasePanel
		);
		return tabContent;
	}

	@Override
	public VBox createParametersTab() {
		VBox tabContent = new VBox(10);
		tabContent.setPadding(new Insets(10));
		tabContent.setMinWidth(960);

		// Create main grid layout
		GridPane mainGrid = new GridPane();
		mainGrid.setHgap(15);
		mainGrid.setVgap(15);
		mainGrid.setPadding(new Insets(10));
		mainGrid.setMinWidth(940);

		// Column constraints
		ColumnConstraints col1 = new ColumnConstraints();
		col1.setPercentWidth(50);
		col1.setMinWidth(450);
		col1.setHgrow(Priority.ALWAYS);
		ColumnConstraints col2 = new ColumnConstraints();
		col2.setPercentWidth(50);
		col2.setMinWidth(450);
		col2.setHgrow(Priority.ALWAYS);
		mainGrid.getColumnConstraints().addAll(col1, col2);

		// Row constraints
		RowConstraints row1 = new RowConstraints();
		row1.setVgrow(Priority.ALWAYS);
		mainGrid.getRowConstraints().add(row1);

		// Left column: Performance Settings, Genome Selection, Modifications
		VBox leftColumn = new VBox(15);
		leftColumn.setMinWidth(430);
		TitledPane workflowSection = createWorkflowSection();
		TitledPane performanceSection = createPerformanceSection();
		TitledPane genomeSelectionSection = createGenomeSelectionSection();
		TitledPane modificationsSection = createModificationsSection();
		leftColumn.getChildren().addAll(workflowSection, performanceSection, genomeSelectionSection, modificationsSection);
		VBox.setVgrow(modificationsSection, Priority.ALWAYS);
		GridPane.setVgrow(leftColumn, Priority.ALWAYS);

		// Right column: Metadata Configuration
		VBox rightColumn = new VBox(15);
		rightColumn.setMinWidth(430);
		rightColumn.getChildren().addAll(metadataPanel);
		GridPane.setVgrow(rightColumn, Priority.ALWAYS);

		mainGrid.add(leftColumn, 0, 0);
		mainGrid.add(rightColumn, 1, 0);

		tabContent.getChildren().add(ParameterSection.createViewport(mainGrid));
		return tabContent;
	}

	private TitledPane createWorkflowSection() {
		VBox section = new VBox(10);

		Label workflowLabel = new Label("Workflow:");
		workflowLabel.setStyle("-fx-font-weight: bold;");
		workflowModeComboBox = new ComboBox<>(
				FXCollections.observableArrayList(WORKFLOW_FOCUS, WORKFLOW_COMPREHENSIVE, WORKFLOW_DEEP, WORKFLOW_DEEP_PER_FILE));
		workflowModeComboBox.setMinWidth(160);
		workflowModeComboBox.setPrefWidth(200);
		workflowModeComboBox.setValue(WORKFLOW_FOCUS);
		workflowModeComboBox.setOnAction(e -> {
			updateRepresentativeClusteringAvailability();
			updateEntrapmentAvailability();
			updateParameterView();
		});
		Tooltip workflowTooltip = new Tooltip(
				"Core / Focus: knee-based genome selection with a compact set of genomes.\n"
				+ "Comprehensive: more inclusive genome set using comprehensive knee mode.\n"
				+ "Deep together: deep iterative workflow over all DIA files together.\n"
				+ "Deep per file: run the full deep-marker iterative workflow independently for each raw file (results under marker/<rawname>/...).");
		workflowModeComboBox.setTooltip(workflowTooltip);
		entrapmentCheckBox = new CheckBox("Paired-shuffle entrapment (evaluation)");
		entrapmentCheckBox.setWrapText(true);
		entrapmentCheckBox.setTooltip(new Tooltip(
				"Deep together only. Same processing as dia-marker-deep-entrapment: each round also searches a shuffled copy of that round's genomes.\n"
						+ "Results go to marker_deep_entrapment/ and do not overwrite a production Deep together run.\n"
						+ "Off by default."));
		entrapmentCheckBox.setOnAction(e -> updateParameterView());

		VBox workflowBox = new VBox(5);
		workflowBox.getChildren().addAll(workflowLabel, workflowModeComboBox, entrapmentCheckBox);
		updateEntrapmentAvailability();

		Label representativeClusteringLabel = new Label("Deep per-file clustering:");
		representativeClusteringLabel.setStyle("-fx-font-weight: bold;");
		representativeClusteringComboBox = new ComboBox<>(FXCollections.observableArrayList(
				REP_CLUSTER_OFF, REP_CLUSTER_CONSERVATIVE, REP_CLUSTER_BALANCED, REP_CLUSTER_AGGRESSIVE,
				REP_CLUSTER_CUSTOM));
		representativeClusteringComboBox.setMinWidth(160);
		representativeClusteringComboBox.setPrefWidth(180);
		representativeClusteringComboBox.setValue(REP_CLUSTER_OFF);
		representativeClusteringComboBox.setOnAction(e -> {
			updateRepresentativeClusteringAvailability();
			updateParameterView();
		});
		representativeClusteringComboBox.setTooltip(new Tooltip(
				"Optional speed-up for large Deep per file cohorts.\n"
				+ "Similar raw files are clustered after Marker4 search; only representative files run iterative rounds.\n"
				+ "Off keeps the original Deep per file strategy."));

		representativeJaccardField = new TextField("0.60");
		representativeJaccardField.setPrefWidth(70);
		representativeJaccardField.setOnAction(e -> updateParameterView());
		representativeJaccardField.focusedProperty().addListener((obs, oldValue, newValue) -> {
			if (!newValue) updateParameterView();
		});
		representativeContainmentField = new TextField("0.78");
		representativeContainmentField.setPrefWidth(70);
		representativeContainmentField.setOnAction(e -> updateParameterView());
		representativeContainmentField.focusedProperty().addListener((obs, oldValue, newValue) -> {
			if (!newValue) updateParameterView();
		});
		HBox customThresholdRow = new HBox(8);
		customThresholdRow.getChildren().addAll(
				new Label("Jaccard:"), representativeJaccardField,
				new Label("Containment:"), representativeContainmentField);
		representativeCustomThresholdBox = new VBox(4);
		representativeCustomThresholdBox.getChildren().add(customThresholdRow);

		representativeDryRunCheckBox = new CheckBox("Grouping dry run");
		representativeDryRunCheckBox.setTooltip(new Tooltip(
				"Run initial Marker searches and write clustering reports, then stop before iterative search."));
		representativeDryRunCheckBox.setOnAction(e -> updateParameterView());

		representativeClusteringBox = new VBox(5);
		representativeClusteringBox.getChildren().addAll(representativeClusteringLabel,
				representativeClusteringComboBox, representativeCustomThresholdBox, representativeDryRunCheckBox);
		updateRepresentativeClusteringAvailability();

		section.getChildren().addAll(workflowBox, representativeClusteringBox);
		return ParameterSection.create("Workflow", section);
	}

	private TitledPane createPerformanceSection() {
		VBox section = new VBox(10);


		Label threadsLabel = new Label("Threads:");
		threadsLabel.setStyle("-fx-font-weight: bold;");
		threadsComboBox = new ComboBox<>();
		threadsComboBox.setMinWidth(100);
		threadsComboBox.setPrefWidth(100);
		threadsComboBox.setMaxWidth(Double.MAX_VALUE);

		ObservableList<Integer> threadOptions = FXCollections.observableArrayList();
		for (int i = 1; i <= Runtime.getRuntime().availableProcessors(); i++) {
			threadOptions.add(i);
		}
		threadsComboBox.setItems(threadOptions);
		threadsComboBox.setValue(4);
		threadsComboBox.valueProperty().addListener((obs, oldVal, newVal) -> updateParameterView());

		VBox threadsBox = new VBox(5);
		threadsBox.getChildren().addAll(threadsLabel, threadsComboBox);

		section.getChildren().addAll(threadsBox);
		return ParameterSection.create("Performance Settings", section);
	}

	private TitledPane createGenomeSelectionSection() {
		VBox section = new VBox(10);


		GridPane paramGrid = new GridPane();
		paramGrid.setHgap(15);
		paramGrid.setVgap(8);

		ColumnConstraints col1 = new ColumnConstraints();
		col1.setMinWidth(150);
		col1.setHgrow(Priority.ALWAYS);
		ColumnConstraints col2 = new ColumnConstraints();
		col2.setMinWidth(150);
		col2.setHgrow(Priority.ALWAYS);
		paramGrid.getColumnConstraints().addAll(col1, col2);

		Label genomeSelectorLabel = new Label("Genome selector:");
		genomeSelectorLabel.setStyle("-fx-font-weight: bold;");
		Label genomeSelectorValueLabel = new Label(MagDbItem.MARKER4);
		genomeSelectorValueLabel.setStyle("-fx-font-weight: normal; -fx-text-fill: #333333;");
		VBox genomeSelectorBox = new VBox(5);
		genomeSelectorBox.getChildren().addAll(genomeSelectorLabel, genomeSelectorValueLabel);

		Label marker4StatusTitleLabel = new Label("Marker4 Status:");
		marker4StatusTitleLabel.setStyle("-fx-font-weight: bold;");
		marker4StatusLabel = new Label();
		marker4StatusLabel.setMinHeight(Label.USE_PREF_SIZE);
		marker4StatusLabel.setStyle("-fx-text-fill: #333333;");
		marker4StatusLabel.setWrapText(true);
		VBox marker4StatusBox = new VBox(5);
		marker4StatusBox.getChildren().addAll(marker4StatusTitleLabel, marker4StatusLabel);

		paramGrid.add(genomeSelectorBox, 0, 0);
		paramGrid.add(marker4StatusBox, 0, 1);
		paramGrid.add(new Region(), 1, 0); // placeholder for alignment

		section.getChildren().addAll(paramGrid);
		return ParameterSection.create("Genome Selection", section);
	}

	private TitledPane createModificationsSection() {
		VBox section = new VBox(10);


		VBox checkboxesBox = new VBox(8);
		
		mExcisionCheckBox = new CheckBox("M excision");
		mExcisionCheckBox.setSelected(true); // Default from DiannParameter
		mExcisionCheckBox.setOnAction(e -> updateParameterView());

		carbamiCCheckBox = new CheckBox("Carbamidomethyl C");
		carbamiCCheckBox.setSelected(true); // Default from DiannParameter
		carbamiCCheckBox.setOnAction(e -> updateParameterView());

		oxidationMCheckBox = new CheckBox("Oxidation M");
		oxidationMCheckBox.setSelected(false); // Default from DiannParameter
		oxidationMCheckBox.setOnAction(e -> updateParameterView());

		acetyProNCheckBox = new CheckBox("Acetylation protein N-term");
		acetyProNCheckBox.setSelected(false); // Default from DiannParameter
		acetyProNCheckBox.setOnAction(e -> updateParameterView());

		phosphoSTYCheckBox = new CheckBox("Phosphorylation STY");
		phosphoSTYCheckBox.setSelected(false); // Default from DiannParameter
		phosphoSTYCheckBox.setOnAction(e -> updateParameterView());

		kGGCheckBox = new CheckBox("K GG");
		kGGCheckBox.setSelected(false); // Default from DiannParameter
		kGGCheckBox.setOnAction(e -> updateParameterView());

		checkboxesBox.getChildren().addAll(
			mExcisionCheckBox,
			carbamiCCheckBox,
			oxidationMCheckBox,
			acetyProNCheckBox,
			phosphoSTYCheckBox,
			kGGCheckBox
		);

		section.getChildren().addAll(checkboxesBox);
		return ParameterSection.create("Modifications", section);
	}

	private Label createSectionTitle(String title) {
		Label label = new Label(title);
		label.setFont(Font.font("System", FontWeight.BOLD, 14));
		label.setTextFill(Color.DARKBLUE);
		return label;
	}

	/**
	 * Update Marker4 status label (sync, no generation). Called when database changes or tab is shown.
	 */
	private void updateMarker4Status() {
		if (marker4StatusLabel == null) {
			return;
		}
		if (magDatabasePanel == null) {
			marker4StatusLabel.setText("—");
			return;
		}
		MagDbItem selectedDb = magDatabasePanel.getSelectedMagDatabase();
		if (selectedDb == null) {
			marker4StatusLabel.setText("Select a MAG database");
			return;
		}
		if (marker4RequiredFilesExist()) {
			marker4StatusLabel.setText("Ready");
			return;
		}
		marker4StatusLabel.setText("Will be generated on first run (may take extra time)");
	}

	/**
	 * Quick check: do Marker4_highpep.fasta and genomes-taxa.tsv exist?
	 */
	private boolean marker4RequiredFilesExist() {
		if (magDatabasePanel == null) {
			return false;
		}
		MagDbItem selectedDb = magDatabasePanel.getSelectedMagDatabase();
		if (selectedDb == null) {
			return false;
		}
		File currentFile = selectedDb.getCurrentFile();
		if (currentFile == null) {
			return false;
		}
		File markersHighpepFile = new File(currentFile, "Marker4_highpep.fasta");
		File genomeTaxaFile = new File(currentFile, "genomes-taxa.tsv");
		return markersHighpepFile.exists() && genomeTaxaFile.exists();
	}

	private void updateMetadataTableFromIOPanel() {
		if (ioPanel == null || metadataPanel == null) {
			return;
		}
		ObservableList<MetaPilotIOPanelFX.ExperimentRow> experiments = ioPanel.getExperimentData();
		metadataPanel.updateMetadataFromIOPanel(experiments);
	}

	private void updateRepresentativeClusteringAvailability() {
		if (representativeClusteringComboBox == null) {
			return;
		}
		if (representativeClusteringBox != null) {
			representativeClusteringBox.setVisible(SHOW_REPRESENTATIVE_CLUSTERING_UI);
			representativeClusteringBox.setManaged(SHOW_REPRESENTATIVE_CLUSTERING_UI);
		}
		if (!SHOW_REPRESENTATIVE_CLUSTERING_UI) {
			representativeClusteringComboBox.setValue(REP_CLUSTER_OFF);
			representativeClusteringComboBox.setDisable(true);
			if (representativeCustomThresholdBox != null) {
				representativeCustomThresholdBox.setVisible(false);
				representativeCustomThresholdBox.setManaged(false);
				representativeCustomThresholdBox.setDisable(true);
			}
			if (representativeDryRunCheckBox != null) {
				representativeDryRunCheckBox.setSelected(false);
				representativeDryRunCheckBox.setVisible(false);
				representativeDryRunCheckBox.setManaged(false);
				representativeDryRunCheckBox.setDisable(true);
			}
			return;
		}
		boolean deepPerFile = workflowModeComboBox != null
				&& WORKFLOW_DEEP_PER_FILE.equals(workflowModeComboBox.getValue());
		representativeClusteringComboBox.setDisable(!deepPerFile);
		if (!deepPerFile) {
			representativeClusteringComboBox.setValue(REP_CLUSTER_OFF);
			if (representativeDryRunCheckBox != null) {
				representativeDryRunCheckBox.setSelected(false);
			}
		}
		String preset = representativeClusteringComboBox.getValue();
		boolean clusteringEnabled = deepPerFile && preset != null && !REP_CLUSTER_OFF.equals(preset);
		boolean customEnabled = clusteringEnabled && REP_CLUSTER_CUSTOM.equals(preset);
		if (representativeCustomThresholdBox != null) {
			representativeCustomThresholdBox.setVisible(customEnabled);
			representativeCustomThresholdBox.setManaged(customEnabled);
			representativeCustomThresholdBox.setDisable(!customEnabled);
		}
		if (representativeDryRunCheckBox != null) {
			if (!clusteringEnabled) {
				representativeDryRunCheckBox.setSelected(false);
			}
			representativeDryRunCheckBox.setVisible(clusteringEnabled);
			representativeDryRunCheckBox.setManaged(clusteringEnabled);
			representativeDryRunCheckBox.setDisable(!clusteringEnabled);
		}
	}

	private void updateEntrapmentAvailability() {
		if (entrapmentCheckBox == null) {
			return;
		}
		boolean deepTogether = workflowModeComboBox != null
				&& WORKFLOW_DEEP.equals(workflowModeComboBox.getValue());
		entrapmentCheckBox.setVisible(deepTogether);
		entrapmentCheckBox.setManaged(deepTogether);
		entrapmentCheckBox.setDisable(!deepTogether);
		if (!deepTogether) {
			entrapmentCheckBox.setSelected(false);
		}
	}

	@Override
	protected void updateParameterView() {
		if (consolePanel == null) {
			return;
		}
		updateRepresentativeClusteringAvailability();
		updateEntrapmentAvailability();
		StringBuilder sb = new StringBuilder();
		sb.append("=== Analysis Parameters ===\n\n");

		if (workflowModeComboBox != null && workflowModeComboBox.getValue() != null) {
			sb.append("Workflow: ").append(workflowModeComboBox.getValue()).append("\n");
		}
		if (entrapmentCheckBox != null && entrapmentCheckBox.isVisible()) {
			sb.append("Paired-shuffle entrapment: ").append(entrapmentCheckBox.isSelected()).append("\n");
		}
		sb.append("\n");

		// Database settings
		sb.append("Database Settings:\n");
		if (magDatabasePanel != null) {
			MagDbItem selectedDb = magDatabasePanel.getSelectedMagDatabase();
			if (selectedDb != null) {
				sb.append("  MAG Database: ").append(selectedDb.getCatalogueID()).append("\n");
			}
			MagVersion selectedVersion = magDatabasePanel.getSelectedMagVersion();
			if (selectedVersion != null) {
				sb.append("  MAG Version: ").append(selectedVersion.getCatalogueVersion()).append("\n");
			}
			sb.append("  Append Host DB: ").append(magDatabasePanel.isAppendHostDb()).append("\n");
			if (magDatabasePanel.isAppendHostDb()
					&& magDatabasePanel.getHostDatabase() != null
					&& !magDatabasePanel.getHostDatabase().isEmpty()) {
				sb.append("  Host Database: ").append(magDatabasePanel.getHostDatabase()).append("\n");
			}
		}
		sb.append("\n");

		// Performance settings
		sb.append("Performance Settings:\n");
		if (threadsComboBox != null) {
			sb.append("  Threads: ").append(threadsComboBox.getValue()).append("\n");
		}
		if (SHOW_REPRESENTATIVE_CLUSTERING_UI && representativeClusteringComboBox != null) {
			sb.append("  Deep per-file clustering: ")
					.append(representativeClusteringComboBox.getValue()).append("\n");
			if (REP_CLUSTER_CUSTOM.equals(representativeClusteringComboBox.getValue())) {
				sb.append("  Clustering Jaccard threshold: ")
						.append(representativeJaccardField != null ? representativeJaccardField.getText() : "")
						.append("\n");
				sb.append("  Clustering containment threshold: ")
						.append(representativeContainmentField != null ? representativeContainmentField.getText() : "")
						.append("\n");
			}
			if (representativeDryRunCheckBox != null && representativeDryRunCheckBox.isVisible()) {
				sb.append("  Grouping dry run: ").append(representativeDryRunCheckBox.isSelected()).append("\n");
			}
		}
		sb.append("\n");

		// Genome selection
		sb.append("Genome Selection:\n");
		sb.append("  Selector: ").append(MagDbItem.MARKER4).append("\n");
		if (marker4StatusLabel != null && !marker4StatusLabel.getText().isEmpty()) {
			sb.append("  Marker4 Status: ").append(marker4StatusLabel.getText()).append("\n");
		}
		sb.append("\n");

		// Modifications
		sb.append("Modifications:\n");
		if (mExcisionCheckBox != null) {
			sb.append("  M excision: ").append(mExcisionCheckBox.isSelected()).append("\n");
		}
		if (carbamiCCheckBox != null) {
			sb.append("  Carbamidomethyl C: ").append(carbamiCCheckBox.isSelected()).append("\n");
		}
		if (oxidationMCheckBox != null) {
			sb.append("  Oxidation M: ").append(oxidationMCheckBox.isSelected()).append("\n");
		}
		if (acetyProNCheckBox != null) {
			sb.append("  Acetylation protein N-term: ").append(acetyProNCheckBox.isSelected()).append("\n");
		}
		if (phosphoSTYCheckBox != null) {
			sb.append("  Phosphorylation STY: ").append(phosphoSTYCheckBox.isSelected()).append("\n");
		}
		if (kGGCheckBox != null) {
			sb.append("  K GG: ").append(kGGCheckBox.isSelected()).append("\n");
		}
		sb.append("\n");

		consolePanel.setParameterText(sb.toString());
	}

	@Override
	protected MetaParameterDiaMarker collectParametersFromUI() {
		try {
			consolePanel.appendConsoleText("Collecting parameters from UI...\n");
			LOGGER.info("Starting parameter collection from UI for DIA Marker workflow");

			// Validate inputs
			if (!validateInputsInternal()) {
				consolePanel.appendConsoleText("ERROR: Input validation failed\n");
				LOGGER.error("Input validation failed");
				return null;
			}
			consolePanel.appendConsoleText("✓ Input validation passed\n");

			// Collect experiment data
			ObservableList<MetaPilotIOPanelFX.ExperimentRow> fileRows = ioPanel.getExperimentData();
			if (fileRows.isEmpty()) {
				LOGGER.error("No experiment data found");
				return null;
			}

			String[] rawFiles = new String[fileRows.size()];
			String[] expNames = new String[fileRows.size()];
			int[] fractions = new int[fileRows.size()];
			int[] replicates = new int[fileRows.size()];

			for (int i = 0; i < fileRows.size(); i++) {
				MetaPilotIOPanelFX.ExperimentRow row = fileRows.get(i);
				rawFiles[i] = row.getFilePath();
				expNames[i] = row.getExperimentName();
				try {
					fractions[i] = Integer.parseInt(row.getFraction());
					if (fractions[i] <= 0) {
						fractions[i] = 1;
					}
				} catch (NumberFormatException e) {
					fractions[i] = 1;
				}
				replicates[i] = 1;
			}

			// Create MetaData
			MetaData metaData = new MetaData(rawFiles, expNames, fractions, replicates);

			// Collect basic parameters
			String result = ioPanel.getOutputDirectory();
			MagDbItem selectedMagDb = magDatabasePanel.getSelectedMagDatabase();
			String microDb = selectedMagDb != null ? selectedMagDb.getCatalogueID() : "";
			boolean appendHostDb = magDatabasePanel.isAppendHostDb();
			String hostDB = appendHostDb ? magDatabasePanel.getHostDatabase() : "";

			int threadCount = threadsComboBox.getValue() != null ? threadsComboBox.getValue() : 4;

			// Create base parameter
			uk.ac.quadram.metapilot.task.par.MetaBaseParameter baseParams = 
				new uk.ac.quadram.metapilot.task.par.MetaBaseParameter(
			uk.ac.quadram.metapilot.workflow.WorkflowType.DIA_MARKER.getAnalysisType(),
					metaData, result, microDb, appendHostDb, hostDB, "DIA", 1, threadCount, true);

			// Create parameter object (MetaParameterDiaMarker for DIA Marker workflow)
			MetaParameterDiaMarker params = new MetaParameterDiaMarker(baseParams, selectedMagDb, true, false,
				new String[0], false, new String[0]);

			// Set MAG database
			if (selectedMagDb != null) {
				params.setMagDbItem(selectedMagDb);
				MagVersion selectedVersion = magDatabasePanel.getSelectedMagVersion();
				if (selectedVersion != null) {
					selectedMagDb.setUsedVersion(selectedVersion);
				}
			}

			// Set modification settings (with null checks)
			if (mExcisionCheckBox != null) {
				params.setMExcision(mExcisionCheckBox.isSelected());
			}
			if (carbamiCCheckBox != null) {
				params.setCarbamiC(carbamiCCheckBox.isSelected());
			}
			if (oxidationMCheckBox != null) {
				params.setOxidationM(oxidationMCheckBox.isSelected());
			}
			if (acetyProNCheckBox != null) {
				params.setAcetyProN(acetyProNCheckBox.isSelected());
			}
			if (phosphoSTYCheckBox != null) {
				params.setPhosphoSTY(phosphoSTYCheckBox.isSelected());
			}
			if (kGGCheckBox != null) {
				params.setKGG(kGGCheckBox.isSelected());
			}

			// Set genome selection strategy based on workflow
			String workflow = workflowModeComboBox != null && workflowModeComboBox.getValue() != null
					? workflowModeComboBox.getValue()
					: WORKFLOW_FOCUS;
			params.setMarkerRefinementMode(DiaRefinementMode.fromDisplayLabel(workflow, DiaRefinementMode.CORE));
			applyRepresentativeClusteringFromUI(params);
			applyEntrapmentFromUI(params);

			if (metadataPanel != null) {
				metadataPanel.applyTo(metaData);
			}

			consolePanel.appendConsoleText("✓ Parameters collected successfully\n");
			LOGGER.info("Parameters collected successfully");
			return params;

		} catch (Exception e) {
			LOGGER.error("Failed to collect parameters", e);
			consolePanel.appendConsoleText("ERROR: " + e.getMessage() + "\n");
			return null;
		}
	}

	public MetaParameterDiaMarker collectParametersFromUIForAutoSave() {
		// Lenient version for auto-save - skip validation
		try {
			// Create minimal MetaData even if no experiment data
			MetaData metaData = new MetaData(new String[0], new String[0], new int[0], new int[0]);
			
			if (ioPanel != null) {
				ObservableList<MetaPilotIOPanelFX.ExperimentRow> fileRows = ioPanel.getExperimentData();
				if (!fileRows.isEmpty()) {
					String[] rawFiles = new String[fileRows.size()];
					String[] expNames = new String[fileRows.size()];
					int[] fractions = new int[fileRows.size()];
					int[] replicates = new int[fileRows.size()];

					for (int i = 0; i < fileRows.size(); i++) {
						MetaPilotIOPanelFX.ExperimentRow row = fileRows.get(i);
						rawFiles[i] = row.getFilePath();
						expNames[i] = row.getExperimentName();
						try {
							fractions[i] = Integer.parseInt(row.getFraction());
							if (fractions[i] <= 0) {
								fractions[i] = 1;
							}
						} catch (NumberFormatException e) {
							fractions[i] = 1;
						}
						replicates[i] = 1;
					}
					metaData = new MetaData(rawFiles, expNames, fractions, replicates);
				}
			}

			String result = ioPanel != null ? ioPanel.getOutputDirectory() : "";
			MagDbItem selectedMagDb = magDatabasePanel != null ? magDatabasePanel.getSelectedMagDatabase() : null;
			String microDb = selectedMagDb != null ? selectedMagDb.getCatalogueID() : "";
			boolean appendHostDb = magDatabasePanel != null && magDatabasePanel.isAppendHostDb();
			String hostDB = appendHostDb && magDatabasePanel != null ? magDatabasePanel.getHostDatabase() : "";
			int threadCount = threadsComboBox != null && threadsComboBox.getValue() != null ? threadsComboBox.getValue() : 4;

			uk.ac.quadram.metapilot.task.par.MetaBaseParameter baseParams = 
				new uk.ac.quadram.metapilot.task.par.MetaBaseParameter(
			uk.ac.quadram.metapilot.workflow.WorkflowType.DIA_MARKER.getAnalysisType(),
					metaData, result, microDb, appendHostDb, hostDB, "DIA", 1, threadCount, true);

			MetaParameterDiaMarker params = new MetaParameterDiaMarker(baseParams, selectedMagDb, true, false,
				new String[0], false, new String[0]);

			// Set MAG database
			if (selectedMagDb != null) {
				params.setMagDbItem(selectedMagDb);
				if (magDatabasePanel != null) {
					MagVersion selectedVersion = magDatabasePanel.getSelectedMagVersion();
					if (selectedVersion != null) {
						selectedMagDb.setUsedVersion(selectedVersion);
					}
				}
			}

			// Set modification settings (with null checks)
			if (mExcisionCheckBox != null) {
				params.setMExcision(mExcisionCheckBox.isSelected());
			}
			if (carbamiCCheckBox != null) {
				params.setCarbamiC(carbamiCCheckBox.isSelected());
			}
			if (oxidationMCheckBox != null) {
				params.setOxidationM(oxidationMCheckBox.isSelected());
			}
			if (acetyProNCheckBox != null) {
				params.setAcetyProN(acetyProNCheckBox.isSelected());
			}
			if (phosphoSTYCheckBox != null) {
				params.setPhosphoSTY(phosphoSTYCheckBox.isSelected());
			}
			if (kGGCheckBox != null) {
				params.setKGG(kGGCheckBox.isSelected());
			}

			// Set genome selection strategy based on workflow (best effort for auto-save)
			String workflow = workflowModeComboBox != null && workflowModeComboBox.getValue() != null
					? workflowModeComboBox.getValue()
					: WORKFLOW_FOCUS;
			params.setMarkerRefinementMode(DiaRefinementMode.fromDisplayLabel(workflow, DiaRefinementMode.CORE));
			applyRepresentativeClusteringFromUI(params);
			applyEntrapmentFromUI(params);

			return params;
		} catch (Exception e) {
			LOGGER.warn("Failed to collect parameters for auto-save", e);
			return null;
		}
	}

	private boolean validateInputsInternal() {
		if (ioPanel == null || magDatabasePanel == null) {
			LOGGER.error("UI components not initialized");
			return false;
		}

		String outputDir = ioPanel.getOutputDirectory();
		if (outputDir == null || outputDir.trim().isEmpty()) {
			showError("Validation Error", "Please select an output directory.");
			return false;
		}

		if (!magDatabasePanel.validateSelection()) {
			String errorMsg = magDatabasePanel.getValidationErrorMessage();
			LOGGER.error("MAG database validation failed: {}", errorMsg);
			showError("Validation Error", errorMsg);
			return false;
		}

		MagDbItem selectedMagDb = magDatabasePanel.getSelectedMagDatabase();
		if (selectedMagDb == null) {
			showError("Validation Error", "Please select a MAG database");
			return false;
		}

		if (SHOW_REPRESENTATIVE_CLUSTERING_UI
				&& workflowModeComboBox != null && WORKFLOW_DEEP_PER_FILE.equals(workflowModeComboBox.getValue())
				&& representativeClusteringComboBox != null
				&& REP_CLUSTER_CUSTOM.equals(representativeClusteringComboBox.getValue())) {
			try {
				parseThreshold(representativeJaccardField, "representative peptide Jaccard threshold");
				parseThreshold(representativeContainmentField, "representative peptide containment threshold");
			} catch (IllegalArgumentException e) {
				showError("Validation Error", e.getMessage());
				return false;
			}
		}

		// Marker4 files are generated automatically by the task when missing
		return true;
	}

	@Override
	protected void loadLastUsedParameters() {
		try {
			MetaParameterDiaMarker params = (MetaParameterDiaMarker) parameterIO.loadLastUsedParameters();
			if (params != null) {
				applyParametersToUI(params);
			}
		} catch (Exception e) {
			LOGGER.warn("Failed to load last used parameters", e);
		}
	}

	@Override
	protected void applyParametersToUI(MetaParameterDiaMarker params) {
		if (params == null) {
			return;
		}

		// Apply metadata to IO panel
		if (params.getMetadata() != null && ioPanel != null) {
			MetaData metadata = params.getMetadata();
			String[] rawFiles = metadata.getRawFiles();
			String[] expNames = metadata.getExpNames();
			int[] fractions = metadata.getFractions();

			if (rawFiles != null && expNames != null && rawFiles.length > 0) {
				ObservableList<MetaPilotIOPanelFX.ExperimentRow> experimentRows = FXCollections.observableArrayList();
				for (int i = 0; i < rawFiles.length; i++) {
					if (rawFiles[i] != null && !rawFiles[i].trim().isEmpty()) {
						MetaPilotIOPanelFX.ExperimentRow row = new MetaPilotIOPanelFX.ExperimentRow(rawFiles[i], expNames[i]);
						if (fractions != null && i < fractions.length) {
							row.setFraction(String.valueOf(fractions[i]));
						}
						experimentRows.add(row);
					}
				}
				ioPanel.setExperimentData(experimentRows);
			}
		}

		// Apply output directory
		if (params.getResult() != null && ioPanel != null) {
			ioPanel.setOutputDirectory(params.getResult());
		}

		// Apply thread count
		if (threadsComboBox != null) {
			threadsComboBox.setValue(params.getThreadCount());
		}

		// Apply MAG database
		if (magDatabasePanel != null && params.getMagDbItem() != null) {
			MagDbItem savedMagDb = params.getMagDbItem();
			ObservableList<MagDbItem> availableDbs = magDatabasePanel.getAvailableMagDatabases();
			if (availableDbs != null) {
				for (MagDbItem item : availableDbs) {
					if (item.getCatalogueID().equals(savedMagDb.getCatalogueID())) {
						magDatabasePanel.setSelectedMagDatabase(item);
						if (savedMagDb.getUsedVersion() != null) {
							for (MagVersion version : item.getVersionList()) {
								if (version.getCatalogueVersion().equals(savedMagDb.getUsedVersion().getCatalogueVersion())) {
									magDatabasePanel.setSelectedMagVersion(version);
									break;
								}
							}
						}
						break;
					}
				}
			}
		}
		if (magDatabasePanel != null) {
			magDatabasePanel.setAppendHostDb(params.isAppendHostDb());
			magDatabasePanel.setHostDatabase(params.getHostDB());
		}

		// Apply modification settings
		if (mExcisionCheckBox != null) {
			mExcisionCheckBox.setSelected(params.isMExcision());
		}
		if (carbamiCCheckBox != null) {
			carbamiCCheckBox.setSelected(params.isCarbamiC());
		}
		if (oxidationMCheckBox != null) {
			oxidationMCheckBox.setSelected(params.isOxidationM());
		}
		if (acetyProNCheckBox != null) {
			acetyProNCheckBox.setSelected(params.isAcetyProN());
		}
		if (phosphoSTYCheckBox != null) {
			phosphoSTYCheckBox.setSelected(params.isPhosphoSTY());
		}
		if (kGGCheckBox != null) {
			kGGCheckBox.setSelected(params.isKGG());
		}

		// Map stored genome selection strategy back to workflow combo (best effort)
		if (workflowModeComboBox != null) {
			DiaRefinementMode markerMode = params.getMarkerRefinementMode();
			if (markerMode == DiaRefinementMode.DEEP_TOGETHER) {
				workflowModeComboBox.setValue(WORKFLOW_DEEP);
			} else if (markerMode == DiaRefinementMode.DEEP_PER_FILE) {
				workflowModeComboBox.setValue(WORKFLOW_DEEP_PER_FILE);
			} else if (markerMode == DiaRefinementMode.COMPREHENSIVE) {
				workflowModeComboBox.setValue(WORKFLOW_COMPREHENSIVE);
			} else {
				// Default and legacy fallback: treat anything else as FOCUS
				workflowModeComboBox.setValue(WORKFLOW_FOCUS);
			}
		}
		if (representativeClusteringComboBox != null) {
			representativeClusteringComboBox.setValue(representativeClusteringPresetToDisplay(
					params.getRepresentativeClusteringPreset()));
		}
		if (representativeJaccardField != null) {
			representativeJaccardField.setText(formatThreshold(params.getRepresentativePeptideJaccardThreshold()));
		}
		if (representativeContainmentField != null) {
			representativeContainmentField.setText(formatThreshold(params.getRepresentativePeptideContainmentThreshold()));
		}
		if (representativeDryRunCheckBox != null) {
			representativeDryRunCheckBox.setSelected(params.isRepresentativeClusteringDryRun());
		}
		updateRepresentativeClusteringAvailability();
		if (entrapmentCheckBox != null) {
			entrapmentCheckBox.setSelected(params.isEntrapmentEnabled());
		}
		updateEntrapmentAvailability();
		if (entrapmentCheckBox != null && entrapmentCheckBox.isVisible()) {
			entrapmentCheckBox.setSelected(params.isEntrapmentEnabled());
		}
	}

	private void applyRepresentativeClusteringFromUI(MetaParameterDiaMarker params) {
		if (params == null) {
			return;
		}
		boolean deepPerFile = workflowModeComboBox != null
				&& WORKFLOW_DEEP_PER_FILE.equals(workflowModeComboBox.getValue());
		if (!SHOW_REPRESENTATIVE_CLUSTERING_UI || !deepPerFile || representativeClusteringComboBox == null
				|| representativeClusteringComboBox.getValue() == null) {
			params.setRepresentativeClusteringPreset(MetaParameterDiaMarker.REPRESENTATIVE_CLUSTERING_OFF);
			params.setRepresentativeClusteringDryRun(false);
			return;
		}
		String preset = representativeClusteringDisplayToPreset(representativeClusteringComboBox.getValue());
		params.setRepresentativeClusteringPreset(preset);
		if (MetaParameterDiaMarker.REPRESENTATIVE_CLUSTERING_CUSTOM.equals(preset)) {
			params.setRepresentativePeptideJaccardThreshold(parseThreshold(representativeJaccardField,
					"representative peptide Jaccard threshold"));
			params.setRepresentativePeptideContainmentThreshold(parseThreshold(representativeContainmentField,
					"representative peptide containment threshold"));
		}
		params.setRepresentativeClusteringDryRun(representativeDryRunCheckBox != null
				&& representativeDryRunCheckBox.isVisible()
				&& representativeDryRunCheckBox.isSelected());
	}

	private void applyEntrapmentFromUI(MetaParameterDiaMarker params) {
		if (params == null) {
			return;
		}
		boolean deepTogether = workflowModeComboBox != null
				&& WORKFLOW_DEEP.equals(workflowModeComboBox.getValue());
		params.setEntrapmentEnabled(deepTogether && entrapmentCheckBox != null && entrapmentCheckBox.isSelected());
	}

	private String representativeClusteringDisplayToPreset(String value) {
		if (REP_CLUSTER_CONSERVATIVE.equals(value)) {
			return MetaParameterDiaMarker.REPRESENTATIVE_CLUSTERING_CONSERVATIVE;
		}
		if (REP_CLUSTER_BALANCED.equals(value)) {
			return MetaParameterDiaMarker.REPRESENTATIVE_CLUSTERING_BALANCED;
		}
		if (REP_CLUSTER_AGGRESSIVE.equals(value)) {
			return MetaParameterDiaMarker.REPRESENTATIVE_CLUSTERING_AGGRESSIVE;
		}
		if (REP_CLUSTER_CUSTOM.equals(value)) {
			return MetaParameterDiaMarker.REPRESENTATIVE_CLUSTERING_CUSTOM;
		}
		return MetaParameterDiaMarker.REPRESENTATIVE_CLUSTERING_OFF;
	}

	private String representativeClusteringPresetToDisplay(String preset) {
		if (MetaParameterDiaMarker.REPRESENTATIVE_CLUSTERING_CONSERVATIVE.equalsIgnoreCase(preset)) {
			return REP_CLUSTER_CONSERVATIVE;
		}
		if (MetaParameterDiaMarker.REPRESENTATIVE_CLUSTERING_BALANCED.equalsIgnoreCase(preset)) {
			return REP_CLUSTER_BALANCED;
		}
		if (MetaParameterDiaMarker.REPRESENTATIVE_CLUSTERING_AGGRESSIVE.equalsIgnoreCase(preset)) {
			return REP_CLUSTER_AGGRESSIVE;
		}
		if (MetaParameterDiaMarker.REPRESENTATIVE_CLUSTERING_CUSTOM.equalsIgnoreCase(preset)) {
			return REP_CLUSTER_CUSTOM;
		}
		return REP_CLUSTER_OFF;
	}

	private double parseThreshold(TextField field, String label) {
		if (field == null || field.getText() == null || field.getText().trim().isEmpty()) {
			throw new IllegalArgumentException("Missing " + label);
		}
		double value = Double.parseDouble(field.getText().trim());
		if (value <= 0.0 || value > 1.0) {
			throw new IllegalArgumentException(label + " must be > 0 and <= 1");
		}
		return value;
	}

	private String formatThreshold(double value) {
		return String.format(java.util.Locale.ROOT, "%.2f", value);
	}

	@Override
	protected Class<MetaParameterDiaMarker> getParameterClass() {
		return MetaParameterDiaMarker.class;
	}

	public void saveCurrentParameters() {
		try {
			MetaParameterDiaMarker currentParams = collectParametersFromUI();
			if (currentParams == null) {
				currentParams = collectParametersFromUIForAutoSave();
			}
			if (currentParams != null) {
				parameterIO.saveLastUsedParameters(currentParams);
			}
		} catch (Exception e) {
			LOGGER.warn("Failed to save current parameters", e);
		}
	}

	@Override
	protected void startAnalysis() {
		try {
			consolePanel.appendConsoleText("Starting DIA Marker analysis...\n");
			LOGGER.info("Starting DIA Marker analysis");

			MetaParameterDiaMarker analysisParams = collectParametersFromUI();
			if (analysisParams == null) {
				consolePanel.appendConsoleText("ERROR: Failed to collect parameters\n");
				return;
			}

			this.currentAnalysisParams = analysisParams;

			// Save parameters
			String outputDir = analysisParams.getResult();
			if (outputDir != null && !outputDir.trim().isEmpty()) {
				parameterIO.saveAnalysisParameters(outputDir, analysisParams);
				parameterIO.saveLastUsedParameters(analysisParams);
			}

			// Load resource parameters
			MetaResourceIO resourceIO = new MetaResourceIO();
			MetaResourceParameters mrp = resourceIO.loadParameters();

			// Create task from canonical mode. Persisted legacy values are normalized in the parameter class.
			String workflow = workflowModeComboBox != null && workflowModeComboBox.getValue() != null
					? workflowModeComboBox.getValue()
					: WORKFLOW_FOCUS;
			DiaRefinementMode markerMode = DiaRefinementMode.fromDisplayLabel(workflow, DiaRefinementMode.CORE);
			analysisParams.setMarkerRefinementMode(markerMode);
			if (markerMode == DiaRefinementMode.DEEP_TOGETHER) {
				if (analysisParams.isEntrapmentEnabled()) {
					File entrapmentDir = new File(analysisParams.getResult(),
							DiaMarkerDeepEntrapmentHeadlessRunner.RESULT_FOLDER);
					entrapmentDir.mkdirs();
					analysisParams.setDbSearchResultFile(entrapmentDir);
					consolePanel.appendConsoleText("Evaluation: paired-shuffle entrapment -> "
							+ entrapmentDir.getAbsolutePath() + "\n");
				}
				currentTask = new DiaMarkerDeepTogetherTask(analysisParams, mrp,
						consolePanel.getOverallProgressBar(), consolePanel.getCurrentTaskProgressBar());
				consolePanel.appendConsoleText("Workflow: " + markerMode.getDisplayLabel() + "\n");
			} else if (markerMode == DiaRefinementMode.DEEP_PER_FILE) {
				currentTask = new DiaMarkerDeepPerFileTask(analysisParams, mrp,
						consolePanel.getOverallProgressBar(), consolePanel.getCurrentTaskProgressBar());
				consolePanel.appendConsoleText("Workflow: " + markerMode.getDisplayLabel() + "\n");
			} else if (markerMode == DiaRefinementMode.COMPREHENSIVE) {
				currentTask = new DiaMarkerComprehensiveTask(analysisParams, mrp,
						consolePanel.getOverallProgressBar(), consolePanel.getCurrentTaskProgressBar());
				consolePanel.appendConsoleText("Workflow: " + markerMode.getDisplayLabel() + "\n");
			} else {
				currentTask = new DiaMarkerCoreTask(analysisParams, mrp,
						consolePanel.getOverallProgressBar(), consolePanel.getCurrentTaskProgressBar());
				consolePanel.appendConsoleText("Workflow: " + markerMode.getDisplayLabel() + "\n");
			}

			// Set up task callbacks
			currentTask.messageProperty().addListener((obs, oldMsg, newMsg) -> {
				Platform.runLater(() -> {
					if (newMsg != null) {
						consolePanel.appendConsoleText("Current Task: " + newMsg + "\n");
						consolePanel.updateCurrentTaskFromMessage(newMsg);
					}
				});
			});

			currentTask.setOnSucceeded(e -> {
				Platform.runLater(() -> {
					try {
						Boolean result = currentTask.get();
						if (result != null && result) {
							consolePanel.appendConsoleText("✓ Analysis completed successfully!\n");
							consolePanel.enableViewReportButton();
						} else {
							consolePanel.appendConsoleText("⚠ Analysis completed but may not have been fully successful.\n");
						}
					} catch (Exception ex) {
						consolePanel.appendConsoleText("✗ Analysis failed: " + ex.getMessage() + "\n");
					}
					consolePanel.setButtonsEnabled(true, false, true);
					consolePanel.stopStatusTimer();
				});
			});

			currentTask.setOnFailed(e -> {
				Platform.runLater(() -> {
					consolePanel.appendConsoleText("✗ Analysis failed. Check console output for details.\n");
					consolePanel.setButtonsEnabled(true, false, false);
					consolePanel.stopStatusTimer();
				});
			});

			// Run task
			Thread taskThread = new Thread(currentTask);
			taskThread.setDaemon(true);
			taskThread.start();

			consolePanel.appendConsoleText("✓ Analysis task started successfully.\n");
			LOGGER.info("DIA Marker analysis task started");

		} catch (Exception e) {
			LOGGER.error("Failed to start DIA Marker analysis", e);
			consolePanel.appendConsoleText("ERROR: " + e.getMessage() + "\n");
		}
	}

	@Override
	protected void stopAnalysis() {
		Alert alert = new Alert(AlertType.CONFIRMATION);
		alert.setTitle("Stop Analysis");
		alert.setHeaderText("Are you sure you want to stop the analysis?");
		alert.setContentText("This action cannot be undone.");

		Optional<ButtonType> result = alert.showAndWait();
		if (result.isPresent() && result.get() == ButtonType.OK) {
			consolePanel.appendConsoleText("Stopping DIA Marker analysis...\n");
			if (currentTask != null && currentTask.isRunning()) {
				currentTask.forceStop();
			}
			consolePanel.setButtonsEnabled(true, false, false);
			consolePanel.stopStatusTimer();
		}
	}

	// viewReport() uses base implementation: currentAnalysisParams.getReportDir()
	// (set by the task to the workflow folder, e.g. result/marker_focus/report)

	@Override
	public void loadParametersFromObject(Object parameters) {
		if (!(parameters instanceof MetaParameterDiaMarker)) {
			LOGGER.warn("Invalid parameter type for DIA Marker: {}", 
				parameters != null ? parameters.getClass().getSimpleName() : "null");
			return;
		}

		MetaParameterDiaMarker params = (MetaParameterDiaMarker) parameters;

		// Load basic parameters - use applyParametersToUI which handles this properly
		applyParametersToUI(params);

		if (params.getResult() != null && ioPanel != null) {
			ioPanel.setOutputDirectory(params.getResult());
		}

		// Load MAG database
		if (magDatabasePanel != null && params.getMagDbItem() != null) {
			MagDbItem savedMagDb = params.getMagDbItem();
			ObservableList<MagDbItem> availableDbs = magDatabasePanel.getAvailableMagDatabases();
			if (availableDbs != null) {
				for (MagDbItem item : availableDbs) {
					if (item.getCatalogueID().equals(savedMagDb.getCatalogueID())) {
						magDatabasePanel.setSelectedMagDatabase(item);
						if (savedMagDb.getUsedVersion() != null) {
							for (MagVersion version : item.getVersionList()) {
								if (version.getCatalogueVersion().equals(savedMagDb.getUsedVersion().getCatalogueVersion())) {
									magDatabasePanel.setSelectedMagVersion(version);
									break;
								}
							}
						}
						break;
					}
				}
			}
		}
		if (magDatabasePanel != null) {
			magDatabasePanel.setAppendHostDb(params.isAppendHostDb());
			magDatabasePanel.setHostDatabase(params.getHostDB());
		}

		// Load thread count
		if (threadsComboBox != null) {
			int threadCount = params.getThreadCount();
			if (threadCount > 0) {
				threadsComboBox.setValue(threadCount);
			}
		}

		// Modifications are stored in MetaParameterDiaMarker
		updateParameterView();
		updateMarker4Status();
	}
}
