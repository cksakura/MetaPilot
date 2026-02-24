package uk.ac.quadram.metapilot.task.dia.gui;

import javafx.application.Platform;
import javafx.geometry.Insets;
import javafx.scene.control.*;
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
import uk.ac.quadram.metapilot.task.mag.MarkerGenomeSelector;
import uk.ac.quadram.metapilot.task.dia.par.DiaFocusParameterIO;
import uk.ac.quadram.metapilot.task.dia.par.MetaParameterDiaFocus;
import uk.ac.quadram.metapilot.task.dia.MetaDiaMarkerTask;
import uk.ac.quadram.metapilot.task.dia.MetaDiaDeepMarkerTask;
import uk.ac.quadram.metapilot.task.mag.MagDbItem.MagVersion;
import uk.ac.quadram.metapilot.task.dda.igc.gui.MetaPilotDdaIgcMaxQuantPanel;
import uk.ac.quadram.metapilot.resources.MetaResourceIO;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.File;
import java.util.Optional;

/**
 * DIA Marker workflow panel: marker-based genome selection (Focus or Comprehensive)
 * using Marker4 and MetaDiaMarkerTask.
 */
public class MetaPilotDiaMarkerPanel extends AbstractWorkflowPanel<MetaParameterDiaFocus> {
	private static final Logger LOGGER = LogManager.getLogger(MetaPilotDiaMarkerPanel.class);

	// UI Components
	private MetaPilotIOPanelFX ioPanel;
	private MagDatabaseSettingsPanel magDatabasePanel;
	private MetadataConfigurationPanel metadataPanel;
	// Note: consolePanel is inherited from AbstractWorkflowPanel

	/** Workflow mode: "Marker" (standard) or "Deep Marker" (iterative rounds + combined search). */
	private static final String WORKFLOW_MARKER = "Marker";
	private static final String WORKFLOW_DEEP_MARKER = "Deep Marker";

	// Parameter components
	private ComboBox<String> workflowModeComboBox;
	private ComboBox<Integer> threadsComboBox;
	private ComboBox<String> genomeSelectionStrategyComboBox;
	private Label marker4StatusLabel;

	// Modification checkboxes (from DiannParameter lines 28-33)
	private CheckBox mExcisionCheckBox;
	private CheckBox carbamiCCheckBox;
	private CheckBox oxidationMCheckBox;
	private CheckBox acetyProNCheckBox;
	private CheckBox phosphoSTYCheckBox;
	private CheckBox kGGCheckBox;

	// Task management
	private MetaDiaMarkerTask currentTask;
	private MetaResourceParameters resourceParams;
	private MetaParameterDiaFocus currentAnalysisParams;


	public MetaPilotDiaMarkerPanel(MagDbItem[] magDbItems, MetaResourceParameters resourceParams) {
		super();
		parameterIO = new DiaFocusParameterIO();
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
		metadataPanel = new MetadataConfigurationPanel();
		magDatabasePanel = new MagDatabaseSettingsPanel(false); // No host database for DIA Marker

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
				// Auto-save parameters
				try {
					MetaParameterDiaFocus currentParams = collectParametersFromUI();
					if (currentParams == null) {
						currentParams = collectParametersFromUIForAutoSave();
					}
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
		tabContent.setMinHeight(760);

		// Create main grid layout
		GridPane mainGrid = new GridPane();
		mainGrid.setHgap(15);
		mainGrid.setVgap(15);
		mainGrid.setPadding(new Insets(10));
		mainGrid.setMinWidth(940);
		mainGrid.setMinHeight(740);
		VBox.setVgrow(mainGrid, Priority.ALWAYS);

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
		row1.setMinHeight(400);
		mainGrid.getRowConstraints().add(row1);

		// Left column: Performance Settings, Genome Selection, Modifications
		VBox leftColumn = new VBox(15);
		leftColumn.setMinWidth(430);
		VBox workflowSection = createWorkflowSection();
		VBox performanceSection = createPerformanceSection();
		VBox genomeSelectionSection = createGenomeSelectionSection();
		VBox modificationsSection = createModificationsSection();
		leftColumn.getChildren().addAll(workflowSection, performanceSection, genomeSelectionSection, modificationsSection);
		GridPane.setVgrow(leftColumn, Priority.ALWAYS);

		// Right column: Metadata Configuration
		VBox rightColumn = new VBox(15);
		rightColumn.setMinWidth(430);
		rightColumn.getChildren().addAll(metadataPanel);
		GridPane.setVgrow(rightColumn, Priority.ALWAYS);

		mainGrid.add(leftColumn, 0, 0);
		mainGrid.add(rightColumn, 1, 0);

		tabContent.getChildren().add(mainGrid);
		return tabContent;
	}

	private VBox createWorkflowSection() {
		VBox section = new VBox(10);
		section.setPadding(new Insets(10));
		section.setStyle("-fx-background-color: #f0f8ff; -fx-background-radius: 5;");

		Label titleLabel = createSectionTitle("Workflow");
		Label workflowLabel = new Label("Mode:");
		workflowLabel.setStyle("-fx-font-weight: bold;");
		workflowModeComboBox = new ComboBox<>(
				FXCollections.observableArrayList(WORKFLOW_MARKER, WORKFLOW_DEEP_MARKER));
		workflowModeComboBox.setMinWidth(160);
		workflowModeComboBox.setPrefWidth(200);
		workflowModeComboBox.setValue(WORKFLOW_MARKER);
		workflowModeComboBox.setOnAction(e -> updateParameterView());
		Tooltip workflowTooltip = new Tooltip(
				"Marker: one combined library from selected genomes, then identification.\n"
				+ "Deep Marker: iterative rounds (round1.parquet, round2.parquet, ...) with round libraries (round1-lib.parquet, ...), then combined search and report.");
		workflowModeComboBox.setTooltip(workflowTooltip);
		VBox workflowBox = new VBox(5);
		workflowBox.getChildren().addAll(workflowLabel, workflowModeComboBox);
		section.getChildren().addAll(titleLabel, workflowBox);
		return section;
	}

	private VBox createPerformanceSection() {
		VBox section = new VBox(10);
		section.setPadding(new Insets(10));
		section.setStyle("-fx-background-color: #f0f8ff; -fx-background-radius: 5;");

		Label titleLabel = createSectionTitle("Performance Settings");

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

		section.getChildren().addAll(titleLabel, threadsBox);
		return section;
	}

	private VBox createGenomeSelectionSection() {
		VBox section = new VBox(10);
		section.setPadding(new Insets(10));
		section.setStyle("-fx-background-color: #f0f8ff; -fx-background-radius: 5;");

		Label titleLabel = createSectionTitle("Genome Selection");

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

		Label genomeSelectionStrategyLabel = new Label("Selection Strategy:");
		genomeSelectionStrategyLabel.setStyle("-fx-font-weight: bold;");
		genomeSelectionStrategyComboBox = new ComboBox<>(
				FXCollections.observableArrayList(MarkerGenomeSelector.KNEE_MODE_FOCUS,
						MarkerGenomeSelector.KNEE_MODE_COMPREHENSIVE));
		genomeSelectionStrategyComboBox.setMinWidth(140);
		genomeSelectionStrategyComboBox.setPrefWidth(160);
		genomeSelectionStrategyComboBox.setMaxWidth(Double.MAX_VALUE);
		genomeSelectionStrategyComboBox.setValue(MarkerGenomeSelector.KNEE_MODE_FOCUS);
		genomeSelectionStrategyComboBox.setOnAction(e -> updateParameterView());
		VBox genomeSelectionStrategyBox = new VBox(5);
		genomeSelectionStrategyBox.getChildren().addAll(genomeSelectionStrategyLabel, genomeSelectionStrategyComboBox);

		Label genomeSelectorLabel = new Label("Genome Selection:");
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

		paramGrid.add(genomeSelectionStrategyBox, 0, 0);
		paramGrid.add(genomeSelectorBox, 1, 0);
		paramGrid.add(marker4StatusBox, 0, 1);
		paramGrid.add(new Region(), 1, 1); // placeholder for alignment

		section.getChildren().addAll(titleLabel, paramGrid);
		return section;
	}

	private VBox createModificationsSection() {
		VBox section = new VBox(10);
		section.setPadding(new Insets(10));
		section.setStyle("-fx-background-color: #f0f8ff; -fx-background-radius: 5;");

		Label titleLabel = createSectionTitle("Modifications");

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

		section.getChildren().addAll(titleLabel, checkboxesBox);
		return section;
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

	@Override
	protected void updateParameterView() {
		if (consolePanel == null) {
			return;
		}
		StringBuilder sb = new StringBuilder();
		sb.append("=== Analysis Parameters ===\n\n");

		if (workflowModeComboBox != null && workflowModeComboBox.getValue() != null) {
			sb.append("Workflow: ").append(workflowModeComboBox.getValue()).append("\n\n");
		}

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
		}
		sb.append("\n");

		// Performance settings
		sb.append("Performance Settings:\n");
		if (threadsComboBox != null) {
			sb.append("  Threads: ").append(threadsComboBox.getValue()).append("\n");
		}
		sb.append("\n");

		// Genome selection
		sb.append("Genome Selection:\n");
		if (genomeSelectionStrategyComboBox != null) {
			sb.append("  Selection Strategy: ").append(genomeSelectionStrategyComboBox.getValue()).append("\n");
		}
		sb.append("  Mode: ").append(MagDbItem.MARKER4).append("\n");
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
	protected MetaParameterDiaFocus collectParametersFromUI() {
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

			int threadCount = threadsComboBox.getValue() != null ? threadsComboBox.getValue() : 4;

			// Create base parameter
			uk.ac.quadram.metapilot.task.par.MetaBaseParameter baseParams = 
				new uk.ac.quadram.metapilot.task.par.MetaBaseParameter(
			uk.ac.quadram.metapilot.workflow.WorkflowType.DIA_MARKER.getAnalysisType(),
					metaData, result, microDb, false, "", "DIA", 1, threadCount, true);

			// Create MetaParameterDiaFocus
			MetaParameterDiaFocus params = new MetaParameterDiaFocus(baseParams, selectedMagDb, true, false, 
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
			if (genomeSelectionStrategyComboBox != null && genomeSelectionStrategyComboBox.getValue() != null) {
				params.setGenomeSelectionStrategy(genomeSelectionStrategyComboBox.getValue());
			}

			// Collect metadata
			if (metadataPanel != null) {
				ObservableList<MetaPilotDdaIgcMaxQuantPanel.MetadataRow> metadataRows = metadataPanel.getMetadataTableData();
				java.util.List<String> colNames = metadataPanel.getMetadataColumnNames();
				int metaTypeCount = colNames.size();
				String[][] metaInfo = new String[metadataRows.size()][Math.max(1, metaTypeCount)];
				if (!colNames.isEmpty()) {
					for (int i = 0; i < metadataRows.size(); i++) {
						MetaPilotDdaIgcMaxQuantPanel.MetadataRow row = metadataRows.get(i);
						for (int j = 0; j < colNames.size(); j++) {
							String v = row.getValueByColumnName(colNames.get(j));
							metaInfo[i][j] = v != null ? v : "";
						}
					}
				}
				metaData.setMetaTypeCount(Math.max(1, metaTypeCount));
				metaData.setMetaInfo(metaInfo);
				if (!colNames.isEmpty()) {
					metaData.setMetaTypeTitles(colNames.toArray(new String[0]));
				}
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

	public MetaParameterDiaFocus collectParametersFromUIForAutoSave() {
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
			int threadCount = threadsComboBox != null && threadsComboBox.getValue() != null ? threadsComboBox.getValue() : 4;

			uk.ac.quadram.metapilot.task.par.MetaBaseParameter baseParams = 
				new uk.ac.quadram.metapilot.task.par.MetaBaseParameter(
			uk.ac.quadram.metapilot.workflow.WorkflowType.DIA_MARKER.getAnalysisType(),
					metaData, result, microDb, false, "", "DIA", 1, threadCount, true);

			MetaParameterDiaFocus params = new MetaParameterDiaFocus(baseParams, selectedMagDb, true, false, 
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
			if (genomeSelectionStrategyComboBox != null && genomeSelectionStrategyComboBox.getValue() != null) {
				params.setGenomeSelectionStrategy(genomeSelectionStrategyComboBox.getValue());
			}

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

		// Marker4 files are generated automatically by the task when missing
		return true;
	}

	@Override
	protected void loadLastUsedParameters() {
		try {
			MetaParameterDiaFocus params = (MetaParameterDiaFocus) parameterIO.loadLastUsedParameters();
			if (params != null) {
				applyParametersToUI(params);
			}
		} catch (Exception e) {
			LOGGER.warn("Failed to load last used parameters", e);
		}
	}

	@Override
	protected void applyParametersToUI(MetaParameterDiaFocus params) {
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
		if (genomeSelectionStrategyComboBox != null && params.getGenomeSelectionStrategy() != null) {
			String strategy = params.getGenomeSelectionStrategy();
			if (MarkerGenomeSelector.KNEE_MODE_FOCUS.equals(strategy)
					|| MarkerGenomeSelector.KNEE_MODE_COMPREHENSIVE.equals(strategy)) {
				genomeSelectionStrategyComboBox.setValue(strategy);
			} else {
				genomeSelectionStrategyComboBox.setValue(MarkerGenomeSelector.KNEE_MODE_FOCUS);
			}
		}
	}

	@Override
	protected Class<MetaParameterDiaFocus> getParameterClass() {
		return MetaParameterDiaFocus.class;
	}

	public void saveCurrentParameters() {
		try {
			MetaParameterDiaFocus currentParams = collectParametersFromUI();
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

			MetaParameterDiaFocus analysisParams = collectParametersFromUI();
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

			// Create task - choose Marker (standard) or Deep Marker (iterative rounds + combined)
			boolean deepMarker = WORKFLOW_DEEP_MARKER.equals(workflowModeComboBox != null ? workflowModeComboBox.getValue() : WORKFLOW_MARKER);
			if (deepMarker) {
				currentTask = new MetaDiaDeepMarkerTask(analysisParams, mrp,
						consolePanel.getOverallProgressBar(), consolePanel.getCurrentTaskProgressBar());
				consolePanel.appendConsoleText("Workflow: Deep Marker (iterative rounds + combined search)\n");
			} else {
				currentTask = new MetaDiaMarkerTask(analysisParams, mrp,
						consolePanel.getOverallProgressBar(), consolePanel.getCurrentTaskProgressBar());
				consolePanel.appendConsoleText("Workflow: Marker\n");
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

	@Override
	protected void viewReport() {
		if (currentAnalysisParams != null) {
			String resultDir = currentAnalysisParams.getResult();
			if (resultDir != null && !resultDir.isEmpty()) {
				File reportDir = new File(resultDir, "report");
				if (reportDir.exists()) {
					try {
						java.awt.Desktop.getDesktop().open(reportDir);
					} catch (Exception e) {
						LOGGER.error("Failed to open report directory", e);
						showError("Error", "Failed to open report directory: " + e.getMessage());
					}
				} else {
					showError("Report Not Found", "Report directory does not exist: " + reportDir.getAbsolutePath());
				}
			}
		}
	}

	@Override
	public void loadParametersFromObject(Object parameters) {
		if (!(parameters instanceof MetaParameterDiaFocus)) {
			LOGGER.warn("Invalid parameter type for DIA Marker: {}", 
				parameters != null ? parameters.getClass().getSimpleName() : "null");
			return;
		}

		MetaParameterDiaFocus params = (MetaParameterDiaFocus) parameters;

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

		// Load thread count
		if (threadsComboBox != null) {
			int threadCount = params.getThreadCount();
			if (threadCount > 0) {
				threadsComboBox.setValue(threadCount);
			}
		}

		// Modifications are stored in MetaParameterDiaFocus
		updateParameterView();
		updateMarker4Status();
	}
}
