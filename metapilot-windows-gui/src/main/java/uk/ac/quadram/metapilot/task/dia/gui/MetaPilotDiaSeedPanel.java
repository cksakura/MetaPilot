package uk.ac.quadram.metapilot.task.dia.gui;

import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.control.*;
import javafx.scene.layout.*;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.scene.text.Font;
import javafx.scene.text.FontWeight;
import javafx.scene.paint.Color;
import javafx.application.Platform;
import javafx.scene.control.Alert.AlertType;
import java.util.Optional;

import uk.ac.quadram.metapilot.task.common.gui.MetaPilotIOPanelFX;
import uk.ac.quadram.metapilot.task.common.gui.MagDatabaseSettingsPanel;
import uk.ac.quadram.metapilot.task.common.gui.MetadataConfigurationPanel;
import uk.ac.quadram.metapilot.task.common.gui.AbstractWorkflowPanel;
import uk.ac.quadram.metapilot.task.common.gui.ConsolePanel;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.resources.MetaResourceIO;
import uk.ac.quadram.metapilot.task.dda.igc.gui.MetaPilotDdaIgcMaxQuantPanel;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;
import uk.ac.quadram.metapilot.task.mag.MagDbItem.MagVersion;
import uk.ac.quadram.metapilot.task.par.MetaBaseParameter;
import uk.ac.quadram.metapilot.task.par.MetaData;
import uk.ac.quadram.metapilot.task.dia.par.DiaParameterIO;
import uk.ac.quadram.metapilot.task.dia.par.MetaParameterDia;
import uk.ac.quadram.metapilot.task.dia.par.MetaParameterDiaMarker;
import uk.ac.quadram.metapilot.task.dia.MetaDiaTask;
import uk.ac.quadram.metapilot.task.dia.MetaDiaSeedDirectTask;
import uk.ac.quadram.metapilot.task.dia.MetaDiaSeedIterativeTask;
import uk.ac.quadram.metapilot.task.dia.MetaDiaDeepComprehensiveTask;
import uk.ac.quadram.metapilot.dbSearch.deepDetect.DeepDetectTask;
import uk.ac.quadram.metapilot.dbSearch.deepDetect.DeepDetectParameter;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.File;
import java.time.LocalDateTime;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

/**
 * DIA Seed workflow panel - library-based identification with seed spectral library.
 * Library table is placed in the right column for better view with many items.
 */
public class MetaPilotDiaSeedPanel extends AbstractWorkflowPanel<MetaParameterDia>
		implements DiaLibraryPanelNoModel.ParentPanel {

	private static final Logger LOGGER = LogManager.getLogger(MetaPilotDiaSeedPanel.class);

	private MetaPilotIOPanelFX ioPanel;
	private MagDatabaseSettingsPanel magDatabasePanel;
	private VBox magDbStatusPanel;
	private TextArea predictionConsole;
	private ProgressBar predictionProgressBar;
	private Button predictButton;
	private TextField timeTextField;

	private DeepDetectTask deepDetectTask;
	private ExecutorService taskExecutor;
	private MetaResourceParameters resourceParams;
	private MetaDiaTask currentTask;

	private ComboBox<String> modeSelectionComboBox;
	private ComboBox<Integer> threadsComboBox;
	private CheckBox usePeptideRoundCheckBox;
	private TextField detectabilityThresholdField;

	// Modification checkboxes (same as DIA Marker Parameters tab)
	private CheckBox mExcisionCheckBox;
	private CheckBox carbamiCCheckBox;
	private CheckBox oxidationMCheckBox;
	private CheckBox acetyProNCheckBox;
	private CheckBox phosphoSTYCheckBox;
	private CheckBox kGGCheckBox;

	private DiaLibraryPanelNoModel libraryPanel;
	private MetadataConfigurationPanel metadataConfigurationPanel;

	private VBox leftColumn;
	private VBox rightColumn;
	private VBox libraryControlsPlaceholder;
	private VBox libraryTablePlaceholder;

	public MetaPilotDiaSeedPanel(MagDbItem[] magDbItems, MetaResourceParameters resourceParams) {
		super();
		parameterIO = new DiaParameterIO();
		this.resourceParams = resourceParams;
		initializeSpecificUI();
		setAvailableMagDatabases(magDbItems);
	}

	@Override
	public MetaParameterDia collectParametersFromUI() {
		try {
			if (consolePanel != null) {
				consolePanel.appendConsoleText("Collecting parameters from UI...\n");
			}
			LOGGER.info("Starting parameter collection from UI for DIA Seed workflow");

			if (!validateInputs()) {
				if (consolePanel != null) {
					consolePanel.appendConsoleText("ERROR: Input validation failed\n");
				}
				LOGGER.error("Input validation failed");
				return null;
			}
			if (consolePanel != null) {
				consolePanel.appendConsoleText("✓ Input validation passed\n");
			}

			ObservableList<MetaPilotIOPanelFX.ExperimentRow> fileRows = ioPanel.getExperimentData();
			if (fileRows.isEmpty()) {
				if (consolePanel != null) {
					consolePanel.appendConsoleText("ERROR: No experiment data found\n");
				}
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

			MetaData metaData = new MetaData(rawFiles, expNames, fractions, replicates);

			List<String> colNames = metadataConfigurationPanel != null ? metadataConfigurationPanel.getMetadataColumnNames() : java.util.Collections.emptyList();
			int metaTypeCount = Math.max(1, colNames.size());
			List<MetaPilotDdaIgcMaxQuantPanel.MetadataRow> metadataTableData = metadataConfigurationPanel != null ? metadataConfigurationPanel.getMetadataTableData() : java.util.Collections.emptyList();
			int rowCount = metadataTableData.size();
			String[][] metaInfo = new String[rowCount][metaTypeCount];
			if (!colNames.isEmpty()) {
				for (int i = 0; i < rowCount; i++) {
					MetaPilotDdaIgcMaxQuantPanel.MetadataRow row = metadataTableData.get(i);
					for (int j = 0; j < colNames.size(); j++) {
						String v = row.getValueByColumnName(colNames.get(j));
						metaInfo[i][j] = v != null ? v : "";
					}
				}
			}
			metaData.setMetaTypeCount(metaTypeCount);
			metaData.setMetaInfo(metaInfo);
			if (!colNames.isEmpty()) {
				metaData.setMetaTypeTitles(colNames.toArray(new String[0]));
			}

			String result = ioPanel.getOutputDirectory();
			if (result == null || result.trim().isEmpty()) {
				showError("Validation Error", "Please select an output directory.");
				return null;
			}

			MagDbItem selectedMagDb = magDatabasePanel.getSelectedMagDatabase();
			if (selectedMagDb == null) {
				showError("Validation Error", "Please select a MAG database.");
				return null;
			}

			MagVersion selectedVersion = magDatabasePanel.getSelectedMagVersion();
			if (selectedVersion == null) {
				showError("Validation Error", "Please select a MAG database version.");
				return null;
			}

			selectedMagDb.setUsedVersion(selectedVersion);

			boolean appendHostDb = magDatabasePanel.isAppendHostDb();
			String hostDB = magDatabasePanel.getHostDatabase();
			int threadCount = threadsComboBox != null ? threadsComboBox.getValue() : 4;

			MetaBaseParameter metaBasePar = new MetaBaseParameter(
					uk.ac.quadram.metapilot.workflow.WorkflowType.DIA_SEED.getAnalysisType(),
					metaData, result, selectedMagDb.getCurrentFile().getAbsolutePath(),
					appendHostDb, hostDB, "DIA", 1, threadCount, true);

			String[] library = new String[0];
			if (libraryPanel != null) {
				String[] selectedLibs = libraryPanel.getSelectedLibraryNames();
				if (selectedLibs != null && selectedLibs.length > 0) {
					library = selectedLibs;
				}
			}

			// Workflow: 0 = Direct search, 1 = Focus search, 2 = Comprehensive search, 3 = Deep comprehensive
			int modeIndex = modeSelectionComboBox != null ? modeSelectionComboBox.getSelectionModel().getSelectedIndex() : 0;
			boolean librarySearch = (modeIndex == 0); // Direct search only
			boolean selfModelSearch = false;
			boolean useModel = false;
			String[] models = new String[0];

			MetaParameterDia metaParameterDia = new MetaParameterDia(
					metaBasePar, selectedMagDb, librarySearch, selfModelSearch, library, useModel, models);
			
			// Persist seed workflow mode so CLI can reproduce GUI behaviour
			if (modeIndex == 0) {
				metaParameterDia.setSeedMode("DIRECT");
			} else if (modeIndex == 1) {
				metaParameterDia.setSeedMode("FOCUS");
			} else if (modeIndex == 2) {
				metaParameterDia.setSeedMode("COMPREHENSIVE");
			} else {
				metaParameterDia.setSeedMode("DEEP");
			}

			// Modifications (used for library generation and DIA-NN search)
			if (mExcisionCheckBox != null) metaParameterDia.setMExcision(mExcisionCheckBox.isSelected());
			if (carbamiCCheckBox != null) metaParameterDia.setCarbamiC(carbamiCCheckBox.isSelected());
			if (oxidationMCheckBox != null) metaParameterDia.setOxidationM(oxidationMCheckBox.isSelected());
			if (acetyProNCheckBox != null) metaParameterDia.setAcetyProN(acetyProNCheckBox.isSelected());
			if (phosphoSTYCheckBox != null) metaParameterDia.setPhosphoSTY(phosphoSTYCheckBox.isSelected());
			if (kGGCheckBox != null) metaParameterDia.setKGG(kGGCheckBox.isSelected());
			boolean peptideEnabledByMode = modeIndex != 0;
			if (usePeptideRoundCheckBox != null) {
				metaParameterDia.setUsePeptideBasedRound(peptideEnabledByMode && usePeptideRoundCheckBox.isSelected());
			}
			metaParameterDia.setPeptideDetectabilityThreshold(parseDetectabilityThresholdFromUI());

			LOGGER.info("Successfully collected parameters from UI for DIA Seed workflow");
			if (consolePanel != null) {
				consolePanel.appendConsoleText("✓ Successfully collected parameters from UI\n");
			}
			return metaParameterDia;

		} catch (Exception e) {
			String errorMsg = "ERROR: Failed to collect parameters: " + e.getMessage();
			if (consolePanel != null) {
				consolePanel.appendConsoleText(errorMsg + "\n");
			}
			LOGGER.error("Error generating MetaParameterDia: " + e.getMessage(), e);
			return null;
		}
	}

	public void setAvailableMagDatabases(MagDbItem[] magDbItems) {
		if (magDatabasePanel != null) {
			magDatabasePanel.setAvailableMagDatabases(magDbItems);
			createLibraryPanel();
			MagDbItem selected = magDatabasePanel.getSelectedMagDatabase();
			if (selected != null) {
				updateLibraryPanelWithMagDb(selected);
			} else if (magDbItems != null && magDbItems.length > 0) {
				updateLibraryPanelWithMagDb(magDbItems[0]);
			}
			updateMagDbStatus();
		}
	}

	private void updateLibraryPanelWithMagDb(MagDbItem magDbItem) {
		if (libraryPanel != null) {
			libraryPanel.setMagDbItem(magDbItem);
		}
	}

	private void createLibraryPanel() {
		if (libraryPanel != null) {
			return;
		}

		loadLastUsedParameters();
		MetaParameterDia initialParams = generateMetaParameterDiaForLibraryGeneration();
		if (initialParams != null) {
			libraryPanel = new DiaLibraryPanelNoModel(initialParams, resourceParams, new File("MAGDB_1_0.json"));
			libraryPanel.setParameterUpdateCallback(() -> updateParameterView());
			libraryPanel.setParentPanel(this);
			// Ensure saved library selections are applied even when the panel is created after load
			try {
				Object loaded = parameterIO.loadLastUsedParameters();
				if (loaded instanceof MetaParameterDia) {
					libraryPanel.updateMetaParameter((MetaParameterDia) loaded);
				}
			} catch (Exception e) {
				LOGGER.debug("Failed to apply saved library selections to library panel: {}", e.getMessage());
			}

			replaceLibraryPlaceholders();
		}
	}

	private void ensureLibraryPanelCreated() {
		if (libraryPanel != null) {
			return;
		}
		if (magDatabasePanel == null || magDatabasePanel.getSelectedMagDatabase() == null) {
			return;
		}
		createLibraryPanel();
	}

	private void replaceLibraryPlaceholders() {
		if (libraryPanel == null || leftColumn == null || rightColumn == null) {
			return;
		}
		if (libraryControlsPlaceholder != null) {
			int idx = leftColumn.getChildren().indexOf(libraryControlsPlaceholder);
			if (idx >= 0) {
				VBox controlsSection = libraryPanel.getControlsSection();
				leftColumn.getChildren().set(idx, controlsSection);
				VBox.setVgrow(controlsSection, Priority.ALWAYS);
			}
		}
		if (libraryTablePlaceholder != null) {
			int idx = rightColumn.getChildren().indexOf(libraryTablePlaceholder);
			if (idx >= 0) {
				rightColumn.getChildren().set(idx, libraryPanel.getTableSection());
			}
		}
	}

	private MetaParameterDia generateMetaParameterDiaForLibraryGeneration() {
		try {
			if (magDatabasePanel == null || magDatabasePanel.getSelectedMagDatabase() == null) {
				return null;
			}
			if (magDatabasePanel.getSelectedMagVersion() == null) {
				return null;
			}

			MetaData metaData = new MetaData(new String[0], new String[0], new int[0], new int[0]);
			String result = ioPanel != null ? ioPanel.getOutputDirectory() : "";
			MagDbItem selectedMagDb = magDatabasePanel.getSelectedMagDatabase();
			MagVersion selectedVersion = magDatabasePanel.getSelectedMagVersion();
			selectedMagDb.setUsedVersion(selectedVersion);

			boolean appendHostDb = magDatabasePanel.isAppendHostDb();
			String hostDB = magDatabasePanel.getHostDatabase();
			int threadCount = threadsComboBox != null ? threadsComboBox.getValue() : 4;

			MetaBaseParameter metaBasePar = new MetaBaseParameter(
					uk.ac.quadram.metapilot.workflow.WorkflowType.DIA_SEED.getAnalysisType(),
					metaData, result,
					selectedMagDb.getCurrentFile() != null ? selectedMagDb.getCurrentFile().getAbsolutePath() : "",
					appendHostDb, hostDB, "DDA", 1, threadCount, true);

			return new MetaParameterDia(metaBasePar, selectedMagDb, false, false, new String[0], false, new String[0]);
		} catch (Exception e) {
			LOGGER.warn("Error generating MetaParameterDia for library generation: " + e.getMessage());
			return null;
		}
	}

	@Override
	protected String getWorkflowName() {
		return "DIA Seed";
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
				magDatabasePanel,
				createSectionTitle("MAG Library Status"),
				magDbStatusPanel
		);
		return tabContent;
	}

	@Override
	public VBox createParametersTab() {
		VBox tabContent = new VBox(10);
		tabContent.setPadding(new Insets(10));
		tabContent.setMinWidth(960);
		tabContent.setMinHeight(840);

		GridPane mainGrid = new GridPane();
		mainGrid.setHgap(15);
		mainGrid.setVgap(15);
		mainGrid.setPadding(new Insets(10));
		mainGrid.setMinWidth(940);
		mainGrid.setMinHeight(740);
		VBox.setVgrow(mainGrid, Priority.ALWAYS);

		ColumnConstraints col1 = new ColumnConstraints();
		col1.setPercentWidth(50);
		col1.setMinWidth(450);
		col1.setHgrow(Priority.ALWAYS);
		ColumnConstraints col2 = new ColumnConstraints();
		col2.setPercentWidth(50);
		col2.setMinWidth(450);
		col2.setHgrow(Priority.ALWAYS);
		mainGrid.getColumnConstraints().addAll(col1, col2);

		RowConstraints row1 = new RowConstraints();
		row1.setVgrow(Priority.ALWAYS);
		row1.setMinHeight(400);
		mainGrid.getRowConstraints().add(row1);

		// Left column: DIA Settings, Modifications, Library controls (Spectral library generation) - same column, modifications above progress bar
		leftColumn = new VBox(15);
		leftColumn.setMinWidth(430);
		VBox diaSettingsSection = createDiaSettingsSection();
		VBox modificationsSection = createModificationsSection();

		libraryControlsPlaceholder = new VBox();
		libraryControlsPlaceholder.setMinHeight(100);
		Label controlsPlaceholderLabel = new Label("Library Generation (available after MAG database is selected)");
		controlsPlaceholderLabel.setStyle("-fx-text-alignment: center; -fx-font-style: italic; -fx-text-fill: #666666;");
		libraryControlsPlaceholder.getChildren().add(controlsPlaceholderLabel);
		libraryControlsPlaceholder.setAlignment(Pos.CENTER);

		leftColumn.getChildren().addAll(diaSettingsSection, modificationsSection, libraryControlsPlaceholder);
		VBox.setVgrow(libraryControlsPlaceholder, Priority.ALWAYS);
		GridPane.setVgrow(leftColumn, Priority.ALWAYS);

		// Right column: Library table (moved here for better view) + Metadata
		rightColumn = new VBox(15);
		rightColumn.setMinWidth(430);

		libraryTablePlaceholder = createPlaceholderTableSection();
		metadataConfigurationPanel = new MetadataConfigurationPanel();
		metadataConfigurationPanel.setParameterUpdateCallback(() -> updateParameterView());

		rightColumn.getChildren().addAll(libraryTablePlaceholder, metadataConfigurationPanel);
		GridPane.setVgrow(rightColumn, Priority.ALWAYS);

		mainGrid.add(leftColumn, 0, 0);
		mainGrid.add(rightColumn, 1, 0);

		tabContent.getChildren().add(mainGrid);
		return tabContent;
	}

	private VBox createPlaceholderTableSection() {
		VBox section = new VBox(10);
		section.setPadding(new Insets(10));
		section.setStyle("-fx-background-color: #f8f8f8; -fx-background-radius: 5;");
		section.setMinHeight(200);
		Label placeholderLabel = new Label("Library Table\n(available after MAG database is selected)");
		placeholderLabel.setStyle("-fx-text-alignment: center; -fx-font-style: italic; -fx-text-fill: #666666;");
		section.getChildren().add(placeholderLabel);
		section.setAlignment(Pos.CENTER);
		return section;
	}

	private VBox createDiaSettingsSection() {
		VBox section = new VBox(10);
		section.setPadding(new Insets(10));
		section.setStyle("-fx-background-color: #f0f8ff; -fx-background-radius: 5;");

		Label titleLabel = createSectionTitle("DIA Analysis Settings");

		GridPane settingsGrid = new GridPane();
		settingsGrid.setHgap(15);
		settingsGrid.setVgap(8);

		ColumnConstraints col1 = new ColumnConstraints();
		col1.setMinWidth(150);
		col1.setHgrow(Priority.ALWAYS);
		ColumnConstraints col2 = new ColumnConstraints();
		col2.setMinWidth(150);
		col2.setHgrow(Priority.ALWAYS);
		settingsGrid.getColumnConstraints().addAll(col1, col2);

		Label modeLabel = new Label("Workflow:");
		modeLabel.setStyle("-fx-font-weight: bold;");
		modeSelectionComboBox = new ComboBox<>(FXCollections.observableArrayList(
				"Direct search",
				"Focus search",
				"Comprehensive search",
				"Deep comprehensive search"));
		modeSelectionComboBox.setValue("Direct search");
		modeSelectionComboBox.setTooltip(new Tooltip(
				"Direct search: identify with selected libraries only (no genome selection).\n" +
				"Focus search: iterative workflow with focus genome selection (fewer, high-confidence genomes).\n" +
				"Comprehensive search: iterative workflow with comprehensive genome selection (more genomes).\n" +
				"Deep comprehensive search: iterative workflow with deep learning-based library refinement."));
		modeSelectionComboBox.setMinWidth(200);
		modeSelectionComboBox.setPrefWidth(300);
		modeSelectionComboBox.setMaxWidth(Double.MAX_VALUE);
		modeSelectionComboBox.setOnAction(e -> updateParameterView());

		VBox modeBox = new VBox(5);
		modeBox.getChildren().addAll(modeLabel, modeSelectionComboBox);

		Label threadsLabel = new Label("Thread Count:");
		threadsLabel.setStyle("-fx-font-weight: bold;");
		threadsComboBox = new ComboBox<>();
		for (int i = 1; i <= Runtime.getRuntime().availableProcessors(); i++) {
			threadsComboBox.getItems().add(i);
		}
		threadsComboBox.setValue(Math.min(4, Runtime.getRuntime().availableProcessors()));
		threadsComboBox.setMinWidth(100);
		threadsComboBox.setPrefWidth(120);
		threadsComboBox.setMaxWidth(Double.MAX_VALUE);
		threadsComboBox.setOnAction(e -> updateParameterView());

		VBox threadsBox = new VBox(5);
		threadsBox.getChildren().addAll(threadsLabel, threadsComboBox);

		settingsGrid.add(modeBox, 0, 0);
		settingsGrid.add(threadsBox, 1, 0);

		VBox deepPeptideBox = new VBox(5);
		Label deepPeptideLabel = new Label("Peptide-based library settings:");
		deepPeptideLabel.setStyle("-fx-font-weight: bold;");
		usePeptideRoundCheckBox = new CheckBox("Use peptide-based rounds");
		usePeptideRoundCheckBox.setSelected(false);
		usePeptideRoundCheckBox.setTooltip(new Tooltip(
				"For Focus/Comprehensive/Deep modes: build iterative microbiome library from predicted peptides."));
		usePeptideRoundCheckBox.setOnAction(e -> updateParameterView());

		HBox thresholdBox = new HBox(8);
		thresholdBox.setAlignment(Pos.CENTER_LEFT);
		Label thresholdLabel = new Label("Detectability threshold:");
		detectabilityThresholdField = new TextField("0.5");
		detectabilityThresholdField.setPrefWidth(90);
		detectabilityThresholdField.setTooltip(new Tooltip(
				"Keep predicted peptides with detectability >= this value."));
		detectabilityThresholdField.textProperty().addListener((obs, ov, nv) -> updateParameterView());
		thresholdBox.getChildren().addAll(thresholdLabel, detectabilityThresholdField);
		deepPeptideBox.getChildren().addAll(deepPeptideLabel, usePeptideRoundCheckBox, thresholdBox);
		updatePeptideBasedControlsState();

		section.getChildren().addAll(titleLabel, settingsGrid, deepPeptideBox);
		return section;
	}

	private VBox createModificationsSection() {
		VBox section = new VBox(10);
		section.setPadding(new Insets(10));
		section.setStyle("-fx-background-color: #f0f8ff; -fx-background-radius: 5;");

		Label titleLabel = createSectionTitle("Modifications");

		VBox checkboxesBox = new VBox(8);

		mExcisionCheckBox = new CheckBox("M excision");
		mExcisionCheckBox.setSelected(true);
		mExcisionCheckBox.setOnAction(e -> updateParameterView());

		carbamiCCheckBox = new CheckBox("Carbamidomethyl C");
		carbamiCCheckBox.setSelected(true);
		carbamiCCheckBox.setOnAction(e -> updateParameterView());

		oxidationMCheckBox = new CheckBox("Oxidation M");
		oxidationMCheckBox.setSelected(false);
		oxidationMCheckBox.setOnAction(e -> updateParameterView());

		acetyProNCheckBox = new CheckBox("Acetylation protein N-term");
		acetyProNCheckBox.setSelected(false);
		acetyProNCheckBox.setOnAction(e -> updateParameterView());

		phosphoSTYCheckBox = new CheckBox("Phosphorylation STY");
		phosphoSTYCheckBox.setSelected(false);
		phosphoSTYCheckBox.setOnAction(e -> updateParameterView());

		kGGCheckBox = new CheckBox("K GG");
		kGGCheckBox.setSelected(false);
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

	private VBox createMagDbStatusPanel() {
		VBox statusPanel = new VBox(10);
		statusPanel.setPadding(new Insets(10));
		statusPanel.setStyle("-fx-background-color: #f8f8f8; -fx-border-color: #cccccc; -fx-border-radius: 5; -fx-border-width: 1;");
		statusPanel.setMinHeight(200);
		statusPanel.setPrefHeight(300);

		// Console output panel for prediction process
		VBox consolePanel = new VBox(5);
		consolePanel.setPadding(new Insets(5));
		consolePanel.setStyle("-fx-background-color: #f0f0f0; -fx-background-radius: 3; -fx-border-color: #cccccc; -fx-border-radius: 3; -fx-border-width: 1;");

		Label consoleTitle = new Label("MAG Database Status Console");
		consoleTitle.setStyle("-fx-font-weight: bold; -fx-font-size: 12px;");

		predictionConsole = new TextArea();
		predictionConsole.setEditable(false);
		predictionConsole.setWrapText(true);
		predictionConsole.setMinHeight(150);
		predictionConsole.setPrefHeight(180);
		predictionConsole.setMaxHeight(200);
		predictionConsole.setStyle("-fx-font-family: 'Courier New', monospace; -fx-font-size: 11px;");
		predictionConsole.setPromptText("Database status information will appear here...");
		predictionConsole.setId("predictionConsole");

		consolePanel.getChildren().addAll(consoleTitle, predictionConsole);

		// Prediction controls panel
		VBox predictionPanel = new VBox(5);
		predictionPanel.setPadding(new Insets(5));
		predictionPanel.setStyle("-fx-background-color: #f0f8ff; -fx-background-radius: 3; -fx-border-color: #b0d4f1; -fx-border-radius: 3; -fx-border-width: 1;");

		Label predictionLabel = new Label("Peptide Prediction Controls");
		predictionLabel.setStyle("-fx-font-weight: bold; -fx-font-size: 12px;");

		HBox predictionControls = new HBox(10);
		predictionControls.setAlignment(Pos.CENTER_LEFT);

		predictButton = new Button("Start Prediction");
		predictButton.setMinWidth(100);
		predictButton.setPrefWidth(120);
		predictButton.setMaxWidth(140);
		predictButton.setDisable(true);
		predictButton.setOnAction(e -> startPrediction());

		predictionProgressBar = new ProgressBar(0);
		predictionProgressBar.setMaxWidth(Double.MAX_VALUE);
		HBox.setHgrow(predictionProgressBar, Priority.ALWAYS);

		timeTextField = new TextField();
		timeTextField.setEditable(false);
		timeTextField.setMinWidth(100);
		timeTextField.setPrefWidth(120);
		timeTextField.setPromptText("Status");

		predictionControls.getChildren().addAll(predictButton, predictionProgressBar, timeTextField);

		predictionPanel.getChildren().addAll(predictionControls);

		statusPanel.getChildren().addAll(consolePanel, predictionPanel);
		return statusPanel;
	}

	private Label createSectionTitle(String title) {
		Label label = new Label(title);
		label.setFont(Font.font("System", FontWeight.BOLD, 14));
		label.setTextFill(Color.DARKBLUE);
		return label;
	}

	private void initializeSpecificUI() {
		setSpacing(10);
		setPadding(new Insets(10));
		setMinWidth(1200);
		setMinHeight(900);
		setPrefWidth(1600);
		setPrefHeight(1000);
		setMaxHeight(Double.MAX_VALUE);

		consolePanel = new ConsolePanel();
		magDatabasePanel = new MagDatabaseSettingsPanel(false);
		magDatabasePanel.setParameterUpdateCallback(() -> {
			syncSelectedMagDbToPanels();
			updateParameterView();
			updateMagDbStatus();
		});

		consolePanel.setStartCallback(this::startAnalysis);
		consolePanel.setStopCallback(this::stopAnalysis);
		consolePanel.setViewReportCallback(this::viewReport);

		magDbStatusPanel = createMagDbStatusPanel();

		TabPane tabPane = new TabPane();
		tabPane.setTabClosingPolicy(TabPane.TabClosingPolicy.UNAVAILABLE);
		tabPane.setMinWidth(1080);
		tabPane.setMinHeight(900);
		VBox.setVgrow(tabPane, Priority.ALWAYS);

		Tab ioTab = new Tab("Input/Output");
		ioTab.setContent(createIOTab());

		Tab paramTab = new Tab("Parameters");
		paramTab.setContent(createParametersTab());

		Tab consoleTab = new Tab("Console");
		consoleTab.setContent(consolePanel);

		tabPane.getTabs().addAll(ioTab, paramTab, consoleTab);

		tabPane.getSelectionModel().selectedItemProperty().addListener((obs, oldTab, newTab) -> {
			if (newTab == paramTab) {
				syncSelectedMagDbToPanels();
				updateMetadataTableFromIOPanel();
				ensureLibraryPanelCreated();
				updateParameterView();
			} else if (newTab == consoleTab) {
				syncSelectedMagDbToPanels();
				updateParameterView();
				// Auto-save without validation or alerts (validation/alert only when user clicks Run)
				try {
					MetaParameterDia currentParams = collectParametersFromUIForAutoSave();
					if (currentParams == null) {
						currentParams = generateMetaParameterDiaForLibraryGeneration();
					}
					if (currentParams != null) {
						parameterIO.saveLastUsedParameters(currentParams);
					}
				} catch (Exception e) {
					LOGGER.warn("Failed to auto-save parameters", e);
				}
			}
		});

		getChildren().add(tabPane);

		Platform.runLater(this::setupWindowCloseHandler);
	}

	/**
	 * Keep Parameters/Console panels aligned with Input/Output MAG DB selection.
	 * Input/Output (magDatabasePanel) is the single source of truth.
	 */
	private void syncSelectedMagDbToPanels() {
		if (magDatabasePanel == null) {
			return;
		}
		MagDbItem selectedMagDb = magDatabasePanel.getSelectedMagDatabase();
		if (selectedMagDb != null) {
			updateLibraryPanelWithMagDb(selectedMagDb);
		}
	}

	private void updateMetadataTableFromIOPanel() {
		if (ioPanel == null || metadataConfigurationPanel == null) {
			return;
		}
		ObservableList<MetaPilotIOPanelFX.ExperimentRow> experiments = ioPanel.getExperimentData();
		metadataConfigurationPanel.updateMetadataFromIOPanel(experiments);
	}

	@Override
	protected void updateParameterView() {
		if (consolePanel == null) {
			return;
		}
		updatePeptideBasedControlsState();

		StringBuilder sb = new StringBuilder();
		sb.append("=== DIA Analysis Parameters ===\n\n");

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
		}
		sb.append("\n");

		sb.append("DIA Settings:\n");
		if (modeSelectionComboBox != null) {
			sb.append("  Workflow: ").append(modeSelectionComboBox.getValue()).append("\n");
		}
		if (threadsComboBox != null) {
			sb.append("  Thread Count: ").append(threadsComboBox.getValue()).append("\n");
		}
		if (usePeptideRoundCheckBox != null) {
			sb.append("  Use peptide rounds: ").append(usePeptideRoundCheckBox.isSelected()).append("\n");
			sb.append("  Detectability threshold: ").append(parseDetectabilityThresholdFromUI()).append("\n");
		}
		sb.append("\n");

		sb.append("Modifications:\n");
		if (mExcisionCheckBox != null) sb.append("  M excision: ").append(mExcisionCheckBox.isSelected()).append("\n");
		if (carbamiCCheckBox != null) sb.append("  Carbamidomethyl C: ").append(carbamiCCheckBox.isSelected()).append("\n");
		if (oxidationMCheckBox != null) sb.append("  Oxidation M: ").append(oxidationMCheckBox.isSelected()).append("\n");
		if (acetyProNCheckBox != null) sb.append("  Acetylation protein N-term: ").append(acetyProNCheckBox.isSelected()).append("\n");
		if (phosphoSTYCheckBox != null) sb.append("  Phosphorylation STY: ").append(phosphoSTYCheckBox.isSelected()).append("\n");
		if (kGGCheckBox != null) sb.append("  K GG: ").append(kGGCheckBox.isSelected()).append("\n");
		sb.append("\n");

		sb.append("Spectral Libraries:\n");
		if (libraryPanel != null) {
			String[] selectedLibs = libraryPanel.getSelectedLibraryNames();
			if (selectedLibs != null && selectedLibs.length > 0) {
				for (String libName : selectedLibs) {
					sb.append("  ").append(libName).append("\n");
				}
			} else {
				sb.append("  No libraries selected\n");
			}
		} else {
			sb.append("  Library panel not available\n");
		}
		sb.append("\n");

		consolePanel.setParameterText(sb.toString());

		if (libraryPanel != null) {
			// Use lenient collect (no validation/alert) so switching to Parameters tab does not pop up alerts
			MetaParameterDia currentMetaParameter = collectParametersFromUIForAutoSave();
			if (currentMetaParameter == null) {
				currentMetaParameter = generateMetaParameterDiaForLibraryGeneration();
			}
			if (currentMetaParameter != null) {
				libraryPanel.updateMetaParameter(currentMetaParameter);
			}
		}
	}

	/**
	 * Collect parameters from UI without validation or alerts. Used for auto-save on tab switch
	 * and for updating the library panel. Validation and alerts only run when user clicks Run.
	 */
	private MetaParameterDia collectParametersFromUIForAutoSave() {
		try {
			if (ioPanel == null || ioPanel.getExperimentData().isEmpty()) {
				return null;
			}
			String result = ioPanel.getOutputDirectory();
			if (result == null || result.trim().isEmpty()) {
				return null;
			}
			MagDbItem selectedMagDb = magDatabasePanel != null ? magDatabasePanel.getSelectedMagDatabase() : null;
			if (selectedMagDb == null) {
				return null;
			}
			MagVersion selectedVersion = magDatabasePanel != null ? magDatabasePanel.getSelectedMagVersion() : null;
			if (selectedVersion == null) {
				return null;
			}
			selectedMagDb.setUsedVersion(selectedVersion);

			ObservableList<MetaPilotIOPanelFX.ExperimentRow> fileRows = ioPanel.getExperimentData();
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
					if (fractions[i] <= 0) fractions[i] = 1;
				} catch (NumberFormatException e) {
					fractions[i] = 1;
				}
				replicates[i] = 1;
			}
			MetaData metaData = new MetaData(rawFiles, expNames, fractions, replicates);
			List<String> colNames = metadataConfigurationPanel != null ? metadataConfigurationPanel.getMetadataColumnNames() : java.util.Collections.emptyList();
			List<MetaPilotDdaIgcMaxQuantPanel.MetadataRow> metadataTableData = metadataConfigurationPanel != null ? metadataConfigurationPanel.getMetadataTableData() : java.util.Collections.emptyList();
			int rowCount = metadataTableData.size();
			int metaTypeCount = Math.max(1, colNames.size());
			String[][] metaInfo = new String[rowCount][metaTypeCount];
			if (!colNames.isEmpty()) {
				for (int i = 0; i < rowCount; i++) {
					MetaPilotDdaIgcMaxQuantPanel.MetadataRow row = metadataTableData.get(i);
					for (int j = 0; j < colNames.size(); j++) {
						String v = row.getValueByColumnName(colNames.get(j));
						metaInfo[i][j] = v != null ? v : "";
					}
				}
			}
			metaData.setMetaTypeCount(metaTypeCount);
			metaData.setMetaInfo(metaInfo);
			if (!colNames.isEmpty()) {
				metaData.setMetaTypeTitles(colNames.toArray(new String[0]));
			}

			boolean appendHostDb = magDatabasePanel != null && magDatabasePanel.isAppendHostDb();
			String hostDB = magDatabasePanel != null ? magDatabasePanel.getHostDatabase() : "";
			int threadCount = threadsComboBox != null && threadsComboBox.getValue() != null ? threadsComboBox.getValue() : 4;
			MetaBaseParameter metaBasePar = new MetaBaseParameter(
					uk.ac.quadram.metapilot.workflow.WorkflowType.DIA_SEED.getAnalysisType(),
					metaData, result, selectedMagDb.getCurrentFile() != null ? selectedMagDb.getCurrentFile().getAbsolutePath() : "",
					appendHostDb, hostDB, "DIA", 1, threadCount, true);

			String[] library = new String[0];
			if (libraryPanel != null) {
				String[] selectedLibs = libraryPanel.getSelectedLibraryNames();
				if (selectedLibs != null && selectedLibs.length > 0) {
					library = selectedLibs;
				}
			}
			int modeIndex = modeSelectionComboBox != null ? modeSelectionComboBox.getSelectionModel().getSelectedIndex() : 0;
			boolean librarySearch = (modeIndex == 0);
			MetaParameterDia metaParameterDia = new MetaParameterDia(
					metaBasePar, selectedMagDb, librarySearch, false, library, false, new String[0]);
			// Keep auto-saved mode consistent with Run-mode parameter collection
			if (modeIndex == 0) {
				metaParameterDia.setSeedMode("DIRECT");
			} else if (modeIndex == 1) {
				metaParameterDia.setSeedMode("FOCUS");
			} else if (modeIndex == 2) {
				metaParameterDia.setSeedMode("COMPREHENSIVE");
			} else {
				metaParameterDia.setSeedMode("DEEP");
			}
			if (mExcisionCheckBox != null) metaParameterDia.setMExcision(mExcisionCheckBox.isSelected());
			if (carbamiCCheckBox != null) metaParameterDia.setCarbamiC(carbamiCCheckBox.isSelected());
			if (oxidationMCheckBox != null) metaParameterDia.setOxidationM(oxidationMCheckBox.isSelected());
			if (acetyProNCheckBox != null) metaParameterDia.setAcetyProN(acetyProNCheckBox.isSelected());
			if (phosphoSTYCheckBox != null) metaParameterDia.setPhosphoSTY(phosphoSTYCheckBox.isSelected());
			if (kGGCheckBox != null) metaParameterDia.setKGG(kGGCheckBox.isSelected());
			boolean peptideEnabledByMode = modeIndex != 0;
			if (usePeptideRoundCheckBox != null) {
				metaParameterDia.setUsePeptideBasedRound(peptideEnabledByMode && usePeptideRoundCheckBox.isSelected());
			}
			metaParameterDia.setPeptideDetectabilityThreshold(parseDetectabilityThresholdFromUI());
			return metaParameterDia;
		} catch (Exception e) {
			LOGGER.debug("Lenient parameter collect failed: " + e.getMessage());
			return null;
		}
	}

	private double parseDetectabilityThresholdFromUI() {
		if (detectabilityThresholdField == null) {
			return 0.5;
		}
		try {
			return Double.parseDouble(detectabilityThresholdField.getText().trim());
		} catch (Exception e) {
			return 0.5;
		}
	}

	private void updatePeptideBasedControlsState() {
		if (usePeptideRoundCheckBox == null || detectabilityThresholdField == null || modeSelectionComboBox == null) {
			return;
		}
		int modeIndex = modeSelectionComboBox.getSelectionModel().getSelectedIndex();
		boolean peptideSectionEnabled = modeIndex != 0; // disable in Direct search
		usePeptideRoundCheckBox.setDisable(!peptideSectionEnabled);
		detectabilityThresholdField
				.setDisable(!(peptideSectionEnabled && usePeptideRoundCheckBox.isSelected()));
	}

	private boolean validateInputs() {
		if (ioPanel == null || ioPanel.getExperimentData().isEmpty()) {
			showError("Validation Error", "Please add at least one experiment.");
			return false;
		}
		String outputDir = ioPanel.getOutputDirectory();
		if (outputDir == null || outputDir.trim().isEmpty()) {
			showError("Validation Error", "Please select an output directory.");
			return false;
		}
		if (magDatabasePanel == null || magDatabasePanel.getSelectedMagDatabase() == null) {
			showError("Validation Error", "Please select a MAG database.");
			return false;
		}
			if (magDatabasePanel.getSelectedMagVersion() == null) {
				showError("Validation Error", "Please select a MAG database version.");
				return false;
			}
			int modeIndex = modeSelectionComboBox != null ? modeSelectionComboBox.getSelectionModel().getSelectedIndex() : 0;
			if (modeIndex != 0 && libraryPanel != null) {
				String[] libs = libraryPanel.getSelectedLibraryNames();
				if (libs == null || libs.length == 0) {
					showError("Validation Error", "For Focus or Comprehensive search, please select at least one spectral library.");
					return false;
				}
			}
			if (modeIndex != 0 && usePeptideRoundCheckBox != null && usePeptideRoundCheckBox.isSelected()) {
				double threshold;
				try {
					threshold = Double.parseDouble(detectabilityThresholdField.getText().trim());
				} catch (Exception e) {
					showError("Validation Error", "Detectability threshold must be a number.");
					return false;
				}
				if (threshold < 0.0 || threshold > 1.0) {
					showError("Validation Error", "Detectability threshold must be between 0 and 1.");
					return false;
				}
			}
			return true;
		}

	@Override
	protected void startAnalysis() {
		try {
			if (consolePanel != null) {
				consolePanel.appendConsoleText("Starting DIA Seed analysis...\n");
			}
			LOGGER.info("Starting DIA Seed analysis");

			MetaParameterDia metaParameterDia = collectParametersFromUI();
			if (metaParameterDia == null) {
				if (consolePanel != null) {
					consolePanel.appendConsoleText("ERROR: Failed to collect parameters from UI.\n");
				}
				return;
			}

			String outputDir = metaParameterDia.getResult();
			if (outputDir != null && !outputDir.trim().isEmpty()) {
				parameterIO.saveAnalysisParameters(outputDir, metaParameterDia);
				parameterIO.saveLastUsedParameters(metaParameterDia);
			}

			MetaResourceIO resourceIO = new MetaResourceIO();
			MetaResourceParameters mrp = resourceIO.loadParameters();
			if (mrp == null) {
				showError("Configuration Error", "Failed to load resource parameters from config file.");
				return;
			}

			consolePanel.clearConsole();
			consolePanel.appendConsoleText("Parameters validated successfully.\n");
			consolePanel.appendConsoleText("Result directory: " + metaParameterDia.getResult() + "\n");

			int modeIndex = modeSelectionComboBox.getSelectionModel().getSelectedIndex();
			if (modeIndex == 0) {
				// Direct seed search (no genome selection, no deep model)
				this.currentAnalysisParams = metaParameterDia;
				consolePanel.appendConsoleText("Workflow: Direct search\n");
				currentTask = new MetaDiaSeedDirectTask(metaParameterDia, mrp,
						consolePanel.getOverallProgressBar(), consolePanel.getCurrentTaskProgressBar());
			} else if (modeIndex == 3) {
				// Deep comprehensive search (seed strategy only)
				consolePanel.appendConsoleText("Workflow: Deep comprehensive search\n");
				MetaParameterDiaMarker focusParams = new MetaParameterDiaMarker(
						metaParameterDia, metaParameterDia.getMagDbItem(),
						false, false, metaParameterDia.getLibrary(), false, new String[0]);
				focusParams.setMarkerMode("DEEP");
				this.currentAnalysisParams = focusParams;
				currentTask = new MetaDiaDeepComprehensiveTask(focusParams, mrp,
						consolePanel.getOverallProgressBar(), consolePanel.getCurrentTaskProgressBar());
			} else {
				// Focus search (1) or Comprehensive search (2) – iterative seed workflow
				consolePanel.appendConsoleText(
						"Workflow: " + (modeIndex == 2 ? "Comprehensive" : "Focus") + " search\n");
				MetaParameterDiaMarker focusParams = new MetaParameterDiaMarker(
						metaParameterDia, metaParameterDia.getMagDbItem(),
						false, false, metaParameterDia.getLibrary(), false, new String[0]);
				if (modeIndex == 2) {
					focusParams.setMarkerMode("COMPREHENSIVE");
				} else {
					focusParams.setMarkerMode("FOCUS");
				}
				this.currentAnalysisParams = focusParams;
				currentTask = new MetaDiaSeedIterativeTask(focusParams, mrp,
						consolePanel.getOverallProgressBar(), consolePanel.getCurrentTaskProgressBar());
			}

			currentTask.setOnSucceeded(e -> Platform.runLater(() -> {
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
			}));

			currentTask.setOnFailed(e -> Platform.runLater(() -> {
				Throwable exception = currentTask.getException();
				String errorMessage = exception != null ? exception.getMessage() : "Unknown error occurred";
				consolePanel.appendConsoleText("✗ Analysis failed: " + errorMessage + "\n");
				consolePanel.setButtonsEnabled(true, false, false);
				consolePanel.stopStatusTimer();
				showError("Analysis Failed", "DIA Seed analysis failed: " + errorMessage);
			}));

			currentTask.setOnCancelled(e -> Platform.runLater(() -> {
				consolePanel.appendConsoleText("DIA Seed analysis was cancelled.\n");
				consolePanel.setButtonsEnabled(true, false, false);
				consolePanel.stopStatusTimer();
			}));

			currentTask.messageProperty().addListener((obs, oldMsg, newMsg) -> {
				Platform.runLater(() -> {
					if (newMsg != null) {
						consolePanel.appendConsoleText("Current Task: " + newMsg + "\n");
						consolePanel.updateCurrentTaskFromMessage(newMsg);
					}
				});
			});

			Thread taskThread = new Thread(currentTask);
			taskThread.setDaemon(true);
			taskThread.start();

			consolePanel.appendConsoleText("✓ Analysis task started successfully.\n");
			LOGGER.info("DIA Seed analysis task started");

		} catch (Exception e) {
			String errorMsg = "ERROR: Failed to start analysis: " + e.getMessage();
			consolePanel.appendConsoleText(errorMsg + "\n");
			LOGGER.error("Failed to start DIA Seed analysis", e);
			showError("Start Error", "Failed to start DIA Seed analysis: " + e.getMessage());
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
			consolePanel.appendConsoleText("Stopping DIA analysis...\n");
			if (currentTask != null && currentTask.isRunning()) {
				currentTask.forceStop();
			}
			consolePanel.setButtonsEnabled(true, false, false);
			consolePanel.stopStatusTimer();
		}
	}

	private void updateMagDbStatus() {
		if (predictionConsole == null) {
			return;
		}

		predictionConsole.clear();
		predictionConsole.appendText("=== MAG Database Status Check ===\n");
		predictionConsole.appendText("Timestamp: " + LocalDateTime.now().toString() + "\n\n");

		if (magDatabasePanel == null) {
			predictionConsole.appendText("ERROR: Database panel not initialized\n");
			predictionConsole.appendText("Please restart the application.\n");
			updateButtonStates(false);
			return;
		}

		MagDbItem magDbItem = magDatabasePanel.getSelectedMagDatabase();
		if (magDbItem == null) {
			predictionConsole.appendText("INFO: No database selected\n");
			predictionConsole.appendText("Please select a MAG database from the dropdown above.\n");
			predictionConsole.appendText("If databases are available, they should appear in the dropdown.\n");
			updateButtonStates(false);
			return;
		}

		MagVersion usedVersion = magDatabasePanel.getSelectedMagVersion();
		if (usedVersion == null) {
			predictionConsole.appendText("INFO: No version selected for database: " + magDbItem.getCatalogueID() + "\n");
			predictionConsole.appendText("Please select a version from the version dropdown.\n");
			updateButtonStates(false);
			return;
		}

		predictionConsole.appendText("Database Information:\n");
		predictionConsole.appendText("  Catalogue ID: " + magDbItem.getCatalogueID() + "\n");
		predictionConsole.appendText("  Version: " + usedVersion.getCatalogueVersion() + "\n");
		predictionConsole.appendText("  Species Count: " + usedVersion.getSpeciesCount() + "\n");
		predictionConsole.appendText("  Repository Path: " + (magDbItem.getRepositoryFile() != null ? magDbItem.getRepositoryFile().getAbsolutePath() : "Not set") + "\n\n");

		int speciesCount = usedVersion.getSpeciesCount();
		File predictedPepFile = magDbItem.getPredictedPepFile();

		if (predictedPepFile == null || !predictedPepFile.exists()) {
			predictionConsole.appendText("STATUS: Predicted peptides folder does not exist\n");
			predictionConsole.appendText("Expected path: " + (predictedPepFile != null ? predictedPepFile.getAbsolutePath() : "null") + "\n");
			predictionConsole.appendText("ACTION: Click 'Start Prediction' to generate predicted peptides\n");
			updateButtonStates(false);
			return;
		}

		File[] files = predictedPepFile.listFiles();
		if (files == null) {
			predictionConsole.appendText("ERROR: Cannot read predicted peptides folder\n");
			predictionConsole.appendText("Path: " + predictedPepFile.getAbsolutePath() + "\n");
			predictionConsole.appendText("ACTION: Check folder permissions\n");
			updateButtonStates(false);
			return;
		}

		int fileCount = files.length;
		predictionConsole.appendText("Predicted Peptides Status:\n");
		predictionConsole.appendText("  Expected files: " + speciesCount + "\n");
		predictionConsole.appendText("  Found files: " + fileCount + "\n");
		predictionConsole.appendText("  Completion: " + String.format("%.1f%%", (fileCount * 100.0 / speciesCount)) + "\n\n");

		if (fileCount == speciesCount) {
			predictionConsole.appendText("STATUS: ✓ All predicted peptides are available\n");
			predictionConsole.appendText("ACTION: Ready for analysis\n");
			updateButtonStates(true);
		} else if (fileCount > 0) {
			predictionConsole.appendText("STATUS: ⚠ Prediction incomplete (" + fileCount + "/" + speciesCount + " files)\n");
			predictionConsole.appendText("ACTION: Click 'Start Prediction' to continue\n");
			updateButtonStates(false);
		} else {
			predictionConsole.appendText("STATUS: ⚠ No predicted peptide files found\n");
			predictionConsole.appendText("ACTION: Click 'Start Prediction' to start prediction\n");
			updateButtonStates(false);
		}
	}

	private void updateButtonStates(boolean isComplete) {
		if (predictionProgressBar == null || timeTextField == null || predictButton == null) {
			return;
		}
		if (isComplete) {
			predictionProgressBar.setProgress(1.0);
			timeTextField.setText("Complete");
			predictButton.setDisable(true);
		} else {
			predictionProgressBar.setProgress(0);
			timeTextField.setText("Not Ready");
			predictButton.setDisable(false);
		}
	}

	private void startPrediction() {
		if (magDatabasePanel == null || resourceParams == null) {
			predictionConsole.appendText("ERROR: Cannot start prediction - missing database or resource parameters\n");
			return;
		}

		MagDbItem magDbItem = magDatabasePanel.getSelectedMagDatabase();
		if (magDbItem == null) {
			predictionConsole.appendText("ERROR: No database selected for prediction\n");
			return;
		}

		File predictedPepFile = magDbItem.getPredictedPepFile();
		if (predictedPepFile != null && predictedPepFile.exists()) {
			File[] files = predictedPepFile.listFiles();
			if (files != null && files.length == magDbItem.getUsedVersion().getSpeciesCount()) {
				predictionConsole.appendText("INFO: Prediction already complete. All files exist.\n");
				return;
			}
		}

		predictButton.setDisable(true);
		predictionProgressBar.setProgress(0);
		timeTextField.setText("Starting...");

		predictionConsole.clear();
		predictionConsole.appendText("=== Starting Peptide Prediction ===\n");
		predictionConsole.appendText("Timestamp: " + LocalDateTime.now().toString() + "\n");
		predictionConsole.appendText("Database: " + magDbItem.getCatalogueID() + " " + magDbItem.getUsedVersion().getCatalogueVersion() + "\n");
		predictionConsole.appendText("Species Count: " + magDbItem.getUsedVersion().getSpeciesCount() + "\n\n");

		MetaBaseParameter metaBaseParameter = new MetaBaseParameter();
		DeepDetectParameter deepDetectParam = new DeepDetectParameter();

		int availableProcessors = Runtime.getRuntime().availableProcessors();
		int threadCount = Math.max(1, (int) (availableProcessors * 0.6));

		predictionConsole.appendText("Thread Configuration:\n");
		predictionConsole.appendText("  Available Processors: " + availableProcessors + "\n");
		predictionConsole.appendText("  Using Threads: " + threadCount + " (60% of available)\n\n");

		deepDetectTask = new DeepDetectTask(metaBaseParameter, resourceParams, predictionProgressBar, deepDetectParam,
				magDbItem, threadCount,
				timeTextField,
				LOGGER);

		deepDetectTask.setOnSucceeded(e -> {
			predictionConsole.appendText("\n=== Prediction Completed ===\n");
			predictionConsole.appendText("Timestamp: " + LocalDateTime.now().toString() + "\n");
			if (deepDetectTask.getValue()) {
				predictionConsole.appendText("STATUS: ✓ All predictions completed successfully\n");
				predictionConsole.appendText("ACTION: Ready for analysis\n");
				updateButtonStates(true);
			} else {
				predictionConsole.appendText("STATUS: ⚠ Prediction completed with some issues\n");
				predictionConsole.appendText("ACTION: Check individual files or restart prediction\n");
				updateButtonStates(false);
			}
			predictButton.setDisable(false);
		});

		deepDetectTask.setOnFailed(e -> {
			predictionConsole.appendText("\n=== Prediction Failed ===\n");
			predictionConsole.appendText("Timestamp: " + LocalDateTime.now().toString() + "\n");
			predictionConsole.appendText("ERROR: " + deepDetectTask.getException().getMessage() + "\n");
			predictionConsole.appendText("ACTION: Check error details and try again\n");
			updateButtonStates(false);
			predictButton.setDisable(false);
		});

		deepDetectTask.messageProperty().addListener((obs, oldMsg, newMsg) -> {
			if (newMsg != null && !newMsg.isEmpty()) {
				predictionConsole.appendText("Progress: " + newMsg + "\n");
			}
		});

		taskExecutor = Executors.newSingleThreadExecutor();
		taskExecutor.submit(deepDetectTask);

		predictionConsole.appendText("DeepDetect task started...\n");
		predictionConsole.appendText("Monitoring progress...\n\n");
	}

	@Override
	protected void loadLastUsedParameters() {
		if (parametersLoadedFromFile) {
			return;
		}
		try {
			MetaParameterDia params = (MetaParameterDia) parameterIO.loadLastUsedParameters();
			if (params != null) {
				applyParametersToUI(params);
			} else {
				applyParametersToUI(new MetaParameterDia());
			}
		} catch (Exception e) {
			LOGGER.warn("Failed to load last used parameters", e);
		}
	}

	@Override
	protected void applyParametersToUI(MetaParameterDia params) {
		if (params == null) {
			return;
		}

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

		if (params.getResult() != null && ioPanel != null) {
			ioPanel.setOutputDirectory(params.getResult());
		}

		if (threadsComboBox != null) {
			threadsComboBox.setValue(params.getThreadCount());
		}

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

		if (modeSelectionComboBox != null) {
			// Prefer explicit seed mode if available
			String seedMode = params.getSeedMode();
			if ("FOCUS".equalsIgnoreCase(seedMode)) {
				modeSelectionComboBox.setValue("Focus search");
			} else if ("COMPREHENSIVE".equalsIgnoreCase(seedMode)) {
				modeSelectionComboBox.setValue("Comprehensive search");
			} else if ("DEEP".equalsIgnoreCase(seedMode)) {
				modeSelectionComboBox.setValue("Deep comprehensive search");
			} else if ("DIRECT".equalsIgnoreCase(seedMode)) {
				modeSelectionComboBox.setValue("Direct search");
			} else {
				// Final fallback based on librarySearch flag
				modeSelectionComboBox.setValue(params.isLibrarySearch() ? "Direct search" : "Focus search");
			}
		}

		// Modifications (MetaParameterDia and MetaParameterDiaMarker both have getters)
		if (mExcisionCheckBox != null) mExcisionCheckBox.setSelected(params.isMExcision());
		if (carbamiCCheckBox != null) carbamiCCheckBox.setSelected(params.isCarbamiC());
		if (oxidationMCheckBox != null) oxidationMCheckBox.setSelected(params.isOxidationM());
		if (acetyProNCheckBox != null) acetyProNCheckBox.setSelected(params.isAcetyProN());
		if (phosphoSTYCheckBox != null) phosphoSTYCheckBox.setSelected(params.isPhosphoSTY());
		if (kGGCheckBox != null) kGGCheckBox.setSelected(params.isKGG());
		if (usePeptideRoundCheckBox != null) usePeptideRoundCheckBox.setSelected(params.isUsePeptideBasedRound());
		if (detectabilityThresholdField != null) {
			detectabilityThresholdField.setText(String.valueOf(params.getPeptideDetectabilityThreshold()));
		}
		updatePeptideBasedControlsState();

		if (libraryPanel != null) {
			libraryPanel.updateMetaParameter(params);
		}
	}

	@Override
	protected Class<MetaParameterDia> getParameterClass() {
		return MetaParameterDia.class;
	}

	public void saveCurrentParameters() {
		try {
			MetaParameterDia currentParams = collectParametersFromUI();
			if (currentParams == null) {
				currentParams = generateMetaParameterDiaForLibraryGeneration();
			}
			if (currentParams != null) {
				parameterIO.saveLastUsedParameters(currentParams);
			}
		} catch (Exception e) {
			LOGGER.warn("Failed to save current parameters", e);
		}
	}

	@Override
	protected void showError(String title, String content) {
		Alert alert = new Alert(AlertType.ERROR);
		alert.setTitle(title);
		alert.setHeaderText(null);
		alert.setContentText(content);
		alert.showAndWait();
	}
}
