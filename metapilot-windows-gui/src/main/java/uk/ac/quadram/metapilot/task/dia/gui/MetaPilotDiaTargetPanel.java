package uk.ac.quadram.metapilot.task.dia.gui;

import javafx.application.Platform;
import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.control.*;
import javafx.scene.control.Alert.AlertType;
import javafx.scene.layout.*;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.scene.text.Font;
import javafx.scene.text.FontWeight;
import javafx.scene.paint.Color;
import javafx.scene.control.cell.PropertyValueFactory;
import javafx.scene.control.cell.CheckBoxTableCell;
import javafx.scene.control.Hyperlink;
import javafx.beans.property.SimpleStringProperty;
import javafx.beans.property.StringProperty;
import javafx.beans.property.SimpleBooleanProperty;
import javafx.beans.property.BooleanProperty;
import javafx.concurrent.Task;
import uk.ac.quadram.metapilot.task.common.gui.MetaPilotIOPanelFX;
import uk.ac.quadram.metapilot.task.common.gui.MagDatabaseSettingsPanel;
import uk.ac.quadram.metapilot.task.common.gui.MetadataConfigurationPanel;
import uk.ac.quadram.metapilot.task.par.MetaData;
import uk.ac.quadram.metapilot.task.par.MetaBaseParameter;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.task.dda.igc.gui.MetaPilotDdaIgcMaxQuantPanel;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;
import uk.ac.quadram.metapilot.task.dia.par.DiaTargetParameterIO;
import uk.ac.quadram.metapilot.task.dia.par.MetaParameterDiaTg;
import uk.ac.quadram.metapilot.task.common.gui.AbstractWorkflowPanel;
import uk.ac.quadram.metapilot.task.common.gui.ConsolePanel;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.HashSet;
import java.util.Optional;
import java.util.Arrays;
import java.io.File;
import java.net.URI;
import java.awt.Desktop;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.lang.management.ManagementFactory;
import com.sun.management.OperatingSystemMXBean;
import java.text.SimpleDateFormat;
import java.util.Date;
import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNTask;
import uk.ac.quadram.metapilot.dbSearch.diann.DiannParameter;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters.ToolType;
import uk.ac.quadram.metapilot.task.dia.MetaDiaTargetTask;

public class MetaPilotDiaTargetPanel extends AbstractWorkflowPanel<MetaParameterDiaTg> {
	private static final Logger LOGGER = LogManager.getLogger(MetaPilotDiaTargetPanel.class);

	// UI Components - Database specific for MAG (like DIA workflows)
	private MetaPilotIOPanelFX ioPanel;
	private MagDatabaseSettingsPanel magDatabasePanel;
	private DiaTargetDBGenerationPanel dbGenerationPanel;

	// DIA Target-specific components
	private ComboBox<Integer> threadsComboBox;

	// Library Generation Panel Components
	private ComboBox<Float> peptideDetectabilityComboBox;
	private ComboBox<Float> thresholdComboBox; // Alias for peptideDetectabilityComboBox
	private TextField libraryPathTextField;
	private TextField libraryFilePathTextField;
	private Button testThresholdButton;
	private Button createLibraryButton;
	private TextArea libraryConsoleOutput;
	private ProgressBar libraryProgressBar;
	private Label libraryProgressLabel;

	// Taxa Selection Panel Components
	private TableView<TaxaRow> taxaTable;
	private ObservableList<TaxaRow> taxaData = FXCollections.observableArrayList();
	private ComboBox<String> taxaLevelComboBox;
	private Button selectAllTaxaButton;
	private Button deselectAllTaxaButton;

	// Function Selection Panel Components
	private TableView<FunctionRow> functionTable;
	private ObservableList<FunctionRow> functionData = FXCollections.observableArrayList();
	private ComboBox<String> functionLevelComboBox;
	private Button selectAllFunctionButton;
	private Button deselectAllFunctionButton;

	// Metadata configuration components
	private MetadataConfigurationPanel metadataPanel;

	// Parameter I/O

	// Task management (following MetaPilotDiaNormalPanel pattern)
	private MetaDiaTargetTask currentTask;
	private MetaResourceParameters resourceParams;

	// Function types for database queries (same as MetaDiaTargetDBTask.FUNCTIONS)
	private static final String[] FUNCTION_TYPES = new String[] { "GOs", "EC", "KEGG_ko", "KEGG_Pathway", "KEGG_Module",
			"KEGG_Reaction", "KEGG_rclass", "BRITE", "KEGG_TC", "CAZy" };

	// Taxonomic levels for database queries
	private static final String[] TAXONOMIC_LEVELS = new String[] { "genome_name", "superkingdom", "phylum", "class",
			"order_name", "family", "genus", "species" };

	// Date formatter for logging
	private SimpleDateFormat format = new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss");

	public MetaPilotDiaTargetPanel(MagDbItem[] magDbItems, MetaResourceParameters resourceParams) {
		// Call parent constructor first
		super();

		// Initialize parameter IO instance
		parameterIO = new DiaTargetParameterIO();
		this.resourceParams = resourceParams;
		initializeSpecificUI();
		// Set MAG databases after UI is initialized
		setAvailableMagDatabases(magDbItems);
		// Load last used parameters after UI and databases are initialized
		loadLastUsedParameters();
	}

	// Load available MAG databases
	public void setAvailableMagDatabases(MagDbItem[] magDbItems) {
		if (magDatabasePanel != null) {
			magDatabasePanel.setAvailableMagDatabases(magDbItems);
			// Update the database generation panel with the selected MAG database
			if (magDbItems != null && magDbItems.length > 0 && dbGenerationPanel != null) {
				MagDbItem firstDb = magDbItems[0];
				MagDbItem.MagVersion firstVersion = firstDb.getVersionList() != null
						&& !firstDb.getVersionList().isEmpty() ? firstDb.getVersionList().get(0) : null;
				dbGenerationPanel.setMagDbItemWithVersion(firstDb, firstVersion);
				// Also set the resource parameters and parent panel
				dbGenerationPanel.setResourceParameters(resourceParams);
				dbGenerationPanel.setParentPanel(this);
			}
		}
	}

	private void initializeSpecificUI() {
		setSpacing(10);
		setPadding(new Insets(10));
		setMinWidth(1200);
		setMinHeight(800);
		setPrefWidth(1600);
		setPrefHeight(900);

		// Create tab pane with responsive sizing
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

		// Add tab change listener for parameter updates and auto-save
		tabPane.getSelectionModel().selectedItemProperty().addListener((obs, oldTab, newTab) -> {
			if (newTab == paramTab) {
				// Update metadata table when switching to Parameters tab
				updateMetadataTableFromIOPanel();
				// Update parameter view
				updateParameterView();
			} else if (newTab == consoleTab) {
				// Update parameter view and auto-save when switching to Console tab
				updateParameterView();

				// Save current parameters when Console tab is selected
				LOGGER.debug("[DIA Target] Auto-save triggered when switching to Console tab");
				try {
					LOGGER.debug("[DIA Target] Attempting to collect parameters for auto-save...");
					MetaParameterDiaTg currentParams = collectParametersFromUI();
					LOGGER.debug("[DIA Target] collectParametersFromUI() for auto-save returned: {}",
							currentParams != null ? "SUCCESS" : "NULL");

					if (currentParams == null) {
						LOGGER.debug("[DIA Target] Full validation failed, trying fallback for auto-save...");
						// If full validation fails (e.g., no raw files), use minimal validation for
						// auto-save
						currentParams = collectParametersFromUIForAutoSave();
						LOGGER.debug("[DIA Target] collectParametersFromUIForAutoSave() returned: {}",
								currentParams != null ? "SUCCESS" : "NULL");
					}

					if (currentParams != null) {
						LOGGER.debug("[DIA Target] Attempting auto-save...");
						boolean saved = parameterIO.saveLastUsedParameters(currentParams);
						LOGGER.debug("[DIA Target] Auto-save result: {}", saved);
						if (saved) {
							LOGGER.info("[DIA Target] Auto-saved parameters when switching to Console tab");
						} else {
							LOGGER.warn("[DIA Target] Failed to save parameters to file during auto-save");
						}
					} else {
						LOGGER.warn("[DIA Target] No parameters to auto-save - both collection methods returned null");
					}
				} catch (Exception e) {
					LOGGER.warn("Failed to auto-save parameters: " + e.getMessage(), e);
				}
			}
		});

		// Add tab pane to main panel
		getChildren().add(tabPane);

		// Initialize parameter view
		updateParameterView();
	}

	@Override
	protected String getWorkflowName() {
		return "DIA Target";
	}

	@Override
	public VBox createIOTab() {
		VBox tabContent = new VBox(10);
		tabContent.setPadding(new Insets(10));
		ioPanel = new MetaPilotIOPanelFX();
		magDatabasePanel = new MagDatabaseSettingsPanel(false); // false = don't show host database options for DIA
		magDatabasePanel.setParameterUpdateCallback(() -> {
			updateParameterView();
			// Update database generation panel when database selection changes
			if (dbGenerationPanel != null) {
				MagDbItem selectedDb = magDatabasePanel.getSelectedMagDatabase();
				MagDbItem.MagVersion selectedVersion = magDatabasePanel.getSelectedMagVersion();
				dbGenerationPanel.setMagDbItemWithVersion(selectedDb, selectedVersion);
				dbGenerationPanel.setResourceParameters(resourceParams);
				dbGenerationPanel.setParentPanel(this);
			}
		});
		dbGenerationPanel = new DiaTargetDBGenerationPanel();
		dbGenerationPanel.setParameterUpdateCallback(() -> updateParameterView());
		// Set initial database if available
		if (magDatabasePanel != null && magDatabasePanel.getSelectedMagDatabase() != null) {
			MagDbItem selectedDb = magDatabasePanel.getSelectedMagDatabase();
			MagDbItem.MagVersion selectedVersion = magDatabasePanel.getSelectedMagVersion();
			dbGenerationPanel.setMagDbItemWithVersion(selectedDb, selectedVersion);
			dbGenerationPanel.setResourceParameters(resourceParams);
			dbGenerationPanel.setParentPanel(this);
		}

		// Initialize console panel
		consolePanel = new ConsolePanel();
		consolePanel.setStartCallback(this::startAnalysis);
		consolePanel.setStopCallback(this::stopAnalysis);
		consolePanel.setViewReportCallback(this::viewReport);
		tabContent.getChildren().addAll(createSectionTitle("Input/Output Settings"), ioPanel,
				createSectionTitle("Database Settings"), magDatabasePanel,
				createSectionTitle("Taxa-Function Database Generation"), dbGenerationPanel);
		return tabContent;
	}

	@Override
	public VBox createParametersTab() {
		VBox tabContent = new VBox(10);
		tabContent.setPadding(new Insets(10));
		tabContent.setMinWidth(960);
		tabContent.setMinHeight(760);

		// Create main grid layout for parameters with 2-column responsive sizing
		GridPane mainGrid = new GridPane();
		mainGrid.setHgap(15);
		mainGrid.setVgap(15);
		mainGrid.setPadding(new Insets(10));
		mainGrid.setMinWidth(940);
		mainGrid.setMinHeight(740);
		VBox.setVgrow(mainGrid, Priority.ALWAYS);

		// Column constraints for 2-column layout
		ColumnConstraints col1 = new ColumnConstraints();
		col1.setPercentWidth(50);
		col1.setMinWidth(450);
		col1.setHgrow(Priority.ALWAYS);
		ColumnConstraints col2 = new ColumnConstraints();
		col2.setPercentWidth(50);
		col2.setMinWidth(450);
		col2.setHgrow(Priority.ALWAYS);
		mainGrid.getColumnConstraints().addAll(col1, col2);

		// Row constraints for better layout control
		RowConstraints row1 = new RowConstraints();
		row1.setMinHeight(60);
		row1.setVgrow(Priority.NEVER);
		RowConstraints row2 = new RowConstraints();
		row2.setVgrow(Priority.ALWAYS);
		row2.setMinHeight(300);
		RowConstraints row3 = new RowConstraints();
		row3.setVgrow(Priority.ALWAYS);
		row3.setMinHeight(200);
		mainGrid.getRowConstraints().addAll(row1, row2, row3);

		// Thread Count Settings (spans both columns)
		VBox threadSection = createThreadCountSection();
		mainGrid.add(threadSection, 0, 0, 2, 1);

		// Row 2: Taxa Selection Panel (left) | Function Selection Panel (right)
		VBox taxaSection = createTaxaSelectionPanel();
		mainGrid.add(taxaSection, 0, 1);
		GridPane.setVgrow(taxaSection, Priority.ALWAYS);

		VBox functionSection = createFunctionSelectionPanel();
		mainGrid.add(functionSection, 1, 1);
		GridPane.setVgrow(functionSection, Priority.ALWAYS);

		// Row 3: Library Generation Panel (left) | Metadata Configuration (right)
		VBox librarySection = createLibraryGenerationPanel();
		mainGrid.add(librarySection, 0, 2);
		GridPane.setVgrow(librarySection, Priority.ALWAYS);

		// Create reusable metadata configuration panel
		metadataPanel = new MetadataConfigurationPanel();
		mainGrid.add(metadataPanel, 1, 2);
		GridPane.setVgrow(metadataPanel, Priority.ALWAYS);

		tabContent.getChildren().add(mainGrid);
		return tabContent;
	}

	private VBox createThreadCountSection() {
		VBox section = new VBox(10);
		section.setPadding(new Insets(10));
		section.setStyle("-fx-background-color: #f0f8ff; -fx-background-radius: 5;");

		Label titleLabel = createSectionTitle("Performance Settings");

		HBox threadBox = new HBox(15);
		threadBox.setAlignment(Pos.CENTER_LEFT);

		Label threadsLabel = new Label("Number of Threads:");
		threadsLabel.setStyle("-fx-font-weight: bold;");
		threadsComboBox = new ComboBox<>();
		for (int i = 1; i <= Runtime.getRuntime().availableProcessors(); i++) {
			threadsComboBox.getItems().add(i);
		}
		threadsComboBox.setValue(Math.min(4, Runtime.getRuntime().availableProcessors()));
		threadsComboBox.setMinWidth(100);
		threadsComboBox.setPrefWidth(120);
		threadsComboBox.setOnAction(e -> updateParameterView());

		threadBox.getChildren().addAll(threadsLabel, threadsComboBox);
		section.getChildren().addAll(titleLabel, threadBox);
		return section;
	}

	private VBox createLibraryGenerationPanel() {
		VBox section = new VBox(10);
		section.setPadding(new Insets(10));
		section.setStyle("-fx-background-color: #ffefd5; -fx-background-radius: 5;");
		section.setMinHeight(180);
		VBox.setVgrow(section, Priority.ALWAYS);

		Label titleLabel = createSectionTitle("Library Generation");

		// Threshold selection and buttons in one row
		HBox thresholdAndButtonsBox = new HBox(15);
		thresholdAndButtonsBox.setAlignment(Pos.CENTER_LEFT);

		// Threshold selection
		HBox thresholdBox = new HBox(10);
		thresholdBox.setAlignment(Pos.CENTER_LEFT);
		Label thresholdLabel = new Label("Detectability Threshold:");
		thresholdLabel.setStyle("-fx-font-weight: bold;");
		peptideDetectabilityComboBox = new ComboBox<>();
		for (int i = 0; i <= 100; i++) {
			float threshold = (float) (i / 100.0);
			peptideDetectabilityComboBox.getItems().add(threshold);
		}
		peptideDetectabilityComboBox.setValue(0.50f); // Default to 0.5

		// Initialize thresholdComboBox to point to the same object
		thresholdComboBox = peptideDetectabilityComboBox;
		peptideDetectabilityComboBox.setMinWidth(80);
		peptideDetectabilityComboBox.setPrefWidth(100);
		peptideDetectabilityComboBox.setOnAction(e -> updateParameterView());
		thresholdBox.getChildren().addAll(thresholdLabel, peptideDetectabilityComboBox);

		// Action buttons
		testThresholdButton = new Button("Get Proper Threshold");
		testThresholdButton.setPrefWidth(150);
		testThresholdButton.setOnAction(e -> {
			LOGGER.debug("[DIA Target] Get Proper Threshold button clicked");
			calculateProperThreshold();
		});

		createLibraryButton = new Button("Create Library");
		createLibraryButton.setPrefWidth(120);
		createLibraryButton.setOnAction(e -> createLibrary());

		// Add threshold and buttons to the same row
		thresholdAndButtonsBox.getChildren().addAll(thresholdBox, testThresholdButton, createLibraryButton);

		// Library name input
		HBox libraryNameBox = new HBox(10);
		libraryNameBox.setAlignment(Pos.CENTER_LEFT);
		Label libraryNameLabel = new Label("Library Name:");
		libraryNameLabel.setStyle("-fx-font-weight: bold;");

		libraryPathTextField = new TextField();
		libraryPathTextField.setPromptText("Enter library name (e.g., my_target_library)");
		libraryPathTextField.setText("target_library");
		libraryPathTextField.setPrefWidth(200);
		libraryPathTextField.setOnAction(e -> updateParameterView());
		libraryNameBox.getChildren().addAll(libraryNameLabel, libraryPathTextField);

		// Library file path display
		HBox libraryPathBox = new HBox(10);
		libraryPathBox.setAlignment(Pos.CENTER_LEFT);
		Label libraryPathLabel = new Label("Library File Path:");
		libraryPathLabel.setStyle("-fx-font-weight: bold;");

		libraryFilePathTextField = new TextField();
		libraryFilePathTextField.setPromptText("Library file path will be shown here after creation");
		libraryFilePathTextField.setPrefWidth(400);
		libraryFilePathTextField.setEditable(false); // Read-only field
		libraryFilePathTextField.setStyle("-fx-background-color: #f0f0f0;"); // Gray background to indicate read-only
		libraryPathBox.getChildren().addAll(libraryPathLabel, libraryFilePathTextField);

		// Progress indicator section
		VBox progressSection = new VBox(8);
		progressSection.setPadding(new Insets(10));
		progressSection.setStyle("-fx-background-color: #f8f8f8; -fx-background-radius: 5;");
		progressSection.setMinHeight(80);

		libraryProgressLabel = new Label("Ready to generate library");
		libraryProgressLabel.setStyle("-fx-font-weight: bold; -fx-font-size: 12px; -fx-text-fill: #333333;");
		libraryProgressLabel.setMaxWidth(Double.MAX_VALUE);

		libraryProgressBar = new ProgressBar(0);
		libraryProgressBar.setMaxWidth(Double.MAX_VALUE);
		libraryProgressBar.setPrefHeight(20);
		libraryProgressBar.setStyle("-fx-accent: #4CAF50;");
		libraryProgressBar.setVisible(false); // Initially hidden

		progressSection.getChildren().addAll(libraryProgressLabel, libraryProgressBar);

		// Console output
		Label consoleLabel = new Label("Progress Output:");
		consoleLabel.setStyle("-fx-font-weight: bold; -fx-font-size: 12px;");
		libraryConsoleOutput = new TextArea();
		libraryConsoleOutput.setEditable(false);
		libraryConsoleOutput.setPrefHeight(200);
		libraryConsoleOutput.setMinHeight(150);
		libraryConsoleOutput.setStyle("-fx-font-family: monospace; -fx-font-size: 12px;");
		VBox.setVgrow(libraryConsoleOutput, Priority.ALWAYS);

		section.getChildren().addAll(titleLabel, thresholdAndButtonsBox, libraryNameBox, libraryPathBox,
				progressSection, consoleLabel, libraryConsoleOutput);
		return section;
	}

	private VBox createTaxaSelectionPanel() {
		VBox section = new VBox(10);
		section.setPadding(new Insets(10));
		section.setStyle("-fx-background-color: #fff8dc; -fx-background-radius: 5;");
		section.setMinHeight(250);
		VBox.setVgrow(section, Priority.ALWAYS);

		Label titleLabel = createSectionTitle("Taxa Selection");

		// Taxa level selection
		HBox levelBox = new HBox(10);
		levelBox.setAlignment(Pos.CENTER_LEFT);
		Label levelLabel = new Label("Taxonomic Level:");
		levelLabel.setStyle("-fx-font-weight: bold;");
		taxaLevelComboBox = new ComboBox<>(FXCollections.observableArrayList("Choose the taxonomic level",
				"Superkingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Genome"));
		taxaLevelComboBox.setValue("Choose the taxonomic level");
		taxaLevelComboBox.setMinWidth(200);
		taxaLevelComboBox.setPrefWidth(300);
		taxaLevelComboBox.setOnAction(e -> loadTaxaData());
		levelBox.getChildren().addAll(levelLabel, taxaLevelComboBox);

		// Taxa table
		taxaTable = new TableView<>(taxaData);
		taxaTable.setMinHeight(150);
		taxaTable.setPrefHeight(200);
		taxaTable.setMaxHeight(Double.MAX_VALUE);
		taxaTable.setEditable(true);
		taxaTable.setColumnResizePolicy(TableView.CONSTRAINED_RESIZE_POLICY_FLEX_LAST_COLUMN);
		VBox.setVgrow(taxaTable, Priority.ALWAYS);

		// Select column (checkbox)
		TableColumn<TaxaRow, Boolean> selectCol = new TableColumn<>("Select");
		selectCol.setCellValueFactory(new PropertyValueFactory<>("selected"));
		selectCol.setCellFactory(CheckBoxTableCell.forTableColumn(selectCol));
		selectCol.setOnEditCommit(e -> updateParameterView());
		selectCol.setMinWidth(60);
		selectCol.setPrefWidth(80);
		selectCol.setSortable(false);

		// Taxa Name column with clickable links
		TableColumn<TaxaRow, String> nameCol = new TableColumn<>("Taxonomic Name");
		nameCol.setCellValueFactory(new PropertyValueFactory<>("name"));
		nameCol.setCellFactory(column -> {
			return new TableCell<TaxaRow, String>() {
				private final Hyperlink hyperlink = new Hyperlink();

				@Override
				protected void updateItem(String item, boolean empty) {
					super.updateItem(item, empty);

					if (empty || item == null) {
						setGraphic(null);
						setText(null);
					} else {
						hyperlink.setText(item);
						hyperlink.setOnAction(e -> {
							// Only create genome links for genome level
							String selectedLevel = taxaLevelComboBox.getValue();
							if ("Genome".equals(selectedLevel)) {
								String url = getGenomeURL(item);
								if (url != null) {
									openURL(url);
								}
							}
						});
						setGraphic(hyperlink);
						setText(null);
					}
				}
			};
		});
		nameCol.setMinWidth(200);
		nameCol.setPrefWidth(300);
		nameCol.setSortable(true);

		taxaTable.getColumns().clear();
		taxaTable.getColumns().add(selectCol);
		taxaTable.getColumns().add(nameCol);

		// Control buttons
		HBox controlBox = new HBox(10);
		controlBox.setAlignment(Pos.CENTER_LEFT);

		selectAllTaxaButton = new Button("Select All");
		selectAllTaxaButton.setOnAction(e -> selectAllTaxa(true));

		deselectAllTaxaButton = new Button("Deselect All");
		deselectAllTaxaButton.setOnAction(e -> selectAllTaxa(false));

		controlBox.getChildren().addAll(selectAllTaxaButton, deselectAllTaxaButton);

		section.getChildren().addAll(titleLabel, levelBox, taxaTable, controlBox);
		return section;
	}

	private VBox createFunctionSelectionPanel() {
		VBox section = new VBox(10);
		section.setPadding(new Insets(10));
		section.setStyle("-fx-background-color: #f0fff0; -fx-background-radius: 5;");
		section.setMinHeight(250);
		VBox.setVgrow(section, Priority.ALWAYS);

		Label titleLabel = createSectionTitle("Function Selection");

		// Function level selection
		HBox levelBox = new HBox(10);
		levelBox.setAlignment(Pos.CENTER_LEFT);
		Label levelLabel = new Label("Functional Level:");
		levelLabel.setStyle("-fx-font-weight: bold;");
		ObservableList<String> functionTypes = FXCollections.observableArrayList("Choose the functional level");
		functionTypes.addAll(FXCollections.observableArrayList(FUNCTION_TYPES));
		functionLevelComboBox = new ComboBox<>(functionTypes);
		functionLevelComboBox.setValue("Choose the functional level");
		functionLevelComboBox.setMinWidth(200);
		functionLevelComboBox.setPrefWidth(300);
		functionLevelComboBox.setOnAction(e -> loadFunctionData());
		levelBox.getChildren().addAll(levelLabel, functionLevelComboBox);

		// Function table
		functionTable = new TableView<>(functionData);
		functionTable.setMinHeight(150);
		functionTable.setPrefHeight(200);
		functionTable.setMaxHeight(Double.MAX_VALUE);
		functionTable.setEditable(true);
		functionTable.setColumnResizePolicy(TableView.CONSTRAINED_RESIZE_POLICY_FLEX_LAST_COLUMN);
		VBox.setVgrow(functionTable, Priority.ALWAYS);

		// Select column (checkbox)
		TableColumn<FunctionRow, Boolean> selectCol = new TableColumn<>("Select");
		selectCol.setCellValueFactory(new PropertyValueFactory<>("selected"));
		selectCol.setCellFactory(CheckBoxTableCell.forTableColumn(selectCol));
		selectCol.setOnEditCommit(e -> updateParameterView());
		selectCol.setMinWidth(60);
		selectCol.setPrefWidth(80);
		selectCol.setSortable(false);

		// Function Name column with clickable links
		TableColumn<FunctionRow, String> nameCol = new TableColumn<>("Function Name");
		nameCol.setCellValueFactory(new PropertyValueFactory<>("name"));
		nameCol.setCellFactory(column -> {
			return new TableCell<FunctionRow, String>() {
				private final Hyperlink hyperlink = new Hyperlink();

				@Override
				protected void updateItem(String item, boolean empty) {
					super.updateItem(item, empty);

					if (empty || item == null) {
						setGraphic(null);
						setText(null);
					} else {
						hyperlink.setText(item);
						hyperlink.setOnAction(e -> {
							String selectedFunctionType = functionLevelComboBox.getValue();
							if (selectedFunctionType != null
									&& !selectedFunctionType.equals("Choose the functional level")) {
								String url = getFunctionURL(selectedFunctionType, item);
								if (url != null) {
									openURL(url);
								}
							}
						});
						setGraphic(hyperlink);
						setText(null);
					}
				}
			};
		});
		nameCol.setMinWidth(200);
		nameCol.setPrefWidth(300);
		nameCol.setSortable(true);

		functionTable.getColumns().clear();
		functionTable.getColumns().add(selectCol);
		functionTable.getColumns().add(nameCol);

		// Control buttons
		HBox controlBox = new HBox(10);
		controlBox.setAlignment(Pos.CENTER_LEFT);

		selectAllFunctionButton = new Button("Select All");
		selectAllFunctionButton.setOnAction(e -> selectAllFunctions(true));

		deselectAllFunctionButton = new Button("Deselect All");
		deselectAllFunctionButton.setOnAction(e -> selectAllFunctions(false));

		controlBox.getChildren().addAll(selectAllFunctionButton, deselectAllFunctionButton);

		section.getChildren().addAll(titleLabel, levelBox, functionTable, controlBox);
		return section;
	}

	// Taxa Row class
	public static class TaxaRow {
		private final SimpleBooleanProperty selected;
		private final SimpleStringProperty name;

		public TaxaRow(boolean selected, String name) {
			this.selected = new SimpleBooleanProperty(selected);
			this.name = new SimpleStringProperty(name);
		}

		public boolean getSelected() {
			return selected.get();
		}

		public boolean isSelected() {
			return selected.get();
		} // Alias for getSelected()

		public void setSelected(boolean value) {
			selected.set(value);
		}

		public BooleanProperty selectedProperty() {
			return selected;
		}

		public String getName() {
			return name.get();
		}

		public void setName(String value) {
			name.set(value);
		}

		public StringProperty nameProperty() {
			return name;
		}
	}

	// Function Row class
	public static class FunctionRow {
		private final SimpleBooleanProperty selected;
		private final SimpleStringProperty name;

		public FunctionRow(boolean selected, String name) {
			this.selected = new SimpleBooleanProperty(selected);
			this.name = new SimpleStringProperty(name);
		}

		public boolean getSelected() {
			return selected.get();
		}

		public boolean isSelected() {
			return selected.get();
		} // Alias for getSelected()

		public void setSelected(boolean value) {
			selected.set(value);
		}

		public BooleanProperty selectedProperty() {
			return selected;
		}

		public String getName() {
			return name.get();
		}

		public void setName(String value) {
			name.set(value);
		}

		public StringProperty nameProperty() {
			return name;
		}
	}

	private void loadTaxaData() {
		String selectedLevel = taxaLevelComboBox.getValue();
		taxaData.clear();

		if (selectedLevel == null || selectedLevel.equals("Choose the taxonomic level")) {
			return;
		}

		// Check if database generation panel has the database available
		if (dbGenerationPanel == null || dbGenerationPanel.getCurrentMagDbItem() == null) {
			System.err.println("DEBUG: No MAG database available for taxa data loading");
			return;
		}

		MagDbItem magDbItem = dbGenerationPanel.getCurrentMagDbItem();
		File currentFile = magDbItem.getCurrentFile();
		if (currentFile == null) {
			System.err.println("DEBUG: No current file available for MAG database");
			return;
		}

		File sqlDbFile = new File(currentFile, "pro_func.db");
		if (!sqlDbFile.exists()) {
			System.err.println("DEBUG: pro_func.db not found at: " + sqlDbFile.getAbsolutePath());
			return;
		}

		// Map UI level names to database column names
		String[] levelNames = { "Superkingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Genome" };
		int levelIndex = -1;
		for (int i = 0; i < levelNames.length; i++) {
			if (levelNames[i].equals(selectedLevel)) {
				levelIndex = i;
				break;
			}
		}

		if (levelIndex == -1) {
			System.err.println("DEBUG: Invalid taxonomic level: " + selectedLevel);
			return;
		}

		// Query the database
		String url = "jdbc:sqlite:" + sqlDbFile.getAbsolutePath().replaceAll("\\\\", "/");
		HashSet<String> taxaSet = new HashSet<>();

		try (Connection conn = DriverManager.getConnection(url)) {
			String columnName;
			String query;

			if (levelIndex < TAXONOMIC_LEVELS.length - 1) {
				// Map UI level to database column (TAXONOMIC_LEVELS starts with genome_name, so
				// we need to offset by 1)
				columnName = TAXONOMIC_LEVELS[levelIndex + 1];
				query = "SELECT DISTINCT " + columnName + " FROM taxonomy WHERE " + columnName + " IS NOT NULL AND "
						+ columnName + " != ''";
			} else {
				// Genome level (last item in levelNames array)
				columnName = "genome_name";
				query = "SELECT DISTINCT genome_name FROM taxonomy";
			}

			try (PreparedStatement pstmt = conn.prepareStatement(query); ResultSet rs = pstmt.executeQuery()) {
				while (rs.next()) {
					String taxon = rs.getString(columnName);
					if (taxon != null && !taxon.trim().isEmpty()) {
						taxaSet.add(taxon);
					}
				}
			}
		} catch (SQLException e) {
			System.err.println("ERROR: Failed to query taxa data from database: " + e.getMessage());
			e.printStackTrace();
			return;
		}

		// Convert to sorted list and add to taxaData
		String[] sortedTaxa = taxaSet.toArray(new String[0]);
		Arrays.sort(sortedTaxa);

		for (String taxon : sortedTaxa) {
			taxaData.add(new TaxaRow(false, taxon));
		}

		LOGGER.debug("[DIA Target] Loaded {} taxa from database for level: {}", taxaData.size(), selectedLevel);
		updateParameterView();
	}

	private void loadFunctionData() {
		String selectedLevel = functionLevelComboBox.getValue();
		functionData.clear();

		if (selectedLevel == null || selectedLevel.equals("Choose the functional level")) {
			return;
		}

		// Check if database generation panel has the database available
		if (dbGenerationPanel == null || dbGenerationPanel.getCurrentMagDbItem() == null) {
			System.err.println("DEBUG: No MAG database available for function data loading");
			return;
		}

		MagDbItem magDbItem = dbGenerationPanel.getCurrentMagDbItem();
		File currentFile = magDbItem.getCurrentFile();
		if (currentFile == null) {
			System.err.println("DEBUG: No current file available for MAG database");
			return;
		}

		File sqlDbFile = new File(currentFile, "pro_func.db");
		if (!sqlDbFile.exists()) {
			System.err.println("DEBUG: pro_func.db not found at: " + sqlDbFile.getAbsolutePath());
			return;
		}

		// Query the database for the selected function type
		String url = "jdbc:sqlite:" + sqlDbFile.getAbsolutePath().replaceAll("\\\\", "/");
		HashSet<String> functionSet = new HashSet<>();

		try (Connection conn = DriverManager.getConnection(url)) {
			String query = "SELECT DISTINCT func_name FROM " + selectedLevel
					+ " WHERE func_name IS NOT NULL AND func_name != ''";

			try (PreparedStatement pstmt = conn.prepareStatement(query); ResultSet rs = pstmt.executeQuery()) {
				while (rs.next()) {
					String function = rs.getString("func_name");
					if (function != null && !function.trim().isEmpty()) {
						functionSet.add(function);
					}
				}
			}
		} catch (SQLException e) {
			System.err.println("ERROR: Failed to query function data from database: " + e.getMessage());
			e.printStackTrace();
			return;
		}

		// Convert to sorted list and add to functionData
		String[] sortedFunctions = functionSet.toArray(new String[0]);
		Arrays.sort(sortedFunctions);

		for (String function : sortedFunctions) {
			functionData.add(new FunctionRow(false, function));
		}

		LOGGER.debug("[DIA Target] Loaded {} functions from database for type: {}", functionData.size(), selectedLevel);
		updateParameterView();
	}

	/**
	 * Generate URL for function based on function type and name Based on
	 * MetaPilotMainPanelDiaTg.getFunURL method
	 */
	private String getFunctionURL(String functionType, String functionName) {
		switch (functionType) {
		case "GOs":
			return "https://www.ebi.ac.uk/QuickGO/term/" + functionName;
		case "EC":
			return "https://www.genome.jp/dbget-bin/www_bget?ec:" + functionName;
		case "KEGG_ko":
			return "https://www.genome.jp/dbget-bin/www_bget?" + functionName;
		case "KEGG_Pathway":
			return "https://www.genome.jp/dbget-bin/www_bget?path:" + functionName;
		case "KEGG_Module":
			return "https://www.genome.jp/module/ath_" + functionName;
		case "KEGG_Reaction":
			return "https://www.genome.jp/dbget-bin/www_bget?rn:" + functionName;
		case "KEGG_rclass":
			return "https://www.genome.jp/entry/" + functionName;
		case "BRITE":
			return "https://www.genome.jp/entry/" + functionName;
		case "KEGG_TC":
			return "https://www.genome.jp/entry/" + functionName;
		case "CAZy":
			return "https://www.cazy.org/" + functionName + ".html";
		default:
			return null;
		}
	}

	/**
	 * Generate URL for genome based on genome name Based on MetaPilotMainPanelDiaTg
	 * genome URL pattern
	 */
	private String getGenomeURL(String genomeName) {
		return "https://www.ebi.ac.uk/metagenomics/genomes/" + genomeName + "#overview";
	}

	/**
	 * Open URL in default browser
	 */
	private void openURL(String url) {
		try {
			if (Desktop.isDesktopSupported()) {
				Desktop.getDesktop().browse(new URI(url));
			} else {
				System.err.println("Desktop not supported, cannot open URL: " + url);
			}
		} catch (Exception e) {
			LOGGER.error("Failed to open URL: {} - {}", url, e.getMessage(), e);
		}
	}

	private void selectAllTaxa(boolean select) {
		for (TaxaRow row : taxaData) {
			row.setSelected(select);
		}
		taxaTable.refresh();
		updateParameterView();
	}

	private void selectAllFunctions(boolean select) {
		for (FunctionRow row : functionData) {
			row.setSelected(select);
		}
		functionTable.refresh();
		updateParameterView();
	}

	private Label createSectionTitle(String title) {
		Label label = new Label(title);
		label.setFont(Font.font("System", FontWeight.BOLD, 14));
		label.setTextFill(Color.DARKBLUE);
		return label;
	}

	/**
	 * Show progress indicator for library operations
	 */
	private void showLibraryProgress(String message) {
		if (libraryProgressBar != null && libraryProgressLabel != null) {
			libraryProgressBar.setVisible(true);
			libraryProgressBar.setProgress(-1); // Indeterminate progress (spinning)
			libraryProgressLabel.setText(message);
			libraryProgressLabel.setStyle("-fx-font-weight: bold; -fx-font-size: 12px; -fx-text-fill: #0066CC;");
		}
	}

	/**
	 * Hide progress indicator and show completion message
	 */
	private void hideLibraryProgress(String message) {
		if (libraryProgressBar != null && libraryProgressLabel != null) {
			libraryProgressBar.setVisible(false);
			libraryProgressLabel.setText(message);
			libraryProgressLabel.setStyle("-fx-font-weight: bold; -fx-font-size: 12px; -fx-text-fill: #009900;");
		}
	}

	/**
	 * Reset progress indicator to ready state
	 */
	private void resetLibraryProgress() {
		if (libraryProgressBar != null && libraryProgressLabel != null) {
			libraryProgressBar.setVisible(false);
			libraryProgressLabel.setText("Ready to generate library");
			libraryProgressLabel.setStyle("-fx-font-weight: bold; -fx-font-size: 12px; -fx-text-fill: #333333;");
		}
	}

	private void updateMetadataTableFromIOPanel() {
		if (ioPanel == null || metadataPanel == null) {
			return;
		}

		ObservableList<MetaPilotIOPanelFX.ExperimentRow> experiments = ioPanel.getExperimentData();
		// Update the reusable metadata panel with experiment data
		metadataPanel.updateMetadataFromIOPanel(experiments);
	}

	@Override
	protected void updateParameterView() {
		if (consolePanel == null) {
			return;
		}

		StringBuilder sb = new StringBuilder();
		sb.append("=== DIA Target Analysis Parameters ===\n\n");

		// Database settings
		sb.append("Database Settings:\n");
		if (magDatabasePanel != null) {
			MagDbItem selectedDb = magDatabasePanel.getSelectedMagDatabase();
			if (selectedDb != null) {
				sb.append("  MAG Database: ").append(selectedDb.getCatalogueID()).append("\n");
			} else {
				sb.append("  MAG Database: Not selected\n");
			}
			MagDbItem.MagVersion selectedVersion = magDatabasePanel.getSelectedMagVersion();
			if (selectedVersion != null) {
				sb.append("  MAG Version: ").append(selectedVersion.getCatalogueVersion()).append("\n");
			}
		}
		sb.append("\n");

		// Performance settings
		sb.append("Performance Settings:\n");
		if (threadsComboBox != null) {
			sb.append("  Thread Count: ").append(threadsComboBox.getValue()).append("\n");
		}
		sb.append("\n");

		// Library generation settings
		sb.append("Library Generation:\n");
		if (peptideDetectabilityComboBox != null) {
			sb.append("  Peptide Detectability Threshold: ").append(peptideDetectabilityComboBox.getValue())
					.append("\n");
		}
		if (libraryFilePathTextField != null && !libraryFilePathTextField.getText().isEmpty()) {
			sb.append("  Library Path: ").append(libraryFilePathTextField.getText()).append("\n");
		}
		sb.append("\n");

		// Taxa selection
		sb.append("Taxa Selection:\n");
		if (taxaLevelComboBox != null) {
			sb.append("  Taxonomic Level: ").append(taxaLevelComboBox.getValue()).append("\n");
		}
		if (taxaData != null) {
			int selectedCount = 0;
			for (TaxaRow row : taxaData) {
				if (row.getSelected()) {
					selectedCount++;
				}
			}
			sb.append("  Selected Taxa: ").append(selectedCount).append(" of ").append(taxaData.size()).append("\n");
		}
		sb.append("\n");

		// Function selection
		sb.append("Function Selection:\n");
		if (functionLevelComboBox != null) {
			sb.append("  Functional Level: ").append(functionLevelComboBox.getValue()).append("\n");
		}
		if (functionData != null) {
			int selectedCount = 0;
			for (FunctionRow row : functionData) {
				if (row.getSelected()) {
					selectedCount++;
				}
			}
			sb.append("  Selected Functions: ").append(selectedCount).append(" of ").append(functionData.size())
					.append("\n");
		}
		sb.append("\n");

		consolePanel.setParameterText(sb.toString());
	}

	@Override
	protected void startAnalysis() {
		try {
			consolePanel.appendConsoleText("Starting DIA Target analysis...\n");
			LOGGER.info("startAnalysis called - DIA Target workflow");

			// 1. Validate inputs (custom validation for DIA Target)
			String resultDir = ioPanel.getOutputDirectory();
			if (resultDir == null || resultDir.trim().isEmpty()) {
				showError("Configuration Error",
						"Output directory not set. Please set output directory in the Input/Output tab.");
				return;
			}

			// Check if library file path is set
			String libraryFilePath = libraryFilePathTextField != null ? libraryFilePathTextField.getText().trim() : "";
			if (libraryFilePath.isEmpty()) {
				showError("Library Error",
						"No library file path found. Please create a target library first using the 'Create Library' button.");
				return;
			}

			File finalLibFile = new File(libraryFilePath);
			if (!finalLibFile.exists()) {
				showError("Library Error", "Target library not found at: " + libraryFilePath
						+ "\nPlease create a target library first using the 'Create Library' button.");
				return;
			}

			// 2. Collect parameters from UI
			consolePanel.appendConsoleText("Collecting parameters from UI...\n");
			MetaParameterDiaTg targetParams = collectParametersFromUI();
			if (targetParams == null) {
				consolePanel.appendConsoleText("ERROR: Failed to collect parameters from UI.\n");
				LOGGER.error("Failed to collect parameters from UI in startAnalysis");
				showError("Parameter Error", "Failed to collect parameters from UI. Please check your settings.");
				return;
			}

			consolePanel.appendConsoleText("✓ Successfully collected parameters from UI\n");
			LOGGER.info("Successfully collected parameters from UI for DIA Target workflow");

			// 3. Store analysis parameters for viewReport
			currentAnalysisParams = targetParams;

			targetParams.setLibraryPath(libraryFilePath);

			// Set the result directory to the library directory for DIA Target workflow
			targetParams.setDbSearchResultFile(finalLibFile.getParentFile());

			// 4. Save parameters to file
			String outputDir = targetParams.getResult();
			if (outputDir != null && !outputDir.trim().isEmpty()) {
				consolePanel.appendConsoleText("Saving parameters to output directory...\n");
				
				// Save analysis-specific parameters
				boolean analysisParamsSaved = parameterIO.saveAnalysisParameters(outputDir, targetParams);
				if (!analysisParamsSaved) {
					consolePanel.appendConsoleText("WARNING: Failed to save analysis parameters to result directory\n");
					LOGGER.warn("Failed to save analysis parameters to {}/analysis_parameters.json", outputDir);
				} else {
					consolePanel.appendConsoleText("✓ Analysis parameters saved to result directory\n");
					LOGGER.info("Saved analysis parameters to {}/analysis_parameters.json", outputDir);
				}
				
				// Save as last used parameters
				boolean lastParamsSaved = parameterIO.saveLastUsedParameters(targetParams);
				if (!lastParamsSaved) {
					consolePanel.appendConsoleText("WARNING: Failed to save last used parameters\n");
					LOGGER.warn("Failed to save last used parameters");
				} else {
					consolePanel.appendConsoleText("✓ Last used parameters saved\n");
				}
				
				// Also save to preferences for quick defaults
				boolean prefsSaved = parameterIO.savePreferences(targetParams);
				if (prefsSaved) {
					consolePanel.appendConsoleText("✓ Preferences saved\n");
				} else {
					consolePanel.appendConsoleText("WARNING: Failed to save preferences\n");
				}
			} else {
				consolePanel.appendConsoleText("WARNING: No output directory specified, skipping parameter save\n");
				LOGGER.warn("No output directory specified, skipping parameter save");
			}

			// 5. Create MetaResourceParameters (load from config file)
			consolePanel.appendConsoleText("Loading resource parameters from config...\n");
			consolePanel.appendConsoleText("✓ Resource parameters loaded from config file\n");

			// 6. Create progress bars and status label
			consolePanel.appendConsoleText("Setting up progress tracking...\n");
			ProgressBar overallBar = new ProgressBar();
			ProgressBar currentBar = new ProgressBar();
			Label statusLabel = new Label("Ready to start analysis...");
			consolePanel.appendConsoleText("✓ Progress bars and status label initialized\n");

			// 7. Status label initialized (not stored in DIA Target parameters)
			consolePanel.appendConsoleText("✓ Status label initialized\n");

			consolePanel.appendConsoleText("Library search result folder: " + targetParams.getDbSearchResultFile());

			// Update UI
			consolePanel.clearConsole();
			consolePanel.appendConsoleText("Using target library: " + finalLibFile.getAbsolutePath() + "\n");
			consolePanel.appendConsoleText("Library search mode: ENABLED\n");

			// Reset progress bars
			consolePanel.getOverallProgressBar().setProgress(0);
			consolePanel.getCurrentTaskProgressBar().setProgress(0);
			consolePanel.setCurrentTaskLabel("Current Task: Initializing...");

			// Read target library info file to get genome and protein filtering information
			TargetLibraryInfo libraryInfo = readTargetLibraryInfo(finalLibFile.getParentFile(),
					libraryPathTextField != null ? libraryPathTextField.getText().trim() : "target_library");

			if (libraryInfo == null) {
				showError("Library Error", "Failed to read target library information. Please recreate the library.");
				return;
			}

			// 8. Create and configure task
			consolePanel.appendConsoleText("Creating DIA Target task...\n");
			currentTask = new MetaDiaTargetTask(targetParams, resourceParams,
					consolePanel.getOverallProgressBar(), consolePanel.getCurrentTaskProgressBar());
			consolePanel.appendConsoleText("✓ Task created successfully\n");

			// Set target library information for peptide mapping
			currentTask.setTargetLibraryInfo(libraryInfo);

			// 9. Set up task event handlers
			consolePanel.appendConsoleText("Setting up task callbacks...\n");
			currentTask.setOnSucceeded(e -> {
				Platform.runLater(() -> {
					consolePanel.appendConsoleText("✓ Analysis completed successfully!\n");
					consolePanel.appendConsoleText("Results saved to: " + resultDir + "\n");
					consolePanel.enableViewReportButton();
					consolePanel.setButtonsEnabled(true, false, true);
					consolePanel.stopStatusTimer(); // Stop the timer when task completes
				});
			});

			currentTask.setOnFailed(e -> {
				Platform.runLater(() -> {
					Throwable exception = currentTask.getException();
					String errorMessage = exception != null ? exception.getMessage() : "Unknown error occurred";
					consolePanel.appendConsoleText("✗ Analysis failed: " + errorMessage + "\n");
					consolePanel.setButtonsEnabled(true, false, false);
					consolePanel.stopStatusTimer(); // Stop the timer when task fails
					showError("Analysis Failed", "DIA Target analysis failed: " + errorMessage);
				});
			});

			currentTask.setOnCancelled(e -> {
				Platform.runLater(() -> {
					consolePanel.appendConsoleText("DIA Target analysis was cancelled.\n");
					consolePanel.setButtonsEnabled(true, false, false);
					consolePanel.stopStatusTimer(); // Stop the timer when task is cancelled
				});
			});

			consolePanel.appendConsoleText("✓ Task callbacks configured\n");

			// 10. Run task in background thread
			consolePanel.appendConsoleText("Starting task execution...\n");
			Thread taskThread = new Thread(currentTask);
			taskThread.setDaemon(true);
			taskThread.start();

			consolePanel.appendConsoleText("✓ Analysis task started successfully.\n");
			LOGGER.info("DIA Target analysis task started successfully");
			LOGGER.info("DIA Target analysis started with library: " + finalLibFile.getAbsolutePath());

		} catch (Exception e) {
			String errorMsg = "ERROR: Failed to start analysis: " + e.getMessage();
			consolePanel.appendConsoleText(errorMsg + "\n");
			LOGGER.error("Failed to start DIA Target analysis", e);
			showError("Start Error", "Failed to start DIA Target analysis: " + e.getMessage());
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
			consolePanel.appendConsoleText("Stopping DIA Target analysis...\n");
			LOGGER.info("stopAnalysis called - user confirmed stop");

			// Force stop the running task if it exists
			if (currentTask != null && currentTask.isRunning()) {
				consolePanel.appendConsoleText("Force stopping DIA Target task...\n");
				currentTask.forceStop();
				consolePanel.appendConsoleText("DIA Target task stop requested.\n");
			} else {
				consolePanel.appendConsoleText("No running DIA Target task to stop.\n");
			}
			
			// Unbind progress bars to prevent memory leaks
			consolePanel.getOverallProgressBar().progressProperty().unbind();
			consolePanel.getCurrentTaskProgressBar().progressProperty().unbind();
			
			consolePanel.setButtonsEnabled(true, false, false);
			consolePanel.stopStatusTimer();
			consolePanel.getOverallProgressBar().setProgress(0);
			consolePanel.getCurrentTaskProgressBar().setProgress(0);
			consolePanel.setCurrentTaskLabel("Current Task: None");
		} else {
			LOGGER.info("stopAnalysis cancelled by user");
		}
	}

	/**
	 * Collect parameters from UI for DIA Target analysis
	 */
	public MetaParameterDiaTg collectParametersFromUI() {
		try {
			consolePanel.appendConsoleText("Collecting parameters from UI...\n");
			LOGGER.info("Starting parameter collection from UI for DIA Target workflow");

			// Validate inputs first
			consolePanel.appendConsoleText("Validating UI inputs...\n");
			if (!validateInputs()) {
				consolePanel.appendConsoleText("ERROR: Input validation failed\n");
				LOGGER.error("Input validation failed - cannot generate MetaParameterDiaTg");
				return null;
			}
			consolePanel.appendConsoleText("✓ Input validation passed\n");
			LOGGER.info("Input validation passed");

			// Collect experiment data from ioPanel
			consolePanel.appendConsoleText("Collecting experiment data from IO panel...\n");
			ObservableList<MetaPilotIOPanelFX.ExperimentRow> fileRows = ioPanel.getExperimentData();
			consolePanel.appendConsoleText("Found " + fileRows.size() + " experiment rows\n");
			if (fileRows.isEmpty()) {
				consolePanel.appendConsoleText("ERROR: No experiment data found in IO panel\n");
				LOGGER.error("No experiment data found in IO panel");
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
				// Get fraction from the IO panel, default to 1 if invalid
				try {
					fractions[i] = Integer.parseInt(row.getFraction());
					if (fractions[i] <= 0) {
						fractions[i] = 1; // Ensure positive fraction
					}
				} catch (NumberFormatException e) {
					fractions[i] = 1; // Default to 1 if not a valid integer
				}
				replicates[i] = 1; // Default replicate
			}

			// Create MetaData object
			MetaData metaData = new MetaData(rawFiles, expNames, fractions, replicates);

			// Collect metadata information from the reusable metadata panel
			ObservableList<MetaPilotDdaIgcMaxQuantPanel.MetadataRow> metadataRows = metadataPanel
					.getMetadataTableData();
			String[][] metaInfo = new String[metadataRows.size()][3];
			for (int i = 0; i < metadataRows.size(); i++) {
				MetaPilotDdaIgcMaxQuantPanel.MetadataRow row = metadataRows.get(i);
				metaInfo[i][0] = row.getAge() != null ? row.getAge() : "";
				metaInfo[i][1] = row.getGender() != null ? row.getGender() : "";
				metaInfo[i][2] = row.getCondition() != null ? row.getCondition() : "";
			}
			metaData.setMetaInfo(metaInfo);

			// Get result directory
			String result = ioPanel.getOutputDirectory();
			if (result == null || result.trim().isEmpty()) {
				System.err.println("Result directory not set");
				return null;
			}

			// Get MAG database information
			MagDbItem selectedMagDb = magDatabasePanel.getSelectedMagDatabase();
			if (selectedMagDb == null) {
				System.err.println("No MAG database selected");
				return null;
			}

			// Get thread count
			int threadCount = threadsComboBox != null ? threadsComboBox.getValue() : 4;

			// Create MetaBaseParameter
			MetaBaseParameter metaBasePar = new MetaBaseParameter(
					uk.ac.quadram.metapilot.MetaPilotApp.WorkflowType.DIA_TARGET.getAnalysisType(), // Workflow type
					metaData, result,
					selectedMagDb.getCurrentFile() != null ? selectedMagDb.getCurrentFile().getAbsolutePath() : "", // MicroDb
																													// path
					false, // appendHostDb - not used in DIA Target
					"", // hostDB - not used in DIA Target
					"DIA", // MS2 scan mode (default for DIA)
					1, // Core count
					threadCount, true // Meta workflow
			);

			// Create MetaParameterDiaTarget
			MetaParameterDiaTg metaParameterDiaTarget = new MetaParameterDiaTg(metaBasePar, selectedMagDb, true, false,
					new String[0], false, new String[0]);

			// Set DIA Target specific parameters
			metaParameterDiaTarget.setTargetLibraryName(
					libraryPathTextField != null ? libraryPathTextField.getText().trim() : "target_library");

			// Collect selected taxa
			if (taxaData != null) {
				java.util.List<String> selectedTaxaList = new java.util.ArrayList<>();
				for (TaxaRow row : taxaData) {
					if (row.getSelected()) {
						selectedTaxaList.add(row.getName());
					}
				}
				metaParameterDiaTarget.setSelectedTaxa(selectedTaxaList.toArray(new String[0]));
			}

			// Collect selected functions
			if (functionData != null) {
				java.util.List<String> selectedFunctionsList = new java.util.ArrayList<>();
				for (FunctionRow row : functionData) {
					if (row.getSelected()) {
						selectedFunctionsList.add(row.getName());
					}
				}
				metaParameterDiaTarget.setSelectedFunctions(selectedFunctionsList.toArray(new String[0]));
			}

			// Save taxonomic and functional levels
			if (taxaLevelComboBox != null) {
				metaParameterDiaTarget.setTaxonomicLevel(taxaLevelComboBox.getValue());
			}
			if (functionLevelComboBox != null) {
				metaParameterDiaTarget.setFunctionalLevel(functionLevelComboBox.getValue());
			}

			// Note: metadata info is not stored in MetaParameterDiaTg

			consolePanel.appendConsoleText("✓ Successfully collected parameters from UI\n");
			LOGGER.info("Successfully collected parameters from UI for DIA Target workflow");
			return metaParameterDiaTarget;

		} catch (Exception e) {
			String errorMsg = "ERROR: Failed to collect parameters from UI: " + e.getMessage();
			consolePanel.appendConsoleText(errorMsg + "\n");
			LOGGER.error("Error generating MetaParameterDiaTg: " + e.getMessage(), e);
			return null;
		}
	}

	/**
	 * Collect parameters from UI for auto-save (lenient validation)
	 */
	public MetaParameterDiaTg collectParametersFromUIForAutoSave() {
		try {
			// Create minimal MetaData even if no experiment data
			MetaData metaData = new MetaData(new String[0], new String[0], new int[0], new int[0]);

			// Get result directory (allow empty for auto-save)
			String result = ioPanel != null ? ioPanel.getOutputDirectory() : "";
			if (result == null)
				result = "";

			// Get MAG database information (allow null for auto-save)
			MagDbItem selectedMagDb = magDatabasePanel != null ? magDatabasePanel.getSelectedMagDatabase() : null;

			// Get thread count
			int threadCount = threadsComboBox != null ? threadsComboBox.getValue() : 4;

			// Create MetaBaseParameter with minimal data
			MetaBaseParameter metaBasePar = new MetaBaseParameter(
					uk.ac.quadram.metapilot.MetaPilotApp.WorkflowType.DIA_TARGET.getAnalysisType(), // Workflow type
					metaData, result,
					selectedMagDb != null && selectedMagDb.getCurrentFile() != null
							? selectedMagDb.getCurrentFile().getAbsolutePath()
							: "", // MicroDb path
					false, // appendHostDb - not used in DIA Target
					"", // hostDB - not used in DIA Target
					"DIA", // MS2 scan mode (default for DIA)
					1, // Core count
					threadCount, true // Meta workflow
			);

			// Create MetaParameterDiaTarget
			MetaParameterDiaTg metaParameterDiaTarget = new MetaParameterDiaTg(metaBasePar, selectedMagDb, false, false,
					new String[0], false, new String[0]);

			// Set DIA Target specific parameters (use defaults if UI components are null)
			metaParameterDiaTarget.setTargetLibraryName(
					libraryPathTextField != null ? libraryPathTextField.getText().trim() : "target_library");

			// Collect selected taxa (allow empty for auto-save)
			if (taxaData != null) {
				java.util.List<String> selectedTaxaList = new java.util.ArrayList<>();
				for (TaxaRow row : taxaData) {
					if (row.getSelected()) {
						selectedTaxaList.add(row.getName());
					}
				}
				metaParameterDiaTarget.setSelectedTaxa(selectedTaxaList.toArray(new String[0]));
			}

			// Collect selected functions (allow empty for auto-save)
			if (functionData != null) {
				java.util.List<String> selectedFunctionsList = new java.util.ArrayList<>();
				for (FunctionRow row : functionData) {
					if (row.getSelected()) {
						selectedFunctionsList.add(row.getName());
					}
				}
				metaParameterDiaTarget.setSelectedFunctions(selectedFunctionsList.toArray(new String[0]));
			}

			// Set metadata info (allow empty for auto-save)
			if (metadataPanel != null) {
				ObservableList<MetaPilotDdaIgcMaxQuantPanel.MetadataRow> metadataRows = metadataPanel
						.getMetadataTableData();
				String[][] metaInfo = new String[metadataRows.size()][3];
				for (int i = 0; i < metadataRows.size(); i++) {
					MetaPilotDdaIgcMaxQuantPanel.MetadataRow row = metadataRows.get(i);
					metaInfo[i][0] = row.getAge() != null ? row.getAge() : "";
					metaInfo[i][1] = row.getGender() != null ? row.getGender() : "";
					metaInfo[i][2] = row.getCondition() != null ? row.getCondition() : "";
				}
				// Note: metadata info is not stored in MetaParameterDiaTg
			}

			// Save taxonomic and functional levels (allow null for auto-save)
			if (taxaLevelComboBox != null) {
				metaParameterDiaTarget.setTaxonomicLevel(taxaLevelComboBox.getValue());
			}
			if (functionLevelComboBox != null) {
				metaParameterDiaTarget.setFunctionalLevel(functionLevelComboBox.getValue());
			}

			return metaParameterDiaTarget;

		} catch (Exception e) {
			LOGGER.error("Error generating MetaParameterDiaTg for auto-save: " + e.getMessage(), e);
			return null;
		}
	}

	/**
	 * Validate inputs before generating parameters
	 */
	private boolean validateInputs() {
		LOGGER.info("Starting input validation for DIA Target workflow");
		
		// Check if IO panel has experiment data
		if (ioPanel == null) {
			LOGGER.error("IO panel is null");
			showError("Validation Error", "IO panel is not initialized.");
			return false;
		}
		if (ioPanel.getExperimentData().isEmpty()) {
			LOGGER.error("No experiment data found in IO panel");
			showError("Validation Error", "Please add at least one experiment.");
			return false;
		}
		LOGGER.info("Experiment data validation passed: {} experiments found", ioPanel.getExperimentData().size());

		// Check if result directory is set
		String outputDir = ioPanel.getOutputDirectory();
		if (outputDir == null || outputDir.trim().isEmpty()) {
			LOGGER.error("Result directory not set");
			showError("Validation Error", "Please select an output directory.");
			return false;
		}
		LOGGER.info("Output directory validation passed: {}", outputDir);

		// Check if MAG database is selected
		if (magDatabasePanel == null || magDatabasePanel.getSelectedMagDatabase() == null) {
			LOGGER.error("No MAG database selected");
			showError("Validation Error", "Please select a MAG database.");
			return false;
		}
		LOGGER.info("MAG database validation passed: {}", magDatabasePanel.getSelectedMagDatabase().getCatalogueID());

		LOGGER.info("All validations passed for DIA Target workflow");
		return true;
	}

	/**
	 * Load last used parameters and apply them to the UI
	 */
	@Override
	protected void loadLastUsedParameters() {
		try {
			consolePanel.appendConsoleText("Loading last used parameters...\n");
			LOGGER.info("loadLastUsedParameters() called for DIA Target panel");
			MetaParameterDiaTg params = (MetaParameterDiaTg) parameterIO.loadLastUsedParameters();
			if (params != null) {
				consolePanel.appendConsoleText("✓ Successfully loaded last used parameters\n");
				LOGGER.info("Successfully loaded parameters from file");
				applyParametersToUI(params);
				LOGGER.info("Applied parameters to UI successfully");
			} else {
				consolePanel.appendConsoleText("No previous parameters found, creating defaults\n");
				LOGGER.info("No parameters loaded from file, creating defaults");
				MetaParameterDiaTg defaultParams = new MetaParameterDiaTg();
				applyParametersToUI(defaultParams);
				consolePanel.appendConsoleText("✓ Default parameters loaded into UI\n");
				LOGGER.info("Applied default parameters to UI successfully");
			}
		} catch (Exception e) {
			String errorMsg = "WARNING: Failed to load last used parameters: " + e.getMessage();
			consolePanel.appendConsoleText(errorMsg + "\n");
			LOGGER.warn("Failed to load last used DIA Target parameters: " + e.getMessage(), e);
			// Continue with default UI state
		}
	}

	/**
	 * Apply function selections to the current function data
	 * 
	 * @param selectedFunctionNames Array of function names to select
	 * @return true if at least one function was found and selected, false otherwise
	 */
	private boolean applyFunctionSelections(String[] selectedFunctionNames) {
		if (selectedFunctionNames == null || selectedFunctionNames.length == 0) {
			return false;
		}

		final FunctionRow[] firstSelectedFunction = { null };
		int foundCount = 0;

		for (String functionName : selectedFunctionNames) {
			for (FunctionRow row : functionData) {
				if (row.getName().equals(functionName)) {
					row.setSelected(true);
					foundCount++;
					if (firstSelectedFunction[0] == null) {
						firstSelectedFunction[0] = row;
					}
					break;
				}
			}
		}

		if (foundCount > 0) {
			if (functionTable != null) {
				functionTable.refresh();
				// Scroll to first selected function
				if (firstSelectedFunction[0] != null) {
					Platform.runLater(() -> {
						functionTable.getSelectionModel().select(firstSelectedFunction[0]);
						functionTable.scrollTo(firstSelectedFunction[0]);
					});
				}
			}
			System.out.println(
					"DEBUG: Applied " + foundCount + " of " + selectedFunctionNames.length + " function selections");
			return true;
		}

		return false;
	}

	/**
	 * Apply taxa selections to the current taxa data
	 * 
	 * @param selectedTaxaNames Array of taxon names to select
	 * @return true if at least one taxon was found and selected, false otherwise
	 */
	private boolean applyTaxaSelections(String[] selectedTaxaNames) {
		if (selectedTaxaNames == null || selectedTaxaNames.length == 0) {
			return false;
		}

		final TaxaRow[] firstSelectedTaxa = { null };
		int foundCount = 0;

		for (String taxonName : selectedTaxaNames) {
			for (TaxaRow row : taxaData) {
				if (row.getName().equals(taxonName)) {
					row.setSelected(true);
					foundCount++;
					if (firstSelectedTaxa[0] == null) {
						firstSelectedTaxa[0] = row;
					}
					break;
				}
			}
		}

		if (foundCount > 0) {
			if (taxaTable != null) {
				taxaTable.refresh();
				// Scroll to first selected taxa
				if (firstSelectedTaxa[0] != null) {
					Platform.runLater(() -> {
						taxaTable.getSelectionModel().select(firstSelectedTaxa[0]);
						taxaTable.scrollTo(firstSelectedTaxa[0]);
					});
				}
			}
			LOGGER.debug("[DIA Target] Applied {} of {} taxa selections", foundCount, selectedTaxaNames.length);
			return true;
		}

		return false;
	}

	/**
	 * Apply parameters to the UI components
	 */
	@Override
	protected void applyParametersToUI(MetaParameterDiaTg params) {
		try {
			consolePanel.appendConsoleText("Applying parameters to UI...\n");
			if (params == null) {
				consolePanel.appendConsoleText("WARNING: No parameters to apply\n");
				LOGGER.warn("applyParametersToUI called with null params");
				return;
			}

			LOGGER.info("applyParametersToUI called with valid params");

			// Apply metadata to IO panel if available
			if (params.getMetadata() != null && ioPanel != null) {
				MetaData metadata = params.getMetadata();
				String[] rawFiles = metadata.getRawFiles();
				String[] expNames = metadata.getExpNames();

				if (rawFiles != null && expNames != null && rawFiles.length > 0) {
					ObservableList<MetaPilotIOPanelFX.ExperimentRow> experimentRows = FXCollections
							.observableArrayList();
					for (int i = 0; i < rawFiles.length; i++) {
						if (rawFiles[i] != null && !rawFiles[i].trim().isEmpty()) {
							experimentRows.add(new MetaPilotIOPanelFX.ExperimentRow(rawFiles[i], expNames[i]));
						}
					}
					ioPanel.setExperimentData(experimentRows);
					consolePanel.appendConsoleText("Applied " + rawFiles.length + " experiment rows to IO panel\n");
					LOGGER.debug("[DIA Target] Applied metadata with {} raw files to IO panel", rawFiles.length);
				}
			}

			// Apply output directory if available
			if (params.getResult() != null && !params.getResult().trim().isEmpty()) {
				if (ioPanel != null) {
					ioPanel.setOutputDirectory(params.getResult());
					consolePanel.appendConsoleText("Applied output directory: " + params.getResult() + "\n");
					LOGGER.debug("[DIA Target] Applied output directory: {}", params.getResult());
				}
			}

			// Apply thread count
			if (threadsComboBox != null) {
				threadsComboBox.setValue(params.getThreadCount());
			}

			// Check if the saved database matches the currently selected database
			boolean databaseMatches = false;
			if (params.getMagDbItem() != null && magDatabasePanel != null) {
				MagDbItem savedMagDbItem = params.getMagDbItem();
				MagDbItem currentMagDbItem = magDatabasePanel.getSelectedMagDatabase();

				// Compare database identifiers
				if (currentMagDbItem != null
						&& savedMagDbItem.getCatalogueID().equals(currentMagDbItem.getCatalogueID())
						&& savedMagDbItem.getIdentifier().equals(currentMagDbItem.getIdentifier())) {
					databaseMatches = true;
					System.out.println(
							"DEBUG: Saved database matches current database: " + savedMagDbItem.getCatalogueID());
				} else {
					LOGGER.debug("[DIA Target] Saved database differs from current database");
					System.out
							.println("  Saved: " + (savedMagDbItem != null ? savedMagDbItem.getCatalogueID() : "null"));
					System.out.println(
							"  Current: " + (currentMagDbItem != null ? currentMagDbItem.getCatalogueID() : "null"));
				}
			}

			// Apply DIA Target specific parameters
			// Note: peptideDetectabilityComboBox value is not restored from parameters
			if (libraryPathTextField != null) {
				libraryPathTextField.setText(
						params.getTargetLibraryName() != null ? params.getTargetLibraryName() : "target_library");
			}
			if (libraryFilePathTextField != null) {
				libraryFilePathTextField.setText(""); // Clear the path field
			}
			// Library path checkbox is no longer used - library name is always editable
			// Note: taxonomic and functional level selections are not restored from
			// parameters

			// Apply taxonomic level first
			String savedTaxonomicLevel = params.getTaxonomicLevel();
			if (savedTaxonomicLevel != null && !savedTaxonomicLevel.equals("Choose the taxonomic level")) {
				if (taxaLevelComboBox != null) {
					taxaLevelComboBox.setValue(savedTaxonomicLevel);
					LOGGER.debug("[DIA Target] Restored taxonomic level: {}", savedTaxonomicLevel);
				}
			} else {
				LOGGER.debug("[DIA Target] No taxonomic level in saved parameters, will try to infer from taxa data");
			}

			// Apply functional level first
			String savedFunctionalLevel = params.getFunctionalLevel();
			if (savedFunctionalLevel != null && !savedFunctionalLevel.equals("Choose the functional level")) {
				if (functionLevelComboBox != null) {
					functionLevelComboBox.setValue(savedFunctionalLevel);
					LOGGER.debug("[DIA Target] Restored functional level: {}", savedFunctionalLevel);
				}
			} else {
				System.out.println(
						"DEBUG: No functional level in saved parameters, will try to infer from function data");
			}

			// Apply selected taxa only if database matches
			if (databaseMatches && params.getSelectedTaxa() != null && params.getSelectedTaxa().length > 0) {
				LOGGER.debug("[DIA Target] Applying saved taxa selections (database matches)");

				boolean taxaApplied = false;

				// Try to load taxa data for the restored taxonomic level first
				if (savedTaxonomicLevel != null && !savedTaxonomicLevel.equals("Choose the taxonomic level")) {
					loadTaxaData(); // Load taxa data for the restored level
					taxaApplied = applyTaxaSelections(params.getSelectedTaxa());
					if (taxaApplied) {
						LOGGER.debug("[DIA Target] Successfully applied taxa selections using saved taxonomic level: {}",
								savedTaxonomicLevel);
					}
				}

				// If that didn't work (or no level was saved), try to find the correct level by
				// testing each one
				if (!taxaApplied) {
					LOGGER.debug("[DIA Target] Trying to find correct taxonomic level for saved taxa...");
					String[] taxonomicLevels = { "Genus", "Species", "Phylum", "Class", "Order", "Family" };

					for (String level : taxonomicLevels) {
						if (taxaLevelComboBox != null) {
							taxaLevelComboBox.setValue(level);
							loadTaxaData();
							taxaApplied = applyTaxaSelections(params.getSelectedTaxa());
							if (taxaApplied) {
								LOGGER.debug("[DIA Target] Found correct taxonomic level: " + level);
								break;
							}
						}
					}
				}

				if (!taxaApplied) {
					LOGGER.debug("[DIA Target] Could not find taxonomic level that contains the saved taxa");
				}
			} else if (!databaseMatches) {
				LOGGER.debug("[DIA Target] Skipping taxa selections (database differs)");
			}

			// Apply selected functions only if database matches
			if (databaseMatches && params.getSelectedFunctions() != null && params.getSelectedFunctions().length > 0) {
				LOGGER.debug("[DIA Target] Applying saved function selections (database matches)");

				boolean functionsApplied = false;

				// Try to load function data for the restored functional level first
				if (savedFunctionalLevel != null && !savedFunctionalLevel.equals("Choose the functional level")) {
					loadFunctionData(); // Load function data for the restored level
					functionsApplied = applyFunctionSelections(params.getSelectedFunctions());
					if (functionsApplied) {
						System.out.println(
								"DEBUG: Successfully applied function selections using saved functional level: "
										+ savedFunctionalLevel);
					}
				}

				// If that didn't work (or no level was saved), try to find the correct level by
				// testing each one
				if (!functionsApplied) {
					LOGGER.debug("[DIA Target] Trying to find correct functional level for saved functions...");
					String[] functionalLevels = { "EC", "KO", "COG", "CAZy", "PFAM" };

					for (String level : functionalLevels) {
						if (functionLevelComboBox != null) {
							functionLevelComboBox.setValue(level);
							loadFunctionData();
							functionsApplied = applyFunctionSelections(params.getSelectedFunctions());
							if (functionsApplied) {
								LOGGER.debug("[DIA Target] Found correct functional level: " + level);
								break;
							}
						}
					}
				}

				if (!functionsApplied) {
					LOGGER.debug("[DIA Target] Could not find functional level that contains the saved functions");
				}
			} else if (!databaseMatches) {
				LOGGER.debug("[DIA Target] Skipping function selections (database differs)");
			}

			// Update parameter view after applying parameters
			Platform.runLater(() -> updateParameterView());

			consolePanel.appendConsoleText("✓ Successfully applied parameters to UI\n");
			LOGGER.info("Successfully applied parameters to UI for DIA Target workflow");

		} catch (Exception e) {
			String errorMsg = "ERROR: Failed to apply parameters to UI: " + e.getMessage();
			consolePanel.appendConsoleText(errorMsg + "\n");
			LOGGER.error("Failed to apply parameters to UI for DIA Target workflow", e);
		}
	}

	/**
	 * Calculate proper threshold for peptide detectability Based on
	 * MetaPilotMainPanelDiaTg testButton implementation
	 */
	private void calculateProperThreshold() {
		LOGGER.debug("[DIA Target] calculateProperThreshold() method called");

		// Update progress output
		if (libraryConsoleOutput != null) {
			libraryConsoleOutput.appendText("\n=== Starting Threshold Calculation ===\n");
			libraryConsoleOutput.appendText("Initializing threshold calculation...\n");
		}

		// Check if database generation panel has the database available
		if (dbGenerationPanel == null || dbGenerationPanel.getCurrentMagDbItem() == null) {
			LOGGER.debug("[DIA Target] Database generation panel or MagDbItem is null");
			if (libraryConsoleOutput != null) {
				libraryConsoleOutput.appendText("ERROR: No MAG database available for threshold calculation\n");
			}
			showError("Database Error", "No MAG database available for threshold calculation");
			return;
		}

		LOGGER.debug("[DIA Target] Database generation panel and MagDbItem are available");

		MagDbItem magDbItem = dbGenerationPanel.getCurrentMagDbItem();
		File currentFile = magDbItem.getCurrentFile();
		if (currentFile == null) {
			LOGGER.debug("[DIA Target] No current file available for MAG database");
			if (libraryConsoleOutput != null) {
				libraryConsoleOutput.appendText("ERROR: No current file available for MAG database\n");
			}
			showError("Database Error", "No current file available for MAG database");
			return;
		}

		LOGGER.debug("[DIA Target] Current file available: " + currentFile.getAbsolutePath());

		if (libraryConsoleOutput != null) {
			libraryConsoleOutput.appendText("Database path: " + currentFile.getAbsolutePath() + "\n");
		}

		File sqlDbFile = new File(currentFile, "pro_func.db");
		if (!sqlDbFile.exists()) {
			LOGGER.debug("[DIA Target] SQL functional database not found: " + sqlDbFile.getAbsolutePath());
			if (libraryConsoleOutput != null) {
				libraryConsoleOutput
						.appendText("ERROR: SQL functional database not found: " + sqlDbFile.getAbsolutePath() + "\n");
			}
			showError("Database Error", "Please generate the SQL functional database first.");
			return;
		}

		LOGGER.debug("[DIA Target] SQL functional database found: " + sqlDbFile.getAbsolutePath());

		if (libraryConsoleOutput != null) {
			libraryConsoleOutput.appendText("SQL database found: " + sqlDbFile.getAbsolutePath() + "\n");
		}

		// Validate selections - allow only taxa OR function selection
		LOGGER.debug("[DIA Target] Starting validation of user selections");
		if (libraryConsoleOutput != null) {
			libraryConsoleOutput.appendText("Validating user selections...\n");
		}

		String selectedTaxaLevel = taxaLevelComboBox.getValue();
		String selectedFunctionType = functionLevelComboBox.getValue();

		LOGGER.debug("[DIA Target] Selected taxonomic level: " + selectedTaxaLevel);
		LOGGER.debug("[DIA Target] Selected function type: " + selectedFunctionType);

		// Get selected taxa and functions
		HashSet<String> selectedTaxa = getSelectedTaxa();
		HashSet<String> selectedFunctions = getSelectedFunctions();

		LOGGER.debug("[DIA Target] Selected taxa count: " + (selectedTaxa != null ? selectedTaxa.size() : "null"));
		System.out.println(
				"DEBUG: Selected functions count: " + (selectedFunctions != null ? selectedFunctions.size() : "null"));

		// Check if at least one selection is made
		boolean hasTaxaSelection = (selectedTaxaLevel != null && !selectedTaxaLevel.equals("Choose the taxonomic level")
				&& selectedTaxa != null && !selectedTaxa.isEmpty());
		boolean hasFunctionSelection = (selectedFunctionType != null
				&& !selectedFunctionType.equals("Choose the functional level") && selectedFunctions != null
				&& !selectedFunctions.isEmpty());

		LOGGER.debug("[DIA Target] hasTaxaSelection: " + hasTaxaSelection);
		LOGGER.debug("[DIA Target] hasFunctionSelection: " + hasFunctionSelection);

		// Check if both taxonomic level and functional level are "Choose level" (no
		// filtering at all)
		boolean noTaxaLevelSelected = (selectedTaxaLevel == null
				|| selectedTaxaLevel.equals("Choose the taxonomic level"));
		boolean noFunctionLevelSelected = (selectedFunctionType == null
				|| selectedFunctionType.equals("Choose the functional level"));

		if (noTaxaLevelSelected && noFunctionLevelSelected) {
			LOGGER.debug("[DIA Target] No taxonomic or functional level selected - returning");
			if (libraryConsoleOutput != null) {
				libraryConsoleOutput.appendText("ERROR: No taxonomic or functional level selected\n");
			}
			showError("Selection Error",
					"Please select either a taxonomic level or functional level (or both) for threshold calculation.");
			return;
		}

		// Check for invalid taxa selection (taxonomic level chosen but no taxa
		// selected)
		if (selectedTaxaLevel != null && !selectedTaxaLevel.equals("Choose the taxonomic level")
				&& (selectedTaxa == null || selectedTaxa.isEmpty())) {
			if (libraryConsoleOutput != null) {
				libraryConsoleOutput.appendText(
						"ERROR: Taxonomic level '" + selectedTaxaLevel + "' is selected but no taxa are chosen\n");
			}
			showError("Selection Error", "You have selected a taxonomic level (" + selectedTaxaLevel
					+ ") but no taxa are selected. Please select taxa or change the taxonomic level to 'Choose the taxonomic level'.");
			return;
		}

		// Check for invalid function selection (function level chosen but no functions
		// selected)
		if (selectedFunctionType != null && !selectedFunctionType.equals("Choose the functional level")
				&& (selectedFunctions == null || selectedFunctions.isEmpty())) {
			if (libraryConsoleOutput != null) {
				libraryConsoleOutput.appendText("ERROR: Function level '" + selectedFunctionType
						+ "' is selected but no functions are chosen\n");
			}
			showError("Selection Error", "You have selected a functional level (" + selectedFunctionType
					+ ") but no functions are selected. Please select functions or change the functional level to 'Choose the functional level'.");
			return;
		}

		LOGGER.debug("[DIA Target] Validation passed, proceeding with threshold calculation");

		if (libraryConsoleOutput != null) {
			libraryConsoleOutput.appendText("Selected taxonomic level: " + selectedTaxaLevel + "\n");
			libraryConsoleOutput.appendText("Selected function type: " + selectedFunctionType + "\n");
			libraryConsoleOutput.appendText("Selected taxa: " + selectedTaxa.size() + " items\n");
			libraryConsoleOutput.appendText("Selected functions: " + selectedFunctions.size() + " items\n");
			libraryConsoleOutput.appendText("Starting threshold calculation...\n");
		}

		// Disable buttons during calculation
		testThresholdButton.setDisable(true);
		createLibraryButton.setDisable(true);

		// Show progress indicator
		showLibraryProgress("Calculating proper threshold...");

		LOGGER.debug("[DIA Target] Creating threshold calculation task");

		// Run calculation in background task
		Task<Void> task = new Task<Void>() {
			@Override
			protected Void call() throws Exception {
				LOGGER.debug("[DIA Target] Task call() method started");
				try {
					// Get genomes from selected taxa
					LOGGER.debug("[DIA Target] Calling getGenomesFromTaxa");
					HashSet<String> taxaGenomeSet = getGenomesFromTaxa(selectedTaxa, selectedTaxaLevel, sqlDbFile);
					LOGGER.debug("[DIA Target] getGenomesFromTaxa returned " + taxaGenomeSet.size() + " genomes");

					// Get proteins from selected functions
					LOGGER.debug("[DIA Target] Calling getProteinsFromFunctions");
					HashSet<String> funcProSet = getProteinsFromFunctions(selectedFunctions, selectedFunctionType,
							taxaGenomeSet, sqlDbFile);
					LOGGER.debug("[DIA Target] getProteinsFromFunctions returned " + funcProSet.size() + " proteins");

					// Check if predicted files exist
					LOGGER.debug("[DIA Target] Checking if predicted files exist");
					File predictedFile = new File(currentFile, "predicted");
					LOGGER.debug("[DIA Target] Predicted file path: " + predictedFile.getAbsolutePath());
					LOGGER.debug("[DIA Target] Predicted file exists: " + predictedFile.exists());
					if (predictedFile.exists()) {
						File[] files = predictedFile.listFiles();
						LOGGER.debug("[DIA Target] Predicted files count: " + (files != null ? files.length : 0));
					}

					if (!predictedFile.exists() || predictedFile.listFiles().length == 0) {
						LOGGER.debug("[DIA Target] Predicted files not found or empty - returning");
						Platform.runLater(
								() -> showError("Data Error", "Please predict the peptide detectability first."));
						return null;
					}

					LOGGER.debug("[DIA Target] Predicted files found, proceeding with threshold calculation");

					// Process predicted files to calculate threshold
					LOGGER.debug("[DIA Target] Calling calculateThresholdFromPredictedFiles");
					calculateThresholdFromPredictedFiles(predictedFile, taxaGenomeSet, funcProSet);
					LOGGER.debug("[DIA Target] calculateThresholdFromPredictedFiles completed");

				} catch (Exception e) {
					Platform.runLater(() -> {
						if (libraryConsoleOutput != null) {
							libraryConsoleOutput
									.appendText("ERROR: Threshold calculation failed - " + e.getMessage() + "\n");
						}
						showError("Calculation Error", "Failed to calculate threshold: " + e.getMessage());
					});
					e.printStackTrace();
				}
				return null;
			}

			@Override
			protected void succeeded() {
				LOGGER.debug("[DIA Target] Threshold calculation Task succeeded");
				Platform.runLater(() -> {
					if (libraryConsoleOutput != null) {
						libraryConsoleOutput.appendText("=== Threshold Calculation Completed ===\n");
						libraryConsoleOutput.appendText("Recommended threshold has been updated!\n");
					}
					testThresholdButton.setDisable(false);
					createLibraryButton.setDisable(false);
					hideLibraryProgress("✓ Threshold calculation completed successfully!");
				});
			}

			@Override
			protected void failed() {
				LOGGER.debug("[DIA Target] Threshold calculation Task failed");
				Platform.runLater(() -> {
					if (libraryConsoleOutput != null) {
						libraryConsoleOutput.appendText("=== Threshold Calculation Failed ===\n");
					}
					testThresholdButton.setDisable(false);
					createLibraryButton.setDisable(false);
					resetLibraryProgress();
				});
			}
		};

		Thread thread = new Thread(task);
		thread.setDaemon(true);
		thread.start();
	}

	/**
	 * Create peptide library Based on MetaPilotMainPanelDiaTg createButton
	 * implementation
	 */
	private void createLibrary() {
		// Update progress output
		if (libraryConsoleOutput != null) {
			libraryConsoleOutput.appendText("\n=== Starting Library Creation ===\n");
			libraryConsoleOutput.appendText("Initializing library creation process...\n");
		}

		// Validate Input/Output panel settings first
		if (ioPanel == null) {
			if (libraryConsoleOutput != null) {
				libraryConsoleOutput.appendText("ERROR: Input/Output panel not available\n");
			}
			showError("Configuration Error", "Input/Output panel not available");
			return;
		}

		// Check if raw files are set
		if (libraryConsoleOutput != null) {
			libraryConsoleOutput.appendText("Validating raw files...\n");
		}

		ObservableList<MetaPilotIOPanelFX.ExperimentRow> experiments = ioPanel.getExperimentData();
		if (experiments == null || experiments.isEmpty()) {
			if (libraryConsoleOutput != null) {
				libraryConsoleOutput.appendText("ERROR: No raw files found\n");
			}
			showError("Input Error", "Raw files are not found. Please set raw files in the Input/Output tab.");
			return;
		}

		// Check if raw files exist
		for (MetaPilotIOPanelFX.ExperimentRow exp : experiments) {
			String rawFile = exp.getFilePath();
			if (rawFile == null || rawFile.trim().isEmpty()) {
				if (libraryConsoleOutput != null) {
					libraryConsoleOutput
							.appendText("ERROR: Empty raw file path for experiment: " + exp.getExperimentName() + "\n");
				}
				showError("Input Error", "Raw file path is empty for experiment: " + exp.getExperimentName());
				return;
			}
			File file = new File(rawFile);
			if (!file.exists()) {
				if (libraryConsoleOutput != null) {
					libraryConsoleOutput.appendText("ERROR: Raw file not found: " + rawFile + "\n");
				}
				showError("Input Error", "Raw file not found: " + rawFile);
				return;
			}
		}

		if (libraryConsoleOutput != null) {
			libraryConsoleOutput.appendText("Found " + experiments.size() + " raw files\n");
		}

		// Check if result folder is set
		if (libraryConsoleOutput != null) {
			libraryConsoleOutput.appendText("Validating output directory...\n");
		}

		String resultDir = ioPanel.getOutputDirectory();
		if (resultDir == null || resultDir.trim().isEmpty()) {
			if (libraryConsoleOutput != null) {
				libraryConsoleOutput.appendText("ERROR: Result folder not set\n");
			}
			showError("Output Error", "Result folder is not set. Please set output directory in the Input/Output tab.");
			return;
		}

		if (libraryConsoleOutput != null) {
			libraryConsoleOutput.appendText("Output directory: " + resultDir + "\n");
		}

		// Check if database generation panel has the database available
		if (dbGenerationPanel == null || dbGenerationPanel.getCurrentMagDbItem() == null) {
			showError("Database Error", "No MAG database available for library creation");
			return;
		}

		MagDbItem magDbItem = dbGenerationPanel.getCurrentMagDbItem();
		File currentFile = magDbItem.getCurrentFile();
		if (currentFile == null) {
			showError("Database Error", "No current file available for MAG database");
			return;
		}

		File sqlDbFile = new File(currentFile, "pro_func.db");
		if (!sqlDbFile.exists()) {
			showError("Database Error", "Please generate the SQL functional database first.");
			return;
		}

		// Validate selections - allow only taxa OR function selection
		String selectedTaxaLevel = taxaLevelComboBox.getValue();
		String selectedFunctionType = functionLevelComboBox.getValue();

		// Get selected taxa and functions
		HashSet<String> selectedTaxa = getSelectedTaxa();
		HashSet<String> selectedFunctions = getSelectedFunctions();

		// Check if at least one selection is made
		boolean hasTaxaSelection = (selectedTaxaLevel != null && !selectedTaxaLevel.equals("Choose the taxonomic level")
				&& selectedTaxa != null && !selectedTaxa.isEmpty());
		boolean hasFunctionSelection = (selectedFunctionType != null
				&& !selectedFunctionType.equals("Choose the functional level") && selectedFunctions != null
				&& !selectedFunctions.isEmpty());

		// Check if both taxonomic level and functional level are "Choose level" (no
		// filtering at all)
		boolean noTaxaLevelSelected = (selectedTaxaLevel == null
				|| selectedTaxaLevel.equals("Choose the taxonomic level"));
		boolean noFunctionLevelSelected = (selectedFunctionType == null
				|| selectedFunctionType.equals("Choose the functional level"));

		if (noTaxaLevelSelected && noFunctionLevelSelected) {
			if (libraryConsoleOutput != null) {
				libraryConsoleOutput.appendText("ERROR: No taxonomic or functional level selected\n");
			}
			showError("Selection Error",
					"Please select either a taxonomic level or functional level (or both) for library creation.");
			return;
		}

		// Check for invalid taxa selection (taxonomic level chosen but no taxa
		// selected)
		if (selectedTaxaLevel != null && !selectedTaxaLevel.equals("Choose the taxonomic level")
				&& (selectedTaxa == null || selectedTaxa.isEmpty())) {
			if (libraryConsoleOutput != null) {
				libraryConsoleOutput.appendText(
						"ERROR: Taxonomic level '" + selectedTaxaLevel + "' is selected but no taxa are chosen\n");
			}
			showError("Selection Error", "You have selected a taxonomic level (" + selectedTaxaLevel
					+ ") but no taxa are selected. Please select taxa or change the taxonomic level to 'Choose the taxonomic level'.");
			return;
		}

		// Check for invalid function selection (function level chosen but no functions
		// selected)
		if (selectedFunctionType != null && !selectedFunctionType.equals("Choose the functional level")
				&& (selectedFunctions == null || selectedFunctions.isEmpty())) {
			if (libraryConsoleOutput != null) {
				libraryConsoleOutput.appendText("ERROR: Function level '" + selectedFunctionType
						+ "' is selected but no functions are chosen\n");
			}
			showError("Selection Error", "You have selected a functional level (" + selectedFunctionType
					+ ") but no functions are selected. Please select functions or change the functional level to 'Choose the functional level'.");
			return;
		}

		// Validate library name
		String libraryName = libraryPathTextField != null ? libraryPathTextField.getText().trim() : "target_library";
		if (libraryName.isEmpty()) {
			showError("Invalid Library Name", "Library name cannot be empty. Please enter a valid library name.");
			return;
		}

		// Check for invalid characters in library name
		if (libraryName.contains("/") || libraryName.contains("\\") || libraryName.contains(":")
				|| libraryName.contains("*") || libraryName.contains("?") || libraryName.contains("\"")
				|| libraryName.contains("<") || libraryName.contains(">") || libraryName.contains("|")) {
			showError("Invalid Library Name",
					"Library name contains invalid characters. Please use only letters, numbers, underscores, and hyphens.");
			return;
		}

		// Check if library already exists and ask for confirmation
		File resultFile = new File(resultDir);
		File libraryDir = new File(resultFile, libraryName);
		File fastaFile = new File(libraryDir, libraryName + ".fasta");
		File infoFile = new File(libraryDir, libraryName + ".info");
		File speclibFile = new File(libraryDir, libraryName + ".predicted.speclib");

		// Check if any library files exist
		boolean fastaExists = fastaFile.exists();
		boolean infoExists = infoFile.exists();
		boolean speclibExists = speclibFile.exists();
		boolean libraryExists = fastaExists || infoExists || speclibExists;

		if (libraryExists) {
			// Build detailed message about which files exist
			StringBuilder message = new StringBuilder();
			message.append("A library named '").append(libraryName).append("' already exists in:\n");
			message.append(libraryDir.getAbsolutePath()).append("\n\n");
			message.append("Existing files:\n");
			if (fastaExists) {
				message.append("✓ ").append(fastaFile.getName()).append("\n");
			}
			if (infoExists) {
				message.append("✓ ").append(infoFile.getName()).append("\n");
			}
			if (speclibExists) {
				message.append("✓ ").append(speclibFile.getName()).append("\n");
			}
			message.append("\nWhat would you like to do?");

			javafx.scene.control.Alert confirm = new javafx.scene.control.Alert(
					javafx.scene.control.Alert.AlertType.CONFIRMATION);
			confirm.setTitle("Library Already Exists");
			confirm.setHeaderText("Library Found");
			confirm.setContentText(message.toString());

			// Create custom buttons
			javafx.scene.control.ButtonType useExistingButton = new javafx.scene.control.ButtonType(
					"Use Existing Library");
			javafx.scene.control.ButtonType overwriteButton = new javafx.scene.control.ButtonType("Overwrite Library");
			javafx.scene.control.ButtonType cancelButton = new javafx.scene.control.ButtonType("Cancel",
					javafx.scene.control.ButtonBar.ButtonData.CANCEL_CLOSE);

			confirm.getButtonTypes().setAll(useExistingButton, overwriteButton, cancelButton);

			java.util.Optional<javafx.scene.control.ButtonType> result = confirm.showAndWait();

			if (result.isEmpty() || result.get() == cancelButton) {
				if (libraryConsoleOutput != null) {
					libraryConsoleOutput.appendText("Library creation cancelled by user.\n");
				}
				return;
			} else if (result.get() == useExistingButton) {
				// Use existing library - check if it's complete
				if (fastaExists && infoExists && speclibExists) {
					// Complete library exists
					if (libraryFilePathTextField != null) {
						libraryFilePathTextField.setText(speclibFile.getAbsolutePath());
					}
					if (libraryConsoleOutput != null) {
						libraryConsoleOutput.appendText("Using existing complete library: " + speclibFile.getAbsolutePath() + "\n");
					}
				} else {
					// Incomplete library - show warning
					StringBuilder warningMessage = new StringBuilder();
					warningMessage.append("The existing library is incomplete. Missing files:\n");
					if (!fastaExists) warningMessage.append("✗ ").append(fastaFile.getName()).append("\n");
					if (!infoExists) warningMessage.append("✗ ").append(infoFile.getName()).append("\n");
					if (!speclibExists) warningMessage.append("✗ ").append(speclibFile.getName()).append("\n");
					warningMessage.append("\nPlease choose 'Overwrite Library' to create a complete library.");

					javafx.scene.control.Alert warning = new javafx.scene.control.Alert(
							javafx.scene.control.Alert.AlertType.WARNING);
					warning.setTitle("Incomplete Library");
					warning.setHeaderText("Library Files Missing");
					warning.setContentText(warningMessage.toString());
					warning.showAndWait();

					if (libraryConsoleOutput != null) {
						libraryConsoleOutput.appendText("Cannot use incomplete library. Please choose 'Overwrite Library'.\n");
					}
					return;
				}
				return;
			}
			// If overwriteButton is selected, continue with library creation
		}

		// Get threshold value
		float threshold = (float) thresholdComboBox.getValue();

		if (libraryConsoleOutput != null) {
			libraryConsoleOutput.appendText("Library name: " + libraryName + "\n");
			libraryConsoleOutput.appendText("Detectability threshold: " + threshold + "\n");
			libraryConsoleOutput.appendText("Starting library creation process...\n");
		}

		// Disable buttons during creation
		testThresholdButton.setDisable(true);
		createLibraryButton.setDisable(true);

		// Show progress indicator
		showLibraryProgress("Creating library...");

		// Run library creation in background task
		Task<Void> task = new Task<Void>() {
			@Override
			protected Void call() throws Exception {
				try {
					Platform.runLater(() -> {
						if (libraryConsoleOutput != null) {
							libraryConsoleOutput.appendText("Processing selected taxa and functions...\n");
						}
					});

					// Get genomes from selected taxa (if any)
					HashSet<String> taxaGenomeSet = getGenomesFromTaxa(selectedTaxa, selectedTaxaLevel, sqlDbFile);

					Platform.runLater(() -> {
						if (libraryConsoleOutput != null) {
							if (taxaGenomeSet.isEmpty()) {
								libraryConsoleOutput.appendText("No taxonomic filtering applied - using all genomes\n");
							} else {
								libraryConsoleOutput
										.appendText("Found " + taxaGenomeSet.size() + " genomes from selected taxa\n");
							}
						}
					});

					// Get proteins from selected functions (if any)
					HashSet<String> funcProSet = getProteinsFromFunctions(selectedFunctions, selectedFunctionType,
							taxaGenomeSet, sqlDbFile);

					Platform.runLater(() -> {
						if (libraryConsoleOutput != null) {
							if (funcProSet.isEmpty()) {
								libraryConsoleOutput
										.appendText("No functional filtering applied - using all proteins\n");
							} else {
								libraryConsoleOutput.appendText(
										"Found " + funcProSet.size() + " proteins from selected functions\n");
							}
						}
					});

					// Check if predicted files exist
					File predictedFile = new File(currentFile, "predicted");
					if (!predictedFile.exists() || predictedFile.listFiles().length == 0) {
						Platform.runLater(() -> {
							if (libraryConsoleOutput != null) {
								libraryConsoleOutput.appendText("ERROR: Predicted peptide files not found\n");
							}
							showError("Data Error", "Please predict the peptide detectability first.");
						});
						return null;
					}

					Platform.runLater(() -> {
						if (libraryConsoleOutput != null) {
							libraryConsoleOutput.appendText("Processing predicted peptide files...\n");
						}
					});

					// Create library files
					String libraryName = libraryPathTextField != null ? libraryPathTextField.getText().trim()
							: "target_library";
					createLibraryFiles(predictedFile, taxaGenomeSet, funcProSet, selectedTaxa, selectedFunctions,
							threshold, resultDir, libraryName);

				} catch (Exception e) {
					Platform.runLater(() -> {
						if (libraryConsoleOutput != null) {
							libraryConsoleOutput
									.appendText("ERROR: Library creation failed - " + e.getMessage() + "\n");
						}
						showError("Library Creation Error", "Failed to create library: " + e.getMessage());
					});
					e.printStackTrace();
				}
				return null;
			}

			@Override
			protected void succeeded() {
				LOGGER.debug("[DIA Target] Library creation Task succeeded");
				Platform.runLater(() -> {
					if (libraryConsoleOutput != null) {
						libraryConsoleOutput.appendText("=== Library Creation Completed Successfully ===\n");
						libraryConsoleOutput.appendText("Target library is ready for analysis!\n");
					}
					testThresholdButton.setDisable(false);
					createLibraryButton.setDisable(false);
					hideLibraryProgress("✓ Library created successfully!");
				});
			}

			@Override
			protected void failed() {
				LOGGER.debug("[DIA Target] Library creation Task failed");
				Platform.runLater(() -> {
					if (libraryConsoleOutput != null) {
						libraryConsoleOutput.appendText("=== Library Creation Failed ===\n");
					}
					testThresholdButton.setDisable(false);
					createLibraryButton.setDisable(false);
					resetLibraryProgress();
				});
			}
		};

		Thread thread = new Thread(task);
		thread.setDaemon(true);
		thread.start();
	}

	/**
	 * Get selected taxa from the taxa table
	 */
	private HashSet<String> getSelectedTaxa() {
		HashSet<String> selectedTaxa = new HashSet<>();
		for (TaxaRow row : taxaData) {
			if (row.isSelected()) {
				selectedTaxa.add(row.getName());
			}
		}
		return selectedTaxa;
	}

	/**
	 * Get selected functions from the function table
	 */
	private HashSet<String> getSelectedFunctions() {
		HashSet<String> selectedFunctions = new HashSet<>();
		for (FunctionRow row : functionData) {
			if (row.isSelected()) {
				selectedFunctions.add(row.getName());
			}
		}
		return selectedFunctions;
	}

	/**
	 * Get genomes from selected taxa using database query
	 */
	private HashSet<String> getGenomesFromTaxa(HashSet<String> selectedTaxa, String selectedTaxaLevel, File sqlDbFile)
			throws SQLException {
		HashSet<String> taxaGenomeSet = new HashSet<>();

		// If taxonomic level is "Choose the taxonomic level", return all genomes (no
		// taxonomic filtering)
		if (selectedTaxaLevel == null || selectedTaxaLevel.equals("Choose the taxonomic level")) {
			Platform.runLater(() -> {
				if (libraryConsoleOutput != null) {
					libraryConsoleOutput
							.appendText("No taxonomic level selected - using all genomes (no taxonomic filtering)\n");
				}
			});
			// Get all genome names from the database (no filtering)
			String url = "jdbc:sqlite:" + sqlDbFile.getAbsolutePath().replaceAll("\\\\", "/");
			try (Connection conn = DriverManager.getConnection(url)) {
				String query = "SELECT DISTINCT genome_name FROM taxonomy";
				try (PreparedStatement pstmt = conn.prepareStatement(query); ResultSet rs = pstmt.executeQuery()) {
					while (rs.next()) {
						String genomeName = rs.getString("genome_name");
						if (genomeName != null && !genomeName.trim().isEmpty()) {
							taxaGenomeSet.add(genomeName);
						}
					}
				}
			} catch (SQLException e) {
				System.err.println("ERROR: Failed to get all genomes from database: " + e.getMessage());
				e.printStackTrace();
			}
			return taxaGenomeSet;
		}

		// If taxonomic level is selected but no taxa are chosen, this is an error
		if (selectedTaxa == null || selectedTaxa.isEmpty()) {
			throw new IllegalArgumentException("Taxonomic level '" + selectedTaxaLevel
					+ "' is selected but no taxa are chosen. Please select taxa or change the taxonomic level to 'Choose the taxonomic level'.");
		}

		String[] levelNames = { "Superkingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Genome" };
		int levelIndex = -1;
		for (int i = 0; i < levelNames.length; i++) {
			if (levelNames[i].equals(selectedTaxaLevel)) {
				levelIndex = i;
				break;
			}
		}

		String url = "jdbc:sqlite:" + sqlDbFile.getAbsolutePath().replaceAll("\\\\", "/");

		if (levelIndex < TAXONOMIC_LEVELS.length - 1) {
			// Map UI level to database column
			String columnName = TAXONOMIC_LEVELS[levelIndex + 1];
			StringBuilder sql = new StringBuilder("SELECT genome_name FROM taxonomy WHERE " + columnName + " IN (");
			boolean first = true;
			for (String taxaName : selectedTaxa) {
				if (!first) {
					sql.append(", ");
				} else {
					first = false;
				}
				sql.append("'").append(taxaName).append("'");
			}
			sql.append(")");

			try (Connection conn = DriverManager.getConnection(url);
					PreparedStatement pstmt = conn.prepareStatement(sql.toString());
					ResultSet rs = pstmt.executeQuery()) {
				while (rs.next()) {
					taxaGenomeSet.add(rs.getString("genome_name"));
				}
			}

			final HashSet<String> finalTaxaGenomeSet1 = taxaGenomeSet;
			Platform.runLater(() -> {
				if (libraryConsoleOutput != null) {
					libraryConsoleOutput.appendText(format.format(new Date()) + ":\t" + selectedTaxa.size()
							+ " taxa are selected in " + columnName + " level + " + finalTaxaGenomeSet1.size()
							+ " genomes will be selected for the following analysis\n");
				}
			});
		} else {
			// Genome level
			taxaGenomeSet = selectedTaxa;
			final HashSet<String> finalTaxaGenomeSet2 = taxaGenomeSet;
			Platform.runLater(() -> {
				if (libraryConsoleOutput != null) {
					libraryConsoleOutput.appendText(format.format(new Date()) + ":\t" + finalTaxaGenomeSet2.size()
							+ " genomes will be selected for the following analysis\n");
				}
			});
		}

		return taxaGenomeSet;
	}

	/**
	 * Get proteins from selected functions using database query
	 */
	private HashSet<String> getProteinsFromFunctions(HashSet<String> selectedFunctions, String selectedFunctionType,
			HashSet<String> taxaGenomeSet, File sqlDbFile) throws SQLException {
		HashSet<String> funcProSet = new HashSet<>();

		// If no functions selected, return all proteins (empty set means no filtering)
		if (selectedFunctions == null || selectedFunctions.isEmpty()) {
			Platform.runLater(() -> {
				if (libraryConsoleOutput != null) {
					libraryConsoleOutput.appendText(
							"No functions selected - will use all proteins from selected taxa (no functional filtering)\n");
				}
			});
			return funcProSet; // Empty set means no filtering
		}

		String url = "jdbc:sqlite:" + sqlDbFile.getAbsolutePath().replaceAll("\\\\", "/");

		StringBuilder sql = new StringBuilder(
				"SELECT protein_name FROM " + selectedFunctionType + " WHERE func_name IN (");
		boolean first = true;
		for (String funcName : selectedFunctions) {
			if (!first) {
				sql.append(", ");
			} else {
				first = false;
			}
			sql.append("'").append(funcName).append("'");
		}
		sql.append(")");

		try (Connection conn = DriverManager.getConnection(url);
				PreparedStatement pstmt = conn.prepareStatement(sql.toString());
				ResultSet rs = pstmt.executeQuery()) {
			while (rs.next()) {
				String proName = rs.getString("protein_name");
				String genomeName = proName.substring(0, proName.indexOf("_"));
				if (taxaGenomeSet.contains(genomeName)) {
					funcProSet.add(proName);
				}
			}
		}

		Platform.runLater(() -> {
			if (libraryConsoleOutput != null) {
				libraryConsoleOutput.appendText(format.format(new Date()) + ":\t" + selectedFunctions.size()
						+ " functions are selected from the " + selectedFunctionType + " database + " + funcProSet.size()
						+ " proteins will be selected for the following analysis\n");
			}
		});

		return funcProSet;
	}

	/**
	 * Calculate threshold from predicted files
	 */
	private void calculateThresholdFromPredictedFiles(File predictedFile, HashSet<String> taxaGenomeSet,
			HashSet<String> funcProSet) throws Exception {
		LOGGER.debug("[DIA Target] calculateThresholdFromPredictedFiles started");
		LOGGER.debug("[DIA Target] taxaGenomeSet size: " + taxaGenomeSet.size());
		LOGGER.debug("[DIA Target] funcProSet size: " + funcProSet.size());

		long totalMemorySize = ((OperatingSystemMXBean) ManagementFactory.getOperatingSystemMXBean())
				.getTotalMemorySize() / 1024 / 1024 / 1024;
		ExecutorService es = Executors.newFixedThreadPool((int) totalMemorySize / 2);

		ConcurrentHashMap<String, Float> pepMap = new ConcurrentHashMap<>();
		for (String genome : taxaGenomeSet) {
			es.submit(() -> {
				try (BufferedReader reader = new BufferedReader(
						new FileReader(new File(predictedFile, genome + ".predicted.tsv")))) {
					String line = reader.readLine(); // Skip header
					while ((line = reader.readLine()) != null) {
						String[] cs = line.split("\t");
						// If funcProSet is empty, include all proteins (no functional filtering)
						// If funcProSet is not empty, only include proteins that match selected
						// functions
						if (funcProSet.isEmpty() || funcProSet.contains(cs[0])) {
							if (pepMap.containsKey(cs[1])) {
								float score = Float.parseFloat(cs[2]);
								if (score > pepMap.get(cs[1])) {
									pepMap.put(cs[1], score);
								}
							} else {
								pepMap.put(cs[1], Float.parseFloat(cs[2]));
							}
						}
					}
				} catch (Exception e) {
					System.err.println("Error processing file for genome " + genome + " - " + e.getMessage());
				}
				return null;
			});
		}
		es.shutdown();
		try {
			es.awaitTermination(2, TimeUnit.HOURS);
		} catch (InterruptedException e) {
			System.err.println("Error waiting for file processing threads: " + e.getMessage());
		}

		// Calculate detectability distribution
		int[] detectabilities = new int[101];
		for (String pep : pepMap.keySet()) {
			float score = pepMap.get(pep);
			int id = (int) (score * 100);
			for (int i = 0; i <= id; i++) {
				detectabilities[i]++;
			}
		}

		// Find optimal threshold
		final int[] lowPepThresCount = { -1 };
		final int[] highPepThresCount = { -1 };
		for (int i = 0; i < detectabilities.length; i++) {
			final int finalI = i;
			Platform.runLater(() -> {
				if (libraryConsoleOutput != null) {
					libraryConsoleOutput.appendText(format.format(new Date()) + ":\t" + "Peptide detectability: "
							+ (double) finalI / 100 + "; peptide count: " + detectabilities[finalI] + "\n");
				}
			});

			if (detectabilities[i] < 3000000 && lowPepThresCount[0] == -1) {
				lowPepThresCount[0] = i;
			}

			if (detectabilities[i] < 1000000 && highPepThresCount[0] == -1) {
				highPepThresCount[0] = i;
			}
		}

		Platform.runLater(() -> {
			if (libraryConsoleOutput != null) {
				libraryConsoleOutput.appendText(format.format(new Date()) + ":\t"
						+ "The proper size of a peptide library is from 1,000,000 to 3,000,000 peptides\n");
				libraryConsoleOutput.appendText(
						format.format(new Date()) + ":\t" + "The recommended peptide detectability threshold is:\n");
			}
		});

		// Update threshold combo box
		Platform.runLater(() -> {
			if (lowPepThresCount[0] == -1) {
				if (highPepThresCount[0] == -1) {
					if (libraryConsoleOutput != null) {
						libraryConsoleOutput.appendText(
								format.format(new Date()) + ":\t" + "Peptide detectability threshold=1.0\n");
					}
					thresholdComboBox.setValue(1.0f);
				}
			} else {
				if (highPepThresCount[0] == -1) {
					if (libraryConsoleOutput != null) {
						libraryConsoleOutput
								.appendText(format.format(new Date()) + ":\t" + "Peptide detectability threshold>"
										+ (double) lowPepThresCount[0] / 100 + " and <=1.0\n");
					}
					thresholdComboBox.setValue((float) lowPepThresCount[0] / 100);
				} else {
					if (highPepThresCount[0] > lowPepThresCount[0]) {
						if (libraryConsoleOutput != null) {
							libraryConsoleOutput.appendText(format.format(new Date()) + ":\t"
									+ "Peptide detectability threshold>=" + (double) lowPepThresCount[0] / 100
									+ " and <=" + (double) highPepThresCount[0] / 100 + "\n");
						}
					} else {
						if (libraryConsoleOutput != null) {
							libraryConsoleOutput.appendText(format.format(new Date()) + ":\t"
									+ "Peptide detectability threshold>=" + (double) lowPepThresCount[0] / 100 + "\n");
						}
					}
					thresholdComboBox.setValue((float) lowPepThresCount[0] / 100);
				}
			}

			// Debug output to show threshold calculation completed
			LOGGER.debug("[DIA Target] Threshold calculation completed successfully");
			if (libraryConsoleOutput != null) {
				libraryConsoleOutput.appendText("=== Threshold Calculation Completed Successfully ===\n");
			}
		});
	}

	/**
	 * Create library files
	 */
	private void createLibraryFiles(File predictedFile, HashSet<String> taxaGenomeSet, HashSet<String> funcProSet,
			HashSet<String> selectedTaxa, HashSet<String> selectedFunctions, float threshold, String resultDir,
			String libraryName) throws Exception {
		long totalMemorySize = ((OperatingSystemMXBean) ManagementFactory.getOperatingSystemMXBean())
				.getTotalMemorySize() / 1024 / 1024 / 1024;
		ExecutorService es = Executors.newFixedThreadPool((int) totalMemorySize / 2);

		// Debug output
		Platform.runLater(() -> {
			if (libraryConsoleOutput != null) {
				libraryConsoleOutput.appendText("DEBUG: createLibraryFiles called with:\n");
				libraryConsoleOutput.appendText("  - taxaGenomeSet size: " + taxaGenomeSet.size() + "\n");
				libraryConsoleOutput.appendText("  - funcProSet size: " + funcProSet.size() + "\n");
				libraryConsoleOutput.appendText("  - threshold: " + threshold + "\n");
				libraryConsoleOutput.appendText("  - predictedFile: " + predictedFile.getAbsolutePath() + "\n");
				libraryConsoleOutput.appendText("  - predictedFile exists: " + predictedFile.exists() + "\n");
				if (predictedFile.exists()) {
					String[] files = predictedFile.list();
					libraryConsoleOutput
							.appendText("  - predicted files count: " + (files != null ? files.length : 0) + "\n");
				}
			}
		});

		HashSet<String> pepSet = new HashSet<>();
		for (String genome : taxaGenomeSet) {
			es.submit(() -> {
				try {
					File genomeFile = new File(predictedFile, genome + ".predicted.tsv");
					if (!genomeFile.exists()) {
						Platform.runLater(() -> {
							if (libraryConsoleOutput != null) {
								libraryConsoleOutput
										.appendText("WARNING: Predicted file not found for genome: " + genome + "\n");
							}
						});
						return null;
					}

					final int[] lineCount = { 0 };
					final int[] peptideCount = { 0 };
					try (BufferedReader reader = new BufferedReader(new FileReader(genomeFile))) {
						String line = reader.readLine(); // Skip header
						while ((line = reader.readLine()) != null) {
							lineCount[0]++;
							String[] cs = line.split("\t");
							// If no functions selected (funcProSet is empty), include all proteins
							// Otherwise, only include proteins that match selected functions
							if (funcProSet.isEmpty() || funcProSet.contains(cs[0])) {
								float score = Float.parseFloat(cs[2]);
								if (score > threshold) {
									pepSet.add(cs[1]);
									peptideCount[0]++;
								}
							}
						}
					}

					Platform.runLater(() -> {
						if (libraryConsoleOutput != null) {
							libraryConsoleOutput.appendText("Processed genome " + genome + ": " + lineCount[0]
									+ " lines + " + peptideCount[0] + " peptides above threshold " + threshold + "\n");
						}
					});

				} catch (Exception e) {
					Platform.runLater(() -> {
						if (libraryConsoleOutput != null) {
							libraryConsoleOutput.appendText(
									"ERROR processing file for genome " + genome + " - " + e.getMessage() + "\n");
						}
					});
				}
				return null;
			});
		}
		es.shutdown();
		try {
			es.awaitTermination(2, TimeUnit.HOURS);
		} catch (InterruptedException e) {
			System.err.println("Error waiting for file processing threads: " + e.getMessage());
		}

		// Debug output for final results
		Platform.runLater(() -> {
			if (libraryConsoleOutput != null) {
				libraryConsoleOutput.appendText("DEBUG: Final results:\n");
				libraryConsoleOutput.appendText("  - Total unique peptides found: " + pepSet.size() + "\n");
				libraryConsoleOutput.appendText("  - Threshold used: " + threshold + "\n");
				libraryConsoleOutput
						.appendText("  - Functional filtering: " + (funcProSet.isEmpty() ? "NONE (using all proteins)"
								: "Applied to " + funcProSet.size() + " proteins") + "\n");
			}
		});

		// Create output directory with library name
		File libraryDir = new File(resultDir, libraryName);
		LOGGER.debug("[DIA Target] Library directory path: " + libraryDir.getAbsolutePath());
		LOGGER.debug("[DIA Target] Library directory exists before creation: " + libraryDir.exists());
		if (!libraryDir.exists()) {
			boolean created = libraryDir.mkdirs();
			LOGGER.debug("[DIA Target] Library directory created: " + created);
			LOGGER.debug("[DIA Target] Library directory exists after creation: " + libraryDir.exists());
		} else {
			LOGGER.debug("[DIA Target] Library directory already exists");
		}

		// Create library files with library name
		File fastaFile = new File(libraryDir, libraryName + ".fasta");
		File infoFile = new File(libraryDir, libraryName + ".info");
		File libFile = new File(libraryDir, libraryName + ".speclib");

		// Write FASTA file
		LOGGER.debug("[DIA Target] Writing FASTA file: " + fastaFile.getAbsolutePath());
		LOGGER.debug("[DIA Target] FASTA file will contain " + pepSet.size() + " peptides");
		try (PrintWriter writer = new PrintWriter(fastaFile)) {
			int count = 0;
			for (String pep : pepSet) {
				writer.println(">pro" + count + " pro" + count);
				writer.println(pep);
				count++;
			}
		}
		LOGGER.debug("[DIA Target] FASTA file written successfully");

		// Write info file
		LOGGER.debug("[DIA Target] Writing info file: " + infoFile.getAbsolutePath());
		try (PrintWriter writer = new PrintWriter(infoFile)) {
			writer.println("Genomes:");
			for (String genome : taxaGenomeSet) {
				writer.println(genome);
			}
			writer.println("Functions:");
			for (String func : selectedFunctions) {
				writer.println(func);
			}
			writer.println("Proteins:");
			for (String pro : funcProSet) {
				writer.println(pro);
			}
			writer.println("Peptide detectability threshold:");
			writer.print(threshold);
		}
		LOGGER.debug("[DIA Target] Info file written successfully");

		// Execute DiaNN task to create spectral library
		Platform.runLater(() -> {
			if (libraryConsoleOutput != null) {
				libraryConsoleOutput.appendText(
						format.format(new Date()) + ":\tExecuting DiaNN task to create spectral library...\n");
			}
		});

		LOGGER.debug("[DIA Target] About to execute DiaNN task");
		LOGGER.debug("[DIA Target] FASTA file: " + fastaFile.getAbsolutePath());
		LOGGER.debug("[DIA Target] Library file: " + libFile.getAbsolutePath());
		LOGGER.debug("[DIA Target] Threads: " + threadsComboBox.getValue());

		try {
			DiaNNTask diaNNTask = new DiaNNTask(resourceParams.getToolPath(ToolType.DIANN));
			DiannParameter diannPar = new DiannParameter();
			diannPar.setThreads(threadsComboBox.getValue()); // Use selected thread count
			diaNNTask.addTask(diannPar, fastaFile.getAbsolutePath(), libFile.getAbsolutePath(), false);

			LOGGER.debug("[DIA Target] DiaNN task created, about to run...");
			diaNNTask.run(2);
			
			if (!libFile.exists()) {
				File tempLibFile = new File(libraryDir, libraryName + ".predicted.speclib");
				if (tempLibFile.exists()) {
					tempLibFile.renameTo(libFile);
				}
			}

			if (libFile.exists()) {
				LOGGER.debug("[DIA Target] DiaNN task completed successfully");
			} else {
				LOGGER.debug("[DIA Target] DiaNN task failed");
			}

			Platform.runLater(() -> {
				if (libraryConsoleOutput != null) {
					libraryConsoleOutput.appendText(format.format(new Date()) + ":\tLibrary created successfully with "
							+ pepSet.size() + " peptides\n");
					libraryConsoleOutput.appendText("FASTA file: " + fastaFile.getAbsolutePath() + "\n");
					libraryConsoleOutput.appendText("Info file: " + infoFile.getAbsolutePath() + "\n");
					libraryConsoleOutput.appendText("Spectral library: " + libFile.getAbsolutePath() + "\n");
				}
			});

			// Update library file path text field if it exists
			Platform.runLater(() -> {
				if (libraryFilePathTextField != null) {
					libraryFilePathTextField.setText(libFile.getAbsolutePath());
				}
			});

		} catch (Exception e) {
			System.err.println("ERROR: Failed to execute DiaNN task: " + e.getMessage());
			System.err.println("DEBUG: Exception type: " + e.getClass().getName());
			e.printStackTrace();
			throw new Exception("DiaNN task execution failed: " + e.getMessage(), e);
		}
	}

	/**
	 * Show error dialog
	 */

	/**
	 * Save current parameters to file
	 */
	public void saveCurrentParameters() {
		LOGGER.debug("[DIA Target] saveCurrentParameters() called");
		try {
			LOGGER.debug("[DIA Target] Attempting to collect parameters from UI...");
			MetaParameterDiaTg currentParams = collectParametersFromUI();
			System.out.println(
					"DEBUG: collectParametersFromUI() returned: " + (currentParams != null ? "SUCCESS" : "NULL"));

			if (currentParams != null) {
				LOGGER.debug("[DIA Target] Parameters collected successfully, attempting to save...");
				LOGGER.debug("[DIA Target] ParameterIO instance: " + (parameterIO != null ? "EXISTS" : "NULL"));
				LOGGER.debug("[DIA Target] Parameter type: " + currentParams.getClass().getSimpleName());
				LOGGER.debug("[DIA Target] Base parameter: " + (currentParams != null ? "EXISTS" : "NULL"));
				LOGGER.debug("[DIA Target] MAG DB item: " + (currentParams.getMagDbItem() != null ? "EXISTS" : "NULL"));
				LOGGER.debug("[DIA Target] Target library name: " + currentParams.getTargetLibraryName());

				boolean saved = parameterIO.saveLastUsedParameters(currentParams);
				LOGGER.debug("[DIA Target] saveLastUsedParameters() returned: " + saved);

				if (saved) {
					System.out.println(
							"DEBUG: Successfully saved current parameters to last_used_dia_target_parameters.json");
				} else {
					System.err.println(
							"WARNING: Failed to save current parameters - saveLastUsedParameters returned false");
				}
			} else {
				System.err.println("WARNING: No parameters to save - collectParametersFromUI returned null");
				LOGGER.debug("[DIA Target] Attempting fallback with collectParametersFromUIForAutoSave...");
				MetaParameterDiaTg fallbackParams = collectParametersFromUIForAutoSave();
				LOGGER.debug("[DIA Target] collectParametersFromUIForAutoSave() returned: "
						+ (fallbackParams != null ? "SUCCESS" : "NULL"));

				if (fallbackParams != null) {
					LOGGER.debug("[DIA Target] Attempting to save fallback parameters...");
					boolean saved = parameterIO.saveLastUsedParameters(fallbackParams);
					LOGGER.debug("[DIA Target] Fallback save returned: " + saved);
					if (saved) {
						LOGGER.debug("[DIA Target] Successfully saved fallback parameters");
					} else {
						System.err.println("WARNING: Failed to save fallback parameters");
					}
				}
			}
		} catch (Exception e) {
			System.err.println("ERROR: Failed to save current parameters: " + e.getMessage());
			System.err.println("DEBUG: Exception type: " + e.getClass().getSimpleName());
			e.printStackTrace();
		}
	}

	/**
	 * Load analysis parameters from a specific result directory
	 */
	public void loadAnalysisParameters(String resultDirectory) {
		try {
			LOGGER.info("Loading DIA Target analysis parameters from: {}", resultDirectory);
			MetaParameterDiaTg params = (MetaParameterDiaTg) parameterIO.loadAnalysisParameters(resultDirectory);
			if (params != null) {
				applyParametersToUI(params);
				LOGGER.info("Loaded DIA Target analysis parameters from {}", resultDirectory);
				consolePanel.appendConsoleText("✓ Loaded DIA Target analysis parameters from " + resultDirectory + "\n");
			} else {
				LOGGER.warn("No parameters found in {}", resultDirectory);
				consolePanel.appendConsoleText("WARNING: No parameters found in " + resultDirectory + "\n");
			}
		} catch (Exception e) {
			LOGGER.error("Failed to load DIA Target analysis parameters from {}", resultDirectory, e);
			consolePanel.appendConsoleText("WARNING: Failed to load parameters from " + resultDirectory + ": " + e.getMessage() + "\n");
		}
	}
	
	/**
	 * Get the parameter class for this workflow
	 */
	@Override
	protected Class<MetaParameterDiaTg> getParameterClass() {
		return MetaParameterDiaTg.class;
	}

	/**
	 * Target library information container
	 */
	public static class TargetLibraryInfo {
		private HashSet<String> genomes;
		private HashSet<String> functions;
		private HashSet<String> proteins;
		private float detectabilityThreshold;

		public TargetLibraryInfo(HashSet<String> genomes, HashSet<String> functions, HashSet<String> proteins,
				float detectabilityThreshold) {
			this.genomes = genomes != null ? genomes : new HashSet<>();
			this.functions = functions != null ? functions : new HashSet<>();
			this.proteins = proteins != null ? proteins : new HashSet<>();
			this.detectabilityThreshold = detectabilityThreshold;
		}

		public HashSet<String> getGenomes() {
			return genomes;
		}

		public HashSet<String> getFunctions() {
			return functions;
		}

		public HashSet<String> getProteins() {
			return proteins;
		}

		public float getDetectabilityThreshold() {
			return detectabilityThreshold;
		}

		public boolean hasGenomeFiltering() {
			return !genomes.isEmpty();
		}

		public boolean hasFunctionFiltering() {
			return !functions.isEmpty();
		}

		public boolean hasProteinFiltering() {
			return !proteins.isEmpty();
		}
	}

	/**
	 * Read target library information from the info file
	 */
	private TargetLibraryInfo readTargetLibraryInfo(File libraryDir, String libraryName) {
		File infoFile = new File(libraryDir, libraryName + ".info");

		if (!infoFile.exists()) {
			LOGGER.error("Target library info file not found: {}", infoFile.getAbsolutePath());
			return null;
		}

		HashSet<String> genomes = new HashSet<>();
		HashSet<String> functions = new HashSet<>();
		HashSet<String> proteins = new HashSet<>();
		float detectabilityThreshold = 0.5f; // Default threshold

		try (BufferedReader reader = new BufferedReader(new FileReader(infoFile))) {
			String line;
			String currentSection = null;

			while ((line = reader.readLine()) != null) {
				line = line.trim();

				if (line.equals("Genomes:")) {
					currentSection = "genomes";
				} else if (line.equals("Functions:")) {
					currentSection = "functions";
				} else if (line.equals("Proteins:")) {
					currentSection = "proteins";
				} else if (line.equals("Peptide detectability threshold:")) {
					currentSection = "threshold";
				} else if (!line.isEmpty()) {
					switch (currentSection) {
					case "genomes":
						genomes.add(line);
						break;
					case "functions":
						functions.add(line);
						break;
					case "proteins":
						proteins.add(line);
						break;
					case "threshold":
						try {
							detectabilityThreshold = Float.parseFloat(line);
						} catch (NumberFormatException e) {
							LOGGER.warn("Invalid threshold value in info file: {}", line);
						}
						break;
					}
				}
			}
		} catch (IOException e) {
			LOGGER.error("Failed to read target library info file: {}", e.getMessage());
			return null;
		}

		LOGGER.info("Read target library info: {} genomes, {} functions, {} proteins, threshold={}", genomes.size(),
				functions.size(), proteins.size(), detectabilityThreshold);

		return new TargetLibraryInfo(genomes, functions, proteins, detectabilityThreshold);
	}

}