package uk.ac.quadram.metapilot.task.dda.mag.gui;

import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.control.*;
import javafx.scene.layout.*;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.scene.text.Font;
import javafx.scene.text.FontWeight;
import javafx.scene.paint.Color;
import javafx.scene.control.Alert.AlertType;
import java.util.Optional;
import uk.ac.quadram.metapilot.task.common.gui.MetaPilotIOPanelFX;
import uk.ac.quadram.metapilot.task.common.gui.ConsolePanel;
import uk.ac.quadram.metapilot.task.common.gui.MetadataConfigurationPanel;
import uk.ac.quadram.metapilot.task.common.gui.IsobaricLabelingPanel;
import uk.ac.quadram.metapilot.task.common.gui.ModificationsPanel;
import uk.ac.quadram.metapilot.task.common.gui.MagDatabaseSettingsPanel;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;
import uk.ac.quadram.metapilot.task.par.MetaData;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.resources.MetaResourceIO;
import uk.ac.quadram.metapilot.task.common.gui.MetaPilotIOPanelFX.ExperimentRow;
import uk.ac.quadram.metapilot.MetaPilotApp;
import javafx.application.Platform;
import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import uk.ac.quadram.metapilot.task.dda.mag.par.DDAMagSageParameters;
import uk.ac.quadram.metapilot.task.dda.mag.task.DDAMagSageTask;
import uk.ac.quadram.metapilot.task.dda.mag.par.DDAMagSageParameterIO;
import uk.ac.quadram.metapilot.task.mag.MagDbItem.MagVersion;
import uk.ac.quadram.metapilot.task.common.gui.AbstractWorkflowPanel;

/**
 * Comprehensive JavaFX panel for DDA MAG Sage workflow configuration Mirrors
 * the functionality of MetaPilotDdaMagPfindPanel with Sage-specific parameters
 */
public class MetaPilotDdaMagSagePanel extends AbstractWorkflowPanel<DDAMagSageParameters> {
	// UI Components - Database specific
	private MetaPilotIOPanelFX ioPanel;
	private MagDatabaseSettingsPanel magDatabasePanel;
	private ToggleGroup quantModeGroup;

	// Reusable UI Components
	private MetadataConfigurationPanel metadataPanel;
	private IsobaricLabelingPanel isobaricPanel;

	// Sage-specific modification tables (same as FragPipe)
	// Reusable modifications panel
	private ModificationsPanel modificationsPanel;

	// Sage-specific parameter components
	private ToggleGroup searchTypeGroup;

	// Closed search parameters (symmetric)
	private Spinner<Double> precursorToleranceSpinner;
	private Spinner<Double> fragmentToleranceSpinner;
	private Label precursorToleranceUnitLabel;
	private Label fragmentToleranceUnitLabel;

	// Open search parameters (asymmetric)
	private Spinner<Double> precursorToleranceMinSpinner;
	private Spinner<Double> precursorToleranceMaxSpinner;
	private Label openPrecursorToleranceUnitLabel;
	private Spinner<Double> openFragmentToleranceSpinner;

	// Parameter panel references for direct access
	private VBox closedSearchParamsPanel;
	private VBox openSearchParamsPanel;
	private ComboBox<String> enzymeComboBox;
	private Spinner<Integer> missedCleavagesSpinner;
	private ComboBox<String> digestionModeComboBox;
	private ComboBox<Integer> threadsComboBox;

	// Data
	private ObservableList<String> enzymes = FXCollections.observableArrayList();

	// Logger
	private static final Logger LOGGER = LogManager.getLogger(MetaPilotDdaMagSagePanel.class);

	public MetaPilotDdaMagSagePanel(MagDbItem[] magDbItems) {
		// Call parent constructor first
		super();

		LOGGER.debug("[DDA MAG Sage] Constructor started");
		try {
			// Initialize parameter IO instance
			parameterIO = new DDAMagSageParameterIO();

			loadEnzymesFromIni();
			LOGGER.debug("[DDA MAG Sage] Loading MAG databases...");
			LOGGER.debug("[DDA MAG Sage] Found {} MAG database items", magDbItems != null ? magDbItems.length : 0);
			initializeSpecificUI();
			setAvailableMagDatabases(magDbItems);
			LOGGER.info("[DDA MAG Sage] Constructor completed successfully");
		} catch (Exception e) {
			System.err.println("ERROR: Exception in MetaPilotDdaMagSagePanel constructor: " + e.getMessage());
			e.printStackTrace();
			throw e;
		}
	}

	private void initializeSpecificUI() {
		System.out.println("DEBUG: initializeUI() started");
		try {
			setSpacing(10);
			setPadding(new Insets(10));
			setMinWidth(1200);
			setMinHeight(800);
			setPrefWidth(1600);
			setPrefHeight(900);

			// Initialize reusable components
			consolePanel = new ConsolePanel();
			consolePanel.setStartCallback(this::startAnalysis);
			consolePanel.setStopCallback(this::stopAnalysis);
			consolePanel.setViewReportCallback(this::viewReport);
			metadataPanel = new MetadataConfigurationPanel();
			isobaricPanel = new IsobaricLabelingPanel();
			magDatabasePanel = new MagDatabaseSettingsPanel();

			// Set up callbacks
			metadataPanel.setParameterUpdateCallback(this::updateParameterView);
			isobaricPanel.setParameterUpdateCallback(() -> {
				updateParameterView();
				updateMetadataTable();
			});
			magDatabasePanel.setParameterUpdateCallback(this::updateParameterView);
			consolePanel.setStartCallback(this::startAnalysis);
			consolePanel.setStopCallback(this::stopAnalysis);
			consolePanel.setViewReportCallback(this::viewReport);

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

			// Add tab change listener to update metadata when Parameters tab is selected
			// and save parameters when Console tab is selected
			tabPane.getSelectionModel().selectedItemProperty().addListener((obs, oldTab, newTab) -> {
				if (newTab == paramTab) {
					updateMetadataTableFromIOPanel();
				} else if (newTab == consoleTab) {
					// Save current parameters when Console tab is selected
					try {
						consolePanel.appendConsoleText("Auto-saving parameters when switching to Console tab...\n");
						DDAMagSageParameters currentParams = collectParametersFromUIForSaving();
						if (currentParams != null) {
							boolean saved = parameterIO.saveLastUsedParameters(currentParams);
							if (saved) {
								consolePanel.appendConsoleText("✓ Parameters auto-saved to config folder\n");
								LOGGER.info(
										"Successfully auto-saved parameters to config folder when Console tab was selected");
							} else {
								consolePanel.appendConsoleText(
										"WARNING: Failed to auto-save parameters to config folder\n");
								LOGGER.warn(
										"Failed to auto-save parameters to config folder when Console tab was selected");
							}

							// Also save preferences
							boolean prefsSaved = parameterIO.savePreferences(currentParams);
							if (prefsSaved) {
								consolePanel.appendConsoleText("✓ Preferences auto-saved\n");
							} else {
								consolePanel.appendConsoleText("WARNING: Failed to auto-save preferences\n");
							}
						} else {
							consolePanel.appendConsoleText("WARNING: No parameters to auto-save\n");
							LOGGER.warn("No parameters to auto-save when Console tab was selected");
						}
					} catch (Exception e) {
						String errorMsg = "ERROR: Failed to auto-save parameters: " + e.getMessage();
						consolePanel.appendConsoleText(errorMsg + "\n");
						LOGGER.warn("Failed to auto-save parameters when Console tab was selected: " + e.getMessage(),
								e);
					}
				}
			});

			getChildren().add(tabPane);

			// Load parameters after UI is fully initialized
			Platform.runLater(() -> {
				loadLastUsedParameters();

				// Add window close event handler to save parameters
				setupWindowCloseHandler();
			});

			System.out.println("DEBUG: initializeUI() completed successfully");

		} catch (Exception e) {
			System.err.println("ERROR: Exception in initializeUI(): " + e.getMessage());
			e.printStackTrace();
			throw e;
		}
	}

	@Override
	protected String getWorkflowName() {
		return "DDA MAG Sage";
	}

	@Override
	public VBox createIOTab() {
		VBox tabContent = new VBox(10);
		tabContent.setPadding(new Insets(10));
		ioPanel = new MetaPilotIOPanelFX();
		tabContent.getChildren().addAll(createSectionTitle("Input/Output Settings"), ioPanel,
				createSectionTitle("Database Settings"), magDatabasePanel);
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

		// Left column: Sage-specific parameters, Enzyme Settings, Additional
		// Parameters, and Modifications Panel
		VBox leftColumn = new VBox(15);
		leftColumn.setMinWidth(430);
		VBox sageParamsSection = createSageParametersSection();
		VBox enzymeSection = createEnzymeSection();
		VBox additionalParamsSection = createAdditionalParametersSection();

		// Create horizontal modifications panel with fixed mods (narrow) and variable
		// mods (wider)
		HBox modificationsPanel = createModificationsPanel();

		leftColumn.getChildren().addAll(sageParamsSection, enzymeSection, additionalParamsSection, modificationsPanel);
		GridPane.setVgrow(leftColumn, Priority.ALWAYS);

		// Right column: Quantitative Mode and Metadata Configuration
		VBox rightColumn = new VBox(15);
		rightColumn.setMinWidth(430);
		VBox quantModeSection = createQuantitativeModeSection();

		rightColumn.getChildren().addAll(quantModeSection, metadataPanel);
		GridPane.setVgrow(rightColumn, Priority.ALWAYS);

		// Add columns to grid
		mainGrid.add(leftColumn, 0, 0);
		mainGrid.add(rightColumn, 1, 0);

		tabContent.getChildren().add(mainGrid);
		return tabContent;
	}

	private VBox createSageParametersSection() {
		VBox section = new VBox(10);
		section.setPadding(new Insets(10));
		section.setStyle("-fx-background-color: #e8f5e8; -fx-background-radius: 5;");

		Label titleLabel = createSectionTitle("Sage-Specific Parameters");

		// Search type selection
		VBox searchTypeBox = new VBox(8);
		Label searchTypeLabel = new Label("Search Type:");
		searchTypeLabel.setStyle("-fx-font-weight: bold;");

		searchTypeGroup = new ToggleGroup();
		RadioButton closedSearchRadio = new RadioButton("Closed Search");
		closedSearchRadio.setToggleGroup(searchTypeGroup);
		closedSearchRadio.setSelected(true); // Default to closed search

		RadioButton openSearchRadio = new RadioButton("Open Search");
		openSearchRadio.setToggleGroup(searchTypeGroup);

		HBox radioBox = new HBox(15);
		radioBox.setAlignment(Pos.CENTER_LEFT);
		radioBox.getChildren().addAll(closedSearchRadio, openSearchRadio);

		searchTypeBox.getChildren().addAll(searchTypeLabel, radioBox);

		// Create tolerance parameters container
		VBox toleranceContainer = new VBox(10);

		// Closed search parameters (symmetric)
		VBox closedSearchParams = createClosedSearchParameters();

		// Open search parameters (asymmetric)
		VBox openSearchParams = createOpenSearchParameters();
		openSearchParams.setVisible(false);
		openSearchParams.setManaged(false);

		toleranceContainer.getChildren().addAll(closedSearchParams, openSearchParams);

		// Set up search type change listeners
		closedSearchRadio.setOnAction(e -> {
			closedSearchParams.setVisible(true);
			closedSearchParams.setManaged(true);
			openSearchParams.setVisible(false);
			openSearchParams.setManaged(false);
			updateParameterView();
		});

		openSearchRadio.setOnAction(e -> {
			closedSearchParams.setVisible(false);
			closedSearchParams.setManaged(false);
			openSearchParams.setVisible(true);
			openSearchParams.setManaged(true);
			updateParameterView();
		});

		section.getChildren().addAll(titleLabel, searchTypeBox, toleranceContainer);
		return section;
	}

	private VBox createClosedSearchParameters() {
		VBox section = new VBox(10);
		this.closedSearchParamsPanel = section; // Store reference for direct access
		section.setPadding(new Insets(5));
		section.setStyle("-fx-background-color: #f0f8ff; -fx-background-radius: 3;");

		Label titleLabel = new Label("Closed Search Parameters");
		titleLabel.setStyle("-fx-font-weight: bold; -fx-text-fill: #2E8B57;");

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

		// Precursor tolerance (symmetric)
		Label precursorLabel = new Label("Precursor Tolerance:");
		precursorLabel.setStyle("-fx-font-weight: bold;");
		precursorToleranceSpinner = new Spinner<>(1.0, 50.0, 20.0, 1.0);
		precursorToleranceSpinner.setMinWidth(100);
		precursorToleranceSpinner.setPrefWidth(100);
		precursorToleranceSpinner.valueProperty().addListener((obs, oldVal, newVal) -> updateParameterView());

		precursorToleranceUnitLabel = new Label("ppm");
		precursorToleranceUnitLabel.setStyle("-fx-font-weight: bold; -fx-text-fill: #2E8B57;");

		HBox precursorBox = new HBox(5);
		precursorBox.setAlignment(Pos.CENTER_LEFT);
		precursorBox.getChildren().addAll(precursorToleranceSpinner, precursorToleranceUnitLabel);

		VBox precursorContainer = new VBox(5);
		precursorContainer.getChildren().addAll(precursorLabel, precursorBox);

		// Fragment tolerance (symmetric)
		Label fragmentLabel = new Label("Fragment Tolerance:");
		fragmentLabel.setStyle("-fx-font-weight: bold;");
		fragmentToleranceSpinner = new Spinner<>(1.0, 50.0, 20.0, 1.0);
		fragmentToleranceSpinner.setMinWidth(100);
		fragmentToleranceSpinner.setPrefWidth(100);
		fragmentToleranceSpinner.valueProperty().addListener((obs, oldVal, newVal) -> updateParameterView());

		fragmentToleranceUnitLabel = new Label("ppm");
		fragmentToleranceUnitLabel.setStyle("-fx-font-weight: bold; -fx-text-fill: #2E8B57;");

		HBox fragmentBox = new HBox(5);
		fragmentBox.setAlignment(Pos.CENTER_LEFT);
		fragmentBox.getChildren().addAll(fragmentToleranceSpinner, fragmentToleranceUnitLabel);

		VBox fragmentContainer = new VBox(5);
		fragmentContainer.getChildren().addAll(fragmentLabel, fragmentBox);

		paramGrid.add(precursorContainer, 0, 0);
		paramGrid.add(fragmentContainer, 1, 0);

		section.getChildren().addAll(titleLabel, paramGrid);
		return section;
	}

	private VBox createOpenSearchParameters() {
		VBox section = new VBox(10);
		this.openSearchParamsPanel = section; // Store reference for direct access
		section.setPadding(new Insets(5));
		section.setStyle("-fx-background-color: #fff8dc; -fx-background-radius: 3;");

		Label titleLabel = new Label("Open Search Parameters");
		titleLabel.setStyle("-fx-font-weight: bold; -fx-text-fill: #B8860B;");

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

		// Precursor tolerance (asymmetric)
		Label precursorLabel = new Label("Precursor Tolerance:");
		precursorLabel.setStyle("-fx-font-weight: bold;");

		HBox precursorRangeBox = new HBox(5);
		precursorRangeBox.setAlignment(Pos.CENTER_LEFT);

		precursorToleranceMinSpinner = new Spinner<>(-500.0, 0.0, -500.0, 10.0);
		precursorToleranceMinSpinner.setMinWidth(80);
		precursorToleranceMinSpinner.setPrefWidth(80);
		precursorToleranceMinSpinner.valueProperty().addListener((obs, oldVal, newVal) -> updateParameterView());

		Label toLabel = new Label("to");
		toLabel.setStyle("-fx-font-weight: bold;");

		precursorToleranceMaxSpinner = new Spinner<>(-0.0, 500.0, 150.0, 10.0);
		precursorToleranceMaxSpinner.setMinWidth(80);
		precursorToleranceMaxSpinner.setPrefWidth(80);
		precursorToleranceMaxSpinner.valueProperty().addListener((obs, oldVal, newVal) -> updateParameterView());

		openPrecursorToleranceUnitLabel = new Label("Da");
		openPrecursorToleranceUnitLabel.setStyle("-fx-font-weight: bold; -fx-text-fill: #B8860B;");

		precursorRangeBox.getChildren().addAll(precursorToleranceMinSpinner, toLabel, precursorToleranceMaxSpinner,
				openPrecursorToleranceUnitLabel);

		VBox precursorContainer = new VBox(5);
		precursorContainer.getChildren().addAll(precursorLabel, precursorRangeBox);

		// Fragment tolerance (symmetric - create separate components for open search)
		Label fragmentLabel = new Label("Fragment Tolerance:");
		fragmentLabel.setStyle("-fx-font-weight: bold;");

		// Create a separate fragment tolerance spinner for open search to avoid
		// component sharing issues
		openFragmentToleranceSpinner = new Spinner<>(1.0, 50.0, 20.0, 1.0);
		openFragmentToleranceSpinner.setMinWidth(100);
		openFragmentToleranceSpinner.setPrefWidth(100);
		openFragmentToleranceSpinner.valueProperty().addListener((obs, oldVal, newVal) -> {
			// Sync with the main fragment tolerance spinner
			if (fragmentToleranceSpinner != null) {
				fragmentToleranceSpinner.getValueFactory().setValue(newVal);
			}
			updateParameterView();
		});

		Label openFragmentToleranceUnitLabel = new Label("ppm");
		openFragmentToleranceUnitLabel.setStyle("-fx-font-weight: bold; -fx-text-fill: #B8860B;");

		HBox fragmentBox = new HBox(5);
		fragmentBox.setAlignment(Pos.CENTER_LEFT);
		fragmentBox.getChildren().addAll(openFragmentToleranceSpinner, openFragmentToleranceUnitLabel);

		VBox fragmentContainer = new VBox(5);
		fragmentContainer.getChildren().addAll(fragmentLabel, fragmentBox);

		paramGrid.add(precursorContainer, 0, 0);
		paramGrid.add(fragmentContainer, 1, 0);

		section.getChildren().addAll(titleLabel, paramGrid);
		return section;
	}

	private VBox createQuantitativeModeSection() {
		VBox section = new VBox(10);
		section.setPadding(new Insets(10));
		section.setStyle("-fx-background-color: #f5f5f5; -fx-background-radius: 5;");
		section.setMinHeight(300);
		VBox.setVgrow(section, Priority.ALWAYS);

		Label titleLabel = createSectionTitle("Quantitative Mode");
		quantModeGroup = new ToggleGroup();

		RadioButton labelFreeRadio = new RadioButton("Label-free");
		labelFreeRadio.setToggleGroup(quantModeGroup);
		labelFreeRadio.setSelected(true);
		RadioButton isobaricRadio = new RadioButton("Isobaric (Not Available)");
		isobaricRadio.setToggleGroup(quantModeGroup);
		isobaricRadio.setDisable(true); // Disable isobaric option since it's not implemented yet

		HBox radioBox = new HBox(15);
		radioBox.setAlignment(Pos.CENTER_LEFT);
		radioBox.getChildren().addAll(labelFreeRadio, isobaricRadio);

		StackPane labelingStackPane = new StackPane();
		labelingStackPane.setMinHeight(200);
		labelingStackPane.setMaxHeight(Double.MAX_VALUE);
		VBox.setVgrow(labelingStackPane, Priority.ALWAYS);

		isobaricPanel.setMinHeight(200);
		isobaricPanel.setPrefHeight(220);
		isobaricPanel.setMaxHeight(Double.MAX_VALUE);
		isobaricPanel.setVisible(false);
		isobaricPanel.setManaged(false);
		labelingStackPane.getChildren().add(isobaricPanel);

		labelFreeRadio.setOnAction(e -> {
			isobaricPanel.setVisible(false);
			isobaricPanel.setManaged(false);
			updateParameterView();
			updateMetadataTable();
		});

		isobaricRadio.setOnAction(e -> {
			isobaricPanel.setVisible(true);
			isobaricPanel.setManaged(true);
			updateParameterView();
			updateMetadataTable();
		});

		section.getChildren().addAll(titleLabel, radioBox, labelingStackPane);
		return section;
	}

	private Label createSectionTitle(String title) {
		Label label = new Label(title);
		label.setFont(Font.font("System", FontWeight.BOLD, 14));
		label.setTextFill(Color.DARKBLUE);
		return label;
	}

	private VBox createEnzymeSection() {
		VBox section = new VBox(10);
		section.setPadding(new Insets(10));
		section.setStyle("-fx-background-color: #f5f5f5; -fx-background-radius: 5;");

		Label titleLabel = createSectionTitle("Enzyme Settings");

		HBox enzymeParamsBox = new HBox(15);
		enzymeParamsBox.setAlignment(Pos.CENTER_LEFT);

		VBox enzymeBox = new VBox(5);
		enzymeBox.setMinWidth(140);
		HBox.setHgrow(enzymeBox, Priority.ALWAYS);
		Label enzymeLabel = new Label("Enzyme:");
		enzymeLabel.setStyle("-fx-font-weight: bold;");
		enzymeBox.getChildren().add(enzymeLabel);

		enzymeComboBox = new ComboBox<>(enzymes);
		enzymeComboBox.setMinWidth(140);
		enzymeComboBox.setPrefWidth(160);
		enzymeComboBox.setMaxWidth(Double.MAX_VALUE);
		if (enzymes.size() > 0) {
			enzymeComboBox.setValue(enzymes.get(0));
		}
		enzymeComboBox.setOnAction(e -> updateParameterView());
		enzymeBox.getChildren().add(enzymeComboBox);

		VBox missedCleavagesBox = new VBox(5);
		missedCleavagesBox.setMinWidth(100);
		Label missedLabel = new Label("Missed Cleavages:");
		missedLabel.setStyle("-fx-font-weight: bold;");
		missedCleavagesSpinner = new Spinner<>(0, 5, 2);
		missedCleavagesSpinner.setMinWidth(100);
		missedCleavagesSpinner.setPrefWidth(100);
		missedCleavagesSpinner.setMaxWidth(Double.MAX_VALUE);
		missedCleavagesSpinner.valueProperty().addListener((obs, oldVal, newVal) -> updateParameterView());
		missedCleavagesBox.getChildren().add(missedLabel);
		missedCleavagesBox.getChildren().add(missedCleavagesSpinner);

		VBox digestionModeBox = new VBox(5);
		digestionModeBox.setMinWidth(140);
		HBox.setHgrow(digestionModeBox, Priority.ALWAYS);
		Label digestionLabel = new Label("Digestion Mode:");
		digestionLabel.setStyle("-fx-font-weight: bold;");
		digestionModeBox.getChildren().add(digestionLabel);

		digestionModeComboBox = new ComboBox<>(FXCollections.observableArrayList("Specific", "Semi free N-term",
				"Semi free C-term", "Semi specific", "Unspecific", "None"));
		digestionModeComboBox.setMinWidth(140);
		digestionModeComboBox.setPrefWidth(160);
		digestionModeComboBox.setMaxWidth(Double.MAX_VALUE);
		digestionModeComboBox.setValue("Specific");
		digestionModeComboBox.setOnAction(e -> updateParameterView());
		digestionModeBox.getChildren().add(digestionModeComboBox);

		enzymeParamsBox.getChildren().addAll(enzymeBox, missedCleavagesBox, digestionModeBox);

		section.getChildren().addAll(titleLabel, enzymeParamsBox);
		return section;
	}

	private VBox createAdditionalParametersSection() {
		VBox section = new VBox(10);
		section.setPadding(new Insets(10));
		section.setStyle("-fx-background-color: #f0f8ff; -fx-background-radius: 5;");

		Label titleLabel = createSectionTitle("Additional Parameters");

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

		Label threadsLabel = new Label("Threads:");
		threadsLabel.setStyle("-fx-font-weight: bold;");
		threadsComboBox = new ComboBox<>();
		threadsComboBox.setMinWidth(100);
		threadsComboBox.setPrefWidth(100);
		threadsComboBox.setMaxWidth(Double.MAX_VALUE);

		// Populate with thread count options (1 to available processors)
		ObservableList<Integer> threadOptions = FXCollections.observableArrayList();
		for (int i = 1; i <= Runtime.getRuntime().availableProcessors(); i++) {
			threadOptions.add(i);
		}
		threadsComboBox.setItems(threadOptions);
		threadsComboBox.setValue(4); // Default value
		threadsComboBox.valueProperty().addListener((obs, oldVal, newVal) -> updateParameterView());

		VBox threadsBox = new VBox(5);
		threadsBox.getChildren().addAll(threadsLabel, threadsComboBox);

		paramGrid.add(threadsBox, 0, 0);

		section.getChildren().addAll(titleLabel, paramGrid);
		return section;
	}

	// Combined modifications panel with fixed mods (narrow) and variable mods
	// (wider) side by side
	private ModificationsPanel createModificationsPanel() {
		modificationsPanel = new ModificationsPanel(this::updateParameterView);
		return modificationsPanel;
	}

	// Load available MAG databases
	public void setAvailableMagDatabases(MagDbItem[] magDbItems) {
		System.out.println("DEBUG: setAvailableMagDatabases called with " + (magDbItems != null ? magDbItems.length : 0)
				+ " items");
		if (magDatabasePanel != null) {
			magDatabasePanel.setAvailableMagDatabases(magDbItems);
			System.out.println("DEBUG: Set MAG databases in reusable panel");
		} else {
			System.out.println("DEBUG: magDatabasePanel is null - UI not initialized yet");
		}
	}

	@Override
	protected void startAnalysis() {
		try {
			consolePanel.appendConsoleText("Starting MAG Sage analysis...\n");
			System.out.println("DEBUG: startAnalysis called - MAG Sage workflow");

			// Collect parameters from UI
			consolePanel.appendConsoleText("Collecting parameters from UI...\n");
			DDAMagSageParameters analysisParams = collectParametersFromUI();
			if (analysisParams == null) {
				consolePanel.appendConsoleText("ERROR: Failed to collect parameters from UI.\n");
				LOGGER.error("Failed to collect parameters from UI in startAnalysis");
				return;
			}

			consolePanel.appendConsoleText("✓ Successfully collected parameters from UI\n");
			LOGGER.info("Successfully collected parameters from UI for MAG Sage workflow");

			// Store analysis parameters for viewReport
			this.currentAnalysisParams = analysisParams;

			// Save parameters for this analysis
			String outputDir = analysisParams.getResult();
			if (outputDir != null && !outputDir.trim().isEmpty()) {
				consolePanel.appendConsoleText("Saving parameters to output directory...\n");

				// Save analysis-specific parameters
				boolean analysisParamsSaved = parameterIO.saveAnalysisParameters(outputDir, analysisParams);
				if (!analysisParamsSaved) {
					consolePanel.appendConsoleText("WARNING: Failed to save analysis parameters to result directory\n");
					LOGGER.warn("Failed to save analysis parameters to {}/analysis_parameters_sage.json", outputDir);
				} else {
					consolePanel.appendConsoleText("✓ Analysis parameters saved to result directory\n");
					LOGGER.info("Saved analysis parameters to {}/analysis_parameters_sage.json", outputDir);
				}

				// Save as last used parameters
				boolean lastParamsSaved = parameterIO.saveLastUsedParameters(analysisParams);
				if (!lastParamsSaved) {
					consolePanel.appendConsoleText("WARNING: Failed to save last used parameters\n");
					LOGGER.warn("Failed to save last used parameters");
				} else {
					consolePanel.appendConsoleText("✓ Last used parameters saved\n");
				}

				// Also save to preferences for quick defaults
				boolean prefsSaved = parameterIO.savePreferences(analysisParams);
				if (prefsSaved) {
					consolePanel.appendConsoleText("✓ Preferences saved\n");
				} else {
					consolePanel.appendConsoleText("WARNING: Failed to save preferences\n");
				}
			} else {
				consolePanel.appendConsoleText("WARNING: No output directory specified, skipping parameter save\n");
				LOGGER.warn("No output directory specified, skipping parameter save");
			}

			// Load resource parameters
			consolePanel.appendConsoleText("Loading resource parameters from config...\n");
			MetaResourceIO resourceIO = new MetaResourceIO();
			MetaResourceParameters mrp = resourceIO.loadParameters();
			consolePanel.appendConsoleText("✓ Resource parameters loaded from config file\n");

			// Set up progress tracking
			consolePanel.appendConsoleText("Setting up progress tracking...\n");
			Label statusLabel = new Label("Ready to start analysis...");
			consolePanel.appendConsoleText("✓ Progress bars and status label initialized\n");

			// Create and run DDAMagSageTask - pass statusLabel directly to task (not through parameters)
			consolePanel.appendConsoleText("Creating MAG Sage task...\n");
			currentTask = new DDAMagSageTask(analysisParams, mrp, consolePanel.getOverallProgressBar(),
					consolePanel.getCurrentTaskProgressBar(), statusLabel);
			consolePanel.appendConsoleText("✓ Task created successfully\n");

			// Set up task callbacks
			consolePanel.appendConsoleText("Setting up task callbacks...\n");
			currentTask.messageProperty().addListener((obs, oldMsg, newMsg) -> {
				Platform.runLater(() -> {
					if (newMsg != null) {
						consolePanel.appendConsoleText("Current Task: " + newMsg + "\n");
						consolePanel.updateCurrentTaskFromMessage(newMsg);
						statusLabel.setText(newMsg);
					}
				});
			});

			currentTask.setOnSucceeded(e -> {
				Platform.runLater(() -> {
					// Check if the task actually succeeded by examining the result
					try {
						Boolean result = currentTask.get();
						if (result != null && result) {
							consolePanel.appendConsoleText("✓ Analysis completed successfully!\n");
							consolePanel.enableViewReportButton();
						} else {
							consolePanel.appendConsoleText(
									"⚠ Analysis completed but may not have been fully successful.\n");
						}
					} catch (Exception ex) {
						consolePanel.appendConsoleText("✗ Analysis failed: " + ex.getMessage() + "\n");
					}
					consolePanel.setButtonsEnabled(true, false, true);
					consolePanel.stopStatusTimer(); // Stop the timer when task completes
				});
			});

			currentTask.setOnFailed(e -> {
				Platform.runLater(() -> {
					consolePanel.appendConsoleText("✗ Analysis failed. Check console output for details.\n");
					consolePanel.setButtonsEnabled(true, false, false);
					consolePanel.stopStatusTimer(); // Stop the timer when task fails
				});
			});

			consolePanel.appendConsoleText("✓ Task callbacks configured\n");

			// Run task in background thread
			consolePanel.appendConsoleText("Starting task execution...\n");
			Thread taskThread = new Thread(currentTask);
			taskThread.setDaemon(true);
			taskThread.start();

			consolePanel.appendConsoleText("✓ Analysis task started successfully.\n");
			LOGGER.info("MAG Sage analysis task started successfully");

		} catch (Exception e) {
			String errorMsg = "ERROR: Failed to start analysis: " + e.getMessage();
			consolePanel.appendConsoleText(errorMsg + "\n");
			LOGGER.error("Failed to start MAG Sage analysis", e);
			e.printStackTrace();
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
			consolePanel.appendConsoleText("Stopping MAG Sage analysis...\n");
			LOGGER.info("stopAnalysis called - user confirmed stop");

			// Force stop the running task if it exists
			if (currentTask != null && currentTask.isRunning()) {
				consolePanel.appendConsoleText("Force stopping MAG Sage task...\n");
				currentTask.forceStop();
				consolePanel.appendConsoleText("MAG Sage task stop requested.\n");
			} else {
				consolePanel.appendConsoleText("No running MAG Sage task to stop.\n");
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

	private void updateMetadataTableFromIOPanel() {
		if (ioPanel == null || metadataPanel == null) {
			System.out.println(
					"DEBUG: updateMetadataTableFromIOPanel called but some components are null - UI not initialized yet");
			return;
		}

		ObservableList<MetaPilotIOPanelFX.ExperimentRow> experiments = ioPanel.getExperimentData();
		metadataPanel.updateMetadataFromIOPanel(experiments);
	}

	private void updateMetadataTable() {
		if (metadataPanel == null || ioPanel == null) {
			System.out.println(
					"DEBUG: updateMetadataTable called but metadata components are null - UI not initialized yet");
			return;
		}

		ObservableList<MetaPilotIOPanelFX.ExperimentRow> fileRows = ioPanel.getExperimentData();
		if (fileRows.isEmpty())
			return;

		if (quantModeGroup != null && quantModeGroup.getSelectedToggle() != null) {
			RadioButton selected = (RadioButton) quantModeGroup.getSelectedToggle();
			String mode = selected.getText();
			boolean hasLabels = isobaricPanel.hasLabels();

			metadataPanel.updateMetadataForQuantMode(mode, hasLabels, fileRows, isobaricPanel.getIsobaricTableData());
		}
	}

	private void loadEnzymesFromIni() {
		try (InputStream in = getClass().getClassLoader().getResourceAsStream("enzyme_default.ini")) {
			if (in != null) {
				BufferedReader reader = new BufferedReader(new InputStreamReader(in));
				Pattern pattern = Pattern.compile("EN\\d+=(.+?) ");
				String line;
				while ((line = reader.readLine()) != null) {
					Matcher m = pattern.matcher(line);
					if (m.find()) {
						String enzymeName = m.group(1).trim();
						enzymes.add(enzymeName);
					}
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	@Override
	protected void updateParameterView() {
		// Update the console panel's parameter view
		if (consolePanel != null) {
			consolePanel.setParameterText(getParameterDisplay());
		}
	}

	private String getParameterDisplay() {
		StringBuilder sb = new StringBuilder();
		sb.append("=== MAG Sage Analysis Parameters ===\n\n");

		// Sage-specific parameters
		sb.append("Sage-Specific Parameters:\n");
		if (searchTypeGroup != null && searchTypeGroup.getSelectedToggle() != null) {
			RadioButton selected = (RadioButton) searchTypeGroup.getSelectedToggle();
			String searchType = selected.getText();
			sb.append("  Search Type: ").append(searchType).append("\n");

			if (searchType.equals("Closed Search")) {
				// Closed search parameters
				if (precursorToleranceSpinner != null && precursorToleranceUnitLabel != null) {
					sb.append("  Precursor Tolerance: ").append(precursorToleranceSpinner.getValue()).append(" ")
							.append(precursorToleranceUnitLabel.getText()).append("\n");
				}
				if (fragmentToleranceSpinner != null && fragmentToleranceUnitLabel != null) {
					sb.append("  Fragment Tolerance: ").append(fragmentToleranceSpinner.getValue()).append(" ")
							.append(fragmentToleranceUnitLabel.getText()).append("\n");
				}
			} else if (searchType.equals("Open Search")) {
				// Open search parameters
				if (precursorToleranceMinSpinner != null && precursorToleranceMaxSpinner != null
						&& openPrecursorToleranceUnitLabel != null) {
					sb.append("  Precursor Tolerance: ").append(precursorToleranceMinSpinner.getValue()).append(" to ")
							.append(precursorToleranceMaxSpinner.getValue()).append(" ")
							.append(openPrecursorToleranceUnitLabel.getText()).append("\n");
				}
				if (fragmentToleranceSpinner != null && fragmentToleranceUnitLabel != null) {
					sb.append("  Fragment Tolerance: ").append(fragmentToleranceSpinner.getValue()).append(" ")
							.append(fragmentToleranceUnitLabel.getText()).append("\n");
				}
			}
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
			if (magDatabasePanel.isAppendHostDb() && magDatabasePanel.getHostDatabase() != null
					&& !magDatabasePanel.getHostDatabase().isEmpty()) {
				sb.append("  Host Database: ").append(magDatabasePanel.getHostDatabase()).append("\n");
			}
		}
		sb.append("\n");

		// Enzyme settings
		sb.append("Enzyme Settings:\n");
		if (enzymeComboBox != null) {
			sb.append("  Enzyme: ").append(enzymeComboBox.getValue()).append("\n");
		}
		if (missedCleavagesSpinner != null) {
			sb.append("  Missed Cleavages: ").append(missedCleavagesSpinner.getValue()).append("\n");
		}
		if (digestionModeComboBox != null) {
			sb.append("  Digestion Mode: ").append(digestionModeComboBox.getValue()).append("\n");
		}
		sb.append("\n");

		// Additional parameters
		sb.append("Additional Parameters:\n");
		if (threadsComboBox != null) {
			sb.append("  Threads: ").append(threadsComboBox.getValue()).append("\n");
		}
		sb.append("\n");

		// Quantitative mode
		sb.append("Quantitative Mode:\n");
		if (quantModeGroup != null && quantModeGroup.getSelectedToggle() != null) {
			RadioButton selected = (RadioButton) quantModeGroup.getSelectedToggle();
			sb.append("  Mode: ").append(selected.getText()).append("\n");

			if (selected.getText().contains("Isobaric") && isobaricPanel != null) {
				sb.append("  Isobaric Type: ").append(isobaricPanel.getSelectedIsobaricType()).append("\n");
			}
		}
		sb.append("\n");

		// Fixed modifications
		sb.append("Fixed Modifications:\n");
		if (modificationsPanel != null) {
			for (ModificationsPanel.FixedModificationRow row : modificationsPanel.getFixedModData()) {
				if (row.getDeltaMass() != 0.0) {
					sb.append("  ").append(row.getSite()).append(": ").append(row.getDeltaMass()).append("\n");
				}
			}
		}
		sb.append("\n");

		// Variable modifications
		sb.append("Variable Modifications:\n");
		if (modificationsPanel != null) {
			for (ModificationsPanel.VariableModificationRow row : modificationsPanel.getVariableModData()) {
				if (row.getUse()) {
					sb.append("  ").append(row.getModName()).append(" (").append(row.getSites()).append("): ")
							.append(row.getDeltaMass()).append("\n");
				}
			}
		}
		sb.append("\n");

		return sb.toString();
	}

	/**
	 * Collect parameters from UI for saving (less strict validation)
	 */
	private DDAMagSageParameters collectParametersFromUIForSaving() {
		try {
			LOGGER.info("Starting parameter collection from UI for saving (MAG Sage workflow)");

			// Create basic parameters with defaults
			MetaData metaData = new MetaData();
			String result = "";
			String microDb = "";
			boolean appendHostDb = false;
			String hostDB = "";

			// Try to get data from UI if available, otherwise use defaults
			if (ioPanel != null) {
				result = ioPanel.getOutputDirectory() != null ? ioPanel.getOutputDirectory() : "";

				// Try to get experiment data
				ObservableList<ExperimentRow> fileRows = ioPanel.getExperimentData();
				if (!fileRows.isEmpty()) {
					String[] rawFiles = new String[fileRows.size()];
					String[] expNames = new String[fileRows.size()];
					int[] fractions = new int[fileRows.size()];
					int[] replicates = new int[fileRows.size()];

					for (int i = 0; i < fileRows.size(); i++) {
						ExperimentRow row = fileRows.get(i);
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
						replicates[i] = 1;
					}
					metaData = new MetaData(rawFiles, expNames, fractions, replicates);
				}
			}

			// Try to get database info if available
			if (magDatabasePanel != null) {
				MagDbItem selectedMagDb = magDatabasePanel.getSelectedMagDatabase();
				if (selectedMagDb != null) {
					microDb = selectedMagDb.getCatalogueID();
				}
				appendHostDb = magDatabasePanel.isAppendHostDb();
				hostDB = magDatabasePanel.getHostDatabase() != null ? magDatabasePanel.getHostDatabase() : "";
			}

			// Get basic parameters with defaults
			String enzyme = enzymeComboBox != null && enzymeComboBox.getValue() != null ? enzymeComboBox.getValue()
					: "Trypsin";
			int missCleavages = missedCleavagesSpinner != null && missedCleavagesSpinner.getValue() != null
					? missedCleavagesSpinner.getValue()
					: 2;
			String digestMode = digestionModeComboBox != null && digestionModeComboBox.getValue() != null
					? digestionModeComboBox.getValue()
					: "Specific";
			int threadCount = threadsComboBox != null && threadsComboBox.getValue() != null ? threadsComboBox.getValue()
					: 4;

			// Create DDAMagSageParameters
			DDAMagSageParameters analysisParams = new DDAMagSageParameters(
					uk.ac.quadram.metapilot.workflow.WorkflowType.DDA_MAG_SAGE.getAnalysisType(), metaData, result, microDb, appendHostDb,
					hostDB, "HCD-FTMS", // Default MS2 scan mode for Sage
					threadCount);

			// Set basic parameters
			analysisParams.setEnzyme(enzyme);
			analysisParams.setMissCleavages(missCleavages);
			analysisParams.setDigestMode(digestMode);

			// Set Sage-specific parameters
			String searchType = "Closed Search"; // Default
			double fragmentTolerance = 20.0;
			String precursorToleranceUnit = "ppm";
			String fragmentToleranceUnit = "ppm";

			if (searchTypeGroup != null && searchTypeGroup.getSelectedToggle() != null) {
				RadioButton selected = (RadioButton) searchTypeGroup.getSelectedToggle();
				searchType = selected.getText();

				if (searchType.equals("Closed Search")) {
					// Closed search: symmetric tolerances
					fragmentTolerance = fragmentToleranceSpinner != null && fragmentToleranceSpinner.getValue() != null
							? fragmentToleranceSpinner.getValue()
							: 20.0;
					precursorToleranceUnit = "ppm";
					fragmentToleranceUnit = "ppm";

					double precursorTolerance = precursorToleranceSpinner != null
							&& precursorToleranceSpinner.getValue() != null ? precursorToleranceSpinner.getValue()
									: 20.0;
					analysisParams.setPrecursorToleranceMin(-precursorTolerance);
					analysisParams.setPrecursorToleranceMax(precursorTolerance);

				} else if (searchType.equals("Open Search")) {
					// Open search: asymmetric tolerances
					double precursorMin = precursorToleranceMinSpinner != null
							&& precursorToleranceMinSpinner.getValue() != null ? precursorToleranceMinSpinner.getValue()
									: -500.0;
					double precursorMax = precursorToleranceMaxSpinner != null
							&& precursorToleranceMaxSpinner.getValue() != null ? precursorToleranceMaxSpinner.getValue()
									: 150.0;
					// Fragment tolerance is symmetric for both open and closed search
					fragmentTolerance = openFragmentToleranceSpinner != null
							&& openFragmentToleranceSpinner.getValue() != null ? openFragmentToleranceSpinner.getValue()
									: 20.0;
					precursorToleranceUnit = "Da";
					fragmentToleranceUnit = "ppm";

					analysisParams.setPrecursorToleranceMin(precursorMin);
					analysisParams.setPrecursorToleranceMax(precursorMax);
				}
			}

			// Don't set the old precursor_tolerance field - it's only for backward
			// compatibility
			analysisParams.setFragmentTolerance(fragmentTolerance);
			analysisParams.setPrecursorToleranceUnit(precursorToleranceUnit);
			analysisParams.setFragmentToleranceUnit(fragmentToleranceUnit);

			// Set quantitative mode
			String quanMode = "Label-Free";
			if (quantModeGroup != null && quantModeGroup.getSelectedToggle() != null) {
				RadioButton selected = (RadioButton) quantModeGroup.getSelectedToggle();
				quanMode = selected.getText();
			}
			analysisParams.setQuantitativeMode(quanMode);

			// Set modifications
			List<String> fixedMods = new ArrayList<>();
			List<String> variableMods = new ArrayList<>();

			// Convert fixed modifications
			for (ModificationsPanel.FixedModificationRow row : modificationsPanel.getFixedModData()) {
				if (row.getDeltaMass() != 0.0) {
					fixedMods.add(row.getSite() + ":" + row.getDeltaMass());
				}
			}

			// Convert variable modifications
			for (ModificationsPanel.VariableModificationRow row : modificationsPanel.getVariableModData()) {
				if (row.getUse()) {
					variableMods.add(row.getModName() + ":" + row.getSites() + ":" + row.getDeltaMass());
				}
			}

			analysisParams.setFixedModifications(fixedMods);
			analysisParams.setVariableModifications(variableMods);

			// Set experiment data if available - convert JavaFX rows to String arrays
			if (ioPanel != null) {
				ObservableList<ExperimentRow> fileRows = ioPanel.getExperimentData();
				List<String[]> experimentArrays = new ArrayList<>();
				for (ExperimentRow row : fileRows) {
					experimentArrays.add(row.toArray());
				}
				analysisParams.setExperimentData(experimentArrays);
			}

			// Set metadata table data if available - convert JavaFX rows to String arrays
			if (metadataPanel != null) {
				List<String[]> metadataArrays = new ArrayList<>();
				for (var row : metadataPanel.getMetadataTableData()) {
					metadataArrays.add(row.toArray());
				}
				analysisParams.setMetadataTableData(metadataArrays);
			}

			LOGGER.info("Successfully collected parameters from UI for saving (MAG Sage workflow)");
			return analysisParams;

		} catch (Exception e) {
			LOGGER.error("Failed to collect parameters from UI for saving: " + e.getMessage(), e);
			return null;
		}
	}

	/**
	 * Collect all parameters from the UI into a DDAMagSageParameters object
	 */
	@Override
	protected DDAMagSageParameters collectParametersFromUI() {
		try {
			consolePanel.appendConsoleText("Collecting parameters from UI...\n");
			LOGGER.info("Starting parameter collection from UI for MAG Sage workflow");

			// Validate inputs first
			consolePanel.appendConsoleText("Validating UI inputs...\n");
			if (!validateInputs()) {
				consolePanel.appendConsoleText("ERROR: Input validation failed\n");
				LOGGER.error("Input validation failed");
				return null;
			}
			consolePanel.appendConsoleText("✓ Input validation passed\n");
			LOGGER.info("Input validation passed");

			// Collect experiment data from ioPanel
			consolePanel.appendConsoleText("Collecting experiment data from IO panel...\n");
			ObservableList<ExperimentRow> fileRows = ioPanel.getExperimentData();
			consolePanel.appendConsoleText("Found " + fileRows.size() + " experiment rows\n");
			LOGGER.info("Found {} experiment rows", fileRows.size());

			if (fileRows.isEmpty()) {
				LOGGER.error("No experiment data found in IO panel");
				return null;
			}

			String[] rawFiles = new String[fileRows.size()];
			String[] expNames = new String[fileRows.size()];
			int[] fractions = new int[fileRows.size()];
			int[] replicates = new int[fileRows.size()];

			for (int i = 0; i < fileRows.size(); i++) {
				ExperimentRow row = fileRows.get(i);
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
				LOGGER.debug("Experiment {}: file={}, name={}, fraction={}", i, rawFiles[i], expNames[i], fractions[i]);
			}

			// Create MetaData object
			MetaData metaData = new MetaData(rawFiles, expNames, fractions, replicates);
			LOGGER.info("MetaData created successfully with {} raw files", rawFiles.length);

			// Collect basic parameters
			String result = ioPanel.getOutputDirectory();
			LOGGER.info("Output directory: {}", result);

			MagDbItem selectedMagDb = magDatabasePanel.getSelectedMagDatabase();
			String microDb = selectedMagDb != null ? selectedMagDb.getCatalogueID() : "";
			LOGGER.info("Selected MAG database: {}", microDb);

			boolean appendHostDb = magDatabasePanel.isAppendHostDb();
			String hostDB = appendHostDb ? magDatabasePanel.getHostDatabase() : "";

			String enzyme = enzymeComboBox.getValue() != null ? enzymeComboBox.getValue() : "Trypsin";
			int missCleavages = missedCleavagesSpinner.getValue() != null ? missedCleavagesSpinner.getValue() : 2;
			String digestMode = digestionModeComboBox.getValue() != null ? digestionModeComboBox.getValue()
					: "Specific";
			int threadCount = threadsComboBox.getValue() != null ? threadsComboBox.getValue() : 4;

			// Collect Sage-specific parameters
			String searchType = "Closed Search"; // Default
			double precursorTolerance = 20.0;
			double fragmentTolerance = 20.0;
			String precursorToleranceUnit = "ppm";
			String fragmentToleranceUnit = "ppm";

			if (searchTypeGroup != null && searchTypeGroup.getSelectedToggle() != null) {
				RadioButton selected = (RadioButton) searchTypeGroup.getSelectedToggle();
				searchType = selected.getText();

				if (searchType.equals("Closed Search")) {
					// Closed search: symmetric tolerances
					precursorTolerance = precursorToleranceSpinner.getValue() != null
							? precursorToleranceSpinner.getValue()
							: 20.0;
					fragmentTolerance = fragmentToleranceSpinner.getValue() != null
							? fragmentToleranceSpinner.getValue()
							: 20.0;
					precursorToleranceUnit = "ppm";
					fragmentToleranceUnit = "ppm";
				} else if (searchType.equals("Open Search")) {
					// Open search: asymmetric precursor tolerance, symmetric fragment tolerance
					double precursorMin = precursorToleranceMinSpinner.getValue() != null
							? precursorToleranceMinSpinner.getValue()
							: -500.0;
					double precursorMax = precursorToleranceMaxSpinner.getValue() != null
							? precursorToleranceMaxSpinner.getValue()
							: 150.0;

					// Store the absolute maximum tolerance for the parameters
					precursorTolerance = Math.max(Math.abs(precursorMin), Math.abs(precursorMax));
					fragmentTolerance = openFragmentToleranceSpinner.getValue() != null
							? openFragmentToleranceSpinner.getValue()
							: 20.0;
					precursorToleranceUnit = "Da";
					fragmentToleranceUnit = "ppm";
				}
			}

			// Create DDAMagSageParameters
			DDAMagSageParameters analysisParams = new DDAMagSageParameters(
					uk.ac.quadram.metapilot.workflow.WorkflowType.DDA_MAG_SAGE.getAnalysisType(), metaData, result, microDb, appendHostDb,
					hostDB, "HCD-FTMS", // Default MS2 scan mode for Sage
					threadCount);

			// Set Sage-specific parameters
			// Don't set the old precursor_tolerance field - it's only for backward
			// compatibility
			analysisParams.setFragmentTolerance(fragmentTolerance);
			analysisParams.setPrecursorToleranceUnit(precursorToleranceUnit);
			analysisParams.setFragmentToleranceUnit(fragmentToleranceUnit);

			// Set min/max tolerance values for open search
			if (searchType.equals("Open Search")) {
				double precursorMin = precursorToleranceMinSpinner.getValue() != null
						? precursorToleranceMinSpinner.getValue()
						: -500.0;
				double precursorMax = precursorToleranceMaxSpinner.getValue() != null
						? precursorToleranceMaxSpinner.getValue()
						: 150.0;

				// Debug logging
				LOGGER.info("Open Search - saving precursor tolerance: min={}, max={}", precursorMin, precursorMax);

				analysisParams.setPrecursorToleranceMin(precursorMin);
				analysisParams.setPrecursorToleranceMax(precursorMax);
			} else {
				// For closed search, use symmetric tolerances
				analysisParams.setPrecursorToleranceMin(-precursorTolerance);
				analysisParams.setPrecursorToleranceMax(precursorTolerance);
			}

			// Set search type (this will be used to determine if it's open search)
			// For now, we'll use the existing isOpenSearch() method based on tolerance
			// units
			// You might want to add a specific field for search type in
			// DDAMagSageParameters

			// Set additional parameters
			analysisParams.setEnzyme(enzyme);
			analysisParams.setMissCleavages(missCleavages);
			analysisParams.setDigestMode(digestMode);

			// Get MAG database version
			String magDbVersion = "";
			if (selectedMagDb != null && selectedMagDb.getVersionList() != null
					&& !selectedMagDb.getVersionList().isEmpty()) {
				magDbVersion = selectedMagDb.getVersionList().get(0).getCatalogueVersion();
			}

			analysisParams.setMagDbCatalogueId(microDb);
			analysisParams.setMagDbVersion(magDbVersion);

			// Set the used version on the MagDbItem
			MagVersion selectedVersion = magDatabasePanel.getSelectedMagVersion();
			if (selectedVersion != null && selectedMagDb != null) {
				selectedMagDb.setUsedVersion(selectedVersion);
			}

			analysisParams.setMagDbItem(selectedMagDb);

			// Get quantitative mode
			String quanMode = "Label-Free";
			if (quantModeGroup != null && quantModeGroup.getSelectedToggle() != null) {
				RadioButton selected = (RadioButton) quantModeGroup.getSelectedToggle();
				quanMode = selected.getText();
			}

			analysisParams.setQuantitativeMode(quanMode);

			// Set modifications (convert from FragPipe format to standard format)
			List<String> fixedMods = new ArrayList<>();
			List<String> variableMods = new ArrayList<>();

			// Convert fixed modifications
			for (ModificationsPanel.FixedModificationRow row : modificationsPanel.getFixedModData()) {
				if (row.getDeltaMass() != 0.0) {
					fixedMods.add(row.getSite() + ":" + row.getDeltaMass());
				}
			}

			// Convert variable modifications
			for (ModificationsPanel.VariableModificationRow row : modificationsPanel.getVariableModData()) {
				if (row.getUse()) {
					variableMods.add(row.getModName() + ":" + row.getSites() + ":" + row.getDeltaMass());
				}
			}

			analysisParams.setFixedModifications(fixedMods);
			analysisParams.setVariableModifications(variableMods);

			// Set experiment data - convert JavaFX rows to String arrays
			List<String[]> experimentArrays = new ArrayList<>();
			for (ExperimentRow row : fileRows) {
				experimentArrays.add(row.toArray());
			}
			analysisParams.setExperimentData(experimentArrays);

			consolePanel.appendConsoleText("✓ Successfully collected parameters from UI\n");
			LOGGER.info("Successfully collected parameters from UI for MAG Sage workflow");
			return analysisParams;

		} catch (Exception e) {
			String errorMsg = "ERROR: Failed to collect parameters from UI: " + e.getMessage();
			consolePanel.appendConsoleText(errorMsg + "\n");
			LOGGER.error("Failed to collect parameters from UI: " + e.getMessage(), e);
			return null;
		}
	}

	/**
	 * Validate inputs before starting analysis
	 */
	private boolean validateInputs() {
		LOGGER.info("Starting input validation for MAG Sage workflow");

		// Check experiment data
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

		// Check output directory
		String outputDir = ioPanel.getOutputDirectory();
		if (outputDir == null || outputDir.trim().isEmpty()) {
			LOGGER.error("Output directory is empty or null");
			showError("Validation Error", "Please select an output directory.");
			return false;
		}

		// Check MAG database selection
		if (magDatabasePanel == null) {
			LOGGER.error("MAG database panel is null");
			showError("Validation Error", "MAG database panel is not initialized.");
			return false;
		}

		if (!magDatabasePanel.validateSelection()) {
			String errorMsg = magDatabasePanel.getValidationErrorMessage();
			LOGGER.error("MAG database validation failed: {}", errorMsg);
			showError("Validation Error", errorMsg);
			return false;
		}

		// Check enzyme selection
		if (enzymeComboBox == null) {
			LOGGER.error("Enzyme combo box is null");
			showError("Validation Error", "Enzyme selection is not initialized.");
			return false;
		}

		String selectedEnzyme = enzymeComboBox.getValue();
		if (selectedEnzyme == null) {
			LOGGER.error("No enzyme selected");
			showError("Validation Error", "Please select an enzyme.");
			return false;
		}

		LOGGER.info("All input validation checks passed successfully");
		return true;
	}

	/**
	 * Show error dialog
	 */

	/**
	 * Load last used parameters and apply them to the UI
	 */
	@Override
	protected void loadLastUsedParameters() {
		try {
			consolePanel.appendConsoleText("Loading last used parameters...\n");
			LOGGER.info("Attempting to load last used parameters for DDA MAG Sage workflow");
			DDAMagSageParameters params = (DDAMagSageParameters) parameterIO.loadLastUsedParameters();
			if (params != null) {
				consolePanel.appendConsoleText("✓ Successfully loaded last used parameters\n");
				LOGGER.info("Successfully loaded parameters from file");
				applyParametersToUI(params);
				LOGGER.info("Applied parameters to UI for DDA MAG Sage workflow");

				// Also apply experiment data to the IO panel if available - convert String arrays back to JavaFX rows
				if (params.getExperimentData() != null && !params.getExperimentData().isEmpty()) {
					List<ExperimentRow> experimentRows = new ArrayList<>();
					for (String[] data : params.getExperimentData()) {
						experimentRows.add(ExperimentRow.fromArray(data));
					}
					ioPanel.setExperimentData(experimentRows);
					consolePanel.appendConsoleText(
							"Applied " + params.getExperimentData().size() + " experiment rows to IO panel\n");
					LOGGER.info("Applied {} experiment rows to IO panel", params.getExperimentData().size());
				}

				// Apply output directory if available
				if (params.getResult() != null && !params.getResult().trim().isEmpty()) {
					ioPanel.setOutputDirectory(params.getResult());
					consolePanel.appendConsoleText("Applied output directory: " + params.getResult() + "\n");
					LOGGER.info("Applied output directory: {}", params.getResult());
				}

				// Update metadata table after loading parameters
				Platform.runLater(() -> {
					updateMetadataTable();
					updateParameterView();
				});
			} else {
				consolePanel.appendConsoleText("No previous parameters found, creating defaults\n");
				LOGGER.info("No parameters loaded from file, creating defaults");
				DDAMagSageParameters defaultParams = new DDAMagSageParameters();
				applyParametersToUI(defaultParams);
				consolePanel.appendConsoleText("✓ Default parameters loaded into UI\n");
				LOGGER.info("Applied default parameters to UI successfully");
			}
		} catch (Exception e) {
			String errorMsg = "WARNING: Failed to load last used parameters: " + e.getMessage();
			consolePanel.appendConsoleText(errorMsg + "\n");
			LOGGER.error("Failed to load last used MAG Sage parameters: " + e.getMessage(), e);
			// Continue with default UI state
		}
	}

	/**
	 * Apply parameters to the UI components
	 */
	@Override
	protected void applyParametersToUI(DDAMagSageParameters params) {
		try {
			consolePanel.appendConsoleText("Applying parameters to UI...\n");
			if (params == null) {
				consolePanel.appendConsoleText("WARNING: No parameters to apply\n");
				LOGGER.warn("No parameters to apply to UI");
				return;
			}

			// Apply metadata to IO panel if available
			if (params.getMetadata() != null) {
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
					LOGGER.info("Applied metadata with {} raw files to IO panel", rawFiles.length);
				}
			}

			// Apply only essential parameters from saved settings
			// Thread count is the only parameter that should be loaded from saved
			// parameters
			if (threadsComboBox != null) {
				threadsComboBox.setValue(params.getThreadCount());
			}

			// Apply Sage-specific parameters
			if (searchTypeGroup != null) {
				// Set search type based on tolerance units
				if ("Da".equals(params.getPrecursorToleranceUnit())) {
					// Open search
					RadioButton openSearchRadio = (RadioButton) searchTypeGroup.getToggles().get(1);
					openSearchRadio.setSelected(true);

					// Show open search parameters panel, hide closed search panel
					if (closedSearchParamsPanel != null) {
						closedSearchParamsPanel.setVisible(false);
						closedSearchParamsPanel.setManaged(false);
					}
					if (openSearchParamsPanel != null) {
						openSearchParamsPanel.setVisible(true);
						openSearchParamsPanel.setManaged(true);
					}

					// Set open search tolerance values using the min/max values from parameters
					if (precursorToleranceMinSpinner != null) {
						precursorToleranceMinSpinner.getValueFactory().setValue(params.getPrecursorToleranceMin());
					}
					if (precursorToleranceMaxSpinner != null) {
						precursorToleranceMaxSpinner.getValueFactory().setValue(params.getPrecursorToleranceMax());
					}
					// Fragment tolerance uses the same spinner for both open and closed search
					if (fragmentToleranceSpinner != null) {
						fragmentToleranceSpinner.getValueFactory().setValue(params.getFragmentTolerance());
					}
					if (precursorToleranceSpinner != null) {
						precursorToleranceSpinner.getValueFactory().setValue(20.0);
					}
				} else {
					// Closed search
					RadioButton closedSearchRadio = (RadioButton) searchTypeGroup.getToggles().get(0);
					closedSearchRadio.setSelected(true);

					// Show closed search parameters panel, hide open search panel
					if (closedSearchParamsPanel != null) {
						closedSearchParamsPanel.setVisible(true);
						closedSearchParamsPanel.setManaged(true);
					}
					if (openSearchParamsPanel != null) {
						openSearchParamsPanel.setVisible(false);
						openSearchParamsPanel.setManaged(false);
					}

					if (precursorToleranceMinSpinner != null) {
						precursorToleranceMinSpinner.getValueFactory().setValue(-500.0);
					}
				if (precursorToleranceMaxSpinner != null) {
					precursorToleranceMaxSpinner.getValueFactory().setValue(150.0);
				}
				// Fragment tolerance uses the same spinner for both open and closed search
				if (fragmentToleranceSpinner != null) {
					fragmentToleranceSpinner.getValueFactory().setValue(params.getFragmentTolerance());
				}
				// Sync the open search fragment tolerance spinner
				if (openFragmentToleranceSpinner != null) {
					openFragmentToleranceSpinner.getValueFactory().setValue(params.getFragmentTolerance());
				}
				if (precursorToleranceSpinner != null) {
					precursorToleranceSpinner.getValueFactory().setValue(params.getPrecursorToleranceMax());
				}
				}
			}

			// Apply database settings using reusable component
			if (magDatabasePanel != null) {
				if (params.getMagDbCatalogueId() != null) {
					// Find and select the matching MAG database
					for (MagDbItem item : magDatabasePanel.getAvailableMagDatabases()) {
						if (item.getCatalogueID().equals(params.getMagDbCatalogueId())) {
							magDatabasePanel.setSelectedMagDatabase(item);

							// Also apply version selection if available
							if (params.getMagDbVersion() != null) {
								for (MagVersion version : magDatabasePanel.getSelectedMagDatabase().getVersionList()) {
									if (version.getCatalogueVersion().equals(params.getMagDbVersion())) {
										magDatabasePanel.setSelectedMagVersion(version);
										LOGGER.info("Applied version selection: {}", version.getCatalogueVersion());
										break;
									}
								}
							}
							break;
						}
					}
				}
				magDatabasePanel.setAppendHostDb(params.isAppendHostDb());
				if (params.getHostDB() != null) {
					magDatabasePanel.setHostDatabase(params.getHostDB());
				}
			}

			// Apply modifications (convert from standard format to Sage format)
			if (params.getFixedModifications() != null && modificationsPanel != null) {
				// Apply saved fixed modifications
				for (String mod : params.getFixedModifications()) {
					if (mod.contains(":")) {
						String[] parts = mod.split(":", 2);
						if (parts.length == 2) {
							String site = parts[0];
							try {
								double mass = Double.parseDouble(parts[1]);
								// Find and update the corresponding row
								for (ModificationsPanel.FixedModificationRow row : modificationsPanel
										.getFixedModData()) {
									if (row.getSite().equals(site)) {
										row.setDeltaMass(mass);
										break;
									}
								}
							} catch (NumberFormatException e) {
								LOGGER.warn("Invalid mass format in fixed modification: {}", mod);
							}
						}
					}
				}
				LOGGER.info("Applied {} fixed modifications", params.getFixedModifications().size());
			}

			if (params.getVariableModifications() != null && modificationsPanel != null) {
				// Apply saved variable modifications
				for (String mod : params.getVariableModifications()) {
					if (mod.contains(":")) {
						String[] parts = mod.split(":", 3);
						if (parts.length == 3) {
							String modName = parts[0];
							String sites = parts[1];
							try {
								double mass = Double.parseDouble(parts[2]);
								// Find and update the corresponding row
								for (ModificationsPanel.VariableModificationRow row : modificationsPanel
										.getVariableModData()) {
									if (row.getModName().equals(modName)) {
										row.setUse(true);
										row.setSites(sites);
										row.setDeltaMass(mass);
										break;
									}
								}
							} catch (NumberFormatException e) {
								LOGGER.warn("Invalid mass format in variable modification: {}", mod);
							}
						}
					}
				}
				LOGGER.info("Applied {} variable modifications", params.getVariableModifications().size());
			}

			// Apply metadata table data if available - convert String arrays back to JavaFX rows
			if (params.getMetadataTableData() != null && !params.getMetadataTableData().isEmpty()
					&& metadataPanel != null) {
				List<uk.ac.quadram.metapilot.task.dda.igc.gui.MetaPilotDdaIgcMaxQuantPanel.MetadataRow> metadataRows = new ArrayList<>();
				for (String[] data : params.getMetadataTableData()) {
					metadataRows.add(uk.ac.quadram.metapilot.task.dda.igc.gui.MetaPilotDdaIgcMaxQuantPanel.MetadataRow.fromArray(data));
				}
				metadataPanel.setMetadataTableData(metadataRows);
				LOGGER.info("Applied {} metadata table rows", params.getMetadataTableData().size());
			}

			// Update parameter view after applying parameters
			Platform.runLater(() -> updateParameterView());

			consolePanel.appendConsoleText("✓ Successfully applied parameters to UI\n");
			LOGGER.info("Successfully applied parameters to UI for MAG Sage workflow");

		} catch (Exception e) {
			String errorMsg = "ERROR: Failed to apply parameters to UI: " + e.getMessage();
			consolePanel.appendConsoleText(errorMsg + "\n");
			LOGGER.error("Failed to apply parameters to UI for MAG Sage workflow", e);
		}
	}

	/**
	 * Load analysis parameters from a specific result directory
	 */
	public void loadAnalysisParameters(String resultDirectory) {
		try {
			LOGGER.info("Loading MAG Sage analysis parameters from: {}", resultDirectory);
			DDAMagSageParameters params = (DDAMagSageParameters) parameterIO.loadAnalysisParameters(resultDirectory);
			if (params != null) {
				applyParametersToUI(params);
				LOGGER.info("Loaded MAG Sage analysis parameters from {}", resultDirectory);
				consolePanel.appendConsoleText("✓ Loaded MAG Sage analysis parameters from " + resultDirectory + "\n");
			} else {
				LOGGER.warn("No parameters found in {}", resultDirectory);
				consolePanel.appendConsoleText("WARNING: No parameters found in " + resultDirectory + "\n");
			}
		} catch (Exception e) {
			LOGGER.error("Failed to load MAG Sage analysis parameters from {}", resultDirectory, e);
			consolePanel.appendConsoleText("WARNING: Failed to load parameters from " + resultDirectory + ": " + e.getMessage() + "\n");
		}
	}

	/**
	 * Get the parameter class for this workflow
	 */
	@Override
	protected Class<DDAMagSageParameters> getParameterClass() {
		return DDAMagSageParameters.class;
	}

}
