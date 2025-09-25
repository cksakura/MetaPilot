package uk.ac.quadram.metapilot.task.dda.igc.gui;

import javafx.application.Platform;
import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.control.*;
import javafx.scene.layout.*;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.scene.text.Font;
import javafx.scene.text.FontWeight;
import javafx.scene.paint.Color;
import javafx.scene.control.cell.PropertyValueFactory;
import javafx.scene.control.cell.TextFieldTableCell;
import javafx.util.converter.DefaultStringConverter;
import javafx.beans.property.SimpleStringProperty;
import javafx.beans.property.StringProperty;
import javafx.scene.control.Alert.AlertType;
import uk.ac.quadram.metapilot.task.common.gui.MetaPilotIOPanelFX;
import uk.ac.quadram.metapilot.task.common.gui.MetaPilotIOPanelFX.ExperimentRow;
import uk.ac.quadram.metapilot.task.common.gui.DdaIgcInputOutputPanel;
import uk.ac.quadram.metapilot.task.common.gui.ConsolePanel;
import uk.ac.quadram.metapilot.task.par.MetaData;
import uk.ac.quadram.metapilot.task.par.MetaConstants;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.MetaPilotApp;
import uk.ac.quadram.metapilot.core.mod.IsobaricTag;
import uk.ac.quadram.metapilot.task.dda.igc.DDAIgcMaxQuantTask;
import uk.ac.quadram.metapilot.task.dda.igc.par.DDAIgcMaxQuantParameters;
import uk.ac.quadram.metapilot.task.dda.igc.par.DDAIgcMaxQuantParameterIO;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Optional;
import java.io.InputStream;
import java.io.InputStreamReader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class MetaPilotDdaIgcMaxQuantPanel extends VBox {
    private static final Logger LOGGER = LogManager.getLogger(MetaPilotDdaIgcMaxQuantPanel.class);
    
    // UI Components
    private DdaIgcInputOutputPanel ioPanel;
    private ConsolePanel consolePanel;
    private ToggleGroup quantModeGroup;
    
    // Parameter Panel Components
    private ComboBox<String> enzymeComboBox;
    private Spinner<Integer> missedCleavagesSpinner;
    private ComboBox<String> digestionModeComboBox;
    private ListView<String> fixedModListView;
    private ListView<String> variableModListView;
    private ComboBox<String> modificationComboBox;
    private Button addFixedModButton;
    private Button addVariableModButton;
    private Button removeFixedModButton;
    private Button removeVariableModButton;
    
    // Additional parameters
    private Spinner<Integer> threadsSpinner;
    private ComboBox<String> ms2ScanModeComboBox;
    private CheckBox generateSampleSpecificDBCheckBox;
    private CheckBox spectraClusteringCheckBox;
    
    // Data
    private ObservableList<String> enzymes;
    private ObservableList<String> modifications;
    private ObservableList<String> fixedModifications;
    private ObservableList<String> variableModifications;
    // Labeling
    private ObservableList<String> labelModifications = FXCollections.observableArrayList();
    private ObservableList<String> isobaricModifications = FXCollections.observableArrayList();
    private ObservableList<String> isobaricBaseTypes = FXCollections.observableArrayList();
    
    // Isotopic labeling components
    private ListView<String> lightLabelListView;
    private ListView<String> mediumLabelListView;
    private ListView<String> heavyLabelListView;
    private ComboBox<String> labelModificationComboBox;
    private ObservableList<String> lightLabels = FXCollections.observableArrayList();
    private ObservableList<String> mediumLabels = FXCollections.observableArrayList();
    private ObservableList<String> heavyLabels = FXCollections.observableArrayList();
    
    // Isobaric labeling components
    private ComboBox<String> isobaricTypeComboBox;
    private TableView<IsobaricLabelRow> isobaricTable;
    private ObservableList<IsobaricLabelRow> isobaricTableData = FXCollections.observableArrayList();
    
    
    // Metadata configuration components
    private TableView<MetadataRow> metadataTable;
    private ObservableList<MetadataRow> metadataTableData = FXCollections.observableArrayList();
    private ComboBox<String> metadataColumnComboBox;
    private Button addMetadataColumnButton;
    private Button removeMetadataColumnButton;
    private Button resetMetadataButton;
    
    
    // Add parameter IO instance
    private DDAIgcMaxQuantParameterIO parameterIO;
    
    // Current running task
    private DDAIgcMaxQuantTask currentTask;
    
    public MetaPilotDdaIgcMaxQuantPanel() {
        // Initialize parameter IO instance
        parameterIO = new DDAIgcMaxQuantParameterIO();
        
        // Initialize unified IO panel
        ioPanel = new DdaIgcInputOutputPanel(this::updateParameterView);
        
        // Initialize console panel
        consolePanel = new ConsolePanel();
        consolePanel.setStartCallback(this::startAnalysis);
        consolePanel.setStopCallback(this::stopAnalysis);
        consolePanel.setViewReportCallback(this::viewReport);
        
        initializeUI();
        
        // Note: Parameters are only saved when starting analysis or explicitly requested
        
        // Load parameters after UI is fully initialized
        Platform.runLater(() -> {
            loadLastUsedParameters();
            
            // Add window close event handler to save parameters
            setupWindowCloseHandler();
        });
    }
    
    private void initializeUI() {
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
        
        Tab consoleTabPane = new Tab("Console");
        consoleTabPane.setContent(consolePanel);
        
        tabPane.getTabs().addAll(ioTab, paramTab, consoleTabPane);
        
        // Add tab change listener to update metadata when Parameters tab is selected
        // and save parameters when Console tab is selected
        tabPane.getSelectionModel().selectedItemProperty().addListener((obs, oldTab, newTab) -> {
            if (newTab == paramTab) {
                updateMetadataTableFromIOPanel();
            } else if (newTab == consoleTabPane) {
                // Save current parameters when Console tab is selected
                try {
                    DDAIgcMaxQuantParameters currentParams = collectParametersFromUI();
                    if (currentParams != null) {
                        parameterIO.saveLastUsedParameters(currentParams);
                        parameterIO.savePreferences(currentParams);
                        LOGGER.info("Saved current parameters when Console tab was selected");
                    }
                } catch (Exception e) {
                    LOGGER.warn("Failed to save parameters when Console tab was selected: " + e.getMessage());
                }
            }
        });
        
        // Add tab pane to main panel
        getChildren().addAll(tabPane);
    }
    
    private int ProcessorCount() {
		int total = 0;
		try {
			String[] args = { "wmic", "cpu", "get", "NumberOfLogicalProcessors" };
			ProcessBuilder pb = new ProcessBuilder(args);
			Process p = pb.start();
			BufferedInputStream in = new BufferedInputStream(p.getInputStream());
			BufferedReader inBr = new BufferedReader(new InputStreamReader(in));
			String lineStr = null;
			L: while ((lineStr = inBr.readLine()) != null) {
				String trimLine = lineStr.trim();
				if (trimLine.length() > 0) {
					int oneCore = 0;
					for (int i = 0; i < trimLine.length(); i++) {
						char aa = trimLine.charAt(i);
						if (aa >= '0' && aa <= '9') {
							oneCore += ((int) (aa - '0')) * Math.pow(10, (trimLine.length() - i - 1));
						} else {
							continue L;
						}
					}
					total += oneCore;
				}
			}
			if (p.waitFor() != 0) {
				if (p.exitValue() == 1)
					System.err.println("false");
			}
			inBr.close();
			in.close();
		} catch (IOException | InterruptedException e) {
			System.err.println("Failed to get logical processors using wmic: " + e.getMessage());
		}

		if (total == 0) {
			total = Runtime.getRuntime().availableProcessors();
		}
		return total;
	}
    
    private VBox createIOTab() {
        VBox tabContent = new VBox(10);
        tabContent.setPadding(new Insets(10));
        
        // Add unified IO Panel (includes both file selection and database settings)
        tabContent.getChildren().add(ioPanel);
        
        return tabContent;
    }
    
    
    
    private VBox createParametersTab() {
        VBox tabContent = new VBox(10);
        tabContent.setPadding(new Insets(10));
        tabContent.setMinWidth(960);
        tabContent.setMinHeight(760);
        
        // Initialize data
        enzymes = FXCollections.observableArrayList();
        modifications = FXCollections.observableArrayList();
        fixedModifications = FXCollections.observableArrayList();
        variableModifications = FXCollections.observableArrayList();
        // Set default modifications
        fixedModifications.add("Carbamidomethyl (C)");
        variableModifications.addAll("Acetyl (Protein N-term)", "Oxidation (M)");
        
        // Load enzymes and modifications from XML
        loadEnzymesAndModifications();
        
        
        // Create main grid layout for parameters with 2-column responsive sizing
        GridPane mainGrid = new GridPane();
        mainGrid.setHgap(15);
        mainGrid.setVgap(15);
        mainGrid.setPadding(new Insets(10));
        mainGrid.setMinWidth(940);
        mainGrid.setMinHeight(740);
        VBox.setVgrow(mainGrid, Priority.ALWAYS);
        
        // Column constraints for 2-column layout (50% + 50%) with minimum widths
        ColumnConstraints col1 = new ColumnConstraints();
        col1.setPercentWidth(50);
        col1.setMinWidth(450);
        col1.setHgrow(Priority.ALWAYS);
        ColumnConstraints col2 = new ColumnConstraints();
        col2.setPercentWidth(50);
        col2.setMinWidth(450);
        col2.setHgrow(Priority.ALWAYS);
        mainGrid.getColumnConstraints().addAll(col1, col2);
        
        // Row constraints for proper height distribution
        RowConstraints row1 = new RowConstraints();
        row1.setVgrow(Priority.ALWAYS);
        row1.setMinHeight(400);
        mainGrid.getRowConstraints().add(row1);
        
        // Left column: Enzyme Settings, Additional Parameters, Taxonomic Analysis, Modifications (stacked)
        VBox leftColumn = new VBox(15);
        leftColumn.setMinWidth(430);
        VBox enzymeSection = createEnzymeSection();
        VBox additionalParamsSection = createAdditionalParametersSection();
        VBox modificationsSection = createModificationsSection();
        leftColumn.getChildren().addAll(enzymeSection, additionalParamsSection, modificationsSection);
        GridPane.setVgrow(leftColumn, Priority.ALWAYS);
        
        // Right column: Quantitative Mode and Metadata Configuration (stacked)
        VBox rightColumn = new VBox(15);
        rightColumn.setMinWidth(430);
        VBox quantModeSection = createQuantitativeModeSection();
        VBox metadataSection = createMetadataConfigurationSection();
        rightColumn.getChildren().addAll(quantModeSection, metadataSection);
        GridPane.setVgrow(rightColumn, Priority.ALWAYS);
        
        // Add columns to grid
        mainGrid.add(leftColumn, 0, 0);
        mainGrid.add(rightColumn, 1, 0);
        
        tabContent.getChildren().add(mainGrid);
        return tabContent;
    }
    
    private VBox createEnzymeSection() {
        VBox section = new VBox(10);
        section.setPadding(new Insets(10));
        section.setStyle("-fx-background-color: #f5f5f5; -fx-background-radius: 5;");
        
        Label titleLabel = createSectionTitle("Enzyme Settings");
        
        // All enzyme parameters in one horizontal line with responsive sizing
        HBox enzymeParamsBox = new HBox(15);
        enzymeParamsBox.setAlignment(Pos.CENTER_LEFT);
        
        // Enzyme Selection
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
        
        // Missed Cleavages
        VBox missedCleavagesBox = new VBox(5);
        missedCleavagesBox.setMinWidth(100);
        
        Label missedLabel = new Label("Missed Cleavages:");
        missedLabel.setStyle("-fx-font-weight: bold;");
        missedCleavagesBox.getChildren().add(missedLabel);
        
        missedCleavagesSpinner = new Spinner<>(0, 5, 2);
        missedCleavagesSpinner.setMinWidth(100);
        missedCleavagesSpinner.setPrefWidth(100);
        missedCleavagesSpinner.setMaxWidth(Double.MAX_VALUE);
        missedCleavagesSpinner.valueProperty().addListener((obs, oldVal, newVal) -> updateParameterView());
        missedCleavagesBox.getChildren().add(missedCleavagesSpinner);
        
        // Digestion Mode
        VBox digestionModeBox = new VBox(5);
        digestionModeBox.setMinWidth(140);
        HBox.setHgrow(digestionModeBox, Priority.ALWAYS);
        
        Label digestionLabel = new Label("Digestion Mode:");
        digestionLabel.setStyle("-fx-font-weight: bold;");
        digestionModeBox.getChildren().add(digestionLabel);
        
        digestionModeComboBox = new ComboBox<>(FXCollections.observableArrayList(
            "Specific", "Semi free N-term", "Semi free C-term", "Semi specific", "Unspecific", "None"));
        digestionModeComboBox.setMinWidth(140);
        digestionModeComboBox.setPrefWidth(160);
        digestionModeComboBox.setMaxWidth(Double.MAX_VALUE);
        digestionModeComboBox.setValue("Specific"); // Set default value
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
        
        // Create a grid for better organization and responsive layout
        GridPane paramGrid = new GridPane();
        paramGrid.setHgap(15);
        paramGrid.setVgap(8);
        
        // Column constraints for responsive sizing
        ColumnConstraints col1 = new ColumnConstraints();
        col1.setMinWidth(150);
        col1.setHgrow(Priority.ALWAYS);
        ColumnConstraints col2 = new ColumnConstraints();
        col2.setMinWidth(150);
        col2.setHgrow(Priority.ALWAYS);
        paramGrid.getColumnConstraints().addAll(col1, col2);
        
        // Threads
        Label threadsLabel = new Label("Threads:");
        threadsLabel.setStyle("-fx-font-weight: bold;");
        threadsSpinner = new Spinner<>(1, ProcessorCount(), 4);
        threadsSpinner.setMinWidth(100);
        threadsSpinner.setPrefWidth(100);
        threadsSpinner.setMaxWidth(Double.MAX_VALUE);
        threadsSpinner.valueProperty().addListener((obs, oldVal, newVal) -> updateParameterView());
        
        VBox threadsBox = new VBox(5);
        threadsBox.getChildren().addAll(threadsLabel, threadsSpinner);
        
        // MS2 Scan Mode
        Label ms2Label = new Label("MS2 Scan Mode:");
        ms2Label.setStyle("-fx-font-weight: bold;");
        ms2ScanModeComboBox = new ComboBox<>(FXCollections.observableArrayList("FTMS", "ITMS"));
        ms2ScanModeComboBox.setValue("FTMS");
        ms2ScanModeComboBox.setMinWidth(100);
        ms2ScanModeComboBox.setPrefWidth(120);
        ms2ScanModeComboBox.setMaxWidth(Double.MAX_VALUE);
        ms2ScanModeComboBox.setOnAction(e -> updateParameterView());
        
        VBox ms2Box = new VBox(5);
        ms2Box.getChildren().addAll(ms2Label, ms2ScanModeComboBox);
        
        // Add to grid
        paramGrid.add(threadsBox, 0, 0);
        paramGrid.add(ms2Box, 1, 0);
        
        // Checkboxes with proper spacing
        VBox checkBoxContainer = new VBox(8);
        
        generateSampleSpecificDBCheckBox = new CheckBox("Generate sample specific database");
        generateSampleSpecificDBCheckBox.setSelected(true);
        generateSampleSpecificDBCheckBox.setOnAction(e -> updateParameterView());
        
        spectraClusteringCheckBox = new CheckBox("Spectra clustering in first search");
        spectraClusteringCheckBox.setSelected(false);
        spectraClusteringCheckBox.setOnAction(e -> updateParameterView());
        
        checkBoxContainer.getChildren().addAll(generateSampleSpecificDBCheckBox, spectraClusteringCheckBox);
        
        section.getChildren().addAll(titleLabel, paramGrid, checkBoxContainer);
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
        
        // Radio buttons for quantitative modes
        RadioButton labelFreeRadio = new RadioButton("Label-free");
        labelFreeRadio.setToggleGroup(quantModeGroup);
        labelFreeRadio.setSelected(true);
        
        RadioButton isotopicRadio = new RadioButton("Isotopic");
        isotopicRadio.setToggleGroup(quantModeGroup);
        
        RadioButton isobaricRadio = new RadioButton("Isobaric");
        isobaricRadio.setToggleGroup(quantModeGroup);
        
        HBox radioBox = new HBox(15);
        radioBox.setAlignment(Pos.CENTER_LEFT);
        radioBox.getChildren().addAll(labelFreeRadio, isotopicRadio, isobaricRadio);
        
        // StackPane for overlapping panels - both panels must occupy same space
        StackPane labelingStackPane = new StackPane();
        labelingStackPane.setMinHeight(200);
        labelingStackPane.setMaxHeight(Double.MAX_VALUE);
        VBox.setVgrow(labelingStackPane, Priority.ALWAYS);
        
        // Create both panels with identical sizing constraints
        VBox isotopicPanel = createCompactIsotopicLabelingPanel();
        VBox isobaricPanel = createCompactIsobaricLabelingPanel();
        
        // Ensure both panels have identical size constraints for proper StackPane behavior
        isotopicPanel.setMinHeight(200);
        isotopicPanel.setPrefHeight(220);
        isotopicPanel.setMaxHeight(Double.MAX_VALUE);
        
        isobaricPanel.setMinHeight(200);
        isobaricPanel.setPrefHeight(220);
        isobaricPanel.setMaxHeight(Double.MAX_VALUE);
        
        // Initially both panels are invisible
        isotopicPanel.setVisible(false);
        isotopicPanel.setManaged(false); // Key: don't take up space when hidden
        isobaricPanel.setVisible(false);
        isobaricPanel.setManaged(false); // Key: don't take up space when hidden
        
        labelingStackPane.getChildren().addAll(isotopicPanel, isobaricPanel);
        
        // Listeners to show/hide labeling panels with proper space management
        labelFreeRadio.setOnAction(e -> {
            isotopicPanel.setVisible(false);
            isotopicPanel.setManaged(false);
            isobaricPanel.setVisible(false);
            isobaricPanel.setManaged(false);
            updateParameterView();
            updateMetadataTable();
        });
        
        isotopicRadio.setOnAction(e -> {
            isotopicPanel.setVisible(true);
            isotopicPanel.setManaged(true);
            isobaricPanel.setVisible(false);
            isobaricPanel.setManaged(false);
            updateParameterView();
            updateMetadataTable();
        });
        
        isobaricRadio.setOnAction(e -> {
            isotopicPanel.setVisible(false);
            isotopicPanel.setManaged(false);
            isobaricPanel.setVisible(true);
            isobaricPanel.setManaged(true);
            updateParameterView();
            updateMetadataTable();
        });
        
        section.getChildren().addAll(titleLabel, radioBox, labelingStackPane);
        return section;
    }
    
    private VBox createCompactIsotopicLabelingPanel() {
        VBox panel = new VBox(8);
        panel.setPadding(new Insets(8));
        panel.setStyle("-fx-background-color: #e8f4f8; -fx-background-radius: 5;");
        
        // Label selection row
        HBox labelSelectionBox = new HBox(8);
        labelSelectionBox.setAlignment(Pos.CENTER_LEFT);
        
        Label selectLabel = new Label("Select label:");
        selectLabel.setMinWidth(120);
        selectLabel.setPrefWidth(120);
        labelModificationComboBox = new ComboBox<>(labelModifications);
        labelModificationComboBox.setMinWidth(120);
        labelModificationComboBox.setPrefWidth(180);
        labelModificationComboBox.setMaxWidth(Double.MAX_VALUE);
        labelModificationComboBox.setPromptText("Select label");
        HBox.setHgrow(labelModificationComboBox, Priority.ALWAYS);
        
        Button addLightBtn = new Button("Add Light");
        addLightBtn.setMinWidth(100);
        addLightBtn.setPrefWidth(120);
        addLightBtn.setOnAction(e -> addLightLabel());
        
        Button addMediumBtn = new Button("Add Medium");
        addMediumBtn.setMinWidth(100);
        addMediumBtn.setPrefWidth(120);
        addMediumBtn.setOnAction(e -> addMediumLabel());
        
        Button addHeavyBtn = new Button("Add Heavy");
        addHeavyBtn.setMinWidth(100);
        addHeavyBtn.setPrefWidth(120);
        addHeavyBtn.setOnAction(e -> addHeavyLabel());
        
        labelSelectionBox.getChildren().addAll(selectLabel, labelModificationComboBox, addLightBtn, addMediumBtn, addHeavyBtn);
        
        // Lists with remove buttons - use VBox.setVgrow to fill available space
        HBox listsBox = new HBox(8);
        listsBox.setAlignment(Pos.TOP_LEFT);
        VBox.setVgrow(listsBox, Priority.ALWAYS); // Allow lists box to grow
        
        // Light labels section
        VBox lightBox = new VBox(5);
        lightBox.setMinWidth(120);
        HBox.setHgrow(lightBox, Priority.ALWAYS);
        
        Label lightLabel = new Label("Light Labels:");
        lightLabel.setStyle("-fx-font-weight: bold;");
        
        lightLabelListView = new ListView<>(lightLabels);
        lightLabelListView.setMinHeight(80);
        lightLabelListView.setPrefHeight(100);
        lightLabelListView.setMaxHeight(Double.MAX_VALUE); // Allow to grow
        VBox.setVgrow(lightLabelListView, Priority.ALWAYS); // Fill available space
        
        Button removeLightBtn = new Button("Remove Selected");
        removeLightBtn.setMinWidth(120);
        removeLightBtn.setPrefWidth(120);
        removeLightBtn.setMaxWidth(Double.MAX_VALUE);
        removeLightBtn.setOnAction(e -> removeLightLabel());
        
        lightBox.getChildren().addAll(lightLabel, lightLabelListView, removeLightBtn);
        
        // Medium labels section
        VBox mediumBox = new VBox(5);
        mediumBox.setMinWidth(120);
        HBox.setHgrow(mediumBox, Priority.ALWAYS);
        
        Label mediumLabel = new Label("Medium Labels:");
        mediumLabel.setStyle("-fx-font-weight: bold;");
        
        mediumLabelListView = new ListView<>(mediumLabels);
        mediumLabelListView.setMinHeight(80);
        mediumLabelListView.setPrefHeight(100);
        mediumLabelListView.setMaxHeight(Double.MAX_VALUE); // Allow to grow
        VBox.setVgrow(mediumLabelListView, Priority.ALWAYS); // Fill available space
        
        Button removeMediumBtn = new Button("Remove Selected");
        removeMediumBtn.setMinWidth(120);
        removeMediumBtn.setPrefWidth(120);
        removeMediumBtn.setMaxWidth(Double.MAX_VALUE);
        removeMediumBtn.setOnAction(e -> removeMediumLabel());
        
        mediumBox.getChildren().addAll(mediumLabel, mediumLabelListView, removeMediumBtn);
        
        // Heavy labels section
        VBox heavyBox = new VBox(5);
        heavyBox.setMinWidth(120);
        HBox.setHgrow(heavyBox, Priority.ALWAYS);
        
        Label heavyLabel = new Label("Heavy Labels:");
        heavyLabel.setStyle("-fx-font-weight: bold;");
        
        heavyLabelListView = new ListView<>(heavyLabels);
        heavyLabelListView.setMinHeight(80);
        heavyLabelListView.setPrefHeight(100);
        heavyLabelListView.setMaxHeight(Double.MAX_VALUE); // Allow to grow
        VBox.setVgrow(heavyLabelListView, Priority.ALWAYS); // Fill available space
        
        Button removeHeavyBtn = new Button("Remove Selected");
        removeHeavyBtn.setMinWidth(120);
        removeHeavyBtn.setPrefWidth(120);
        removeHeavyBtn.setMaxWidth(Double.MAX_VALUE);
        removeHeavyBtn.setOnAction(e -> removeHeavyLabel());
        
        heavyBox.getChildren().addAll(heavyLabel, heavyLabelListView, removeHeavyBtn);
        
        listsBox.getChildren().addAll(lightBox, mediumBox, heavyBox);
        
        panel.getChildren().addAll(labelSelectionBox, listsBox);
        
        return panel;
    }
    
    private VBox createCompactIsobaricLabelingPanel() {
        VBox panel = new VBox(8);
        panel.setPadding(new Insets(8));
        panel.setStyle("-fx-background-color: #f8e8f4; -fx-background-radius: 5;");
        // Remove fixed height constraint to allow filling
        
        // Isobaric type selection
        HBox typeBox = new HBox(10);
        typeBox.setAlignment(Pos.CENTER_LEFT);
        
        Label typeLabel = new Label("Isobaric Type:");
        typeLabel.setStyle("-fx-font-weight: bold;");
        
        isobaricTypeComboBox = new ComboBox<>(isobaricBaseTypes);
        isobaricTypeComboBox.setMinWidth(200);
        isobaricTypeComboBox.setPrefWidth(220);
        isobaricTypeComboBox.setMaxWidth(Double.MAX_VALUE);
        isobaricTypeComboBox.setOnAction(e -> updateIsobaricTable());
        HBox.setHgrow(isobaricTypeComboBox, Priority.ALWAYS);
        
        typeBox.getChildren().addAll(typeLabel, isobaricTypeComboBox);
        
        // Isobaric table - allow it to fill available space
        isobaricTable = new TableView<>(isobaricTableData);
        isobaricTable.setMinHeight(120);
        isobaricTable.setPrefHeight(140);
        isobaricTable.setMaxHeight(Double.MAX_VALUE); // Allow to grow
        isobaricTable.setEditable(true);
        isobaricTable.setColumnResizePolicy(TableView.CONSTRAINED_RESIZE_POLICY);
        VBox.setVgrow(isobaricTable, Priority.ALWAYS); // Fill available space
        
        // Table columns with responsive widths
        TableColumn<IsobaricLabelRow, String> labelCol = new TableColumn<>("Label");
        labelCol.setCellValueFactory(new PropertyValueFactory<>("labelName"));
        labelCol.setMinWidth(120);
        labelCol.setPrefWidth(150);
        
        TableColumn<IsobaricLabelRow, String> minus2Col = new TableColumn<>("-2");
        minus2Col.setCellValueFactory(new PropertyValueFactory<>("minus2"));
        minus2Col.setCellFactory(TextFieldTableCell.forTableColumn(new DefaultStringConverter()));
        minus2Col.setOnEditCommit(e -> e.getRowValue().setMinus2(e.getNewValue()));
        minus2Col.setMinWidth(50);
        minus2Col.setPrefWidth(60);
        
        TableColumn<IsobaricLabelRow, String> minus1Col = new TableColumn<>("-1");
        minus1Col.setCellValueFactory(new PropertyValueFactory<>("minus1"));
        minus1Col.setCellFactory(TextFieldTableCell.forTableColumn(new DefaultStringConverter()));
        minus1Col.setOnEditCommit(e -> e.getRowValue().setMinus1(e.getNewValue()));
        minus1Col.setMinWidth(50);
        minus1Col.setPrefWidth(60);
        
        TableColumn<IsobaricLabelRow, String> plus1Col = new TableColumn<>("+1");
        plus1Col.setCellValueFactory(new PropertyValueFactory<>("plus1"));
        plus1Col.setCellFactory(TextFieldTableCell.forTableColumn(new DefaultStringConverter()));
        plus1Col.setOnEditCommit(e -> e.getRowValue().setPlus1(e.getNewValue()));
        plus1Col.setMinWidth(50);
        plus1Col.setPrefWidth(60);
        
        TableColumn<IsobaricLabelRow, String> plus2Col = new TableColumn<>("+2");
        plus2Col.setCellValueFactory(new PropertyValueFactory<>("plus2"));
        plus2Col.setCellFactory(TextFieldTableCell.forTableColumn(new DefaultStringConverter()));
        plus2Col.setOnEditCommit(e -> e.getRowValue().setPlus2(e.getNewValue()));
        plus2Col.setMinWidth(50);
        plus2Col.setPrefWidth(60);
        
        isobaricTable.getColumns().addAll(labelCol, minus2Col, minus1Col, plus1Col, plus2Col);
        
        panel.getChildren().addAll(typeBox, isobaricTable);
        return panel;
    }
    

    
    private VBox createMetadataConfigurationSection() {
        VBox section = new VBox(10);
        section.setPadding(new Insets(10));
        section.setStyle("-fx-background-color: #fff8f0; -fx-background-radius: 5;");
        section.setMinHeight(200);
        VBox.setVgrow(section, Priority.ALWAYS);
        
        Label titleLabel = createSectionTitle("Metadata Configuration");
        
        // Dynamic content based on quantitative mode
        Label infoLabel = new Label("Metadata will be configured based on selected quantitative mode");
        infoLabel.setStyle("-fx-font-style: italic; -fx-text-fill: #666666;");
        
        // Metadata table
        metadataTable = new TableView<>(metadataTableData);
        metadataTable.setMinHeight(120);
        metadataTable.setPrefHeight(150);
        metadataTable.setMaxHeight(Double.MAX_VALUE);
        metadataTable.setEditable(true);
        metadataTable.setColumnResizePolicy(TableView.CONSTRAINED_RESIZE_POLICY);
        VBox.setVgrow(metadataTable, Priority.ALWAYS);
        
        // Initially create basic columns
        createMetadataTableColumns();
        
        // Control buttons
        HBox controlBox = new HBox(8);
        controlBox.setAlignment(Pos.CENTER_LEFT);
        
        metadataColumnComboBox = new ComboBox<>(FXCollections.observableArrayList(
            "Age", "Gender", "Condition", "Treatment", "Time Point", "BMI", "Custom"));
        metadataColumnComboBox.setMinWidth(100);
        metadataColumnComboBox.setPrefWidth(120);
        metadataColumnComboBox.setPromptText("Add column");
        
        // Custom column TextField (hidden by default)
        TextField customColumnTextField = new TextField();
        customColumnTextField.setPromptText("Enter custom column name");
        customColumnTextField.setMinWidth(120);
        customColumnTextField.setVisible(false);
        customColumnTextField.setManaged(false);
        
        // Show/hide the custom column TextField based on ComboBox selection
        metadataColumnComboBox.setOnAction(e -> {
            String selected = metadataColumnComboBox.getValue();
            boolean isCustom = "Custom".equals(selected);
            customColumnTextField.setVisible(isCustom);
            customColumnTextField.setManaged(isCustom);
            if (isCustom) {
                customColumnTextField.requestFocus();
            }
        });
        
        addMetadataColumnButton = new Button("Add Column");
        addMetadataColumnButton.setMinWidth(120);
        addMetadataColumnButton.setPrefWidth(150);
        addMetadataColumnButton.setOnAction(e -> {
            String columnName = metadataColumnComboBox.getValue();
            if ("Custom".equals(columnName)) {
                String customName = customColumnTextField.getText();
                if (customName != null && !customName.trim().isEmpty()) {
                    addMetadataColumnWithProperty(customName.trim(), getPropertyNameForColumn(customName.trim()));
                    customColumnTextField.clear();
                }
            } else if (columnName != null && !columnName.isEmpty()) {
                addMetadataColumnWithProperty(columnName, getPropertyNameForColumn(columnName));
            }
        });
        
        removeMetadataColumnButton = new Button("Remove Column");
        removeMetadataColumnButton.setMinWidth(120);
        removeMetadataColumnButton.setPrefWidth(150);
        removeMetadataColumnButton.setOnAction(e -> removeMetadataColumn());
        
        resetMetadataButton = new Button("Reset All");
        resetMetadataButton.setMinWidth(120);
        resetMetadataButton.setPrefWidth(150);
        resetMetadataButton.setStyle("-fx-background-color: #ff6b6b; -fx-text-fill: white;");
        resetMetadataButton.setOnAction(e -> resetMetadataTable());
        
        controlBox.getChildren().addAll(metadataColumnComboBox, customColumnTextField, addMetadataColumnButton, removeMetadataColumnButton, resetMetadataButton);
        
        section.getChildren().addAll(titleLabel, infoLabel, metadataTable, controlBox);
        
        // Update metadata table when quantitative mode changes
        updateMetadataTable();
        
        return section;
    }
    
    private void createMetadataTableColumns() {
        metadataTable.getColumns().clear();
        
        // Sample/Experiment column
        TableColumn<MetadataRow, String> sampleCol = new TableColumn<>("Sample");
        sampleCol.setCellValueFactory(new PropertyValueFactory<>("sampleName"));
        sampleCol.setCellFactory(TextFieldTableCell.forTableColumn(new DefaultStringConverter()));
        sampleCol.setOnEditCommit(e -> {
            e.getRowValue().setSampleName(e.getNewValue());
            updateParameterView();
        });
        sampleCol.setMinWidth(100);
        sampleCol.setPrefWidth(120);
        
        metadataTable.getColumns().add(sampleCol);
        
        // Add label type column if using isotopic or isobaric labeling
        if (quantModeGroup.getSelectedToggle() != null) {
            RadioButton selected = (RadioButton) quantModeGroup.getSelectedToggle();
            String mode = selected.getText();
            
            if (mode.contains("Isotopic") || mode.contains("Isobaric")) {
                TableColumn<MetadataRow, String> labelCol = new TableColumn<>("Label Type");
                labelCol.setCellValueFactory(new PropertyValueFactory<>("labelType"));
                labelCol.setMinWidth(120);
                labelCol.setPrefWidth(150);
                metadataTable.getColumns().add(labelCol);
            }
        }
        
        // Add common metadata columns
        addMetadataColumnWithProperty("Age", "age");
        addMetadataColumnWithProperty("Gender", "gender");
        addMetadataColumnWithProperty("Condition", "condition");
    }
    
    private void updateMetadataTable() {
        // Always show label-free rows by default
        metadataTableData.clear();
        List<String> experimentNames = ioPanel.getIoPanel().getExperimentNames();
        if (experimentNames.isEmpty()) return;
        createMetadataTableColumns();
        // Determine quantitative mode
        if (quantModeGroup.getSelectedToggle() != null) {
            RadioButton selected = (RadioButton) quantModeGroup.getSelectedToggle();
            String mode = selected.getText();
            boolean hasLabels = false;
            if (mode.contains("Isotopic")) {
                hasLabels = !lightLabels.isEmpty() || !mediumLabels.isEmpty() || !heavyLabels.isEmpty();
            } else if (mode.contains("Isobaric")) {
                hasLabels = !isobaricTableData.isEmpty();
            }
            if (mode.contains("Isotopic") && hasLabels) {
                for (String expName : experimentNames) {
                    for (String label : lightLabels) {
                        MetadataRow row = new MetadataRow(expName);
                        row.setLabelType("Light - " + label);
                        metadataTableData.add(row);
                    }
                    for (String label : mediumLabels) {
                        MetadataRow row = new MetadataRow(expName);
                        row.setLabelType("Medium - " + label);
                        metadataTableData.add(row);
                    }
                    for (String label : heavyLabels) {
                        MetadataRow row = new MetadataRow(expName);
                        row.setLabelType("Heavy - " + label);
                        metadataTableData.add(row);
                    }
                }
            } else if (mode.contains("Isobaric") && hasLabels) {
                for (String expName : experimentNames) {
                    for (IsobaricLabelRow isoLabel : isobaricTableData) {
                        // Create sample name combining experiment name and label name
                        String sampleName = expName + "_" + isoLabel.getLabelName();
                        MetadataRow row = new MetadataRow(sampleName);
                        row.setLabelType("Isobaric - " + isoLabel.getLabelName());
                        metadataTableData.add(row);
                    }
                }
            } else {
                // Show label-free by default
                for (String expName : experimentNames) {
                    metadataTableData.add(new MetadataRow(expName));
                }
            }
        }
    }

    private void addMetadataColumnWithProperty(String columnName, String propertyName) {
        TableColumn<MetadataRow, String> column = new TableColumn<>(columnName);
        column.setCellValueFactory(new PropertyValueFactory<>(propertyName));
        column.setCellFactory(TextFieldTableCell.forTableColumn(new DefaultStringConverter()));
        column.setOnEditCommit(event -> {
            MetadataRow row = event.getRowValue();
            String newValue = event.getNewValue();
            switch (propertyName) {
                case "age": row.setAge(newValue); break;
                case "gender": row.setGender(newValue); break;
                case "condition": row.setCondition(newValue); break;
                case "treatment": row.setTreatment(newValue); break;
                case "timePoint": row.setTimePoint(newValue); break;
                case "customField1": row.setCustomField1(newValue); break;
                case "customField2": row.setCustomField2(newValue); break;
                case "customField3": row.setCustomField3(newValue); break;
            }
            updateParameterView();
        });
        column.setMinWidth(80);
        column.setPrefWidth(100);
        metadataTable.getColumns().add(column);
    }
    
    private String getPropertyNameForColumn(String columnName) {
        switch (columnName) {
            case "Age": return "age";
            case "Gender": return "gender";
            case "Condition": return "condition";
            case "Treatment": return "treatment";
            case "Time Point": return "timePoint";
            case "BMI": return "customField1";
            default: return "customField" + (metadataTable.getColumns().size());
        }
    }
    
    private void removeMetadataColumn() {
        if (metadataTable.getColumns().size() > 1) { // Keep at least the sample column
            metadataTable.getColumns().remove(metadataTable.getColumns().size() - 1);
        }
    }
    
    private void resetMetadataTable() {
        // Clear all data
        metadataTableData.clear();
        
        // Remove all columns except the Sample column
        metadataTable.getColumns().clear();
        
        // Recreate basic table structure with only Sample column
        createMetadataTableColumns();
        
        // Update parameter view
        updateParameterView();
    }
    
    private VBox createModificationsSection() {
        VBox section = new VBox(10);
        section.setPadding(new Insets(10));
        section.setStyle("-fx-background-color: #f5f5f5; -fx-background-radius: 5;");
        section.setMinHeight(250);
        VBox.setVgrow(section, Priority.ALWAYS);
        
        Label titleLabel = createSectionTitle("Modifications");
        
        GridPane modGrid = new GridPane();
        modGrid.setHgap(15);
        modGrid.setVgap(10);
        
        // Column constraints for responsive sizing
        ColumnConstraints col1 = new ColumnConstraints();
        col1.setMinWidth(200);
        col1.setHgrow(Priority.ALWAYS);
        ColumnConstraints col2 = new ColumnConstraints();
        col2.setMinWidth(200);
        col2.setHgrow(Priority.ALWAYS);
        modGrid.getColumnConstraints().addAll(col1, col2);
        
        // Row constraints for proper height distribution
        RowConstraints row1 = new RowConstraints();
        row1.setVgrow(Priority.NEVER);
        RowConstraints row2 = new RowConstraints();
        row2.setVgrow(Priority.ALWAYS);
        row2.setMinHeight(120);
        RowConstraints row3 = new RowConstraints();
        row3.setVgrow(Priority.NEVER);
        RowConstraints row4 = new RowConstraints();
        row4.setVgrow(Priority.NEVER);
        modGrid.getRowConstraints().addAll(row1, row2, row3, row4);
        
        // Fixed Modifications
        Label fixedModLabel = new Label("Fixed Modifications:");
        fixedModLabel.setStyle("-fx-font-weight: bold;");
        
        fixedModListView = new ListView<>(fixedModifications);
        fixedModListView.setMinHeight(120);
        fixedModListView.setPrefHeight(150);
        fixedModListView.setMaxHeight(Double.MAX_VALUE);
        GridPane.setVgrow(fixedModListView, Priority.ALWAYS);
        
        removeFixedModButton = new Button("Remove Selected");
        removeFixedModButton.setMinWidth(120);
        removeFixedModButton.setPrefWidth(120);
        removeFixedModButton.setMaxWidth(Double.MAX_VALUE);
        removeFixedModButton.setOnAction(e -> removeFixedModification());
        
        // Variable Modifications
        Label variableModLabel = new Label("Variable Modifications:");
        variableModLabel.setStyle("-fx-font-weight: bold;");
        
        variableModListView = new ListView<>(variableModifications);
        variableModListView.setMinHeight(120);
        variableModListView.setPrefHeight(150);
        variableModListView.setMaxHeight(Double.MAX_VALUE);
        GridPane.setVgrow(variableModListView, Priority.ALWAYS);
        
        removeVariableModButton = new Button("Remove Selected");
        removeVariableModButton.setMinWidth(120);
        removeVariableModButton.setPrefWidth(120);
        removeVariableModButton.setMaxWidth(Double.MAX_VALUE);
        removeVariableModButton.setOnAction(e -> removeVariableModification());
        
        // Modification Selection with responsive layout
        VBox modSelectionContainer = new VBox(8);
        
        Label selectionLabel = new Label("Add Modification:");
        selectionLabel.setStyle("-fx-font-weight: bold;");
        
        modificationComboBox = new ComboBox<>(modifications);
        modificationComboBox.setMinWidth(200);
        modificationComboBox.setPrefWidth(250);
        modificationComboBox.setMaxWidth(Double.MAX_VALUE);
        modificationComboBox.setPromptText("Select modification to add");
        
        HBox buttonBox = new HBox(10);
        buttonBox.setAlignment(Pos.CENTER_LEFT);
        
        addFixedModButton = new Button("Add to Fixed");
        addFixedModButton.setMinWidth(150);
        addFixedModButton.setPrefWidth(150);
        addFixedModButton.setOnAction(e -> addFixedModification());
        
        addVariableModButton = new Button("Add to Variable");
        addVariableModButton.setMinWidth(150);
        addVariableModButton.setPrefWidth(150);
        addVariableModButton.setOnAction(e -> addVariableModification());
        
        buttonBox.getChildren().addAll(addFixedModButton, addVariableModButton);
        
        modSelectionContainer.getChildren().addAll(selectionLabel, modificationComboBox, buttonBox);
        
        // Add components to grid
        modGrid.add(fixedModLabel, 0, 0);
        modGrid.add(variableModLabel, 1, 0);
        modGrid.add(fixedModListView, 0, 1);
        modGrid.add(variableModListView, 1, 1);
        modGrid.add(removeFixedModButton, 0, 2);
        modGrid.add(removeVariableModButton, 1, 2);
        modGrid.add(modSelectionContainer, 0, 3, 2, 1); // Span 2 columns
        
        section.getChildren().addAll(titleLabel, modGrid);
        return section;
    }
    
    
    // Inner class for metadata table rows
    public static class MetadataRow {
        private final SimpleStringProperty sampleName;
        private final SimpleStringProperty labelType;
        private final SimpleStringProperty age;
        private final SimpleStringProperty gender;
        private final SimpleStringProperty condition;
        private final SimpleStringProperty treatment;
        private final SimpleStringProperty timePoint;
        private final SimpleStringProperty customField1;
        private final SimpleStringProperty customField2;
        private final SimpleStringProperty customField3;
        
        public MetadataRow(String sampleName) {
            this.sampleName = new SimpleStringProperty(sampleName);
            this.labelType = new SimpleStringProperty("");
            this.age = new SimpleStringProperty("");
            this.gender = new SimpleStringProperty("");
            this.condition = new SimpleStringProperty("");
            this.treatment = new SimpleStringProperty("");
            this.timePoint = new SimpleStringProperty("");
            this.customField1 = new SimpleStringProperty("");
            this.customField2 = new SimpleStringProperty("");
            this.customField3 = new SimpleStringProperty("");
        }
        
        // Getters and setters
        public String getSampleName() { return sampleName.get(); }
        public void setSampleName(String value) { sampleName.set(value); }
        public StringProperty sampleNameProperty() { return sampleName; }
        
        public String getLabelType() { return labelType.get(); }
        public void setLabelType(String value) { labelType.set(value); }
        public StringProperty labelTypeProperty() { return labelType; }
        
        public String getAge() { return age.get(); }
        public void setAge(String value) { age.set(value); }
        public StringProperty ageProperty() { return age; }
        
        public String getGender() { return gender.get(); }
        public void setGender(String value) { gender.set(value); }
        public StringProperty genderProperty() { return gender; }
        
        public String getCondition() { return condition.get(); }
        public void setCondition(String value) { condition.set(value); }
        public StringProperty conditionProperty() { return condition; }
        
        public String getTreatment() { return treatment.get(); }
        public void setTreatment(String value) { treatment.set(value); }
        public StringProperty treatmentProperty() { return treatment; }
        
        public String getTimePoint() { return timePoint.get(); }
        public void setTimePoint(String value) { timePoint.set(value); }
        public StringProperty timePointProperty() { return timePoint; }
        
        public String getCustomField1() { return customField1.get(); }
        public void setCustomField1(String value) { customField1.set(value); }
        public StringProperty customField1Property() { return customField1; }
        
        public String getCustomField2() { return customField2.get(); }
        public void setCustomField2(String value) { customField2.set(value); }
        public StringProperty customField2Property() { return customField2; }
        
        public String getCustomField3() { return customField3.get(); }
        public void setCustomField3(String value) { customField3.set(value); }
        public StringProperty customField3Property() { return customField3; }
    }
    
    private void loadEnzymesAndModifications() {
        try {
            // Load enzymes from XML using dom4j
            InputStream enzymesStream = getClass().getClassLoader().getResourceAsStream("enzymes.xml");
            if (enzymesStream != null) {
                org.dom4j.io.SAXReader reader = new org.dom4j.io.SAXReader();
                org.dom4j.Document doc = reader.read(enzymesStream);
                
                // Use element() method instead of XPath to avoid Jaxen dependency
                org.dom4j.Element root = doc.getRootElement();
                java.util.List<org.dom4j.Element> enzymeElements = root.elements("enzyme");
                for (org.dom4j.Element enzyme : enzymeElements) {
                    String title = enzyme.attributeValue("title");
                    if (title != null && !title.isEmpty()) {
                        enzymes.add(title);
                    }
                }
                enzymesStream.close();
                System.out.println("Loaded " + enzymes.size() + " enzymes from enzymes.xml");
            } else {
                showError("Configuration Error", "Could not find enzymes.xml in resources");
                System.err.println("ERROR: enzymes.xml not found in classpath");
            }
            
            // Load modifications from XML using dom4j
            InputStream modificationsStream = getClass().getClassLoader().getResourceAsStream("modifications.xml");
            if (modificationsStream != null) {
                org.dom4j.io.SAXReader reader = new org.dom4j.io.SAXReader();
                org.dom4j.Document doc = reader.read(modificationsStream);
                
                // Use element() method instead of XPath to avoid Jaxen dependency
                org.dom4j.Element root = doc.getRootElement();
                java.util.List<org.dom4j.Element> modElements = root.elements("modification");
                for (org.dom4j.Element mod : modElements) {
                    String title = mod.attributeValue("title");
                    if (title != null && !title.isEmpty()) {
                        modifications.add(title);
                    }
                    // Add label modifications
                    org.dom4j.Element typeElement = mod.element("type");
                    if (typeElement != null) {
                        String type = typeElement.getTextTrim();
                        if ("Label".equalsIgnoreCase(type)) {
                            if (title != null && !title.isEmpty() && !labelModifications.contains(title)) {
                                labelModifications.add(title);
                            }
                        } else if ("IsobaricLabel".equalsIgnoreCase(type)) {
                            if (title != null && !title.isEmpty() && !isobaricModifications.contains(title)) {
                                isobaricModifications.add(title);
                                // Extract base type (e.g., "TMTpro16plex" from "TMTpro16plex-Nter133C")
                                String baseType = extractIsobaricBaseType(title);
                                if (baseType != null && !isobaricBaseTypes.contains(baseType)) {
                                    isobaricBaseTypes.add(baseType);
                                }
                            }
                        }
                    }
                }
                modificationsStream.close();
                System.out.println("Loaded " + modifications.size() + " modifications from modifications.xml");
                System.out.println("Loaded " + labelModifications.size() + " label modifications");
                System.out.println("Loaded " + isobaricModifications.size() + " isobaric modifications");
                System.out.println("Loaded " + isobaricBaseTypes.size() + " isobaric base types");
            } else {
                showError("Configuration Error", "Could not find modifications.xml in resources");
                System.err.println("ERROR: modifications.xml not found in classpath");
            }
        } catch (Exception e) {
            showError("Configuration Error", "Failed to load enzymes and modifications: " + e.getMessage());
            e.printStackTrace();
        }
    }
    
    private String extractIsobaricBaseType(String fullName) {
        // Extract base type from names like "TMTpro16plex-Nter133C" -> "TMTpro16plex"
        if (fullName.contains("-")) {
            return fullName.substring(0, fullName.indexOf("-"));
        }
        return fullName;
    }
    
    private void addFixedModification() {
        String selected = modificationComboBox.getValue();
        if (selected != null && !fixedModifications.contains(selected)) {
            fixedModifications.add(selected);
            updateParameterView();
        }
    }
    
    private void addVariableModification() {
        String selected = modificationComboBox.getValue();
        if (selected != null && !variableModifications.contains(selected)) {
            variableModifications.add(selected);
            updateParameterView();
        }
    }
    
    private void removeFixedModification() {
        String selected = fixedModListView.getSelectionModel().getSelectedItem();
        if (selected != null) {
            fixedModifications.remove(selected);
            updateParameterView();
        }
    }
    
    private void removeVariableModification() {
        String selected = variableModListView.getSelectionModel().getSelectedItem();
        if (selected != null) {
            variableModifications.remove(selected);
            updateParameterView();
        }
    }
    
    private void updateParameterView() {
        // Update the console panel's parameter view
        if (consolePanel != null) {
            consolePanel.setParameterText(getParameterDisplay());
        }
    }
    
    private String getParameterDisplay() {
        StringBuilder sb = new StringBuilder();
        
        // Enzyme settings
        sb.append("=== MaxQuant Analysis Parameters ===\n\n");
        sb.append("=== ENZYME SETTINGS ===\n");
        sb.append("Enzyme: ").append(enzymeComboBox.getValue()).append("\n");
        sb.append("Missed Cleavages: ").append(missedCleavagesSpinner.getValue()).append("\n");
        sb.append("Digestion Mode: ").append(digestionModeComboBox.getValue()).append("\n\n");
        
        // Additional parameters
        sb.append("=== ADDITIONAL PARAMETERS ===\n");
        sb.append("Threads: ").append(threadsSpinner.getValue()).append("\n");
        sb.append("MS2 Scan Mode: ").append(ms2ScanModeComboBox.getValue()).append("\n");
        sb.append("Generate Sample Specific DB: ").append(generateSampleSpecificDBCheckBox.isSelected() ? "Yes" : "No").append("\n");
        sb.append("Spectra Clustering: ").append(spectraClusteringCheckBox.isSelected() ? "Yes" : "No").append("\n\n");
        
        // Modifications
        sb.append("=== MODIFICATIONS ===\n");
        sb.append("Fixed Modifications (").append(fixedModifications.size()).append("):\n");
        for (String mod : fixedModifications) {
            sb.append("  - ").append(mod).append("\n");
        }
        
        sb.append("\nVariable Modifications (").append(variableModifications.size()).append("):\n");
        for (String mod : variableModifications) {
            sb.append("  - ").append(mod).append("\n");
        }
        
        // Quantitative mode
        sb.append("\n=== QUANTITATIVE MODE ===\n");
        if (quantModeGroup.getSelectedToggle() != null) {
            RadioButton selected = (RadioButton) quantModeGroup.getSelectedToggle();
            sb.append("Mode: ").append(selected.getText()).append("\n");
            
            if (selected.getText().contains("Isotopic")) {
                sb.append("Light Labels (").append(lightLabels.size()).append("): ");
                sb.append(lightLabels.isEmpty() ? "None" : String.join(", ", lightLabels)).append("\n");
                sb.append("Medium Labels (").append(mediumLabels.size()).append("): ");
                sb.append(mediumLabels.isEmpty() ? "None" : String.join(", ", mediumLabels)).append("\n");
                sb.append("Heavy Labels (").append(heavyLabels.size()).append("): ");
                sb.append(heavyLabels.isEmpty() ? "None" : String.join(", ", heavyLabels)).append("\n");
            } else if (selected.getText().contains("Isobaric")) {
                sb.append("Isobaric Type: ").append(isobaricTypeComboBox.getValue() != null ? isobaricTypeComboBox.getValue() : "None").append("\n");
                sb.append("Configured Labels: ").append(isobaricTableData.size()).append("\n");
            }
        }
        
        // Taxonomic analysis
        
        return sb.toString();
    }
    
    private boolean validateInputs() {
        if (ioPanel.getIoPanel().getExperimentData().isEmpty()) {
            showError("Validation Error", "Please add at least one experiment.");
            return false;
        }
        if (ioPanel.getIoPanel().getOutputDirectory().isEmpty()) {
            showError("Validation Error", "Please select an output directory.");
            return false;
        }
        if (ioPanel.getMicrobiomeDBPath().isEmpty()) {
            showError("Validation Error", "Please select a microbiome database file.");
            return false;
        }
        if (ioPanel.isHostDBEnabled() && ioPanel.getHostDBPath().isEmpty()) {
            showError("Validation Error", "Please select a host database file.");
            return false;
        }
        if (enzymeComboBox.getValue() == null) {
            showError("Validation Error", "Please select an enzyme.");
            return false;
        }
        return true;
    }
    
    private void startAnalysis() {
        try {
            // Validate inputs first
            if (!validateInputs()) {
                return;
            }
            
            // Collect ALL parameters from UI into a single object
            DDAIgcMaxQuantParameters analysisParams = collectParametersFromUI();
            
            // Save parameters for this analysis and as last used
            String outputDir = analysisParams.getResult();
            if (outputDir != null && !outputDir.trim().isEmpty()) {
                // Save analysis-specific parameters
                boolean analysisParamsSaved = parameterIO.saveAnalysisParameters(outputDir, analysisParams);
                if (!analysisParamsSaved) {
                    LOGGER.warn("Failed to save analysis parameters to {}/analysis_parameters.json", outputDir);
                    showAlert("Warning", "Failed to save analysis parameters to result directory. Analysis will continue but parameters may not be preserved.");
                } else {
                    LOGGER.info("Saved analysis parameters to {}/analysis_parameters.json", outputDir);
                }
                
                // Save as last used parameters
                boolean lastParamsSaved = parameterIO.saveLastUsedParameters(analysisParams);
                if (!lastParamsSaved) {
                    LOGGER.warn("Failed to save last used parameters");
                }
                
                // Also save to preferences for quick defaults
                parameterIO.savePreferences(analysisParams);
            }
            
            // Create MetaResourceParameters by loading from config file
            uk.ac.quadram.metapilot.resources.MetaResourceIO resourceIO = new uk.ac.quadram.metapilot.resources.MetaResourceIO();
            MetaResourceParameters mrp = resourceIO.loadParameters();

            // Create and run DDAIgcMaxQuantTask with the separate progress bars
			currentTask = new DDAIgcMaxQuantTask(analysisParams, mrp, consolePanel.getOverallProgressBar(),
					consolePanel.getCurrentTaskProgressBar());

            // Bind task message to console panel
            currentTask.messageProperty().addListener((obs, oldMsg, newMsg) -> {
                Platform.runLater(() -> {
                    if (newMsg != null && !newMsg.isEmpty()) {
                        consolePanel.updateCurrentTaskFromMessage(newMsg);
                    }
                });
            });
            
            // Handle task completion
            currentTask.setOnSucceeded(e -> {
                Platform.runLater(() -> {
                    // Unbind progress bars to prevent memory leaks
                    consolePanel.getOverallProgressBar().progressProperty().unbind();
                    consolePanel.getCurrentTaskProgressBar().progressProperty().unbind();
                    
                    consolePanel.setButtonsEnabled(true, false, true);
                    consolePanel.stopStatusTimer();
                    
                    // Check if the task actually succeeded by examining the result
                    try {
                        Boolean result = currentTask.get();
                        if (result != null && result) {
                            showAlert("Success", "Analysis completed successfully.");
                        } else {
                            showAlert("Warning", "Analysis completed but may not have been fully successful.");
                        }
                    } catch (Exception ex) {
                        showAlert("Error", "Analysis failed: " + ex.getMessage());
                    }
                });
            });
            currentTask.setOnFailed(e -> {
                Platform.runLater(() -> {
                    // Unbind progress bars to prevent memory leaks
                    consolePanel.getOverallProgressBar().progressProperty().unbind();
                    consolePanel.getCurrentTaskProgressBar().progressProperty().unbind();
                    
                    consolePanel.setButtonsEnabled(true, false, false);
                    consolePanel.stopStatusTimer();
                    showAlert("Error", "Analysis failed: " + currentTask.getException());
                });
            });
            
            // Start the task in a new thread
            Thread t = new Thread(currentTask);
            t.setDaemon(true);
            t.start();
            
            // Update UI buttons and start status timer
            consolePanel.setButtonsEnabled(false, true, false);
            consolePanel.startStatusTimer();
        } catch (Exception e) {
            showAlert("Error", "Failed to start analysis: " + e.getMessage());
            e.printStackTrace();
        }
    }
    
    private void stopAnalysis() {
        Alert alert = new Alert(AlertType.CONFIRMATION);
        alert.setTitle("Stop Analysis");
        alert.setHeaderText("Are you sure you want to stop the analysis?");
        alert.setContentText("This action cannot be undone.");
        
        Optional<ButtonType> result = alert.showAndWait();
        if (result.isPresent() && result.get() == ButtonType.OK) {
            // Force stop the running task if it exists
            if (currentTask != null && currentTask.isRunning()) {
                consolePanel.appendConsoleText("Stopping MaxQuant analysis...\n");
                currentTask.forceStop();
                consolePanel.appendConsoleText("MaxQuant analysis stop requested.\n");
            } else {
                consolePanel.appendConsoleText("No running MaxQuant analysis to stop.\n");
            }
            
            // Unbind progress bars to prevent memory leaks
            consolePanel.getOverallProgressBar().progressProperty().unbind();
            consolePanel.getCurrentTaskProgressBar().progressProperty().unbind();
            
            consolePanel.setButtonsEnabled(true, false, false);
            consolePanel.stopStatusTimer();
            consolePanel.getOverallProgressBar().setProgress(0);
            consolePanel.getCurrentTaskProgressBar().setProgress(0);
            consolePanel.setCurrentTaskLabel("Current Task: None");
        }
    }
    
    private void viewReport() {
        try {
            // Get the result directory from the IO panel
            String resultDir = ioPanel.getIoPanel().getOutputDirectory();
            if (resultDir == null || resultDir.trim().isEmpty()) {
                showAlert("Error", "No result directory available. Please run analysis first.");
                return;
            }
            
            // Look for the main report file
            java.io.File reportDir = new java.io.File(resultDir);
            java.io.File indexFile = new java.io.File(reportDir, "index.html");
            
            if (!indexFile.exists()) {
                showAlert("Error", "Report not found. Please ensure the analysis has completed successfully.");
                return;
            }
            
            // Open the report in the default browser
            if (java.awt.Desktop.isDesktopSupported()) {
                java.awt.Desktop.getDesktop().browse(indexFile.toURI());
                showAlert("Success", "Report opened in default browser.");
            } else {
                showAlert("Error", "Desktop not supported. Please open manually: " + indexFile.getAbsolutePath());
            }
            
        } catch (Exception e) {
            showAlert("Error", "Failed to open report: " + e.getMessage());
            e.printStackTrace();
        }
    }
    
    private void showError(String title, String content) {
        Platform.runLater(() -> {
            Alert alert = new Alert(AlertType.ERROR);
            alert.setTitle(title);
            alert.setHeaderText(null);
            alert.setContentText(content);
            alert.showAndWait();
        });
    }
    
    private Label createSectionTitle(String title) {
        Label label = new Label(title);
        label.setFont(Font.font("System", FontWeight.BOLD, 14));
        label.setTextFill(Color.DARKBLUE);
        return label;
    }
    
    private void updateIsobaricTable() {
        String selectedType = isobaricTypeComboBox.getValue();
        if (selectedType != null) {
            isobaricTableData.clear();
            List<String> labels = getIsobaricLabelsForType(selectedType);
            for (String label : labels) {
                isobaricTableData.add(new IsobaricLabelRow(label));
            }
            updateMetadataTable();
        }
    }
    
    private List<String> getIsobaricLabelsForType(String baseType) {
        List<String> labels = new ArrayList<>();
        
        switch (baseType) {
            case "iTRAQ4plex":
                labels.add("iTRAQ4plex-114");
                labels.add("iTRAQ4plex-115");
                labels.add("iTRAQ4plex-116");
                labels.add("iTRAQ4plex-117");
                break;
                
            case "iTRAQ8plex":
                labels.add("iTRAQ8plex-113");
                labels.add("iTRAQ8plex-114");
                labels.add("iTRAQ8plex-115");
                labels.add("iTRAQ8plex-116");
                labels.add("iTRAQ8plex-117");
                labels.add("iTRAQ8plex-118");
                labels.add("iTRAQ8plex-119");
                labels.add("iTRAQ8plex-121");
                break;
                
            case "TMT6plex":
                labels.add("TMT6plex-126");
                labels.add("TMT6plex-127");
                labels.add("TMT6plex-128");
                labels.add("TMT6plex-129");
                labels.add("TMT6plex-130");
                labels.add("TMT6plex-131");
                break;
                
            case "TMT10plex":
                labels.add("TMT10plex-126");
                labels.add("TMT10plex-127N");
                labels.add("TMT10plex-127C");
                labels.add("TMT10plex-128N");
                labels.add("TMT10plex-128C");
                labels.add("TMT10plex-129N");
                labels.add("TMT10plex-129C");
                labels.add("TMT10plex-130N");
                labels.add("TMT10plex-130C");
                labels.add("TMT10plex-131");
                break;
                
            case "TMT11plex":
                // Include all TMT10plex labels
                labels.addAll(getIsobaricLabelsForType("TMT10plex"));
                // Add TMT11plex specific labels
                labels.add("TMT11plex-131C");
                break;
                
            case "TMTpro16plex":
                labels.add("TMTpro16plex-126");
                labels.add("TMTpro16plex-127N");
                labels.add("TMTpro16plex-127C");
                labels.add("TMTpro16plex-128N");
                labels.add("TMTpro16plex-128C");
                labels.add("TMTpro16plex-129N");
                labels.add("TMTpro16plex-129C");
                labels.add("TMTpro16plex-130N");
                labels.add("TMTpro16plex-130C");
                labels.add("TMTpro16plex-131N");
                labels.add("TMTpro16plex-131C");
                labels.add("TMTpro16plex-132N");
                labels.add("TMTpro16plex-132C");
                labels.add("TMTpro16plex-133N");
                labels.add("TMTpro16plex-133C");
                labels.add("TMTpro16plex-134N");
                break;
                
            case "TMTpro18plex":
                // Include all TMTpro16plex labels
                labels.addAll(getIsobaricLabelsForType("TMTpro16plex"));
                // Add TMTpro18plex specific labels
                labels.add("TMTpro18plex-134C");
                labels.add("TMTpro18plex-135N");
                break;
                
            default:
                // Fallback: try to find modifications that start with the selected base type
                for (String modification : isobaricModifications) {
                    if (modification.startsWith(baseType + "-")) {
                        labels.add(modification);
                    }
                }
                break;
        }
        
        return labels;
    }
    
    // Methods for isotopic labeling
    private void addLightLabel() {
        String selected = labelModificationComboBox.getValue();
        if (selected != null && !lightLabels.contains(selected)) {
            lightLabels.add(selected);
            updateParameterView();
            updateMetadataTable();
        }
    }
    
    private void addMediumLabel() {
        String selected = labelModificationComboBox.getValue();
        if (selected != null && !mediumLabels.contains(selected)) {
            mediumLabels.add(selected);
            updateParameterView();
            updateMetadataTable();
        }
    }
    
    private void addHeavyLabel() {
        String selected = labelModificationComboBox.getValue();
        if (selected != null && !heavyLabels.contains(selected)) {
            heavyLabels.add(selected);
            updateParameterView();
            updateMetadataTable();
        }
    }
    
    private void removeLightLabel() {
        String selected = lightLabelListView.getSelectionModel().getSelectedItem();
        if (selected != null) {
            lightLabels.remove(selected);
            updateParameterView();
            updateMetadataTable();
        }
    }
    
    private void removeMediumLabel() {
        String selected = mediumLabelListView.getSelectionModel().getSelectedItem();
        if (selected != null) {
            mediumLabels.remove(selected);
            updateParameterView();
            updateMetadataTable();
        }
    }
    
    private void removeHeavyLabel() {
        String selected = heavyLabelListView.getSelectionModel().getSelectedItem();
        if (selected != null) {
            heavyLabels.remove(selected);
            updateParameterView();
            updateMetadataTable();
        }
    }
    
    // Inner class for isobaric labeling table rows
    public static class IsobaricLabelRow {
        private final SimpleStringProperty labelName;
        private final SimpleStringProperty minus2;
        private final SimpleStringProperty minus1;
        private final SimpleStringProperty plus1;
        private final SimpleStringProperty plus2;
        
        public IsobaricLabelRow(String labelName) {
            this.labelName = new SimpleStringProperty(labelName);
            this.minus2 = new SimpleStringProperty("0.0");
            this.minus1 = new SimpleStringProperty("0.0");
            this.plus1 = new SimpleStringProperty("0.0");
            this.plus2 = new SimpleStringProperty("0.0");
        }
        
        public String getLabelName() { return labelName.get(); }
        public void setLabelName(String value) { labelName.set(value); }
        public StringProperty labelNameProperty() { return labelName; }
        
        public String getMinus2() { return minus2.get(); }
        public void setMinus2(String value) { minus2.set(value); }
        public StringProperty minus2Property() { return minus2; }
        
        public String getMinus1() { return minus1.get(); }
        public void setMinus1(String value) { minus1.set(value); }
        public StringProperty minus1Property() { return minus1; }
        
        public String getPlus1() { return plus1.get(); }
        public void setPlus1(String value) { plus1.set(value); }
        public StringProperty plus1Property() { return plus1; }
        
        public String getPlus2() { return plus2.get(); }
        public void setPlus2(String value) { plus2.set(value); }
        public StringProperty plus2Property() { return plus2; }
    }
    
    private void updateMetadataTableFromIOPanel() {
        // Update metadata table when Parameters tab is selected
        updateMetadataTable();
    }
    
    private void showAlert(String title, String message) {
        Alert alert = new Alert(AlertType.INFORMATION);
        alert.setTitle(title);
        alert.setHeaderText(null);
        alert.setContentText(message);
        alert.showAndWait();
    }
    
    
    /**
     * Collect all parameters from the UI into a parameters object
     */
    private DDAIgcMaxQuantParameters collectParametersFromUI() {
        // Collect experiment data from ioPanel
        ObservableList<ExperimentRow> fileRows = ioPanel.getIoPanel().getExperimentData();
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
        }
        MetaData metadata = new MetaData(rawFiles, expNames, fractions, replicates);
        
        // Collect basic parameters with null checks
        String resultDir = ioPanel.getIoPanel().getOutputDirectory();
        String microDb = ioPanel.getMicrobiomeDBPath();
        boolean appendHostDb = ioPanel.isHostDBEnabled();
        String hostDB = appendHostDb ? ioPanel.getHostDBPath() : "";
        String selectedEnzymeStr = enzymeComboBox.getValue();
        if (selectedEnzymeStr == null) {
            selectedEnzymeStr = "Trypsin"; // Default
        }
        int missCleavages = missedCleavagesSpinner.getValue();
        String digestMode = digestionModeComboBox.getValue();
        if (digestMode == null) {
            digestMode = "Specific"; // Default digestion mode
        }
        String ms2ScanMode = ms2ScanModeComboBox.getValue();
        if (ms2ScanMode == null) {
            ms2ScanMode = MetaConstants.HCD_FTMS; // Default
        }
        int threadCount = threadsSpinner.getValue();
        boolean spectraCluster = spectraClusteringCheckBox.isSelected();
        
        // Create DDAIgcMaxQuantParameters using parameterized constructor
		DDAIgcMaxQuantParameters params = new DDAIgcMaxQuantParameters(
				MetaPilotApp.WorkflowType.DDA_IGC_MAXQUANT.getAnalysisType(), // analysisType
				metadata, // MetaData object
				resultDir, // output directory (result)
				microDb, // microbiome database (microDb)
				appendHostDb, // append host DB
				hostDB, // host database path (may be empty)
				ms2ScanMode, // ms2ScanMode
				threadCount // thread count
		);
        
        // Set additional parameters not covered by constructor
        params.setEnzyme(selectedEnzymeStr); // Set enzyme
		params.setMissCleavages(missCleavages);
        params.setDigestMode(digestMode); // Set digestion mode
        params.setSpectraCluster(spectraCluster);

        // Additional parameters
        params.setReduceDb(generateSampleSpecificDBCheckBox.isSelected());
        
        // Modifications
        params.setFixedModifications(new ArrayList<>(fixedModifications));
        params.setVariableModifications(new ArrayList<>(variableModifications));
        
        // Quantitative mode
        if (quantModeGroup.getSelectedToggle() != null) {
            RadioButton selected = (RadioButton) quantModeGroup.getSelectedToggle();
            params.setQuantitativeMode(selected.getText());
        }
        
        // Isotopic labeling
        params.setLightLabels(new ArrayList<>(lightLabels));
        params.setMediumLabels(new ArrayList<>(mediumLabels));
        params.setHeavyLabels(new ArrayList<>(heavyLabels));
        
        // Set isotopic labels (String[][])
        String[][] labels = new String[3][];
        labels[0] = lightLabels.toArray(new String[0]);
        labels[1] = mediumLabels.toArray(new String[0]);
        labels[2] = heavyLabels.toArray(new String[0]);
        params.setIsotopicLabels(labels);
        
        // Isobaric labeling
        if (isobaricTypeComboBox.getValue() != null) {
            params.setIsobaricType(isobaricTypeComboBox.getValue());
            
            // Set isobaric tag (IsobaricTag) from ComboBox using getName()
            String isobaricType = isobaricTypeComboBox.getValue();
            if (isobaricType != null && !isobaricType.isEmpty()) {
                for (IsobaricTag tag : IsobaricTag.values()) {
                    if (tag.getName().equals(isobaricType)) {
                        params.setIsobaricTag(tag);
                        break;
                    }
                }
            } else {
                params.setIsobaricTag(null);
            }
        }
        params.setIsobaricTableData(new ArrayList<>(isobaricTableData));
        
        // Set isobaric correction factors (double[][])
        if (!isobaricTableData.isEmpty()) {
            double[][] isoCorFactor = new double[isobaricTableData.size()][4];
            for (int i = 0; i < isobaricTableData.size(); i++) {
                IsobaricLabelRow row = isobaricTableData.get(i);
                try {
                    isoCorFactor[i][0] = Double.parseDouble(row.getMinus2());
                    isoCorFactor[i][1] = Double.parseDouble(row.getMinus1());
                    isoCorFactor[i][2] = Double.parseDouble(row.getPlus1());
                    isoCorFactor[i][3] = Double.parseDouble(row.getPlus2());
                } catch (NumberFormatException e) {
                    isoCorFactor[i][0] = 0.0;
                    isoCorFactor[i][1] = 0.0;
                    isoCorFactor[i][2] = 0.0;
                    isoCorFactor[i][3] = 0.0;
                }
            }
            params.setIsobaricRefTable(isoCorFactor);
        }
        
        // Taxonomic analysis
        
        // Metadata
        params.setMetadataTableData(new ArrayList<>(metadataTableData));
        
        return params;
    }
    
    /**
     * Load last used parameters and apply them to the UI
     */
    private void loadLastUsedParameters() {
        try {
            DDAIgcMaxQuantParameters params = (DDAIgcMaxQuantParameters) parameterIO.loadLastUsedParameters();
            applyParametersToUI(params);
            LOGGER.info("Loaded last used parameters for DDA IGC MaxQuant workflow");
        } catch (Exception e) {
            LOGGER.warn("Failed to load last used parameters: " + e.getMessage());
            // Continue with default UI state
        }
    }
    
    /**
     * Apply parameters to the UI components
     */
    private void applyParametersToUI(DDAIgcMaxQuantParameters params) {
        if (params == null) return;
        
        // Apply Input/Output settings
        if (params.getResult() != null && ioPanel != null) {
            ioPanel.getIoPanel().setOutputDirectory(params.getResult());
        }
        
        // Set the result directory for logging (this is used by MetaPilotTask.configLog)
        if (params.getResult() != null && !params.getResult().trim().isEmpty()) {
            params.setResult(params.getResult());
        }
        
        // Apply database settings
        if (params.getMicroDb() != null && ioPanel != null) {
            ioPanel.setMicrobiomeDBPath(params.getMicroDb());
        }
        if (ioPanel != null) {
            ioPanel.setHostDBEnabled(params.isAppendHostDb());
        }
        if (params.getHostDB() != null && ioPanel != null) {
            ioPanel.setHostDBPath(params.getHostDB());
        }
        
        // Apply experiment data to IO panel
        if (params.getMetadata() != null && ioPanel != null) {
            MetaData metadata = params.getMetadata();
            String[] rawFiles = metadata.getRawFiles();
            String[] expNames = metadata.getExpNames();
            
            if (rawFiles != null && rawFiles.length > 0) {
                List<MetaPilotIOPanelFX.ExperimentRow> experimentData = new ArrayList<>();
                for (int i = 0; i < rawFiles.length; i++) {
                    String rawFile = rawFiles[i];
                    String expName = (expNames != null && i < expNames.length) ? expNames[i] : "";
                    
                    // Only add rows that have meaningful data (not blank rows)
                    if (rawFile != null && !rawFile.trim().isEmpty()) {
                        experimentData.add(new MetaPilotIOPanelFX.ExperimentRow(rawFile, expName));
                    }
                }
                ioPanel.getIoPanel().setExperimentData(experimentData);
            }
        }
        
        // Apply basic settings that don't depend on file paths
        if (params.getEnzyme() != null && enzymeComboBox != null) {
            enzymeComboBox.setValue(params.getEnzyme());
        }
        if (missedCleavagesSpinner != null) {
            missedCleavagesSpinner.getValueFactory().setValue(params.getMissCleavages());
        }
        if (params.getDigestModeString() != null && digestionModeComboBox != null) {
            digestionModeComboBox.setValue(params.getDigestModeString());
        }
        if (threadsSpinner != null) {
            threadsSpinner.getValueFactory().setValue(params.getThreadCount());
        }
        if (params.getMs2ScanMode() != null && ms2ScanModeComboBox != null) {
            ms2ScanModeComboBox.setValue(params.getMs2ScanMode());
        }
        if (generateSampleSpecificDBCheckBox != null) {
            generateSampleSpecificDBCheckBox.setSelected(params.isReduceDb());
        }
        if (spectraClusteringCheckBox != null) {
            spectraClusteringCheckBox.setSelected(params.isSpectraCluster());
        }
        
        // Apply modifications (clear existing and add from params)
        if (fixedModifications != null && params.getFixedModifications() != null) {
            fixedModifications.clear();
            fixedModifications.addAll(params.getFixedModifications());
        }
        if (variableModifications != null && params.getVariableModifications() != null) {
            variableModifications.clear();
            variableModifications.addAll(params.getVariableModifications());
        }
        
        // Apply quantitative mode settings
        if (params.getQuantitativeMode() != null && quantModeGroup != null) {
            String mode = params.getQuantitativeMode();
            for (Toggle toggle : quantModeGroup.getToggles()) {
                RadioButton radioButton = (RadioButton) toggle;
                if (radioButton.getText().equals(mode)) {
                    radioButton.setSelected(true);
                    break;
                }
            }
        }
        
        // Apply isotopic labeling settings
        if (params.getLightLabels() != null && lightLabels != null) {
            lightLabels.clear();
            lightLabels.addAll(params.getLightLabels());
        }
        if (params.getMediumLabels() != null && mediumLabels != null) {
            mediumLabels.clear();
            mediumLabels.addAll(params.getMediumLabels());
        }
        if (params.getHeavyLabels() != null && heavyLabels != null) {
            heavyLabels.clear();
            heavyLabels.addAll(params.getHeavyLabels());
        }
        
        // Apply isobaric labeling settings
        if (params.getIsobaricType() != null && isobaricTypeComboBox != null) {
            isobaricTypeComboBox.setValue(params.getIsobaricType());
        }
        if (params.getIsobaricTableData() != null && isobaricTableData != null) {
            isobaricTableData.clear();
            isobaricTableData.addAll(params.getIsobaricTableData());
        }
        
        // Apply taxonomic settings
        
        // Apply metadata table data
        if (params.getMetadataTableData() != null && metadataTableData != null) {
            metadataTableData.clear();
            metadataTableData.addAll(params.getMetadataTableData());
        }
        
        // Update parameter view and metadata table after applying parameters
        Platform.runLater(() -> {
            updateParameterView();
            updateMetadataTable();
        });
    }
    
    /**
     * Load parameters from a specific analysis result directory
     */
    public void loadAnalysisParameters(String resultDirectory) {
        try {
            DDAIgcMaxQuantParameters params = (DDAIgcMaxQuantParameters) parameterIO.loadAnalysisParameters(resultDirectory);
            applyParametersToUI(params);
            showAlert("Success", "Loaded analysis parameters from " + resultDirectory);
            LOGGER.info("Loaded analysis parameters from {}", resultDirectory);
        } catch (Exception e) {
            showAlert("Error", "Failed to load analysis parameters: " + e.getMessage());
            LOGGER.error("Failed to load analysis parameters from " + resultDirectory, e);
        }
    }
    
    /**
     * Set up window close handler to save parameters when the window is closed
     */
    private void setupWindowCloseHandler() {
        Platform.runLater(() -> {
            try {
                // Get the current scene and window
                if (getScene() != null) {
                    javafx.stage.Window window = getScene().getWindow();
                    if (window != null) {
                        window.setOnCloseRequest(event -> {
                            try {
                                DDAIgcMaxQuantParameters currentParams = collectParametersFromUI();
                                if (currentParams != null) {
                                    parameterIO.saveLastUsedParameters(currentParams);
                                    parameterIO.savePreferences(currentParams);
                                    LOGGER.info("Saved current parameters when window was closed");
                                }
                            } catch (Exception e) {
                                LOGGER.warn("Failed to save parameters when window was closed: " + e.getMessage(), e);
                            }
                        });
                        LOGGER.info("Window close handler set up successfully");
                    } else {
                        LOGGER.warn("Window is null, cannot set close handler");
                    }
                } else {
                    LOGGER.warn("Scene is null, cannot set close handler");
                }
            } catch (Exception e) {
                LOGGER.error("Failed to set up window close handler", e);
            }
        });
    }
} 