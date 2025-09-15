package uk.ac.quadram.metapilot.task.dda.igc.gui;

import javafx.application.Platform;
import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.control.*;
import javafx.scene.control.Toggle;
import javafx.scene.layout.*;
import javafx.stage.FileChooser;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.scene.text.Font;
import javafx.scene.text.FontWeight;
import javafx.scene.paint.Color;
import javafx.scene.control.cell.PropertyValueFactory;
import javafx.scene.control.cell.TextFieldTableCell;
import javafx.util.converter.DefaultStringConverter;
import javafx.scene.control.Alert.AlertType;
import uk.ac.quadram.metapilot.task.common.gui.MetaPilotIOPanelFX;
import uk.ac.quadram.metapilot.task.common.gui.IsobaricLabelingPanel;
import uk.ac.quadram.metapilot.task.common.gui.TaxonomicAnalysisPanel;
import uk.ac.quadram.metapilot.task.common.gui.DdaIgcInputOutputPanel;
import uk.ac.quadram.metapilot.task.par.MetaData;
import uk.ac.quadram.metapilot.MetaPilotApp;
import uk.ac.quadram.metapilot.core.mod.IsobaricTag;
import uk.ac.quadram.metapilot.task.dda.igc.DDAIgcPFindTask;
import uk.ac.quadram.metapilot.task.dda.igc.par.DDAIgcPFindParameters;
import uk.ac.quadram.metapilot.task.dda.igc.par.DDAIgcPFindParameterIO;
import javafx.animation.KeyFrame;
import javafx.animation.Timeline;
import javafx.util.Duration;
import java.lang.management.ManagementFactory;
import java.lang.management.MemoryMXBean;
import java.lang.management.MemoryUsage;
import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.io.File;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class MetaPilotDdaIgcPFindPanel extends VBox {
    private static final Logger LOGGER = LogManager.getLogger(MetaPilotDdaIgcPFindPanel.class);
    // UI Components (similar to MaxQuant panel)
    private DdaIgcInputOutputPanel ioPanel;
    private TextArea consoleOutput;
    private TextArea parameterView;
    private ProgressBar overallProgressBar;
    private ProgressBar currentTaskProgressBar;
    private Label currentTaskLabel;
    private Button startButton;
    private Button stopButton;
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
    private ObservableList<String> enzymes = FXCollections.observableArrayList();
    private ObservableList<String> modifications = FXCollections.observableArrayList();
    private ObservableList<String> fixedModifications = FXCollections.observableArrayList();
    private ObservableList<String> variableModifications = FXCollections.observableArrayList();
    // Labeling
    private ObservableList<String> isobaricModifications = FXCollections.observableArrayList();
    private ObservableList<String> isobaricBaseTypes = FXCollections.observableArrayList();

    // Isobaric labeling components
    private ComboBox<String> isobaricTypeComboBox;
    private TableView<IsobaricLabelingPanel.IsobaricLabelRow> isobaricTable;
    private ObservableList<IsobaricLabelingPanel.IsobaricLabelRow> isobaricTableData = FXCollections.observableArrayList();

    // Taxonomic analysis panel
    private TaxonomicAnalysisPanel taxonomicAnalysisPanel;

    // Metadata configuration components
    private TableView<MetaPilotDdaIgcMaxQuantPanel.MetadataRow> metadataTable;
    private ObservableList<MetaPilotDdaIgcMaxQuantPanel.MetadataRow> metadataTableData = FXCollections.observableArrayList();
    private ComboBox<String> metadataColumnComboBox;
    private Button addMetadataColumnButton;
    private Button removeMetadataColumnButton;

    // Status panel fields
    private Label statusTimeLabel;
    private Label statusCpuLabel;
    private Label statusMemLabel;
    private Timeline statusTimeline;
    private long taskStartTimeMillis;
    
    // Add parameter IO instance
    private DDAIgcPFindParameterIO parameterIO;
    
    public MetaPilotDdaIgcPFindPanel() {
        // Initialize parameter IO instance
        parameterIO = new DDAIgcPFindParameterIO();
        
        loadEnzymesFromIni();
        loadModificationsFromIni();
        loadIsobaricBaseTypes();
        
        // Initialize unified IO panel
        ioPanel = new DdaIgcInputOutputPanel(this::updateParameterView);
        
        // Initialize taxonomic analysis panel
        taxonomicAnalysisPanel = new TaxonomicAnalysisPanel(this::updateParameterView);
        
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

        Tab consoleTab = new Tab("Console");
        consoleTab.setContent(createConsoleTab());

        tabPane.getTabs().addAll(ioTab, paramTab, consoleTab);

        // Add tab change listener to update metadata when Parameters tab is selected
        // and save parameters when Console tab is selected
        tabPane.getSelectionModel().selectedItemProperty().addListener((obs, oldTab, newTab) -> {
            if (newTab == paramTab) {
                updateMetadataTableFromIOPanel();
            } else if (newTab == consoleTab) {
                // Save current parameters when Console tab is selected
                try {
                    DDAIgcPFindParameters currentParams = collectParametersFromUI();
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
        getChildren().add(tabPane);
    }

    // --- BEGIN: Methods adapted from MetaPilotDdaIgcMaxQuantPanel ---
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
        // Data already loaded in constructor
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
        leftColumn.getChildren().addAll(enzymeSection, additionalParamsSection, taxonomicAnalysisPanel, modificationsSection);
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
        RadioButton isobaricRadio = new RadioButton("Isobaric");
        isobaricRadio.setToggleGroup(quantModeGroup);
        HBox radioBox = new HBox(15);
        radioBox.setAlignment(Pos.CENTER_LEFT);
        radioBox.getChildren().addAll(labelFreeRadio, isobaricRadio);
        // StackPane for overlapping panels
        StackPane labelingStackPane = new StackPane();
        labelingStackPane.setMinHeight(200);
        labelingStackPane.setMaxHeight(Double.MAX_VALUE);
        VBox.setVgrow(labelingStackPane, Priority.ALWAYS);
        VBox isobaricPanel = createCompactIsobaricLabelingPanel();
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
    private VBox createCompactIsobaricLabelingPanel() {
        VBox panel = new VBox(8);
        panel.setPadding(new Insets(8));
        panel.setStyle("-fx-background-color: #f8e8f4; -fx-background-radius: 5;");
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
        isobaricTable = new TableView<>(isobaricTableData);
        isobaricTable.setMinHeight(120);
        isobaricTable.setPrefHeight(140);
        isobaricTable.setMaxHeight(Double.MAX_VALUE);
        isobaricTable.setEditable(true);
        isobaricTable.setColumnResizePolicy(TableView.CONSTRAINED_RESIZE_POLICY);
        VBox.setVgrow(isobaricTable, Priority.ALWAYS);
        TableColumn<IsobaricLabelingPanel.IsobaricLabelRow, String> labelCol = new TableColumn<>("Label");
        labelCol.setCellValueFactory(new PropertyValueFactory<>("labelName"));
        labelCol.setMinWidth(120);
        labelCol.setPrefWidth(150);
        TableColumn<IsobaricLabelingPanel.IsobaricLabelRow, String> minus2Col = new TableColumn<>("-2");
        minus2Col.setCellValueFactory(new PropertyValueFactory<>("minus2"));
        minus2Col.setCellFactory(TextFieldTableCell.forTableColumn(new DefaultStringConverter()));
        minus2Col.setOnEditCommit(e -> e.getRowValue().setMinus2(e.getNewValue()));
        minus2Col.setMinWidth(50);
        minus2Col.setPrefWidth(60);
        TableColumn<IsobaricLabelingPanel.IsobaricLabelRow, String> minus1Col = new TableColumn<>("-1");
        minus1Col.setCellValueFactory(new PropertyValueFactory<>("minus1"));
        minus1Col.setCellFactory(TextFieldTableCell.forTableColumn(new DefaultStringConverter()));
        minus1Col.setOnEditCommit(e -> e.getRowValue().setMinus1(e.getNewValue()));
        minus1Col.setMinWidth(50);
        minus1Col.setPrefWidth(60);
        TableColumn<IsobaricLabelingPanel.IsobaricLabelRow, String> plus1Col = new TableColumn<>("+1");
        plus1Col.setCellValueFactory(new PropertyValueFactory<>("plus1"));
        plus1Col.setCellFactory(TextFieldTableCell.forTableColumn(new DefaultStringConverter()));
        plus1Col.setOnEditCommit(e -> e.getRowValue().setPlus1(e.getNewValue()));
        plus1Col.setMinWidth(50);
        plus1Col.setPrefWidth(60);
        TableColumn<IsobaricLabelingPanel.IsobaricLabelRow, String> plus2Col = new TableColumn<>("+2");
        plus2Col.setCellValueFactory(new PropertyValueFactory<>("plus2"));
        plus2Col.setCellFactory(TextFieldTableCell.forTableColumn(new DefaultStringConverter()));
        plus2Col.setOnEditCommit(e -> e.getRowValue().setPlus2(e.getNewValue()));
        plus2Col.setMinWidth(50);
        plus2Col.setPrefWidth(60);
        isobaricTable.getColumns().addAll(labelCol, minus2Col, minus1Col, plus1Col, plus2Col);
        panel.getChildren().addAll(typeBox, isobaricTable);
        return panel;
    }
    private VBox createConsoleTab() {
        VBox tabContent = new VBox(10);
        tabContent.setPadding(new Insets(10));
        // Parameter View
        VBox paramViewBox = new VBox(5);
        paramViewBox.getChildren().add(new Label("Parameter View:"));
        parameterView = new TextArea();
        parameterView.setEditable(false);
        parameterView.setPrefHeight(200); // Increased height
        parameterView.setStyle("-fx-font-family: monospace;");
        VBox.setVgrow(parameterView, Priority.ALWAYS); // Allow to expand
        paramViewBox.getChildren().add(parameterView);
        VBox.setVgrow(paramViewBox, Priority.ALWAYS);
        // Progress Bars
        VBox progressBox = new VBox(5);
        currentTaskLabel = new Label("Current Task: None");
        overallProgressBar = new ProgressBar(0);
        overallProgressBar.setMaxWidth(Double.MAX_VALUE);
        currentTaskProgressBar = new ProgressBar(0);
        currentTaskProgressBar.setMaxWidth(Double.MAX_VALUE);
        progressBox.getChildren().addAll(
            currentTaskLabel,
            new Label("Overall Progress:"),
            overallProgressBar,
            new Label("Current Task Progress:"),
            currentTaskProgressBar
        );
        // Console Output
        VBox consoleBox = new VBox(5);
        consoleBox.getChildren().add(new Label("Console Output:"));
        consoleOutput = new TextArea();
        consoleOutput.setEditable(false);
        consoleOutput.setPrefHeight(250); // Increased height
        consoleOutput.setStyle("-fx-font-family: monospace;");
        VBox.setVgrow(consoleOutput, Priority.ALWAYS); // Allow to expand
        consoleBox.getChildren().add(consoleOutput);
        VBox.setVgrow(consoleBox, Priority.ALWAYS);
        // Control Buttons
        HBox controlBox = new HBox(10);
        startButton = new Button("Start Analysis");
        stopButton = new Button("Stop");
        stopButton.setDisable(true);
        // TODO: Implement startAnalysis() and stopAnalysis() for pFind workflow
        startButton.setOnAction(e -> { startAnalysis(); });
        stopButton.setOnAction(e -> stopAnalysis());
        
        // Debug button for testing metadata flow
//        Button testButton = new Button("Test Metadata");
//       testButton.setOnAction(e -> testMetadataFlow());
//        testButton.setStyle("-fx-background-color: #ffeb3b; -fx-text-fill: #000;");
        
        controlBox.getChildren().addAll(startButton, stopButton);
        // Status Panel (bottom)
        HBox statusPanel = new HBox(20);
        statusPanel.setPadding(new Insets(5, 0, 0, 0));
        statusPanel.setAlignment(Pos.CENTER_LEFT);
        statusTimeLabel = new Label("Time: 00:00:00");
        statusCpuLabel = new Label("CPU: 0.0%");
        statusMemLabel = new Label("Memory: 0 / 0 MB");
        statusPanel.getChildren().addAll(statusTimeLabel, statusCpuLabel, statusMemLabel);
        // Layout: parameterView and consoleOutput fill available space, statusPanel at bottom
        tabContent.getChildren().addAll(
            paramViewBox,
            progressBox,
            consoleBox,
            controlBox
        );
        VBox.setVgrow(paramViewBox, Priority.ALWAYS);
        VBox.setVgrow(consoleBox, Priority.ALWAYS);
        VBox.setVgrow(tabContent, Priority.ALWAYS);
        // Always keep status panel at the bottom
        tabContent.getChildren().add(statusPanel);
        return tabContent;
    }
    private void updateMetadataTableFromIOPanel() {
        // ... (copy/adapt from MaxQuant panel)
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
    private void loadModificationsFromIni() {
        try (InputStream in = getClass().getClassLoader().getResourceAsStream("modification_default.ini")) {
            if (in != null) {
                BufferedReader reader = new BufferedReader(new InputStreamReader(in));
                Pattern pattern = Pattern.compile("name\\d+=(.+?) ");
                String line;
                while ((line = reader.readLine()) != null) {
                    Matcher m = pattern.matcher(line);
                    if (m.find()) {
                        String modName = m.group(1).trim();
                        modifications.add(modName);
                    }
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
    private void loadIsobaricBaseTypes() {
        IsobaricTag[] tags = IsobaricTag.getPFindTags();
        for (IsobaricTag tag : tags) {
            isobaricBaseTypes.add(tag.getName());
        }
    }
    private Label createSectionTitle(String title) {
        Label label = new Label(title);
        label.setFont(Font.font("System", FontWeight.BOLD, 14));
        label.setTextFill(Color.DARKBLUE);
        return label;
    }
    // --- END: Methods adapted from MetaPilotDdaIgcMaxQuantPanel ---

    // --- BEGIN: Adapted methods from MaxQuant panel ---
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
        missedCleavagesBox.getChildren().add(missedLabel);
        missedCleavagesSpinner = new Spinner<>(0, 5, 2);
        missedCleavagesSpinner.setMinWidth(100);
        missedCleavagesSpinner.setPrefWidth(100);
        missedCleavagesSpinner.setMaxWidth(Double.MAX_VALUE);
        missedCleavagesSpinner.valueProperty().addListener((obs, oldVal, newVal) -> updateParameterView());
        missedCleavagesBox.getChildren().add(missedCleavagesSpinner);
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
        threadsSpinner = new Spinner<>(1, Runtime.getRuntime().availableProcessors(), 4);
        threadsSpinner.setMinWidth(100);
        threadsSpinner.setPrefWidth(100);
        threadsSpinner.setMaxWidth(Double.MAX_VALUE);
        threadsSpinner.valueProperty().addListener((obs, oldVal, newVal) -> updateParameterView());
        VBox threadsBox = new VBox(5);
        threadsBox.getChildren().addAll(threadsLabel, threadsSpinner);
        Label ms2Label = new Label("MS2 Scan Mode:");
        ms2Label.setStyle("-fx-font-weight: bold;");
        ms2ScanModeComboBox = new ComboBox<>(FXCollections.observableArrayList("HCD-FTMS", "HCD-ITMS", "CID-FTMS", "CID-ITMS"));
        ms2ScanModeComboBox.setValue("HCD-FTMS");
        ms2ScanModeComboBox.setMinWidth(100);
        ms2ScanModeComboBox.setPrefWidth(120);
        ms2ScanModeComboBox.setMaxWidth(Double.MAX_VALUE);
        ms2ScanModeComboBox.setOnAction(e -> updateParameterView());
        VBox ms2Box = new VBox(5);
        ms2Box.getChildren().addAll(ms2Label, ms2ScanModeComboBox);
        paramGrid.add(threadsBox, 0, 0);
        paramGrid.add(ms2Box, 1, 0);
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
        ColumnConstraints col1 = new ColumnConstraints();
        col1.setMinWidth(200);
        col1.setHgrow(Priority.ALWAYS);
        ColumnConstraints col2 = new ColumnConstraints();
        col2.setMinWidth(200);
        col2.setHgrow(Priority.ALWAYS);
        modGrid.getColumnConstraints().addAll(col1, col2);
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
        modGrid.add(fixedModLabel, 0, 0);
        modGrid.add(variableModLabel, 1, 0);
        modGrid.add(fixedModListView, 0, 1);
        modGrid.add(variableModListView, 1, 1);
        modGrid.add(removeFixedModButton, 0, 2);
        modGrid.add(removeVariableModButton, 1, 2);
        modGrid.add(modSelectionContainer, 0, 3, 2, 1);
        section.getChildren().addAll(titleLabel, modGrid);
        return section;
    }
    private VBox createMetadataConfigurationSection() {
        VBox section = new VBox(10);
        section.setPadding(new Insets(10));
        section.setStyle("-fx-background-color: #fff8f0; -fx-background-radius: 5;");
        section.setMinHeight(200);
        VBox.setVgrow(section, Priority.ALWAYS);
        Label titleLabel = createSectionTitle("Metadata Configuration");
        Label infoLabel = new Label("Metadata will be configured based on selected quantitative mode");
        infoLabel.setStyle("-fx-font-style: italic; -fx-text-fill: #666666;");
        metadataTable = new TableView<>(metadataTableData);
        metadataTable.setMinHeight(120);
        metadataTable.setPrefHeight(150);
        metadataTable.setMaxHeight(Double.MAX_VALUE);
        metadataTable.setEditable(true);
        metadataTable.setColumnResizePolicy(TableView.CONSTRAINED_RESIZE_POLICY);
        VBox.setVgrow(metadataTable, Priority.ALWAYS);
        createMetadataTableColumns();
        HBox controlBox = new HBox(8);
        controlBox.setAlignment(Pos.CENTER_LEFT);
        metadataColumnComboBox = new ComboBox<>(FXCollections.observableArrayList(
            "Age", "Gender", "Condition", "Treatment", "Time Point", "BMI", "Custom"));
        metadataColumnComboBox.setMinWidth(100);
        metadataColumnComboBox.setPrefWidth(120);
        metadataColumnComboBox.setPromptText("Add column");
        TextField customColumnTextField = new TextField();
        customColumnTextField.setPromptText("Enter custom column name");
        customColumnTextField.setMinWidth(120);
        customColumnTextField.setVisible(false);
        customColumnTextField.setManaged(false);
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
        controlBox.getChildren().addAll(metadataColumnComboBox, customColumnTextField, addMetadataColumnButton, removeMetadataColumnButton);
        section.getChildren().addAll(titleLabel, infoLabel, metadataTable, controlBox);
        updateMetadataTable();
        return section;
    }
    private void updateParameterView() {
        StringBuilder sb = new StringBuilder();
        // Enzyme settings
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
            if (selected.getText().contains("Isobaric")) {
                sb.append("Isobaric Type: ").append(isobaricTypeComboBox.getValue() != null ? isobaricTypeComboBox.getValue() : "None").append("\n");
                sb.append("Configured Labels: ").append(isobaricTableData.size()).append("\n");
            }
        }
        // Taxonomic analysis
        sb.append("\n").append(taxonomicAnalysisPanel.getTaxonomicSummary());
        parameterView.setText(sb.toString());
    }
    private void updateMetadataTable() {
        metadataTableData.clear();
        ObservableList<MetaPilotIOPanelFX.ExperimentRow> fileRows = ioPanel.getIoPanel().getExperimentData();
        if (fileRows.isEmpty()) return;
        createMetadataTableColumns();
        if (quantModeGroup.getSelectedToggle() != null) {
            RadioButton selected = (RadioButton) quantModeGroup.getSelectedToggle();
            String mode = selected.getText();
            boolean hasLabels = false;
            if (mode.contains("Isobaric")) {
                hasLabels = !isobaricTableData.isEmpty();
            }
            if (mode.contains("Isobaric") && hasLabels) {
                for (MetaPilotIOPanelFX.ExperimentRow row : fileRows) {
                    for (IsobaricLabelingPanel.IsobaricLabelRow isoLabel : isobaricTableData) {
                        MetaPilotDdaIgcMaxQuantPanel.MetadataRow metadataRow = new MetaPilotDdaIgcMaxQuantPanel.MetadataRow(row.getExperimentName());
                        metadataRow.setLabelType("Isobaric - " + isoLabel.getLabelName());
                        metadataTableData.add(metadataRow);
                    }
                }
            } else {
                for (MetaPilotIOPanelFX.ExperimentRow row : fileRows) {
                    metadataTableData.add(new MetaPilotDdaIgcMaxQuantPanel.MetadataRow(row.getExperimentName()));
                }
            }
        }
    }
    private void updateIsobaricTable() {
        String selectedType = isobaricTypeComboBox.getValue();
        isobaricTableData.clear();
        
        if (selectedType != null) {
            // Find the IsobaricTag that matches the selected type
            IsobaricTag selectedTag = null;
            IsobaricTag[] tags = IsobaricTag.getPFindTags();
            for (IsobaricTag tag : tags) {
                if (tag.getName().equals(selectedType)) {
                    selectedTag = tag;
                    break;
                }
            }
            
            if (selectedTag != null) {
                // Use the getReportTitle() method to get the proper labels
                String[] reportTitles = selectedTag.getReportTitle();
                for (String title : reportTitles) {
                    isobaricTableData.add(new IsobaricLabelingPanel.IsobaricLabelRow(title));
                }
            }
            updateMetadataTable();
        }
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
    // --- END: Adapted methods from MaxQuant panel ---

    private void createMetadataTableColumns() {
        metadataTable.getColumns().clear();
        TableColumn<MetaPilotDdaIgcMaxQuantPanel.MetadataRow, String> sampleCol = new TableColumn<>("Sample");
        sampleCol.setCellValueFactory(new PropertyValueFactory<>("sampleName"));
        sampleCol.setMinWidth(100);
        sampleCol.setPrefWidth(120);
        metadataTable.getColumns().add(sampleCol);
        if (quantModeGroup.getSelectedToggle() != null) {
            RadioButton selected = (RadioButton) quantModeGroup.getSelectedToggle();
            String mode = selected.getText();
            if (mode.contains("Isobaric")) {
                TableColumn<MetaPilotDdaIgcMaxQuantPanel.MetadataRow, String> labelCol = new TableColumn<>("Label Type");
                labelCol.setCellValueFactory(new PropertyValueFactory<>("labelType"));
                labelCol.setMinWidth(120);
                labelCol.setPrefWidth(150);
                metadataTable.getColumns().add(labelCol);
            }
        }
        addMetadataColumnWithProperty("Age", "age");
        addMetadataColumnWithProperty("Gender", "gender");
        addMetadataColumnWithProperty("Condition", "condition");
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
        if (metadataTable.getColumns().size() > 1) {
            metadataTable.getColumns().remove(metadataTable.getColumns().size() - 1);
        }
    }
    
    private void addMetadataColumnWithProperty(String columnName, String propertyName) {
        TableColumn<MetaPilotDdaIgcMaxQuantPanel.MetadataRow, String> column = new TableColumn<>(columnName);
        column.setCellValueFactory(new PropertyValueFactory<>(propertyName));
        column.setCellFactory(TextFieldTableCell.forTableColumn(new DefaultStringConverter()));
        column.setOnEditCommit(event -> {
            MetaPilotDdaIgcMaxQuantPanel.MetadataRow row = event.getRowValue();
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
        });
        column.setMinWidth(80);
        column.setPrefWidth(100);
        metadataTable.getColumns().add(column);
    }

    private void startAnalysis() {
        // Disable start, enable stop
        startButton.setDisable(true);
        stopButton.setDisable(false);
        consoleOutput.clear();
        parameterView.setText("Starting analysis...\n" + parameterView.getText());

        // Gather parameters from UI
        ObservableList<MetaPilotIOPanelFX.ExperimentRow> fileRows = ioPanel.getIoPanel().getExperimentData();
        if (fileRows.isEmpty()) {
            showAlert("Error", "No files selected for analysis.");
            return;
        }
        
        LOGGER.info("Starting pFind analysis with {} experiment rows", fileRows.size());
        for (int i = 0; i < fileRows.size(); i++) {
            MetaPilotIOPanelFX.ExperimentRow row = fileRows.get(i);
            LOGGER.info("Experiment row {}: file={}, name={}", i, row.getFilePath(), row.getExperimentName());
        }
        
        // Save parameters for this analysis and as last used
        DDAIgcPFindParameters analysisParams = collectParametersFromUI();
        
        // Debug logging for the collected parameters
        LOGGER.info("Collected parameters object: {}", analysisParams);
        if (analysisParams.getMetadata() != null) {
            MetaData metadata = analysisParams.getMetadata();
            String[] rawFiles = metadata.getRawFiles();
            String[] expNames = metadata.getExpNames();
            LOGGER.info("Parameters metadata contains {} raw files", rawFiles != null ? rawFiles.length : 0);
            if (rawFiles != null) {
                for (int i = 0; i < rawFiles.length; i++) {
                    LOGGER.info("Parameters raw file {}: {} -> {}", i, rawFiles[i], expNames != null && i < expNames.length ? expNames[i] : "N/A");
                }
            }
        } else {
            LOGGER.warn("Parameters metadata is null!");
        }
        
        String outputDir = ioPanel.getIoPanel().getOutputDirectory();
        if (outputDir != null && !outputDir.trim().isEmpty()) {
            // Save analysis-specific parameters
            parameterIO.saveAnalysisParameters(outputDir, analysisParams);
            
            // Save as last used parameters
            parameterIO.saveLastUsedParameters(analysisParams);
            
            // Also save to preferences for quick defaults
            parameterIO.savePreferences(analysisParams);
            
            LOGGER.info("Saved pFind analysis parameters to {}/analysis_parameters_pfind.json", outputDir);
        }
        
        // Note: We now use DDAIgcPFindParameters directly instead of MetaParameterPFind
        // The analysisParams object already contains all the necessary metadata and parameters

        // Note: We now use DDAIgcPFindParameters directly instead of MetaParameterPFind

        // Load resource parameters
        uk.ac.quadram.metapilot.resources.MetaResourceParameters mrp = (new uk.ac.quadram.metapilot.resources.MetaResourceIO()).loadParameters();

        // Progress bars
        ProgressBar overallBar = new ProgressBar();
        ProgressBar currentBar = new ProgressBar();

        // Create and run the DDAIgcPFindTask
        LOGGER.info("Creating DDAIgcPFindTask with analysisParams: {}", analysisParams);
        DDAIgcPFindTask task = new DDAIgcPFindTask(analysisParams, mrp, overallBar, currentBar);
        // Bind progress and messages
        overallProgressBar.progressProperty().bind(overallBar.progressProperty());
        currentTaskProgressBar.progressProperty().bind(currentBar.progressProperty());
        task.messageProperty().addListener((obs, oldMsg, newMsg) -> {
            Platform.runLater(() -> {
                if (newMsg != null && !newMsg.isEmpty()) {
                    consoleOutput.appendText(newMsg + "\n");
                }
            });
        });
        // Status panel timer
        taskStartTimeMillis = System.currentTimeMillis();
        if (statusTimeline != null) statusTimeline.stop();
        statusTimeline = new Timeline(new KeyFrame(Duration.seconds(1), e -> updateStatusPanel()));
        statusTimeline.setCycleCount(Timeline.INDEFINITE);
        statusTimeline.play();
        // On succeeded/failed/cancelled
        task.setOnSucceeded(e -> {
            statusTimeline.stop();
            startButton.setDisable(false);
            stopButton.setDisable(true);
            updateStatusPanel();
            
            // Check if the task actually succeeded by examining the result
            try {
                Boolean result = task.get();
                if (result != null && result) {
                    consoleOutput.appendText("\nAnalysis completed successfully.\n");
                } else {
                    consoleOutput.appendText("\nAnalysis completed but may not have been fully successful.\n");
                }
            } catch (Exception ex) {
                consoleOutput.appendText("\nAnalysis failed: " + ex.getMessage() + "\n");
            }
        });
        task.setOnFailed(e -> {
            statusTimeline.stop();
            startButton.setDisable(false);
            stopButton.setDisable(true);
            updateStatusPanel();
            consoleOutput.appendText("\nAnalysis failed: " + task.getException() + "\n");
        });
        task.setOnCancelled(e -> {
            statusTimeline.stop();
            startButton.setDisable(false);
            stopButton.setDisable(true);
            updateStatusPanel();
            consoleOutput.appendText("\nAnalysis cancelled.\n");
        });
        // Run the task
        Thread t = new Thread(task);
        t.setDaemon(true);
        t.start();
    }

    private void updateStatusPanel() {
        // Elapsed time
        long elapsedMillis = System.currentTimeMillis() - taskStartTimeMillis;
        long seconds = elapsedMillis / 1000;
        long hh = seconds / 3600;
        long mm = (seconds % 3600) / 60;
        long ss = seconds % 60;
        String timeStr = String.format("Time: %02d:%02d:%02d", hh, mm, ss);
        // CPU usage (system-wide)
        double cpuLoad = 0.0;
        try {
            java.lang.management.OperatingSystemMXBean osBean = ManagementFactory.getOperatingSystemMXBean();
            if (osBean instanceof com.sun.management.OperatingSystemMXBean) {
                cpuLoad = ((com.sun.management.OperatingSystemMXBean) osBean).getSystemCpuLoad();
                if (cpuLoad >= 0) cpuLoad *= 100.0;
                else cpuLoad = 0.0;
            }
        } catch (Throwable t) { cpuLoad = 0.0; }
        String cpuStr = String.format("CPU: %.1f%%", cpuLoad);
        // Memory usage (heap)
        long usedMem = 0, maxMem = 0;
        try {
            MemoryMXBean memBean = ManagementFactory.getMemoryMXBean();
            MemoryUsage heap = memBean.getHeapMemoryUsage();
            usedMem = heap.getUsed() / (1024 * 1024);
            maxMem = heap.getMax() / (1024 * 1024);
        } catch (Throwable t) { usedMem = 0; maxMem = 0; }
        String memStr = String.format("Memory: %d / %d MB", usedMem, maxMem);
        Platform.runLater(() -> {
            statusTimeLabel.setText(timeStr);
            statusCpuLabel.setText(cpuStr);
            statusMemLabel.setText(memStr);
        });
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
    private DDAIgcPFindParameters collectParametersFromUI() {
        // Collect required parameters for DDAIgcPFindParameters constructor
        String result = ioPanel.getIoPanel().getOutputDirectory();
        String microDb = ioPanel.getMicrobiomeDBPath();
        boolean appendHostDb = ioPanel.isHostDBEnabled();
        String hostDb = ioPanel.getHostDBPath();
        String enzyme = enzymeComboBox.getValue();
        int missCleavages = missedCleavagesSpinner.getValue();
        String digestMode = digestionModeComboBox.getValue();
        int threadCount = threadsSpinner.getValue();
        String ms2ScanMode = ms2ScanModeComboBox.getValue();
        boolean reduceDb = generateSampleSpecificDBCheckBox.isSelected();
        boolean spectraCluster = spectraClusteringCheckBox.isSelected();

        // Create MetaData object from experiment data
        ObservableList<MetaPilotIOPanelFX.ExperimentRow> fileRows = ioPanel.getIoPanel().getExperimentData();
        LOGGER.info("collectParametersFromUI: Found {} experiment rows", fileRows.size());
        
        MetaData metaData;
        if (!fileRows.isEmpty()) {
            String[] rawFiles = new String[fileRows.size()];
            String[] expNames = new String[fileRows.size()];
            int[] fractions = new int[fileRows.size()];
            int[] replicates = new int[fileRows.size()];
            
            for (int i = 0; i < fileRows.size(); i++) {
                MetaPilotIOPanelFX.ExperimentRow row = fileRows.get(i);
                rawFiles[i] = row.getFilePath();
                expNames[i] = row.getExperimentName();
                fractions[i] = 0;
                replicates[i] = 1;
                LOGGER.info("collectParametersFromUI: Row {} - file: {}, name: {}", i, rawFiles[i], expNames[i]);
            }
            
            LOGGER.info("collectParametersFromUI: Creating MetaData with {} raw files", rawFiles.length);
            metaData = new MetaData(rawFiles, expNames, fractions, replicates);
            LOGGER.info("collectParametersFromUI: MetaData created: {}", metaData);
            LOGGER.info("collectParametersFromUI: MetaData.getRawFiles(): {}", metaData.getRawFiles());
        } else {
            LOGGER.warn("No experiment data available, creating empty MetaData");
            metaData = new MetaData();
        }

        // Create DDAIgcPFindParameters using the correct constructor
        DDAIgcPFindParameters params = new DDAIgcPFindParameters(
            MetaPilotApp.WorkflowType.DDA_IGC_PFIND.getAnalysisType(),
            metaData,
            result,
            microDb,
            appendHostDb,
            hostDb,
            ms2ScanMode,
            threadCount
        );
        
        // Set additional parameters that aren't in the constructor
        params.setExperimentData(new ArrayList<>(fileRows));
        
        // Enzyme settings
        params.setEnzyme(enzyme);
        params.setMissCleavages(missCleavages);
        params.setDigestMode(digestMode);
        
        // Additional parameters
        params.setReduceDb(reduceDb);
        params.setSpectraCluster(spectraCluster);
        
        // Modifications
        params.setFixedModifications(new ArrayList<>(fixedModifications));
        params.setVariableModifications(new ArrayList<>(variableModifications));
        
        // Quantitative mode
        if (quantModeGroup.getSelectedToggle() != null) {
            RadioButton selected = (RadioButton) quantModeGroup.getSelectedToggle();
            params.setQuantitativeMode(selected.getText());
        }
        
        // Isobaric labeling
        if (isobaricTypeComboBox.getValue() != null) {
            params.setIsobaricType(isobaricTypeComboBox.getValue());
        }
        params.setIsobaricTableData(new ArrayList<>(isobaricTableData));
        
        // Taxonomic analysis
        params.setEukaryota(taxonomicAnalysisPanel.isEukaryota());
        params.setBacteria(taxonomicAnalysisPanel.isBacteria());
        params.setArchaea(taxonomicAnalysisPanel.isArchaea());
        params.setViruses(taxonomicAnalysisPanel.isViruses());
        params.setViroids(taxonomicAnalysisPanel.isViroids());
        params.setUnclassified(taxonomicAnalysisPanel.isUnclassified());
        params.setOther(taxonomicAnalysisPanel.isOther());
        params.setIgnoreRank(taxonomicAnalysisPanel.getIgnoreRank());
        params.setIgnoreEnvironmentalSamples(taxonomicAnalysisPanel.isIgnoreEnvironmentalSamples());
        
        // Metadata table data
        params.setMetadataTableData(new ArrayList<>(metadataTableData));
        
        return params;
    }
    
    /**
     * Load last used parameters and apply them to the UI
     */
    private void loadLastUsedParameters() {
        try {
            DDAIgcPFindParameters params = (DDAIgcPFindParameters) parameterIO.loadLastUsedParameters();
            if (params != null) {
                applyParametersToUI(params);
                LOGGER.info("Loaded last used parameters for DDA IGC pFind workflow");
                
                // Also apply experiment data to the IO panel if available
                if (params.getExperimentData() != null && !params.getExperimentData().isEmpty()) {
                    ioPanel.getIoPanel().setExperimentData(params.getExperimentData());
                    LOGGER.info("Applied {} experiment rows to IO panel", params.getExperimentData().size());
                }
                
                // Apply output directory if available
                if (params.getResult() != null && !params.getResult().trim().isEmpty()) {
                    ioPanel.getIoPanel().setOutputDirectory(params.getResult());
                    LOGGER.info("Applied output directory: {}", params.getResult());
                }
                
                // Update metadata table after loading parameters
                Platform.runLater(() -> {
                    updateMetadataTable();
                    updateParameterView();
                });
            } else {
                LOGGER.warn("No parameters loaded, using defaults");
            }
        } catch (Exception e) {
            LOGGER.warn("Failed to load last used pFind parameters: " + e.getMessage(), e);
            // Continue with default UI state
        }
    }
    
    /**
     * Apply parameters to the UI components
     */
    private void applyParametersToUI(DDAIgcPFindParameters params) {
        if (params == null) return;
        
        // Apply metadata to IO panel if available
        if (params.getMetadata() != null) {
            MetaData metadata = params.getMetadata();
            String[] rawFiles = metadata.getRawFiles();
            String[] expNames = metadata.getExpNames();
            
            if (rawFiles != null && expNames != null && rawFiles.length > 0) {
                ObservableList<MetaPilotIOPanelFX.ExperimentRow> experimentRows = FXCollections.observableArrayList();
                for (int i = 0; i < rawFiles.length; i++) {
                    if (rawFiles[i] != null && !rawFiles[i].trim().isEmpty()) {
                        experimentRows.add(new MetaPilotIOPanelFX.ExperimentRow(rawFiles[i], expNames[i]));
                    }
                }
                ioPanel.getIoPanel().setExperimentData(experimentRows);
                LOGGER.info("Applied metadata with {} raw files to IO panel", rawFiles.length);
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
            if (mode.contains("Isobaric")) {
                // Find and select the isobaric radio button
                for (Toggle toggle : quantModeGroup.getToggles()) {
                    if (toggle instanceof RadioButton && ((RadioButton) toggle).getText().contains("Isobaric")) {
                        quantModeGroup.selectToggle(toggle);
                        break;
                    }
                }
            } else {
                // Find and select the label-free radio button
                for (Toggle toggle : quantModeGroup.getToggles()) {
                    if (toggle instanceof RadioButton && ((RadioButton) toggle).getText().contains("Label-free")) {
                        quantModeGroup.selectToggle(toggle);
                        break;
                    }
                }
            }
        }
        
        // Apply isobaric type if available
        if (params.getIsobaricType() != null && isobaricTypeComboBox != null) {
            isobaricTypeComboBox.setValue(params.getIsobaricType());
        }
        
        // Apply isobaric table data if available
        if (params.getIsobaricTableData() != null && !params.getIsobaricTableData().isEmpty()) {
            isobaricTableData.clear();
            isobaricTableData.addAll(params.getIsobaricTableData());
        }
        
        // Apply taxonomic settings
        if (taxonomicAnalysisPanel != null) {
            taxonomicAnalysisPanel.setEukaryota(params.isEukaryota());
            taxonomicAnalysisPanel.setBacteria(params.isBacteria());
            taxonomicAnalysisPanel.setArchaea(params.isArchaea());
            taxonomicAnalysisPanel.setViruses(params.isViruses());
            taxonomicAnalysisPanel.setViroids(params.isViroids());
            taxonomicAnalysisPanel.setUnclassified(params.isUnclassified());
            taxonomicAnalysisPanel.setOther(params.isOther());
            if (params.getIgnoreRank() != null) {
                taxonomicAnalysisPanel.setIgnoreRank(params.getIgnoreRank());
            }
            taxonomicAnalysisPanel.setIgnoreEnvironmentalSamples(params.isIgnoreEnvironmentalSamples());
        }
        
        // Update parameter view after applying parameters
        Platform.runLater(() -> updateParameterView());
    }
    
    /**
     * Load parameters from a specific analysis result directory
     */
    public void loadAnalysisParameters(String resultDirectory) {
        try {
            DDAIgcPFindParameters params = (DDAIgcPFindParameters) parameterIO.loadAnalysisParameters(resultDirectory);
            applyParametersToUI(params);
            showAlert("Success", "Loaded pFind analysis parameters from " + resultDirectory);
            LOGGER.info("Loaded pFind analysis parameters from {}", resultDirectory);
        } catch (Exception e) {
            showAlert("Error", "Failed to load pFind analysis parameters: " + e.getMessage());
            LOGGER.error("Failed to load pFind analysis parameters from " + resultDirectory, e);
        }
    }
    
    /**
     * Test method to verify metadata flow - can be called for debugging
     */
    public void testMetadataFlow() {
        LOGGER.info("=== Testing Metadata Flow ===");
        
        // Test 1: Check current experiment data
        ObservableList<MetaPilotIOPanelFX.ExperimentRow> fileRows = ioPanel.getIoPanel().getExperimentData();
        LOGGER.info("Test 1: Current experiment data has {} rows", fileRows.size());
        
        // Test 2: Collect parameters and check metadata
        DDAIgcPFindParameters testParams = collectParametersFromUI();
        LOGGER.info("Test 2: Collected parameters object: {}", testParams);
        
        if (testParams.getMetadata() != null) {
            MetaData metadata = testParams.getMetadata();
            String[] rawFiles = metadata.getRawFiles();
            LOGGER.info("Test 2: Metadata contains {} raw files", rawFiles != null ? rawFiles.length : 0);
        } else {
            LOGGER.warn("Test 2: Metadata is null in collected parameters!");
        }
        
        // Test 3: Create a task and check if metadata is preserved
        try {
            uk.ac.quadram.metapilot.resources.MetaResourceParameters mrp = (new uk.ac.quadram.metapilot.resources.MetaResourceIO()).loadParameters();
            ProgressBar testBar = new ProgressBar();
            
            LOGGER.info("Test 3: Creating test task...");
            DDAIgcPFindTask testTask = new DDAIgcPFindTask(testParams, mrp, testBar, testBar);
            LOGGER.info("Test 3: Test task created successfully");
            
        } catch (Exception e) {
            LOGGER.error("Test 3: Failed to create test task: " + e.getMessage(), e);
        }
        
        LOGGER.info("=== Metadata Flow Test Complete ===");
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
                                DDAIgcPFindParameters currentParams = collectParametersFromUI();
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
    
    private void stopAnalysis() {
        startButton.setDisable(false);
        stopButton.setDisable(true);
        consoleOutput.appendText("Stopping pFind analysis...\n");
        System.out.println("DEBUG: stopAnalysis called");
    }
} 