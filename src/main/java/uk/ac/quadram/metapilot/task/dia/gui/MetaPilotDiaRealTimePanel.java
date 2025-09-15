package uk.ac.quadram.metapilot.task.dia.gui;

import javafx.application.Platform;
import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.control.*;
import javafx.scene.layout.*;
import javafx.stage.DirectoryChooser;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.scene.text.Font;
import javafx.scene.text.FontWeight;
import javafx.scene.paint.Color;
import javafx.scene.control.cell.PropertyValueFactory;
import javafx.beans.property.SimpleStringProperty;
import javafx.beans.property.StringProperty;
import javafx.concurrent.Task;
import uk.ac.quadram.metapilot.task.common.gui.MagDatabaseSettingsPanel;
import uk.ac.quadram.metapilot.task.par.MetaData;
import uk.ac.quadram.metapilot.task.par.MetaBaseParameter;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;
import uk.ac.quadram.metapilot.task.mag.MagDbItem.MagVersion;
import uk.ac.quadram.metapilot.task.mag.par.MetaParameterMag;
import uk.ac.quadram.metapilot.task.dia.par.DiaParameterIO;
import uk.ac.quadram.metapilot.task.dia.par.MetaParameterDia;
import java.io.File;
import java.io.IOException;
import java.nio.channels.FileChannel;
import java.nio.channels.FileLock;
import java.nio.file.*;
import java.nio.file.StandardOpenOption;

public class MetaPilotDiaRealTimePanel extends VBox {
    // UI Components - Database
    private MagDatabaseSettingsPanel magDatabasePanel;
    private TextArea consoleOutput;
    private TextArea parameterView;
    private ProgressBar overallProgressBar;
    private ProgressBar currentTaskProgressBar;
    private Label currentTaskLabel;

    // Real-time monitoring components
    private TextField monitorFolderField;
    private Button browseMonitorFolderButton;
    private TextField resultFolderField;
    private Button browseResultFolderButton;
    private ListView<String> monitorFilesList;
    private TableView<ProcessingRow> resultTable;
    private ObservableList<ProcessingRow> resultData = FXCollections.observableArrayList();
    
    // Control panel components
    private Button startMonitorButton;
    private Button stopMonitorButton;
    private Button stopAfterCountButton;
    private Button resetButton;
    private ComboBox<Integer> fileCountComboBox;
    private CheckBox performAnalysisCheckBox;

    // DIA settings
    private ComboBox<Integer> threadsComboBox;

    // Separate panels for library and model generation
    private DiaLibraryGenerationPanel libraryGenerationPanel;
    private DiaModelGenerationPanel modelGenerationPanel;

    // Status panel fields
    private Label statusTimeLabel;
    private Label statusCpuLabel;
    private Label statusMemLabel;

    // Parameter I/O
    private DiaParameterIO parameterIO;
    private MetaParameterDia currentAnalysisParameters;
    
    // Resource parameters
    private MetaResourceParameters resourceParams;

    // Monitoring status
    private boolean isMonitoring = false;
    private int processedFileCount = 0;
    private int targetFileCount = 0;

    public MetaPilotDiaRealTimePanel() {
        // Initialize parameter IO instance
        parameterIO = new DiaParameterIO();
        initializeUI();
    }

    public MetaPilotDiaRealTimePanel(MetaParameterMag metaParameterMag) {
        // Initialize parameter IO instance
        parameterIO = new DiaParameterIO();
        
        initializeUI();
        
        // Initialize MAG databases from parameter after UI is initialized
        if (metaParameterMag != null) {
            MagDbItem[] magDbItems = metaParameterMag.getAvailableMagDbItem();
            if (magDbItems != null) {
                setAvailableMagDatabases(magDbItems);
            }
        }
    }
    
    public MetaPilotDiaRealTimePanel(MagDbItem[] magDbItems, MetaResourceParameters resourceParams) {
        // Initialize parameter IO instance
        parameterIO = new DiaParameterIO();
        this.resourceParams = resourceParams;
        initializeUI();
        // Set MAG databases after UI is initialized
        setAvailableMagDatabases(magDbItems);
        // Load last used parameters after UI and databases are initialized
        loadLastUsedParameters();
    }

    // Load available MAG databases
    public void setAvailableMagDatabases(MagDbItem[] magDbItems) {
        if (magDatabasePanel != null) {
            magDatabasePanel.setAvailableMagDatabases(magDbItems);
            // Update the panels with the selected MAG database
            if (magDbItems != null && magDbItems.length > 0) {
                updatePanelsWithMagDb(magDbItems[0]);
            }
        } else {
            System.out.println("DEBUG: magDatabasePanel is null - UI not initialized yet");
        }
    }

    private void initializeUI() {
        setSpacing(10);
        setPadding(new Insets(10));
        setMinWidth(1200);
        setMinHeight(800);
        setPrefWidth(1600);
        setPrefHeight(900);

        // Initialize reusable components
        magDatabasePanel = new MagDatabaseSettingsPanel();
        magDatabasePanel.setParameterUpdateCallback(() -> {
            updateParameterView();
            // Update panels when database selection changes
            MagDbItem selectedMagDb = magDatabasePanel.getSelectedMagDatabase();
            if (selectedMagDb != null) {
                updatePanelsWithMagDb(selectedMagDb);
            }
        });

        // Create tab pane with responsive sizing
        TabPane tabPane = new TabPane();
        tabPane.setTabClosingPolicy(TabPane.TabClosingPolicy.UNAVAILABLE);
        tabPane.setMinWidth(980);
        tabPane.setMinHeight(780);
        VBox.setVgrow(tabPane, Priority.ALWAYS);

        // Create tabs: Parameters, Monitor (removed I/O tab for real-time workflow)
        Tab parametersTab = new Tab("Parameters");
        parametersTab.setContent(createParametersTab());

        Tab monitorTab = new Tab("Monitor");
        monitorTab.setContent(createMonitorTab());

        tabPane.getTabs().addAll(parametersTab, monitorTab);

        // Add tab change listener for parameter updates and auto-save
        tabPane.getSelectionModel().selectedItemProperty().addListener((obs, oldTab, newTab) -> {
            if (newTab == parametersTab) {
                // Update parameter view when switching to Parameters tab
                updateParameterView();
            } else if (newTab == monitorTab) {
                // Update parameter view and auto-save when switching to Monitor tab
                updateParameterView();
                
                // Save current parameters when Monitor tab is selected
                try {
                    MetaParameterDia currentParams = collectParametersFromUI();
                    if (currentParams != null) {
                        boolean saved = parameterIO.saveLastUsedParameters(currentParams);
                        if (saved) {
                            System.out.println("DEBUG: Auto-saved parameters when switching to Monitor tab");
                        }
                    }
                } catch (Exception e) {
                    System.err.println("WARNING: Failed to auto-save parameters: " + e.getMessage());
                }
            }
        });

        // Add tab pane to main panel
        getChildren().add(tabPane);
        
        // Initialize parameter view
        updateParameterView();
    }

    // Create IO tab (same as DIA Normal Panel)

    // Processing Row class for result table
    public static class ProcessingRow {
        private final SimpleStringProperty fileName;
        private final SimpleStringProperty round1;
        private final SimpleStringProperty model;
        private final SimpleStringProperty round2;
        private final SimpleStringProperty round3;
        private final SimpleStringProperty round4;
        
        public ProcessingRow(String fileName) {
            this.fileName = new SimpleStringProperty(fileName);
            this.round1 = new SimpleStringProperty("Pending");
            this.model = new SimpleStringProperty("Pending");
            this.round2 = new SimpleStringProperty("Pending");
            this.round3 = new SimpleStringProperty("Pending");
            this.round4 = new SimpleStringProperty("Pending");
        }
        
        public String getFileName() { return fileName.get(); }
        public void setFileName(String value) { fileName.set(value); }
        public StringProperty fileNameProperty() { return fileName; }
        
        public String getRound1() { return round1.get(); }
        public void setRound1(String value) { round1.set(value); }
        public StringProperty round1Property() { return round1; }
        
        public String getModel() { return model.get(); }
        public void setModel(String value) { model.set(value); }
        public StringProperty modelProperty() { return model; }
        
        public String getRound2() { return round2.get(); }
        public void setRound2(String value) { round2.set(value); }
        public StringProperty round2Property() { return round2; }
        
        public String getRound3() { return round3.get(); }
        public void setRound3(String value) { round3.set(value); }
        public StringProperty round3Property() { return round3; }
        
        public String getRound4() { return round4.get(); }
        public void setRound4(String value) { round4.set(value); }
        public StringProperty round4Property() { return round4; }
    }

    private VBox createParametersTab() {
        VBox tabContent = new VBox(10);
        tabContent.setPadding(new Insets(10));
        tabContent.setMinWidth(960);
        tabContent.setMinHeight(760);

        // Add database settings at the top
        tabContent.getChildren().addAll(
            createSectionTitle("Database Settings"),
            magDatabasePanel
        );

        // Create main grid layout for parameters with 2-column responsive sizing
        GridPane mainGrid = new GridPane();
        mainGrid.setHgap(15);
        mainGrid.setVgap(15);
        mainGrid.setPadding(new Insets(10));
        mainGrid.setMinWidth(940);
        mainGrid.setMinHeight(600); // Reduced height since database settings are now above
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

        // Row constraints - 2 rows: settings and panels
        RowConstraints row1 = new RowConstraints();
        row1.setMinHeight(60);
        row1.setVgrow(Priority.NEVER);
        RowConstraints row2 = new RowConstraints();
        row2.setVgrow(Priority.ALWAYS);
        row2.setMinHeight(400);
        mainGrid.getRowConstraints().addAll(row1, row2);

        // Thread Count Settings (spans both columns)
        VBox settingsSection = createThreadCountSection();
        mainGrid.add(settingsSection, 0, 0, 2, 1);

        // Row 2: Library Generation Panel (left) | Model Generation Panel (right)
        libraryGenerationPanel = new DiaLibraryGenerationPanel(consoleOutput);
        libraryGenerationPanel.setParameterUpdateCallback(() -> updateParameterView());
        mainGrid.add(libraryGenerationPanel, 0, 1);
        GridPane.setVgrow(libraryGenerationPanel, Priority.ALWAYS);

        modelGenerationPanel = new DiaModelGenerationPanel(consoleOutput);
        modelGenerationPanel.setParameterUpdateCallback(() -> updateParameterView());
        mainGrid.add(modelGenerationPanel, 1, 1);
        GridPane.setVgrow(modelGenerationPanel, Priority.ALWAYS);
        
        tabContent.getChildren().add(mainGrid);
        return tabContent;
    }

    private VBox createMonitorTab() {
        VBox tabContent = new VBox(10);
        tabContent.setPadding(new Insets(10));
        tabContent.setMinWidth(960);
        tabContent.setMinHeight(760);

        // Create main grid layout with 2 columns
        GridPane mainGrid = new GridPane();
        mainGrid.setHgap(15);
        mainGrid.setVgap(15);
        mainGrid.setPadding(new Insets(10));
        mainGrid.setMinWidth(940);
        mainGrid.setMinHeight(740);
        VBox.setVgrow(mainGrid, Priority.ALWAYS);

        // Column constraints for 2-column layout
        ColumnConstraints col1 = new ColumnConstraints();
        col1.setPercentWidth(60); // Monitor controls on left
        col1.setMinWidth(550);
        col1.setHgrow(Priority.ALWAYS);
        ColumnConstraints col2 = new ColumnConstraints();
        col2.setPercentWidth(40); // Console on right
        col2.setMinWidth(350);
        col2.setHgrow(Priority.ALWAYS);
        mainGrid.getColumnConstraints().addAll(col1, col2);

        // Row constraint
        RowConstraints row1 = new RowConstraints();
        row1.setVgrow(Priority.ALWAYS);
        row1.setMinHeight(720);
        mainGrid.getRowConstraints().add(row1);

        // Left side: Monitor controls
        VBox monitorControls = createMonitorControlsPanel();
        mainGrid.add(monitorControls, 0, 0);
        GridPane.setVgrow(monitorControls, Priority.ALWAYS);

        // Right side: Console content
        VBox consolePanel = createConsolePanel();
        mainGrid.add(consolePanel, 1, 0);
        GridPane.setVgrow(consolePanel, Priority.ALWAYS);

        tabContent.getChildren().add(mainGrid);
        return tabContent;
    }

    private VBox createMonitorControlsPanel() {
        VBox controlsPanel = new VBox(15);
        controlsPanel.setPadding(new Insets(0));
        VBox.setVgrow(controlsPanel, Priority.ALWAYS);

        // Create sections for monitor controls
        VBox folderSection = createFolderSection();
        VBox monitorSection = createMonitorSection();
        VBox resultSection = createResultSection();
        VBox controlSection = createControlSection();

        // Adjust heights for the 4-section layout
        folderSection.setMinHeight(80);
        folderSection.setMaxHeight(100);
        
        monitorSection.setMinHeight(180);
        monitorSection.setPrefHeight(220);
        VBox.setVgrow(monitorSection, Priority.ALWAYS);
        
        resultSection.setMinHeight(180);
        resultSection.setPrefHeight(220);
        VBox.setVgrow(resultSection, Priority.ALWAYS);
        
        controlSection.setMinHeight(100);
        controlSection.setMaxHeight(120);

        controlsPanel.getChildren().addAll(folderSection, monitorSection, resultSection, controlSection);
        return controlsPanel;
    }

    private VBox createConsolePanel() {
        VBox consolePanel = new VBox(10);
        consolePanel.setPadding(new Insets(10));
        consolePanel.setStyle("-fx-background-color: #f8f8f8; -fx-background-radius: 5;");
        VBox.setVgrow(consolePanel, Priority.ALWAYS);

        Label consoleTitle = createSectionTitle("Console & Progress");

        // Parameter View
        VBox paramViewBox = new VBox(5);
        Label paramLabel = new Label("Parameter View:");
        paramLabel.setStyle("-fx-font-weight: bold; -fx-font-size: 11px;");
        parameterView = new TextArea();
        parameterView.setEditable(false);
        parameterView.setPrefHeight(150);
        parameterView.setMaxHeight(180);
        parameterView.setStyle("-fx-font-family: monospace; -fx-font-size: 10px;");
        paramViewBox.getChildren().addAll(paramLabel, parameterView);

        // Progress Bars
        VBox progressBox = new VBox(5);
        currentTaskLabel = new Label("Current Task: None");
        currentTaskLabel.setStyle("-fx-font-size: 11px;");
        overallProgressBar = new ProgressBar(0);
        overallProgressBar.setMaxWidth(Double.MAX_VALUE);
        currentTaskProgressBar = new ProgressBar(0);
        currentTaskProgressBar.setMaxWidth(Double.MAX_VALUE);
        
        Label overallLabel = new Label("Overall Progress:");
        overallLabel.setStyle("-fx-font-size: 11px;");
        Label currentLabel = new Label("Current Task Progress:");
        currentLabel.setStyle("-fx-font-size: 11px;");
        
        progressBox.getChildren().addAll(
            currentTaskLabel, overallLabel, overallProgressBar, currentLabel, currentTaskProgressBar
        );

        // Console Output
        VBox consoleBox = new VBox(5);
        Label consoleLabel = new Label("Console Output:");
        consoleLabel.setStyle("-fx-font-weight: bold; -fx-font-size: 11px;");
        consoleOutput = new TextArea();
        consoleOutput.setEditable(false);
        consoleOutput.setPrefHeight(200);
        consoleOutput.setMaxHeight(Double.MAX_VALUE);
        consoleOutput.setStyle("-fx-font-family: monospace; -fx-font-size: 10px;");
        VBox.setVgrow(consoleOutput, Priority.ALWAYS);
        consoleBox.getChildren().addAll(consoleLabel, consoleOutput);
        VBox.setVgrow(consoleBox, Priority.ALWAYS);

        // Status Panel
        HBox statusPanel = new HBox(15);
        statusPanel.setPadding(new Insets(5, 0, 0, 0));
        statusPanel.setAlignment(Pos.CENTER_LEFT);
        statusTimeLabel = new Label("Time: 00:00:00");
        statusCpuLabel = new Label("CPU: 0.0%");
        statusMemLabel = new Label("Memory: 0 / 0 MB");
        statusTimeLabel.setStyle("-fx-font-size: 10px;");
        statusCpuLabel.setStyle("-fx-font-size: 10px;");
        statusMemLabel.setStyle("-fx-font-size: 10px;");
        statusPanel.getChildren().addAll(statusTimeLabel, statusCpuLabel, statusMemLabel);

        consolePanel.getChildren().addAll(
            consoleTitle, paramViewBox, progressBox, consoleBox, statusPanel
        );
        
        VBox.setVgrow(consoleBox, Priority.ALWAYS);
        return consolePanel;
    }

    private VBox createFolderSection() {
        VBox section = new VBox(10);
        section.setPadding(new Insets(10));
        section.setStyle("-fx-background-color: #f0f8ff; -fx-background-radius: 5;");

        Label titleLabel = createSectionTitle("Folder Configuration");

        GridPane folderGrid = new GridPane();
        folderGrid.setHgap(15);
        folderGrid.setVgap(10);

        // Monitor folder
        Label monitorLabel = new Label("Monitor Folder:");
        monitorLabel.setStyle("-fx-font-weight: bold;");
        folderGrid.add(monitorLabel, 0, 0);

        monitorFolderField = new TextField();
        monitorFolderField.setPromptText("Select folder to monitor for new files");
        monitorFolderField.setPrefWidth(400);
        folderGrid.add(monitorFolderField, 1, 0);

        browseMonitorFolderButton = new Button("Browse");
        browseMonitorFolderButton.setOnAction(e -> browseMonitorFolder());
        folderGrid.add(browseMonitorFolderButton, 2, 0);

        // Result folder
        Label resultLabel = new Label("Result Folder:");
        resultLabel.setStyle("-fx-font-weight: bold;");
        folderGrid.add(resultLabel, 0, 1);

        resultFolderField = new TextField();
        resultFolderField.setPromptText("Select folder to save analysis results");
        resultFolderField.setPrefWidth(400);
        folderGrid.add(resultFolderField, 1, 1);

        browseResultFolderButton = new Button("Browse");
        browseResultFolderButton.setOnAction(e -> browseResultFolder());
        folderGrid.add(browseResultFolderButton, 2, 1);

        section.getChildren().addAll(titleLabel, folderGrid);
        return section;
    }

    private VBox createMonitorSection() {
        VBox section = new VBox(10);
        section.setPadding(new Insets(10));
        section.setStyle("-fx-background-color: #fff8dc; -fx-background-radius: 5;");
        VBox.setVgrow(section, Priority.ALWAYS);

        Label titleLabel = createSectionTitle("Monitor Folder Files");

        monitorFilesList = new ListView<>();
        monitorFilesList.setMinHeight(150);
        monitorFilesList.setPrefHeight(200);
        monitorFilesList.setMaxHeight(Double.MAX_VALUE);
        VBox.setVgrow(monitorFilesList, Priority.ALWAYS);

        section.getChildren().addAll(titleLabel, monitorFilesList);
        return section;
    }

    private VBox createResultSection() {
        VBox section = new VBox(10);
        section.setPadding(new Insets(10));
        section.setStyle("-fx-background-color: #f0fff0; -fx-background-radius: 5;");
        VBox.setVgrow(section, Priority.ALWAYS);

        Label titleLabel = createSectionTitle("Processing Results");

        resultTable = new TableView<>(resultData);
        resultTable.setMinHeight(150);
        resultTable.setPrefHeight(200);
        resultTable.setMaxHeight(Double.MAX_VALUE);
        resultTable.setColumnResizePolicy(TableView.CONSTRAINED_RESIZE_POLICY_FLEX_LAST_COLUMN);
        VBox.setVgrow(resultTable, Priority.ALWAYS);

        // File Name column
        TableColumn<ProcessingRow, String> fileCol = new TableColumn<>("File");
        fileCol.setCellValueFactory(new PropertyValueFactory<>("fileName"));
        fileCol.setMinWidth(150);
        fileCol.setPrefWidth(200);

        // Round 1 column
        TableColumn<ProcessingRow, String> round1Col = new TableColumn<>("Round 1");
        round1Col.setCellValueFactory(new PropertyValueFactory<>("round1"));
        round1Col.setMinWidth(80);
        round1Col.setPrefWidth(100);

        // Model column
        TableColumn<ProcessingRow, String> modelCol = new TableColumn<>("Model");
        modelCol.setCellValueFactory(new PropertyValueFactory<>("model"));
        modelCol.setMinWidth(80);
        modelCol.setPrefWidth(100);

        // Round 2 column
        TableColumn<ProcessingRow, String> round2Col = new TableColumn<>("Round 2");
        round2Col.setCellValueFactory(new PropertyValueFactory<>("round2"));
        round2Col.setMinWidth(80);
        round2Col.setPrefWidth(100);

        // Round 3 column
        TableColumn<ProcessingRow, String> round3Col = new TableColumn<>("Round 3");
        round3Col.setCellValueFactory(new PropertyValueFactory<>("round3"));
        round3Col.setMinWidth(80);
        round3Col.setPrefWidth(100);

        // Round 4 column
        TableColumn<ProcessingRow, String> round4Col = new TableColumn<>("Round 4");
        round4Col.setCellValueFactory(new PropertyValueFactory<>("round4"));
        round4Col.setMinWidth(80);
        round4Col.setPrefWidth(100);

        resultTable.getColumns().clear();
        resultTable.getColumns().add(fileCol);
        resultTable.getColumns().add(round1Col);
        resultTable.getColumns().add(modelCol);
        resultTable.getColumns().add(round2Col);
        resultTable.getColumns().add(round3Col);
        resultTable.getColumns().add(round4Col);

        section.getChildren().addAll(titleLabel, resultTable);
        return section;
    }

    private VBox createControlSection() {
        VBox section = new VBox(10);
        section.setPadding(new Insets(10));
        section.setStyle("-fx-background-color: #ffefd5; -fx-background-radius: 5;");

        Label titleLabel = createSectionTitle("Control Panel");

        GridPane controlGrid = new GridPane();
        controlGrid.setHgap(15);
        controlGrid.setVgap(10);

        // Row 1: Start button and combine checkbox
        startMonitorButton = new Button("Start Monitor");
        startMonitorButton.setPrefWidth(120);
        startMonitorButton.setOnAction(e -> startMonitoring());
        controlGrid.add(startMonitorButton, 0, 0);

        performAnalysisCheckBox = new CheckBox("Perform taxonomic and functional analysis after monitoring is completed");
        controlGrid.add(performAnalysisCheckBox, 1, 0, 3, 1);

        // Row 2: Stop button and stop after count
        stopMonitorButton = new Button("Stop Monitor");
        stopMonitorButton.setPrefWidth(120);
        stopMonitorButton.setDisable(true);
        stopMonitorButton.setOnAction(e -> stopMonitoring());
        controlGrid.add(stopMonitorButton, 0, 1);

        stopAfterCountButton = new Button("Stop after files reached:");
        stopAfterCountButton.setPrefWidth(150);
        stopAfterCountButton.setDisable(true);
        stopAfterCountButton.setOnAction(e -> stopAfterCount());
        controlGrid.add(stopAfterCountButton, 1, 1);

        fileCountComboBox = new ComboBox<>();
        for (int i = 1; i < 1000; i++) {
            fileCountComboBox.getItems().add(i);
        }
        fileCountComboBox.setValue(10);
        fileCountComboBox.setPrefWidth(80);
        controlGrid.add(fileCountComboBox, 2, 1);

        resetButton = new Button("Reset");
        resetButton.setPrefWidth(80);
        resetButton.setDisable(true);
        resetButton.setOnAction(e -> resetMonitoring());
        controlGrid.add(resetButton, 3, 1);

        section.getChildren().addAll(titleLabel, controlGrid);
        return section;
    }

    private void browseMonitorFolder() {
        DirectoryChooser chooser = new DirectoryChooser();
        chooser.setTitle("Select Monitor Folder");
        File directory = chooser.showDialog(getScene().getWindow());
        if (directory != null) {
            monitorFolderField.setText(directory.getAbsolutePath());
            updateParameterView();
        }
    }

    private void browseResultFolder() {
        DirectoryChooser chooser = new DirectoryChooser();
        chooser.setTitle("Select Result Folder");
        File directory = chooser.showDialog(getScene().getWindow());
        if (directory != null) {
            resultFolderField.setText(directory.getAbsolutePath());
            updateParameterView();
        }
    }

    private void startMonitoring() {
        if (!validateFolders()) {
            return;
        }

        isMonitoring = true;
        processedFileCount = 0;
        
        startMonitorButton.setDisable(true);
        stopMonitorButton.setDisable(false);
        stopAfterCountButton.setDisable(false);
        resetButton.setDisable(false);

        // Clear previous results
        monitorFilesList.getItems().clear();
        resultData.clear();

        consoleOutput.appendText("Starting real-time monitoring...\n");
        consoleOutput.appendText("Monitor folder: " + monitorFolderField.getText() + "\n");
        consoleOutput.appendText("Result folder: " + resultFolderField.getText() + "\n");

        // Start monitoring task
        startMonitoringTask();
    }

    private void stopMonitoring() {
        isMonitoring = false;
        
        stopMonitorButton.setDisable(true);
        
        consoleOutput.appendText("Stopping monitoring after current files are processed...\n");
        
        // Create stop signal file
        createStopSignalFile(false);
    }

    private void stopAfterCount() {
        targetFileCount = fileCountComboBox.getValue();
        
        stopAfterCountButton.setDisable(true);
        
        consoleOutput.appendText("Will stop monitoring after " + targetFileCount + " files are processed...\n");
        
        // Create stop signal file with count
        createStopSignalFile(true);
    }

    private void resetMonitoring() {
        isMonitoring = false;
        processedFileCount = 0;
        targetFileCount = 0;
        
        startMonitorButton.setDisable(false);
        stopMonitorButton.setDisable(true);
        stopAfterCountButton.setDisable(true);
        resetButton.setDisable(true);
        
        fileCountComboBox.setValue(10);
        
        consoleOutput.appendText("Monitoring reset.\n");
    }

    private boolean validateFolders() {
        if (monitorFolderField.getText() == null || monitorFolderField.getText().trim().isEmpty()) {
            showError("Validation Error", "Please select a monitor folder.");
            return false;
        }
        
        if (resultFolderField.getText() == null || resultFolderField.getText().trim().isEmpty()) {
            showError("Validation Error", "Please select a result folder.");
            return false;
        }
        
        File monitorDir = new File(monitorFolderField.getText());
        if (!monitorDir.exists() || !monitorDir.isDirectory()) {
            showError("Validation Error", "Monitor folder does not exist or is not a directory.");
            return false;
        }
        
        File resultDir = new File(resultFolderField.getText());
        if (!resultDir.exists()) {
            resultDir.mkdirs();
        }
        
        return true;
    }

    private void createStopSignalFile(boolean useCount) {
        // Implementation for creating stop signal file
        // This would write a control file to the monitor folder
        consoleOutput.appendText("Stop signal created.\n");
    }

    private void startMonitoringTask() {
        // Create background task for real file monitoring using WatchService
        Task<Void> monitorTask = new Task<Void>() {
            @Override
            protected Void call() throws Exception {
                try {
                    // Set up WatchService to monitor the folder for new files
                    WatchService watchService = FileSystems.getDefault().newWatchService();
                    Path monitorPath = Paths.get(monitorFolderField.getText());
                    monitorPath.register(watchService, StandardWatchEventKinds.ENTRY_CREATE);
                    
                    consoleOutput.appendText("File monitoring started. Watching folder: " + monitorPath + "\n");
                    
                    // Process any existing files in the folder first
                    File[] existingFiles = monitorPath.toFile().listFiles();
                    if (existingFiles != null && existingFiles.length > 0) {
                        consoleOutput.appendText("Found " + existingFiles.length + " existing files, processing...\n");
                        for (File file : existingFiles) {
                            if (isValidRawFile(file.getName())) {
                                processNewFile(file);
                            }
                        }
                    }
                    
                    // Watch for new files
                    while (isMonitoring && !isCancelled()) {
                        WatchKey key;
                        try {
                            key = watchService.take();
                        } catch (InterruptedException e) {
                            Thread.currentThread().interrupt();
                            break;
                        }
                        
                        for (WatchEvent<?> event : key.pollEvents()) {
                            WatchEvent.Kind<?> kind = event.kind();
                            if (kind == StandardWatchEventKinds.OVERFLOW) {
                                continue;
                            }
                            
                            @SuppressWarnings("unchecked")
                            WatchEvent<Path> ev = (WatchEvent<Path>) event;
                            Path filename = ev.context();
                            Path filePath = monitorPath.resolve(filename);
                            
                            if (Files.isRegularFile(filePath) && isValidRawFile(filename.toString())) {
                                consoleOutput.appendText("New file detected: " + filename + "\n");
                                
                                // Wait for file to be completely written before processing
                                waitForFileCompletionAndProcess(filePath);
                            }
                        }
                        
                        boolean valid = key.reset();
                        if (!valid) {
                            break;
                        }
                        
                        // Check stop conditions
                        if (targetFileCount > 0 && processedFileCount >= targetFileCount) {
                            Platform.runLater(() -> {
                                consoleOutput.appendText("Target file count reached. Stopping monitoring.\n");
                                resetMonitoring();
                            });
                            break;
                        }
                    }
                    
                    watchService.close();
                    consoleOutput.appendText("File monitoring stopped.\n");
                    
                } catch (IOException e) {
                    consoleOutput.appendText("Error in file monitoring: " + e.getMessage() + "\n");
                    e.printStackTrace();
                }
                
                return null;
            }
        };
        
        new Thread(monitorTask).start();
    }
    
    /**
     * Check if a file is a valid raw file for processing
     */
    private boolean isValidRawFile(String fileName) {
        String lowerName = fileName.toLowerCase();
        return lowerName.endsWith(".raw") || lowerName.endsWith(".d") || 
               lowerName.endsWith(".dia") || lowerName.endsWith(".d");
    }
    
    /**
     * Process a new file that was detected
     */
    private void processNewFile(File file) {
        String fileName = file.getName();
        String baseFileName = getBaseFileName(fileName);
        
        Platform.runLater(() -> {
            monitorFilesList.getItems().add(baseFileName);
            ProcessingRow newRow = new ProcessingRow(baseFileName);
            resultData.add(newRow);
            consoleOutput.appendText("Processing file: " + baseFileName + "\n");
            processFileInBackground(newRow);
        });
        
        processedFileCount++;
    }
    
    /**
     * Get the base file name without extension
     */
    private String getBaseFileName(String fileName) {
        String lowerName = fileName.toLowerCase();
        if (lowerName.endsWith(".raw")) {
            return fileName.substring(0, fileName.length() - 4);
        } else if (lowerName.endsWith(".d")) {
            return fileName.substring(0, fileName.length() - 2);
        } else if (lowerName.endsWith(".dia")) {
            if (lowerName.endsWith(".raw.dia")) {
                return fileName.substring(0, fileName.length() - 8);
            } else if (lowerName.endsWith(".raw.d")) {
                return fileName.substring(0, fileName.length() - 6);
            } else {
                return fileName.substring(0, fileName.length() - 4);
            }
        }
        return fileName;
    }
    
    /**
     * Wait for file to be completely written before processing
     * Based on MetaDiaMonitorTask.waitForFileCompletionAndProcess()
     */
    private void waitForFileCompletionAndProcess(Path filePath) {
        try {
            long previousSize = -1;
            while (true) {
                // Check the file size at regular intervals
                long currentSize = Files.size(filePath);
                if (currentSize == previousSize) {
                    // Try to acquire an exclusive lock to see if any other process is still writing
                    try (FileChannel fileChannel = FileChannel.open(filePath, StandardOpenOption.WRITE)) {
                        FileLock lock = fileChannel.tryLock();
                        if (lock != null) {
                            lock.release(); // Release the lock immediately
                            break; // File is not in use
                        }
                    } catch (Exception e) {
                        // The file is still in use or another process is writing to it
                        consoleOutput.appendText("File " + filePath.getFileName() + " is still being written, waiting...\n");
                    }
                }
                previousSize = currentSize;
                Thread.sleep(1000); // Wait for 1 second before checking again
            }
            
            // File is ready for processing
            processNewFile(filePath.toFile());
            
        } catch (IOException | InterruptedException e) {
            consoleOutput.appendText("Error waiting for file completion: " + e.getMessage() + "\n");
            e.printStackTrace();
            Thread.currentThread().interrupt();
        }
    }

    private void processFileInBackground(ProcessingRow row) {
        // Mock processing of a file
        Task<Void> processTask = new Task<Void>() {
            @Override
            protected Void call() throws Exception {
                // Round 1
                Platform.runLater(() -> {
                    row.setRound1("Running");
                    consoleOutput.appendText("Starting Round 1 for " + row.getFileName() + "\n");
                });
                Thread.sleep(3000);
                Platform.runLater(() -> row.setRound1("Done"));
                
                // Model
                Platform.runLater(() -> {
                    row.setModel("Running");
                    consoleOutput.appendText("Running Model for " + row.getFileName() + "\n");
                });
                Thread.sleep(2000);
                Platform.runLater(() -> row.setModel("Done"));
                
                // Round 2
                Platform.runLater(() -> {
                    row.setRound2("Running");
                    consoleOutput.appendText("Starting Round 2 for " + row.getFileName() + "\n");
                });
                Thread.sleep(2500);
                Platform.runLater(() -> row.setRound2("Done"));
                
                // Round 3
                Platform.runLater(() -> {
                    row.setRound3("Running");
                    consoleOutput.appendText("Starting Round 3 for " + row.getFileName() + "\n");
                });
                Thread.sleep(2000);
                Platform.runLater(() -> row.setRound3("Done"));
                
                // Round 4
                Platform.runLater(() -> {
                    row.setRound4("Running");
                    consoleOutput.appendText("Starting Round 4 for " + row.getFileName() + "\n");
                });
                Thread.sleep(1500);
                Platform.runLater(() -> {
                    row.setRound4("Done");
                    consoleOutput.appendText("Processing complete for " + row.getFileName() + "\n");
                });
                
                return null;
            }
        };
        
        new Thread(processTask).start();
    }

    private VBox createThreadCountSection() {
        VBox section = new VBox(10);
        section.setPadding(new Insets(10));
        section.setStyle("-fx-background-color: #f0f8ff; -fx-background-radius: 5;");

        Label titleLabel = createSectionTitle("Performance Settings");

        GridPane settingsGrid = new GridPane();
        settingsGrid.setHgap(20);
        settingsGrid.setVgap(10);
        settingsGrid.setAlignment(Pos.CENTER_LEFT);

        // Threads selection only (database settings are in I/O tab)
        Label threadsLabel = new Label("Number of Threads:");
        threadsLabel.setStyle("-fx-font-weight: bold;");
        settingsGrid.add(threadsLabel, 0, 0);

        threadsComboBox = new ComboBox<>();
        for (int i = 1; i <= Runtime.getRuntime().availableProcessors(); i++) {
            threadsComboBox.getItems().add(i);
        }
        threadsComboBox.setValue(Math.min(4, Runtime.getRuntime().availableProcessors()));
        threadsComboBox.setMinWidth(100);
        threadsComboBox.setPrefWidth(120);
        threadsComboBox.setOnAction(e -> updateParameterView());
        settingsGrid.add(threadsComboBox, 1, 0);

        section.getChildren().addAll(titleLabel, settingsGrid);
        return section;
    }

    private void updateParameterView() {
        if (parameterView == null) {
            return;
        }
        
        StringBuilder sb = new StringBuilder();
        sb.append("=== DIA Real-time Analysis Parameters ===\n\n");

        // Real-time monitoring settings (no input/output needed)
        sb.append("Real-time Monitoring Settings:\n");
        sb.append("  Monitor Folder: ").append(monitorFolderField != null ? monitorFolderField.getText() : "Not set").append("\n");
        sb.append("  Result Folder: ").append(resultFolderField != null ? resultFolderField.getText() : "Not set").append("\n");
        sb.append("\n");

        // Database settings
        sb.append("Database Settings:\n");
        if (magDatabasePanel != null) {
            MagDbItem selectedDb = magDatabasePanel.getSelectedMagDatabase();
            if (selectedDb != null) {
                sb.append("  MAG Database: ").append(selectedDb.getCatalogueID()).append("\n");
            } else {
                sb.append("  MAG Database: Not selected\n");
            }
            MagVersion selectedVersion = magDatabasePanel.getSelectedMagVersion();
            if (selectedVersion != null) {
                sb.append("  MAG Version: ").append(selectedVersion.getCatalogueVersion()).append("\n");
            }
            sb.append("  Append Host DB: ").append(magDatabasePanel.isAppendHostDb()).append("\n");
            if (magDatabasePanel.isAppendHostDb() && magDatabasePanel.getHostDatabase() != null && !magDatabasePanel.getHostDatabase().isEmpty()) {
                sb.append("  Host Database: ").append(magDatabasePanel.getHostDatabase()).append("\n");
            }
        } else {
            sb.append("  MAG Database: Database panel not available\n");
        }
        sb.append("\n");

        // Monitoring settings
        sb.append("Monitoring Settings:\n");
        sb.append("  Monitor Folder: ").append(monitorFolderField != null ? monitorFolderField.getText() : "").append("\n");
        sb.append("  Result Folder: ").append(resultFolderField != null ? resultFolderField.getText() : "").append("\n");
        sb.append("  Perform Analysis: ").append(performAnalysisCheckBox != null ? performAnalysisCheckBox.isSelected() : false).append("\n");
        sb.append("\n");

        // Performance settings
        sb.append("Performance Settings:\n");
        if (threadsComboBox != null) {
            sb.append("  Thread Count: ").append(threadsComboBox.getValue()).append("\n");
        }
        sb.append("\n");

        // Library configuration
        sb.append("Spectral Libraries:\n");
        if (libraryGenerationPanel != null) {
            String[] selectedLibs = libraryGenerationPanel.getSelectedLibraryNames();
            if (selectedLibs.length > 0) {
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

        // Model configuration
        sb.append("ML Models (Optional):\n");
        if (modelGenerationPanel != null) {
            String[] selectedModels = modelGenerationPanel.getSelectedModelNames();
            if (selectedModels.length > 0) {
                for (String modelName : selectedModels) {
                    sb.append("  ").append(modelName).append("\n");
                }
            } else {
                sb.append("  No models selected (models are optional)\n");
            }
        } else {
            sb.append("  Model panel not available\n");
        }
        sb.append("\n");

        parameterView.setText(sb.toString());
    }

    private Label createSectionTitle(String title) {
        Label label = new Label(title);
        label.setFont(Font.font("System", FontWeight.BOLD, 14));
        label.setTextFill(Color.DARKBLUE);
        return label;
    }

    private void showError(String title, String content) {
        Alert alert = new Alert(Alert.AlertType.ERROR);
        alert.setTitle(title);
        alert.setHeaderText(null);
        alert.setContentText(content);
        alert.showAndWait();
    }

    private void updatePanelsWithMagDb(MagDbItem magDbItem) {
        if (libraryGenerationPanel != null) {
            libraryGenerationPanel.setMagDbItem(magDbItem);
        }
        if (modelGenerationPanel != null) {
            modelGenerationPanel.setMagDbItem(magDbItem);
        }
    }

    /**
     * Collect parameters from UI for real-time monitoring
     */
    private MetaParameterDia collectParametersFromUI() {
        try {
            // For real-time monitoring, we may not have raw files yet
            // So we create minimal parameters for configuration
            
            // Create minimal MetaData (empty experiments for real-time monitoring)
            MetaData metaData = new MetaData(new String[0], new String[0], new int[0], new int[0]);
            
            // Get result directory from monitor settings
            String result = resultFolderField != null ? resultFolderField.getText() : "";
            
            // Get MAG database information
            MagDbItem selectedMagDb = null;
            if (magDatabasePanel != null) {
                selectedMagDb = magDatabasePanel.getSelectedMagDatabase();
            }
            if (selectedMagDb == null) {
                System.err.println("No MAG database selected");
                return null;
            }
            
            MagVersion selectedVersion = null;
            if (magDatabasePanel != null) {
                selectedVersion = magDatabasePanel.getSelectedMagVersion();
            }
            if (selectedVersion == null) {
                System.err.println("No MAG database version selected");
                return null;
            }
            
            // Set the version on the MagDbItem
            selectedMagDb.setUsedVersion(selectedVersion);
            
            // Get host database settings
            boolean appendHostDb = magDatabasePanel != null ? magDatabasePanel.isAppendHostDb() : false;
            String hostDB = magDatabasePanel != null ? magDatabasePanel.getHostDatabase() : "";
            
            // Get thread count
            int threadCount = threadsComboBox != null ? threadsComboBox.getValue() : 4;
            
            // Create MetaBaseParameter
            MetaBaseParameter metaBasePar = new MetaBaseParameter(
                uk.ac.quadram.metapilot.task.MetaLabWorkflowType.DiaNNMAGRT.name(), // Workflow type
                metaData,
                result,
                selectedMagDb.getCurrentFile() != null ? selectedMagDb.getCurrentFile().getAbsolutePath() : "", // MicroDb path
                appendHostDb,
                hostDB,
                "DIA", // MS2 scan mode (default for DIA)
                1, // Core count
                threadCount,
                true // Meta workflow
            );
            
            // Get library information
            String[] library = new String[0];
            if (libraryGenerationPanel != null) {
                String[] selectedLibs = libraryGenerationPanel.getSelectedLibraryNames();
                if (selectedLibs != null && selectedLibs.length > 0) {
                    library = selectedLibs;
                }
            }
            
            // Get model information
            String[] models = new String[0];
            if (modelGenerationPanel != null) {
                String[] selectedModels = modelGenerationPanel.getSelectedModelNames();
                if (selectedModels != null && selectedModels.length > 0) {
                    models = selectedModels;
                }
            }
            
            // For real-time monitoring, assume library search mode (first mode)
            boolean librarySearch = true;
            boolean selfModelSearch = false;
            boolean useModel = models.length > 0;
            
            // Create MetaParameterDia
            MetaParameterDia metaParameterDia = new MetaParameterDia(
                metaBasePar,
                selectedMagDb,
                librarySearch,
                selfModelSearch,
                library,
                useModel,
                models
            );
            
            return metaParameterDia;
            
        } catch (Exception e) {
            System.err.println("Error generating MetaParameterDia: " + e.getMessage());
            e.printStackTrace();
            return null;
        }
    }

    /**
     * Load last used parameters and apply them to the UI
     */
    private void loadLastUsedParameters() {
        System.out.println("DEBUG: loadLastUsedParameters() called for DIA Real-time panel");
        try {
            MetaParameterDia params = (MetaParameterDia) parameterIO.loadLastUsedParameters();
            if (params != null) {
                System.out.println("DEBUG: Successfully loaded parameters from file");
                applyParametersToUI(params);
                System.out.println("DEBUG: Applied parameters to UI successfully");
            } else {
                System.out.println("DEBUG: No parameters loaded from file, using defaults");
            }
        } catch (Exception e) {
            System.err.println("WARNING: Failed to load last used DIA parameters: " + e.getMessage());
            e.printStackTrace();
            // Continue with default UI state
        }
    }
    
    /**
     * Apply parameters to the UI components
     */
    private void applyParametersToUI(MetaParameterDia params) {
        if (params == null) {
            System.out.println("DEBUG: applyParametersToUI called with null params");
            return;
        }
        
        System.out.println("DEBUG: applyParametersToUI called with valid params");
        
        // Apply result directory to monitor settings if available
        if (params.getResult() != null && !params.getResult().trim().isEmpty()) {
            if (resultFolderField != null) {
                resultFolderField.setText(params.getResult());
                System.out.println("DEBUG: Applied result directory: " + params.getResult());
            }
        }
        
        // Apply thread count
        if (threadsComboBox != null) {
            threadsComboBox.setValue(params.getThreadCount());
        }
        
        // Apply MAG database settings
        if (params.getMagDbItem() != null && magDatabasePanel != null) {
            MagDbItem magDbItem = params.getMagDbItem();
            // Find and select the matching MAG database
            for (MagDbItem item : magDatabasePanel.getAvailableMagDatabases()) {
                if (item.getCatalogueID().equals(magDbItem.getCatalogueID()) && 
                    item.getIdentifier().equals(magDbItem.getIdentifier())) {
                    magDatabasePanel.setSelectedMagDatabase(item);
                    
                    // Also apply version selection if available
                    if (magDbItem.getUsedVersion() != null) {
                        for (MagVersion version : item.getVersionList()) {
                            if (version.getCatalogueVersion().equals(magDbItem.getUsedVersion().getCatalogueVersion())) {
                                magDatabasePanel.setSelectedMagVersion(version);
                                System.out.println("DEBUG: Applied version selection: " + version.getCatalogueVersion());
                                break;
                            }
                        }
                    }
                    break;
                }
            }
            magDatabasePanel.setAppendHostDb(params.isAppendHostDb());
            if (params.getHostDB() != null) {
                magDatabasePanel.setHostDatabase(params.getHostDB());
            }
        }
        
        // Update parameter view after applying parameters
        Platform.runLater(() -> updateParameterView());
    }
    
    /**
     * Save current parameters to file
     */
    public void saveCurrentParameters() {
        try {
            MetaParameterDia currentParams = collectParametersFromUI();
            if (currentParams != null) {
                boolean saved = parameterIO.saveLastUsedParameters(currentParams);
                if (saved) {
                    System.out.println("DEBUG: Successfully saved current parameters to last_used_dia_parameters.json");
                } else {
                    System.err.println("WARNING: Failed to save current parameters");
                }
            } else {
                System.err.println("WARNING: No parameters to save - validation failed");
            }
        } catch (Exception e) {
            System.err.println("ERROR: Failed to save current parameters: " + e.getMessage());
            e.printStackTrace();
        }
    }
} 