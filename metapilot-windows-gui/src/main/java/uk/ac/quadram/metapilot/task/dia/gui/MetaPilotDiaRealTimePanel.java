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
// import javafx.beans.property.SimpleStringProperty; // Not used
// import javafx.beans.property.StringProperty; // Not used
import javafx.concurrent.Task;
import uk.ac.quadram.metapilot.task.common.gui.MagDatabaseSettingsPanel;
import uk.ac.quadram.metapilot.task.common.gui.ConsolePanel;
import uk.ac.quadram.metapilot.task.par.MetaData;
import uk.ac.quadram.metapilot.task.par.MetaBaseParameter;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;
import uk.ac.quadram.metapilot.task.mag.MagDbItem.MagVersion;
import uk.ac.quadram.metapilot.task.dia.par.DiaRTParameterIO;
import uk.ac.quadram.metapilot.task.dia.par.MetaParameterDiaRT;
import uk.ac.quadram.metapilot.task.dia.MetaDiaCombineTask;
import uk.ac.quadram.metapilot.task.dia.MetaDiaMonitorTask;

import uk.ac.quadram.metapilot.task.dia.TableRowData;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.io.File;

public class MetaPilotDiaRealTimePanel extends VBox {
    // UI Components - Database
    private MagDatabaseSettingsPanel magDatabasePanel;
    private ConsolePanel consolePanel;
    private TextArea consoleOutput;

    // Real-time monitoring components
    private TextField monitorFolderField;
    private Button browseMonitorFolderButton;
    private TextField resultFolderField;
    private Button browseResultFolderButton;
    private ListView<String> monitorFilesList;
    private TableView<TableRowData> resultTable;
    
    // Control panel components
    private Button startMonitorButton;
    private Button stopMonitorButton;
    private Button stopAfterCountButton;
    private Button stopCombineButton;
    private Button resetButton;
    private ComboBox<Integer> fileCountComboBox;
    private CheckBox performAnalysisCheckBox;

    // DIA settings
    private ComboBox<Integer> threadsComboBox;
    /** Real-time workflow mode: 3 marker + 4 seed modes. */
    private ComboBox<String> rtWorkflowModeComboBox;

    // DIA-NN modification checkboxes (six standard modifications)
    private CheckBox mExcisionCheckBox;
    private CheckBox carbamiCCheckBox;
    private CheckBox oxidationMCheckBox;
    private CheckBox acetyProNCheckBox;
    private CheckBox phosphoSTYCheckBox;
    private CheckBox kGGCheckBox;

    // Separate panels for library and model generation
    private DiaLibraryGenerationPanel libraryGenerationPanel;
    private DiaModelGenerationPanel modelGenerationPanel;


    // Parameter I/O
    private MetaParameterDiaRT currentAnalysisParameters;
    
    // Real-time task management
    private MetaDiaMonitorTask monitorTask;
    private MetaDiaCombineTask combineTask;
    private ObservableList<String> javafxRawListItems;
    private ObservableList<TableRowData> javafxResultTableItems;
    private HashMap<String, String[]> rawResultFileMap;
    
    // Resource parameters
    private MetaResourceParameters resourceParams;

    // Monitoring status
    private boolean isMonitoring = false;
    private int targetFileCount = 0;
    private boolean performAnalysisAfterComplete = false;
    private DiaRTParameterIO parameterIO;

    public MetaPilotDiaRealTimePanel() {
        // Initialize parameter IO instance
        parameterIO = new DiaRTParameterIO();
        
        // Initialize JavaFX models
        javafxRawListItems = FXCollections.observableArrayList();
        javafxResultTableItems = FXCollections.observableArrayList();
        
        // Initialize console panel
        consolePanel = new ConsolePanel();
        // Don't set start/stop callbacks for real-time workflow - use custom buttons instead
        consolePanel.setViewReportCallback(this::viewReport);
        // Hide the Start Analysis button since we use Start Monitor instead
        consolePanel.setButtonsEnabled(false, false, true); // Disable start and stop, enable view report
        // Explicitly hide progress indicator since we control it manually in real-time workflow
        consolePanel.hideProgressIndicator();
        
        initializeUI();
    }

    public MetaPilotDiaRealTimePanel(MagDbItem[] magDbItems, MetaResourceParameters resourceParams) {
        // Initialize parameter IO instance
        parameterIO = new DiaRTParameterIO();
        this.resourceParams = resourceParams;
        
        // Initialize JavaFX models
        javafxRawListItems = FXCollections.observableArrayList();
        javafxResultTableItems = FXCollections.observableArrayList();
        
        // Initialize console panel
        consolePanel = new ConsolePanel();
        // Don't set start/stop callbacks for real-time workflow - use custom buttons instead
        consolePanel.setViewReportCallback(this::viewReport);
        // Hide the Start Analysis button since we use Start Monitor instead
        consolePanel.setButtonsEnabled(false, false, true); // Disable start and stop, enable view report
        // Explicitly hide progress indicator since we control it manually in real-time workflow
        consolePanel.hideProgressIndicator();
        
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
            // Sync library table to the current combobox selection (after selection is established)
            Platform.runLater(() -> {
                MagDbItem selected = magDatabasePanel.getSelectedMagDatabase();
                if (selected != null) {
                    updatePanelsWithMagDb(selected);
                } else if (magDbItems != null && magDbItems.length > 0) {
                    updatePanelsWithMagDb(magDbItems[0]);
                }
            });
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
        magDatabasePanel = new MagDatabaseSettingsPanel(false);
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
                    MetaParameterDiaRT currentParams = collectParametersFromUIForRealTime();
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

    /**
     * Real-time workflow selection section.
     * Provides unified choices for marker and seed workflows that will be applied
     * to each raw file during monitoring.
     */
    private VBox createWorkflowSelectionSection() {
        VBox section = new VBox(10);
        section.setPadding(new Insets(10));
        section.setStyle("-fx-background-color: #f0f8ff; -fx-background-radius: 5;");

        Label titleLabel = createSectionTitle("Real-time Workflow");

        Label workflowLabel = new Label("Workflow mode:");
        workflowLabel.setStyle("-fx-font-weight: bold;");

        rtWorkflowModeComboBox = new ComboBox<>(
                FXCollections.observableArrayList(
                        "Marker - Focus",
                        "Marker - Comprehensive",
                        "Marker - Deep",
                        "Seed - Direct search",
                        "Seed - Focus search",
                        "Seed - Comprehensive search",
                        "Seed - Deep comprehensive search"));
        rtWorkflowModeComboBox.setMinWidth(260);
        rtWorkflowModeComboBox.setPrefWidth(320);
        // Default to a simple seed-direct mode
        rtWorkflowModeComboBox.setValue("Seed - Direct search");
        rtWorkflowModeComboBox.setOnAction(e -> updateParameterView());

        Tooltip workflowTooltip = new Tooltip(
                "Select which DIA workflow to apply to each new raw file:\n"
                + " - Marker modes use Marker4-based genome selection (Focus / Comprehensive / Deep).\n"
                + " - Seed modes use seed library workflows (Direct / Focus / Comprehensive / Deep comprehensive).\n"
                + "Per-file searches run in no-quant mode; a final combined quant search is run after monitoring stops.");
        rtWorkflowModeComboBox.setTooltip(workflowTooltip);

        VBox inner = new VBox(5);
        inner.getChildren().addAll(workflowLabel, rtWorkflowModeComboBox);

        section.getChildren().addAll(titleLabel, inner);
        return section;
    }

    // Create IO tab (same as DIA Normal Panel)

    // ProcessingRow class removed - now using TableRowData from MetaDiaMonitorTask

    private VBox createParametersTab() {
        VBox tabContent = new VBox(10);
        tabContent.setPadding(new Insets(10));
        tabContent.setMinWidth(960);
        tabContent.setMinHeight(760);

        // 2-column grid like seed panel (no metadata)
        GridPane mainGrid = new GridPane();
        mainGrid.setHgap(15);
        mainGrid.setVgap(15);
        mainGrid.setPadding(new Insets(10));
        mainGrid.setMinWidth(940);
        mainGrid.setMinHeight(600);
        VBox.setVgrow(mainGrid, Priority.ALWAYS);

        ColumnConstraints col1 = new ColumnConstraints();
        col1.setPercentWidth(50);
        col1.setMinWidth(430);
        col1.setHgrow(Priority.ALWAYS);
        ColumnConstraints col2 = new ColumnConstraints();
        col2.setPercentWidth(50);
        col2.setMinWidth(430);
        col2.setHgrow(Priority.ALWAYS);
        mainGrid.getColumnConstraints().addAll(col1, col2);

        RowConstraints row1 = new RowConstraints();
        row1.setVgrow(Priority.ALWAYS);
        row1.setMinHeight(400);
        mainGrid.getRowConstraints().add(row1);

        // Left column: Database, Workflow & threads, Library controls
        VBox leftColumn = new VBox(15);
        leftColumn.setMinWidth(430);
        leftColumn.getChildren().addAll(
            createSectionTitle("Database Settings"),
            magDatabasePanel,
            createWorkflowAndThreadSection()
        );
        libraryGenerationPanel = new DiaLibraryGenerationPanel(consoleOutput, true);
        libraryGenerationPanel.setParameterUpdateCallback(() -> updateParameterView());
        leftColumn.getChildren().add(libraryGenerationPanel.getControlsSection());
        VBox.setVgrow(libraryGenerationPanel.getControlsSection(), Priority.ALWAYS);
        GridPane.setVgrow(leftColumn, Priority.ALWAYS);

        // Right column: Library table only
        VBox rightColumn = new VBox(10);
        rightColumn.setMinWidth(430);
        rightColumn.getChildren().add(libraryGenerationPanel.getTableSection());
        VBox.setVgrow(libraryGenerationPanel.getTableSection(), Priority.ALWAYS);
        GridPane.setVgrow(rightColumn, Priority.ALWAYS);

        mainGrid.add(leftColumn, 0, 0);
        mainGrid.add(rightColumn, 1, 0);

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
        
        monitorSection.setMinHeight(150);
        monitorSection.setPrefHeight(220);
        VBox.setVgrow(monitorSection, Priority.ALWAYS);
        
        resultSection.setMinHeight(180);
        resultSection.setPrefHeight(220);
        VBox.setVgrow(resultSection, Priority.ALWAYS);
        
        controlSection.setMinHeight(150);
        controlSection.setMaxHeight(180);

        controlsPanel.getChildren().addAll(folderSection, monitorSection, resultSection, controlSection);
        return controlsPanel;
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

        // Create ListView and connect to JavaFX model
        monitorFilesList = new ListView<>();
        monitorFilesList.setItems(javafxRawListItems);
        monitorFilesList.setMinHeight(150);
        monitorFilesList.setPrefHeight(180);
        monitorFilesList.setMaxHeight(300);
        VBox.setVgrow(monitorFilesList, Priority.ALWAYS);

        // Add components only if they're not null
        if (monitorFilesList != null) {
        section.getChildren().addAll(titleLabel, monitorFilesList);
            System.out.println("DEBUG: Successfully added monitorFilesList to VBox");
        } else {
            // Fallback to label if ListView is null
            Label fallbackLabel = new Label("Monitor Files List (Not Available)");
            fallbackLabel.setStyle("-fx-background-color: #f0f0f0; -fx-border-color: #ccc; -fx-padding: 10;");
            section.getChildren().addAll(titleLabel, fallbackLabel);
            System.out.println("DEBUG: monitorFilesList is null, using fallback label");
        }
        return section;
    }

    private VBox createResultSection() {
        VBox section = new VBox(10);
        section.setPadding(new Insets(10));
        section.setStyle("-fx-background-color: #f0fff0; -fx-background-radius: 5;");
        VBox.setVgrow(section, Priority.ALWAYS);

        Label titleLabel = createSectionTitle("Processing Results");

        // Restore TableView for the working sync solution
        resultTable = new TableView<>(javafxResultTableItems);
        resultTable.setMinHeight(150);
        resultTable.setPrefHeight(200);
        resultTable.setMaxHeight(Double.MAX_VALUE);
        resultTable.setColumnResizePolicy(TableView.CONSTRAINED_RESIZE_POLICY_FLEX_LAST_COLUMN);
        VBox.setVgrow(resultTable, Priority.ALWAYS);

        // File Name column
        TableColumn<TableRowData, String> fileCol = new TableColumn<>("File");
        fileCol.setCellValueFactory(new PropertyValueFactory<>("fileName"));
        fileCol.setMinWidth(180);
        fileCol.setPrefWidth(240);

        TableColumn<TableRowData, String> statusCol = new TableColumn<>("Status");
        statusCol.setCellValueFactory(new PropertyValueFactory<>("status"));
        statusCol.setMinWidth(100);
        statusCol.setPrefWidth(120);

        TableColumn<TableRowData, String> detailsCol = new TableColumn<>("Details");
        detailsCol.setCellValueFactory(new PropertyValueFactory<>("details"));
        detailsCol.setMinWidth(260);
        detailsCol.setPrefWidth(420);

        resultTable.getColumns().clear();
        resultTable.getColumns().add(fileCol);
        resultTable.getColumns().add(statusCol);
        resultTable.getColumns().add(detailsCol);

        // Add components only if they're not null
        if (resultTable != null) {
        section.getChildren().addAll(titleLabel, resultTable);
            System.out.println("DEBUG: Successfully added resultTable to VBox");
        } else {
            // Fallback to label if TableView is null
            Label fallbackLabel = new Label("Processing Results Table (Not Available)");
            fallbackLabel.setStyle("-fx-background-color: #f0f0f0; -fx-border-color: #ccc; -fx-padding: 10;");
            section.getChildren().addAll(titleLabel, fallbackLabel);
            System.out.println("DEBUG: resultTable is null, using fallback label");
        }
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

        // Row 3: Stop combine button
        stopCombineButton = new Button("Stop Combined Analysis");
        stopCombineButton.setPrefWidth(180);
        stopCombineButton.setDisable(true);
        stopCombineButton.setOnAction(e -> stopCombinedAnalysis());
        controlGrid.add(stopCombineButton, 0, 2, 2, 1);

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
            
            // If we're in the middle of starting monitoring and user changed folder,
            // automatically restart the monitoring process with the new folder
            if (isMonitoring) {
                consolePanel.appendConsoleText("Monitor folder changed to: " + directory.getAbsolutePath() + "\n");
                // The monitoring will continue with the new folder
            }
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

        // Check if monitor folder contains existing raw files
        File monitorFolder = new File(monitorFolderField.getText());
        List<String> existingRawFiles = findExistingRawFiles(monitorFolder);
        
        if (!existingRawFiles.isEmpty()) {
            // Ask user what to do with existing files
            String message = "The monitor folder contains " + existingRawFiles.size() + " existing raw file(s):\n\n";
            for (int i = 0; i < Math.min(existingRawFiles.size(), 5); i++) {
                message += "• " + existingRawFiles.get(i) + "\n";
            }
            if (existingRawFiles.size() > 5) {
                message += "• ... and " + (existingRawFiles.size() - 5) + " more files\n";
            }
            message += "\nWhat would you like to do?";
            
            Alert choiceAlert = new Alert(Alert.AlertType.CONFIRMATION);
            choiceAlert.setTitle("Existing Files Found");
            choiceAlert.setHeaderText("Monitor folder contains existing raw files");
            choiceAlert.setContentText(message);
            
            ButtonType continueButton = new ButtonType("Continue with existing files");
            ButtonType changeFolderButton = new ButtonType("Change to another folder");
            ButtonType cancelButton = new ButtonType("Cancel", ButtonBar.ButtonData.CANCEL_CLOSE);
            
            choiceAlert.getButtonTypes().setAll(continueButton, changeFolderButton, cancelButton);
            
            java.util.Optional<ButtonType> result = choiceAlert.showAndWait();
            
            if (result.get() == cancelButton) {
                return; // User cancelled
            } else if (result.get() == changeFolderButton) {
                browseMonitorFolder(); // Let user select different folder
                // After user selects new folder, check again for existing files
                File newMonitorFolder = new File(monitorFolderField.getText());
                List<String> newExistingFiles = findExistingRawFiles(newMonitorFolder);
                if (!newExistingFiles.isEmpty()) {
                    consolePanel.appendConsoleText("New folder also contains " + newExistingFiles.size() + " existing raw file(s), will process them first.\n");
                }
                // Continue with the new folder (don't return, let it proceed)
            }
            // If user chose continueButton, proceed with existing files
        }

        isMonitoring = true;
        
        // Control console panel state - like other workflows do with Start Analysis
        consolePanel.setButtonsEnabled(false, true, true); // Disable start, enable stop, enable view report
        consolePanel.startStatusTimer(); // Start the status timer
        
        startMonitorButton.setDisable(true);
        stopMonitorButton.setDisable(false);
        stopAfterCountButton.setDisable(false);
        resetButton.setDisable(false);

        // Clear previous results - with null checks
        if (monitorFilesList != null && monitorFilesList.getItems() != null) {
            monitorFilesList.getItems().clear();
        }
        if (javafxResultTableItems != null) {
            javafxResultTableItems.clear();
        }

        consolePanel.appendConsoleText("Starting real-time monitoring...\n");
        consolePanel.appendConsoleText("Monitor folder: " + monitorFolderField.getText() + "\n");
        consolePanel.appendConsoleText("Result folder: " + resultFolderField.getText() + "\n");
        
        if (!existingRawFiles.isEmpty()) {
            consolePanel.appendConsoleText("Found " + existingRawFiles.size() + " existing raw file(s), will process them first.\n");
        }

        // Store analysis preference
        performAnalysisAfterComplete = performAnalysisCheckBox.isSelected();
        
        // Start real monitoring task
        startRealMonitoringTask();
    }

    private void stopMonitoring() {
        isMonitoring = false;
        
        stopMonitorButton.setDisable(true);
        
        // Control console panel state - like other workflows do with Stop Analysis
        consolePanel.setButtonsEnabled(false, false, true); // Disable both start and stop, enable view report
        consolePanel.stopStatusTimer(); // Stop the status timer
        
        consolePanel.appendConsoleText("Stopping monitoring after current files are processed...\n");
        
        // Force stop the monitor task if it's running
        if (monitorTask != null && monitorTask.isRunning()) {
            consolePanel.appendConsoleText("Force stopping monitoring task...\n");
            monitorTask.forceStop();
            consolePanel.appendConsoleText("Monitoring task stop requested.\n");
        }
        
        // Signal the monitor task to stop
        if (monitorTask != null && currentAnalysisParameters != null) {
            currentAnalysisParameters.setStopNow(true);
            consolePanel.appendConsoleText("Stop signal sent to monitoring task.\n");
        }
    }

    private void stopAfterCount() {
        targetFileCount = fileCountComboBox.getValue();
        
        stopAfterCountButton.setDisable(true);
        
        consolePanel.appendConsoleText("Will stop monitoring after " + targetFileCount + " files are processed...\n");
        
        // Signal the monitor task to stop after count
        if (monitorTask != null && currentAnalysisParameters != null) {
            currentAnalysisParameters.setStopComplete(true);
            currentAnalysisParameters.setTaskCount(targetFileCount);
            consolePanel.appendConsoleText("Stop-after-count signal sent to monitoring task.\n");
        }
    }
    
    private void stopCombinedAnalysis() {
        consolePanel.appendConsoleText("Stopping combined analysis...\n");
        
        // Force stop the combine task if it's running
        if (combineTask != null && combineTask.isRunning()) {
            consolePanel.appendConsoleText("Force stopping combined analysis task...\n");
            combineTask.forceStop();
            consolePanel.appendConsoleText("Combined analysis task stop requested.\n");
        }
    }

    private void resetMonitoring() {
        isMonitoring = false;
        targetFileCount = 0;
        
        startMonitorButton.setDisable(false);
        stopMonitorButton.setDisable(true);
        stopAfterCountButton.setDisable(true);
        stopCombineButton.setDisable(true);
        resetButton.setDisable(true);
        
        fileCountComboBox.setValue(10);
        
        // Reset console panel state - like other workflows do when reset
        consolePanel.setButtonsEnabled(true, false, true); // Enable start, disable stop, enable view report
        consolePanel.stopStatusTimer(); // Stop the status timer
        
        consolePanel.appendConsoleText("Monitoring reset.\n");
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
    
    /**
     * Find existing raw files in the monitor folder
     * Supports .raw, .RAW, .d, .D, .dia, .DIA file extensions
     */
    private List<String> findExistingRawFiles(File monitorFolder) {
        List<String> rawFiles = new ArrayList<>();
        
        if (monitorFolder == null || !monitorFolder.exists() || !monitorFolder.isDirectory()) {
            return rawFiles;
        }
        
        File[] files = monitorFolder.listFiles();
        if (files == null) {
            return rawFiles;
        }
        
        for (File file : files) {
            if (file.isFile()) {
                String fileName = file.getName().toLowerCase();
                if (fileName.endsWith(".raw") || fileName.endsWith(".d") || 
                    fileName.endsWith(".dia") || fileName.endsWith(".mzml")) {
                    rawFiles.add(file.getName());
                }
            }
        }
        
        return rawFiles;
    }

    /**
     * Start combined analysis after monitoring is complete
     */
    private void startCombinedAnalysis() {
        try {
            // Prepare raw result file map from monitoring results
            rawResultFileMap = new HashMap<>();
            
            // Collect results from each processed file
            File resultFolder = new File(resultFolderField.getText());
            if (resultFolder.exists()) {
                File[] subDirs = resultFolder.listFiles(File::isDirectory);
                if (subDirs != null) {
                    for (File subDir : subDirs) {
                        if (!subDir.getName().equals("combined")) {
                            File finalParquet = new File(subDir, subDir.getName() + ".parquet");
                            if (finalParquet.exists()) {
                                String baseName = subDir.getName();
                                String rawPath = resolveRawPathForBaseName(baseName);
                                rawResultFileMap.put(baseName, new String[]{
                                    rawPath,
                                    subDir.getAbsolutePath()
                                });
                            }
                        }
                    }
                }
            }
            
            if (rawResultFileMap.isEmpty()) {
                consolePanel.appendConsoleText("No processed files found for combined analysis.\n");
                resetMonitoring();
                return;
            }
            
            consolePanel.appendConsoleText("Found " + rawResultFileMap.size() + " processed files for combined analysis.\n");
            
            // Create MetaParameterDiaRT for combined analysis
            MetaParameterDiaRT combineParams = collectParametersFromUIForRealTime();
            if (combineParams == null) {
                showError("Parameter Error", "Failed to collect parameters for combined analysis.");
                resetMonitoring();
                return;
            }
            
            // Create and start MetaDiaCombineTask
            combineTask = new MetaDiaCombineTask(
                combineParams,
                resourceParams,
                consolePanel.getOverallProgressBar(),
                consolePanel.getCurrentTaskProgressBar(),
                rawResultFileMap
            );
            
            // Enable stop combine button
            stopCombineButton.setDisable(false);
            
            // Set up task callbacks
            combineTask.setOnSucceeded(e -> {
                Platform.runLater(() -> {
                    consolePanel.appendConsoleText("Combined analysis completed successfully!\n");
                    consolePanel.appendConsoleText("Combined results saved to: " + resultFolder.getAbsolutePath() + "/combined\n");
                    // Enable view report button when combined analysis completes
                    consolePanel.enableViewReportButton();
                    stopCombineButton.setDisable(true);
                    resetMonitoring();
                });
            });
            
            combineTask.setOnFailed(e -> {
                Platform.runLater(() -> {
                    Throwable exception = combineTask.getException();
                    consolePanel.appendConsoleText("Combined analysis failed: " + (exception != null ? exception.getMessage() : "Unknown error") + "\n");
                    showError("Analysis Error", "Combined analysis failed: " + (exception != null ? exception.getMessage() : "Unknown error"));
                    stopCombineButton.setDisable(true);
                    resetMonitoring();
                });
            });
            
            combineTask.setOnCancelled(e -> {
                Platform.runLater(() -> {
                    consolePanel.appendConsoleText("Combined analysis was cancelled.\n");
                    stopCombineButton.setDisable(true);
                    resetMonitoring();
                });
            });
            
            // Set up progress monitoring for combine task
            combineTask.messageProperty().addListener((obs, oldMsg, newMsg) -> {
                if (newMsg != null && !newMsg.isEmpty()) {
                    Platform.runLater(() -> {
                        consolePanel.appendConsoleText("Combine: " + newMsg + "\n");
                        consolePanel.updateCurrentTaskFromMessage(newMsg);
                    });
                }
            });
            
            // Start the combine task properly as a JavaFX Task
            new Thread(combineTask).start();
            
        } catch (Exception e) {
            consolePanel.appendConsoleText("ERROR: Failed to start combined analysis - " + e.getMessage() + "\n");
            showError("Analysis Error", "Failed to start combined analysis: " + e.getMessage());
            resetMonitoring();
            e.printStackTrace();
        }
    }

    private void startRealMonitoringTask() {
        try {
            // 1. Collect parameters for real-time monitoring
            consolePanel.appendConsoleText("Collecting parameters from UI...\n");
            MetaParameterDiaRT rtParams = collectParametersFromUIForRealTime();
            if (rtParams == null) {
                consolePanel.appendConsoleText("ERROR: Failed to collect parameters from UI.\n");
                showError("Parameter Error", "Failed to collect parameters. Please check your settings.");
                resetMonitoring();
                return;
            }
            consolePanel.appendConsoleText("✓ Successfully collected parameters from UI\n");
            
            // Set monitoring-specific parameters
            rtParams.setMonitorFolder(new File(monitorFolderField.getText()));
            rtParams.setResultFolder(new File(resultFolderField.getText()));
            rtParams.setCombineResult(performAnalysisAfterComplete);
            rtParams.setStopNow(false);
            rtParams.setStopComplete(targetFileCount > 0);
            rtParams.setTaskCount(targetFileCount > 0 ? targetFileCount : Integer.MAX_VALUE);
            
            // Store current parameters for use in stop control
            currentAnalysisParameters = rtParams;
            
            // 2. Save parameters to file
            String resultDir = rtParams.getResultFolder().getAbsolutePath();
            if (resultDir != null && !resultDir.trim().isEmpty()) {
                consolePanel.appendConsoleText("Saving parameters...\n");
                
                // Save as last used parameters
                boolean lastParamsSaved = parameterIO.saveLastUsedParameters(rtParams);
                if (lastParamsSaved) {
                    consolePanel.appendConsoleText("✓ Last used parameters saved\n");
                } else {
                    consolePanel.appendConsoleText("WARNING: Failed to save last used parameters\n");
                }
                
                // Save to preferences
                boolean prefsSaved = parameterIO.savePreferences(rtParams);
                if (prefsSaved) {
                    consolePanel.appendConsoleText("✓ Preferences saved\n");
                } else {
                    consolePanel.appendConsoleText("WARNING: Failed to save preferences\n");
                }
            }
            
            // 3. Display monitoring configuration
            consolePanel.appendConsoleText("Setting up monitoring configuration...\n");
            consolePanel.appendConsoleText("Monitor folder: " + rtParams.getMonitorFolder().getAbsolutePath() + "\n");
            consolePanel.appendConsoleText("Result folder: " + rtParams.getResultFolder().getAbsolutePath() + "\n");
            consolePanel.appendConsoleText("Perform analysis after completion: " + performAnalysisAfterComplete + "\n");
            consolePanel.appendConsoleText("✓ Monitoring configuration ready\n");
            
            // 4. Create and start MetaDiaMonitorTask
            consolePanel.appendConsoleText("Creating DIA Real-time monitoring task...\n");
            monitorTask = new MetaDiaMonitorTask(
                rtParams,
                resourceParams,
                consolePanel.getOverallProgressBar(),
                consolePanel.getCurrentTaskProgressBar(),
                rtParams.getMonitorFolder(),
                rtParams.getResultFolder(),
                javafxResultTableItems,
                javafxRawListItems,
                consolePanel
            );
            consolePanel.appendConsoleText("✓ Task created successfully\n");
            
            // 5. Set up task monitoring callbacks
            consolePanel.appendConsoleText("Setting up task callbacks...\n");
            monitorTask.setOnRunning(e -> {
                Platform.runLater(() -> {
                    consolePanel.appendConsoleText("Real-time monitoring task started.\n");
                });
            });
            
            monitorTask.setOnSucceeded(e -> {
                Platform.runLater(() -> {
                    consolePanel.appendConsoleText("Real-time monitoring completed successfully!\n");
                    
                    if (performAnalysisAfterComplete) {
                        consolePanel.appendConsoleText("Starting combined analysis...\n");
                        startCombinedAnalysis();
                    } else {
                        consolePanel.appendConsoleText("Monitoring completed. No combined analysis requested.\n");
                        // Enable view report button when monitoring completes
                        consolePanel.enableViewReportButton();
                        resetMonitoring();
                    }
                });
            });
            
            monitorTask.setOnFailed(e -> {
                Platform.runLater(() -> {
                    Throwable exception = monitorTask.getException();
                    consolePanel.appendConsoleText("Real-time monitoring failed: " + (exception != null ? exception.getMessage() : "Unknown error") + "\n");
                    showError("Monitoring Error", "Real-time monitoring failed: " + (exception != null ? exception.getMessage() : "Unknown error"));
                    resetMonitoring();
                });
            });
            
            monitorTask.setOnCancelled(e -> {
                Platform.runLater(() -> {
                    consolePanel.appendConsoleText("Real-time monitoring was cancelled.\n");
                    resetMonitoring();
                });
            });
            
            // Set up progress monitoring to capture task messages
            monitorTask.messageProperty().addListener((obs, oldMsg, newMsg) -> {
                if (newMsg != null && !newMsg.isEmpty()) {
                    Platform.runLater(() -> {
                        consolePanel.appendConsoleText("Monitor: " + newMsg + "\n");
                        consolePanel.updateCurrentTaskFromMessage(newMsg);
                    });
                }
            });
            
            consolePanel.appendConsoleText("✓ Task callbacks configured\n");
            
            // Completely disable all sync tasks to isolate threading issues
            // startProgressSyncTask();
            // Sync mechanism removed - using direct JavaFX model updates
            
            // Manual refresh button removed - UI updates automatically now
            
            // Add a debug button to show Swing model contents
            addDebugButton();
            
            // 6. Start the monitoring task
            consolePanel.appendConsoleText("Starting task execution...\n");
            new Thread(monitorTask).start();
            consolePanel.appendConsoleText("✓ Monitoring task started successfully.\n");
            
            // Enable the working safer sync task
            // startSafeProgressSyncTask(); // Removed - using direct JavaFX model updates
            
        } catch (Exception e) {
            consolePanel.appendConsoleText("ERROR: Failed to start real-time monitoring - " + e.getMessage() + "\n");
            showError("Start Error", "Failed to start real-time monitoring: " + e.getMessage());
            resetMonitoring();
            e.printStackTrace();
        }
    }

    /**
     * Resolve the original monitored file for a base name across supported formats.
     */
    private String resolveRawPathForBaseName(String baseName) {
        if (monitorFolderField == null || monitorFolderField.getText() == null) {
            return baseName;
        }
        File monitorDir = new File(monitorFolderField.getText());
        File[] candidates = monitorDir.listFiles();
        if (candidates != null) {
            for (File candidate : candidates) {
                String name = candidate.getName();
                String lower = name.toLowerCase();
                if (!(lower.endsWith(".raw") || lower.endsWith(".d") || lower.endsWith(".dia") || lower.endsWith(".mzml"))) {
                    continue;
                }
                String normalized = name;
                if (lower.endsWith(".raw")) {
                    normalized = name.substring(0, name.length() - 4);
                } else if (lower.endsWith(".d")) {
                    normalized = name.substring(0, name.length() - 2);
                } else if (lower.endsWith(".dia")) {
                    normalized = name.substring(0, name.length() - 4);
                } else if (lower.endsWith(".mzml")) {
                    normalized = name.substring(0, name.length() - 5);
                }
                if (baseName.equalsIgnoreCase(normalized)) {
                    return candidate.getAbsolutePath();
                }
            }
        }
        return new File(monitorDir, baseName + ".raw").getAbsolutePath();
    }
    
    

    /**
     * Start progress synchronization between Swing components (used by MetaDiaMonitorTask)
     * and JavaFX UI components (this panel)
     */
    private void startProgressSyncTask() {
        Task<Void> syncTask = new Task<Void>() {
            @Override
            protected Void call() throws Exception {
                while (isMonitoring && !isCancelled()) {
                    try {
                        // Platform.runLater(() -> syncProgressFromSwingToJavaFX()); // Removed - using direct JavaFX model updates
                        Thread.sleep(200); // Update every 200ms for more responsive UI
                    } catch (InterruptedException e) {
                        Thread.currentThread().interrupt();
                        break;
                    }
                }
                return null;
            }
        };
        
        new Thread(syncTask).start();
    }
    
    // startSafeProgressSyncTask method removed - using direct JavaFX model updates
    
    // syncProgressFromSwingToJavaFX method removed - using direct JavaFX model updates
    
    // safeSyncProgressFromSwingToJavaFX method removed - using direct JavaFX model updates
    
    // addManualRefreshButton method removed - UI updates automatically now
    
    /**
     * Add a debug button to show Swing model contents
     */
    private void addDebugButton() {
        try {
            Button debugButton = new Button("🐛 Debug Swing Models");
            debugButton.setStyle("-fx-background-color: #FF9800; -fx-text-fill: white; -fx-font-weight: bold;");
            debugButton.setOnAction(e -> {
                try {
                    consolePanel.appendConsoleText("=== DEBUG: Swing Model Contents ===\n");
                    
                    // Check monitor files list
                    consolePanel.appendConsoleText("Monitor Files List Size: " + monitorFilesList.getItems().size() + "\n");
                    for (int i = 0; i < monitorFilesList.getItems().size(); i++) {
                        consolePanel.appendConsoleText("  [" + i + "] " + monitorFilesList.getItems().get(i) + "\n");
                    }
                    
                    // Check result table
                    consolePanel.appendConsoleText("Result Table Size: " + resultTable.getItems().size() + "\n");
                    for (int i = 0; i < resultTable.getItems().size(); i++) {
                        TableRowData row = resultTable.getItems().get(i);
                        consolePanel.appendConsoleText("  [" + i + "] " + row.getFileName() + " - Status: " + row.getStatus() + "\n");
                    }
                    
                    // Check Swing models from monitor task
                    if (monitorTask != null) {
                        consolePanel.appendConsoleText("Monitor Task: " + (monitorTask.isRunning() ? "RUNNING" : "STOPPED") + "\n");
                        // Try to access the Swing models if possible
                        try {
                            // Note: We can't directly access monitorTask's Swing models from here
                            // They are encapsulated within the task for thread safety
                            consolePanel.appendConsoleText("Monitor Task Swing Models: Encapsulated (cannot access directly)\n");
                            consolePanel.appendConsoleText("Note: Swing models are updated correctly (see console output)\n");
                            consolePanel.appendConsoleText("Issue: JavaFX sync is disabled to prevent crashes\n");
                        } catch (Exception ex) {
                            consolePanel.appendConsoleText("Monitor Task Swing Models: Error (" + ex.getMessage() + ")\n");
                        }
                    } else {
                        consolePanel.appendConsoleText("Monitor Task: NULL\n");
                    }
                    
                    consolePanel.appendConsoleText("=== END DEBUG ===\n");
                } catch (Exception ex) {
                    consolePanel.appendConsoleText("DEBUG: Error in debug button: " + ex.getMessage() + "\n");
                }
            });
            
            consolePanel.appendConsoleText("DEBUG: Debug button created\n");
        } catch (Exception e) {
            consolePanel.appendConsoleText("DEBUG: Error creating debug button: " + e.getMessage() + "\n");
        }
    }
    
    /**
     * Schedule a one-time sync after a delay to show initial data
     */
    private void scheduleOneTimeSync() {
        try {
            consolePanel.appendConsoleText("DEBUG: Scheduling one-time sync in 3 seconds...\n");
            
            // Use a simple timer for one-time sync
            new Thread(() -> {
                try {
                    Thread.sleep(3000); // Wait 3 seconds
                    
                    Platform.runLater(() -> {
                        try {
                            consolePanel.appendConsoleText("DEBUG: Performing one-time sync...\n");
                            // syncProgressFromSwingToJavaFX(); // Removed - using direct JavaFX model updates
                            consolePanel.appendConsoleText("DEBUG: One-time sync completed\n");
                        } catch (Exception e) {
                            consolePanel.appendConsoleText("DEBUG: One-time sync error: " + e.getMessage() + "\n");
                        }
                    });
                } catch (InterruptedException e) {
                    consolePanel.appendConsoleText("DEBUG: One-time sync interrupted\n");
                }
            }).start();
            
        } catch (Exception e) {
            consolePanel.appendConsoleText("DEBUG: Error scheduling one-time sync: " + e.getMessage() + "\n");
        }
    }
    
    /**
     * Safely append text to console to prevent TextArea IndexOutOfBoundsException
     */
    private void safeAppendConsoleText(String text) {
        try {
            if (consolePanel != null && text != null) {
                consolePanel.appendConsoleText(text);
            }
        } catch (Exception e) {
            // If console panel fails, just log to system out
            System.out.println("Console panel error: " + e.getMessage());
            System.out.println("Original message: " + text);
        }
    }

    private VBox createWorkflowAndThreadSection() {
        VBox section = new VBox(10);
        section.setPadding(new Insets(10));
        section.setStyle("-fx-background-color: #f0f8ff; -fx-background-radius: 5;");

        Label titleLabel = createSectionTitle("Workflow, Performance & Modifications");

        GridPane settingsGrid = new GridPane();
        settingsGrid.setHgap(20);
        settingsGrid.setVgap(10);
        settingsGrid.setAlignment(Pos.CENTER_LEFT);

        ColumnConstraints col1 = new ColumnConstraints();
        col1.setMinWidth(260);
        col1.setHgrow(Priority.ALWAYS);
        ColumnConstraints col2 = new ColumnConstraints();
        col2.setMinWidth(200);
        col2.setHgrow(Priority.ALWAYS);
        settingsGrid.getColumnConstraints().addAll(col1, col2);

        // Workflow selection: 3 marker + 4 seed modes
        Label workflowLabel = new Label("Workflow mode:");
        workflowLabel.setStyle("-fx-font-weight: bold;");
        rtWorkflowModeComboBox = new ComboBox<>(FXCollections.observableArrayList(
                "Marker – Focus",
                "Marker – Comprehensive",
                "Marker – Deep",
                "Seed – Direct",
                "Seed – Focus",
                "Seed – Comprehensive",
                "Seed – Deep"));
        rtWorkflowModeComboBox.setMinWidth(220);
        rtWorkflowModeComboBox.setPrefWidth(260);
        rtWorkflowModeComboBox.setMaxWidth(Double.MAX_VALUE);
        rtWorkflowModeComboBox.setValue("Seed – Direct");
        rtWorkflowModeComboBox.setOnAction(e -> updateParameterView());

        VBox workflowBox = new VBox(5);
        workflowBox.getChildren().addAll(workflowLabel, rtWorkflowModeComboBox);
        settingsGrid.add(workflowBox, 0, 0);

        // Threads selection (database settings are in Database Settings section)
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

        VBox threadsBox = new VBox(5);
        threadsBox.getChildren().addAll(threadsLabel, threadsComboBox);
        settingsGrid.add(threadsBox, 1, 0);

        // Modification row: all six DIA-NN modifications in one line
        HBox modsBox = new HBox(10);
        modsBox.setAlignment(Pos.CENTER_LEFT);

        Label modsLabel = new Label("Modifications:");
        modsLabel.setStyle("-fx-font-weight: bold;");

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

        modsBox.getChildren().addAll(
                modsLabel,
                mExcisionCheckBox,
                carbamiCCheckBox,
                oxidationMCheckBox,
                acetyProNCheckBox,
                phosphoSTYCheckBox,
                kGGCheckBox);

        // Place modification row spanning both columns underneath workflow/threads
        settingsGrid.add(modsBox, 0, 1, 2, 1);

        section.getChildren().addAll(titleLabel, settingsGrid);
        return section;
    }

    private void updateParameterView() {
        if (consolePanel == null) {
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

        // Performance and workflow settings
        sb.append("Performance Settings:\n");
        if (threadsComboBox != null) {
            sb.append("  Thread Count: ").append(threadsComboBox.getValue()).append("\n");
        }
        if (rtWorkflowModeComboBox != null && rtWorkflowModeComboBox.getValue() != null) {
            sb.append("  Workflow Mode: ").append(rtWorkflowModeComboBox.getValue()).append("\n");
        }
        sb.append("\n");

        // DIA-NN modification settings
        sb.append("Modification Settings:\n");
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

        // Update the console panel's parameter view
        if (consolePanel != null) {
            consolePanel.setParameterText(sb.toString());
        }
    }

    private Label createSectionTitle(String title) {
        Label label = new Label(title);
        // Use default font instead of "System" to avoid PrismTextLayout issues
        label.setFont(Font.font(null, FontWeight.BOLD, 14));
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
    private MetaParameterDiaRT collectParametersFromUIForRealTime() {
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
					uk.ac.quadram.metapilot.workflow.WorkflowType.DIA_REALTIME.getAnalysisType(), // Workflow type
					metaData, result,
					selectedMagDb.getCurrentFile() != null ? selectedMagDb.getCurrentFile().getAbsolutePath() : "", // MicroDb
																													// path
					appendHostDb, hostDB, "DIA", // MS2 scan mode (default for DIA)
					1, // Core count
					threadCount, true // Meta workflow
			);
            
            // Get library information (used by seed modes)
            String[] library = new String[0];
            if (libraryGenerationPanel != null) {
                String[] selectedLibs = libraryGenerationPanel.getSelectedLibraryNames();
                if (selectedLibs != null && selectedLibs.length > 0) {
                    library = selectedLibs;
                }
            }
            
            // Workflow selection: map UI label to internal mode code and DIA flags
            String uiMode = (rtWorkflowModeComboBox != null && rtWorkflowModeComboBox.getValue() != null)
                    ? rtWorkflowModeComboBox.getValue()
                    : "Seed \u2013 Direct";
            String rtModeCode;
            boolean librarySearch;
            boolean selfModelSearch = false;
            boolean useModel = false;
            String[] models = new String[0];

            switch (uiMode) {
                case "Marker \u2013 Focus":
                case "Marker - Focus":
                    rtModeCode = "MARKER_FOCUS";
                    librarySearch = false;
                    break;
                case "Marker \u2013 Comprehensive":
                case "Marker - Comprehensive":
                    rtModeCode = "MARKER_COMPREHENSIVE";
                    librarySearch = false;
                    break;
                case "Marker \u2013 Deep":
                case "Marker - Deep":
                    rtModeCode = "MARKER_DEEP";
                    librarySearch = false;
                    break;
                case "Seed \u2013 Focus":
                case "Seed - Focus search":
                    rtModeCode = "SEED_FOCUS";
                    librarySearch = true;
                    break;
                case "Seed \u2013 Comprehensive":
                case "Seed - Comprehensive search":
                    rtModeCode = "SEED_COMPREHENSIVE";
                    librarySearch = true;
                    break;
                case "Seed \u2013 Deep":
                case "Seed - Deep comprehensive search":
                    rtModeCode = "SEED_DEEP";
                    librarySearch = true;
                    break;
                case "Seed \u2013 Direct":
                case "Seed - Direct search":
                default:
                    rtModeCode = "SEED_DIRECT";
                    librarySearch = true;
                    break;
            }
            
            // Create MetaParameterDiaRT (extends MetaParameterDia)
            MetaParameterDiaRT metaParameterDiaRT = new MetaParameterDiaRT(
                metaBasePar,
                selectedMagDb,
                librarySearch,
                selfModelSearch,
                library,
                useModel,
                models
            );
            metaParameterDiaRT.setRtWorkflowMode(rtModeCode);

            // DIA-NN modification settings
            if (mExcisionCheckBox != null) {
                metaParameterDiaRT.setMExcision(mExcisionCheckBox.isSelected());
            }
            if (carbamiCCheckBox != null) {
                metaParameterDiaRT.setCarbamiC(carbamiCCheckBox.isSelected());
            }
            if (oxidationMCheckBox != null) {
                metaParameterDiaRT.setOxidationM(oxidationMCheckBox.isSelected());
            }
            if (acetyProNCheckBox != null) {
                metaParameterDiaRT.setAcetyProN(acetyProNCheckBox.isSelected());
            }
            if (phosphoSTYCheckBox != null) {
                metaParameterDiaRT.setPhosphoSTY(phosphoSTYCheckBox.isSelected());
            }
            if (kGGCheckBox != null) {
                metaParameterDiaRT.setKGG(kGGCheckBox.isSelected());
            }
            
            // Set real-time specific parameters
            if (monitorFolderField != null && monitorFolderField.getText() != null && !monitorFolderField.getText().trim().isEmpty()) {
                metaParameterDiaRT.setMonitorFolder(new File(monitorFolderField.getText()));
            }
            
            if (resultFolderField != null && resultFolderField.getText() != null && !resultFolderField.getText().trim().isEmpty()) {
                metaParameterDiaRT.setResultFolder(new File(resultFolderField.getText()));
            }
            
            if (performAnalysisCheckBox != null) {
                metaParameterDiaRT.setCombineResult(performAnalysisCheckBox.isSelected());
            }
            
            return metaParameterDiaRT;
            
        } catch (Exception e) {
            System.err.println("Error generating MetaParameterDiaRT: " + e.getMessage());
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
            MetaParameterDiaRT params = (MetaParameterDiaRT) parameterIO.loadLastUsedParameters();
            if (params != null) {
                System.out.println("DEBUG: Successfully loaded parameters from file");
                applyParametersToUI(params);
                System.out.println("DEBUG: Applied parameters to UI successfully");
            } else {
                System.out.println("DEBUG: No parameters loaded from file, using defaults");
            }
        } catch (Exception e) {
            System.err.println("WARNING: Failed to load last used DIA Real-time parameters: " + e.getMessage());
            e.printStackTrace();
            // Continue with default UI state
        }
    }
    
    /**
     * Apply parameters to the UI components
     */
    public void applyParametersToUI(MetaParameterDiaRT params) {
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
        
        // Apply real-time specific parameters
        if (params.getMonitorFolder() != null) {
            if (monitorFolderField != null) {
                monitorFolderField.setText(params.getMonitorFolder().getAbsolutePath());
                System.out.println("DEBUG: Applied monitor folder: " + params.getMonitorFolder().getAbsolutePath());
            }
        }
        
        if (params.getResultFolder() != null) {
            if (resultFolderField != null) {
                resultFolderField.setText(params.getResultFolder().getAbsolutePath());
                System.out.println("DEBUG: Applied result folder: " + params.getResultFolder().getAbsolutePath());
            }
        }
        
        if (performAnalysisCheckBox != null) {
            performAnalysisCheckBox.setSelected(params.isCombineResult());
            System.out.println("DEBUG: Applied combine result setting: " + params.isCombineResult());
        }
        
        // Apply thread count
        if (threadsComboBox != null) {
            threadsComboBox.setValue(params.getThreadCount());
        }
        
        // Apply workflow mode selection
        if (rtWorkflowModeComboBox != null) {
            String mode = params.getRtWorkflowMode();
            String label;
            if (mode == null || mode.isEmpty()) {
                label = "Seed \u2013 Direct";
            } else {
                switch (mode) {
                    case "MARKER_FOCUS":
                        label = "Marker \u2013 Focus";
                        break;
                    case "MARKER_COMPREHENSIVE":
                        label = "Marker \u2013 Comprehensive";
                        break;
                    case "MARKER_DEEP":
                        label = "Marker \u2013 Deep";
                        break;
                    case "SEED_FOCUS":
                        label = "Seed \u2013 Focus";
                        break;
                    case "SEED_COMPREHENSIVE":
                        label = "Seed \u2013 Comprehensive";
                        break;
                    case "SEED_DEEP":
                        label = "Seed \u2013 Deep";
                        break;
                    case "SEED_DIRECT":
                    default:
                        label = "Seed \u2013 Direct";
                        break;
                }
            }
            rtWorkflowModeComboBox.setValue(label);
        }

        // Apply DIA-NN modification settings
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
            // Refresh library table to show libraries for the restored MAG database
            updatePanelsWithMagDb(magDatabasePanel.getSelectedMagDatabase());
        }
        
        // Update parameter view after applying parameters
        Platform.runLater(() -> updateParameterView());
    }
    
    /**
     * Save current parameters to file
     */
    public void saveCurrentParameters() {
        try {
            MetaParameterDiaRT currentParams = collectParametersFromUIForRealTime();
            if (currentParams != null) {
                boolean saved = parameterIO.saveLastUsedParameters(currentParams);
                if (saved) {
                    System.out.println("DEBUG: Successfully saved current parameters to last_used_dia_rt_parameters.json");
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
    
    private void viewReport() {
        if (consolePanel != null) {
            consolePanel.appendConsoleText("Opening report in default browser...\n");
        }
        
        // Get the result folder
        String resultDir = resultFolderField.getText();
        if (resultDir == null || resultDir.trim().isEmpty()) {
            if (consolePanel != null) {
                consolePanel.appendConsoleText("ERROR: No result directory specified.\n");
            }
            return;
        }
        
        java.io.File reportDir = new java.io.File(resultDir);
        if (!reportDir.exists()) {
            if (consolePanel != null) {
                consolePanel.appendConsoleText("ERROR: Result directory does not exist: " + resultDir + "\n");
            }
            return;
        }
        
        java.io.File indexFile = new java.io.File(reportDir, "index.html");
        if (!indexFile.exists()) {
            if (consolePanel != null) {
                consolePanel.appendConsoleText("ERROR: Report index.html not found in " + resultDir + "\n");
            }
            return;
        }
        
        try {
            if (java.awt.Desktop.isDesktopSupported()) {
                java.awt.Desktop.getDesktop().browse(indexFile.toURI());
                if (consolePanel != null) {
                    consolePanel.appendConsoleText("✓ Report opened in default browser\n");
                }
            } else {
                if (consolePanel != null) {
                    consolePanel.appendConsoleText("ERROR: Desktop not supported. Please open manually: " + indexFile.getAbsolutePath() + "\n");
                }
            }
        } catch (Exception e) {
            if (consolePanel != null) {
                consolePanel.appendConsoleText("ERROR: Failed to open report: " + e.getMessage() + "\n");
            }
        }
    }
    
} 