package uk.ac.quadram.metapilot.task.dia.gui;

import javafx.application.Platform;
import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.control.*;
import javafx.scene.layout.*;
import javafx.stage.FileChooser;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.scene.text.Font;
import javafx.scene.text.FontWeight;
import javafx.scene.paint.Color;
import javafx.scene.control.cell.PropertyValueFactory;
import javafx.scene.control.cell.TextFieldTableCell;
import javafx.scene.control.cell.CheckBoxTableCell;
import javafx.util.converter.DefaultStringConverter;
import javafx.beans.property.SimpleStringProperty;
import javafx.beans.property.StringProperty;
import javafx.beans.property.SimpleBooleanProperty;
import javafx.beans.property.BooleanProperty;
import javafx.concurrent.Task;
import javafx.scene.control.Alert.AlertType;
import uk.ac.quadram.metapilot.task.common.gui.MetaLabIOPanelFX;
import uk.ac.quadram.metapilot.task.common.gui.MetaLabIOPanelFX.ExperimentRow;
import uk.ac.quadram.metapilot.task.par.MetaData;
import uk.ac.quadram.metapilot.task.par.MetaBaseParameter;
import uk.ac.quadram.metapilot.task.par.MetaConstants;
import uk.ac.quadram.metapilot.task.MetaLabWorkflowType;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.task.dda.igc.gui.MetaPilotDdaIgcMaxQuantPanel;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;
import uk.ac.quadram.metapilot.task.mag.par.MetaParameterMag;
import javafx.animation.KeyFrame;
import javafx.animation.Timeline;
import javafx.util.Duration;
import java.lang.management.ManagementFactory;
import java.lang.management.OperatingSystemMXBean;
import java.lang.management.MemoryMXBean;
import java.lang.management.MemoryUsage;
import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Optional;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.io.File;
import uk.ac.quadram.metapilot.gui.JavaFXProgressBar;
import java.util.HashMap;

public class MetaPilotDiaTargetPanel extends VBox {
    // UI Components - Database specific for MAG (like DIA workflows)
    private MetaLabIOPanelFX ioPanel;
    private ComboBox<MagDbItem> magDatabaseComboBox;
    private TextArea consoleOutput;
    private TextArea parameterView;
    private ProgressBar overallProgressBar;
    private ProgressBar currentTaskProgressBar;
    private Label currentTaskLabel;
    private Button startButton;
    private Button stopButton;

    // DIA Target-specific components
    private ComboBox<Integer> threadsComboBox;

    // Library Generation Panel Components
    private ComboBox<Float> peptideDetectabilityComboBox;
    private CheckBox libraryPathCheckBox;
    private TextField libraryPathTextField;
    private Button testThresholdButton;
    private Button createLibraryButton;
    private TextArea libraryConsoleOutput;

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

    // MAG database items
    private ObservableList<MagDbItem> availableMagDatabases = FXCollections.observableArrayList();

    public MetaPilotDiaTargetPanel() {
        initializeUI();
    }

    public MetaPilotDiaTargetPanel(MetaParameterMag metaParameterMag) {
        // Initialize MAG databases from parameter
        if (metaParameterMag != null) {
            MagDbItem[] magDbItems = metaParameterMag.getAvailableMagDbItem();
            if (magDbItems != null) {
                for (MagDbItem item : magDbItems) {
                    availableMagDatabases.add(item);
                }
            }
        }
        initializeUI();
    }

    // Load available MAG databases
    public void setAvailableMagDatabases(MagDbItem[] magDbItems) {
        availableMagDatabases.clear();
        if (magDbItems != null) {
            for (MagDbItem item : magDbItems) {
                availableMagDatabases.add(item);
            }
        }
        if (magDatabaseComboBox != null && !availableMagDatabases.isEmpty()) {
            magDatabaseComboBox.setValue(availableMagDatabases.get(0));
        }
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
        tabPane.getSelectionModel().selectedItemProperty().addListener((obs, oldTab, newTab) -> {
            if (newTab == paramTab) {
                updateMetadataTableFromIOPanel();
            }
        });

        // Add tab pane to main panel
        getChildren().add(tabPane);
    }

    // Create IO tab (same as MAG workflows)
    private VBox createIOTab() {
        VBox tabContent = new VBox(10);
        tabContent.setPadding(new Insets(10));
        ioPanel = new MetaLabIOPanelFX();
        VBox databaseSection = createDatabaseSection();
        tabContent.getChildren().addAll(
            createSectionTitle("Input/Output Settings"),
            ioPanel,
            createSectionTitle("Database Settings"),
            databaseSection
        );
        return tabContent;
    }

    private VBox createDatabaseSection() {
        VBox section = new VBox(10);
        section.setPadding(new Insets(10));
        section.setStyle("-fx-background-color: #f5f5f5; -fx-background-radius: 5;");
        
        // MAG Database ComboBox (not file browser)
        HBox magDBBox = new HBox(10);
        magDatabaseComboBox = new ComboBox<>(availableMagDatabases);
        magDatabaseComboBox.setPromptText("Select MAG database");
        magDatabaseComboBox.setPrefWidth(400);
        magDatabaseComboBox.setOnAction(e -> updateParameterView());
        magDBBox.getChildren().addAll(new Label("MAG Database:"), magDatabaseComboBox);
        
        section.getChildren().add(magDBBox);
        return section;
    }

    private VBox createParametersTab() {
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

        VBox metadataSection = createMetadataConfigurationSection();
        mainGrid.add(metadataSection, 1, 2);
        GridPane.setVgrow(metadataSection, Priority.ALWAYS);
        
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
        peptideDetectabilityComboBox.setMinWidth(80);
        peptideDetectabilityComboBox.setPrefWidth(100);
        peptideDetectabilityComboBox.setOnAction(e -> updateParameterView());
        thresholdBox.getChildren().addAll(thresholdLabel, peptideDetectabilityComboBox);

        // Action buttons
        HBox buttonBox = new HBox(10);
        buttonBox.setAlignment(Pos.CENTER_LEFT);
        testThresholdButton = new Button("Get Proper Threshold");
        testThresholdButton.setPrefWidth(150);
        testThresholdButton.setOnAction(e -> calculateProperThreshold());
        
        createLibraryButton = new Button("Create Library");
        createLibraryButton.setPrefWidth(120);
        createLibraryButton.setOnAction(e -> createLibrary());
        buttonBox.getChildren().addAll(testThresholdButton, createLibraryButton);

        // Library path
        HBox pathBox = new HBox(10);
        pathBox.setAlignment(Pos.CENTER_LEFT);
        libraryPathCheckBox = new CheckBox("Library Path:");
        libraryPathCheckBox.setStyle("-fx-font-weight: bold;");
        libraryPathCheckBox.setSelected(false);
        libraryPathCheckBox.setDisable(true);
        
        libraryPathTextField = new TextField();
        libraryPathTextField.setPromptText("Library will be created here");
        libraryPathTextField.setText("selectedPep.predicted.speclib");
        libraryPathTextField.setDisable(true);
        libraryPathTextField.setPrefWidth(200);
        pathBox.getChildren().addAll(libraryPathCheckBox, libraryPathTextField);

        // Console output
        Label consoleLabel = new Label("Progress Output:");
        consoleLabel.setStyle("-fx-font-weight: bold; -fx-font-size: 11px;");
        libraryConsoleOutput = new TextArea();
        libraryConsoleOutput.setEditable(false);
        libraryConsoleOutput.setPrefHeight(80);
        libraryConsoleOutput.setMaxHeight(120);
        libraryConsoleOutput.setStyle("-fx-font-family: monospace; -fx-font-size: 10px;");
        VBox.setVgrow(libraryConsoleOutput, Priority.ALWAYS);

        section.getChildren().addAll(titleLabel, thresholdBox, buttonBox, pathBox, consoleLabel, libraryConsoleOutput);
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
        taxaLevelComboBox = new ComboBox<>(FXCollections.observableArrayList(
            "Choose the taxonomic level", "Superkingdom", "Phylum", "Class", 
            "Order", "Family", "Genus", "Species", "Genome"));
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
        taxaTable.setColumnResizePolicy(TableView.CONSTRAINED_RESIZE_POLICY);
        VBox.setVgrow(taxaTable, Priority.ALWAYS);

        // Select column (checkbox)
        TableColumn<TaxaRow, Boolean> selectCol = new TableColumn<>("Select");
        selectCol.setCellValueFactory(new PropertyValueFactory<>("selected"));
        selectCol.setCellFactory(CheckBoxTableCell.forTableColumn(selectCol));
        selectCol.setOnEditCommit(e -> updateParameterView());
        selectCol.setMinWidth(60);
        selectCol.setPrefWidth(80);
        selectCol.setSortable(false);

        // Taxa Name column
        TableColumn<TaxaRow, String> nameCol = new TableColumn<>("Taxonomic Name");
        nameCol.setCellValueFactory(new PropertyValueFactory<>("name"));
        nameCol.setMinWidth(200);
        nameCol.setPrefWidth(300);
        nameCol.setSortable(true);

        taxaTable.getColumns().addAll(selectCol, nameCol);

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
        functionLevelComboBox = new ComboBox<>(FXCollections.observableArrayList(
            "Choose the functional level", "COG", "KEGG", "Pfam", "GO"));
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
        functionTable.setColumnResizePolicy(TableView.CONSTRAINED_RESIZE_POLICY);
        VBox.setVgrow(functionTable, Priority.ALWAYS);

        // Select column (checkbox)
        TableColumn<FunctionRow, Boolean> selectCol = new TableColumn<>("Select");
        selectCol.setCellValueFactory(new PropertyValueFactory<>("selected"));
        selectCol.setCellFactory(CheckBoxTableCell.forTableColumn(selectCol));
        selectCol.setOnEditCommit(e -> updateParameterView());
        selectCol.setMinWidth(60);
        selectCol.setPrefWidth(80);
        selectCol.setSortable(false);

        // Function Name column
        TableColumn<FunctionRow, String> nameCol = new TableColumn<>("Function Name");
        nameCol.setCellValueFactory(new PropertyValueFactory<>("name"));
        nameCol.setMinWidth(200);
        nameCol.setPrefWidth(300);
        nameCol.setSortable(true);

        functionTable.getColumns().addAll(selectCol, nameCol);

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
        
        public boolean getSelected() { return selected.get(); }
        public void setSelected(boolean value) { selected.set(value); }
        public BooleanProperty selectedProperty() { return selected; }
        
        public String getName() { return name.get(); }
        public void setName(String value) { name.set(value); }
        public StringProperty nameProperty() { return name; }
    }

    // Function Row class
    public static class FunctionRow {
        private final SimpleBooleanProperty selected;
        private final SimpleStringProperty name;
        
        public FunctionRow(boolean selected, String name) {
            this.selected = new SimpleBooleanProperty(selected);
            this.name = new SimpleStringProperty(name);
        }
        
        public boolean getSelected() { return selected.get(); }
        public void setSelected(boolean value) { selected.set(value); }
        public BooleanProperty selectedProperty() { return selected; }
        
        public String getName() { return name.get(); }
        public void setName(String value) { name.set(value); }
        public StringProperty nameProperty() { return name; }
    }

    private void loadTaxaData() {
        String selectedLevel = taxaLevelComboBox.getValue();
        taxaData.clear();
        
        if (selectedLevel == null || selectedLevel.equals("Choose the taxonomic level")) {
            return;
        }
        
        // Mock data for demonstration - in real implementation, this would query the database
        String[] mockTaxa = {"Bacteria", "Archaea", "Eukaryota", "Viruses"};
        for (String taxon : mockTaxa) {
            taxaData.add(new TaxaRow(false, taxon));
        }
        
        updateParameterView();
    }

    private void loadFunctionData() {
        String selectedLevel = functionLevelComboBox.getValue();
        functionData.clear();
        
        if (selectedLevel == null || selectedLevel.equals("Choose the functional level")) {
            return;
        }
        
        // Mock data for demonstration - in real implementation, this would query the database
        String[] mockFunctions = {"Energy metabolism", "Protein synthesis", "DNA replication", "Cell division"};
        for (String function : mockFunctions) {
            functionData.add(new FunctionRow(false, function));
        }
        
        updateParameterView();
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

    // Add remaining missing methods
    private VBox createConsoleTab() {
        VBox tabContent = new VBox(10);
        tabContent.setPadding(new Insets(10));

        // Parameter View
        VBox paramViewBox = new VBox(5);
        paramViewBox.getChildren().add(new Label("Parameter View:"));
        parameterView = new TextArea();
        parameterView.setEditable(false);
        parameterView.setPrefHeight(200);
        parameterView.setStyle("-fx-font-family: monospace;");
        VBox.setVgrow(parameterView, Priority.ALWAYS);
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
        consoleOutput.setPrefHeight(250);
        consoleOutput.setStyle("-fx-font-family: monospace;");
        VBox.setVgrow(consoleOutput, Priority.ALWAYS);
        consoleBox.getChildren().add(consoleOutput);
        VBox.setVgrow(consoleBox, Priority.ALWAYS);

        // Control Buttons
        HBox controlBox = new HBox(10);
        startButton = new Button("Start Analysis");
        stopButton = new Button("Stop");
        stopButton.setDisable(true);
        startButton.setOnAction(e -> startAnalysis());
        stopButton.setOnAction(e -> stopAnalysis());
        controlBox.getChildren().addAll(startButton, stopButton);

        // Status Panel
        HBox statusPanel = new HBox(20);
        statusPanel.setPadding(new Insets(5, 0, 0, 0));
        statusPanel.setAlignment(Pos.CENTER_LEFT);
        statusTimeLabel = new Label("Time: 00:00:00");
        statusCpuLabel = new Label("CPU: 0.0%");
        statusMemLabel = new Label("Memory: 0 / 0 MB");
        statusPanel.getChildren().addAll(statusTimeLabel, statusCpuLabel, statusMemLabel);

        tabContent.getChildren().addAll(
            paramViewBox,
            progressBox,
            consoleBox,
            controlBox,
            statusPanel
        );
        VBox.setVgrow(paramViewBox, Priority.ALWAYS);
        VBox.setVgrow(consoleBox, Priority.ALWAYS);
        VBox.setVgrow(tabContent, Priority.ALWAYS);

        return tabContent;
    }

    private void updateMetadataTableFromIOPanel() {
        if (ioPanel == null || metadataTable == null) {
            return;
        }
        
        ObservableList<MetaLabIOPanelFX.ExperimentRow> experiments = ioPanel.getExperimentData();
        metadataTableData.clear();
        for (MetaLabIOPanelFX.ExperimentRow exp : experiments) {
            metadataTableData.add(new MetaPilotDdaIgcMaxQuantPanel.MetadataRow(exp.getExperimentName()));
        }
    }

    private VBox createMetadataConfigurationSection() {
        VBox section = new VBox(10);
        section.setPadding(new Insets(10));
        section.setStyle("-fx-background-color: #fff8f0; -fx-background-radius: 5;");
        section.setMinHeight(180);
        VBox.setVgrow(section, Priority.ALWAYS);

        Label titleLabel = createSectionTitle("Metadata Configuration");
        Label infoLabel = new Label("Configure metadata for DIA Target analysis samples");
        infoLabel.setStyle("-fx-font-style: italic; -fx-text-fill: #666666; -fx-font-size: 11px;");

        metadataTable = new TableView<>(metadataTableData);
        metadataTable.setMinHeight(80);
        metadataTable.setPrefHeight(100);
        metadataTable.setMaxHeight(Double.MAX_VALUE);
        metadataTable.setEditable(true);
        metadataTable.setColumnResizePolicy(TableView.CONSTRAINED_RESIZE_POLICY);
        VBox.setVgrow(metadataTable, Priority.ALWAYS);

        createMetadataTableColumns();

        // Compact control layout
        VBox controlContainer = new VBox(5);

        // Row 1: Column selection
        HBox row1 = new HBox(5);
        row1.setAlignment(Pos.CENTER_LEFT);

        metadataColumnComboBox = new ComboBox<>(FXCollections.observableArrayList(
            "Age", "Gender", "Condition", "Treatment", "Time Point", "BMI", "Custom"));
        metadataColumnComboBox.setMinWidth(80);
        metadataColumnComboBox.setPrefWidth(100);
        metadataColumnComboBox.setPromptText("Add column");

        TextField customColumnTextField = new TextField();
        customColumnTextField.setPromptText("Custom name");
        customColumnTextField.setMinWidth(80);
        customColumnTextField.setPrefWidth(100);
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

        row1.getChildren().addAll(metadataColumnComboBox, customColumnTextField);

        // Row 2: Action buttons
        HBox row2 = new HBox(5);
        row2.setAlignment(Pos.CENTER_LEFT);

        addMetadataColumnButton = new Button("Add Column");
        addMetadataColumnButton.setPrefWidth(100);
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
        removeMetadataColumnButton.setPrefWidth(100);
        removeMetadataColumnButton.setOnAction(e -> removeMetadataColumn());

        row2.getChildren().addAll(addMetadataColumnButton, removeMetadataColumnButton);

        controlContainer.getChildren().addAll(row1, row2);

        section.getChildren().addAll(titleLabel, infoLabel, metadataTable, controlContainer);
        return section;
    }

    private void createMetadataTableColumns() {
        metadataTable.getColumns().clear();
        TableColumn<MetaPilotDdaIgcMaxQuantPanel.MetadataRow, String> sampleCol = new TableColumn<>("Sample");
        sampleCol.setCellValueFactory(new PropertyValueFactory<>("sampleName"));
        sampleCol.setMinWidth(120);
        sampleCol.setPrefWidth(150);
        metadataTable.getColumns().add(sampleCol);
    }

    private void updateParameterView() {
        if (parameterView == null) {
            return;
        }
        
        StringBuilder sb = new StringBuilder();
        sb.append("=== DIA Target Analysis Parameters ===\n\n");

        // Database settings
        sb.append("Database Settings:\n");
        sb.append("  MAG Database: ").append(magDatabaseComboBox != null ? magDatabaseComboBox.getValue() != null ? magDatabaseComboBox.getValue().getCatalogueID() : "" : "").append("\n");
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
            sb.append("  Peptide Detectability Threshold: ").append(peptideDetectabilityComboBox.getValue()).append("\n");
        }
        if (libraryPathTextField != null && !libraryPathTextField.getText().isEmpty()) {
            sb.append("  Library Path: ").append(libraryPathTextField.getText()).append("\n");
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
            sb.append("  Selected Functions: ").append(selectedCount).append(" of ").append(functionData.size()).append("\n");
        }
        sb.append("\n");

        parameterView.setText(sb.toString());
    }

    private String getPropertyNameForColumn(String columnName) {
        switch (columnName) {
            case "Age": return "age";
            case "Gender": return "gender";
            case "Condition": return "condition";
            case "Treatment": return "treatment";
            case "Time Point": return "timePoint";
            case "BMI": return "customField1";
            default: return "customField1";
        }
    }

    private void addMetadataColumnWithProperty(String columnName, String propertyName) {
        TableColumn<MetaPilotDdaIgcMaxQuantPanel.MetadataRow, String> newCol = new TableColumn<>(columnName);
        newCol.setCellValueFactory(data -> {
            switch (propertyName) {
                case "age": return new SimpleStringProperty(data.getValue().getAge());
                case "gender": return new SimpleStringProperty(data.getValue().getGender());
                case "condition": return new SimpleStringProperty(data.getValue().getCondition());
                case "treatment": return new SimpleStringProperty(data.getValue().getTreatment());
                case "timePoint": return new SimpleStringProperty(data.getValue().getTimePoint());
                case "customField1": return new SimpleStringProperty(data.getValue().getCustomField1());
                case "customField2": return new SimpleStringProperty(data.getValue().getCustomField2());
                case "customField3": return new SimpleStringProperty(data.getValue().getCustomField3());
                default: return new SimpleStringProperty("");
            }
        });
        newCol.setCellFactory(TextFieldTableCell.forTableColumn());
        newCol.setOnEditCommit(e -> {
            switch (propertyName) {
                case "age": e.getRowValue().setAge(e.getNewValue()); break;
                case "gender": e.getRowValue().setGender(e.getNewValue()); break;
                case "condition": e.getRowValue().setCondition(e.getNewValue()); break;
                case "treatment": e.getRowValue().setTreatment(e.getNewValue()); break;
                case "timePoint": e.getRowValue().setTimePoint(e.getNewValue()); break;
                case "customField1": e.getRowValue().setCustomField1(e.getNewValue()); break;
                case "customField2": e.getRowValue().setCustomField2(e.getNewValue()); break;
                case "customField3": e.getRowValue().setCustomField3(e.getNewValue()); break;
            }
        });
        newCol.setMinWidth(100);
        newCol.setPrefWidth(120);
        metadataTable.getColumns().add(newCol);
    }

    private void removeMetadataColumn() {
        TableColumn<MetaPilotDdaIgcMaxQuantPanel.MetadataRow, ?> selectedColumn = 
            metadataTable.getSelectionModel().getSelectedCells().isEmpty() ? null :
            metadataTable.getColumns().get(metadataTable.getSelectionModel().getSelectedCells().get(0).getColumn());
        
        if (selectedColumn != null && !selectedColumn.getText().equals("Sample")) {
            metadataTable.getColumns().remove(selectedColumn);
        }
    }

    private void startAnalysis() {
        startButton.setDisable(true);
        stopButton.setDisable(false);
        consoleOutput.clear();
        consoleOutput.appendText("Starting DIA Target analysis...\n");
        System.out.println("DEBUG: startAnalysis called - DIA Target workflow");
    }

    private void stopAnalysis() {
        startButton.setDisable(false);
        stopButton.setDisable(true);
        consoleOutput.appendText("Stopping DIA Target analysis...\n");
        System.out.println("DEBUG: stopAnalysis called");
    }

    private void calculateProperThreshold() {
        // Validate inputs
        if (!validateSelections()) {
            return;
        }
        
        // Disable buttons during calculation
        testThresholdButton.setDisable(true);
        createLibraryButton.setDisable(true);
        libraryConsoleOutput.clear();
        libraryConsoleOutput.appendText("=== Calculating Proper Peptide Detectability Threshold ===\n");
        libraryConsoleOutput.appendText("Analyzing selected taxa and functions...\n");
        
        // Create background task for threshold calculation
        Task<Float> task = new Task<Float>() {
            @Override
            protected Float call() throws Exception {
                // Mock calculation - in real implementation, this would query the database
                // and calculate optimal threshold based on selected taxa and functions
                
                Platform.runLater(() -> {
                    libraryConsoleOutput.appendText("Querying MAG database for selected taxa...\n");
                });
                Thread.sleep(1000);
                
                Platform.runLater(() -> {
                    libraryConsoleOutput.appendText("Processing functional annotations...\n");
                });
                Thread.sleep(1000);
                
                Platform.runLater(() -> {
                    libraryConsoleOutput.appendText("Calculating peptide detectability distributions...\n");
                });
                Thread.sleep(1000);
                
                // Mock recommendation: suggest 0.7 as proper threshold
                Platform.runLater(() -> {
                    libraryConsoleOutput.appendText("Analysis complete!\n");
                    libraryConsoleOutput.appendText("Recommended threshold: 0.70\n");
                    libraryConsoleOutput.appendText("The proper size of a peptide library is from 1,000,000 to 3,000,000 peptides\n");
                    libraryConsoleOutput.appendText("=== Threshold Calculation Complete ===\n\n");
                });
                
                return 0.70f;
            }
            
            @Override
            protected void succeeded() {
                Float recommendedThreshold = getValue();
                if (recommendedThreshold != null) {
                    peptideDetectabilityComboBox.setValue(recommendedThreshold);
                    updateParameterView();
                }
                testThresholdButton.setDisable(false);
                createLibraryButton.setDisable(false);
            }
            
            @Override
            protected void failed() {
                libraryConsoleOutput.appendText("ERROR: Failed to calculate threshold - " + getException().getMessage() + "\n");
                testThresholdButton.setDisable(false);
                createLibraryButton.setDisable(false);
            }
        };
        
        new Thread(task).start();
    }
    
    private void createLibrary() {
        // Validate inputs
        if (!validateSelections()) {
            return;
        }
        
        // Check if library file already exists
        if (libraryPathCheckBox.isSelected()) {
            Alert confirm = new Alert(Alert.AlertType.CONFIRMATION);
            confirm.setTitle("Overwrite Library");
            confirm.setHeaderText("Library Already Exists");
            confirm.setContentText("A library already exists at " + libraryPathTextField.getText() + 
                                 ". Do you want to overwrite it?");
            
            Optional<ButtonType> result = confirm.showAndWait();
            if (result.isEmpty() || result.get() != ButtonType.OK) {
                return;
            }
        }
        
        // Disable buttons during library creation
        testThresholdButton.setDisable(true);
        createLibraryButton.setDisable(true);
        libraryConsoleOutput.clear();
        libraryConsoleOutput.appendText("=== Creating Targeted Peptide Library ===\n");
        
        // Create background task for library creation
        Task<Boolean> task = new Task<Boolean>() {
            @Override
            protected Boolean call() throws Exception {
                // Mock library creation - in real implementation, this would:
                // 1. Query database for selected taxa/functions
                // 2. Filter peptides by detectability threshold
                // 3. Generate spectral library file
                
                Platform.runLater(() -> {
                    libraryConsoleOutput.appendText("Step 1: Querying selected taxa and functions...\n");
                });
                Thread.sleep(1500);
                
                Platform.runLater(() -> {
                    libraryConsoleOutput.appendText("Found " + getSelectedTaxaCount() + " taxonomic groups\n");
                    libraryConsoleOutput.appendText("Found " + getSelectedFunctionCount() + " functional groups\n");
                    libraryConsoleOutput.appendText("Step 2: Retrieving protein sequences...\n");
                });
                Thread.sleep(1500);
                
                Platform.runLater(() -> {
                    libraryConsoleOutput.appendText("Step 3: Generating peptide predictions...\n");
                });
                Thread.sleep(1500);
                
                Platform.runLater(() -> {
                    libraryConsoleOutput.appendText("Step 4: Filtering by detectability threshold (" + 
                                           peptideDetectabilityComboBox.getValue() + ")...\n");
                });
                Thread.sleep(1000);
                
                Platform.runLater(() -> {
                    libraryConsoleOutput.appendText("Step 5: Creating spectral library format...\n");
                });
                Thread.sleep(1000);
                
                Platform.runLater(() -> {
                    libraryConsoleOutput.appendText("Library created successfully: " + libraryPathTextField.getText() + "\n");
                    libraryConsoleOutput.appendText("=== Library Generation Complete ===\n");
                });
                
                return true;
            }
            
            @Override
            protected void succeeded() {
                if (getValue()) {
                    libraryPathCheckBox.setSelected(true);
                    Alert success = new Alert(Alert.AlertType.INFORMATION);
                    success.setTitle("Library Created");
                    success.setHeaderText("Success");
                    success.setContentText("Targeted peptide library has been created successfully!");
                    success.showAndWait();
                }
                testThresholdButton.setDisable(false);
                createLibraryButton.setDisable(false);
            }
            
            @Override
            protected void failed() {
                libraryConsoleOutput.appendText("ERROR: Failed to create library - " + getException().getMessage() + "\n");
                testThresholdButton.setDisable(false);
                createLibraryButton.setDisable(false);
            }
        };
        
        new Thread(task).start();
    }
    
    private boolean validateSelections() {
        // Check if taxa are selected
        if (taxaLevelComboBox.getValue() == null || 
            taxaLevelComboBox.getValue().equals("Choose the taxonomic level")) {
            showValidationError("Please specify a taxonomic level and select taxa.");
            return false;
        }
        
        int selectedTaxa = 0;
        for (TaxaRow row : taxaData) {
            if (row.getSelected()) {
                selectedTaxa++;
            }
        }
        if (selectedTaxa == 0) {
            showValidationError("Please select at least one taxonomic group.");
            return false;
        }
        
        // Check if functions are selected
        if (functionLevelComboBox.getValue() == null || 
            functionLevelComboBox.getValue().equals("Choose the functional level")) {
            showValidationError("Please specify a functional level and select functions.");
            return false;
        }
        
        int selectedFunctions = 0;
        for (FunctionRow row : functionData) {
            if (row.getSelected()) {
                selectedFunctions++;
            }
        }
        if (selectedFunctions == 0) {
            showValidationError("Please select at least one functional group.");
            return false;
        }
        
        return true;
    }
    
    private void showValidationError(String message) {
        Alert alert = new Alert(Alert.AlertType.ERROR);
        alert.setTitle("Validation Error");
        alert.setHeaderText("Invalid Selection");
        alert.setContentText(message);
        alert.showAndWait();
    }

    private int getSelectedTaxaCount() {
        int count = 0;
        for (TaxaRow row : taxaData) {
            if (row.getSelected()) {
                count++;
            }
        }
        return count;
    }

    private int getSelectedFunctionCount() {
        int count = 0;
        for (FunctionRow row : functionData) {
            if (row.getSelected()) {
                count++;
            }
        }
        return count;
    }
} 