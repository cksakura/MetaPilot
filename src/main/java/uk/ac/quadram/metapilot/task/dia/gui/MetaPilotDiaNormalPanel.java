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
import javafx.util.converter.DoubleStringConverter;
import javafx.util.converter.IntegerStringConverter;
import javafx.beans.property.SimpleStringProperty;
import javafx.beans.property.StringProperty;
import javafx.beans.property.SimpleBooleanProperty;
import javafx.beans.property.BooleanProperty;
import javafx.beans.property.SimpleDoubleProperty;
import javafx.beans.property.DoubleProperty;
import javafx.beans.property.SimpleIntegerProperty;
import javafx.beans.property.IntegerProperty;
import javafx.concurrent.Task;
import javafx.scene.control.Alert.AlertType;
import uk.ac.quadram.metapilot.task.common.gui.MetaPilotIOPanelFX;
import uk.ac.quadram.metapilot.task.common.gui.MetaPilotIOPanelFX.ExperimentRow;
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

public class MetaPilotDiaNormalPanel extends VBox {
    // UI Components - Database specific for MAG (like DIA workflows)
    private MetaPilotIOPanelFX ioPanel;
    private ComboBox<MagDbItem> magDatabaseComboBox;
    private TextArea consoleOutput;
    private TextArea parameterView;
    private ProgressBar overallProgressBar;
    private ProgressBar currentTaskProgressBar;
    private Label currentTaskLabel;
    private Button startButton;
    private Button stopButton;

    // DIA-specific mode selection
    private ComboBox<String> modeSelectionComboBox;
    private ComboBox<Integer> threadsComboBox;

    // Separate panels for library and model generation
    private DiaLibraryGenerationPanel libraryGenerationPanel;
    private DiaModelGenerationPanel modelGenerationPanel;

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

    public MetaPilotDiaNormalPanel() {
        initializeUI();
    }

    public MetaPilotDiaNormalPanel(MetaParameterMag metaParameterMag) {
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
            // Update the panels with the selected MAG database
            updatePanelsWithMagDb(availableMagDatabases.get(0));
        }
    }

    private void updatePanelsWithMagDb(MagDbItem magDbItem) {
        if (libraryGenerationPanel != null) {
            libraryGenerationPanel.setMagDbItem(magDbItem);
        }
        if (modelGenerationPanel != null) {
            modelGenerationPanel.setMagDbItem(magDbItem);
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
        ioPanel = new MetaPilotIOPanelFX();
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
        section.setStyle("-fx-background-color: #f8f8f8; -fx-border-color: #cccccc; -fx-border-radius: 5; -fx-border-width: 1;");
        
        // MAG Database ComboBox (not file browser)
        HBox magDBBox = new HBox(10);
        magDatabaseComboBox = new ComboBox<>(availableMagDatabases);
        magDatabaseComboBox.setPromptText("Select MAG database");
        magDatabaseComboBox.setPrefWidth(400);
        magDatabaseComboBox.setOnAction(e -> updateParameterView());
        magDBBox.getChildren().addAll(new Label("MAG Database:"), magDatabaseComboBox);
        
        section.getChildren().addAll(magDBBox);
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

        // Row constraints
        RowConstraints row1 = new RowConstraints();
        row1.setVgrow(Priority.ALWAYS);
        row1.setMinHeight(400);
        mainGrid.getRowConstraints().add(row1);

        // Left column: DIA Mode Selection, Thread Settings, Library Panel
        VBox leftColumn = new VBox(15);
        leftColumn.setMinWidth(430);
        VBox diaSettingsSection = createDiaSettingsSection();
        
        // Use the new reusable library generation panel
        libraryGenerationPanel = new DiaLibraryGenerationPanel();
        libraryGenerationPanel.setParameterUpdateCallback(() -> updateParameterView());
        
        leftColumn.getChildren().addAll(diaSettingsSection, libraryGenerationPanel);
        GridPane.setVgrow(leftColumn, Priority.ALWAYS);

        // Right column: Model Panel and Metadata Configuration
        VBox rightColumn = new VBox(15);
        rightColumn.setMinWidth(430);
        
        // Use the new reusable model generation panel
        modelGenerationPanel = new DiaModelGenerationPanel();
        modelGenerationPanel.setParameterUpdateCallback(() -> updateParameterView());
        
        VBox metadataSection = createMetadataConfigurationSection();
        rightColumn.getChildren().addAll(modelGenerationPanel, metadataSection);
        GridPane.setVgrow(rightColumn, Priority.ALWAYS);

        // Add columns to grid
        mainGrid.add(leftColumn, 0, 0);
        mainGrid.add(rightColumn, 1, 0);
        
        tabContent.getChildren().add(mainGrid);
        return tabContent;
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

        // Mode Selection
        Label modeLabel = new Label("Analysis Mode:");
        modeLabel.setStyle("-fx-font-weight: bold;");
        modeSelectionComboBox = new ComboBox<>(FXCollections.observableArrayList(
            "Directly identification (no iterative search)",
            "Iterative search", 
            "Deep learning mode"));
        modeSelectionComboBox.setValue("Directly identification (no iterative search)");
        modeSelectionComboBox.setMinWidth(200);
        modeSelectionComboBox.setPrefWidth(300);
        modeSelectionComboBox.setMaxWidth(Double.MAX_VALUE);
        modeSelectionComboBox.setOnAction(e -> updateParameterView());

        VBox modeBox = new VBox(5);
        modeBox.getChildren().addAll(modeLabel, modeSelectionComboBox);

        // Thread Count
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

        section.getChildren().addAll(titleLabel, settingsGrid);
        return section;
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
        
        ObservableList<MetaPilotIOPanelFX.ExperimentRow> experiments = ioPanel.getExperimentData();
        metadataTableData.clear();
        for (MetaPilotIOPanelFX.ExperimentRow exp : experiments) {
            metadataTableData.add(new MetaPilotDdaIgcMaxQuantPanel.MetadataRow(exp.getExperimentName()));
        }
    }

    private VBox createMetadataConfigurationSection() {
        VBox section = new VBox(10);
        section.setPadding(new Insets(10));
        section.setStyle("-fx-background-color: #fff8f0; -fx-background-radius: 5;");
        section.setMinHeight(150);
        VBox.setVgrow(section, Priority.ALWAYS);

        Label titleLabel = createSectionTitle("Metadata Configuration");
        Label infoLabel = new Label("Configure metadata for DIA analysis samples");
        infoLabel.setStyle("-fx-font-style: italic; -fx-text-fill: #666666;");

        metadataTable = new TableView<>(metadataTableData);
        metadataTable.setMinHeight(100);
        metadataTable.setPrefHeight(120);
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

        controlBox.getChildren().addAll(
            metadataColumnComboBox, 
            customColumnTextField, 
            addMetadataColumnButton, 
            removeMetadataColumnButton
        );

        section.getChildren().addAll(titleLabel, infoLabel, metadataTable, controlBox);
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
        sb.append("=== DIA Analysis Parameters ===\n\n");

        // Database settings
        sb.append("Database Settings:\n");
        sb.append("  MAG Database: ").append(magDatabaseComboBox != null ? magDatabaseComboBox.getValue() != null ? magDatabaseComboBox.getValue().getCatalogueID() : "" : "").append("\n");
        sb.append("\n");

        // DIA settings
        sb.append("DIA Settings:\n");
        if (modeSelectionComboBox != null) {
            sb.append("  Analysis Mode: ").append(modeSelectionComboBox.getValue()).append("\n");
        }
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
        }
        sb.append("\n");

        // Model configuration
        sb.append("ML Models:\n");
        if (modelGenerationPanel != null) {
            String[] selectedModels = modelGenerationPanel.getSelectedModelNames();
            if (selectedModels.length > 0) {
                for (String modelName : selectedModels) {
                    sb.append("  ").append(modelName).append("\n");
                }
            } else {
                sb.append("  No models selected (models are optional)\n");
            }
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
        consoleOutput.appendText("Starting DIA analysis...\n");
        System.out.println("DEBUG: startAnalysis called - DIA Normal workflow");
    }

    private void stopAnalysis() {
        startButton.setDisable(false);
        stopButton.setDisable(true);
        consoleOutput.appendText("Stopping DIA analysis...\n");
        System.out.println("DEBUG: stopAnalysis called");
    }
} 