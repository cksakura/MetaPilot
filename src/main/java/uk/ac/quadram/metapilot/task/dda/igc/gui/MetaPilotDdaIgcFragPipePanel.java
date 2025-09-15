package uk.ac.quadram.metapilot.task.dda.igc.gui;

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
import javafx.scene.control.cell.PropertyValueFactory;
import javafx.scene.control.cell.TextFieldTableCell;
import javafx.beans.property.SimpleStringProperty;
import uk.ac.quadram.metapilot.task.common.gui.MetaPilotIOPanelFX;
import uk.ac.quadram.metapilot.task.common.gui.TaxonomicAnalysisPanel;
import uk.ac.quadram.metapilot.task.common.gui.DdaIgcInputOutputPanel;
import uk.ac.quadram.metapilot.task.common.gui.ModificationsPanel;
import uk.ac.quadram.metapilot.task.dda.igc.gui.MetaPilotDdaIgcMaxQuantPanel;
import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.IOException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.io.File;

public class MetaPilotDdaIgcFragPipePanel extends VBox {
    // UI Components (similar to pFind panel)
    private DdaIgcInputOutputPanel ioPanel;
    private TextArea consoleOutput;
    private TextArea parameterView;
    private ProgressBar overallProgressBar;
    private ProgressBar currentTaskProgressBar;
    private Label currentTaskLabel;
    private Button startButton;
    private Button stopButton;

    // Parameter Panel Components
    private ComboBox<String> enzymeComboBox;
    private Spinner<Integer> missedCleavagesSpinner;
    private ComboBox<String> digestionModeComboBox;

    // Reusable modifications panel
    private ModificationsPanel modificationsPanel;

    // Additional parameters
    private Spinner<Integer> threadsSpinner;
    private ComboBox<String> workflowComboBox;
    private ComboBox<String> msDataTypeComboBox;
    private CheckBox generateSampleSpecificDBCheckBox;
    private CheckBox spectraClusteringCheckBox;

    // Data
    private ObservableList<String> enzymes = FXCollections.observableArrayList();
    private ObservableList<String> workflows = FXCollections.observableArrayList();


    // Taxonomic analysis panel (unified component)
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

    public MetaPilotDdaIgcFragPipePanel() {
        loadEnzymesFromIni();
        loadWorkflowsFromFragPipe();
        
        // Initialize unified IO panel
        ioPanel = new DdaIgcInputOutputPanel(this::updateParameterView);
        
        // Initialize taxonomic analysis panel
        taxonomicAnalysisPanel = new TaxonomicAnalysisPanel(this::updateParameterView);
        
        initializeUI();
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

    // Create IO tab (similar to pFind)
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

        // Left column: Enzyme Settings, Additional Parameters, Taxonomic Analysis, and Modifications Panel
        VBox leftColumn = new VBox(15);
        leftColumn.setMinWidth(430);
        VBox enzymeSection = createEnzymeSection();
        VBox additionalParamsSection = createAdditionalParametersSection();
        
        // Create horizontal modifications panel with fixed mods (narrow) and variable mods (wider)
        modificationsPanel = new ModificationsPanel(this::updateParameterView);
        
        leftColumn.getChildren().addAll(enzymeSection, additionalParamsSection, taxonomicAnalysisPanel, modificationsPanel);
        GridPane.setVgrow(leftColumn, Priority.ALWAYS);

        // Right column: Metadata Configuration only
        VBox rightColumn = new VBox(15);
        rightColumn.setMinWidth(430);
        VBox metadataSection = createMetadataConfigurationSection();
        rightColumn.getChildren().addAll(metadataSection);
        GridPane.setVgrow(rightColumn, Priority.ALWAYS);

        // Add columns to grid
        mainGrid.add(leftColumn, 0, 0);
        mainGrid.add(rightColumn, 1, 0);
        
        tabContent.getChildren().add(mainGrid);
        return tabContent;
    }
    
    // Combined modifications panel with fixed mods (narrow) and variable mods (wider) side by side

    // Add missing utility methods
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

    private void loadWorkflowsFromFragPipe() {
        try {
            // Get FragPipe path from MetaResourceParameters
            uk.ac.quadram.metapilot.resources.MetaResourceIO resourceIO = new uk.ac.quadram.metapilot.resources.MetaResourceIO();
            uk.ac.quadram.metapilot.resources.MetaResourceParameters mrp = resourceIO.loadParameters();
            String fragpipePath = mrp.getToolPath(uk.ac.quadram.metapilot.resources.MetaResourceParameters.ToolType.FRAGPIPE);

            if (fragpipePath != null && !fragpipePath.trim().isEmpty()) {
                File fragpipeFile = new File(fragpipePath);
                if (fragpipeFile.exists()) {
                    // Get parent directory (bin)
                    File binDir = fragpipeFile.getParentFile();
                    if (binDir != null) {
                        // Get workflows directory (same level as bin)
                        File workflowsDir = new File(binDir.getParentFile(), "workflows");
                        if (workflowsDir.exists() && workflowsDir.isDirectory()) {
                            // List all .workflow files
                            File[] workflowFiles = workflowsDir
                                    .listFiles((dir, name) -> name.toLowerCase().endsWith(".workflow"));
                            if (workflowFiles != null) {
                                for (File workflowFile : workflowFiles) {
                                    String fileName = workflowFile.getName();
                                    // Remove .workflow extension
                                    String workflowName = fileName.substring(0, fileName.lastIndexOf('.'));
                                    workflows.add(workflowName);
                                }
                                System.out.println("Loaded " + workflows.size() + " workflows from FragPipe: " + workflowsDir.getAbsolutePath());
                            }
                        } else {
                            System.out.println("FragPipe workflows directory not found: " + workflowsDir.getAbsolutePath());
                        }
                    }
                } else {
                    System.out.println("FragPipe executable not found: " + fragpipePath);
                }
            } else {
                System.out.println("FragPipe path not configured in MetaResourceParameters");
            }
        } catch (Exception e) {
            System.out.println("Failed to load workflows from FragPipe: " + e.getMessage());
            e.printStackTrace();
        }

        // Add default workflow if none found
        if (workflows.isEmpty()) {
            workflows.add("Standard");
            System.out.println("No workflows found, using default: Standard");
        }
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

        Label workflowLabel = new Label("Workflow:");
        workflowLabel.setStyle("-fx-font-weight: bold;");
        workflowComboBox = new ComboBox<>(workflows);
        workflowComboBox.setMinWidth(100);
        workflowComboBox.setPrefWidth(120);
        workflowComboBox.setMaxWidth(Double.MAX_VALUE);
        if (workflows.size() > 0) {
            workflowComboBox.setValue(workflows.get(0));
        }
        workflowComboBox.setOnAction(e -> updateParameterView());

        VBox workflowBox = new VBox(5);
        workflowBox.getChildren().addAll(workflowLabel, workflowComboBox);

        // Add MS data type option
        Label msDataTypeLabel = new Label("MS data type:");
        msDataTypeLabel.setStyle("-fx-font-weight: bold;");
        msDataTypeComboBox = new ComboBox<>(FXCollections.observableArrayList("Regular MS", "IM-MS (ion mobility, timsTOF only)"));
        msDataTypeComboBox.setValue("Regular MS");
        msDataTypeComboBox.setMinWidth(100);
        msDataTypeComboBox.setPrefWidth(200);
        msDataTypeComboBox.setMaxWidth(Double.MAX_VALUE);
        msDataTypeComboBox.setOnAction(e -> updateParameterView());

        VBox msDataTypeBox = new VBox(5);
        msDataTypeBox.getChildren().addAll(msDataTypeLabel, msDataTypeComboBox);

        paramGrid.add(threadsBox, 0, 0);
        paramGrid.add(workflowBox, 1, 0);
        paramGrid.add(msDataTypeBox, 0, 1, 2, 1); // Span both columns

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
        
        ObservableList<MetaPilotIOPanelFX.ExperimentRow> experiments = ioPanel.getIoPanel().getExperimentData();
        metadataTableData.clear();
        for (MetaPilotIOPanelFX.ExperimentRow exp : experiments) {
            metadataTableData.add(new MetaPilotDdaIgcMaxQuantPanel.MetadataRow(exp.getExperimentName()));
        }
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
        metadataTable.setColumnResizePolicy(TableView.UNCONSTRAINED_RESIZE_POLICY);
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
        sb.append("=== FragPipe Analysis Parameters ===\n\n");

        // Database settings
        sb.append("Database Settings:\n");
        sb.append("  Microbiome DB: ").append(ioPanel != null ? ioPanel.getMicrobiomeDBPath() : "").append("\n");
        sb.append("  Append Host DB: ").append(ioPanel != null ? ioPanel.isHostDBEnabled() : false).append("\n");
        if (ioPanel != null && ioPanel.isHostDBEnabled() && !ioPanel.getHostDBPath().isEmpty()) {
            sb.append("  Host Database: ").append(ioPanel.getHostDBPath()).append("\n");
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

        // Additional parameters
        sb.append("Additional Parameters:\n");
        if (threadsSpinner != null) {
            sb.append("  Threads: ").append(threadsSpinner.getValue()).append("\n");
        }
        if (workflowComboBox != null) {
            sb.append("  Workflow: ").append(workflowComboBox.getValue()).append("\n");
        }
        if (msDataTypeComboBox != null) {
            sb.append("  MS data type: ").append(msDataTypeComboBox.getValue()).append("\n");
        }
        if (generateSampleSpecificDBCheckBox != null) {
            sb.append("  Generate Sample-specific DB: ").append(generateSampleSpecificDBCheckBox.isSelected()).append("\n");
        }
        if (spectraClusteringCheckBox != null) {
            sb.append("  Spectra Clustering: ").append(spectraClusteringCheckBox.isSelected()).append("\n");
        }
        sb.append("\n");

        // Note: Quantitative mode is determined by the selected workflow
        sb.append("Quantitative Mode: Determined by workflow selection\n\n");

        // Taxonomic analysis
        sb.append("\n").append(taxonomicAnalysisPanel.getTaxonomicSummary());

        parameterView.setText(sb.toString());
    }

    private void updateMetadataTable() {
        if (metadataTable == null || ioPanel == null) {
            return;
        }
        
        metadataTableData.clear();
        ObservableList<MetaPilotIOPanelFX.ExperimentRow> fileRows = ioPanel.getIoPanel().getExperimentData();
        if (fileRows.isEmpty()) return;

        createMetadataTableColumns();

        // Simple metadata table for label-free mode (quantitative mode determined by workflow)
        for (MetaPilotIOPanelFX.ExperimentRow row : fileRows) {
            metadataTableData.add(new MetaPilotDdaIgcMaxQuantPanel.MetadataRow(row.getExperimentName()));
        }
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
        consoleOutput.appendText("Starting FragPipe analysis...\n");
        System.out.println("DEBUG: startAnalysis called - FragPipe workflow");
    }

    private void stopAnalysis() {
        startButton.setDisable(false);
        stopButton.setDisable(true);
        consoleOutput.appendText("Stopping FragPipe analysis...\n");
        System.out.println("DEBUG: stopAnalysis called");
    }
    
    private boolean validateInputs() {
        if (ioPanel.getIoPanel().getExperimentData().isEmpty()) {
            showAlert("Validation Error", "Please add at least one experiment.");
            return false;
        }
        if (ioPanel.getIoPanel().getOutputDirectory().isEmpty()) {
            showAlert("Validation Error", "Please select an output directory.");
            return false;
        }
        if (ioPanel.getMicrobiomeDBPath().isEmpty()) {
            showAlert("Validation Error", "Please select a microbiome database file.");
            return false;
        }
        if (ioPanel.isHostDBEnabled() && ioPanel.getHostDBPath().isEmpty()) {
            showAlert("Validation Error", "Please select a host database file.");
            return false;
        }
        if (enzymeComboBox.getValue() == null) {
            showAlert("Validation Error", "Please select an enzyme.");
            return false;
        }
        return true;
    }
    
    private void showAlert(String title, String message) {
        Alert alert = new Alert(AlertType.INFORMATION);
        alert.setTitle(title);
        alert.setHeaderText(null);
        alert.setContentText(message);
        alert.showAndWait();
    }
    
    private void showError(String title, String content) {
        Alert alert = new Alert(AlertType.ERROR);
        alert.setTitle(title);
        alert.setHeaderText(null);
        alert.setContentText(content);
        alert.showAndWait();
    }
    
} 