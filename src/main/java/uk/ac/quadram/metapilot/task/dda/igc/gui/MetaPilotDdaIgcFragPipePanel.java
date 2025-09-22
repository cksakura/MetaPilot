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
import javafx.scene.control.Alert.AlertType;
import javafx.scene.control.cell.PropertyValueFactory;
import javafx.scene.control.cell.TextFieldTableCell;
import javafx.beans.property.SimpleStringProperty;
import uk.ac.quadram.metapilot.task.common.gui.MetaPilotIOPanelFX;
import uk.ac.quadram.metapilot.task.common.gui.TaxonomicAnalysisPanel;
import uk.ac.quadram.metapilot.task.common.gui.DdaIgcInputOutputPanel;
import uk.ac.quadram.metapilot.task.common.gui.ModificationsPanel;
import uk.ac.quadram.metapilot.task.common.gui.ConsolePanel;
import uk.ac.quadram.metapilot.task.dda.igc.DDAIgcFragPipeTask;
import uk.ac.quadram.metapilot.task.dda.igc.par.DDAIgcFragPipeParameters;
import uk.ac.quadram.metapilot.task.dda.igc.par.DDAIgcFragPipeParameterIO;
import uk.ac.quadram.metapilot.task.par.MetaData;
import uk.ac.quadram.metapilot.MetaPilotApp;
import java.util.ArrayList;
import java.util.List;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.IOException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.io.File;

public class MetaPilotDdaIgcFragPipePanel extends VBox {
    private static final Logger LOGGER = LogManager.getLogger(MetaPilotDdaIgcFragPipePanel.class);
    
    // UI Components (similar to pFind panel)
    private DdaIgcInputOutputPanel ioPanel;
    private ConsolePanel consolePanel;

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
    private Button resetMetadataButton;

    // Add parameter IO instance
    private DDAIgcFragPipeParameterIO parameterIO;


    public MetaPilotDdaIgcFragPipePanel() {
        // Initialize parameter IO instance
        parameterIO = new DDAIgcFragPipeParameterIO();
        
        loadEnzymesFromIni();
        loadWorkflowsFromFragPipe();
        
        // Initialize unified IO panel
        ioPanel = new DdaIgcInputOutputPanel(this::updateParameterView);
        
        // Initialize taxonomic analysis panel
        taxonomicAnalysisPanel = new TaxonomicAnalysisPanel(this::updateParameterView);
        
        // Initialize console panel
        consolePanel = new ConsolePanel();
        consolePanel.setStartCallback(this::startAnalysis);
        consolePanel.setStopCallback(this::stopAnalysis);
        consolePanel.setViewReportCallback(this::viewReport);
        
        initializeUI();
        
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
                    DDAIgcFragPipeParameters currentParams = collectParametersFromUI();
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
			// Set default workflow to LFQ-MBR if it exists, otherwise use the first one
			if (workflows.contains("LFQ-MBR")) {
				workflowComboBox.setValue("LFQ-MBR");
				System.out.println("Set default workflow to: LFQ-MBR");
			} else if (!workflows.isEmpty()) {
				workflowComboBox.setValue(workflows.get(0));
				System.out.println("Set default workflow to: " + workflows.get(0));
			}
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
        
        resetMetadataButton = new Button("Reset All");
        resetMetadataButton.setMinWidth(120);
        resetMetadataButton.setPrefWidth(150);
        resetMetadataButton.setStyle("-fx-background-color: #ff6b6b; -fx-text-fill: white;");
        resetMetadataButton.setOnAction(e -> resetMetadataTable());

        controlBox.getChildren().addAll(
            metadataColumnComboBox, 
            customColumnTextField, 
            addMetadataColumnButton, 
            removeMetadataColumnButton,
            resetMetadataButton
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
        if (consolePanel == null) {
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

        // Update the console panel's parameter view
        if (consolePanel != null) {
            consolePanel.setParameterText(sb.toString());
        }
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

    private void startAnalysis() {
        // Clear console and show starting message
        if (consolePanel != null) {
            consolePanel.clearConsole();
            consolePanel.appendConsoleText("Starting FragPipe analysis...\n");
        }

        // Gather parameters from UI
        ObservableList<MetaPilotIOPanelFX.ExperimentRow> fileRows = ioPanel.getIoPanel().getExperimentData();
        if (fileRows.isEmpty()) {
            showAlert("Error", "No files selected for analysis.");
            return;
        }
        
        // Save parameters for this analysis
        DDAIgcFragPipeParameters analysisParams = collectParametersFromUI();
        
        String outputDir = ioPanel.getIoPanel().getOutputDirectory();
        if (outputDir != null && !outputDir.trim().isEmpty()) {
            // Save analysis-specific parameters
            boolean analysisParamsSaved = parameterIO.saveAnalysisParameters(outputDir, analysisParams);
            if (!analysisParamsSaved) {
                LOGGER.warn("Failed to save analysis parameters to {}/analysis_parameters_fragpipe.json", outputDir);
                showAlert("Warning", "Failed to save analysis parameters to result directory. Analysis will continue but parameters may not be preserved.");
            } else {
                LOGGER.info("Saved analysis parameters to {}/analysis_parameters_fragpipe.json", outputDir);
            }
            
            // Save as last used parameters
            parameterIO.saveLastUsedParameters(analysisParams);
            
            // Also save to preferences for quick defaults
            parameterIO.savePreferences(analysisParams);
        }

        // Load resource parameters
        uk.ac.quadram.metapilot.resources.MetaResourceParameters mrp = (new uk.ac.quadram.metapilot.resources.MetaResourceIO()).loadParameters();

        // Create and run the DDAIgcFragPipeTask
        DDAIgcFragPipeTask task = new DDAIgcFragPipeTask(analysisParams, mrp, consolePanel.getOverallProgressBar(),
                consolePanel.getCurrentTaskProgressBar());
        
        // Bind progress and messages
        if (consolePanel != null) {
            task.messageProperty().addListener((obs, oldMsg, newMsg) -> {
                javafx.application.Platform.runLater(() -> {
                    if (newMsg != null && !newMsg.isEmpty()) {
                        consolePanel.appendConsoleText(newMsg + "\n");
                        consolePanel.updateCurrentTaskFromMessage(newMsg);
                    }
                });
            });
        }
        
        // Start status timer in console panel
        if (consolePanel != null) {
            consolePanel.startStatusTimer();
        }
        
        // On succeeded/failed/cancelled
        task.setOnSucceeded(e -> {
            if (consolePanel != null) {
                consolePanel.setButtonsEnabled(true, false, true);
                consolePanel.stopStatusTimer();
            }
            
            // Check if the task actually succeeded by examining the result
            try {
                Boolean result = task.get();
                if (result != null && result) {
                    if (consolePanel != null) {
                        consolePanel.appendConsoleText("\nAnalysis completed successfully.\n");
                        consolePanel.enableViewReportButton();
                    }
                } else {
                    if (consolePanel != null) {
                        consolePanel.appendConsoleText("\nAnalysis completed but may not have been fully successful.\n");
                    }
                }
            } catch (Exception ex) {
                if (consolePanel != null) {
                    consolePanel.appendConsoleText("\nAnalysis failed: " + ex.getMessage() + "\n");
                }
            }
        });
        
        task.setOnFailed(e -> {
            if (consolePanel != null) {
                consolePanel.setButtonsEnabled(true, false, false);
                consolePanel.stopStatusTimer();
                consolePanel.appendConsoleText("\nAnalysis failed: " + task.getException() + "\n");
            }
        });
        
        task.setOnCancelled(e -> {
            if (consolePanel != null) {
                consolePanel.setButtonsEnabled(true, false, false);
                consolePanel.stopStatusTimer();
                consolePanel.appendConsoleText("\nAnalysis cancelled.\n");
            }
        });
        
        // Run the task
        Thread t = new Thread(task);
        t.setDaemon(true);
        t.start();
    }

    private void stopAnalysis() {
        if (consolePanel != null) {
            consolePanel.appendConsoleText("Stopping FragPipe analysis...\n");
        }
        System.out.println("DEBUG: stopAnalysis called");
    }
    
    private void viewReport() {
        if (consolePanel != null) {
            consolePanel.appendConsoleText("Opening report in default browser...\n");
        }
        
        // Get the output directory from the IO panel
        String outputDir = ioPanel.getIoPanel().getOutputDirectory();
        if (outputDir == null || outputDir.trim().isEmpty()) {
            if (consolePanel != null) {
                consolePanel.appendConsoleText("ERROR: No output directory specified.\n");
            }
            return;
        }
        
        java.io.File reportDir = new java.io.File(outputDir);
        if (!reportDir.exists()) {
            if (consolePanel != null) {
                consolePanel.appendConsoleText("ERROR: Report directory does not exist: " + outputDir + "\n");
            }
            return;
        }
        
        java.io.File indexFile = new java.io.File(reportDir, "index.html");
        if (!indexFile.exists()) {
            if (consolePanel != null) {
                consolePanel.appendConsoleText("ERROR: Report index.html not found in " + outputDir + "\n");
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
    
    /**
     * Collect all parameters from the UI into a parameters object
     */
    private DDAIgcFragPipeParameters collectParametersFromUI() {
        // Collect required parameters for DDAIgcFragPipeParameters constructor
        String result = ioPanel.getIoPanel().getOutputDirectory();
        String microDb = ioPanel.getMicrobiomeDBPath();
        boolean appendHostDb = ioPanel.isHostDBEnabled();
        String hostDb = ioPanel.getHostDBPath();
        String enzyme = enzymeComboBox.getValue();
        int missCleavages = missedCleavagesSpinner.getValue();
        String digestMode = digestionModeComboBox.getValue();
        int threadCount = threadsSpinner.getValue();
        String workflow = workflowComboBox.getValue();
        String msDataType = msDataTypeComboBox.getValue();
        boolean reduceDb = generateSampleSpecificDBCheckBox.isSelected();
        boolean spectraCluster = spectraClusteringCheckBox.isSelected();

        // Create MetaData object from experiment data
        ObservableList<MetaPilotIOPanelFX.ExperimentRow> fileRows = ioPanel.getIoPanel().getExperimentData();
        
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
        } else {
            metaData = new MetaData();
        }

        // Create DDAIgcFragPipeParameters using the correct constructor
        DDAIgcFragPipeParameters params = new DDAIgcFragPipeParameters(
            MetaPilotApp.WorkflowType.DDA_IGC_FRAGPIPE.getAnalysisType(),
            metaData,
            result,
            microDb,
            appendHostDb,
            hostDb,
            workflow,
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
        
        // Modifications - convert to String lists
        if (modificationsPanel != null) {
            List<String> fixedMods = new ArrayList<>();
            for (ModificationsPanel.FixedModificationRow row : modificationsPanel.getFixedModData()) {
                if (row.getDeltaMass() != 0.0) {
                    fixedMods.add(row.getSite() + ":" + row.getDeltaMass());
                }
            }
            params.setFixedModifications(fixedMods);
            
            List<String> variableMods = new ArrayList<>();
            for (ModificationsPanel.VariableModificationRow row : modificationsPanel.getVariableModData()) {
                if (row.getUse()) {
                    variableMods.add(row.getModName() + " (" + row.getSites() + "):" + row.getDeltaMass());
                }
            }
            params.setVariableModifications(variableMods);
        }
        
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
            DDAIgcFragPipeParameters params = (DDAIgcFragPipeParameters) parameterIO.loadLastUsedParameters();
            if (params != null) {
                applyParametersToUI(params);
                
                // Also apply experiment data to the IO panel if available
                if (params.getExperimentData() != null && !params.getExperimentData().isEmpty()) {
                    ioPanel.getIoPanel().setExperimentData(params.getExperimentData());
                }
                
                // Apply output directory if available
                if (params.getResult() != null && !params.getResult().trim().isEmpty()) {
                    ioPanel.getIoPanel().setOutputDirectory(params.getResult());
                }
                
                // Update metadata table after loading parameters
                Platform.runLater(() -> {
                    updateMetadataTable();
                    updateParameterView();
                });
            }
        } catch (Exception e) {
            // Continue with default UI state
        }
    }
    
    /**
     * Apply parameters to the UI components
     */
    private void applyParametersToUI(DDAIgcFragPipeParameters params) {
        if (params == null) return;
        
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
        // Note: Workflow and MS data type are FragPipe-specific and may not have getters
        // These would need to be added to DDAIgcFragPipeParameters if needed
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
        
        // Apply modifications - would need to parse string lists back to modification rows
        // This is complex and may not be needed for basic functionality
        
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
            DDAIgcFragPipeParameters params = (DDAIgcFragPipeParameters) parameterIO.loadAnalysisParameters(resultDirectory);
            applyParametersToUI(params);
            showAlert("Success", "Loaded FragPipe analysis parameters from " + resultDirectory);
            LOGGER.info("Loaded FragPipe analysis parameters from {}", resultDirectory);
        } catch (Exception e) {
            showAlert("Error", "Failed to load FragPipe analysis parameters: " + e.getMessage());
            LOGGER.error("Failed to load FragPipe analysis parameters from " + resultDirectory, e);
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
                                DDAIgcFragPipeParameters currentParams = collectParametersFromUI();
                                if (currentParams != null) {
                                    parameterIO.saveLastUsedParameters(currentParams);
                                    parameterIO.savePreferences(currentParams);
                                }
                            } catch (Exception e) {
                                // Ignore errors during window close
                            }
                        });
                    }
                }
            } catch (Exception e) {
                // Ignore errors during setup
            }
        });
    }
    
} 