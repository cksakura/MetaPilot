package uk.ac.quadram.metapilot.task.dda.igc.gui;

import javafx.application.Platform;
import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.control.*;
import javafx.scene.control.Toggle;
import javafx.scene.layout.*;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.scene.text.Font;
import javafx.scene.text.FontWeight;
import javafx.scene.paint.Color;
import javafx.scene.control.cell.PropertyValueFactory;
import javafx.scene.control.cell.TextFieldTableCell;
import javafx.util.converter.DefaultStringConverter;
import uk.ac.quadram.metapilot.task.common.gui.MetaPilotIOPanelFX;
import uk.ac.quadram.metapilot.task.common.gui.IsobaricLabelingPanel;
import uk.ac.quadram.metapilot.task.common.gui.DdaIgcInputOutputPanel;
import uk.ac.quadram.metapilot.task.par.MetaData;
import uk.ac.quadram.metapilot.MetaPilotApp;
import uk.ac.quadram.metapilot.task.dda.igc.DDAIgcPFindTask;
import uk.ac.quadram.metapilot.task.dda.igc.par.DDAIgcPFindParameters;
import uk.ac.quadram.metapilot.task.dda.igc.par.DDAIgcPFindParameterIO;
import uk.ac.quadram.metapilot.resources.MetaResourceIO;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.task.common.gui.AbstractWorkflowPanel;
import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class MetaPilotDdaIgcPFindPanel extends AbstractWorkflowPanel<DDAIgcPFindParameters> {
    private static final Logger LOGGER = LogManager.getLogger(MetaPilotDdaIgcPFindPanel.class);
    
    
    // UI Components (similar to MaxQuant panel)
    private DdaIgcInputOutputPanel ioPanel;
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
    private ComboBox<Integer> threadsComboBox;
    private ComboBox<String> ms2ScanModeComboBox;
    private CheckBox generateSampleSpecificDBCheckBox;
    private CheckBox spectraClusteringCheckBox;

    // Data
    private ObservableList<String> enzymes = FXCollections.observableArrayList();
    private ObservableList<String> modifications = FXCollections.observableArrayList();
    private ObservableList<String> fixedModifications = FXCollections.observableArrayList();
    private ObservableList<String> variableModifications = FXCollections.observableArrayList();
    // Isobaric labeling panel (reusable component)
    private IsobaricLabelingPanel isobaricLabelingPanel;


    // Metadata configuration components
    private TableView<MetaPilotDdaIgcMaxQuantPanel.MetadataRow> metadataTable;
    private ObservableList<MetaPilotDdaIgcMaxQuantPanel.MetadataRow> metadataTableData = FXCollections.observableArrayList();
    private ComboBox<String> metadataColumnComboBox;
    private Button addMetadataColumnButton;
    private Button removeMetadataColumnButton;
    private Button resetMetadataButton;

    public MetaPilotDdaIgcPFindPanel() {
        // Call parent constructor first
        super();
        
        // Initialize parameter IO instance
        parameterIO = new DDAIgcPFindParameterIO();
        
        loadEnzymesFromIni();
        loadModificationsFromIni();
        
        // Initialize isobaric labeling panel
        isobaricLabelingPanel = new IsobaricLabelingPanel();
        isobaricLabelingPanel.setParameterUpdateCallback(this::updateParameterView);
        
        // Initialize unified IO panel
        ioPanel = new DdaIgcInputOutputPanel(this::updateParameterView);
        
        // Initialize specific UI components
        initializeSpecificUI();
        
        // Load parameters after UI is fully initialized
        Platform.runLater(() -> {
            loadLastUsedParameters();
        });
    }

    private void initializeSpecificUI() {
        // Create main tab pane using parent method
        TabPane tabPane = createMainTabPane();
        
        // Replace the default IO tab content with IGC-specific content
        Tab ioTab = tabPane.getTabs().get(0);
        ioTab.setContent(createIOTab());
        
        // Replace the default Parameters tab content with IGC-specific content  
        Tab paramTab = tabPane.getTabs().get(1);
        paramTab.setContent(createParametersTab());
        
        // Add tab pane to main panel
        getChildren().add(tabPane);
    }
    
    // Override parent method to provide IGC-specific behavior
    @Override
    protected void onParametersTabSelected() {
        updateMetadataTableFromIOPanel();
    }

    // --- BEGIN: Methods adapted from MetaPilotDdaIgcMaxQuantPanel ---
    @Override
    public VBox createIOTab() {
        VBox tabContent = new VBox(10);
        tabContent.setPadding(new Insets(10));
        
        // Add unified IO Panel (includes both file selection and database settings)
        tabContent.getChildren().add(ioPanel);
        
        return tabContent;
    }
    @Override
    public VBox createParametersTab() {
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
        
        // Use the reusable isobaric labeling panel
        isobaricLabelingPanel.setVisible(false);
        isobaricLabelingPanel.setManaged(false);
        labelingStackPane.getChildren().add(isobaricLabelingPanel);
        labelFreeRadio.setOnAction(e -> {
            isobaricLabelingPanel.setVisible(false);
            isobaricLabelingPanel.setManaged(false);
            updateParameterView();
            updateMetadataTable();
        });
        isobaricRadio.setOnAction(e -> {
            isobaricLabelingPanel.setVisible(true);
            isobaricLabelingPanel.setManaged(true);
            updateParameterView();
            updateMetadataTable();
        });
        section.getChildren().addAll(titleLabel, radioBox, labelingStackPane);
        return section;
    }
    private void updateMetadataTableFromIOPanel() {
        // Update metadata table when Parameters tab is selected
        updateMetadataTable();
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
        
        resetMetadataButton = new Button("Reset All");
        resetMetadataButton.setMinWidth(120);
        resetMetadataButton.setPrefWidth(150);
        resetMetadataButton.setStyle("-fx-background-color: #ff6b6b; -fx-text-fill: white;");
        resetMetadataButton.setOnAction(e -> resetMetadataTable());
        
        controlBox.getChildren().addAll(metadataColumnComboBox, customColumnTextField, addMetadataColumnButton, removeMetadataColumnButton, resetMetadataButton);
        section.getChildren().addAll(titleLabel, infoLabel, metadataTable, controlBox);
        updateMetadataTable();
        return section;
    }
    @Override
    protected void updateParameterView() {
        StringBuilder sb = new StringBuilder();
        // Enzyme settings
        sb.append("=== ENZYME SETTINGS ===\n");
        sb.append("Enzyme: ").append(enzymeComboBox.getValue()).append("\n");
        sb.append("Missed Cleavages: ").append(missedCleavagesSpinner.getValue()).append("\n");
        sb.append("Digestion Mode: ").append(digestionModeComboBox.getValue()).append("\n\n");
        // Additional parameters
        sb.append("=== ADDITIONAL PARAMETERS ===\n");
        sb.append("Threads: ").append(threadsComboBox.getValue()).append("\n");
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
                sb.append("Isobaric Type: ").append(isobaricLabelingPanel.getSelectedIsobaricType() != null ? isobaricLabelingPanel.getSelectedIsobaricType() : "None").append("\n");
                sb.append("Configured Labels: ").append(isobaricLabelingPanel.getIsobaricTableData().size()).append("\n");
            }
        }
        // Taxonomic analysis
        
        // Update the console panel's parameter view
        if (consolePanel != null) {
            consolePanel.setParameterText(sb.toString());
        }
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
                hasLabels = isobaricLabelingPanel.hasLabels();
            }
            if (mode.contains("Isobaric") && hasLabels) {
                for (MetaPilotIOPanelFX.ExperimentRow row : fileRows) {
                    for (IsobaricLabelingPanel.IsobaricLabelRow isoLabel : isobaricLabelingPanel.getIsobaricTableData()) {
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

    private boolean validateInputs() {
        if (ioPanel.getIoPanel().getExperimentData().isEmpty()) {
            showAlert("Validation Error", "Please select at least one raw file for analysis.");
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

    @Override
    protected void startAnalysis() {
        try {
            consolePanel.appendConsoleText("Starting IGC pFind analysis...\n");
            System.out.println("DEBUG: startAnalysis called - IGC pFind workflow");
            
            // Validate inputs first
            if (!validateInputs()) {
                return;
            }
            
            // Collect parameters from UI
            consolePanel.appendConsoleText("Collecting parameters from UI...\n");
            DDAIgcPFindParameters analysisParams = collectParametersFromUI();
            if (analysisParams == null) {
                consolePanel.appendConsoleText("ERROR: Failed to collect parameters from UI.\n");
                LOGGER.error("Failed to collect parameters from UI in startAnalysis");
                return;
            }
            
            consolePanel.appendConsoleText("✓ Successfully collected parameters from UI\n");
            LOGGER.info("Successfully collected parameters from UI for IGC pFind workflow");
            
            // Store analysis parameters for viewReport
            this.currentAnalysisParams = analysisParams;
        
            // Save parameters to file
            String outputDir = analysisParams.getResult();
            if (outputDir != null && !outputDir.trim().isEmpty()) {
                consolePanel.appendConsoleText("Saving parameters to file...\n");
                boolean paramsSaved = parameterIO.saveAnalysisParameters(outputDir, analysisParams);
                if (!paramsSaved) {
                    consolePanel.appendConsoleText("WARNING: Failed to save parameters to file. Analysis will continue.\n");
                    LOGGER.warn("Failed to save analysis parameters to file: {}", outputDir);
                } else {
                    consolePanel.appendConsoleText("✓ Parameters saved to: " + outputDir + "/analysis_parameters_pfind.json\n");
                    LOGGER.info("Successfully saved analysis parameters to: {}", outputDir);
                }
                
                // Also save as last used parameters
                boolean lastParamsSaved = parameterIO.saveLastUsedParameters(analysisParams);
                if (lastParamsSaved) {
                    consolePanel.appendConsoleText("✓ Last used parameters saved\n");
                    LOGGER.info("Successfully saved last used parameters");
                } else {
                    consolePanel.appendConsoleText("WARNING: Failed to save last used parameters\n");
                    LOGGER.warn("Failed to save last used parameters");
                }
                
                // Save to preferences
                boolean prefsSaved = parameterIO.savePreferences(analysisParams);
                if (prefsSaved) {
                    consolePanel.appendConsoleText("✓ Preferences saved\n");
                    LOGGER.info("Successfully saved preferences");
                } else {
                    consolePanel.appendConsoleText("WARNING: Failed to save preferences\n");
                    LOGGER.warn("Failed to save preferences");
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
            
            // Create and run the DDAIgcPFindTask
            consolePanel.appendConsoleText("Creating IGC pFind task...\n");
            currentTask = new DDAIgcPFindTask(analysisParams, mrp, consolePanel.getOverallProgressBar(),
                    consolePanel.getCurrentTaskProgressBar());
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
            
            // Set up task completion callback to enable View Report button
            currentTask.setOnSucceeded(e -> {
                Platform.runLater(() -> {
                    // Check if the task actually succeeded by examining the result
                    try {
                        Boolean result = currentTask.get();
                        if (result != null && result) {
                            consolePanel.appendConsoleText("✓ Analysis completed successfully!\n");
                            consolePanel.enableViewReportButton();
                        } else {
                            consolePanel.appendConsoleText("⚠ Analysis completed but may not have been fully successful.\n");
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
            LOGGER.info("IGC pFind analysis task started successfully");
            
        } catch (Exception e) {
            String errorMsg = "ERROR: Failed to start analysis: " + e.getMessage();
            consolePanel.appendConsoleText(errorMsg + "\n");
            LOGGER.error("Failed to start IGC pFind analysis", e);
            e.printStackTrace();
        }
    }


    
    /**
     * Collect all parameters from the UI into a parameters object
     */
    @Override
    protected DDAIgcPFindParameters collectParametersFromUI() {
        try {
            consolePanel.appendConsoleText("Collecting parameters from UI...\n");
            
            // Collect required parameters for DDAIgcPFindParameters constructor
            String result = ioPanel.getIoPanel().getOutputDirectory();
            String microDb = ioPanel.getMicrobiomeDBPath();
            boolean appendHostDb = ioPanel.isHostDBEnabled();
            String hostDb = ioPanel.getHostDBPath();
            String enzyme = enzymeComboBox.getValue();
            int missCleavages = missedCleavagesSpinner.getValue();
            String digestMode = digestionModeComboBox.getValue();
            int threadCount = threadsComboBox.getValue() != null ? threadsComboBox.getValue() : 4;
            String ms2ScanMode = ms2ScanModeComboBox.getValue();
            boolean reduceDb = generateSampleSpecificDBCheckBox.isSelected();
            boolean spectraCluster = spectraClusteringCheckBox.isSelected();

            // Create MetaData object from experiment data
            ObservableList<MetaPilotIOPanelFX.ExperimentRow> fileRows = ioPanel.getIoPanel().getExperimentData();
            consolePanel.appendConsoleText("Found " + fileRows.size() + " experiment rows\n");
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
                LOGGER.info("collectParametersFromUI: Row {} - file: {}, name: {}, fraction: {}", i, rawFiles[i], expNames[i], fractions[i]);
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
        if (isobaricLabelingPanel.getSelectedIsobaricType() != null) {
            params.setIsobaricType(isobaricLabelingPanel.getSelectedIsobaricType());
        }
        params.setIsobaricTableData(new ArrayList<>(isobaricLabelingPanel.getIsobaricTableData()));
        
        // Taxonomic analysis
        
        // Metadata table data
        params.setMetadataTableData(new ArrayList<>(metadataTableData));
        
        consolePanel.appendConsoleText("✓ Successfully collected parameters from UI\n");
        LOGGER.info("Successfully collected parameters from UI for IGC pFind workflow");
        return params;
        
        } catch (Exception e) {
            String errorMsg = "ERROR: Failed to collect parameters from UI: " + e.getMessage();
            consolePanel.appendConsoleText(errorMsg + "\n");
            LOGGER.error("Failed to collect parameters from UI for IGC pFind workflow", e);
            return null;
        }
    }
    
    /**
     * Load last used parameters and apply them to the UI
     */
    @Override
    protected void loadLastUsedParameters() {
        try {
            consolePanel.appendConsoleText("Loading last used parameters...\n");
            DDAIgcPFindParameters params = (DDAIgcPFindParameters) parameterIO.loadLastUsedParameters();
            if (params != null) {
                applyParametersToUI(params);
                consolePanel.appendConsoleText("✓ Successfully loaded last used parameters\n");
                LOGGER.info("Loaded last used parameters for DDA IGC pFind workflow");
                
                // Also apply experiment data to the IO panel if available
                if (params.getExperimentData() != null && !params.getExperimentData().isEmpty()) {
                    ioPanel.getIoPanel().setExperimentData(params.getExperimentData());
                    consolePanel.appendConsoleText("Applied " + params.getExperimentData().size() + " experiment rows to IO panel\n");
                    LOGGER.info("Applied {} experiment rows to IO panel", params.getExperimentData().size());
                }
                
                // Apply output directory if available
                if (params.getResult() != null && !params.getResult().trim().isEmpty()) {
                    ioPanel.getIoPanel().setOutputDirectory(params.getResult());
                    consolePanel.appendConsoleText("Applied output directory: " + params.getResult() + "\n");
                    LOGGER.info("Applied output directory: {}", params.getResult());
                }
                
                // Update metadata table after loading parameters
                Platform.runLater(() -> {
                    updateMetadataTable();
                    updateParameterView();
                });
            } else {
                consolePanel.appendConsoleText("No previous parameters found, using defaults\n");
                LOGGER.warn("No parameters loaded, using defaults");
            }
        } catch (Exception e) {
            String errorMsg = "WARNING: Failed to load last used parameters: " + e.getMessage();
            consolePanel.appendConsoleText(errorMsg + "\n");
            LOGGER.warn("Failed to load last used pFind parameters: " + e.getMessage(), e);
            // Continue with default UI state
        }
    }
    
    /**
     * Apply parameters to the UI components
     */
    @Override
    protected void applyParametersToUI(DDAIgcPFindParameters params) {
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
        if (threadsComboBox != null) {
            threadsComboBox.setValue(params.getThreadCount());
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
        if (params.getIsobaricType() != null && isobaricLabelingPanel != null) {
            isobaricLabelingPanel.setSelectedIsobaricType(params.getIsobaricType());
        }
        
        // Apply isobaric table data if available
        if (params.getIsobaricTableData() != null && !params.getIsobaricTableData().isEmpty() && isobaricLabelingPanel != null) {
            isobaricLabelingPanel.setIsobaricTableData(params.getIsobaricTableData());
        }
        
        // Apply taxonomic settings
        
        // Update parameter view after applying parameters
        Platform.runLater(() -> updateParameterView());
        
        consolePanel.appendConsoleText("✓ Successfully applied parameters to UI\n");
        LOGGER.info("Successfully applied parameters to UI for IGC pFind workflow");
        
        } catch (Exception e) {
            String errorMsg = "ERROR: Failed to apply parameters to UI: " + e.getMessage();
            consolePanel.appendConsoleText(errorMsg + "\n");
            LOGGER.error("Failed to apply parameters to UI for IGC pFind workflow", e);
        }
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
     * Get the parameter class for this workflow
     */
    @Override
    protected Class<DDAIgcPFindParameters> getParameterClass() {
        return DDAIgcPFindParameters.class;
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
    
    
    @Override
    protected void stopAnalysis() {
        if (consolePanel != null) {
            consolePanel.appendConsoleText("Stopping pFind analysis...\n");
        }
        System.out.println("DEBUG: stopAnalysis called");
    }
    
    @Override
    protected String getWorkflowName() {
        return "DDA IGC pFind";
    }
} 