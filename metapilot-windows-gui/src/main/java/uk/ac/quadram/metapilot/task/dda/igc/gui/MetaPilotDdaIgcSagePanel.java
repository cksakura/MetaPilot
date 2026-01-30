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
import java.util.Optional;
import javafx.scene.control.cell.PropertyValueFactory;
import javafx.scene.control.cell.TextFieldTableCell;
import javafx.util.converter.DefaultStringConverter;
import javafx.beans.property.SimpleStringProperty;
import uk.ac.quadram.metapilot.task.common.gui.MetaPilotIOPanelFX;
import uk.ac.quadram.metapilot.task.common.gui.IsobaricLabelingPanel;
import uk.ac.quadram.metapilot.task.common.gui.DdaIgcInputOutputPanel;
import uk.ac.quadram.metapilot.task.common.gui.ModificationsPanel;
import uk.ac.quadram.metapilot.task.dda.igc.DDAIgcSageTask;
import uk.ac.quadram.metapilot.task.dda.igc.par.DDAIgcSageParameters;
import uk.ac.quadram.metapilot.task.dda.igc.par.DDAIgcSageParameterIO;
import uk.ac.quadram.metapilot.resources.MetaResourceIO;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.task.common.gui.AbstractWorkflowPanel;
import uk.ac.quadram.metapilot.task.par.MetaData;
import java.util.ArrayList;
import java.util.List;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import uk.ac.quadram.metapilot.core.mod.IsobaricTag;
import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.IOException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class MetaPilotDdaIgcSagePanel extends AbstractWorkflowPanel<DDAIgcSageParameters> {
    private static final Logger LOGGER = LogManager.getLogger(MetaPilotDdaIgcSagePanel.class);
    
    
    // UI Components (similar to pFind panel)
    private DdaIgcInputOutputPanel ioPanel;
    private ToggleGroup quantModeGroup;

    // Parameter Panel Components
    private ComboBox<String> enzymeComboBox;
    private Spinner<Integer> missedCleavagesSpinner;
    private ComboBox<String> digestionModeComboBox;

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

    // Additional parameters
    private ComboBox<Integer> threadsComboBox;
    private ComboBox<String> ms2ScanModeComboBox;
    private CheckBox generateSampleSpecificDBCheckBox;
    private CheckBox spectraClusteringCheckBox;

    // Data
    private ObservableList<String> enzymes = FXCollections.observableArrayList();
    private ObservableList<String> isobaricBaseTypes = FXCollections.observableArrayList();

    // Isobaric labeling components
    private ComboBox<String> isobaricTypeComboBox;
    private TableView<IsobaricLabelingPanel.IsobaricLabelRow> isobaricTable;
    private ObservableList<IsobaricLabelingPanel.IsobaricLabelRow> isobaricTableData = FXCollections.observableArrayList();


    // Metadata configuration components
    private TableView<MetaPilotDdaIgcMaxQuantPanel.MetadataRow> metadataTable;
    private ObservableList<MetaPilotDdaIgcMaxQuantPanel.MetadataRow> metadataTableData = FXCollections.observableArrayList();
    private ComboBox<String> metadataColumnComboBox;
    private Button addMetadataColumnButton;
    private Button removeMetadataColumnButton;
    private Button resetMetadataButton;



    public MetaPilotDdaIgcSagePanel() {
        // Call parent constructor first
        super();
        
        // Initialize parameter IO instance
        parameterIO = new DDAIgcSageParameterIO();
        
        loadEnzymesFromIni();
        loadIsobaricBaseTypes();
        
        // Initialize unified IO panel
        ioPanel = new DdaIgcInputOutputPanel(this::updateParameterView);
        
        // Initialize specific UI components
        initializeSpecificUI();
        
        // Load parameters after UI is fully initialized
        Platform.runLater(() -> {
            loadLastUsedParameters();
            
            // Add window close event handler to save parameters
            setupWindowCloseHandler();
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

    // Create IO tab (similar to pFind)
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

        // Left column: Sage Parameters, Enzyme Settings, Additional Parameters, and Modifications Panel
        VBox leftColumn = new VBox(15);
        leftColumn.setMinWidth(430);
        VBox sageParamsSection = createSageParametersSection();
        VBox enzymeSection = createEnzymeSection();
        VBox additionalParamsSection = createAdditionalParametersSection();
        
        // Create horizontal modifications panel with fixed mods (narrow) and variable mods (wider)
        modificationsPanel = new ModificationsPanel(this::updateParameterView);
        
        leftColumn.getChildren().addAll(sageParamsSection, enzymeSection, additionalParamsSection, modificationsPanel);
        GridPane.setVgrow(leftColumn, Priority.ALWAYS);

        // Right column: Quantitative Mode and Metadata Configuration (keeping original layout)
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
        section.setPadding(new Insets(5));
        section.setStyle("-fx-background-color: #f0f8f0; -fx-background-radius: 3;");

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

        // Create a separate fragment tolerance spinner for open search to avoid component sharing issues
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

        if (quantModeGroup != null && quantModeGroup.getSelectedToggle() != null) {
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
    }

    @Override
    protected void updateParameterView() {
        if (consolePanel == null) {
            return;
        }
        
        StringBuilder sb = new StringBuilder();
        sb.append("=== Sage Analysis Parameters ===\n\n");

        // Database settings
        sb.append("Database Settings:\n");
        sb.append("  Microbiome DB: ").append(ioPanel != null ? ioPanel.getMicrobiomeDBPath() : "").append("\n");
        sb.append("  Append Host DB: ").append(ioPanel != null ? ioPanel.isHostDBEnabled() : false).append("\n");
        if (ioPanel != null && ioPanel.isHostDBEnabled() && !ioPanel.getHostDBPath().isEmpty()) {
            sb.append("  Host Database: ").append(ioPanel.getHostDBPath()).append("\n");
        }
        sb.append("\n");

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
                if (openFragmentToleranceSpinner != null) {
                    sb.append("  Fragment Tolerance: ").append(openFragmentToleranceSpinner.getValue())
                            .append(" ppm\n");
                }
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
        if (threadsComboBox != null) {
            sb.append("  Threads: ").append(threadsComboBox.getValue()).append("\n");
        }
        if (ms2ScanModeComboBox != null) {
            sb.append("  MS2 Scan Mode: ").append(ms2ScanModeComboBox.getValue()).append("\n");
        }
        if (generateSampleSpecificDBCheckBox != null) {
            sb.append("  Generate Sample-specific DB: ").append(generateSampleSpecificDBCheckBox.isSelected()).append("\n");
        }
        if (spectraClusteringCheckBox != null) {
            sb.append("  Spectra Clustering: ").append(spectraClusteringCheckBox.isSelected()).append("\n");
        }
        sb.append("\n");

        // Quantitative mode
        sb.append("Quantitative Mode:\n");
        if (quantModeGroup != null && quantModeGroup.getSelectedToggle() != null) {
            RadioButton selected = (RadioButton) quantModeGroup.getSelectedToggle();
            sb.append("  Mode: ").append(selected.getText()).append("\n");

            if (selected.getText().contains("Isobaric") && isobaricTypeComboBox != null) {
                sb.append("  Isobaric Type: ").append(isobaricTypeComboBox.getValue()).append("\n");
                sb.append("  Configured Labels: ").append(isobaricTableData.size()).append("\n");
            }
        }
        sb.append("\n");

        // Taxonomic analysis

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

        if (quantModeGroup != null && quantModeGroup.getSelectedToggle() != null) {
            RadioButton selected = (RadioButton) quantModeGroup.getSelectedToggle();
            String mode = selected.getText();
            boolean hasLabels = !isobaricTableData.isEmpty();

            if (mode.contains("Isobaric") && hasLabels) {
                for (MetaPilotIOPanelFX.ExperimentRow row : fileRows) {
                    for (IsobaricLabelingPanel.IsobaricLabelRow isoLabel : isobaricTableData) {
                        MetaPilotDdaIgcMaxQuantPanel.MetadataRow metadataRow = 
                            new MetaPilotDdaIgcMaxQuantPanel.MetadataRow(row.getExperimentName());
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

    private boolean validateInputs() {
        LOGGER.info("Starting input validation for IGC Sage workflow");
        
        if (ioPanel.getIoPanel().getExperimentData().isEmpty()) {
            LOGGER.error("No experiment data found in IO panel");
            showError("Validation Error", "Please select at least one raw file for analysis.");
            return false;
        }
        LOGGER.info("Experiment data validation passed: {} experiments found", ioPanel.getIoPanel().getExperimentData().size());
        
        if (ioPanel.getIoPanel().getOutputDirectory().isEmpty()) {
            LOGGER.error("Output directory not set");
            showError("Validation Error", "Please select an output directory.");
            return false;
        }
        LOGGER.info("Output directory validation passed: {}", ioPanel.getIoPanel().getOutputDirectory());
        
        if (ioPanel.getMicrobiomeDBPath().isEmpty()) {
            LOGGER.error("Microbiome database not selected");
            showError("Validation Error", "Please select a microbiome database file.");
            return false;
        }
        LOGGER.info("Microbiome database validation passed: {}", ioPanel.getMicrobiomeDBPath());
        
        if (ioPanel.isHostDBEnabled() && ioPanel.getHostDBPath().isEmpty()) {
            LOGGER.error("Host database enabled but not selected");
            showError("Validation Error", "Please select a host database file.");
            return false;
        }
        if (enzymeComboBox.getValue() == null) {
            LOGGER.error("No enzyme selected");
            showError("Validation Error", "Please select an enzyme.");
            return false;
        }
        LOGGER.info("Enzyme validation passed: {}", enzymeComboBox.getValue());
        
        LOGGER.info("All validations passed for IGC Sage workflow");
        return true;
    }

    @Override
    protected void startAnalysis() {
        try {
            consolePanel.appendConsoleText("Starting IGC Sage analysis...\n");
            LOGGER.debug("[IGC Sage] startAnalysis called");
            
            // Validate inputs first
            if (!validateInputs()) {
                return;
            }
            
            // Collect parameters from UI
            consolePanel.appendConsoleText("Collecting parameters from UI...\n");
            DDAIgcSageParameters analysisParams = collectParametersFromUI();
            if (analysisParams == null) {
                consolePanel.appendConsoleText("ERROR: Failed to collect parameters from UI.\n");
                LOGGER.error("Failed to collect parameters from UI in startAnalysis");
                return;
            }
            
            consolePanel.appendConsoleText("✓ Successfully collected parameters from UI\n");
            LOGGER.info("Successfully collected parameters from UI for IGC Sage workflow");
            
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
                    consolePanel.appendConsoleText("✓ Parameters saved to: " + outputDir + "/analysis_parameters_sage.json\n");
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
            
            // Create and run the DDAIgcSageTask
            consolePanel.appendConsoleText("Creating IGC Sage task...\n");
            currentTask = new DDAIgcSageTask(analysisParams, mrp, consolePanel.getOverallProgressBar(),
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
            LOGGER.info("IGC Sage analysis task started successfully");
            
        } catch (Exception e) {
            String errorMsg = "ERROR: Failed to start analysis: " + e.getMessage();
            consolePanel.appendConsoleText(errorMsg + "\n");
            LOGGER.error("Failed to start IGC Sage analysis", e);
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
            if (consolePanel != null) {
                consolePanel.appendConsoleText("Stopping Sage analysis...\n");
            }
            LOGGER.info("stopAnalysis called - user confirmed stop");
            
            // Force stop the running task if it exists
            if (currentTask != null && currentTask.isRunning()) {
                consolePanel.appendConsoleText("Force stopping Sage task...\n");
                currentTask.forceStop();
                consolePanel.appendConsoleText("Sage task stop requested.\n");
            } else {
                consolePanel.appendConsoleText("No running Sage task to stop.\n");
            }
            
            // Unbind progress bars to prevent memory leaks
            if (consolePanel != null) {
                consolePanel.getOverallProgressBar().progressProperty().unbind();
                consolePanel.getCurrentTaskProgressBar().progressProperty().unbind();
                
                consolePanel.setButtonsEnabled(true, false, false);
                consolePanel.stopStatusTimer();
                consolePanel.getOverallProgressBar().setProgress(0);
                consolePanel.getCurrentTaskProgressBar().setProgress(0);
                consolePanel.setCurrentTaskLabel("Current Task: None");
            }
        } else {
            LOGGER.info("stopAnalysis cancelled by user");
        }
    }
    
    @Override
    protected String getWorkflowName() {
        return "DDA IGC Sage";
    }
    
    
    
    
    /**
     * Collect all parameters from the UI into a parameters object
     */
    @Override
    protected DDAIgcSageParameters collectParametersFromUI() {
        try {
            consolePanel.appendConsoleText("Collecting parameters from UI...\n");
            
            // Collect required parameters for DDAIgcSageParameters constructor
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

        // Create DDAIgcSageParameters using the correct constructor
        DDAIgcSageParameters params = new DDAIgcSageParameters(
            uk.ac.quadram.metapilot.workflow.WorkflowType.DDA_IGC_SAGE.getAnalysisType(),
            metaData,
            result,
            microDb,
            appendHostDb,
            hostDb,
            ms2ScanMode,
            threadCount
        );
        
        // Set additional parameters that aren't in the constructor - Convert ExperimentRow to String[]
        List<String[]> experimentDataArray = new ArrayList<>();
        for (MetaPilotIOPanelFX.ExperimentRow row : fileRows) {
            experimentDataArray.add(row.toArray());
        }
        params.setExperimentData(experimentDataArray);
        
        // Enzyme settings
        params.setEnzyme(enzyme);
        params.setMissCleavages(missCleavages);
        params.setDigestMode(digestMode);
        
        // Additional parameters
        params.setReduceDb(reduceDb);
        params.setSpectraCluster(spectraCluster);
        
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
                params.setPrecursorToleranceMin(-precursorTolerance);
                params.setPrecursorToleranceMax(precursorTolerance);

            } else if (searchType.equals("Open Search")) {
                // Open search: asymmetric tolerances
                double precursorMin = precursorToleranceMinSpinner != null
                        && precursorToleranceMinSpinner.getValue() != null ? precursorToleranceMinSpinner.getValue()
                                : -500.0;
                double precursorMax = precursorToleranceMaxSpinner != null
                        && precursorToleranceMaxSpinner.getValue() != null ? precursorToleranceMaxSpinner.getValue()
                                : 150.0;
                // Fragment tolerance is symmetric for both open and closed search
                fragmentTolerance = openFragmentToleranceSpinner != null && openFragmentToleranceSpinner.getValue() != null
                        ? openFragmentToleranceSpinner.getValue() : 20.0;
                precursorToleranceUnit = "Da";
                fragmentToleranceUnit = "ppm";

                params.setPrecursorToleranceMin(precursorMin);
                params.setPrecursorToleranceMax(precursorMax);
            }
        }

        // Don't set the old precursor_tolerance field - it's only for backward compatibility
        params.setFragmentTolerance(fragmentTolerance);
        params.setPrecursorToleranceUnit(precursorToleranceUnit);
        params.setFragmentToleranceUnit(fragmentToleranceUnit);
        
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
                    variableMods.add(row.getModName() + ":" + row.getSites() + ":" + row.getDeltaMass());
                }
            }
            params.setVariableModifications(variableMods);
        }
        
        // Quantitative mode
        if (quantModeGroup.getSelectedToggle() != null) {
            RadioButton selected = (RadioButton) quantModeGroup.getSelectedToggle();
            params.setQuantitativeMode(selected.getText());
        }
        
        // Isobaric labeling
        if (isobaricTypeComboBox.getValue() != null) {
            params.setIsobaricType(isobaricTypeComboBox.getValue());
        }
        // Convert JavaFX IsobaricLabelRow objects to String arrays
        List<String[]> isobaricDataArray = new ArrayList<>();
        for (IsobaricLabelingPanel.IsobaricLabelRow row : isobaricTableData) {
            isobaricDataArray.add(row.toArray());
        }
        params.setIsobaricTableData(isobaricDataArray);
        
        // Taxonomic analysis - these methods may not exist in DDAIgcSageParameters
        // They would need to be added to the parameter class if needed
        
        // Metadata table data - Convert JavaFX MetadataRow objects to String arrays
        List<String[]> metadataDataArray = new ArrayList<>();
        for (MetaPilotDdaIgcMaxQuantPanel.MetadataRow row : metadataTableData) {
            metadataDataArray.add(row.toArray());
        }
        params.setMetadataTableData(metadataDataArray);
        
        consolePanel.appendConsoleText("✓ Successfully collected parameters from UI\n");
        LOGGER.info("Successfully collected parameters from UI for IGC Sage workflow");
        return params;
        
        } catch (Exception e) {
            String errorMsg = "ERROR: Failed to collect parameters from UI: " + e.getMessage();
            consolePanel.appendConsoleText(errorMsg + "\n");
            LOGGER.error("Failed to collect parameters from UI for IGC Sage workflow", e);
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
            DDAIgcSageParameters params = (DDAIgcSageParameters) parameterIO.loadLastUsedParameters();
            if (params != null) {
                applyParametersToUI(params);
                consolePanel.appendConsoleText("✓ Successfully loaded last used parameters\n");
                LOGGER.info("Loaded last used parameters for DDA IGC Sage workflow");
                
                // Also apply experiment data to the IO panel if available - Convert String[] to ExperimentRow
                if (params.getExperimentData() != null && !params.getExperimentData().isEmpty()) {
                    List<MetaPilotIOPanelFX.ExperimentRow> experimentRows = new ArrayList<>();
                    for (String[] data : params.getExperimentData()) {
                        experimentRows.add(MetaPilotIOPanelFX.ExperimentRow.fromArray(data));
                    }
                    ioPanel.getIoPanel().setExperimentData(experimentRows);
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
                consolePanel.appendConsoleText("No previous parameters found, creating defaults\n");
                LOGGER.info("No parameters loaded from file, creating defaults");
                DDAIgcSageParameters defaultParams = new DDAIgcSageParameters();
                applyParametersToUI(defaultParams);
                consolePanel.appendConsoleText("✓ Default parameters loaded into UI\n");
                LOGGER.info("Applied default parameters to UI successfully");
            }
        } catch (Exception e) {
            String errorMsg = "WARNING: Failed to load last used parameters: " + e.getMessage();
            consolePanel.appendConsoleText(errorMsg + "\n");
            LOGGER.warn("Failed to load last used parameters: " + e.getMessage());
            // Continue with default UI state
        }
    }
    
    /**
     * Apply parameters to the UI components
     */
    @Override
    protected void applyParametersToUI(DDAIgcSageParameters params) {
        try {
            consolePanel.appendConsoleText("Applying parameters to UI...\n");
            if (params == null) {
                consolePanel.appendConsoleText("WARNING: No parameters to apply\n");
                LOGGER.warn("No parameters to apply to UI");
                return;
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
        
        // Apply Sage-specific parameters
        if (params.getPrecursorToleranceMin() != 0.0 || params.getPrecursorToleranceMax() != 0.0) {
            // Determine if this is open search (asymmetric) or closed search (symmetric)
            double min = params.getPrecursorToleranceMin();
            double max = params.getPrecursorToleranceMax();
            
            if (Math.abs(min + max) < 0.1) {
                // Symmetric (closed search) - min and max are opposite values
                double tolerance = Math.abs(max);
                if (searchTypeGroup != null) {
                    // Select closed search radio button
                    RadioButton closedSearchRadio = (RadioButton) searchTypeGroup.getToggles().get(0);
                    closedSearchRadio.setSelected(true);
                }
                if (precursorToleranceSpinner != null) {
                    precursorToleranceSpinner.getValueFactory().setValue(tolerance);
                }
            } else {
                // Asymmetric (open search)
                if (searchTypeGroup != null) {
                    // Select open search radio button
                    RadioButton openSearchRadio = (RadioButton) searchTypeGroup.getToggles().get(1);
                    openSearchRadio.setSelected(true);
                }
                if (precursorToleranceMinSpinner != null) {
                    precursorToleranceMinSpinner.getValueFactory().setValue(min);
                }
                if (precursorToleranceMaxSpinner != null) {
                    precursorToleranceMaxSpinner.getValueFactory().setValue(max);
                }
            }
        }
        
        // Apply fragment tolerance
        if (params.getFragmentTolerance() != 0.0) {
            if (fragmentToleranceSpinner != null) {
                fragmentToleranceSpinner.getValueFactory().setValue(params.getFragmentTolerance());
            }
            // Sync the open search fragment tolerance spinner
            if (openFragmentToleranceSpinner != null) {
                openFragmentToleranceSpinner.getValueFactory().setValue(params.getFragmentTolerance());
            }
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
        
        // Apply modifications - parse string lists back to modification rows
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
                            for (ModificationsPanel.FixedModificationRow row : modificationsPanel.getFixedModData()) {
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
                            for (ModificationsPanel.VariableModificationRow row : modificationsPanel.getVariableModData()) {
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
        
        // Apply isobaric table data if available - Convert String arrays back to JavaFX IsobaricLabelRow objects
        if (params.getIsobaricTableData() != null && !params.getIsobaricTableData().isEmpty()) {
            isobaricTableData.clear();
            for (String[] data : params.getIsobaricTableData()) {
                isobaricTableData.add(IsobaricLabelingPanel.IsobaricLabelRow.fromArray(data));
            }
        }
        
        // Apply taxonomic settings - these methods may not exist in DDAIgcSageParameters
        // They would need to be added to the parameter class if needed
        
        // Update parameter view after applying parameters
        Platform.runLater(() -> updateParameterView());
        
        consolePanel.appendConsoleText("✓ Successfully applied parameters to UI\n");
        LOGGER.info("Successfully applied parameters to UI for IGC Sage workflow");
        
        } catch (Exception e) {
            String errorMsg = "ERROR: Failed to apply parameters to UI: " + e.getMessage();
            consolePanel.appendConsoleText(errorMsg + "\n");
            LOGGER.error("Failed to apply parameters to UI for IGC Sage workflow", e);
        }
    }
    
    /**
     * Load parameters from a specific analysis result directory
     */
    public void loadAnalysisParameters(String resultDirectory) {
        try {
            DDAIgcSageParameters params = (DDAIgcSageParameters) parameterIO.loadAnalysisParameters(resultDirectory);
            applyParametersToUI(params);
            LOGGER.info("Loaded Sage analysis parameters from {}", resultDirectory);
            consolePanel.appendConsoleText("✓ Loaded Sage analysis parameters from " + resultDirectory + "\n");
        } catch (Exception e) {
            LOGGER.error("Failed to load Sage analysis parameters from " + resultDirectory, e);
            consolePanel.appendConsoleText("WARNING: Failed to load Sage analysis parameters: " + e.getMessage() + "\n");
        }
    }
    
    /**
     * Get the parameter class for this workflow
     */
    @Override
    protected Class<DDAIgcSageParameters> getParameterClass() {
        return DDAIgcSageParameters.class;
    }
    
    /**
     * Set up window close handler to save parameters when the window is closed
     */
} 