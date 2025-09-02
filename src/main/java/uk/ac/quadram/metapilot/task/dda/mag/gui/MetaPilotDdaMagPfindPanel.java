package uk.ac.quadram.metapilot.task.dda.mag.gui;

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
import uk.ac.quadram.metapilot.task.common.gui.MetaPilotIOPanelFX;
import uk.ac.quadram.metapilot.task.common.gui.ConsolePanel;
import uk.ac.quadram.metapilot.task.common.gui.MetadataConfigurationPanel;
import uk.ac.quadram.metapilot.task.common.gui.IsobaricLabelingPanel;
import uk.ac.quadram.metapilot.task.common.gui.ModificationsPanel;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;
import uk.ac.quadram.metapilot.task.mag.par.MetaParameterMag;
import uk.ac.quadram.metapilot.task.dda.igc.gui.MetaPilotDdaIgcMaxQuantPanel;
import uk.ac.quadram.metapilot.task.dda.mag.DDAMagPfindTask;
import uk.ac.quadram.metapilot.task.par.MetaBaseParameter;
import uk.ac.quadram.metapilot.task.par.MetaData;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.gui.JavaFXProgressBar;
import uk.ac.quadram.metapilot.task.common.gui.MetaPilotIOPanelFX.ExperimentRow;
import uk.ac.quadram.metapilot.dbSearch.MaxQuant.MaxquantModification;
import javafx.application.Platform;
import javafx.concurrent.Task;
import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.io.File;

public class MetaPilotDdaMagPfindPanel extends VBox {
    // UI Components - Database specific
    private MetaPilotIOPanelFX ioPanel;
    private ComboBox<MagDbItem> magDatabaseComboBox;
    private TextField hostDBField;
    private CheckBox appendHostDBCheckBox;
    private Button browseHostDBButton;
    private ToggleGroup quantModeGroup;

    // Reusable UI Components
    private ConsolePanel consolePanel;
    private MetadataConfigurationPanel metadataPanel;
    private IsobaricLabelingPanel isobaricPanel;
    private ModificationsPanel modificationsPanel;

    // Parameter Panel Components
    private ComboBox<String> enzymeComboBox;
    private Spinner<Integer> missedCleavagesSpinner;
    private ComboBox<String> digestionModeComboBox;

    // Additional parameters
    private Spinner<Integer> threadsSpinner;
    private ComboBox<String> ms2ScanModeComboBox;
    private CheckBox generateSampleSpecificDBCheckBox;
    private CheckBox spectraClusteringCheckBox;

    // Data
    private ObservableList<String> enzymes = FXCollections.observableArrayList();

    // MAG database items
    private ObservableList<MagDbItem> availableMagDatabases = FXCollections.observableArrayList();

    public MetaPilotDdaMagPfindPanel(MagDbItem[] magDbItems) {
        System.out.println("DEBUG: MetaPilotDdaMagPfindPanel() constructor started");
        try {
            loadEnzymesFromIni();
            // Initialize MAG databases from parameter
            System.out.println("DEBUG: Loading MAG databases...");
            System.out.println("DEBUG: Found " + (magDbItems != null ? magDbItems.length : 0) + " MAG database items");
            initializeUI();
            // Now set MAG databases after UI is initialized
            setAvailableMagDatabases(magDbItems);
            System.out.println("DEBUG: MetaPilotDdaMagPfindPanel(MetaParameterMag) constructor completed successfully");
        } catch (Exception e) {
            System.err.println("ERROR: Exception in MetaPilotDdaMagPfindPanel constructor: " + e.getMessage());
            e.printStackTrace();
            throw e;
        }
    }

    public MetaPilotDdaMagPfindPanel(MetaParameterMag metaParameterMag) {
        System.out.println("DEBUG: MetaPilotDdaMagPfindPanel(MetaParameterMag) constructor started");
        try {
            loadEnzymesFromIni();
            // Initialize MAG databases from parameter
            if (metaParameterMag != null) {
                System.out.println("DEBUG: MetaParameterMag provided, loading MAG databases...");
                MagDbItem[] magDbItems = metaParameterMag.getAvailableMagDbItem();
                System.out.println("DEBUG: Found " + (magDbItems != null ? magDbItems.length : 0) + " MAG database items");
            } else {
                System.out.println("DEBUG: MetaParameterMag is null, no MAG databases to load");
            }
            initializeUI();
            // Now set MAG databases after UI is initialized
            if (metaParameterMag != null) {
                MagDbItem[] magDbItems = metaParameterMag.getAvailableMagDbItem();
                setAvailableMagDatabases(magDbItems);
            }
            System.out.println("DEBUG: MetaPilotDdaMagPfindPanel(MetaParameterMag) constructor completed successfully");
        } catch (Exception e) {
            System.err.println("ERROR: Exception in MetaPilotDdaMagPfindPanel(MetaParameterMag) constructor: " + e.getMessage());
            e.printStackTrace();
            throw e;
        }
    }

    // Test method to create sample MAG database items for testing
    public static MagDbItem[] createTestMagDbItems() {
        List<MagDbItem> testItems = new ArrayList<>();
        
        // Create test MAG database items similar to the real ones
        MagDbItem chickenGut = new MagDbItem("chicken-gut", "Chicken Gut v1.0");
        testItems.add(chickenGut);
        
        MagDbItem humanGut = new MagDbItem("human-gut", "Human Gut v2.0.1");
        testItems.add(humanGut);
        
        MagDbItem mouseGut = new MagDbItem("mouse-gut", "Mouse Gut v1.0");
        testItems.add(mouseGut);
        
        MagDbItem pigGut = new MagDbItem("pig-gut", "Pig Gut v1.0");
        testItems.add(pigGut);
        
        return testItems.toArray(new MagDbItem[0]);
    }

    private void initializeUI() {
        System.out.println("DEBUG: initializeUI() started");
        try {
            System.out.println("DEBUG: Setting basic panel properties...");
            setSpacing(10);
            setPadding(new Insets(10));
            setMinWidth(1200);
            setMinHeight(800);
            setPrefWidth(1600);
            setPrefHeight(900);

            System.out.println("DEBUG: Creating reusable UI components...");
            // Initialize reusable components
            consolePanel = new ConsolePanel();
            metadataPanel = new MetadataConfigurationPanel();
            isobaricPanel = new IsobaricLabelingPanel();
            modificationsPanel = new ModificationsPanel();
            
            // Set up callbacks for parameter updates
            metadataPanel.setParameterUpdateCallback(this::updateParameterView);
            isobaricPanel.setParameterUpdateCallback(() -> {
                updateParameterView();
                updateMetadataTable(); // Also update metadata when isobaric settings change
            });
            modificationsPanel.setParameterUpdateCallback(this::updateParameterView);
            consolePanel.setStartCallback(this::startAnalysis);
            consolePanel.setStopCallback(this::stopAnalysis);

            System.out.println("DEBUG: Creating tab pane...");
            // Create tab pane with responsive sizing
            TabPane tabPane = new TabPane();
            tabPane.setTabClosingPolicy(TabPane.TabClosingPolicy.UNAVAILABLE);
            tabPane.setMinWidth(980);
            tabPane.setMinHeight(780);
            VBox.setVgrow(tabPane, Priority.ALWAYS);

            System.out.println("DEBUG: Creating IO tab...");
            // Create tabs
            Tab ioTab = new Tab("Input/Output");
            ioTab.setContent(createIOTab());

            System.out.println("DEBUG: Creating Parameters tab...");
            Tab paramTab = new Tab("Parameters");
            paramTab.setContent(createParametersTab());

            System.out.println("DEBUG: Creating Console tab...");
            Tab consoleTab = new Tab("Console");
            consoleTab.setContent(consolePanel);

            System.out.println("DEBUG: Adding tabs to tab pane...");
            tabPane.getTabs().addAll(ioTab, paramTab, consoleTab);

            System.out.println("DEBUG: Adding tab change listener...");
            // Add tab change listener to update metadata when Parameters tab is selected
            tabPane.getSelectionModel().selectedItemProperty().addListener((obs, oldTab, newTab) -> {
                if (newTab == paramTab) {
                    updateMetadataTableFromIOPanel();
                }
            });

            System.out.println("DEBUG: Adding tab pane to main panel...");
            // Add tab pane to main panel
            getChildren().add(tabPane);
            
            System.out.println("DEBUG: initializeUI() completed successfully");
        } catch (Exception e) {
            System.err.println("ERROR: Exception in initializeUI(): " + e.getMessage());
            e.printStackTrace();
            throw e;
        }
    }

    private VBox createIOTab() {
        System.out.println("DEBUG: createIOTab() started");
        try {
            VBox tabContent = new VBox(10);
            tabContent.setPadding(new Insets(10));
            System.out.println("DEBUG: Creating MetaLabIOPanelFX...");
            ioPanel = new MetaPilotIOPanelFX();
            System.out.println("DEBUG: Creating database section...");
            VBox databaseSection = createDatabaseSection();
            System.out.println("DEBUG: Adding components to IO tab...");
            tabContent.getChildren().addAll(
                createSectionTitle("Input/Output Settings"),
                ioPanel,
                createSectionTitle("Database Settings"),
                databaseSection
            );
            System.out.println("DEBUG: createIOTab() completed successfully");
            return tabContent;
        } catch (Exception e) {
            System.err.println("ERROR: Exception in createIOTab(): " + e.getMessage());
            e.printStackTrace();
            throw e;
        }
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
        
        HBox appendHostBox = new HBox(10);
        appendHostDBCheckBox = new CheckBox("Append Host Database");
        appendHostDBCheckBox.setOnAction(e -> updateHostDBVisibility());
        appendHostBox.getChildren().add(appendHostDBCheckBox);
        
        HBox hostDBBox = new HBox(10);
        hostDBField = new TextField();
        hostDBField.setPromptText("Select host database file");
        hostDBField.setPrefWidth(400);
        hostDBField.setDisable(true);
        browseHostDBButton = new Button("Browse");
        browseHostDBButton.setDisable(true);
        browseHostDBButton.setOnAction(e -> browseHostDatabase());
        hostDBBox.getChildren().addAll(new Label("Host Database:"), hostDBField, browseHostDBButton);
        
        section.getChildren().addAll(magDBBox, appendHostBox, hostDBBox);
        return section;
    }

    private void updateHostDBVisibility() {
        boolean selected = appendHostDBCheckBox.isSelected();
        hostDBField.setDisable(!selected);
        browseHostDBButton.setDisable(!selected);
    }

    private void browseHostDatabase() {
        FileChooser chooser = new FileChooser();
        chooser.setTitle("Select Host Database File");
        chooser.getExtensionFilters().add(
            new FileChooser.ExtensionFilter("FASTA Files", "*.fasta", "*.fa")
        );
        File file = chooser.showOpenDialog(getScene().getWindow());
        if (file != null) {
            hostDBField.setText(file.getAbsolutePath());
        }
    }

    // Load available MAG databases
    public void setAvailableMagDatabases(MagDbItem[] magDbItems) {
        System.out.println("DEBUG: setAvailableMagDatabases called with " + (magDbItems != null ? magDbItems.length : 0) + " items");
        availableMagDatabases.clear();
        if (magDbItems != null) {
            for (MagDbItem item : magDbItems) {
                System.out.println("DEBUG: Adding MAG database: " + item.getCatalogueID() + " - " + item.getIdentifier());
                availableMagDatabases.add(item);
            }
        }
        System.out.println("DEBUG: Available MAG databases size: " + availableMagDatabases.size());
        // Add null check to prevent NPE
        if (magDatabaseComboBox != null && !availableMagDatabases.isEmpty()) {
            magDatabaseComboBox.setValue(availableMagDatabases.get(0));
            System.out.println("DEBUG: Set default MAG database to: " + availableMagDatabases.get(0).getCatalogueID());
        } else if (magDatabaseComboBox == null) {
            System.out.println("DEBUG: magDatabaseComboBox is null - UI not initialized yet");
        } else {
            System.out.println("DEBUG: No MAG databases available to set as default");
        }
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
        
        // Left column: Enzyme Settings, Additional Parameters, Modifications (stacked) - NO taxonomic analysis for MAG
        VBox leftColumn = new VBox(15);
        leftColumn.setMinWidth(430);
        VBox enzymeSection = createEnzymeSection();
        VBox additionalParamsSection = createAdditionalParametersSection();
        
        leftColumn.getChildren().addAll(enzymeSection, additionalParamsSection, modificationsPanel);
        GridPane.setVgrow(leftColumn, Priority.ALWAYS);
        
        // Right column: Quantitative Mode and Metadata Configuration (stacked)
        VBox rightColumn = new VBox(15);
        rightColumn.setMinWidth(430);
        VBox quantModeSection = createQuantitativeModeSection();
        
        rightColumn.getChildren().addAll(quantModeSection, metadataPanel);
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

        section.getChildren().addAll(titleLabel, paramGrid);
        return section;
    }
    
    private void startAnalysis() {
        try {
            consolePanel.appendConsoleText("Starting MAG pFind analysis...\n");
            System.out.println("DEBUG: startAnalysis called - MAG pFind workflow");
            
            // 1. Validate inputs
            ObservableList<ExperimentRow> fileRows = ioPanel.getExperimentData();
            if (fileRows.isEmpty()) {
                consolePanel.appendConsoleText("ERROR: No files selected for analysis.\n");
                return;
            }
            
            MagDbItem selectedMagDb = magDatabaseComboBox.getValue();
            if (selectedMagDb == null) {
                consolePanel.appendConsoleText("ERROR: No MAG database selected.\n");
                return;
            }
            
            String outputDir = ioPanel.getOutputDirectory();
            if (outputDir == null || outputDir.isEmpty()) {
                consolePanel.appendConsoleText("ERROR: No output directory specified.\n");
                return;
            }
            
            // 2. Collect file and experiment information
            String[] rawFiles = new String[fileRows.size()];
            String[] expNames = new String[fileRows.size()];
            int[] fractions = new int[fileRows.size()];
            int[] replicates = new int[fileRows.size()];
            
            for (int i = 0; i < fileRows.size(); i++) {
                ExperimentRow row = fileRows.get(i);
                rawFiles[i] = row.getFilePath();
                expNames[i] = row.getExperimentName();
                fractions[i] = 1; // Default fraction
                replicates[i] = 1; // Default replicate
            }
            
            // 3. Collect isobaric labeling information
            String[] labelTitles = new String[0];
            if (quantModeGroup.getSelectedToggle() != null) {
                RadioButton selected = (RadioButton) quantModeGroup.getSelectedToggle();
                if (selected.getText().contains("Isobaric") && isobaricPanel.hasLabels()) {
                    ObservableList<IsobaricLabelingPanel.IsobaricLabelRow> isoLabels = isobaricPanel.getIsobaricTableData();
                    labelTitles = new String[isoLabels.size()];
                    for (int i = 0; i < isoLabels.size(); i++) {
                        labelTitles[i] = isoLabels.get(i).getLabelName();
                    }
                }
            }
            
            // 4. Collect metadata information
            ObservableList<MetaPilotDdaIgcMaxQuantPanel.MetadataRow> metadataRows = metadataPanel.getMetadataTableData();
            int metaTypeCount = metadataPanel.getMetadataColumnCount();
            String[][] metaInfo = new String[metadataRows.size()][Math.max(3, metaTypeCount)]; // At least 3 columns
            for (int i = 0; i < metadataRows.size(); i++) {
                MetaPilotDdaIgcMaxQuantPanel.MetadataRow row = metadataRows.get(i);
                metaInfo[i][0] = row.getAge() != null ? row.getAge() : "";
                metaInfo[i][1] = row.getGender() != null ? row.getGender() : "";
                metaInfo[i][2] = row.getCondition() != null ? row.getCondition() : "";
                // Add more metadata fields as needed
            }
            
            // 5. Create MetaData object
            MetaData metaData = new MetaData(rawFiles, expNames, fractions, replicates, 
                                           new String[0]); // Empty label titles initially
            metaData.setLabelTitle(labelTitles);
            metaData.setMetaInfo(metaInfo);
            
            // 6. Collect parameters for MetaBaseParameter constructor
            String workflowType = "MAG pFind";
            String microDb = selectedMagDb.getCatalogueID();
            boolean appendHostDb = appendHostDBCheckBox.isSelected();
            String hostDB = appendHostDb ? hostDBField.getText() : "";
            String result = outputDir;
            int coreCount = threadsSpinner.getValue() != null ? threadsSpinner.getValue() : 4;
            boolean isMetaWorkflow = true;
            
            // 7. Create MetaBaseParameter with required constructor parameters
            MetaBaseParameter metaBasePar = new MetaBaseParameter(workflowType, metaData, result, 
                                                                  microDb, appendHostDb, hostDB, 
                                                                  "", coreCount, coreCount, isMetaWorkflow);
            
            // 8. Collect parameters for MetaParameterMag constructor
            String[] fixMods = new String[0];  // TODO: Get from modifications panel
            String[] variMods = new String[0]; // TODO: Get from modifications panel
            String enzyme = enzymeComboBox.getValue() != null ? enzymeComboBox.getValue() : "Trypsin";
            int missCleavages = missedCleavagesSpinner.getValue() != null ? missedCleavagesSpinner.getValue() : 2;
            int digestMode = 1; // Default digest mode
            String quanMode = "Label-Free";
            if (quantModeGroup.getSelectedToggle() != null) {
                RadioButton selected = (RadioButton) quantModeGroup.getSelectedToggle();
                quanMode = selected.getText();
            }
            
            // 9. Create MetaParameterMag with proper constructor
            String magDbVersion = "";
            if (selectedMagDb.getVersionList() != null && !selectedMagDb.getVersionList().isEmpty()) {
                // Get the latest version from the version list
                magDbVersion = selectedMagDb.getVersionList().get(0).getCatalogueVersion();
            }
            MetaParameterMag metaParameterMag = new MetaParameterMag(metaBasePar, fixMods, variMods, enzyme,
                    missCleavages, digestMode, quanMode, new MaxquantModification[0], null, false, 
                    new double[0][0], selectedMagDb.getCatalogueID(), magDbVersion);
            
            // 10. Create MetaResourceParameters (may need to be loaded from config)
            MetaResourceParameters mrp = new MetaResourceParameters();
            // TODO: Load MRP from configuration if needed
            
            // 11. Create progress bars
            JavaFXProgressBar overallBar = new JavaFXProgressBar(0, 100, 0);
            JavaFXProgressBar currentBar = new JavaFXProgressBar(0, 100, 0);
            
            // 12. Create and run DDAMagPfindTask
            DDAMagPfindTask task = new DDAMagPfindTask(metaBasePar, mrp, overallBar, currentBar);
            
            // 13. Set up task callbacks
            task.messageProperty().addListener((obs, oldMsg, newMsg) -> {
                Platform.runLater(() -> {
                    if (newMsg != null) {
                        consolePanel.appendConsoleText("Current Task: " + newMsg + "\n");
                    }
                });
            });
            
            task.setOnSucceeded(e -> {
                Platform.runLater(() -> {
                    // Check if the task actually succeeded by examining the result
                    try {
                        Boolean taskResult = task.get();
                        if (taskResult != null && taskResult) {
                            consolePanel.appendConsoleText("Analysis completed successfully!\n");
                        } else {
                            consolePanel.appendConsoleText("Analysis completed but may not have been fully successful.\n");
                        }
                    } catch (Exception ex) {
                        consolePanel.appendConsoleText("Analysis failed: " + ex.getMessage() + "\n");
                    }
                });
            });
            
            task.setOnFailed(e -> {
                Platform.runLater(() -> {
                    Throwable exception = task.getException();
                    consolePanel.appendConsoleText("Analysis failed: " + 
                        (exception != null ? exception.getMessage() : "Unknown error") + "\n");
                });
            });
            
            // 14. Run task in background thread
            Thread taskThread = new Thread(task);
            taskThread.setDaemon(true);
            taskThread.start();
            
            consolePanel.appendConsoleText("Analysis task started successfully.\n");
            
        } catch (Exception e) {
            consolePanel.appendConsoleText("ERROR: Failed to start analysis: " + e.getMessage() + "\n");
            e.printStackTrace();
        }
    }

    private void stopAnalysis() {
        consolePanel.appendConsoleText("Stopping MAG pFind analysis...\n");
        System.out.println("DEBUG: stopAnalysis called");
    }

    private void updateMetadataTableFromIOPanel() {
        if (ioPanel == null || metadataPanel == null) {
            System.out.println("DEBUG: updateMetadataTableFromIOPanel called but some components are null - UI not initialized yet");
            return;
        }
        
        ObservableList<MetaPilotIOPanelFX.ExperimentRow> experiments = ioPanel.getExperimentData();
        metadataPanel.updateMetadataFromIOPanel(experiments);
    }

    private void updateMetadataTable() {
        if (metadataPanel == null || ioPanel == null) {
            System.out.println("DEBUG: updateMetadataTable called but metadata components are null - UI not initialized yet");
            return;
        }
        
        ObservableList<MetaPilotIOPanelFX.ExperimentRow> fileRows = ioPanel.getExperimentData();
        if (fileRows.isEmpty()) return;
        
        if (quantModeGroup != null && quantModeGroup.getSelectedToggle() != null) {
            RadioButton selected = (RadioButton) quantModeGroup.getSelectedToggle();
            String mode = selected.getText();
            boolean hasLabels = isobaricPanel.hasLabels();
            
            metadataPanel.updateMetadataForQuantMode(mode, hasLabels, fileRows, isobaricPanel.getIsobaricTableData());
        }
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

    private void updateParameterView() {
        if (consolePanel == null) {
            System.out.println("DEBUG: updateParameterView called but consolePanel is null - UI not initialized yet");
            return;
        }
        
        StringBuilder sb = new StringBuilder();
        sb.append("=== Analysis Parameters ===\n\n");
        
        // Database settings
        sb.append("Database Settings:\n");
        MagDbItem selectedDb = magDatabaseComboBox != null ? magDatabaseComboBox.getValue() : null;
        if (selectedDb != null) {
            sb.append("  MAG Database: ").append(selectedDb.getCatalogueID()).append("\n");
        }
        sb.append("  Append Host DB: ").append(appendHostDBCheckBox != null ? appendHostDBCheckBox.isSelected() : false).append("\n");
        if (appendHostDBCheckBox != null && appendHostDBCheckBox.isSelected() && hostDBField != null && !hostDBField.getText().isEmpty()) {
            sb.append("  Host Database: ").append(hostDBField.getText()).append("\n");
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
        
        // Additional parameters
        sb.append("Additional Parameters:\n");
        if (threadsSpinner != null) {
            sb.append("  Threads: ").append(threadsSpinner.getValue()).append("\n");
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
            
            if (selected.getText().contains("Isobaric") && isobaricPanel != null) {
                sb.append("  Isobaric Type: ").append(isobaricPanel.getSelectedIsobaricType()).append("\n");
            }
        }
        sb.append("\n");
        
        // Modifications
        sb.append("Modifications:\n");
        if (modificationsPanel != null) {
            sb.append("  Fixed: ").append(modificationsPanel.getFixedModifications().toString()).append("\n");
            sb.append("  Variable: ").append(modificationsPanel.getVariableModifications().toString()).append("\n");
        }
        sb.append("\n");
        
        consolePanel.setParameterText(sb.toString());
    }

    // Load test MAG databases for debugging
    private void loadTestMagDatabases() {
        System.out.println("DEBUG: Loading test MAG databases...");
        MagDbItem[] testItems = createTestMagDbItems();
        setAvailableMagDatabases(testItems);
    }
} 