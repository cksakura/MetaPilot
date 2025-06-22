package uk.ac.quadram.metapilot.task.dda.igc.gui;

import javafx.application.Platform;
import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.control.*;
import javafx.scene.layout.*;
import javafx.stage.DirectoryChooser;
import javafx.stage.FileChooser;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.scene.text.Font;
import javafx.scene.text.FontWeight;
import javafx.scene.paint.Color;
import javafx.scene.control.cell.PropertyValueFactory;
import javafx.scene.control.cell.TextFieldTableCell;
import javafx.util.converter.DefaultStringConverter;
import javafx.beans.property.SimpleStringProperty;
import javafx.beans.property.StringProperty;
import javafx.concurrent.Task;
import javafx.scene.control.Alert.AlertType;
import uk.ac.quadram.metapilot.task.common.gui.MetaLabIOPanelFX;
import uk.ac.quadram.metapilot.task.common.gui.MetaLabIOPanelFX.ExperimentRow;
import uk.ac.quadram.metapilot.task.par.MetaData;
import uk.ac.quadram.metapilot.task.par.MetaBaseParameter;
import uk.ac.quadram.metapilot.task.par.MetaConstants;
import uk.ac.quadram.metapilot.task.MetaLabWorkflowType;
import uk.ac.quadram.metapilot.task.v2.par.MetaParameterMQ;
import uk.ac.quadram.metapilot.task.v2.par.MetaMaxQuantPar;
import uk.ac.quadram.metapilot.task.v2.par.MetaOpenPar;
import uk.ac.quadram.metapilot.task.par.MetaPep2TaxaPar;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.dbSearch.MaxQuant.MaxquantModification;
import uk.ac.quadram.metapilot.core.mod.IsobaricTag;
import uk.ac.quadram.metapilot.task.dda.igc.DDAIgcMaxQuantTask;
import javafx.animation.KeyFrame;
import javafx.animation.Timeline;
import javafx.util.Duration;
import java.lang.management.ManagementFactory;
import java.lang.management.OperatingSystemMXBean;
import java.lang.management.MemoryMXBean;
import java.lang.management.MemoryUsage;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Optional;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class MetaPilotDdaIgcMaxQuantPanel extends VBox {
    // UI Components
    private MetaLabIOPanelFX ioPanel;
    private TextField microbiomeDBField;
    private TextField hostDBField;
    private CheckBox appendHostDBCheckBox;
    private Button browseMicrobiomeDBButton;
    private Button browseHostDBButton;
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
    private TextArea enzymeSpecificityArea;
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
    private ObservableList<String> enzymes;
    private ObservableList<String> modifications;
    private ObservableList<String> fixedModifications;
    private ObservableList<String> variableModifications;
    // Labeling
    private ObservableList<String> labelModifications = FXCollections.observableArrayList();
    private ObservableList<String> isobaricModifications = FXCollections.observableArrayList();
    private ObservableList<String> isobaricBaseTypes = FXCollections.observableArrayList();
    
    // Isotopic labeling components
    private ListView<String> lightLabelListView;
    private ListView<String> mediumLabelListView;
    private ListView<String> heavyLabelListView;
    private ComboBox<String> labelModificationComboBox;
    private Button addLightLabelButton;
    private Button addMediumLabelButton;
    private Button addHeavyLabelButton;
    private Button removeLightLabelButton;
    private Button removeMediumLabelButton;
    private Button removeHeavyLabelButton;
    private ObservableList<String> lightLabels = FXCollections.observableArrayList();
    private ObservableList<String> mediumLabels = FXCollections.observableArrayList();
    private ObservableList<String> heavyLabels = FXCollections.observableArrayList();
    
    // Isobaric labeling components
    private ComboBox<String> isobaricTypeComboBox;
    private TableView<IsobaricLabelRow> isobaricTable;
    private ObservableList<IsobaricLabelRow> isobaricTableData = FXCollections.observableArrayList();
    
    // Taxonomic analysis components
    private CheckBox eukaryotaCheckBox;
    private CheckBox bacteriaCheckBox;
    private CheckBox archaeaCheckBox;
    private CheckBox virusesCheckBox;
    private CheckBox viroidsCheckBox;
    private CheckBox unclassifiedCheckBox;
    private CheckBox otherCheckBox;
    private ComboBox<String> ignoreRankComboBox;
    private CheckBox environmentalSamplesCheckBox;
    
    // Metadata configuration components
    private TableView<MetadataRow> metadataTable;
    private ObservableList<MetadataRow> metadataTableData = FXCollections.observableArrayList();
    private ComboBox<String> metadataColumnComboBox;
    private Button addMetadataColumnButton;
    private Button removeMetadataColumnButton;
    
    // Add fields for status panel and timer
    private Label statusTimeLabel;
    private Label statusCpuLabel;
    private Label statusMemLabel;
    private Timeline statusTimeline;
    private long taskStartTimeMillis;
    
    public MetaPilotDdaIgcMaxQuantPanel() {
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
    
    private int ProcessorCount() {
		int total = 0;
		try {
			String[] args = { "wmic", "cpu", "get", "NumberOfLogicalProcessors" };
			ProcessBuilder pb = new ProcessBuilder(args);
			Process p = pb.start();
			BufferedInputStream in = new BufferedInputStream(p.getInputStream());
			BufferedReader inBr = new BufferedReader(new InputStreamReader(in));
			String lineStr = null;
			L: while ((lineStr = inBr.readLine()) != null) {
				String trimLine = lineStr.trim();
				if (trimLine.length() > 0) {
					int oneCore = 0;
					for (int i = 0; i < trimLine.length(); i++) {
						char aa = trimLine.charAt(i);
						if (aa >= '0' && aa <= '9') {
							oneCore += ((int) (aa - '0')) * Math.pow(10, (trimLine.length() - i - 1));
						} else {
							continue L;
						}
					}
					total += oneCore;
				}
			}
			if (p.waitFor() != 0) {
				if (p.exitValue() == 1)
					System.err.println("false");
			}
			inBr.close();
			in.close();
		} catch (IOException | InterruptedException e) {
			System.err.println("Failed to get logical processors using wmic: " + e.getMessage());
		}

		if (total == 0) {
			total = Runtime.getRuntime().availableProcessors();
		}
		return total;
	}
    
    private VBox createIOTab() {
        VBox tabContent = new VBox(10);
        tabContent.setPadding(new Insets(10));
        
        // Add IO Panel
        ioPanel = new MetaLabIOPanelFX();
        
        // Add Database Panel
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
        
        // Microbiome Database
        HBox microbiomeDBBox = new HBox(10);
        microbiomeDBField = new TextField();
        microbiomeDBField.setPromptText("Select microbiome database file");
        microbiomeDBField.setPrefWidth(400);
        browseMicrobiomeDBButton = new Button("Browse");
        browseMicrobiomeDBButton.setOnAction(e -> browseMicrobiomeDatabase());
        microbiomeDBBox.getChildren().addAll(new Label("Microbiome Database:"), microbiomeDBField, browseMicrobiomeDBButton);
        
        // Append Host Database
        HBox appendHostBox = new HBox(10);
        appendHostDBCheckBox = new CheckBox("Append Host Database");
        appendHostDBCheckBox.setOnAction(e -> updateHostDBVisibility());
        appendHostBox.getChildren().add(appendHostDBCheckBox);
        
        // Host Database
        HBox hostDBBox = new HBox(10);
        hostDBField = new TextField();
        hostDBField.setPromptText("Select host database file");
        hostDBField.setPrefWidth(400);
        hostDBField.setDisable(true);
        browseHostDBButton = new Button("Browse");
        browseHostDBButton.setDisable(true);
        browseHostDBButton.setOnAction(e -> browseHostDatabase());
        hostDBBox.getChildren().addAll(new Label("Host Database:"), hostDBField, browseHostDBButton);
        
        section.getChildren().addAll(microbiomeDBBox, appendHostBox, hostDBBox);
        return section;
    }
    
    private void updateHostDBVisibility() {
        boolean selected = appendHostDBCheckBox.isSelected();
        hostDBField.setDisable(!selected);
        browseHostDBButton.setDisable(!selected);
    }
    
    private void browseMicrobiomeDatabase() {
        FileChooser chooser = new FileChooser();
        chooser.setTitle("Select Microbiome Database File");
        chooser.getExtensionFilters().add(
            new FileChooser.ExtensionFilter("FASTA Files", "*.fasta", "*.fa")
        );
        File file = chooser.showOpenDialog(getScene().getWindow());
        if (file != null) {
            microbiomeDBField.setText(file.getAbsolutePath());
        }
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
    
    private VBox createParametersTab() {
        VBox tabContent = new VBox(10);
        tabContent.setPadding(new Insets(10));
        tabContent.setMinWidth(960);
        tabContent.setMinHeight(760);
        
        // Initialize data
        enzymes = FXCollections.observableArrayList();
        modifications = FXCollections.observableArrayList();
        fixedModifications = FXCollections.observableArrayList();
        variableModifications = FXCollections.observableArrayList();
        // Set default modifications
        fixedModifications.add("Carbamidomethyl (C)");
        variableModifications.addAll("Acetyl (Protein N-term)", "Oxidation (M)");
        
        // Load enzymes and modifications from XML
        loadEnzymesAndModifications();
        
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
        VBox taxonomicSection = createTaxonomicAnalysisSection();
        VBox modificationsSection = createModificationsSection();
        leftColumn.getChildren().addAll(enzymeSection, additionalParamsSection, taxonomicSection, modificationsSection);
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
    
    private VBox createEnzymeSection() {
        VBox section = new VBox(10);
        section.setPadding(new Insets(10));
        section.setStyle("-fx-background-color: #f5f5f5; -fx-background-radius: 5;");
        
        Label titleLabel = createSectionTitle("Enzyme Settings");
        
        // All enzyme parameters in one horizontal line with responsive sizing
        HBox enzymeParamsBox = new HBox(15);
        enzymeParamsBox.setAlignment(Pos.CENTER_LEFT);
        
        // Enzyme Selection
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
        
        // Missed Cleavages
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
        
        // Digestion Mode
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
        digestionModeComboBox.setValue("Specific"); // Set default value
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
        
        // Create a grid for better organization and responsive layout
        GridPane paramGrid = new GridPane();
        paramGrid.setHgap(15);
        paramGrid.setVgap(8);
        
        // Column constraints for responsive sizing
        ColumnConstraints col1 = new ColumnConstraints();
        col1.setMinWidth(150);
        col1.setHgrow(Priority.ALWAYS);
        ColumnConstraints col2 = new ColumnConstraints();
        col2.setMinWidth(150);
        col2.setHgrow(Priority.ALWAYS);
        paramGrid.getColumnConstraints().addAll(col1, col2);
        
        // Threads
        Label threadsLabel = new Label("Threads:");
        threadsLabel.setStyle("-fx-font-weight: bold;");
        threadsSpinner = new Spinner<>(1, ProcessorCount(), 4);
        threadsSpinner.setMinWidth(100);
        threadsSpinner.setPrefWidth(100);
        threadsSpinner.setMaxWidth(Double.MAX_VALUE);
        threadsSpinner.valueProperty().addListener((obs, oldVal, newVal) -> updateParameterView());
        
        VBox threadsBox = new VBox(5);
        threadsBox.getChildren().addAll(threadsLabel, threadsSpinner);
        
        // MS2 Scan Mode
        Label ms2Label = new Label("MS2 Scan Mode:");
        ms2Label.setStyle("-fx-font-weight: bold;");
        ms2ScanModeComboBox = new ComboBox<>(FXCollections.observableArrayList("FTMS", "ITMS"));
        ms2ScanModeComboBox.setValue("FTMS");
        ms2ScanModeComboBox.setMinWidth(100);
        ms2ScanModeComboBox.setPrefWidth(120);
        ms2ScanModeComboBox.setMaxWidth(Double.MAX_VALUE);
        ms2ScanModeComboBox.setOnAction(e -> updateParameterView());
        
        VBox ms2Box = new VBox(5);
        ms2Box.getChildren().addAll(ms2Label, ms2ScanModeComboBox);
        
        // Add to grid
        paramGrid.add(threadsBox, 0, 0);
        paramGrid.add(ms2Box, 1, 0);
        
        // Checkboxes with proper spacing
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
        
        RadioButton isotopicRadio = new RadioButton("Isotopic");
        isotopicRadio.setToggleGroup(quantModeGroup);
        
        RadioButton isobaricRadio = new RadioButton("Isobaric");
        isobaricRadio.setToggleGroup(quantModeGroup);
        
        HBox radioBox = new HBox(15);
        radioBox.setAlignment(Pos.CENTER_LEFT);
        radioBox.getChildren().addAll(labelFreeRadio, isotopicRadio, isobaricRadio);
        
        // StackPane for overlapping panels - both panels must occupy same space
        StackPane labelingStackPane = new StackPane();
        labelingStackPane.setMinHeight(200);
        labelingStackPane.setMaxHeight(Double.MAX_VALUE);
        VBox.setVgrow(labelingStackPane, Priority.ALWAYS);
        
        // Create both panels with identical sizing constraints
        VBox isotopicPanel = createCompactIsotopicLabelingPanel();
        VBox isobaricPanel = createCompactIsobaricLabelingPanel();
        
        // Ensure both panels have identical size constraints for proper StackPane behavior
        isotopicPanel.setMinHeight(200);
        isotopicPanel.setPrefHeight(220);
        isotopicPanel.setMaxHeight(Double.MAX_VALUE);
        
        isobaricPanel.setMinHeight(200);
        isobaricPanel.setPrefHeight(220);
        isobaricPanel.setMaxHeight(Double.MAX_VALUE);
        
        // Initially both panels are invisible
        isotopicPanel.setVisible(false);
        isotopicPanel.setManaged(false); // Key: don't take up space when hidden
        isobaricPanel.setVisible(false);
        isobaricPanel.setManaged(false); // Key: don't take up space when hidden
        
        labelingStackPane.getChildren().addAll(isotopicPanel, isobaricPanel);
        
        // Listeners to show/hide labeling panels with proper space management
        labelFreeRadio.setOnAction(e -> {
            isotopicPanel.setVisible(false);
            isotopicPanel.setManaged(false);
            isobaricPanel.setVisible(false);
            isobaricPanel.setManaged(false);
            updateParameterView();
            updateMetadataTable();
        });
        
        isotopicRadio.setOnAction(e -> {
            isotopicPanel.setVisible(true);
            isotopicPanel.setManaged(true);
            isobaricPanel.setVisible(false);
            isobaricPanel.setManaged(false);
            updateParameterView();
            updateMetadataTable();
        });
        
        isobaricRadio.setOnAction(e -> {
            isotopicPanel.setVisible(false);
            isotopicPanel.setManaged(false);
            isobaricPanel.setVisible(true);
            isobaricPanel.setManaged(true);
            updateParameterView();
            updateMetadataTable();
        });
        
        section.getChildren().addAll(titleLabel, radioBox, labelingStackPane);
        return section;
    }
    
    private VBox createCompactIsotopicLabelingPanel() {
        VBox panel = new VBox(8);
        panel.setPadding(new Insets(8));
        panel.setStyle("-fx-background-color: #e8f4f8; -fx-background-radius: 5;");
        
        // Label selection row
        HBox labelSelectionBox = new HBox(8);
        labelSelectionBox.setAlignment(Pos.CENTER_LEFT);
        
        Label selectLabel = new Label("Select label:");
        selectLabel.setMinWidth(120);
        selectLabel.setPrefWidth(120);
        labelModificationComboBox = new ComboBox<>(labelModifications);
        labelModificationComboBox.setMinWidth(120);
        labelModificationComboBox.setPrefWidth(180);
        labelModificationComboBox.setMaxWidth(Double.MAX_VALUE);
        labelModificationComboBox.setPromptText("Select label");
        HBox.setHgrow(labelModificationComboBox, Priority.ALWAYS);
        
        Button addLightBtn = new Button("Add Light");
        addLightBtn.setMinWidth(100);
        addLightBtn.setPrefWidth(120);
        addLightBtn.setOnAction(e -> addLightLabel());
        
        Button addMediumBtn = new Button("Add Medium");
        addMediumBtn.setMinWidth(100);
        addMediumBtn.setPrefWidth(120);
        addMediumBtn.setOnAction(e -> addMediumLabel());
        
        Button addHeavyBtn = new Button("Add Heavy");
        addHeavyBtn.setMinWidth(100);
        addHeavyBtn.setPrefWidth(120);
        addHeavyBtn.setOnAction(e -> addHeavyLabel());
        
        labelSelectionBox.getChildren().addAll(selectLabel, labelModificationComboBox, addLightBtn, addMediumBtn, addHeavyBtn);
        
        // Lists with remove buttons - use VBox.setVgrow to fill available space
        HBox listsBox = new HBox(8);
        listsBox.setAlignment(Pos.TOP_LEFT);
        VBox.setVgrow(listsBox, Priority.ALWAYS); // Allow lists box to grow
        
        // Light labels section
        VBox lightBox = new VBox(5);
        lightBox.setMinWidth(120);
        HBox.setHgrow(lightBox, Priority.ALWAYS);
        
        Label lightLabel = new Label("Light Labels:");
        lightLabel.setStyle("-fx-font-weight: bold;");
        
        lightLabelListView = new ListView<>(lightLabels);
        lightLabelListView.setMinHeight(80);
        lightLabelListView.setPrefHeight(100);
        lightLabelListView.setMaxHeight(Double.MAX_VALUE); // Allow to grow
        VBox.setVgrow(lightLabelListView, Priority.ALWAYS); // Fill available space
        
        Button removeLightBtn = new Button("Remove Selected");
        removeLightBtn.setMinWidth(120);
        removeLightBtn.setPrefWidth(120);
        removeLightBtn.setMaxWidth(Double.MAX_VALUE);
        removeLightBtn.setOnAction(e -> removeLightLabel());
        
        lightBox.getChildren().addAll(lightLabel, lightLabelListView, removeLightBtn);
        
        // Medium labels section
        VBox mediumBox = new VBox(5);
        mediumBox.setMinWidth(120);
        HBox.setHgrow(mediumBox, Priority.ALWAYS);
        
        Label mediumLabel = new Label("Medium Labels:");
        mediumLabel.setStyle("-fx-font-weight: bold;");
        
        mediumLabelListView = new ListView<>(mediumLabels);
        mediumLabelListView.setMinHeight(80);
        mediumLabelListView.setPrefHeight(100);
        mediumLabelListView.setMaxHeight(Double.MAX_VALUE); // Allow to grow
        VBox.setVgrow(mediumLabelListView, Priority.ALWAYS); // Fill available space
        
        Button removeMediumBtn = new Button("Remove Selected");
        removeMediumBtn.setMinWidth(120);
        removeMediumBtn.setPrefWidth(120);
        removeMediumBtn.setMaxWidth(Double.MAX_VALUE);
        removeMediumBtn.setOnAction(e -> removeMediumLabel());
        
        mediumBox.getChildren().addAll(mediumLabel, mediumLabelListView, removeMediumBtn);
        
        // Heavy labels section
        VBox heavyBox = new VBox(5);
        heavyBox.setMinWidth(120);
        HBox.setHgrow(heavyBox, Priority.ALWAYS);
        
        Label heavyLabel = new Label("Heavy Labels:");
        heavyLabel.setStyle("-fx-font-weight: bold;");
        
        heavyLabelListView = new ListView<>(heavyLabels);
        heavyLabelListView.setMinHeight(80);
        heavyLabelListView.setPrefHeight(100);
        heavyLabelListView.setMaxHeight(Double.MAX_VALUE); // Allow to grow
        VBox.setVgrow(heavyLabelListView, Priority.ALWAYS); // Fill available space
        
        Button removeHeavyBtn = new Button("Remove Selected");
        removeHeavyBtn.setMinWidth(120);
        removeHeavyBtn.setPrefWidth(120);
        removeHeavyBtn.setMaxWidth(Double.MAX_VALUE);
        removeHeavyBtn.setOnAction(e -> removeHeavyLabel());
        
        heavyBox.getChildren().addAll(heavyLabel, heavyLabelListView, removeHeavyBtn);
        
        listsBox.getChildren().addAll(lightBox, mediumBox, heavyBox);
        
        panel.getChildren().addAll(labelSelectionBox, listsBox);
        
        return panel;
    }
    
    private VBox createCompactIsobaricLabelingPanel() {
        VBox panel = new VBox(8);
        panel.setPadding(new Insets(8));
        panel.setStyle("-fx-background-color: #f8e8f4; -fx-background-radius: 5;");
        // Remove fixed height constraint to allow filling
        
        // Isobaric type selection
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
        
        // Isobaric table - allow it to fill available space
        isobaricTable = new TableView<>(isobaricTableData);
        isobaricTable.setMinHeight(120);
        isobaricTable.setPrefHeight(140);
        isobaricTable.setMaxHeight(Double.MAX_VALUE); // Allow to grow
        isobaricTable.setEditable(true);
        isobaricTable.setColumnResizePolicy(TableView.CONSTRAINED_RESIZE_POLICY);
        VBox.setVgrow(isobaricTable, Priority.ALWAYS); // Fill available space
        
        // Table columns with responsive widths
        TableColumn<IsobaricLabelRow, String> labelCol = new TableColumn<>("Label");
        labelCol.setCellValueFactory(new PropertyValueFactory<>("labelName"));
        labelCol.setMinWidth(120);
        labelCol.setPrefWidth(150);
        
        TableColumn<IsobaricLabelRow, String> minus2Col = new TableColumn<>("-2");
        minus2Col.setCellValueFactory(new PropertyValueFactory<>("minus2"));
        minus2Col.setCellFactory(TextFieldTableCell.forTableColumn(new DefaultStringConverter()));
        minus2Col.setOnEditCommit(e -> e.getRowValue().setMinus2(e.getNewValue()));
        minus2Col.setMinWidth(50);
        minus2Col.setPrefWidth(60);
        
        TableColumn<IsobaricLabelRow, String> minus1Col = new TableColumn<>("-1");
        minus1Col.setCellValueFactory(new PropertyValueFactory<>("minus1"));
        minus1Col.setCellFactory(TextFieldTableCell.forTableColumn(new DefaultStringConverter()));
        minus1Col.setOnEditCommit(e -> e.getRowValue().setMinus1(e.getNewValue()));
        minus1Col.setMinWidth(50);
        minus1Col.setPrefWidth(60);
        
        TableColumn<IsobaricLabelRow, String> plus1Col = new TableColumn<>("+1");
        plus1Col.setCellValueFactory(new PropertyValueFactory<>("plus1"));
        plus1Col.setCellFactory(TextFieldTableCell.forTableColumn(new DefaultStringConverter()));
        plus1Col.setOnEditCommit(e -> e.getRowValue().setPlus1(e.getNewValue()));
        plus1Col.setMinWidth(50);
        plus1Col.setPrefWidth(60);
        
        TableColumn<IsobaricLabelRow, String> plus2Col = new TableColumn<>("+2");
        plus2Col.setCellValueFactory(new PropertyValueFactory<>("plus2"));
        plus2Col.setCellFactory(TextFieldTableCell.forTableColumn(new DefaultStringConverter()));
        plus2Col.setOnEditCommit(e -> e.getRowValue().setPlus2(e.getNewValue()));
        plus2Col.setMinWidth(50);
        plus2Col.setPrefWidth(60);
        
        isobaricTable.getColumns().addAll(labelCol, minus2Col, minus1Col, plus1Col, plus2Col);
        
        panel.getChildren().addAll(typeBox, isobaricTable);
        return panel;
    }
    
    private VBox createTaxonomicAnalysisSection() {
        VBox section = new VBox(8);
        section.setPadding(new Insets(8));
        section.setStyle("-fx-background-color: #f0fff0; -fx-background-radius: 5;");
        
        Label titleLabel = createSectionTitle("Taxonomic Analysis");
        
        // Taxa selection (horizontal compact grid - more columns, fewer rows)
        Label taxaLabel = new Label("Taxa from which top nodes showed in result:");
        taxaLabel.setStyle("-fx-font-weight: bold; -fx-font-size: 11px;");
        
        GridPane taxaGrid = new GridPane();
        taxaGrid.setHgap(8);
        taxaGrid.setVgap(5);
        
        eukaryotaCheckBox = new CheckBox("Eukaryota");
        bacteriaCheckBox = new CheckBox("Bacteria");
        archaeaCheckBox = new CheckBox("Archaea");
        virusesCheckBox = new CheckBox("Viruses");
        viroidsCheckBox = new CheckBox("Viroids");
        unclassifiedCheckBox = new CheckBox("Unclassified");
        otherCheckBox = new CheckBox("Other");
        
        // Arrange in 4 columns, 2 rows for better width utilization
        taxaGrid.add(eukaryotaCheckBox, 0, 0);
        taxaGrid.add(bacteriaCheckBox, 1, 0);
        taxaGrid.add(archaeaCheckBox, 2, 0);
        taxaGrid.add(virusesCheckBox, 3, 0);
        taxaGrid.add(viroidsCheckBox, 0, 1);
        taxaGrid.add(unclassifiedCheckBox, 1, 1);
        taxaGrid.add(otherCheckBox, 2, 1);
        
        // Ignore rank - compact horizontal layout
        HBox rankBox = new HBox(8);
        rankBox.setAlignment(Pos.CENTER_LEFT);
        Label rankLabel = new Label("Ignore blanks below rank:");
        rankLabel.setStyle("-fx-font-weight: bold; -fx-font-size: 11px;");
        ignoreRankComboBox = new ComboBox<>(FXCollections.observableArrayList(
            "Superkingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"));
        ignoreRankComboBox.setValue("Species");
        ignoreRankComboBox.setMinWidth(80);
        ignoreRankComboBox.setPrefWidth(100);
        ignoreRankComboBox.setMaxWidth(Double.MAX_VALUE);
        ignoreRankComboBox.setOnAction(e -> updateParameterView());
        rankBox.getChildren().addAll(rankLabel, ignoreRankComboBox);
        
        environmentalSamplesCheckBox = new CheckBox("Ignore environmental samples");
        environmentalSamplesCheckBox.setOnAction(e -> updateParameterView());
        
        section.getChildren().addAll(titleLabel, taxaLabel, taxaGrid, rankBox, environmentalSamplesCheckBox);
        
        return section;
    }
    
    private VBox createMetadataConfigurationSection() {
        VBox section = new VBox(10);
        section.setPadding(new Insets(10));
        section.setStyle("-fx-background-color: #fff8f0; -fx-background-radius: 5;");
        section.setMinHeight(200);
        VBox.setVgrow(section, Priority.ALWAYS);
        
        Label titleLabel = createSectionTitle("Metadata Configuration");
        
        // Dynamic content based on quantitative mode
        Label infoLabel = new Label("Metadata will be configured based on selected quantitative mode");
        infoLabel.setStyle("-fx-font-style: italic; -fx-text-fill: #666666;");
        
        // Metadata table
        metadataTable = new TableView<>(metadataTableData);
        metadataTable.setMinHeight(120);
        metadataTable.setPrefHeight(150);
        metadataTable.setMaxHeight(Double.MAX_VALUE);
        metadataTable.setEditable(true);
        metadataTable.setColumnResizePolicy(TableView.CONSTRAINED_RESIZE_POLICY);
        VBox.setVgrow(metadataTable, Priority.ALWAYS);
        
        // Initially create basic columns
        createMetadataTableColumns();
        
        // Control buttons
        HBox controlBox = new HBox(8);
        controlBox.setAlignment(Pos.CENTER_LEFT);
        
        metadataColumnComboBox = new ComboBox<>(FXCollections.observableArrayList(
            "Age", "Gender", "Condition", "Treatment", "Time Point", "BMI", "Custom"));
        metadataColumnComboBox.setMinWidth(100);
        metadataColumnComboBox.setPrefWidth(120);
        metadataColumnComboBox.setPromptText("Add column");
        
        // Custom column TextField (hidden by default)
        TextField customColumnTextField = new TextField();
        customColumnTextField.setPromptText("Enter custom column name");
        customColumnTextField.setMinWidth(120);
        customColumnTextField.setVisible(false);
        customColumnTextField.setManaged(false);
        
        // Show/hide the custom column TextField based on ComboBox selection
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
        
        // Update metadata table when quantitative mode changes
        updateMetadataTable();
        
        return section;
    }
    
    private void createMetadataTableColumns() {
        metadataTable.getColumns().clear();
        
        // Sample/Experiment column
        TableColumn<MetadataRow, String> sampleCol = new TableColumn<>("Sample");
        sampleCol.setCellValueFactory(new PropertyValueFactory<>("sampleName"));
        sampleCol.setMinWidth(100);
        sampleCol.setPrefWidth(120);
        
        metadataTable.getColumns().add(sampleCol);
        
        // Add label type column if using isotopic or isobaric labeling
        if (quantModeGroup.getSelectedToggle() != null) {
            RadioButton selected = (RadioButton) quantModeGroup.getSelectedToggle();
            String mode = selected.getText();
            
            if (mode.contains("Isotopic") || mode.contains("Isobaric")) {
                TableColumn<MetadataRow, String> labelCol = new TableColumn<>("Label Type");
                labelCol.setCellValueFactory(new PropertyValueFactory<>("labelType"));
                labelCol.setMinWidth(120);
                labelCol.setPrefWidth(150);
                metadataTable.getColumns().add(labelCol);
            }
        }
        
        // Add common metadata columns
        addMetadataColumnWithProperty("Age", "age");
        addMetadataColumnWithProperty("Gender", "gender");
        addMetadataColumnWithProperty("Condition", "condition");
    }
    
    private void updateMetadataTable() {
        // Always show label-free rows by default
        metadataTableData.clear();
        List<String> experimentNames = ioPanel.getExperimentNames();
        if (experimentNames.isEmpty()) return;
        createMetadataTableColumns();
        // Determine quantitative mode
        if (quantModeGroup.getSelectedToggle() != null) {
            RadioButton selected = (RadioButton) quantModeGroup.getSelectedToggle();
            String mode = selected.getText();
            boolean hasLabels = false;
            if (mode.contains("Isotopic")) {
                hasLabels = !lightLabels.isEmpty() || !mediumLabels.isEmpty() || !heavyLabels.isEmpty();
            } else if (mode.contains("Isobaric")) {
                hasLabels = !isobaricTableData.isEmpty();
            }
            if (mode.contains("Isotopic") && hasLabels) {
                for (String expName : experimentNames) {
                    for (String label : lightLabels) {
                        MetadataRow row = new MetadataRow(expName);
                        row.setLabelType("Light - " + label);
                        metadataTableData.add(row);
                    }
                    for (String label : mediumLabels) {
                        MetadataRow row = new MetadataRow(expName);
                        row.setLabelType("Medium - " + label);
                        metadataTableData.add(row);
                    }
                    for (String label : heavyLabels) {
                        MetadataRow row = new MetadataRow(expName);
                        row.setLabelType("Heavy - " + label);
                        metadataTableData.add(row);
                    }
                }
            } else if (mode.contains("Isobaric") && hasLabels) {
                for (String expName : experimentNames) {
                    for (IsobaricLabelRow isoLabel : isobaricTableData) {
                        MetadataRow row = new MetadataRow(expName);
                        row.setLabelType("Isobaric - " + isoLabel.getLabelName());
                        metadataTableData.add(row);
                    }
                }
            } else {
                // Show label-free by default
                for (String expName : experimentNames) {
                    metadataTableData.add(new MetadataRow(expName));
                }
            }
        }
    }
    
    private void addMetadataColumn() {
        String columnName = metadataColumnComboBox.getValue();
        if (columnName != null && !columnName.isEmpty()) {
            if ("Custom".equals(columnName)) {
                // Could add a dialog for custom column name
                columnName = "Custom" + (metadataTable.getColumns().size());
            }
            
            String propertyName = getPropertyNameForColumn(columnName);
            addMetadataColumnWithProperty(columnName, propertyName);
        }
    }
    
    private void addMetadataColumnWithProperty(String columnName, String propertyName) {
        TableColumn<MetadataRow, String> column = new TableColumn<>(columnName);
        column.setCellValueFactory(new PropertyValueFactory<>(propertyName));
        column.setCellFactory(TextFieldTableCell.forTableColumn(new DefaultStringConverter()));
        column.setOnEditCommit(event -> {
            MetadataRow row = event.getRowValue();
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
        if (metadataTable.getColumns().size() > 1) { // Keep at least the sample column
            metadataTable.getColumns().remove(metadataTable.getColumns().size() - 1);
        }
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
        
        // Column constraints for responsive sizing
        ColumnConstraints col1 = new ColumnConstraints();
        col1.setMinWidth(200);
        col1.setHgrow(Priority.ALWAYS);
        ColumnConstraints col2 = new ColumnConstraints();
        col2.setMinWidth(200);
        col2.setHgrow(Priority.ALWAYS);
        modGrid.getColumnConstraints().addAll(col1, col2);
        
        // Row constraints for proper height distribution
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
        
        // Fixed Modifications
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
        
        // Variable Modifications
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
        
        // Modification Selection with responsive layout
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
        
        // Add components to grid
        modGrid.add(fixedModLabel, 0, 0);
        modGrid.add(variableModLabel, 1, 0);
        modGrid.add(fixedModListView, 0, 1);
        modGrid.add(variableModListView, 1, 1);
        modGrid.add(removeFixedModButton, 0, 2);
        modGrid.add(removeVariableModButton, 1, 2);
        modGrid.add(modSelectionContainer, 0, 3, 2, 1); // Span 2 columns
        
        section.getChildren().addAll(titleLabel, modGrid);
        return section;
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
        startButton.setOnAction(e -> startAnalysis());
        stopButton.setOnAction(e -> stopAnalysis());
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
    
    // Inner class for metadata table rows
    public static class MetadataRow {
        private final SimpleStringProperty sampleName;
        private final SimpleStringProperty labelType;
        private final SimpleStringProperty age;
        private final SimpleStringProperty gender;
        private final SimpleStringProperty condition;
        private final SimpleStringProperty treatment;
        private final SimpleStringProperty timePoint;
        private final SimpleStringProperty customField1;
        private final SimpleStringProperty customField2;
        private final SimpleStringProperty customField3;
        
        public MetadataRow(String sampleName) {
            this.sampleName = new SimpleStringProperty(sampleName);
            this.labelType = new SimpleStringProperty("");
            this.age = new SimpleStringProperty("");
            this.gender = new SimpleStringProperty("");
            this.condition = new SimpleStringProperty("");
            this.treatment = new SimpleStringProperty("");
            this.timePoint = new SimpleStringProperty("");
            this.customField1 = new SimpleStringProperty("");
            this.customField2 = new SimpleStringProperty("");
            this.customField3 = new SimpleStringProperty("");
        }
        
        // Getters and setters
        public String getSampleName() { return sampleName.get(); }
        public void setSampleName(String value) { sampleName.set(value); }
        public StringProperty sampleNameProperty() { return sampleName; }
        
        public String getLabelType() { return labelType.get(); }
        public void setLabelType(String value) { labelType.set(value); }
        public StringProperty labelTypeProperty() { return labelType; }
        
        public String getAge() { return age.get(); }
        public void setAge(String value) { age.set(value); }
        public StringProperty ageProperty() { return age; }
        
        public String getGender() { return gender.get(); }
        public void setGender(String value) { gender.set(value); }
        public StringProperty genderProperty() { return gender; }
        
        public String getCondition() { return condition.get(); }
        public void setCondition(String value) { condition.set(value); }
        public StringProperty conditionProperty() { return condition; }
        
        public String getTreatment() { return treatment.get(); }
        public void setTreatment(String value) { treatment.set(value); }
        public StringProperty treatmentProperty() { return treatment; }
        
        public String getTimePoint() { return timePoint.get(); }
        public void setTimePoint(String value) { timePoint.set(value); }
        public StringProperty timePointProperty() { return timePoint; }
        
        public String getCustomField1() { return customField1.get(); }
        public void setCustomField1(String value) { customField1.set(value); }
        public StringProperty customField1Property() { return customField1; }
        
        public String getCustomField2() { return customField2.get(); }
        public void setCustomField2(String value) { customField2.set(value); }
        public StringProperty customField2Property() { return customField2; }
        
        public String getCustomField3() { return customField3.get(); }
        public void setCustomField3(String value) { customField3.set(value); }
        public StringProperty customField3Property() { return customField3; }
    }
    
    private void loadEnzymesAndModifications() {
        try {
            // Load enzymes from XML using proper resource loading
            InputStream enzymesStream = getClass().getClassLoader().getResourceAsStream("enzymes.xml");
            if (enzymesStream != null) {
                DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
                DocumentBuilder builder = factory.newDocumentBuilder();
                Document doc = builder.parse(enzymesStream);
                
                NodeList enzymeNodes = doc.getElementsByTagName("enzyme");
                for (int i = 0; i < enzymeNodes.getLength(); i++) {
                    Element enzyme = (Element) enzymeNodes.item(i);
                    String title = enzyme.getAttribute("title");
                    if (title != null && !title.isEmpty()) {
                        enzymes.add(title);
                    }
                }
                enzymesStream.close();
                System.out.println("Loaded " + enzymes.size() + " enzymes from enzymes.xml");
            } else {
                showError("Configuration Error", "Could not find enzymes.xml in resources");
                System.err.println("ERROR: enzymes.xml not found in classpath");
            }
            
            // Load modifications from XML using proper resource loading
            InputStream modificationsStream = getClass().getClassLoader().getResourceAsStream("modifications.xml");
            if (modificationsStream != null) {
                DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
                DocumentBuilder builder = factory.newDocumentBuilder();
                Document doc = builder.parse(modificationsStream);
                
                NodeList modNodes = doc.getElementsByTagName("modification");
                for (int i = 0; i < modNodes.getLength(); i++) {
                    Element mod = (Element) modNodes.item(i);
                    String title = mod.getAttribute("title");
                    if (title != null && !title.isEmpty()) {
                        modifications.add(title);
                    }
                    // Add label modifications
                    NodeList typeNodes = mod.getElementsByTagName("type");
                    if (typeNodes.getLength() > 0) {
                        String type = typeNodes.item(0).getTextContent();
                        if ("Label".equalsIgnoreCase(type)) {
                            if (title != null && !title.isEmpty() && !labelModifications.contains(title)) {
                                labelModifications.add(title);
                            }
                        } else if ("IsobaricLabel".equalsIgnoreCase(type)) {
                            if (title != null && !title.isEmpty() && !isobaricModifications.contains(title)) {
                                isobaricModifications.add(title);
                                // Extract base type (e.g., "TMTpro16plex" from "TMTpro16plex-Nter133C")
                                String baseType = extractIsobaricBaseType(title);
                                if (baseType != null && !isobaricBaseTypes.contains(baseType)) {
                                    isobaricBaseTypes.add(baseType);
                                }
                            }
                        }
                    }
                }
                modificationsStream.close();
                System.out.println("Loaded " + modifications.size() + " modifications from modifications.xml");
                System.out.println("Loaded " + labelModifications.size() + " label modifications");
                System.out.println("Loaded " + isobaricModifications.size() + " isobaric modifications");
                System.out.println("Loaded " + isobaricBaseTypes.size() + " isobaric base types");
            } else {
                showError("Configuration Error", "Could not find modifications.xml in resources");
                System.err.println("ERROR: modifications.xml not found in classpath");
            }
        } catch (Exception e) {
            showError("Configuration Error", "Failed to load enzymes and modifications: " + e.getMessage());
            e.printStackTrace();
        }
    }
    
    private String extractIsobaricBaseType(String fullName) {
        // Extract base type from names like "TMTpro16plex-Nter133C" -> "TMTpro16plex"
        if (fullName.contains("-")) {
            return fullName.substring(0, fullName.indexOf("-"));
        }
        return fullName;
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
            
            if (selected.getText().contains("Isotopic")) {
                sb.append("Light Labels (").append(lightLabels.size()).append("): ");
                sb.append(lightLabels.isEmpty() ? "None" : String.join(", ", lightLabels)).append("\n");
                sb.append("Medium Labels (").append(mediumLabels.size()).append("): ");
                sb.append(mediumLabels.isEmpty() ? "None" : String.join(", ", mediumLabels)).append("\n");
                sb.append("Heavy Labels (").append(heavyLabels.size()).append("): ");
                sb.append(heavyLabels.isEmpty() ? "None" : String.join(", ", heavyLabels)).append("\n");
            } else if (selected.getText().contains("Isobaric")) {
                sb.append("Isobaric Type: ").append(isobaricTypeComboBox.getValue() != null ? isobaricTypeComboBox.getValue() : "None").append("\n");
                sb.append("Configured Labels: ").append(isobaricTableData.size()).append("\n");
            }
        }
        
        // Taxonomic analysis
        sb.append("\n=== TAXONOMIC ANALYSIS ===\n");
        sb.append("Eukaryota: ").append(eukaryotaCheckBox.isSelected() ? "Yes" : "No").append("\n");
        sb.append("Bacteria: ").append(bacteriaCheckBox.isSelected() ? "Yes" : "No").append("\n");
        sb.append("Archaea: ").append(archaeaCheckBox.isSelected() ? "Yes" : "No").append("\n");
        sb.append("Viruses: ").append(virusesCheckBox.isSelected() ? "Yes" : "No").append("\n");
        sb.append("Viroids: ").append(viroidsCheckBox.isSelected() ? "Yes" : "No").append("\n");
        sb.append("Unclassified: ").append(unclassifiedCheckBox.isSelected() ? "Yes" : "No").append("\n");
        sb.append("Other: ").append(otherCheckBox.isSelected() ? "Yes" : "No").append("\n");
        sb.append("Ignore Rank: ").append(ignoreRankComboBox.getValue()).append("\n");
        sb.append("Ignore Environmental Samples: ").append(environmentalSamplesCheckBox.isSelected() ? "Yes" : "No").append("\n");
        
        parameterView.setText(sb.toString());
    }
    
    private boolean validateInputs() {
        if (ioPanel.getExperimentData().isEmpty()) {
            showError("Validation Error", "Please add at least one experiment.");
            return false;
        }
        if (ioPanel.getOutputDirectory().isEmpty()) {
            showError("Validation Error", "Please select an output directory.");
            return false;
        }
        if (microbiomeDBField.getText().isEmpty()) {
            showError("Validation Error", "Please select a microbiome database file.");
            return false;
        }
        if (appendHostDBCheckBox.isSelected() && hostDBField.getText().isEmpty()) {
            showError("Validation Error", "Please select a host database file.");
            return false;
        }
        if (enzymeComboBox.getValue() == null) {
            showError("Validation Error", "Please select an enzyme.");
            return false;
        }
        return true;
    }
    
    private void startAnalysis() {
        try {
            // Collect parameters from UI
            ObservableList<ExperimentRow> fileRows = ioPanel.getExperimentData();
            if (fileRows.isEmpty()) {
                showAlert("Error", "No files selected for analysis.");
                return;
            }
            String[] rawFiles = new String[fileRows.size()];
            String[] expNames = new String[fileRows.size()];
            int[] fractions = new int[fileRows.size()];
            int[] replicates = new int[fileRows.size()];
            for (int i = 0; i < fileRows.size(); i++) {
                ExperimentRow row = fileRows.get(i);
                rawFiles[i] = row.getFilePath();
                expNames[i] = row.getExperimentName();
                fractions[i] = 0; // Default fraction
                replicates[i] = 1; // Default replicate
            }
            MetaData metadata = new MetaData(rawFiles, expNames, fractions, replicates);
            String resultDir = ioPanel.getOutputDirectory();
            if (resultDir == null || resultDir.trim().isEmpty()) {
                showAlert("Error", "Please select a result directory.");
                return;
            }
            String microDb = microbiomeDBField.getText();
            if (microDb == null || microDb.trim().isEmpty()) {
                showAlert("Error", "Please select a database file.");
                return;
            }
            boolean appendHostDb = appendHostDBCheckBox.isSelected();
            String hostDB = appendHostDb ? hostDBField.getText() : "";
            String selectedEnzymeStr = enzymeComboBox.getValue();
            if (selectedEnzymeStr == null) {
                selectedEnzymeStr = "Trypsin"; // Default
            }
            int missCleavages = missedCleavagesSpinner.getValue();
            String digestMode = digestionModeComboBox.getValue();
            // 1. Create MetaBaseParameter
            String workflowType = MetaLabWorkflowType.MaxQuantWorkflow.name();
            String ms2ScanMode = MetaConstants.HCD_FTMS; // Default
            int coreCount = Runtime.getRuntime().availableProcessors();
            int threadCount = Math.max(1, coreCount - 1); // Leave one core free
            boolean isMetaWorkflow = true;
            MetaBaseParameter metaBasePar = new MetaBaseParameter(workflowType, metadata, resultDir, microDb, 
                    appendHostDb, hostDB, ms2ScanMode, coreCount, threadCount, isMetaWorkflow);
            // 2. Set modifications (String[])
            metaBasePar.setFixMods(fixedModifications.toArray(new String[0]));
            metaBasePar.setVariMods(variableModifications.toArray(new String[0]));
            // 3. Set enzyme and digestion mode (String)
            metaBasePar.setEnzyme(selectedEnzymeStr);
            metaBasePar.setDigestMode(digestMode);
            metaBasePar.setMissCleavages(missCleavages);
            // 4. Set isotopic labels (String[][])
            String[][] labels = new String[3][];
            labels[0] = lightLabels.toArray(new String[0]);
            labels[1] = mediumLabels.toArray(new String[0]);
            labels[2] = heavyLabels.toArray(new String[0]);
            metaBasePar.setIsotopicLabels(labels);
            // 5. Set isobaric tag (IsobaricTag) from ComboBox using getName()
            String isobaricType = isobaricTypeComboBox.getValue();
            if (isobaricType != null && !isobaricType.isEmpty()) {
                for (IsobaricTag tag : IsobaricTag.values()) {
                    if (tag.getName().equals(isobaricType)) {
                        metaBasePar.setIsobaricTag(tag);
                        break;
                    }
                }
            } else {
                metaBasePar.setIsobaricTag(null);
            }
            // 6. Set isobaric correction factors (double[][])
            if (!isobaricTableData.isEmpty()) {
                double[][] isoCorFactor = new double[isobaricTableData.size()][4];
                for (int i = 0; i < isobaricTableData.size(); i++) {
                    IsobaricLabelRow row = isobaricTableData.get(i);
                    try {
                        isoCorFactor[i][0] = Double.parseDouble(row.getMinus2());
                        isoCorFactor[i][1] = Double.parseDouble(row.getMinus1());
                        isoCorFactor[i][2] = Double.parseDouble(row.getPlus1());
                        isoCorFactor[i][3] = Double.parseDouble(row.getPlus2());
                    } catch (NumberFormatException e) {
                        isoCorFactor[i][0] = 0.0;
                        isoCorFactor[i][1] = 0.0;
                        isoCorFactor[i][2] = 0.0;
                        isoCorFactor[i][3] = 0.0;
                    }
                }
                metaBasePar.setIsobaricRefTable(isoCorFactor);
            }
            // 7. Create MetaResourceParameters
            MetaResourceParameters mrp = new MetaResourceParameters();
            // 8. Create JavaFXProgressBar for overall and current task
            uk.ac.quadram.metapilot.gui.JavaFXProgressBar overallBar = new uk.ac.quadram.metapilot.gui.JavaFXProgressBar(0, 100, 0);
            uk.ac.quadram.metapilot.gui.JavaFXProgressBar currentBar = new uk.ac.quadram.metapilot.gui.JavaFXProgressBar(0, 100, 0);
            // 9. Bind progress bars to UI
            overallProgressBar.progressProperty().bind(overallBar.getProgressBar().progressProperty());
            currentTaskProgressBar.progressProperty().bind(currentBar.getProgressBar().progressProperty());
            // 10. Create and run DDAIgcMaxQuantTask
            DDAIgcMaxQuantTask task = new DDAIgcMaxQuantTask(metaBasePar, mrp, overallBar, currentBar);
            // 11. Bind task message to currentTaskLabel
            task.messageProperty().addListener((obs, oldMsg, newMsg) -> {
                Platform.runLater(() -> currentTaskLabel.setText("Current Task: " + newMsg));
            });
            // 12. Handle task completion
            task.setOnSucceeded(e -> {
                Platform.runLater(() -> {
                    startButton.setDisable(false);
                    stopButton.setDisable(true);
                    if (statusTimeline != null) statusTimeline.stop();
                    showAlert("Success", "Analysis completed successfully.");
                });
            });
            task.setOnFailed(e -> {
                Platform.runLater(() -> {
                    startButton.setDisable(false);
                    stopButton.setDisable(true);
                    if (statusTimeline != null) statusTimeline.stop();
                    showAlert("Error", "Analysis failed: " + task.getException());
                });
            });
            // 13. Start the task in a new thread
            Thread t = new Thread(task);
            t.setDaemon(true);
            t.start();
            // 14. Update UI buttons
            startButton.setDisable(true);
            stopButton.setDisable(false);
            // 15. Start status panel timer
            taskStartTimeMillis = System.currentTimeMillis();
            if (statusTimeline != null) statusTimeline.stop();
            statusTimeline = new Timeline(new KeyFrame(Duration.seconds(1), ev -> updateStatusPanel()));
            statusTimeline.setCycleCount(Timeline.INDEFINITE);
            statusTimeline.play();
        } catch (Exception e) {
            showAlert("Error", "Failed to start analysis: " + e.getMessage());
            e.printStackTrace();
        }
    }
    
    private void stopAnalysis() {
        Alert alert = new Alert(AlertType.CONFIRMATION);
        alert.setTitle("Stop Analysis");
        alert.setHeaderText("Are you sure you want to stop the analysis?");
        alert.setContentText("This action cannot be undone.");
        
        Optional<ButtonType> result = alert.showAndWait();
        if (result.isPresent() && result.get() == ButtonType.OK) {
            startButton.setDisable(false);
            stopButton.setDisable(true);
            overallProgressBar.progressProperty().unbind();
            overallProgressBar.setProgress(0);
            currentTaskProgressBar.setProgress(0);
            currentTaskLabel.setText("Current Task: None");
        }
    }
    
    private void showError(String title, String content) {
        Platform.runLater(() -> {
            Alert alert = new Alert(AlertType.ERROR);
            alert.setTitle(title);
            alert.setHeaderText(null);
            alert.setContentText(content);
            alert.showAndWait();
        });
    }
    
    private Label createSectionTitle(String title) {
        Label label = new Label(title);
        label.setFont(Font.font("System", FontWeight.BOLD, 14));
        label.setTextFill(Color.DARKBLUE);
        return label;
    }
    
    private void updateIsobaricTable() {
        String selectedType = isobaricTypeComboBox.getValue();
        if (selectedType != null) {
            isobaricTableData.clear();
            List<String> labels = getIsobaricLabelsForType(selectedType);
            for (String label : labels) {
                isobaricTableData.add(new IsobaricLabelRow(label));
            }
            updateMetadataTable();
        }
    }
    
    private List<String> getIsobaricLabelsForType(String baseType) {
        List<String> labels = new ArrayList<>();
        
        switch (baseType) {
            case "iTRAQ4plex":
                labels.add("iTRAQ4plex-114");
                labels.add("iTRAQ4plex-115");
                labels.add("iTRAQ4plex-116");
                labels.add("iTRAQ4plex-117");
                break;
                
            case "iTRAQ8plex":
                labels.add("iTRAQ8plex-113");
                labels.add("iTRAQ8plex-114");
                labels.add("iTRAQ8plex-115");
                labels.add("iTRAQ8plex-116");
                labels.add("iTRAQ8plex-117");
                labels.add("iTRAQ8plex-118");
                labels.add("iTRAQ8plex-119");
                labels.add("iTRAQ8plex-121");
                break;
                
            case "TMT6plex":
                labels.add("TMT6plex-126");
                labels.add("TMT6plex-127");
                labels.add("TMT6plex-128");
                labels.add("TMT6plex-129");
                labels.add("TMT6plex-130");
                labels.add("TMT6plex-131");
                break;
                
            case "TMT10plex":
                labels.add("TMT10plex-126");
                labels.add("TMT10plex-127N");
                labels.add("TMT10plex-127C");
                labels.add("TMT10plex-128N");
                labels.add("TMT10plex-128C");
                labels.add("TMT10plex-129N");
                labels.add("TMT10plex-129C");
                labels.add("TMT10plex-130N");
                labels.add("TMT10plex-130C");
                labels.add("TMT10plex-131");
                break;
                
            case "TMT11plex":
                // Include all TMT10plex labels
                labels.addAll(getIsobaricLabelsForType("TMT10plex"));
                // Add TMT11plex specific labels
                labels.add("TMT11plex-131C");
                break;
                
            case "TMTpro16plex":
                labels.add("TMTpro16plex-126");
                labels.add("TMTpro16plex-127N");
                labels.add("TMTpro16plex-127C");
                labels.add("TMTpro16plex-128N");
                labels.add("TMTpro16plex-128C");
                labels.add("TMTpro16plex-129N");
                labels.add("TMTpro16plex-129C");
                labels.add("TMTpro16plex-130N");
                labels.add("TMTpro16plex-130C");
                labels.add("TMTpro16plex-131N");
                labels.add("TMTpro16plex-131C");
                labels.add("TMTpro16plex-132N");
                labels.add("TMTpro16plex-132C");
                labels.add("TMTpro16plex-133N");
                labels.add("TMTpro16plex-133C");
                labels.add("TMTpro16plex-134N");
                break;
                
            case "TMTpro18plex":
                // Include all TMTpro16plex labels
                labels.addAll(getIsobaricLabelsForType("TMTpro16plex"));
                // Add TMTpro18plex specific labels
                labels.add("TMTpro18plex-134C");
                labels.add("TMTpro18plex-135N");
                break;
                
            default:
                // Fallback: try to find modifications that start with the selected base type
                for (String modification : isobaricModifications) {
                    if (modification.startsWith(baseType + "-")) {
                        labels.add(modification);
                    }
                }
                break;
        }
        
        return labels;
    }
    
    // Methods for isotopic labeling
    private void addLightLabel() {
        String selected = labelModificationComboBox.getValue();
        if (selected != null && !lightLabels.contains(selected)) {
            lightLabels.add(selected);
            updateParameterView();
            updateMetadataTable();
        }
    }
    
    private void addMediumLabel() {
        String selected = labelModificationComboBox.getValue();
        if (selected != null && !mediumLabels.contains(selected)) {
            mediumLabels.add(selected);
            updateParameterView();
            updateMetadataTable();
        }
    }
    
    private void addHeavyLabel() {
        String selected = labelModificationComboBox.getValue();
        if (selected != null && !heavyLabels.contains(selected)) {
            heavyLabels.add(selected);
            updateParameterView();
            updateMetadataTable();
        }
    }
    
    private void removeLightLabel() {
        String selected = lightLabelListView.getSelectionModel().getSelectedItem();
        if (selected != null) {
            lightLabels.remove(selected);
            updateParameterView();
            updateMetadataTable();
        }
    }
    
    private void removeMediumLabel() {
        String selected = mediumLabelListView.getSelectionModel().getSelectedItem();
        if (selected != null) {
            mediumLabels.remove(selected);
            updateParameterView();
            updateMetadataTable();
        }
    }
    
    private void removeHeavyLabel() {
        String selected = heavyLabelListView.getSelectionModel().getSelectedItem();
        if (selected != null) {
            heavyLabels.remove(selected);
            updateParameterView();
            updateMetadataTable();
        }
    }
    
    // Inner class for isobaric labeling table rows
    public static class IsobaricLabelRow {
        private final SimpleStringProperty labelName;
        private final SimpleStringProperty minus2;
        private final SimpleStringProperty minus1;
        private final SimpleStringProperty plus1;
        private final SimpleStringProperty plus2;
        
        public IsobaricLabelRow(String labelName) {
            this.labelName = new SimpleStringProperty(labelName);
            this.minus2 = new SimpleStringProperty("0.0");
            this.minus1 = new SimpleStringProperty("0.0");
            this.plus1 = new SimpleStringProperty("0.0");
            this.plus2 = new SimpleStringProperty("0.0");
        }
        
        public String getLabelName() { return labelName.get(); }
        public void setLabelName(String value) { labelName.set(value); }
        public StringProperty labelNameProperty() { return labelName; }
        
        public String getMinus2() { return minus2.get(); }
        public void setMinus2(String value) { minus2.set(value); }
        public StringProperty minus2Property() { return minus2; }
        
        public String getMinus1() { return minus1.get(); }
        public void setMinus1(String value) { minus1.set(value); }
        public StringProperty minus1Property() { return minus1; }
        
        public String getPlus1() { return plus1.get(); }
        public void setPlus1(String value) { plus1.set(value); }
        public StringProperty plus1Property() { return plus1; }
        
        public String getPlus2() { return plus2.get(); }
        public void setPlus2(String value) { plus2.set(value); }
        public StringProperty plus2Property() { return plus2; }
    }
    
    private void updateMetadataTableFromIOPanel() {
        // Update metadata table when Parameters tab is selected
        updateMetadataTable();
    }
    
    private void showAlert(String title, String message) {
        Alert alert = new Alert(AlertType.INFORMATION);
        alert.setTitle(title);
        alert.setHeaderText(null);
        alert.setContentText(message);
        alert.showAndWait();
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
            OperatingSystemMXBean osBean = ManagementFactory.getOperatingSystemMXBean();
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
} 