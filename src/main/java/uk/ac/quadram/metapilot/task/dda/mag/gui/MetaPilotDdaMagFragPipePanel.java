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
import javafx.scene.control.cell.PropertyValueFactory;
import javafx.scene.control.cell.TextFieldTableCell;
import javafx.scene.control.cell.CheckBoxTableCell;
import javafx.util.converter.DefaultStringConverter;
import javafx.util.converter.DoubleStringConverter;
import javafx.beans.property.SimpleStringProperty;
import javafx.beans.property.StringProperty;
import javafx.beans.property.SimpleBooleanProperty;
import javafx.beans.property.BooleanProperty;
import javafx.beans.property.SimpleDoubleProperty;
import javafx.beans.property.DoubleProperty;
import uk.ac.quadram.metapilot.task.common.gui.MetaPilotIOPanelFX;
import uk.ac.quadram.metapilot.task.common.gui.IsobaricLabelingPanel;
import uk.ac.quadram.metapilot.core.mod.IsobaricTag;
import uk.ac.quadram.metapilot.task.dda.igc.gui.MetaPilotDdaIgcMaxQuantPanel;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;
import uk.ac.quadram.metapilot.task.mag.par.MetaParameterMag;
import javafx.animation.Timeline;
import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.IOException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.io.File;

public class MetaPilotDdaMagFragPipePanel extends VBox {
    // UI Components - Database specific for MAG
    private MetaPilotIOPanelFX ioPanel;
    private ComboBox<MagDbItem> magDatabaseComboBox;
    private TextField hostDBField;
    private CheckBox appendHostDBCheckBox;
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

    // FragPipe-specific modification tables
    private TableView<FixedModificationRow> fixedModTable;
    private ObservableList<FixedModificationRow> fixedModData = FXCollections.observableArrayList();
    private Button restoreDefaultButton;
    
    private TableView<VariableModificationRow> variableModTable;
    private ObservableList<VariableModificationRow> variableModData = FXCollections.observableArrayList();

    // Additional parameters
    private Spinner<Integer> threadsSpinner;
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

    // Status panel fields
    private Label statusTimeLabel;
    private Label statusCpuLabel;
    private Label statusMemLabel;
    private Timeline statusTimeline;
    private long taskStartTimeMillis;

    // MAG database items
    private ObservableList<MagDbItem> availableMagDatabases = FXCollections.observableArrayList();

	public MetaPilotDdaMagFragPipePanel(MagDbItem[] magDbItems) {
		loadEnzymesFromIni();
		loadIsobaricBaseTypes();
		initializeModificationData();
		if (magDbItems != null) {
			for (MagDbItem item : magDbItems) {
				availableMagDatabases.add(item);
			}
		}
		initializeUI();
	}

    public MetaPilotDdaMagFragPipePanel(MetaParameterMag metaParameterMag) {
        loadEnzymesFromIni();
        loadIsobaricBaseTypes();
        initializeModificationData();
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

    // Create IO tab (similar to pFind)
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

        // Left column: Enzyme Settings, Additional Parameters, and Modifications Panel (NO taxonomic analysis for MAG)
        VBox leftColumn = new VBox(15);
        leftColumn.setMinWidth(430);
        VBox enzymeSection = createEnzymeSection();
        VBox additionalParamsSection = createAdditionalParametersSection();
        
        // Create horizontal modifications panel with fixed mods (narrow) and variable mods (wider)
        HBox modificationsPanel = createModificationsPanel();
        
        leftColumn.getChildren().addAll(enzymeSection, additionalParamsSection, modificationsPanel);
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
    private HBox createModificationsPanel() {
        HBox modificationsContainer = new HBox(10);
        modificationsContainer.setMinHeight(200);
        HBox.setHgrow(modificationsContainer, Priority.ALWAYS);
        
        // Fixed modifications section (narrow - only 2 columns)
        VBox fixedModSection = createFixedModificationsSection();
        fixedModSection.setMinWidth(180);
        fixedModSection.setPrefWidth(200);
        fixedModSection.setMaxWidth(220);
        
        // Variable modifications section (wider - 4 columns)
        VBox variableModSection = createVariableModificationsSection();
        HBox.setHgrow(variableModSection, Priority.ALWAYS);
        
        modificationsContainer.getChildren().addAll(fixedModSection, variableModSection);
        return modificationsContainer;
    }

    // FragPipe-specific Fixed Modifications Section
    private VBox createFixedModificationsSection() {
        VBox section = new VBox(10);
        section.setPadding(new Insets(10));
        section.setStyle("-fx-background-color: #fff8dc; -fx-background-radius: 5;");
        section.setMinHeight(200);
        VBox.setVgrow(section, Priority.ALWAYS);

        Label titleLabel = createSectionTitle("Fixed Modifications");
        
        // Create table - narrower since only 2 columns
        fixedModTable = new TableView<>(fixedModData);
        fixedModTable.setMinHeight(150);
        fixedModTable.setPrefHeight(180);
        fixedModTable.setMaxHeight(Double.MAX_VALUE);
        fixedModTable.setEditable(true);
        fixedModTable.setColumnResizePolicy(TableView.CONSTRAINED_RESIZE_POLICY);
        VBox.setVgrow(fixedModTable, Priority.ALWAYS);

        // Site column (read-only) - narrower
        TableColumn<FixedModificationRow, String> siteCol = new TableColumn<>("Site");
        siteCol.setCellValueFactory(new PropertyValueFactory<>("site"));
        siteCol.setMinWidth(60);
        siteCol.setPrefWidth(80);
        siteCol.setSortable(false);

        // Delta Mass column (editable) - narrower
        TableColumn<FixedModificationRow, Double> deltaMassCol = new TableColumn<>("Delta Mass");
        deltaMassCol.setCellValueFactory(new PropertyValueFactory<>("deltaMass"));
        deltaMassCol.setCellFactory(TextFieldTableCell.forTableColumn(new DoubleStringConverter()));
        deltaMassCol.setOnEditCommit(e -> {
            e.getRowValue().setDeltaMass(e.getNewValue());
            updateParameterView();
        });
        deltaMassCol.setMinWidth(80);
        deltaMassCol.setPrefWidth(100);
        deltaMassCol.setSortable(false);

        fixedModTable.getColumns().addAll(siteCol, deltaMassCol);

        // Restore default button
        restoreDefaultButton = new Button("Restore Default");
        restoreDefaultButton.setOnAction(e -> restoreDefaultFixedModifications());
        
        HBox buttonBox = new HBox();
        buttonBox.setAlignment(Pos.CENTER_RIGHT);
        buttonBox.getChildren().add(restoreDefaultButton);

        section.getChildren().addAll(titleLabel, fixedModTable, buttonBox);
        return section;
    }

    // FragPipe-specific Variable Modifications Section
    private VBox createVariableModificationsSection() {
        VBox section = new VBox(10);
        section.setPadding(new Insets(10));
        section.setStyle("-fx-background-color: #f0fff0; -fx-background-radius: 5;");
        section.setMinHeight(200);  // Reduced from 300 to match fixed mods
        VBox.setVgrow(section, Priority.ALWAYS);

        Label titleLabel = createSectionTitle("Variable Modifications");
        
        // Create table
        variableModTable = new TableView<>(variableModData);
        variableModTable.setMinHeight(150);  // Reduced to match fixed mods
        variableModTable.setPrefHeight(180);  // Reduced to match fixed mods
        variableModTable.setMaxHeight(Double.MAX_VALUE);
        variableModTable.setEditable(true);
        variableModTable.setColumnResizePolicy(TableView.CONSTRAINED_RESIZE_POLICY);
        VBox.setVgrow(variableModTable, Priority.ALWAYS);

        // Use column (checkbox) - smaller
        TableColumn<VariableModificationRow, Boolean> useCol = new TableColumn<>("Use");
        useCol.setCellValueFactory(new PropertyValueFactory<>("use"));
        useCol.setCellFactory(CheckBoxTableCell.forTableColumn(useCol));
        useCol.setOnEditCommit(e -> updateParameterView());
        useCol.setMinWidth(40);
        useCol.setPrefWidth(50);
        useCol.setSortable(false);

        // Mod Name column (editable) - adjusted width
        TableColumn<VariableModificationRow, String> modNameCol = new TableColumn<>("Modification");
        modNameCol.setCellValueFactory(new PropertyValueFactory<>("modName"));
        modNameCol.setCellFactory(TextFieldTableCell.forTableColumn());
        modNameCol.setOnEditCommit(e -> {
            e.getRowValue().setModName(e.getNewValue());
            updateParameterView();
        });
        modNameCol.setMinWidth(100);
        modNameCol.setPrefWidth(120);
        modNameCol.setSortable(false);

        // Sites column (editable) - smaller
        TableColumn<VariableModificationRow, String> sitesCol = new TableColumn<>("Sites");
        sitesCol.setCellValueFactory(new PropertyValueFactory<>("sites"));
        sitesCol.setCellFactory(TextFieldTableCell.forTableColumn());
        sitesCol.setOnEditCommit(e -> {
            e.getRowValue().setSites(e.getNewValue());
            updateParameterView();
        });
        sitesCol.setMinWidth(60);
        sitesCol.setPrefWidth(70);
        sitesCol.setSortable(false);

        // Delta Mass column (editable) - adjusted width
        TableColumn<VariableModificationRow, Double> deltaMassCol = new TableColumn<>("Delta Mass");
        deltaMassCol.setCellValueFactory(new PropertyValueFactory<>("deltaMass"));
        deltaMassCol.setCellFactory(TextFieldTableCell.forTableColumn(new DoubleStringConverter()));
        deltaMassCol.setOnEditCommit(e -> {
            e.getRowValue().setDeltaMass(e.getNewValue());
            updateParameterView();
        });
        deltaMassCol.setMinWidth(80);
        deltaMassCol.setPrefWidth(90);
        deltaMassCol.setSortable(false);

        variableModTable.getColumns().addAll(useCol, modNameCol, sitesCol, deltaMassCol);

        section.getChildren().addAll(titleLabel, variableModTable);
        return section;
    }

    private void restoreDefaultFixedModifications() {
        for (FixedModificationRow row : fixedModData) {
            if (row.getSite().equals("C")) {
                row.setDeltaMass(57.02146);  // Carbamidomethyl
            } else {
                row.setDeltaMass(0.0);
            }
        }
        fixedModTable.refresh();
        updateParameterView();
    }

    // Fixed Modification Row class for FragPipe
    public static class FixedModificationRow {
        private final SimpleStringProperty site;
        private final SimpleDoubleProperty deltaMass;
        
        public FixedModificationRow(String site, double deltaMass) {
            this.site = new SimpleStringProperty(site);
            this.deltaMass = new SimpleDoubleProperty(deltaMass);
        }
        
        public String getSite() { return site.get(); }
        public void setSite(String value) { site.set(value); }
        public StringProperty siteProperty() { return site; }
        
        public double getDeltaMass() { return deltaMass.get(); }
        public void setDeltaMass(double value) { deltaMass.set(value); }
        public DoubleProperty deltaMassProperty() { return deltaMass; }
    }

    // Variable Modification Row class for FragPipe
    public static class VariableModificationRow {
        private final SimpleBooleanProperty use;
        private final SimpleStringProperty modName;
        private final SimpleStringProperty sites;
        private final SimpleDoubleProperty deltaMass;
        
        public VariableModificationRow(boolean use, String modName, String sites, double deltaMass) {
            this.use = new SimpleBooleanProperty(use);
            this.modName = new SimpleStringProperty(modName);
            this.sites = new SimpleStringProperty(sites);
            this.deltaMass = new SimpleDoubleProperty(deltaMass);
        }
        
        public boolean getUse() { return use.get(); }
        public void setUse(boolean value) { use.set(value); }
        public BooleanProperty useProperty() { return use; }
        
        public String getModName() { return modName.get(); }
        public void setModName(String value) { modName.set(value); }
        public StringProperty modNameProperty() { return modName; }
        
        public String getSites() { return sites.get(); }
        public void setSites(String value) { sites.set(value); }
        public StringProperty sitesProperty() { return sites; }
        
        public double getDeltaMass() { return deltaMass.get(); }
        public void setDeltaMass(double value) { deltaMass.set(value); }
        public DoubleProperty deltaMassProperty() { return deltaMass; }
    }

    private void initializeModificationData() {
        // Initialize fixed modifications (all amino acids + termini)
        String[] allSites = {
            "peptide c-term", "peptide n-term", "protein c-term", "protein n-term",
            "A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"
        };
        
        for (String site : allSites) {
            double defaultMass = site.equals("C") ? 57.02146 : 0.0; // Carbamidomethyl on Cysteine
            fixedModData.add(new FixedModificationRow(site, defaultMass));
        }

        // Initialize variable modifications (15 rows)
        variableModData.addAll(
            // Common modifications
            new VariableModificationRow(false, "Oxidation", "M", 15.9949),
            new VariableModificationRow(false, "Phosphorylation", "STY", 79.9663),
            new VariableModificationRow(false, "Acetylation", "protein n-term", 42.0106),
            new VariableModificationRow(false, "Deamidation", "NQ", 0.9840),
            new VariableModificationRow(false, "Methylation", "KR", 14.0157),
            new VariableModificationRow(false, "Ubiquitination", "K", 114.0429),
            new VariableModificationRow(false, "Nitrosylation", "C", 28.9898),
            new VariableModificationRow(false, "Citrullination", "R", 0.9840),
            // User-defined slots
            new VariableModificationRow(false, "User Mod 1", "", 0.0),
            new VariableModificationRow(false, "User Mod 2", "", 0.0),
            new VariableModificationRow(false, "User Mod 3", "", 0.0),
            new VariableModificationRow(false, "User Mod 4", "", 0.0),
            new VariableModificationRow(false, "User Mod 5", "", 0.0),
            new VariableModificationRow(false, "User Mod 6", "", 0.0),
            new VariableModificationRow(false, "User Mod 7", "", 0.0)
        );
    }

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

    private void updateParameterView() {
        if (parameterView == null) {
            return;
        }
        
        StringBuilder sb = new StringBuilder();
        sb.append("=== FragPipe Analysis Parameters ===\n\n");

        // Database settings
        sb.append("Database Settings:\n");
        sb.append("  MAG Database: ").append(magDatabaseComboBox != null ? magDatabaseComboBox.getValue() != null ? magDatabaseComboBox.getValue().getCatalogueID() : "" : "").append("\n");
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

        // Fixed modifications
        sb.append("Fixed Modifications:\n");
        if (fixedModData != null) {
            for (FixedModificationRow row : fixedModData) {
                if (row.getDeltaMass() != 0.0) {
                    sb.append("  ").append(row.getSite()).append(": ").append(row.getDeltaMass()).append("\n");
                }
            }
        }
        sb.append("\n");

        // Variable modifications
        sb.append("Variable Modifications:\n");
        if (variableModData != null) {
            for (VariableModificationRow row : variableModData) {
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

        parameterView.setText(sb.toString());
    }

    private void updateMetadataTable() {
        if (metadataTable == null || ioPanel == null) {
            return;
        }
        
        metadataTableData.clear();
        ObservableList<MetaPilotIOPanelFX.ExperimentRow> fileRows = ioPanel.getExperimentData();
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
} 