package uk.ac.quadram.metapilot.gui.dda.mag;

import javafx.geometry.Insets;
import javafx.scene.control.*;
import javafx.scene.layout.*;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.beans.property.*;
import uk.ac.quadram.metapilot.gui.common.IOPanel;
import uk.ac.quadram.metapilot.gui.common.MetadataPanel;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

public class MetaPilotDdaMagPFindPanel extends VBox {
    private TabPane tabPane;
    private Tab ioTab;
    private Tab parameterTab;
    private Tab consoleTab;
    private ComboBox<String> enzymeComboBox;
    private TableView<FixedModification> fixedModTable;
    private TableView<VariableModification> variableModTable;
    private ToggleGroup quantModeGroup;
    private VBox isobaricPanel;
    private ComboBox<String> isobaricTypeComboBox;
    private TableView<String> isobaricTable;

    public MetaPilotDdaMagPFindPanel() {
        setPadding(new Insets(10));
        setSpacing(10);

        // Create tab pane
        tabPane = new TabPane();

        // Create tabs
        ioTab = new Tab("Input/Output");
        parameterTab = new Tab("Parameters");
        consoleTab = new Tab("Console");

        // Set up IO tab with reusable IOPanel
        ioTab.setContent(new IOPanel());

        // Set up parameter tab with MetadataPanel and other parameter sections
        VBox parameterPanel = new VBox(10);
        parameterPanel.setPadding(new Insets(10));
        
        // Metadata panel (reusable)
        MetadataPanel metadataPanel = new MetadataPanel();
        parameterPanel.getChildren().add(metadataPanel);

        // Enzyme section
        TitledPane enzymePane = new TitledPane("Enzyme", createEnzymePanel());
        enzymePane.setExpanded(true);
        parameterPanel.getChildren().add(enzymePane);

        // Modifications section
        HBox modificationsBox = new HBox(10);
        modificationsBox.setPadding(new Insets(5));
        
        // Fixed modifications
        TitledPane fixedModPane = new TitledPane("Fixed Modifications", createFixedModificationsPanel());
        fixedModPane.setExpanded(true);
        
        // Variable modifications
        TitledPane variableModPane = new TitledPane("Variable Modifications", createVariableModificationsPanel());
        variableModPane.setExpanded(true);
        
        modificationsBox.getChildren().addAll(fixedModPane, variableModPane);
        parameterPanel.getChildren().add(modificationsBox);

        // Quantitative section
        TitledPane quantitativePane = new TitledPane("Quantitative", createQuantitativePanel());
        quantitativePane.setExpanded(true);
        parameterPanel.getChildren().add(quantitativePane);

        // MAG specific parameters
        TitledPane magParamsPane = new TitledPane("MAG Parameters", createMagParametersPanel());
        magParamsPane.setExpanded(true);
        parameterPanel.getChildren().add(magParamsPane);

        parameterTab.setContent(parameterPanel);

        // Set up console tab
        consoleTab.setContent(createConsolePanel());

        // Add tabs to tab pane
        tabPane.getTabs().addAll(ioTab, parameterTab, consoleTab);

        // Add tab pane to main layout
        getChildren().add(tabPane);

        // Load enzymes and modifications
        loadEnzymesFromIni();
        initializeModifications();
    }

    private VBox createEnzymePanel() {
        VBox panel = new VBox(5);
        panel.setPadding(new Insets(5));

        // Enzyme selection
        HBox enzymeBox = new HBox(10);
        Label enzymeLabel = new Label("Enzyme:");
        enzymeComboBox = new ComboBox<>();
        enzymeBox.getChildren().addAll(enzymeLabel, enzymeComboBox);

        // Missed cleavages
        HBox missedCleavagesBox = new HBox(10);
        Label missedCleavagesLabel = new Label("Missed Cleavages:");
        Spinner<Integer> missedCleavagesSpinner = new Spinner<>(0, 10, 2);
        missedCleavagesBox.getChildren().addAll(missedCleavagesLabel, missedCleavagesSpinner);

        panel.getChildren().addAll(enzymeBox, missedCleavagesBox);
        return panel;
    }

    private VBox createFixedModificationsPanel() {
        VBox panel = new VBox(5);
        panel.setPadding(new Insets(5));

        // Create table for fixed modifications
        fixedModTable = new TableView<>();
        
        // Site column
        TableColumn<FixedModification, String> siteColumn = new TableColumn<>("Site");
        siteColumn.setCellValueFactory(cellData -> cellData.getValue().siteProperty());
        
        // Delta mass column
        TableColumn<FixedModification, Double> deltaMassColumn = new TableColumn<>("Delta Mass");
        deltaMassColumn.setCellValueFactory(cellData -> cellData.getValue().deltaMassProperty().asObject());
        
        fixedModTable.getColumns().addAll(siteColumn, deltaMassColumn);
        
        // Add default sites
        ObservableList<FixedModification> fixedMods = FXCollections.observableArrayList(
            new FixedModification("Peptide C-term", 0.0),
            new FixedModification("Peptide N-term", 0.0),
            new FixedModification("Protein C-term", 0.0),
            new FixedModification("Protein N-term", 0.0),
            new FixedModification("C", 57.02146),
            new FixedModification("D", 0.0),
            new FixedModification("E", 0.0),
            new FixedModification("F", 0.0),
            new FixedModification("G", 0.0),
            new FixedModification("H", 0.0),
            new FixedModification("I", 0.0),
            new FixedModification("K", 0.0),
            new FixedModification("L", 0.0),
            new FixedModification("M", 0.0),
            new FixedModification("N", 0.0),
            new FixedModification("P", 0.0),
            new FixedModification("Q", 0.0),
            new FixedModification("R", 0.0),
            new FixedModification("S", 0.0),
            new FixedModification("T", 0.0),
            new FixedModification("V", 0.0),
            new FixedModification("W", 0.0),
            new FixedModification("Y", 0.0)
        );
        
        fixedModTable.setItems(fixedMods);
        
        // Add restore defaults button
        Button restoreDefaultsButton = new Button("Restore Defaults");
        restoreDefaultsButton.setOnAction(e -> {
            for (FixedModification mod : fixedMods) {
                if (mod.getSite().equals("C")) {
                    mod.setDeltaMass(57.02146);
                } else {
                    mod.setDeltaMass(0.0);
                }
            }
        });
        
        panel.getChildren().addAll(fixedModTable, restoreDefaultsButton);
        return panel;
    }

    private VBox createVariableModificationsPanel() {
        VBox panel = new VBox(5);
        panel.setPadding(new Insets(5));

        // Create table for variable modifications
        variableModTable = new TableView<>();
        
        // Mod name column
        TableColumn<VariableModification, String> modNameColumn = new TableColumn<>("Modification");
        modNameColumn.setCellValueFactory(cellData -> cellData.getValue().modNameProperty());
        
        // Sites column
        TableColumn<VariableModification, String> sitesColumn = new TableColumn<>("Sites");
        sitesColumn.setCellValueFactory(cellData -> cellData.getValue().sitesProperty());
        
        // Delta mass column
        TableColumn<VariableModification, Double> deltaMassColumn = new TableColumn<>("Delta Mass");
        deltaMassColumn.setCellValueFactory(cellData -> cellData.getValue().deltaMassProperty().asObject());
        
        // Use column
        TableColumn<VariableModification, Boolean> useColumn = new TableColumn<>("Use");
        useColumn.setCellValueFactory(cellData -> cellData.getValue().useProperty());
        
        variableModTable.getColumns().addAll(modNameColumn, sitesColumn, deltaMassColumn, useColumn);
        
        // Add default modifications
        ObservableList<VariableModification> variableMods = FXCollections.observableArrayList(
            new VariableModification("Oxidation", "M", 15.9949, true),
            new VariableModification("Phosphorylation", "STY", 79.9663, false),
            new VariableModification("Acetylation", "K", 42.0106, false),
            new VariableModification("Methylation", "KR", 14.0157, false),
            new VariableModification("GlyGly", "K", 114.0429, false),
            new VariableModification("Carbamidomethyl", "C", 57.0215, false),
            new VariableModification("Deamidation", "NQ", 0.9840, false),
            new VariableModification("Pyro-glu", "QE", -17.0265, false),
            new VariableModification("Pyro-carbamidomethyl", "C", 39.9949, false),
            new VariableModification("Sulfo", "STY", 79.9568, false),
            new VariableModification("Gln->pyro-Glu", "Q", -17.0265, false),
            new VariableModification("Asn->pyro-Glu", "N", -17.0265, false),
            new VariableModification("Met->Hsl", "M", -48.0034, false),
            new VariableModification("Met->Hse", "M", -30.0106, false),
            new VariableModification("Met->Hsl", "M", -48.0034, false)
        );
        
        variableModTable.setItems(variableMods);
        
        panel.getChildren().add(variableModTable);
        return panel;
    }

    private VBox createQuantitativePanel() {
        VBox panel = new VBox(5);
        panel.setPadding(new Insets(5));

        // Quantitative mode selection
        quantModeGroup = new ToggleGroup();
        RadioButton labelFreeButton = new RadioButton("Label Free");
        RadioButton isobaricButton = new RadioButton("Isobaric");

        labelFreeButton.setToggleGroup(quantModeGroup);
        isobaricButton.setToggleGroup(quantModeGroup);

        labelFreeButton.setSelected(true);

        // Isobaric panel
        isobaricPanel = new VBox(5);
        isobaricPanel.setPadding(new Insets(5));
        isobaricPanel.setVisible(false);

        // Isobaric type selection
        HBox isobaricTypeBox = new HBox(10);
        Label isobaricTypeLabel = new Label("Isobaric Type:");
        isobaricTypeComboBox = new ComboBox<>();
        isobaricTypeComboBox.getItems().addAll("iTRAQ4plex", "iTRAQ8plex", "TMT6plex", "TMT10plex", "TMT11plex", "TMTpro16plex", "TMTpro18plex");
        isobaricTypeBox.getChildren().addAll(isobaricTypeLabel, isobaricTypeComboBox);

        // Isobaric table
        isobaricTable = new TableView<>();
        TableColumn<String, String> channelColumn = new TableColumn<>("Channel");
        TableColumn<String, String> correctionColumn = new TableColumn<>("Correction");
        isobaricTable.getColumns().addAll(channelColumn, correctionColumn);

        isobaricPanel.getChildren().addAll(isobaricTypeBox, isobaricTable);

        // Add listeners
        isobaricButton.selectedProperty().addListener((obs, oldVal, newVal) -> {
            isobaricPanel.setVisible(newVal);
        });

        panel.getChildren().addAll(labelFreeButton, isobaricButton, isobaricPanel);
        return panel;
    }

    private VBox createMagParametersPanel() {
        VBox panel = new VBox(5);
        panel.setPadding(new Insets(5));

        // MAG specific parameters
        HBox minPeptidesBox = new HBox(10);
        Label minPeptidesLabel = new Label("Minimum Peptides per Protein:");
        Spinner<Integer> minPeptidesSpinner = new Spinner<>(1, 10, 2);
        minPeptidesBox.getChildren().addAll(minPeptidesLabel, minPeptidesSpinner);

        HBox minLengthBox = new HBox(10);
        Label minLengthLabel = new Label("Minimum Peptide Length:");
        Spinner<Integer> minLengthSpinner = new Spinner<>(5, 50, 7);
        minLengthBox.getChildren().addAll(minLengthLabel, minLengthSpinner);

        HBox maxLengthBox = new HBox(10);
        Label maxLengthLabel = new Label("Maximum Peptide Length:");
        Spinner<Integer> maxLengthSpinner = new Spinner<>(5, 50, 30);
        maxLengthBox.getChildren().addAll(maxLengthLabel, maxLengthSpinner);

        // MAG specific features
        HBox matchBetweenRunsBox = new HBox(10);
        Label matchBetweenRunsLabel = new Label("Match Between Runs:");
        ComboBox<String> matchBetweenRunsComboBox = new ComboBox<>();
        matchBetweenRunsComboBox.getItems().addAll("Yes", "No");
        matchBetweenRunsComboBox.setValue("Yes");
        matchBetweenRunsBox.getChildren().addAll(matchBetweenRunsLabel, matchBetweenRunsComboBox);

        HBox iBAQBox = new HBox(10);
        Label iBAQLabel = new Label("iBAQ Quantification:");
        ComboBox<String> iBAQComboBox = new ComboBox<>();
        iBAQComboBox.getItems().addAll("Yes", "No");
        iBAQComboBox.setValue("Yes");
        iBAQBox.getChildren().addAll(iBAQLabel, iBAQComboBox);

        panel.getChildren().addAll(
            minPeptidesBox,
            minLengthBox,
            maxLengthBox,
            matchBetweenRunsBox,
            iBAQBox
        );
        return panel;
    }

    private VBox createConsolePanel() {
        VBox panel = new VBox(10);
        panel.setPadding(new Insets(10));

        // Parameter view
        TitledPane parameterViewPane = new TitledPane("Parameter View", createParameterViewPanel());
        parameterViewPane.setExpanded(true);

        // Progress section
        TitledPane progressPane = new TitledPane("Progress", createProgressPanel());
        progressPane.setExpanded(true);

        // Console output
        TitledPane consolePane = new TitledPane("Console Output", createConsoleOutputPanel());
        consolePane.setExpanded(true);

        panel.getChildren().addAll(parameterViewPane, progressPane, consolePane);
        return panel;
    }

    private VBox createParameterViewPanel() {
        VBox panel = new VBox(5);
        panel.setPadding(new Insets(5));

        // Parameter view text area
        TextArea parameterViewArea = new TextArea();
        parameterViewArea.setEditable(false);
        parameterViewArea.setPrefRowCount(5);

        panel.getChildren().add(parameterViewArea);
        return panel;
    }

    private VBox createProgressPanel() {
        VBox panel = new VBox(5);
        panel.setPadding(new Insets(5));

        // Progress bars
        ProgressBar mainProgressBar = new ProgressBar(0);
        ProgressBar subProgressBar = new ProgressBar(0);

        panel.getChildren().addAll(mainProgressBar, subProgressBar);
        return panel;
    }

    private VBox createConsoleOutputPanel() {
        VBox panel = new VBox(5);
        panel.setPadding(new Insets(5));

        // Console output text area
        TextArea consoleOutputArea = new TextArea();
        consoleOutputArea.setEditable(false);
        consoleOutputArea.setPrefRowCount(10);

        panel.getChildren().add(consoleOutputArea);
        return panel;
    }

    private void loadEnzymesFromIni() {
        try {
            List<String> enzymes = new ArrayList<>();
            String iniPath = "config/enzyme_default.ini";
            if (Files.exists(Paths.get(iniPath))) {
                try (BufferedReader reader = new BufferedReader(new FileReader(iniPath))) {
                    String line;
                    while ((line = reader.readLine()) != null) {
                        if (line.startsWith("name=")) {
                            String enzymeName = line.substring(5).trim();
                            enzymes.add(enzymeName);
                        }
                    }
                }
            }
            enzymeComboBox.getItems().addAll(enzymes);
            if (!enzymes.isEmpty()) {
                enzymeComboBox.setValue(enzymes.get(0));
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private void initializeModifications() {
        // Fixed modifications are initialized in createFixedModificationsPanel
        // Variable modifications are initialized in createVariableModificationsPanel
    }

    // Helper classes for modifications
    public static class FixedModification {
        private final StringProperty site;
        private final DoubleProperty deltaMass;

        public FixedModification(String site, double deltaMass) {
            this.site = new SimpleStringProperty(site);
            this.deltaMass = new SimpleDoubleProperty(deltaMass);
        }

        public StringProperty siteProperty() { return site; }
        public DoubleProperty deltaMassProperty() { return deltaMass; }
        public String getSite() { return site.get(); }
        public double getDeltaMass() { return deltaMass.get(); }
        public void setDeltaMass(double value) { deltaMass.set(value); }
    }

    public static class VariableModification {
        private final StringProperty modName;
        private final StringProperty sites;
        private final DoubleProperty deltaMass;
        private final BooleanProperty use;

        public VariableModification(String modName, String sites, double deltaMass, boolean use) {
            this.modName = new SimpleStringProperty(modName);
            this.sites = new SimpleStringProperty(sites);
            this.deltaMass = new SimpleDoubleProperty(deltaMass);
            this.use = new SimpleBooleanProperty(use);
        }

        public StringProperty modNameProperty() { return modName; }
        public StringProperty sitesProperty() { return sites; }
        public DoubleProperty deltaMassProperty() { return deltaMass; }
        public BooleanProperty useProperty() { return use; }
        public String getModName() { return modName.get(); }
        public String getSites() { return sites.get(); }
        public double getDeltaMass() { return deltaMass.get(); }
        public boolean isUse() { return use.get(); }
        public void setUse(boolean value) { use.set(value); }
    }
} 