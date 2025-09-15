package uk.ac.quadram.metapilot.task.common.gui;

import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.control.*;
import javafx.scene.layout.*;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.scene.text.Font;
import javafx.scene.text.FontWeight;
import javafx.scene.paint.Color;
import javafx.scene.control.cell.PropertyValueFactory;
import javafx.scene.control.cell.TextFieldTableCell;
import javafx.scene.control.cell.CheckBoxTableCell;
import javafx.util.converter.DoubleStringConverter;
import javafx.beans.property.SimpleStringProperty;
import javafx.beans.property.StringProperty;
import javafx.beans.property.SimpleBooleanProperty;
import javafx.beans.property.BooleanProperty;
import javafx.beans.property.SimpleDoubleProperty;
import javafx.beans.property.DoubleProperty;

/**
 * Reusable modifications panel for proteomics workflows
 * 
 * This panel provides fixed and variable modification tables that can be used
 * across different workflow panels (MAG/IGC with pFind/FragPipe/Sage).
 */
public class ModificationsPanel extends HBox {
    
    // Fixed modifications table and data
    private TableView<FixedModificationRow> fixedModTable;
    private ObservableList<FixedModificationRow> fixedModData = FXCollections.observableArrayList();
    private Button restoreDefaultButton;
    
    // Variable modifications table and data
    private TableView<VariableModificationRow> variableModTable;
    private ObservableList<VariableModificationRow> variableModData = FXCollections.observableArrayList();
    
    // Callback for parameter updates
    private Runnable parameterUpdateCallback;
    
    public ModificationsPanel() {
        this(null);
    }
    
    public ModificationsPanel(Runnable parameterUpdateCallback) {
        this.parameterUpdateCallback = parameterUpdateCallback;
        initializeModificationData();
        initializeUI();
    }
    
    private void initializeUI() {
        setSpacing(10);
        setMinHeight(200);
        HBox.setHgrow(this, Priority.ALWAYS);
        
        // Fixed modifications section (narrow - only 2 columns)
        VBox fixedModSection = createFixedModificationsSection();
        fixedModSection.setMinWidth(250);
        fixedModSection.setPrefWidth(270);
        fixedModSection.setMaxWidth(300);
        
        // Variable modifications section (wider - 4 columns)
        VBox variableModSection = createVariableModificationsSection();
        HBox.setHgrow(variableModSection, Priority.ALWAYS);
        
        getChildren().addAll(fixedModSection, variableModSection);
    }
    
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
            fixedModTable.refresh(); // Refresh to show the change
            updateParameterView();
        });
        deltaMassCol.setMinWidth(80);
        deltaMassCol.setPrefWidth(100);
        deltaMassCol.setSortable(false);
        deltaMassCol.setEditable(true);

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
    
    private VBox createVariableModificationsSection() {
        VBox section = new VBox(10);
        section.setPadding(new Insets(10));
        section.setStyle("-fx-background-color: #f0fff0; -fx-background-radius: 5;");
        section.setMinHeight(200);
        VBox.setVgrow(section, Priority.ALWAYS);

        Label titleLabel = createSectionTitle("Variable Modifications");
        
        // Create table
        variableModTable = new TableView<>(variableModData);
        variableModTable.setMinHeight(150);
        variableModTable.setPrefHeight(180);
        variableModTable.setMaxHeight(Double.MAX_VALUE);
        variableModTable.setEditable(true);
        variableModTable.setColumnResizePolicy(TableView.CONSTRAINED_RESIZE_POLICY);
        VBox.setVgrow(variableModTable, Priority.ALWAYS);

        // Use column (checkbox) - smaller
        TableColumn<VariableModificationRow, Boolean> useCol = new TableColumn<>("Use");
        useCol.setCellValueFactory(new PropertyValueFactory<>("use"));
        useCol.setCellFactory(CheckBoxTableCell.forTableColumn(useCol));
        useCol.setOnEditCommit(e -> {
            e.getRowValue().setUse(e.getNewValue());
            variableModTable.refresh();
            updateParameterView();
        });
        useCol.setMinWidth(40);
        useCol.setPrefWidth(50);
        useCol.setSortable(false);
        useCol.setEditable(true);

        // Mod Name column (editable only for User Mod rows) - adjusted width
        TableColumn<VariableModificationRow, String> modNameCol = new TableColumn<>("Modification");
        modNameCol.setCellValueFactory(new PropertyValueFactory<>("modName"));
        modNameCol.setCellFactory(column -> new TextFieldTableCell<VariableModificationRow, String>() {
            @Override
            public void startEdit() {
                VariableModificationRow row = getTableRow().getItem();
                if (row != null && row.getModName().startsWith("User Mod")) {
                    super.startEdit();
                }
            }
        });
        modNameCol.setOnEditCommit(e -> {
            e.getRowValue().setModName(e.getNewValue());
            variableModTable.refresh();
            updateParameterView();
        });
        modNameCol.setMinWidth(100);
        modNameCol.setPrefWidth(120);
        modNameCol.setSortable(false);
        modNameCol.setEditable(true);

        // Sites column (editable only for User Mod rows) - smaller
        TableColumn<VariableModificationRow, String> sitesCol = new TableColumn<>("Sites");
        sitesCol.setCellValueFactory(new PropertyValueFactory<>("sites"));
        sitesCol.setCellFactory(column -> new TextFieldTableCell<VariableModificationRow, String>() {
            @Override
            public void startEdit() {
                VariableModificationRow row = getTableRow().getItem();
                if (row != null && row.getModName().startsWith("User Mod")) {
                    super.startEdit();
                }
            }
        });
        sitesCol.setOnEditCommit(e -> {
            e.getRowValue().setSites(e.getNewValue());
            variableModTable.refresh();
            updateParameterView();
        });
        sitesCol.setMinWidth(60);
        sitesCol.setPrefWidth(70);
        sitesCol.setSortable(false);
        sitesCol.setEditable(true);

        // Delta Mass column (editable only for User Mod rows) - adjusted width
        TableColumn<VariableModificationRow, Double> deltaMassCol = new TableColumn<>("Delta Mass");
        deltaMassCol.setCellValueFactory(new PropertyValueFactory<>("deltaMass"));
        deltaMassCol.setCellFactory(column -> new TextFieldTableCell<VariableModificationRow, Double>(new DoubleStringConverter()) {
            @Override
            public void startEdit() {
                VariableModificationRow row = getTableRow().getItem();
                if (row != null && row.getModName().startsWith("User Mod")) {
                    super.startEdit();
                }
            }
        });
        deltaMassCol.setOnEditCommit(e -> {
            e.getRowValue().setDeltaMass(e.getNewValue());
            variableModTable.refresh();
            updateParameterView();
        });
        deltaMassCol.setMinWidth(80);
        deltaMassCol.setPrefWidth(90);
        deltaMassCol.setSortable(false);
        deltaMassCol.setEditable(true);

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
    
    private void updateParameterView() {
        if (parameterUpdateCallback != null) {
            parameterUpdateCallback.run();
        }
    }
    
    private Label createSectionTitle(String title) {
        Label label = new Label(title);
        label.setFont(Font.font("System", FontWeight.BOLD, 14));
        label.setTextFill(Color.DARKBLUE);
        return label;
    }
    
    private void initializeModificationData() {
        // Clear existing data first to prevent duplicates
        fixedModData.clear();
        variableModData.clear();
        
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
    
    // Getters for the modification data
    public ObservableList<FixedModificationRow> getFixedModData() {
        return fixedModData;
    }
    
    public ObservableList<VariableModificationRow> getVariableModData() {
        return variableModData;
    }
    
    // Setter for parameter update callback
    public void setParameterUpdateCallback(Runnable callback) {
        this.parameterUpdateCallback = callback;
    }
    
    // Fixed Modification Row class
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

    // Variable Modification Row class
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
}