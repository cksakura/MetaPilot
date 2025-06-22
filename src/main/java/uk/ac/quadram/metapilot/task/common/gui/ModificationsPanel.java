package uk.ac.quadram.metapilot.task.common.gui;

import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.control.*;
import javafx.scene.layout.*;
import javafx.scene.paint.Color;
import javafx.scene.text.Font;
import javafx.scene.text.FontWeight;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.IOException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Reusable modifications panel that can be used across different MetaPilot workflows.
 * Handles fixed and variable modification management with loading from INI files.
 */
public class ModificationsPanel extends VBox {
    // UI Components
    private ListView<String> fixedModListView;
    private ListView<String> variableModListView;
    private ComboBox<String> modificationComboBox;
    private Button addFixedModButton;
    private Button addVariableModButton;
    private Button removeFixedModButton;
    private Button removeVariableModButton;
    
    // Data
    private ObservableList<String> modifications = FXCollections.observableArrayList();
    private ObservableList<String> fixedModifications = FXCollections.observableArrayList();
    private ObservableList<String> variableModifications = FXCollections.observableArrayList();
    
    // Callback for parameter updates
    private Runnable parameterUpdateCallback;
    
    public ModificationsPanel() {
        loadModificationsFromIni();
        initializeUI();
    }
    
    private void initializeUI() {
        setSpacing(10);
        setPadding(new Insets(10));
        setStyle("-fx-background-color: #f5f5f5; -fx-background-radius: 5;");
        setMinHeight(250);
        VBox.setVgrow(this, Priority.ALWAYS);
        
        Label titleLabel = createSectionTitle("Modifications");
        
        // Create grid layout for modifications
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
        
        // Fixed modifications
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
        
        // Variable modifications
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
        
        // Modification selection
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
        modGrid.add(modSelectionContainer, 0, 3, 2, 1);
        
        getChildren().addAll(titleLabel, modGrid);
    }
    
    private Label createSectionTitle(String title) {
        Label label = new Label(title);
        label.setFont(Font.font("System", FontWeight.BOLD, 14));
        label.setTextFill(Color.DARKBLUE);
        return label;
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
            // Add some default modifications if file loading fails
            modifications.addAll(
                "Oxidation (M)",
                "Carbamidomethyl (C)",
                "Acetyl (Protein N-term)",
                "Phospho (STY)",
                "Deamidated (NQ)"
            );
        }
    }
    
    private void addFixedModification() {
        String selectedMod = modificationComboBox.getValue();
        if (selectedMod != null && !fixedModifications.contains(selectedMod)) {
            fixedModifications.add(selectedMod);
            notifyParameterUpdate();
        }
    }
    
    private void addVariableModification() {
        String selectedMod = modificationComboBox.getValue();
        if (selectedMod != null && !variableModifications.contains(selectedMod)) {
            variableModifications.add(selectedMod);
            notifyParameterUpdate();
        }
    }
    
    private void removeFixedModification() {
        String selectedMod = fixedModListView.getSelectionModel().getSelectedItem();
        if (selectedMod != null) {
            fixedModifications.remove(selectedMod);
            notifyParameterUpdate();
        }
    }
    
    private void removeVariableModification() {
        String selectedMod = variableModListView.getSelectionModel().getSelectedItem();
        if (selectedMod != null) {
            variableModifications.remove(selectedMod);
            notifyParameterUpdate();
        }
    }
    
    public void setParameterUpdateCallback(Runnable callback) {
        this.parameterUpdateCallback = callback;
    }
    
    private void notifyParameterUpdate() {
        if (parameterUpdateCallback != null) {
            parameterUpdateCallback.run();
        }
    }
    
    // Getters for external access
    public ObservableList<String> getAvailableModifications() {
        return modifications;
    }
    
    public ObservableList<String> getFixedModifications() {
        return fixedModifications;
    }
    
    public ObservableList<String> getVariableModifications() {
        return variableModifications;
    }
    
    // Add common modifications programmatically
    public void addCommonFixedModifications() {
        if (!fixedModifications.contains("Carbamidomethyl (C)")) {
            fixedModifications.add("Carbamidomethyl (C)");
        }
        notifyParameterUpdate();
    }
    
    public void addCommonVariableModifications() {
        String[] common = {"Oxidation (M)", "Acetyl (Protein N-term)", "Deamidated (NQ)"};
        for (String mod : common) {
            if (!variableModifications.contains(mod)) {
                variableModifications.add(mod);
            }
        }
        notifyParameterUpdate();
    }
    
    public void clearAllModifications() {
        fixedModifications.clear();
        variableModifications.clear();
        notifyParameterUpdate();
    }
} 