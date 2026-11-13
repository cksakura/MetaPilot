package uk.ac.quadram.metapilot.task.common.gui;

import javafx.geometry.Pos;
import javafx.scene.control.*;
import javafx.scene.layout.*;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.IOException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Reusable pFind modifications panel for proteomics workflows
 * 
 * This panel provides fixed and variable modification ListViews that can be used
 * across different pFind workflow panels (MAG/IGC with pFind).
 */
public class PFindModificationPanel extends VBox {
    
    // pFind modification components (ListView-based, not table-based)
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
    
    public PFindModificationPanel() {
        this(null);
    }
    
    public PFindModificationPanel(Runnable parameterUpdateCallback) {
        this.parameterUpdateCallback = parameterUpdateCallback;
        loadModificationsFromIni();
        initializeUI();
    }
    
    private void initializeUI() {
        VBox.setVgrow(this, Priority.ALWAYS);
        
        GridPane modGrid = new GridPane();
        modGrid.setHgap(15);
        modGrid.setVgap(10);
        
        // Column constraints
        ColumnConstraints col1 = new ColumnConstraints();
        col1.setMinWidth(200);
        col1.setHgrow(Priority.ALWAYS);
        ColumnConstraints col2 = new ColumnConstraints();
        col2.setMinWidth(200);
        col2.setHgrow(Priority.ALWAYS);
        modGrid.getColumnConstraints().addAll(col1, col2);
        
        // Row constraints
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
        
        // Modification selection and add buttons
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
        
        TitledPane section = ParameterSection.create("Modifications", modGrid);
        VBox.setVgrow(section, Priority.ALWAYS);
        getChildren().add(section);
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
    
    private void addFixedModification() {
        String selected = modificationComboBox.getValue();
        if (selected != null && !fixedModifications.contains(selected)) {
            fixedModifications.add(selected);
            if (parameterUpdateCallback != null) {
                parameterUpdateCallback.run();
            }
        }
    }
    
    private void addVariableModification() {
        String selected = modificationComboBox.getValue();
        if (selected != null && !variableModifications.contains(selected)) {
            variableModifications.add(selected);
            if (parameterUpdateCallback != null) {
                parameterUpdateCallback.run();
            }
        }
    }
    
    private void removeFixedModification() {
        String selected = fixedModListView.getSelectionModel().getSelectedItem();
        if (selected != null) {
            fixedModifications.remove(selected);
            if (parameterUpdateCallback != null) {
                parameterUpdateCallback.run();
            }
        }
    }
    
    private void removeVariableModification() {
        String selected = variableModListView.getSelectionModel().getSelectedItem();
        if (selected != null) {
            variableModifications.remove(selected);
            if (parameterUpdateCallback != null) {
                parameterUpdateCallback.run();
            }
        }
    }
    
    
    // Getters for the modification data
    public ObservableList<String> getFixedModifications() {
        return fixedModifications;
    }
    
    public ObservableList<String> getVariableModifications() {
        return variableModifications;
    }
    
    // Setters for loading saved modifications
    public void setFixedModifications(java.util.List<String> modifications) {
        fixedModifications.clear();
        if (modifications != null) {
            fixedModifications.addAll(modifications);
        }
    }
    
    public void setVariableModifications(java.util.List<String> modifications) {
        variableModifications.clear();
        if (modifications != null) {
            variableModifications.addAll(modifications);
        }
    }
}
