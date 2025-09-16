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
import javafx.scene.control.cell.PropertyValueFactory;
import javafx.scene.control.cell.TextFieldTableCell;
import javafx.util.converter.DefaultStringConverter;
import javafx.beans.property.SimpleStringProperty;
import uk.ac.quadram.metapilot.task.dda.igc.gui.MetaPilotDdaIgcMaxQuantPanel;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Reusable metadata configuration panel that can be used across different MetaPilot workflows.
 * Handles metadata table configuration with dynamic column management.
 */
public class MetadataConfigurationPanel extends VBox {
    // UI Components
    private TableView<MetaPilotDdaIgcMaxQuantPanel.MetadataRow> metadataTable;
    private ObservableList<MetaPilotDdaIgcMaxQuantPanel.MetadataRow> metadataTableData = FXCollections.observableArrayList();
    private ComboBox<String> metadataColumnComboBox;
    private TextField customColumnTextField;
    private Button addMetadataColumnButton;
    private Button removeMetadataColumnButton;
    private Button resetMetadataButton;
    
    // Callback for parameter updates
    private Runnable parameterUpdateCallback;
    
    public MetadataConfigurationPanel() {
        initializeUI();
    }
    
    private void initializeUI() {
        setSpacing(10);
        setPadding(new Insets(10));
        setStyle("-fx-background-color: #fff8f0; -fx-background-radius: 5;");
        setMinHeight(240);
        VBox.setVgrow(this, Priority.ALWAYS);
        
        Label titleLabel = createSectionTitle("Metadata Configuration");
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
        
        createBasicMetadataTableColumns();
        
        // Control panel
        HBox controlBox = new HBox(8);
        controlBox.setAlignment(Pos.CENTER_LEFT);
        
        metadataColumnComboBox = new ComboBox<>(FXCollections.observableArrayList(
            "Age", "Gender", "Condition", "Treatment", "Time Point", "BMI", "Custom"));
        metadataColumnComboBox.setMinWidth(100);
        metadataColumnComboBox.setPrefWidth(120);
        metadataColumnComboBox.setPromptText("Add column");
        
        customColumnTextField = new TextField();
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
        addMetadataColumnButton.setOnAction(e -> addMetadataColumn());
        
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
        
        getChildren().addAll(titleLabel, infoLabel, metadataTable, controlBox);
    }
    
    private Label createSectionTitle(String title) {
        Label label = new Label(title);
        label.setFont(Font.font("System", FontWeight.BOLD, 14));
        label.setTextFill(Color.DARKBLUE);
        return label;
    }
    
    private void createBasicMetadataTableColumns() {
        metadataTable.getColumns().clear();
        TableColumn<MetaPilotDdaIgcMaxQuantPanel.MetadataRow, String> sampleCol = new TableColumn<>("Sample");
        sampleCol.setCellValueFactory(cellData -> cellData.getValue().sampleNameProperty());
        sampleCol.setCellFactory(column -> {
            TextFieldTableCell<MetaPilotDdaIgcMaxQuantPanel.MetadataRow, String> cell = 
                new TextFieldTableCell<MetaPilotDdaIgcMaxQuantPanel.MetadataRow, String>(new DefaultStringConverter()) {
                @Override
                public void startEdit() {
                    super.startEdit();
                    TextField textField = (TextField) getGraphic();
                    if (textField != null) {
                        // Commit edit when focus is lost (like Swing)
                        textField.focusedProperty().addListener((obs, wasFocused, isNowFocused) -> {
                            if (wasFocused && !isNowFocused && isEditing()) {
                                commitEdit(textField.getText());
                            }
                        });
                    }
                }
            };
            return cell;
        });
        sampleCol.setOnEditCommit(e -> {
            e.getRowValue().setSampleName(e.getNewValue());
            notifyParameterUpdate();
        });
        sampleCol.setMinWidth(120);
        sampleCol.setPrefWidth(150);
        metadataTable.getColumns().add(sampleCol);
    }
    
    public void setParameterUpdateCallback(Runnable callback) {
        this.parameterUpdateCallback = callback;
    }
    
    public void updateMetadataFromIOPanel(ObservableList<MetaPilotIOPanelFX.ExperimentRow> experiments) {
        // Check if experiment names have changed before clearing existing data
        boolean experimentsChanged = false;
        
        // Check if the number of experiments changed
        if (metadataTableData.size() != experiments.size()) {
            experimentsChanged = true;
        } else {
            // Check if any experiment names changed
            for (int i = 0; i < experiments.size(); i++) {
                String expName = experiments.get(i).getExperimentName();
                if (i >= metadataTableData.size() || !metadataTableData.get(i).getSampleName().equals(expName)) {
                    experimentsChanged = true;
                    break;
                }
            }
        }
        
        // Only recreate if experiments changed
        if (experimentsChanged) {
            // Save existing metadata before clearing
            Map<String, MetaPilotDdaIgcMaxQuantPanel.MetadataRow> existingData = new HashMap<>();
            for (MetaPilotDdaIgcMaxQuantPanel.MetadataRow row : metadataTableData) {
                existingData.put(row.getSampleName(), row);
            }
            
            metadataTableData.clear();
            for (MetaPilotIOPanelFX.ExperimentRow exp : experiments) {
                String expName = exp.getExperimentName();
                MetaPilotDdaIgcMaxQuantPanel.MetadataRow metadataRow;
                
                // Try to restore existing data
                if (existingData.containsKey(expName)) {
                    metadataRow = existingData.get(expName);
                } else {
                    metadataRow = new MetaPilotDdaIgcMaxQuantPanel.MetadataRow(expName);
                }
                metadataTableData.add(metadataRow);
            }
        }
        // If experiments haven't changed, preserve existing metadata
    }
    
    public void updateMetadataForQuantMode(String quantMode, boolean hasLabels, 
                                          ObservableList<MetaPilotIOPanelFX.ExperimentRow> fileRows,
                                          ObservableList<IsobaricLabelingPanel.IsobaricLabelRow> isobaricTableData) {
        if (fileRows.isEmpty()) return;
        
        // Calculate what the new structure should be
        List<String> expectedSampleNames = new ArrayList<>();
        if (quantMode.contains("Isobaric") && hasLabels) {
            for (MetaPilotIOPanelFX.ExperimentRow row : fileRows) {
                for (IsobaricLabelingPanel.IsobaricLabelRow isoLabel : isobaricTableData) {
                    String sampleName = row.getExperimentName() + "_" + isoLabel.getLabelName();
                    expectedSampleNames.add(sampleName);
                }
            }
        } else {
            for (MetaPilotIOPanelFX.ExperimentRow row : fileRows) {
                expectedSampleNames.add(row.getExperimentName());
            }
        }
        
        // Check if the current structure matches the expected structure
        boolean structureChanged = false;
        if (metadataTableData.size() != expectedSampleNames.size()) {
            structureChanged = true;
        } else {
            for (int i = 0; i < expectedSampleNames.size(); i++) {
                if (!metadataTableData.get(i).getSampleName().equals(expectedSampleNames.get(i))) {
                    structureChanged = true;
                    break;
                }
            }
        }
        
        // Only recreate if structure changed, otherwise preserve existing data
        if (structureChanged) {
            // Save existing metadata values before clearing
            Map<String, MetaPilotDdaIgcMaxQuantPanel.MetadataRow> existingData = new HashMap<>();
            for (MetaPilotDdaIgcMaxQuantPanel.MetadataRow row : metadataTableData) {
                existingData.put(row.getSampleName(), row);
            }
            
            metadataTableData.clear();
            createBasicMetadataTableColumns();
            
            // Create new rows and restore data where possible
            if (quantMode.contains("Isobaric") && hasLabels) {
                for (MetaPilotIOPanelFX.ExperimentRow row : fileRows) {
                    for (IsobaricLabelingPanel.IsobaricLabelRow isoLabel : isobaricTableData) {
                        String sampleName = row.getExperimentName() + "_" + isoLabel.getLabelName();
                        MetaPilotDdaIgcMaxQuantPanel.MetadataRow metadataRow;
                        
                        // Try to restore existing data
                        if (existingData.containsKey(sampleName)) {
                            metadataRow = existingData.get(sampleName);
                        } else {
                            metadataRow = new MetaPilotDdaIgcMaxQuantPanel.MetadataRow(sampleName);
                            metadataRow.setLabelType("Isobaric - " + isoLabel.getLabelName());
                        }
                        metadataTableData.add(metadataRow);
                    }
                }
            } else {
                for (MetaPilotIOPanelFX.ExperimentRow row : fileRows) {
                    String sampleName = row.getExperimentName();
                    MetaPilotDdaIgcMaxQuantPanel.MetadataRow metadataRow;
                    
                    // Try to restore existing data
                    if (existingData.containsKey(sampleName)) {
                        metadataRow = existingData.get(sampleName);
                    } else {
                        metadataRow = new MetaPilotDdaIgcMaxQuantPanel.MetadataRow(sampleName);
                    }
                    metadataTableData.add(metadataRow);
                }
            }
        }
        // If structure hasn't changed, do nothing to preserve user edits
    }
    
    private void addMetadataColumn() {
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
    }
    
    private void addMetadataColumnWithProperty(String columnName, String propertyName) {
        TableColumn<MetaPilotDdaIgcMaxQuantPanel.MetadataRow, String> newCol = new TableColumn<>(columnName);
        newCol.setCellValueFactory(cellData -> {
            switch (propertyName) {
                case "age": return cellData.getValue().ageProperty();
                case "gender": return cellData.getValue().genderProperty();
                case "condition": return cellData.getValue().conditionProperty();
                case "treatment": return cellData.getValue().treatmentProperty();
                case "timePoint": return cellData.getValue().timePointProperty();
                case "customField1": return cellData.getValue().customField1Property();
                case "customField2": return cellData.getValue().customField2Property();
                case "customField3": return cellData.getValue().customField3Property();
                default: return new SimpleStringProperty("");
            }
        });
        newCol.setCellFactory(column -> {
            TextFieldTableCell<MetaPilotDdaIgcMaxQuantPanel.MetadataRow, String> cell = 
                new TextFieldTableCell<MetaPilotDdaIgcMaxQuantPanel.MetadataRow, String>(new DefaultStringConverter()) {
                @Override
                public void startEdit() {
                    super.startEdit();
                    TextField textField = (TextField) getGraphic();
                    if (textField != null) {
                        // Commit edit when focus is lost (like Swing)
                        textField.focusedProperty().addListener((obs, wasFocused, isNowFocused) -> {
                            if (wasFocused && !isNowFocused && isEditing()) {
                                commitEdit(textField.getText());
                            }
                        });
                    }
                }
            };
            return cell;
        });
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
            notifyParameterUpdate();
        });
        newCol.setMinWidth(100);
        newCol.setPrefWidth(120);
        metadataTable.getColumns().add(newCol);
        notifyParameterUpdate();
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
    
    private void removeMetadataColumn() {
        TableColumn<MetaPilotDdaIgcMaxQuantPanel.MetadataRow, ?> selectedColumn = 
            metadataTable.getSelectionModel().getSelectedCells().isEmpty() ? null :
            metadataTable.getColumns().get(metadataTable.getSelectionModel().getSelectedCells().get(0).getColumn());
        
        if (selectedColumn != null && !selectedColumn.getText().equals("Sample")) {
            metadataTable.getColumns().remove(selectedColumn);
            notifyParameterUpdate();
        }
    }
    
    private void notifyParameterUpdate() {
        if (parameterUpdateCallback != null) {
            parameterUpdateCallback.run();
        }
    }
    
    public TableView<MetaPilotDdaIgcMaxQuantPanel.MetadataRow> getMetadataTable() {
        return metadataTable;
    }
    
    public ObservableList<MetaPilotDdaIgcMaxQuantPanel.MetadataRow> getMetadataTableData() {
        return metadataTableData;
    }
    
    public int getMetadataColumnCount() {
        // Return number of metadata columns (excluding Sample and Label Type columns)
        return Math.max(0, metadataTable.getColumns().size() - 1); // Subtract 1 for Sample column
    }
    
    // Setter for parameter loading
    public void setMetadataTableData(List<MetaPilotDdaIgcMaxQuantPanel.MetadataRow> data) {
        metadataTableData.clear();
        if (data != null) {
            metadataTableData.addAll(data);
        }
        notifyParameterUpdate();
    }
    
    /**
     * Resets the metadata table to its default state (only Sample column)
     */
    private void resetMetadataTable() {
        // Clear all data
        metadataTableData.clear();
        
        // Remove all columns except the Sample column
        metadataTable.getColumns().clear();
        
        // Recreate basic table structure with only Sample column
        createBasicMetadataTableColumns();
        
        // Notify that parameters have changed
        notifyParameterUpdate();
    }
} 