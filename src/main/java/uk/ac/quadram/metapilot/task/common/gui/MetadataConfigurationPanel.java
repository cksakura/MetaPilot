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
import javafx.beans.property.SimpleStringProperty;
import uk.ac.quadram.metapilot.task.dda.igc.gui.MetaPilotDdaIgcMaxQuantPanel;

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
    
    // Callback for parameter updates
    private Runnable parameterUpdateCallback;
    
    public MetadataConfigurationPanel() {
        initializeUI();
    }
    
    private void initializeUI() {
        setSpacing(10);
        setPadding(new Insets(10));
        setStyle("-fx-background-color: #fff8f0; -fx-background-radius: 5;");
        setMinHeight(200);
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
        
        controlBox.getChildren().addAll(
            metadataColumnComboBox, 
            customColumnTextField, 
            addMetadataColumnButton, 
            removeMetadataColumnButton
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
        sampleCol.setCellValueFactory(new PropertyValueFactory<>("sampleName"));
        sampleCol.setMinWidth(120);
        sampleCol.setPrefWidth(150);
        metadataTable.getColumns().add(sampleCol);
    }
    
    public void setParameterUpdateCallback(Runnable callback) {
        this.parameterUpdateCallback = callback;
    }
    
    public void updateMetadataFromIOPanel(ObservableList<MetaLabIOPanelFX.ExperimentRow> experiments) {
        metadataTableData.clear();
        for (MetaLabIOPanelFX.ExperimentRow exp : experiments) {
            metadataTableData.add(new MetaPilotDdaIgcMaxQuantPanel.MetadataRow(exp.getExperimentName()));
        }
    }
    
    public void updateMetadataForQuantMode(String quantMode, boolean hasLabels, 
                                          ObservableList<MetaLabIOPanelFX.ExperimentRow> fileRows,
                                          ObservableList<IsobaricLabelingPanel.IsobaricLabelRow> isobaricTableData) {
        metadataTableData.clear();
        if (fileRows.isEmpty()) return;
        
        createBasicMetadataTableColumns();
        
        if (quantMode.contains("Isobaric") && hasLabels) {
            for (MetaLabIOPanelFX.ExperimentRow row : fileRows) {
                for (IsobaricLabelingPanel.IsobaricLabelRow isoLabel : isobaricTableData) {
                    MetaPilotDdaIgcMaxQuantPanel.MetadataRow metadataRow = 
                        new MetaPilotDdaIgcMaxQuantPanel.MetadataRow(row.getExperimentName());
                    metadataRow.setLabelType("Isobaric - " + isoLabel.getLabelName());
                    metadataTableData.add(metadataRow);
                }
            }
        } else {
            for (MetaLabIOPanelFX.ExperimentRow row : fileRows) {
                metadataTableData.add(new MetaPilotDdaIgcMaxQuantPanel.MetadataRow(row.getExperimentName()));
            }
        }
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
} 