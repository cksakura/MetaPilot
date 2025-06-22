package uk.ac.quadram.metapilot.task.common.gui;

import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.control.*;
import javafx.scene.layout.*;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.scene.control.cell.PropertyValueFactory;
import javafx.scene.control.cell.TextFieldTableCell;
import javafx.util.converter.DefaultStringConverter;
import javafx.beans.property.SimpleStringProperty;
import javafx.beans.property.StringProperty;
import uk.ac.quadram.metapilot.core.mod.IsobaricTag;

/**
 * Reusable isobaric labeling panel that can be used across different MetaPilot workflows.
 * Handles isobaric type selection and correction factor configuration.
 */
public class IsobaricLabelingPanel extends VBox {
    // UI Components
    private ComboBox<String> isobaricTypeComboBox;
    private TableView<IsobaricLabelRow> isobaricTable;
    private ObservableList<IsobaricLabelRow> isobaricTableData = FXCollections.observableArrayList();
    private ObservableList<String> isobaricBaseTypes = FXCollections.observableArrayList();
    
    // Callback for parameter updates
    private Runnable parameterUpdateCallback;
    
    public IsobaricLabelingPanel() {
        initializeUI();
        loadIsobaricBaseTypes();
    }
    
    private void initializeUI() {
        setSpacing(8);
        setPadding(new Insets(8));
        setStyle("-fx-background-color: #f8e8f4; -fx-background-radius: 5;");
        setMinHeight(200);
        setPrefHeight(220);
        setMaxHeight(Double.MAX_VALUE);
        
        // Isobaric type selection
        HBox typeBox = new HBox(10);
        typeBox.setAlignment(Pos.CENTER_LEFT);
        Label typeLabel = new Label("Isobaric Type:");
        typeLabel.setStyle("-fx-font-weight: bold;");
        
        isobaricTypeComboBox = new ComboBox<>(isobaricBaseTypes);
        isobaricTypeComboBox.setMinWidth(200);
        isobaricTypeComboBox.setPrefWidth(220);
        isobaricTypeComboBox.setMaxWidth(Double.MAX_VALUE);
        isobaricTypeComboBox.setOnAction(e -> {
            updateIsobaricTable();
            notifyParameterUpdate();
        });
        HBox.setHgrow(isobaricTypeComboBox, Priority.ALWAYS);
        
        typeBox.getChildren().addAll(typeLabel, isobaricTypeComboBox);
        
        // Isobaric correction table
        isobaricTable = new TableView<>(isobaricTableData);
        isobaricTable.setMinHeight(120);
        isobaricTable.setPrefHeight(140);
        isobaricTable.setMaxHeight(Double.MAX_VALUE);
        isobaricTable.setEditable(true);
        isobaricTable.setColumnResizePolicy(TableView.CONSTRAINED_RESIZE_POLICY);
        VBox.setVgrow(isobaricTable, Priority.ALWAYS);
        
        createIsobaricTableColumns();
        
        getChildren().addAll(typeBox, isobaricTable);
    }
    
    private void createIsobaricTableColumns() {
        isobaricTable.getColumns().clear();
        
        TableColumn<IsobaricLabelRow, String> labelCol = new TableColumn<>("Label");
        labelCol.setCellValueFactory(new PropertyValueFactory<>("labelName"));
        labelCol.setMinWidth(120);
        labelCol.setPrefWidth(150);
        
        TableColumn<IsobaricLabelRow, String> minus2Col = new TableColumn<>("-2");
        minus2Col.setCellValueFactory(new PropertyValueFactory<>("minus2"));
        minus2Col.setCellFactory(TextFieldTableCell.forTableColumn(new DefaultStringConverter()));
        minus2Col.setOnEditCommit(e -> {
            e.getRowValue().setMinus2(e.getNewValue());
            notifyParameterUpdate();
        });
        minus2Col.setMinWidth(50);
        minus2Col.setPrefWidth(60);
        
        TableColumn<IsobaricLabelRow, String> minus1Col = new TableColumn<>("-1");
        minus1Col.setCellValueFactory(new PropertyValueFactory<>("minus1"));
        minus1Col.setCellFactory(TextFieldTableCell.forTableColumn(new DefaultStringConverter()));
        minus1Col.setOnEditCommit(e -> {
            e.getRowValue().setMinus1(e.getNewValue());
            notifyParameterUpdate();
        });
        minus1Col.setMinWidth(50);
        minus1Col.setPrefWidth(60);
        
        TableColumn<IsobaricLabelRow, String> plus1Col = new TableColumn<>("+1");
        plus1Col.setCellValueFactory(new PropertyValueFactory<>("plus1"));
        plus1Col.setCellFactory(TextFieldTableCell.forTableColumn(new DefaultStringConverter()));
        plus1Col.setOnEditCommit(e -> {
            e.getRowValue().setPlus1(e.getNewValue());
            notifyParameterUpdate();
        });
        plus1Col.setMinWidth(50);
        plus1Col.setPrefWidth(60);
        
        TableColumn<IsobaricLabelRow, String> plus2Col = new TableColumn<>("+2");
        plus2Col.setCellValueFactory(new PropertyValueFactory<>("plus2"));
        plus2Col.setCellFactory(TextFieldTableCell.forTableColumn(new DefaultStringConverter()));
        plus2Col.setOnEditCommit(e -> {
            e.getRowValue().setPlus2(e.getNewValue());
            notifyParameterUpdate();
        });
        plus2Col.setMinWidth(50);
        plus2Col.setPrefWidth(60);
        
        isobaricTable.getColumns().addAll(labelCol, minus2Col, minus1Col, plus1Col, plus2Col);
    }
    
    private void loadIsobaricBaseTypes() {
        IsobaricTag[] tags = IsobaricTag.getPFindTags();
        for (IsobaricTag tag : tags) {
            isobaricBaseTypes.add(tag.getName());
        }
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
                    isobaricTableData.add(new IsobaricLabelRow(title));
                }
            }
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
    public String getSelectedIsobaricType() {
        return isobaricTypeComboBox.getValue();
    }
    
    public ObservableList<IsobaricLabelRow> getIsobaricTableData() {
        return isobaricTableData;
    }
    
    public boolean hasLabels() {
        return !isobaricTableData.isEmpty();
    }
    
    /**
     * Inner class for isobaric labeling table rows.
     * Moved from MetaPilotDdaIgcMaxQuantPanel to make it reusable.
     */
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
} 