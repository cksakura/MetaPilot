package uk.ac.quadram.metapilot.gui.dia.model;

import javafx.geometry.Insets;
import javafx.scene.control.*;
import javafx.scene.layout.*;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.beans.property.*;
import javafx.stage.FileChooser;
import java.io.File;

public class ModelGenerationPanel extends VBox {
    private TableView<ModelItem> modelTable;
    private Button checkAvailableButton;
    private Button deleteButton;
    private Label hintLabel;

    public ModelGenerationPanel() {
        setPadding(new Insets(10));
        setSpacing(10);

        // Input panel
        TitledPane inputPane = new TitledPane("Model Generation", createInputPanel());
        inputPane.setExpanded(true);

        // Model table
        TitledPane modelTablePane = new TitledPane("Available Models", createModelTablePanel());
        modelTablePane.setExpanded(true);

        // Hint label
        hintLabel = new Label("Note: Model is optional in the data processing workflow");
        hintLabel.setStyle("-fx-text-fill: gray;");

        getChildren().addAll(inputPane, modelTablePane, hintLabel);
    }

    private VBox createInputPanel() {
        VBox panel = new VBox(10);
        panel.setPadding(new Insets(5));

        HBox fileBox = new HBox(10);
        Label fileLabel = new Label("Peptide Intensity Table:");
        TextField fileField = new TextField();
        Button browseButton = new Button("Browse");
        browseButton.setOnAction(e -> {
            FileChooser fileChooser = new FileChooser();
            fileChooser.setTitle("Select Peptide Intensity Table");
            File file = fileChooser.showOpenDialog(getScene().getWindow());
            if (file != null) {
                fileField.setText(file.getAbsolutePath());
            }
        });
        fileBox.getChildren().addAll(fileLabel, fileField, browseButton);

        panel.getChildren().add(fileBox);
        return panel;
    }

    private VBox createModelTablePanel() {
        VBox panel = new VBox(10);
        panel.setPadding(new Insets(5));

        // Create table
        modelTable = new TableView<>();
        
        // Name column
        TableColumn<ModelItem, String> nameColumn = new TableColumn<>("Name");
        nameColumn.setCellValueFactory(cellData -> cellData.getValue().nameProperty());
        
        // Type column
        TableColumn<ModelItem, String> typeColumn = new TableColumn<>("Type");
        typeColumn.setCellValueFactory(cellData -> cellData.getValue().typeProperty());
        
        // Path column
        TableColumn<ModelItem, String> pathColumn = new TableColumn<>("Path");
        pathColumn.setCellValueFactory(cellData -> cellData.getValue().pathProperty());
        
        modelTable.getColumns().addAll(nameColumn, typeColumn, pathColumn);

        // Buttons
        HBox buttonBox = new HBox(10);
        checkAvailableButton = new Button("Check Available");
        deleteButton = new Button("Delete");
        
        checkAvailableButton.setOnAction(e -> checkAvailableModels());
        deleteButton.setOnAction(e -> deleteSelectedModel());
        
        buttonBox.getChildren().addAll(checkAvailableButton, deleteButton);

        panel.getChildren().addAll(modelTable, buttonBox);
        return panel;
    }

    private void checkAvailableModels() {
        // TODO: Implement model detection logic
    }

    private void deleteSelectedModel() {
        ModelItem selectedItem = modelTable.getSelectionModel().getSelectedItem();
        if (selectedItem != null) {
            Alert alert = new Alert(Alert.AlertType.CONFIRMATION);
            alert.setTitle("Delete Model");
            alert.setHeaderText("Delete Model");
            alert.setContentText("How would you like to delete this model?");

            ButtonType deleteFromTableButton = new ButtonType("Delete from Table Only");
            ButtonType deleteFromDiskButton = new ButtonType("Delete from Disk");
            ButtonType cancelButton = ButtonType.CANCEL;

            alert.getButtonTypes().setAll(deleteFromTableButton, deleteFromDiskButton, cancelButton);

            alert.showAndWait().ifPresent(response -> {
                if (response == deleteFromTableButton) {
                    modelTable.getItems().remove(selectedItem);
                } else if (response == deleteFromDiskButton) {
                    // TODO: Implement file deletion logic
                    modelTable.getItems().remove(selectedItem);
                }
            });
        }
    }

    // Helper class for model items
    public static class ModelItem {
        private final StringProperty name;
        private final StringProperty type;
        private final StringProperty path;

        public ModelItem(String name, String type, String path) {
            this.name = new SimpleStringProperty(name);
            this.type = new SimpleStringProperty(type);
            this.path = new SimpleStringProperty(path);
        }

        public StringProperty nameProperty() { return name; }
        public StringProperty typeProperty() { return type; }
        public StringProperty pathProperty() { return path; }
    }
} 