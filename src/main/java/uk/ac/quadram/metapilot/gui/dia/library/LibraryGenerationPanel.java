package uk.ac.quadram.metapilot.gui.dia.library;

import javafx.geometry.Insets;
import javafx.scene.control.*;
import javafx.scene.layout.*;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.beans.property.*;
import javafx.stage.FileChooser;
import java.io.File;

public class LibraryGenerationPanel extends VBox {
    private ToggleGroup inputTypeGroup;
    private VBox libraryInputPanel;
    private VBox peptideListInputPanel;
    private VBox proteinDbInputPanel;
    private TableView<LibraryItem> libraryTable;
    private Button checkAvailableButton;
    private Button deleteButton;

    public LibraryGenerationPanel() {
        setPadding(new Insets(10));
        setSpacing(10);

        // Input type selection
        TitledPane inputTypePane = new TitledPane("Input Type", createInputTypePanel());
        inputTypePane.setExpanded(true);

        // Library table
        TitledPane libraryTablePane = new TitledPane("Available Libraries", createLibraryTablePanel());
        libraryTablePane.setExpanded(true);

        getChildren().addAll(inputTypePane, libraryTablePane);
    }

    private VBox createInputTypePanel() {
        VBox panel = new VBox(10);
        panel.setPadding(new Insets(5));

        // Input type selection
        inputTypeGroup = new ToggleGroup();
        RadioButton libraryButton = new RadioButton("Library");
        RadioButton peptideListButton = new RadioButton("Peptide List");
        RadioButton proteinDbButton = new RadioButton("Protein Database");

        libraryButton.setToggleGroup(inputTypeGroup);
        peptideListButton.setToggleGroup(inputTypeGroup);
        proteinDbButton.setToggleGroup(inputTypeGroup);

        libraryButton.setSelected(true);

        // Input panels
        libraryInputPanel = createLibraryInputPanel();
        peptideListInputPanel = createPeptideListInputPanel();
        proteinDbInputPanel = createProteinDbInputPanel();

        // Show/hide panels based on selection
        libraryButton.selectedProperty().addListener((obs, oldVal, newVal) -> {
            libraryInputPanel.setVisible(newVal);
            peptideListInputPanel.setVisible(false);
            proteinDbInputPanel.setVisible(false);
        });

        peptideListButton.selectedProperty().addListener((obs, oldVal, newVal) -> {
            libraryInputPanel.setVisible(false);
            peptideListInputPanel.setVisible(newVal);
            proteinDbInputPanel.setVisible(false);
        });

        proteinDbButton.selectedProperty().addListener((obs, oldVal, newVal) -> {
            libraryInputPanel.setVisible(false);
            peptideListInputPanel.setVisible(false);
            proteinDbInputPanel.setVisible(newVal);
        });

        panel.getChildren().addAll(
            libraryButton,
            libraryInputPanel,
            peptideListButton,
            peptideListInputPanel,
            proteinDbButton,
            proteinDbInputPanel
        );

        return panel;
    }

    private VBox createLibraryInputPanel() {
        VBox panel = new VBox(10);
        panel.setPadding(new Insets(5));

        HBox fileBox = new HBox(10);
        Label fileLabel = new Label("Library File:");
        TextField fileField = new TextField();
        Button browseButton = new Button("Browse");
        browseButton.setOnAction(e -> {
            FileChooser fileChooser = new FileChooser();
            fileChooser.setTitle("Select Library File");
            File file = fileChooser.showOpenDialog(getScene().getWindow());
            if (file != null) {
                fileField.setText(file.getAbsolutePath());
            }
        });
        fileBox.getChildren().addAll(fileLabel, fileField, browseButton);

        panel.getChildren().add(fileBox);
        return panel;
    }

    private VBox createPeptideListInputPanel() {
        VBox panel = new VBox(10);
        panel.setPadding(new Insets(5));

        HBox fileBox = new HBox(10);
        Label fileLabel = new Label("Peptide List File:");
        TextField fileField = new TextField();
        Button browseButton = new Button("Browse");
        browseButton.setOnAction(e -> {
            FileChooser fileChooser = new FileChooser();
            fileChooser.setTitle("Select Peptide List File");
            File file = fileChooser.showOpenDialog(getScene().getWindow());
            if (file != null) {
                fileField.setText(file.getAbsolutePath());
            }
        });
        fileBox.getChildren().addAll(fileLabel, fileField, browseButton);

        panel.getChildren().add(fileBox);
        return panel;
    }

    private VBox createProteinDbInputPanel() {
        VBox panel = new VBox(10);
        panel.setPadding(new Insets(5));

        // File selection
        HBox fileBox = new HBox(10);
        Label fileLabel = new Label("Protein Database File:");
        TextField fileField = new TextField();
        Button browseButton = new Button("Browse");
        browseButton.setOnAction(e -> {
            FileChooser fileChooser = new FileChooser();
            fileChooser.setTitle("Select Protein Database File");
            File file = fileChooser.showOpenDialog(getScene().getWindow());
            if (file != null) {
                fileField.setText(file.getAbsolutePath());
            }
        });
        fileBox.getChildren().addAll(fileLabel, fileField, browseButton);

        // Enzyme parameters
        HBox enzymeBox = new HBox(10);
        Label enzymeLabel = new Label("Enzyme:");
        ComboBox<String> enzymeComboBox = new ComboBox<>();
        enzymeComboBox.getItems().addAll("Trypsin", "Lys-C", "Arg-C", "Glu-C", "Chymotrypsin");
        enzymeComboBox.setValue("Trypsin");
        enzymeBox.getChildren().addAll(enzymeLabel, enzymeComboBox);

        HBox missedCleavagesBox = new HBox(10);
        Label missedCleavagesLabel = new Label("Missed Cleavages:");
        Spinner<Integer> missedCleavagesSpinner = new Spinner<>(0, 10, 2);
        missedCleavagesBox.getChildren().addAll(missedCleavagesLabel, missedCleavagesSpinner);

        // Keep protein names option
        HBox keepNamesBox = new HBox(10);
        Label keepNamesLabel = new Label("Keep Protein Names in Library:");
        ComboBox<String> keepNamesComboBox = new ComboBox<>();
        keepNamesComboBox.getItems().addAll("Yes", "No");
        keepNamesComboBox.setValue("No");
        keepNamesBox.getChildren().addAll(keepNamesLabel, keepNamesComboBox);

        panel.getChildren().addAll(fileBox, enzymeBox, missedCleavagesBox, keepNamesBox);
        return panel;
    }

    private VBox createLibraryTablePanel() {
        VBox panel = new VBox(10);
        panel.setPadding(new Insets(5));

        // Create table
        libraryTable = new TableView<>();
        
        // Name column
        TableColumn<LibraryItem, String> nameColumn = new TableColumn<>("Name");
        nameColumn.setCellValueFactory(cellData -> cellData.getValue().nameProperty());
        
        // Type column
        TableColumn<LibraryItem, String> typeColumn = new TableColumn<>("Type");
        typeColumn.setCellValueFactory(cellData -> cellData.getValue().typeProperty());
        
        // Path column
        TableColumn<LibraryItem, String> pathColumn = new TableColumn<>("Path");
        pathColumn.setCellValueFactory(cellData -> cellData.getValue().pathProperty());
        
        libraryTable.getColumns().addAll(nameColumn, typeColumn, pathColumn);

        // Buttons
        HBox buttonBox = new HBox(10);
        checkAvailableButton = new Button("Check Available");
        deleteButton = new Button("Delete");
        
        checkAvailableButton.setOnAction(e -> checkAvailableLibraries());
        deleteButton.setOnAction(e -> deleteSelectedLibrary());
        
        buttonBox.getChildren().addAll(checkAvailableButton, deleteButton);

        panel.getChildren().addAll(libraryTable, buttonBox);
        return panel;
    }

    private void checkAvailableLibraries() {
        // TODO: Implement library detection logic
    }

    private void deleteSelectedLibrary() {
        LibraryItem selectedItem = libraryTable.getSelectionModel().getSelectedItem();
        if (selectedItem != null) {
            Alert alert = new Alert(Alert.AlertType.CONFIRMATION);
            alert.setTitle("Delete Library");
            alert.setHeaderText("Delete Library");
            alert.setContentText("How would you like to delete this library?");

            ButtonType deleteFromTableButton = new ButtonType("Delete from Table Only");
            ButtonType deleteFromDiskButton = new ButtonType("Delete from Disk");
            ButtonType cancelButton = ButtonType.CANCEL;

            alert.getButtonTypes().setAll(deleteFromTableButton, deleteFromDiskButton, cancelButton);

            alert.showAndWait().ifPresent(response -> {
                if (response == deleteFromTableButton) {
                    libraryTable.getItems().remove(selectedItem);
                } else if (response == deleteFromDiskButton) {
                    // TODO: Implement file deletion logic
                    libraryTable.getItems().remove(selectedItem);
                }
            });
        }
    }

    // Helper class for library items
    public static class LibraryItem {
        private final StringProperty name;
        private final StringProperty type;
        private final StringProperty path;

        public LibraryItem(String name, String type, String path) {
            this.name = new SimpleStringProperty(name);
            this.type = new SimpleStringProperty(type);
            this.path = new SimpleStringProperty(path);
        }

        public StringProperty nameProperty() { return name; }
        public StringProperty typeProperty() { return type; }
        public StringProperty pathProperty() { return path; }
    }
} 