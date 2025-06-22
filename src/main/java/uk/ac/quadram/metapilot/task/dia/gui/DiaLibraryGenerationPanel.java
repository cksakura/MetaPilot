package uk.ac.quadram.metapilot.task.dia.gui;

import javafx.application.Platform;
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
import javafx.scene.control.cell.CheckBoxTableCell;
import javafx.beans.property.SimpleStringProperty;
import javafx.beans.property.StringProperty;
import javafx.beans.property.SimpleBooleanProperty;
import javafx.beans.property.BooleanProperty;
import javafx.concurrent.Task;
import javafx.scene.control.Alert.AlertType;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;
import uk.ac.quadram.metapilot.task.mag.MagDbItem.MagPepLib;
import java.io.File;
import java.util.Optional;
import javafx.animation.KeyFrame;
import javafx.animation.Timeline;
import javafx.util.Duration;

public class DiaLibraryGenerationPanel extends VBox {
    
    // Library table components
    private TableView<LibraryRow> libraryTable;
    private ObservableList<LibraryRow> libraryData = FXCollections.observableArrayList();
    
    // Generation components
    private ToggleGroup inputTypeGroup;
    private TextField inputFileField;
    private ComboBox<String> enzymeComboBox;
    private ComboBox<Integer> missCleavageComboBox;
    private CheckBox keepProteinNamesCheckBox;
    private TextField libraryNameField;
    private VBox enzymeBox;
    
    // Current MAG database reference
    private MagDbItem currentMagDbItem;
    
    // Console output for feedback
    private TextArea consoleOutput;
    
    // Internal progress tracking
    private ProgressBar libraryProgressBar;
    private Label progressLabel;
    
    // Callback interface for parameter updates
    public interface ParameterUpdateCallback {
        void onParametersChanged();
    }
    
    private ParameterUpdateCallback parameterCallback;

    public DiaLibraryGenerationPanel() {
        initializeUI();
    }
    
    public DiaLibraryGenerationPanel(TextArea consoleOutput) {
        this.consoleOutput = consoleOutput;
        initializeUI();
    }
    
    public void setParameterUpdateCallback(ParameterUpdateCallback callback) {
        this.parameterCallback = callback;
    }
    
    public void setConsoleOutput(TextArea consoleOutput) {
        this.consoleOutput = consoleOutput;
    }
    
    public void setMagDbItem(MagDbItem magDbItem) {
        this.currentMagDbItem = magDbItem;
        loadAvailableLibraries();
    }

    private void initializeUI() {
        setSpacing(10);
        setPadding(new Insets(10));
        setStyle("-fx-background-color: #fff8dc; -fx-background-radius: 5;");
        setMinHeight(400);
        VBox.setVgrow(this, Priority.ALWAYS);

        Label titleLabel = createSectionTitle("Spectral Library Generation");
        
        // Input type selection
        VBox inputTypeBox = createInputTypeSelection();
        
        // File selection
        HBox fileBox = createFileSelection();
        
        // Enzyme parameters (for protein database option)
        enzymeBox = createEnzymeParameters();
        
        // Generation controls
        HBox generationBox = createGenerationControls();
        
        // Available libraries table
        VBox tableSection = createLibraryTable();
        
        // Progress and console section
        VBox progressSection = createProgressSection();
        
        getChildren().addAll(titleLabel, inputTypeBox, fileBox, enzymeBox, 
            generationBox, tableSection, progressSection);
    }
    
    private VBox createInputTypeSelection() {
        VBox inputTypeBox = new VBox(8);
        Label inputLabel = new Label("Input Type for Library Generation:");
        inputLabel.setStyle("-fx-font-weight: bold;");
        
        inputTypeGroup = new ToggleGroup();
        RadioButton libraryRadio = new RadioButton("Existing Library File");
        RadioButton peptideListRadio = new RadioButton("Peptide List File");
        RadioButton proteinDbRadio = new RadioButton("Protein Database File");
        
        libraryRadio.setToggleGroup(inputTypeGroup);
        peptideListRadio.setToggleGroup(inputTypeGroup);
        proteinDbRadio.setToggleGroup(inputTypeGroup);
        libraryRadio.setSelected(true); // Default selection
        
        // Show/hide enzyme parameters based on input type
        inputTypeGroup.selectedToggleProperty().addListener((obs, oldToggle, newToggle) -> {
            boolean showEnzyme = newToggle == proteinDbRadio;
            enzymeBox.setVisible(showEnzyme);
            enzymeBox.setManaged(showEnzyme);
        });
        
        inputTypeBox.getChildren().addAll(inputLabel, libraryRadio, peptideListRadio, proteinDbRadio);
        return inputTypeBox;
    }
    
    private HBox createFileSelection() {
        HBox fileBox = new HBox(10);
        fileBox.setAlignment(Pos.CENTER_LEFT);
        Label fileLabel = new Label("Input File:");
        fileLabel.setStyle("-fx-font-weight: bold;");
        inputFileField = new TextField();
        inputFileField.setPromptText("Select input file");
        inputFileField.setPrefWidth(250);
        inputFileField.setEditable(false);
        Button browseFileButton = new Button("Browse");
        browseFileButton.setOnAction(e -> browseInputFile());
        fileBox.getChildren().addAll(fileLabel, inputFileField, browseFileButton);
        return fileBox;
    }
    
    private VBox createEnzymeParameters() {
        VBox enzymeBox = new VBox(8);
        enzymeBox.setStyle("-fx-background-color: #f0f0f0; -fx-padding: 10; -fx-background-radius: 3;");
        Label enzymeTitle = new Label("Enzyme Parameters (for Protein Database):");
        enzymeTitle.setStyle("-fx-font-weight: bold;");
        
        GridPane enzymeGrid = new GridPane();
        enzymeGrid.setHgap(10);
        enzymeGrid.setVgap(5);
        
        // Enzyme name
        Label enzymeNameLabel = new Label("Enzyme:");
        enzymeComboBox = new ComboBox<>(FXCollections.observableArrayList(
            "Trypsin", "Trypsin/P", "Lys-C", "Arg-C", "Asp-N", "Glu-C", "Pepsin"));
        enzymeComboBox.setValue("Trypsin");
        enzymeComboBox.setPrefWidth(120);
        
        // Miss cleavage
        Label missCleavageLabel = new Label("Max Miss Cleavages:");
        missCleavageComboBox = new ComboBox<>();
        for (int i = 0; i <= 5; i++) {
            missCleavageComboBox.getItems().add(i);
        }
        missCleavageComboBox.setValue(2);
        missCleavageComboBox.setPrefWidth(80);
        
        enzymeGrid.add(enzymeNameLabel, 0, 0);
        enzymeGrid.add(enzymeComboBox, 1, 0);
        enzymeGrid.add(missCleavageLabel, 0, 1);
        enzymeGrid.add(missCleavageComboBox, 1, 1);
        
        // Keep protein names option
        keepProteinNamesCheckBox = new CheckBox("Keep protein names in library");
        keepProteinNamesCheckBox.setSelected(true);
        Label keepNamesHint = new Label("Uncheck if proteins are from DNA sequencing without annotation (will match to MAG database)");
        keepNamesHint.setStyle("-fx-font-size: 10px; -fx-text-fill: #666666;");
        keepNamesHint.setWrapText(true);
        
        enzymeBox.getChildren().addAll(enzymeTitle, enzymeGrid, keepProteinNamesCheckBox, keepNamesHint);
        
        // Initially hide enzyme parameters
        enzymeBox.setVisible(false);
        enzymeBox.setManaged(false);
        
        return enzymeBox;
    }
    
    private HBox createGenerationControls() {
        HBox generationBox = new HBox(10);
        generationBox.setAlignment(Pos.CENTER_LEFT);
        libraryNameField = new TextField();
        libraryNameField.setPromptText("Library name");
        libraryNameField.setPrefWidth(150);
        
        Button generateLibraryButton = new Button("Generate Library");
        generateLibraryButton.setPrefWidth(130);
        generateLibraryButton.setOnAction(e -> generateLibrary());
        
        generationBox.getChildren().addAll(
            new Label("Library Name:"), libraryNameField, generateLibraryButton);
        return generationBox;
    }
    
    private VBox createLibraryTable() {
        VBox tableSection = new VBox(10);
        
        // Header with controls
        HBox headerBox = new HBox(10);
        headerBox.setAlignment(Pos.CENTER_LEFT);
        Label availableLabel = new Label("Available Libraries (select for identification):");
        availableLabel.setStyle("-fx-font-weight: bold;");
        
        Button checkAvailableButton = new Button("Check Available");
        checkAvailableButton.setTooltip(new Tooltip("Auto-detect libraries in MAG database folder"));
        checkAvailableButton.setOnAction(e -> checkAvailableLibraries());
        
        Button deleteButton = new Button("Delete Selected");
        deleteButton.setOnAction(e -> deleteSelectedLibrary());
        
        headerBox.getChildren().addAll(availableLabel, checkAvailableButton, deleteButton);
        
        // Table
        libraryTable = new TableView<>(libraryData);
        libraryTable.setMinHeight(120);
        libraryTable.setPrefHeight(150);
        libraryTable.setMaxHeight(Double.MAX_VALUE);
        libraryTable.setEditable(true);
        libraryTable.setColumnResizePolicy(TableView.CONSTRAINED_RESIZE_POLICY);
        VBox.setVgrow(libraryTable, Priority.ALWAYS);

        // Use column (checkbox)
        TableColumn<LibraryRow, Boolean> useCol = new TableColumn<>("Use");
        useCol.setCellValueFactory(new PropertyValueFactory<>("use"));
        useCol.setCellFactory(CheckBoxTableCell.forTableColumn(useCol));
        useCol.setOnEditCommit(e -> notifyParameterUpdate());
        useCol.setMinWidth(40);
        useCol.setPrefWidth(50);

        // Library Name column
        TableColumn<LibraryRow, String> nameCol = new TableColumn<>("Library Name");
        nameCol.setCellValueFactory(new PropertyValueFactory<>("libName"));
        nameCol.setMinWidth(120);
        nameCol.setPrefWidth(150);

        // Peptide Count column
        TableColumn<LibraryRow, String> pepCountCol = new TableColumn<>("Peptide Count");
        pepCountCol.setCellValueFactory(new PropertyValueFactory<>("pepCount"));
        pepCountCol.setMinWidth(80);
        pepCountCol.setPrefWidth(100);
        
        // Protease column
        TableColumn<LibraryRow, String> proteaseCol = new TableColumn<>("Protease");
        proteaseCol.setCellValueFactory(new PropertyValueFactory<>("protease"));
        proteaseCol.setMinWidth(80);
        proteaseCol.setPrefWidth(100);
        
        // Miss Cleavages column  
        TableColumn<LibraryRow, String> missCol = new TableColumn<>("Miss Cleavages");
        missCol.setCellValueFactory(new PropertyValueFactory<>("miss"));
        missCol.setMinWidth(80);
        missCol.setPrefWidth(100);

        libraryTable.getColumns().addAll(useCol, nameCol, pepCountCol, proteaseCol, missCol);
        
        tableSection.getChildren().addAll(headerBox, libraryTable);
        return tableSection;
    }
    
    private VBox createProgressSection() {
        VBox progressSection = new VBox(8);
        progressSection.setPadding(new Insets(10));
        progressSection.setStyle("-fx-background-color: #f8f8f8; -fx-background-radius: 5;");
        progressSection.setMinHeight(120);
        progressSection.setMaxHeight(150);

//        Label progressTitle = new Label("Generation Progress:");
//        progressTitle.setStyle("-fx-font-weight: bold;");

//        progressLabel = new Label("Ready to generate library");
//        progressLabel.setStyle("-fx-font-size: 11px;");

        libraryProgressBar = new ProgressBar(0);
        libraryProgressBar.setMaxWidth(Double.MAX_VALUE);
        libraryProgressBar.setPrefHeight(20);

        Label consoleTitle = new Label("Console Output:");
        consoleTitle.setStyle("-fx-font-weight: bold; -fx-font-size: 11px;");

        // Create internal console if external one not provided
        if (consoleOutput == null) {
            consoleOutput = new TextArea();
        }
        consoleOutput.setEditable(false);
        consoleOutput.setPrefHeight(200);
        consoleOutput.setMaxHeight(300);
        consoleOutput.setStyle("-fx-font-family: monospace; -fx-font-size: 10px;");
        VBox.setVgrow(consoleOutput, Priority.ALWAYS);

        progressSection.getChildren().addAll(
            libraryProgressBar, consoleTitle, consoleOutput
        );
        
        return progressSection;
    }
    
    // Library Row class based on MagPepLib structure
    public static class LibraryRow {
        private final SimpleBooleanProperty use;
        private final SimpleStringProperty libName;
        private final SimpleStringProperty pepCount;
        private final SimpleStringProperty protease;
        private final SimpleStringProperty miss;
        private MagPepLib magPepLib; // Reference to original MagPepLib
        
        public LibraryRow(boolean use, MagPepLib magPepLib) {
            this.use = new SimpleBooleanProperty(use);
            this.libName = new SimpleStringProperty(magPepLib.getLibName());
            this.pepCount = new SimpleStringProperty(String.valueOf(magPepLib.getPepCount()));
            this.protease = new SimpleStringProperty(magPepLib.getProtease());
            this.miss = new SimpleStringProperty(String.valueOf(magPepLib.getMiss()));
            this.magPepLib = magPepLib;
        }
        
        public boolean getUse() { return use.get(); }
        public void setUse(boolean value) { use.set(value); }
        public BooleanProperty useProperty() { return use; }
        
        public String getLibName() { return libName.get(); }
        public StringProperty libNameProperty() { return libName; }
        
        public String getPepCount() { return pepCount.get(); }
        public StringProperty pepCountProperty() { return pepCount; }
        
        public String getProtease() { return protease.get(); }
        public StringProperty proteaseProperty() { return protease; }
        
        public String getMiss() { return miss.get(); }
        public StringProperty missProperty() { return miss; }
        
        public MagPepLib getMagPepLib() { return magPepLib; }
    }
    
    private void browseInputFile() {
        javafx.stage.FileChooser chooser = new javafx.stage.FileChooser();
        
        RadioButton selected = (RadioButton) inputTypeGroup.getSelectedToggle();
        if (selected.getText().contains("Library")) {
            chooser.setTitle("Select Library File");
            chooser.getExtensionFilters().add(new javafx.stage.FileChooser.ExtensionFilter("Library Files", "*.speclib", "*.sptxt", "*.tsv"));
        } else if (selected.getText().contains("Peptide")) {
            chooser.setTitle("Select Peptide List File");
            chooser.getExtensionFilters().add(new javafx.stage.FileChooser.ExtensionFilter("Peptide Files", "*.txt", "*.tsv", "*.csv"));
        } else if (selected.getText().contains("Protein")) {
            chooser.setTitle("Select Protein Database File");
            chooser.getExtensionFilters().add(new javafx.stage.FileChooser.ExtensionFilter("Protein Database Files", "*.fasta", "*.fa"));
        }
        
        File file = chooser.showOpenDialog(getScene().getWindow());
        if (file != null) {
            inputFileField.setText(file.getAbsolutePath());
        }
    }
    
    private void generateLibrary() {
        String inputFile = inputFileField.getText();
        String libraryName = libraryNameField.getText();
        
        if (inputFile == null || inputFile.trim().isEmpty()) {
            showError("Input Required", "Please select an input file.");
            return;
        }
        
        if (libraryName == null || libraryName.trim().isEmpty()) {
            showError("Name Required", "Please enter a library name.");
            return;
        }
        
        // Mock library generation process
        Task<Void> generateTask = new Task<Void>() {
            @Override
            protected Void call() throws Exception {
                Platform.runLater(() -> {
                    libraryProgressBar.setProgress(0);
                    progressLabel.setText("Starting library generation...");
                    if (consoleOutput != null) {
                        consoleOutput.appendText("Starting library generation: " + libraryName + "\n");
                    }
                });
                
                RadioButton selected = (RadioButton) inputTypeGroup.getSelectedToggle();
                String sourceType = selected.getText();
                
                // Step 1: Input validation and parsing
                Platform.runLater(() -> {
                    libraryProgressBar.setProgress(0.1);
                    progressLabel.setText("Validating input file...");
                    if (consoleOutput != null) {
                        consoleOutput.appendText("Validating input file: " + inputFile + "\n");
                        consoleOutput.appendText("Source type: " + sourceType + "\n");
                    }
                });
                Thread.sleep(1000);
                
                if (sourceType.contains("Protein")) {
                    Platform.runLater(() -> {
                        libraryProgressBar.setProgress(0.2);
                        progressLabel.setText("Configuring enzyme parameters...");
                        if (consoleOutput != null) {
                            consoleOutput.appendText("Using enzyme: " + enzymeComboBox.getValue() + "\n");
                            consoleOutput.appendText("Max miss cleavages: " + missCleavageComboBox.getValue() + "\n");
                            consoleOutput.appendText("Keep protein names: " + keepProteinNamesCheckBox.isSelected() + "\n");
                        }
                    });
                    Thread.sleep(500);
                }
                
                // Step 2: Processing input
                Platform.runLater(() -> {
                    libraryProgressBar.setProgress(0.4);
                    progressLabel.setText("Processing input data...");
                    if (consoleOutput != null) {
                        consoleOutput.appendText("Reading and processing input data...\n");
                    }
                });
                Thread.sleep(1500);
                
                // Step 3: Generating library
                Platform.runLater(() -> {
                    libraryProgressBar.setProgress(0.7);
                    progressLabel.setText("Generating spectral library...");
                    if (consoleOutput != null) {
                        consoleOutput.appendText("Building spectral library database...\n");
                    }
                });
                Thread.sleep(1000);
                
                // Step 4: Finalizing
                Platform.runLater(() -> {
                    libraryProgressBar.setProgress(0.9);
                    progressLabel.setText("Finalizing and saving...");
                    if (consoleOutput != null) {
                        consoleOutput.appendText("Writing library to disk...\n");
                    }
                });
                Thread.sleep(500);
                
                Platform.runLater(() -> {
                    libraryProgressBar.setProgress(1.0);
                    progressLabel.setText("Library generation completed successfully!");
                    if (consoleOutput != null) {
                        consoleOutput.appendText("Library generation completed: " + libraryName + "\n");
                        consoleOutput.appendText("Adding library to MAGDB configuration...\n");
                    }
                    
                    // Add to MAG database (mock implementation)
                    if (currentMagDbItem != null) {
                        // Mock values - in real implementation, these would be calculated
                        int pepCount = (int)(Math.random() * 5000) + 1000;
                        String protease = enzymeComboBox.getValue();
                        int miss = missCleavageComboBox.getValue();
                        
                        currentMagDbItem.addPepLib(libraryName, pepCount, protease, miss);
                        loadAvailableLibraries();
                        
                        if (consoleOutput != null) {
                            consoleOutput.appendText("Library added with " + pepCount + " peptides\n");
                        }
                    }
                    
                    notifyParameterUpdate();
                    
                    // Clear input fields and reset progress
                    inputFileField.clear();
                    libraryNameField.clear();
                    
                    // Reset progress after a delay
                    Timeline resetTimeline = new Timeline(new KeyFrame(Duration.seconds(3), e -> {
                        libraryProgressBar.setProgress(0);
                        progressLabel.setText("Ready to generate library");
                    }));
                    resetTimeline.play();
                });
                
                return null;
            }
        };
        
        new Thread(generateTask).start();
    }
    
    private void checkAvailableLibraries() {
        if (currentMagDbItem == null) {
            showError("No Database", "Please select a MAG database first.");
            return;
        }
        
        Task<Void> checkTask = new Task<Void>() {
            @Override
            protected Void call() throws Exception {
                Platform.runLater(() -> {
                    if (consoleOutput != null) {
                        consoleOutput.appendText("Checking for available libraries in MAG database folder...\n");
                    }
                });
                
                // Mock implementation - in real implementation, this would scan the libraries folder
                // and call currentMagDbItem.addPepLibFromFile(file) for each found library
                Thread.sleep(1000);
                
                Platform.runLater(() -> {
                    if (consoleOutput != null) {
                        consoleOutput.appendText("Library scan completed.\n");
                    }
                    loadAvailableLibraries();
                });
                
                return null;
            }
        };
        
        new Thread(checkTask).start();
    }
    
    private void deleteSelectedLibrary() {
        LibraryRow selected = libraryTable.getSelectionModel().getSelectedItem();
        if (selected == null) {
            showError("No Selection", "Please select a library to delete.");
            return;
        }
        
        // Create confirmation dialog with three options
        Alert alert = new Alert(Alert.AlertType.CONFIRMATION);
        alert.setTitle("Delete Library");
        alert.setHeaderText("Delete library: " + selected.getLibName());
        alert.setContentText("Choose delete option:");
        
        ButtonType tableOnlyButton = new ButtonType("Remove from table only");
        ButtonType deleteFilesButton = new ButtonType("Delete from table and disk files");
        ButtonType cancelButton = new ButtonType("Cancel", ButtonBar.ButtonData.CANCEL_CLOSE);
        
        alert.getButtonTypes().setAll(tableOnlyButton, deleteFilesButton, cancelButton);
        
        Optional<ButtonType> result = alert.showAndWait();
        if (result.isPresent()) {
            if (result.get() == tableOnlyButton) {
                // Remove from table and MAG database only
                if (currentMagDbItem != null) {
                    currentMagDbItem.removePepLib(selected.getMagPepLib());
                }
                libraryData.remove(selected);
                notifyParameterUpdate();
                
                if (consoleOutput != null) {
                    consoleOutput.appendText("Library removed from table: " + selected.getLibName() + "\n");
                }
            } else if (result.get() == deleteFilesButton) {
                // Remove from table, MAG database, and delete files
                if (currentMagDbItem != null) {
                    currentMagDbItem.removePepLib(selected.getMagPepLib());
                }
                // Delete files from disk
                selected.getMagPepLib().delete();
                libraryData.remove(selected);
                notifyParameterUpdate();
                
                if (consoleOutput != null) {
                    consoleOutput.appendText("Library and files deleted: " + selected.getLibName() + "\n");
                }
            }
            // Cancel - do nothing
        }
    }
    
    private void loadAvailableLibraries() {
        libraryData.clear();
        if (currentMagDbItem != null) {
            MagPepLib[] magPepLibs = currentMagDbItem.getMagPepLibs();
            for (MagPepLib pepLib : magPepLibs) {
                libraryData.add(new LibraryRow(false, pepLib));
            }
        }
    }
    
    private Label createSectionTitle(String title) {
        Label label = new Label(title);
        label.setFont(Font.font("System", FontWeight.BOLD, 14));
        label.setTextFill(Color.DARKBLUE);
        return label;
    }
    
    private void showError(String title, String content) {
        Alert alert = new Alert(Alert.AlertType.ERROR);
        alert.setTitle(title);
        alert.setHeaderText(null);
        alert.setContentText(content);
        alert.showAndWait();
    }
    
    private void notifyParameterUpdate() {
        if (parameterCallback != null) {
            parameterCallback.onParametersChanged();
        }
    }
    
    // Public methods for accessing selected libraries
    public MagPepLib[] getSelectedLibraries() {
        return libraryData.stream()
            .filter(LibraryRow::getUse)
            .map(LibraryRow::getMagPepLib)
            .toArray(MagPepLib[]::new);
    }
    
    public String[] getSelectedLibraryNames() {
        return libraryData.stream()
            .filter(LibraryRow::getUse)
            .map(LibraryRow::getLibName)
            .toArray(String[]::new);
    }
} 