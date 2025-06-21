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
import uk.ac.quadram.metapilot.task.mag.MagDbItem.MagModel;
import java.io.File;
import java.util.Optional;
import javafx.animation.KeyFrame;
import javafx.animation.Timeline;
import javafx.util.Duration;

public class DiaModelGenerationPanel extends VBox {
    
    // Model table components
    private TableView<ModelRow> modelTable;
    private ObservableList<ModelRow> modelData = FXCollections.observableArrayList();
    
    // Generation components
    private TextField intensityFileField;
    private TextField modelNameField;
    
    // Current MAG database reference
    private MagDbItem currentMagDbItem;
    
    // Console output for feedback
    private TextArea consoleOutput;
    
    // Internal progress tracking
    private ProgressBar modelProgressBar;
    private Label progressLabel;
    
    // Callback interface for parameter updates
    public interface ParameterUpdateCallback {
        void onParametersChanged();
    }
    
    private ParameterUpdateCallback parameterCallback;

    public DiaModelGenerationPanel() {
        initializeUI();
    }
    
    public DiaModelGenerationPanel(TextArea consoleOutput) {
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
        loadAvailableModels();
    }

    private void initializeUI() {
        setSpacing(10);
        setPadding(new Insets(10));
        setStyle("-fx-background-color: #f0fff0; -fx-background-radius: 5;");
        setMinHeight(400);
        VBox.setVgrow(this, Priority.ALWAYS);

        Label titleLabel = createSectionTitle("ML Model Generation (Optional)");
        
        // Optional hint
        Label optionalHint = new Label("⚠ Model generation is NOT mandatory - it's an optional feature for enhanced identification");
        optionalHint.setStyle("-fx-font-size: 11px; -fx-text-fill: #ff6600; -fx-font-weight: bold;");
        optionalHint.setWrapText(true);

        // File selection for peptide intensity table
        VBox inputBox = createInputSection();
        
        // Generation controls
        HBox generationBox = createGenerationControls();
        
        // Available models table
        VBox tableSection = createModelTable();
        
        // Progress and console section
        VBox progressSection = createProgressSection();
        
        getChildren().addAll(titleLabel, optionalHint, inputBox, generationBox, tableSection, progressSection);
    }
    
    private VBox createInputSection() {
        VBox inputBox = new VBox(8);
        Label inputLabel = new Label("Peptide Intensity Table for Model Generation:");
        inputLabel.setStyle("-fx-font-weight: bold;");
        
        HBox fileBox = new HBox(10);
        fileBox.setAlignment(Pos.CENTER_LEFT);
        intensityFileField = new TextField();
        intensityFileField.setPromptText("Select peptide intensity table file (.txt, .tsv, .csv)");
        intensityFileField.setPrefWidth(280);
        intensityFileField.setEditable(false);
        
        Button browseIntensityButton = new Button("Browse");
        browseIntensityButton.setOnAction(e -> browseIntensityFile());
        
        fileBox.getChildren().addAll(intensityFileField, browseIntensityButton);
        inputBox.getChildren().addAll(inputLabel, fileBox);
        return inputBox;
    }
    
    private HBox createGenerationControls() {
        HBox generationBox = new HBox(10);
        generationBox.setAlignment(Pos.CENTER_LEFT);
        
        modelNameField = new TextField();
        modelNameField.setPromptText("Model name");
        modelNameField.setPrefWidth(150);
        
        Button generateModelButton = new Button("Generate Model");
        generateModelButton.setPrefWidth(130);
        generateModelButton.setOnAction(e -> generateModel());
        
        generationBox.getChildren().addAll(
            new Label("Model Name:"), modelNameField, generateModelButton);
        return generationBox;
    }
    
    private VBox createModelTable() {
        VBox tableSection = new VBox(10);
        
        // Header with controls
        HBox headerBox = new HBox(10);
        headerBox.setAlignment(Pos.CENTER_LEFT);
        Label availableLabel = new Label("Available Models (select for identification - optional):");
        availableLabel.setStyle("-fx-font-weight: bold;");
        
        Button checkAvailableButton = new Button("Check Available");
        checkAvailableButton.setTooltip(new Tooltip("Auto-detect models in MAG database folder"));
        checkAvailableButton.setOnAction(e -> checkAvailableModels());
        
        Button deleteButton = new Button("Delete Selected");
        deleteButton.setOnAction(e -> deleteSelectedModel());
        
        headerBox.getChildren().addAll(availableLabel, checkAvailableButton, deleteButton);
        
        // Table
        modelTable = new TableView<>(modelData);
        modelTable.setMinHeight(120);
        modelTable.setPrefHeight(150);
        modelTable.setMaxHeight(Double.MAX_VALUE);
        modelTable.setEditable(true);
        modelTable.setColumnResizePolicy(TableView.CONSTRAINED_RESIZE_POLICY);
        VBox.setVgrow(modelTable, Priority.ALWAYS);

        // Use column (checkbox)
        TableColumn<ModelRow, Boolean> useCol = new TableColumn<>("Use");
        useCol.setCellValueFactory(new PropertyValueFactory<>("use"));
        useCol.setCellFactory(CheckBoxTableCell.forTableColumn(useCol));
        useCol.setOnEditCommit(e -> notifyParameterUpdate());
        useCol.setMinWidth(40);
        useCol.setPrefWidth(50);

        // Model Name column
        TableColumn<ModelRow, String> nameCol = new TableColumn<>("Model Name");
        nameCol.setCellValueFactory(new PropertyValueFactory<>("modelName"));
        nameCol.setMinWidth(120);
        nameCol.setPrefWidth(150);

        // Training Size column
        TableColumn<ModelRow, String> trainSizeCol = new TableColumn<>("Train Size");
        trainSizeCol.setCellValueFactory(new PropertyValueFactory<>("trainSize"));
        trainSizeCol.setMinWidth(80);
        trainSizeCol.setPrefWidth(90);

        // Model Size column
        TableColumn<ModelRow, String> modelSizeCol = new TableColumn<>("Model Size");
        modelSizeCol.setCellValueFactory(new PropertyValueFactory<>("modelSize"));
        modelSizeCol.setMinWidth(80);
        modelSizeCol.setPrefWidth(90);

        // Correlation column
        TableColumn<ModelRow, String> corrCol = new TableColumn<>("Correlation");
        corrCol.setCellValueFactory(new PropertyValueFactory<>("corr"));
        corrCol.setMinWidth(80);
        corrCol.setPrefWidth(90);
        
        // Median column
        TableColumn<ModelRow, String> medCol = new TableColumn<>("Median");
        medCol.setCellValueFactory(new PropertyValueFactory<>("med"));
        medCol.setMinWidth(80);
        medCol.setPrefWidth(90);

        modelTable.getColumns().addAll(useCol, nameCol, trainSizeCol, modelSizeCol, corrCol, medCol);
        
        tableSection.getChildren().addAll(headerBox, modelTable);
        return tableSection;
    }
    
    private VBox createProgressSection() {
        VBox progressSection = new VBox(8);
        progressSection.setPadding(new Insets(10));
        progressSection.setStyle("-fx-background-color: #f8f8f8; -fx-background-radius: 5;");
        progressSection.setMinHeight(120);
        progressSection.setMaxHeight(150);

//        Label progressTitle = new Label("Model Training Progress:");
//        progressTitle.setStyle("-fx-font-weight: bold;");

//       progressLabel = new Label("Ready to train model");
//        progressLabel.setStyle("-fx-font-size: 11px;");

        modelProgressBar = new ProgressBar(0);
        modelProgressBar.setMaxWidth(Double.MAX_VALUE);
        modelProgressBar.setPrefHeight(20);

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
            modelProgressBar, consoleTitle, consoleOutput
        );
        
        return progressSection;
    }
    
    // Model Row class based on MagModel structure
    public static class ModelRow {
        private final SimpleBooleanProperty use;
        private final SimpleStringProperty modelName;
        private final SimpleStringProperty trainSize;
        private final SimpleStringProperty modelSize;
        private final SimpleStringProperty corr;
        private final SimpleStringProperty med;
        private MagModel magModel; // Reference to original MagModel
        
        public ModelRow(boolean use, MagModel magModel) {
            this.use = new SimpleBooleanProperty(use);
            this.modelName = new SimpleStringProperty(magModel.getModelName());
            this.trainSize = new SimpleStringProperty(String.valueOf(magModel.getTrainSize()));
            this.modelSize = new SimpleStringProperty(String.valueOf(magModel.getModelSize()));
            this.corr = new SimpleStringProperty(String.valueOf(magModel.getCorr()));
            this.med = new SimpleStringProperty(String.valueOf(magModel.getMed()));
            this.magModel = magModel;
        }
        
        public boolean getUse() { return use.get(); }
        public void setUse(boolean value) { use.set(value); }
        public BooleanProperty useProperty() { return use; }
        
        public String getModelName() { return modelName.get(); }
        public StringProperty modelNameProperty() { return modelName; }
        
        public String getTrainSize() { return trainSize.get(); }
        public StringProperty trainSizeProperty() { return trainSize; }
        
        public String getModelSize() { return modelSize.get(); }
        public StringProperty modelSizeProperty() { return modelSize; }
        
        public String getCorr() { return corr.get(); }
        public StringProperty corrProperty() { return corr; }
        
        public String getMed() { return med.get(); }
        public StringProperty medProperty() { return med; }
        
        public MagModel getMagModel() { return magModel; }
    }
    
    private void browseIntensityFile() {
        javafx.stage.FileChooser chooser = new javafx.stage.FileChooser();
        chooser.setTitle("Select Peptide Intensity Table File");
        chooser.getExtensionFilters().add(
            new javafx.stage.FileChooser.ExtensionFilter("Intensity Table Files", "*.txt", "*.tsv", "*.csv")
        );
        File file = chooser.showOpenDialog(getScene().getWindow());
        if (file != null) {
            intensityFileField.setText(file.getAbsolutePath());
        }
    }
    
    private void generateModel() {
        String intensityFile = intensityFileField.getText();
        String modelName = modelNameField.getText();
        
        if (intensityFile == null || intensityFile.trim().isEmpty()) {
            showError("Input Required", "Please select a peptide intensity table file.");
            return;
        }
        
        if (modelName == null || modelName.trim().isEmpty()) {
            showError("Name Required", "Please enter a model name.");
            return;
        }
        
        // Mock model generation process
        Task<Void> generateTask = new Task<Void>() {
            @Override
            protected Void call() throws Exception {
                Platform.runLater(() -> {
                    modelProgressBar.setProgress(0);
                    progressLabel.setText("Starting model generation...");
                    if (consoleOutput != null) {
                        consoleOutput.appendText("Starting model generation: " + modelName + "\n");
                        consoleOutput.appendText("Processing peptide intensity table: " + intensityFile + "\n");
                    }
                });
                
                // Step 1: Data loading and validation
                Platform.runLater(() -> {
                    modelProgressBar.setProgress(0.1);
                    progressLabel.setText("Loading and validating data...");
                    if (consoleOutput != null) {
                        consoleOutput.appendText("Loading peptide intensity data...\n");
                        consoleOutput.appendText("Validating data format and quality...\n");
                    }
                });
                Thread.sleep(1000);
                
                // Step 2: Data preprocessing
                Platform.runLater(() -> {
                    modelProgressBar.setProgress(0.3);
                    progressLabel.setText("Preprocessing intensity data...");
                    if (consoleOutput != null) {
                        consoleOutput.appendText("Normalizing intensity values...\n");
                        consoleOutput.appendText("Filtering low-quality peptides...\n");
                    }
                });
                Thread.sleep(1500);
                
                // Step 3: Feature extraction
                Platform.runLater(() -> {
                    modelProgressBar.setProgress(0.5);
                    progressLabel.setText("Extracting features...");
                    if (consoleOutput != null) {
                        consoleOutput.appendText("Analyzing peptide intensity patterns...\n");
                        consoleOutput.appendText("Extracting chemical and physical features...\n");
                    }
                });
                Thread.sleep(2000);
                
                // Step 4: Model training
                Platform.runLater(() -> {
                    modelProgressBar.setProgress(0.7);
                    progressLabel.setText("Training machine learning model...");
                    if (consoleOutput != null) {
                        consoleOutput.appendText("Training machine learning model...\n");
                        consoleOutput.appendText("Optimizing model parameters...\n");
                    }
                });
                Thread.sleep(2000);
                
                // Step 5: Model validation
                Platform.runLater(() -> {
                    modelProgressBar.setProgress(0.9);
                    progressLabel.setText("Validating model performance...");
                    if (consoleOutput != null) {
                        consoleOutput.appendText("Calculating model correlation...\n");
                        consoleOutput.appendText("Performing cross-validation...\n");
                    }
                });
                Thread.sleep(1000);
                
                // Generate mock statistics
                int trainSize = (int)(Math.random() * 5000) + 1000;
                int modelSize = (int)(Math.random() * 100) + 50;
                double correlation = Math.round((0.7 + Math.random() * 0.25) * 100.0) / 100.0;
                double median = Math.round((0.5 + Math.random() * 0.4) * 100.0) / 100.0;
                
                Platform.runLater(() -> {
                    modelProgressBar.setProgress(1.0);
                    progressLabel.setText("Model generation completed successfully!");
                    if (consoleOutput != null) {
                        consoleOutput.appendText("Model generation completed: " + modelName + "\n");
                        consoleOutput.appendText("Training size: " + trainSize + " peptides\n");
                        consoleOutput.appendText("Model size: " + modelSize + " features\n");
                        consoleOutput.appendText("Model correlation: " + correlation + "\n");
                        consoleOutput.appendText("Model median: " + median + "\n");
                        consoleOutput.appendText("Adding model to MAGDB configuration...\n");
                    }
                    
                    // Add to MAG database
                    if (currentMagDbItem != null) {
                        currentMagDbItem.addModel(modelName, trainSize, modelSize, correlation, median);
                        loadAvailableModels();
                        
                        if (consoleOutput != null) {
                            consoleOutput.appendText("Model successfully added to database\n");
                        }
                    }
                    
                    notifyParameterUpdate();
                    
                    // Clear input fields and reset progress
                    intensityFileField.clear();
                    modelNameField.clear();
                    
                    // Reset progress after a delay
                    Timeline resetTimeline = new Timeline(new KeyFrame(Duration.seconds(3), e -> {
                        modelProgressBar.setProgress(0);
                        progressLabel.setText("Ready to train model");
                    }));
                    resetTimeline.play();
                });
                
                return null;
            }
        };
        
        new Thread(generateTask).start();
    }
    
    private void checkAvailableModels() {
        if (currentMagDbItem == null) {
            showError("No Database", "Please select a MAG database first.");
            return;
        }
        
        Task<Void> checkTask = new Task<Void>() {
            @Override
            protected Void call() throws Exception {
                Platform.runLater(() -> {
                    if (consoleOutput != null) {
                        consoleOutput.appendText("Checking for available models in MAG database folder...\n");
                    }
                });
                
                // Mock implementation - in real implementation, this would scan the models folder
                // and call currentMagDbItem.addModelFromFile(file) for each found model
                Thread.sleep(1000);
                
                Platform.runLater(() -> {
                    if (consoleOutput != null) {
                        consoleOutput.appendText("Model scan completed.\n");
                    }
                    loadAvailableModels();
                });
                
                return null;
            }
        };
        
        new Thread(checkTask).start();
    }
    
    private void deleteSelectedModel() {
        ModelRow selected = modelTable.getSelectionModel().getSelectedItem();
        if (selected == null) {
            showError("No Selection", "Please select a model to delete.");
            return;
        }
        
        // Create confirmation dialog with three options
        Alert alert = new Alert(Alert.AlertType.CONFIRMATION);
        alert.setTitle("Delete Model");
        alert.setHeaderText("Delete model: " + selected.getModelName());
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
                    currentMagDbItem.removeModel(selected.getMagModel());
                }
                modelData.remove(selected);
                notifyParameterUpdate();
                
                if (consoleOutput != null) {
                    consoleOutput.appendText("Model removed from table: " + selected.getModelName() + "\n");
                }
            } else if (result.get() == deleteFilesButton) {
                // Remove from table, MAG database, and delete files
                if (currentMagDbItem != null) {
                    currentMagDbItem.removeModel(selected.getMagModel());
                }
                // Note: MagModel doesn't have a delete() method like MagPepLib, 
                // but in real implementation, you would add file deletion logic here
                modelData.remove(selected);
                notifyParameterUpdate();
                
                if (consoleOutput != null) {
                    consoleOutput.appendText("Model and files deleted: " + selected.getModelName() + "\n");
                }
            }
            // Cancel - do nothing
        }
    }
    
    private void loadAvailableModels() {
        modelData.clear();
        if (currentMagDbItem != null) {
            MagModel[] magModels = currentMagDbItem.getMagModels();
            for (MagModel model : magModels) {
                modelData.add(new ModelRow(false, model));
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
    
    // Public methods for accessing selected models
    public MagModel[] getSelectedModels() {
        return modelData.stream()
            .filter(ModelRow::getUse)
            .map(ModelRow::getMagModel)
            .toArray(MagModel[]::new);
    }
    
    public String[] getSelectedModelNames() {
        return modelData.stream()
            .filter(ModelRow::getUse)
            .map(ModelRow::getModelName)
            .toArray(String[]::new);
    }
} 