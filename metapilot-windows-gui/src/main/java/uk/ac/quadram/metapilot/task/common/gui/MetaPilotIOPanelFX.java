package uk.ac.quadram.metapilot.task.common.gui;

import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.control.*;
import javafx.scene.control.cell.TextFieldTableCell;
import javafx.scene.layout.HBox;
import javafx.scene.layout.VBox;
import javafx.stage.DirectoryChooser;
import javafx.stage.FileChooser;
import javafx.util.converter.DefaultStringConverter;
import javafx.beans.property.SimpleStringProperty;
import javafx.beans.property.StringProperty;
import javafx.beans.property.SimpleBooleanProperty;
import javafx.beans.property.BooleanProperty;

import java.io.File;
import java.util.List;
import java.util.ArrayList;
import java.util.HashSet;

public class MetaPilotIOPanelFX extends VBox {
    private final TableView<ExperimentRow> experimentTable;
    private final ObservableList<ExperimentRow> experimentData;
    private final TextField outputDirField;
    private final Button browseRawFilesButton;
    private final Button removeExperimentButton;
    private final Button removeAllExperimentsButton;
    private final Button browseOutputDirButton;
    
    // Experiment naming components
    private ComboBox<String> separatorComboBox;
    private TextField customSeparatorField;
    private TableView<PatternRow> patternTable;
    private ObservableList<PatternRow> patternData;
    private Button generateNamesButton;
    private Button resetNamesButton;
    
    // Pattern management
    private String currentPattern = "";
    private boolean patternTableUpdate = false;
    private HashSet<String> fileSet;

    @SuppressWarnings("unchecked")
    public MetaPilotIOPanelFX() {
        setSpacing(10);
        setPadding(new Insets(10));
        setStyle("-fx-border-color: #cccccc; -fx-border-radius: 5; -fx-border-width: 1;");

        fileSet = new HashSet<>();

        // Initialize pattern data FIRST
        patternData = FXCollections.observableArrayList();

        // Raw Files Section
        Label rawFilesLabel = new Label("Raw Files:");
        rawFilesLabel.setStyle("-fx-font-weight: bold; -fx-font-size: 14px;");
        
        browseRawFilesButton = new Button("Select Raw Files");
        browseRawFilesButton.setOnAction(e -> selectRawFiles());

        experimentData = FXCollections.observableArrayList();
        experimentTable = new TableView<>(experimentData);
        experimentTable.setEditable(true);
        experimentTable.setPrefHeight(250);
        
        TableColumn<ExperimentRow, String> fileCol = new TableColumn<>("File Path");
        fileCol.setCellValueFactory(cellData -> cellData.getValue().filePathProperty());
        fileCol.setPrefWidth(400);
        fileCol.setCellFactory(TextFieldTableCell.forTableColumn(new DefaultStringConverter()));
        fileCol.setOnEditCommit(event -> {
            event.getRowValue().setFilePath(event.getNewValue());
        });
        
        TableColumn<ExperimentRow, String> expCol = new TableColumn<>("Experiment Name");
        expCol.setCellValueFactory(cellData -> cellData.getValue().experimentNameProperty());
        expCol.setCellFactory(column -> {
            TextFieldTableCell<ExperimentRow, String> cell = new TextFieldTableCell<ExperimentRow, String>(new DefaultStringConverter()) {
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
        expCol.setOnEditCommit(event -> {
            ExperimentRow row = event.getRowValue();
            row.setExperimentName(event.getNewValue());
            
            // Update pattern table if this is the first row
            if (!experimentData.isEmpty() && experimentData.get(0) == event.getRowValue()) {
                updatePatternTable();
            }
        });
        expCol.setPrefWidth(200);
        
        TableColumn<ExperimentRow, String> fractionCol = new TableColumn<>("Fraction");
        fractionCol.setCellValueFactory(cellData -> cellData.getValue().fractionProperty());
        fractionCol.setCellFactory(column -> {
            TextFieldTableCell<ExperimentRow, String> cell = new TextFieldTableCell<ExperimentRow, String>(new DefaultStringConverter()) {
                @Override
                public void startEdit() {
                    super.startEdit();
                    TextField textField = (TextField) getGraphic();
                    if (textField != null) {
                        // Only allow integer input
                        textField.textProperty().addListener((obs, oldVal, newVal) -> {
                            if (newVal != null && !newVal.matches("\\d*")) {
                                textField.setText(oldVal);
                            }
                        });
                        
                        // Commit edit when focus is lost (like Swing)
                        textField.focusedProperty().addListener((obs, wasFocused, isNowFocused) -> {
                            if (wasFocused && !isNowFocused && isEditing()) {
                                String text = textField.getText();
                                // Ensure it's a valid integer, default to 1 if empty or invalid
                                if (text == null || text.trim().isEmpty()) {
                                    text = "1";
                                } else {
                                    try {
                                        int value = Integer.parseInt(text);
                                        if (value <= 0) {
                                            text = "1"; // Ensure positive integer
                                        }
                                    } catch (NumberFormatException e) {
                                        text = "1"; // Default to 1 if not a valid integer
                                    }
                                }
                                textField.setText(text);
                                commitEdit(text);
                            }
                        });
                    }
                }
            };
            return cell;
        });
        fractionCol.setOnEditCommit(event -> {
            ExperimentRow row = event.getRowValue();
            String newValue = event.getNewValue();
            
            // Validate and set the fraction value
            if (newValue == null || newValue.trim().isEmpty()) {
                newValue = "1";
            } else {
                try {
                    int value = Integer.parseInt(newValue);
                    if (value <= 0) {
                        newValue = "1"; // Ensure positive integer
                    }
                } catch (NumberFormatException e) {
                    newValue = "1"; // Default to 1 if not a valid integer
                }
            }
            row.setFraction(newValue);
        });
        fractionCol.setPrefWidth(80);
        
        experimentTable.getColumns().setAll(fileCol, expCol, fractionCol);

        HBox expButtonBox = new HBox(10);
        removeExperimentButton = new Button("Remove Selected");
        removeExperimentButton.setOnAction(e -> removeSelectedExperiment());
        removeExperimentButton.setTooltip(new Tooltip("Remove the selected experiment row"));
        removeAllExperimentsButton = new Button("Remove All");
        removeAllExperimentsButton.setOnAction(e -> removeAllExperiments());
        removeAllExperimentsButton.setTooltip(new Tooltip("Remove all raw files"));
        expButtonBox.getChildren().addAll(removeExperimentButton, removeAllExperimentsButton);

        // Create left panel for Raw Files
        VBox leftPanel = new VBox(10);
        leftPanel.setPadding(new Insets(10));
        leftPanel.setStyle("-fx-background-color: #f8f8f8; -fx-border-color: #cccccc; -fx-border-radius: 5; -fx-border-width: 1;");
        leftPanel.getChildren().addAll(
            rawFilesLabel,
            browseRawFilesButton,
            experimentTable,
            expButtonBox
        );
        
        // Experiment Name Generation Section
        VBox nameGenerationSection = createNameGenerationSection();

        // Create horizontal layout for the two main panels
        HBox mainPanelsBox = new HBox(20);
        mainPanelsBox.getChildren().addAll(leftPanel, nameGenerationSection);
        
        // Make both panels grow equally to fill available space
        HBox.setHgrow(leftPanel, javafx.scene.layout.Priority.ALWAYS);
        HBox.setHgrow(nameGenerationSection, javafx.scene.layout.Priority.ALWAYS);

        // Output Directory Section
        Label outputDirLabel = new Label("Output Directory:");
        outputDirLabel.setStyle("-fx-font-weight: bold; -fx-font-size: 14px;");
        
        outputDirField = new TextField();
        outputDirField.setPrefWidth(400);
        browseOutputDirButton = new Button("Browse");
        browseOutputDirButton.setOnAction(e -> selectOutputDirectory());
        HBox outputDirBox = new HBox(10, outputDirField, browseOutputDirButton);

        getChildren().addAll(
            mainPanelsBox,
            outputDirLabel,
            outputDirBox
        );
        
        updatePatternTable(); // Initial call
    }
    
    @SuppressWarnings("unchecked")
    private VBox createNameGenerationSection() {
        VBox section = new VBox(10);
        section.setPadding(new Insets(10));
        section.setStyle("-fx-background-color: #f8f8f8; -fx-border-color: #cccccc; -fx-border-radius: 5; -fx-border-width: 1;");
        
        Label titleLabel = new Label("Experiment Name Generation");
        titleLabel.setStyle("-fx-font-weight: bold; -fx-font-size: 14px;");
        
        // Separator selection
        HBox separatorBox = new HBox(10);
        separatorBox.setAlignment(Pos.CENTER_LEFT);
        
        Label separatorLabel = new Label("Split Symbol:");
        separatorLabel.setStyle("-fx-font-weight: bold;");
        
        separatorComboBox = new ComboBox<>(FXCollections.observableArrayList("_", "-", "_ and -", "Other"));
        separatorComboBox.setValue("_");
        separatorComboBox.setPrefWidth(120);
        separatorComboBox.setOnAction(e -> {
            String selected = separatorComboBox.getValue();
            customSeparatorField.setDisable(!"Other".equals(selected));
            updatePatternTable();
        });
        
        customSeparatorField = new TextField();
        customSeparatorField.setPromptText("Enter custom separator");
        customSeparatorField.setPrefWidth(150);
        customSeparatorField.setDisable(true);
        customSeparatorField.textProperty().addListener((obs, oldVal, newVal) -> {
            if (!customSeparatorField.isDisabled()) {
                updatePatternTable();
            }
        });
        
        separatorBox.getChildren().addAll(separatorLabel, separatorComboBox, customSeparatorField);
        
        // Pattern table
        Label patternLabel = new Label("Select parts to include in experiment name:");
        patternLabel.setStyle("-fx-font-weight: bold;");
        
        patternTable = new TableView<>(patternData);
        patternTable.setEditable(true);
        patternTable.setPrefHeight(180);
        
        TableColumn<PatternRow, Boolean> includeCol = new TableColumn<>("Keep in experiment name");
        includeCol.setCellValueFactory(cellData -> cellData.getValue().includeProperty());
        includeCol.setCellFactory(column -> {
            CheckBox checkBox = new CheckBox();
            TableCell<PatternRow, Boolean> cell = new TableCell<PatternRow, Boolean>() {
                @Override
                protected void updateItem(Boolean item, boolean empty) {
                    super.updateItem(item, empty);
                    if (empty || getTableRow() == null || getTableRow().getItem() == null) {
                        setGraphic(null);
                    } else {
                        PatternRow row = getTableRow().getItem();
                        checkBox.setSelected(row.isInclude());
                        checkBox.setOnAction(e -> row.setInclude(checkBox.isSelected()));
                        setGraphic(checkBox);
                    }
                }
            };
            return cell;
        });
        includeCol.setPrefWidth(200);
        
        TableColumn<PatternRow, String> contentCol = new TableColumn<>("Content");
        contentCol.setCellValueFactory(cellData -> cellData.getValue().contentProperty());
        contentCol.setPrefWidth(200);
        
        patternTable.getColumns().addAll(includeCol, contentCol);
        
        // Control buttons
        HBox buttonBox = new HBox(10);
        buttonBox.setAlignment(Pos.CENTER_LEFT);
        
        generateNamesButton = new Button("Generate Names");
        generateNamesButton.setOnAction(e -> generateExp());
        
        resetNamesButton = new Button("Reset");
        resetNamesButton.setOnAction(e -> restoreExpNames());
        
        buttonBox.getChildren().addAll(generateNamesButton, resetNamesButton);
        
        section.getChildren().addAll(titleLabel, separatorBox, patternLabel, patternTable, buttonBox);
        return section;
    }
    
    private void updatePatternTable() {
        System.out.println("=== updatePatternTable called ===");
        System.out.println("patternTableUpdate flag: " + patternTableUpdate);
        System.out.println("experimentData.isEmpty(): " + experimentData.isEmpty());
        
        if (patternTableUpdate) {
            System.out.println("Returning early due to patternTableUpdate flag");
            return;
        }
        
        synchronized (this) {
            patternTableUpdate = true;
            
            patternData.clear();
            System.out.println("Cleared patternData");
            
            if (experimentData.isEmpty()) {
                System.out.println("experimentData is empty, returning");
                patternTableUpdate = false;
                return;
            }
            
            // Get experiment name from first row - use filename if experiment name is empty
            String firstExpName = experimentData.get(0).getExperimentName();
            System.out.println("First experiment name: '" + firstExpName + "'");
            
            if (firstExpName == null || firstExpName.trim().isEmpty()) {
                // If experiment name is empty, use filename from file path
                String filePath = experimentData.get(0).getFilePath();
                System.out.println("Experiment name empty, using file path: '" + filePath + "'");
                
                if (filePath != null && !filePath.isEmpty()) {
                    String fileName = new File(filePath).getName();
                    System.out.println("File name: '" + fileName + "'");
                    
                    if (fileName.contains(".")) {
                        firstExpName = fileName.substring(0, fileName.lastIndexOf("."));
                    } else {
                        firstExpName = fileName;
                    }
                    System.out.println("Final firstExpName: '" + firstExpName + "'");
                } else {
                    System.out.println("File path is null or empty, returning");
                    patternTableUpdate = false;
                    return;
                }
            }
            
            // Determine pattern
            String symbol = separatorComboBox.getValue();
            System.out.println("Selected separator: '" + symbol + "'");
            
            if ("Other".equals(symbol)) {
                currentPattern = customSeparatorField.getText();
                System.out.println("Using custom pattern: '" + currentPattern + "'");
                if (currentPattern == null || currentPattern.trim().isEmpty()) {
                    System.out.println("Custom pattern is empty, returning");
                    patternTableUpdate = false;
                    return;
                }
            } else if ("_ and -".equals(symbol)) {
                currentPattern = "[_-]";
                System.out.println("Using _ and - pattern: '" + currentPattern + "'");
            } else {
                currentPattern = symbol;
                System.out.println("Using symbol pattern: '" + currentPattern + "'");
            }
            
            if (currentPattern != null && currentPattern.length() > 0) {
                String[] parts = firstExpName.split(currentPattern);
                System.out.println("Split '" + firstExpName + "' with pattern '" + currentPattern + "' into " + parts.length + " parts:");
                
                for (int i = 0; i < parts.length; i++) {
                    System.out.println("  Part " + i + ": '" + parts[i] + "'");
                    PatternRow row = new PatternRow(false, parts[i]);
                    patternData.add(row);
                }
                System.out.println("Added " + patternData.size() + " rows to patternData");
            } else {
                System.out.println("currentPattern is null or empty");
            }
            
            patternTableUpdate = false;
            System.out.println("=== updatePatternTable finished ===");
        }
    }
    
    private void generateExp() {
        if (patternTableUpdate || experimentData.isEmpty()) {
            return;
        }
        
        // Get pattern selection - which parts are checked
        List<Integer> selectedIndices = new ArrayList<>();
        for (int i = 0; i < patternData.size(); i++) {
            if (patternData.get(i).isInclude()) {
                selectedIndices.add(i);
            }
        }
        
        if (selectedIndices.isEmpty()) {
            return; // No parts selected, don't change anything
        }
        
        // Apply to all rows
        for (ExperimentRow row : experimentData) {
            String filePath = row.getFilePath();
            if (filePath != null && !filePath.isEmpty()) {
                // Get filename without extension
                String fileName = new File(filePath).getName();
                if (fileName.contains(".")) {
                    fileName = fileName.substring(0, fileName.lastIndexOf("."));
                }
                
                // Split the filename using current pattern
                String[] parts = fileName.split(currentPattern);
                
                // Build new experiment name from selected parts
                StringBuilder newExpName = new StringBuilder();
                boolean first = true;
                for (int index : selectedIndices) {
                    if (index < parts.length) {
                        if (!first) {
                            newExpName.append("_");
                        }
                        newExpName.append(parts[index]);
                        first = false;
                    }
                }
                
                if (newExpName.length() > 0) {
                    row.setExperimentName(newExpName.toString());
                }
            }
        }
        
        experimentTable.refresh();
    }
    
    private void restoreExpNames() {
        String[] originalExpNames = getExpFromRawName();
        for (int i = 0; i < originalExpNames.length && i < experimentData.size(); i++) {
            experimentData.get(i).setExperimentName(originalExpNames[i]);
        }
        experimentTable.refresh();
        updatePatternTable();
    }
    
    private String[] getExpFromRawName() {
        List<String> expNames = new ArrayList<>();
        for (ExperimentRow row : experimentData) {
            String path = row.getFilePath();
            if (path != null && !path.isEmpty()) {
                String name = new File(path).getName();
                if (name.contains(".")) {
                    name = name.substring(0, name.lastIndexOf("."));
                }
                expNames.add(name);
            } else {
                expNames.add("");
            }
        }
        return expNames.toArray(new String[0]);
    }

    private void selectRawFiles() {
        // Ask user whether to select files or .d folders
        Alert choiceAlert = new Alert(Alert.AlertType.CONFIRMATION);
        choiceAlert.setTitle("Select Raw Data");
        choiceAlert.setHeaderText("What would you like to select?");
        choiceAlert.setContentText("Choose your option:");
        
        ButtonType filesButton = new ButtonType("Files (.raw, .mzML, .wiff)");
        ButtonType foldersButton = new ButtonType(".d Folders");
        ButtonType cancelButton = new ButtonType("Cancel", ButtonBar.ButtonData.CANCEL_CLOSE);
        
        choiceAlert.getButtonTypes().setAll(filesButton, foldersButton, cancelButton);
        
        choiceAlert.showAndWait().ifPresent(response -> {
            if (response == filesButton) {
                selectFiles();
            } else if (response == foldersButton) {
                selectDFolders();
            }
        });
    }
    
    private void selectFiles() {
        FileChooser fileChooser = new FileChooser();
        fileChooser.setTitle("Select Raw Files");
        fileChooser.getExtensionFilters().addAll(
            new FileChooser.ExtensionFilter("Raw Files", "*.raw", "*.mzML", "*.wiff"),
            new FileChooser.ExtensionFilter("All Files", "*.*")
        );
        List<File> files = fileChooser.showOpenMultipleDialog(getScene().getWindow());
        if (files != null) {
            for (File file : files) {
                addRawFileOrFolder(file);
            }
            updatePatternTable();
        }
    }
    
    private void selectDFolders() {
        // Step 1: Ask user to select parent directory
        DirectoryChooser directoryChooser = new DirectoryChooser();
        directoryChooser.setTitle("Select Parent Directory Containing .d Folders");
        File parentDir = directoryChooser.showDialog(getScene().getWindow());
        
        if (parentDir == null) {
            return; // User cancelled
        }
        
        // Step 2: Find all .d folders in the parent directory
        File[] allFiles = parentDir.listFiles();
        List<File> dFolders = new ArrayList<>();
        
        if (allFiles != null) {
            for (File file : allFiles) {
                if (file.isDirectory() && file.getName().toLowerCase().endsWith(".d")) {
                    dFolders.add(file);
                }
            }
        }
        
        // Step 3: Check if any .d folders were found
        if (dFolders.isEmpty()) {
            Alert noFoldersAlert = new Alert(Alert.AlertType.INFORMATION);
            noFoldersAlert.setTitle("No .d Folders Found");
            noFoldersAlert.setHeaderText("No .d folders in selected directory");
            noFoldersAlert.setContentText("No folders ending with '.d' were found in:\n" + parentDir.getAbsolutePath());
            noFoldersAlert.showAndWait();
            return;
        }
        
        // Step 4: Show dialog with checkboxes for all .d folders
        Dialog<List<File>> selectionDialog = new Dialog<>();
        selectionDialog.setTitle("Select .d Folders");
        selectionDialog.setHeaderText("Found " + dFolders.size() + " .d folder(s) in:\n" + parentDir.getAbsolutePath());
        
        // Create content
        VBox content = new VBox(10);
        content.setPadding(new Insets(10));
        
        Label instructionLabel = new Label("Select the .d folders you want to add:");
        instructionLabel.setStyle("-fx-font-weight: bold;");
        
        // Create scrollable list of checkboxes
        VBox checkBoxContainer = new VBox(5);
        List<CheckBox> checkBoxes = new ArrayList<>();
        
        for (File dFolder : dFolders) {
            CheckBox checkBox = new CheckBox(dFolder.getName());
            checkBox.setUserData(dFolder);
            checkBoxes.add(checkBox);
            checkBoxContainer.getChildren().add(checkBox);
        }
        
        ScrollPane scrollPane = new ScrollPane(checkBoxContainer);
        scrollPane.setFitToWidth(true);
        scrollPane.setPrefHeight(Math.min(300, dFolders.size() * 30 + 20));
        scrollPane.setStyle("-fx-border-color: #cccccc; -fx-border-radius: 3;");
        
        // Select All / Deselect All buttons
        HBox selectButtonBox = new HBox(10);
        Button selectAllButton = new Button("Select All");
        Button deselectAllButton = new Button("Deselect All");
        
        selectAllButton.setOnAction(e -> {
            for (CheckBox cb : checkBoxes) {
                cb.setSelected(true);
            }
        });
        
        deselectAllButton.setOnAction(e -> {
            for (CheckBox cb : checkBoxes) {
                cb.setSelected(false);
            }
        });
        
        selectButtonBox.getChildren().addAll(selectAllButton, deselectAllButton);
        
        content.getChildren().addAll(instructionLabel, scrollPane, selectButtonBox);
        selectionDialog.getDialogPane().setContent(content);
        
        // Add OK and Cancel buttons
        ButtonType okButton = new ButtonType("Add Selected", ButtonBar.ButtonData.OK_DONE);
        ButtonType cancelButton = new ButtonType("Cancel", ButtonBar.ButtonData.CANCEL_CLOSE);
        selectionDialog.getDialogPane().getButtonTypes().addAll(okButton, cancelButton);
        
        // Set result converter
        selectionDialog.setResultConverter(buttonType -> {
            if (buttonType == okButton) {
                List<File> selectedFolders = new ArrayList<>();
                for (CheckBox cb : checkBoxes) {
                    if (cb.isSelected()) {
                        selectedFolders.add((File) cb.getUserData());
                    }
                }
                return selectedFolders;
            }
            return null;
        });
        
        // Show dialog and process results
        selectionDialog.showAndWait().ifPresent(selectedFolders -> {
            if (!selectedFolders.isEmpty()) {
                for (File folder : selectedFolders) {
                    addRawFileOrFolder(folder);
                }
                updatePatternTable();
                
                // Show confirmation
                Alert confirmAlert = new Alert(Alert.AlertType.INFORMATION);
                confirmAlert.setTitle("Folders Added");
                confirmAlert.setHeaderText("Successfully added " + selectedFolders.size() + " .d folder(s)");
                confirmAlert.showAndWait();
            }
        });
    }
    
    private void addRawFileOrFolder(File file) {
        if (!fileSet.contains(file.getAbsolutePath())) {
            String path = file.getAbsolutePath();
            fileSet.add(path);
            
            // Set experiment name to filename without extension (like Swing)
            String name = file.getName();
            if (name.contains(".")) {
                name = name.substring(0, name.lastIndexOf("."));
            }
            
            experimentData.add(new ExperimentRow(path, name));
            
            // Set output directory to parent + "\\MetaPilot"
            outputDirField.setText(file.getParentFile().getAbsolutePath() + "\\MetaPilot");
        }
    }

    private void removeSelectedExperiment() {
        ExperimentRow selected = experimentTable.getSelectionModel().getSelectedItem();
        if (selected != null) {
            String path = selected.getFilePath();
            if (path != null && !path.isEmpty()) {
                fileSet.remove(path);
            }
            experimentData.remove(selected);
            updatePatternTable();
        }
    }

    private void removeAllExperiments() {
        experimentData.clear();
        fileSet.clear();
        updatePatternTable();
    }

    private void selectOutputDirectory() {
        DirectoryChooser directoryChooser = new DirectoryChooser();
        directoryChooser.setTitle("Select Output Directory");
        File dir = directoryChooser.showDialog(getScene().getWindow());
        if (dir != null) {
            outputDirField.setText(dir.getAbsolutePath());
        }
    }

    public ObservableList<ExperimentRow> getExperimentData() {
        return experimentData;
    }

    public String getOutputDirectory() {
        return outputDirField.getText();
    }
    
    public void setOutputDirectory(String outputDirectory) {
        if (outputDirectory != null) {
            outputDirField.setText(outputDirectory);
        }
    }
    
    public void setExperimentData(List<ExperimentRow> experimentData) {
        if (experimentData != null) {
            this.experimentData.clear();
            this.experimentData.addAll(experimentData);
            
            // Update file set for duplicate checking
            fileSet.clear();
            for (ExperimentRow row : experimentData) {
                String path = row.getFilePath();
                if (path != null && !path.isEmpty()) {
                    fileSet.add(path);
                }
            }
            
            updatePatternTable();
        }
    }
    
    // Method to get raw file names for metadata panel updates
    public List<String> getRawFileNames() {
        List<String> fileNames = new ArrayList<>();
        for (ExperimentRow row : experimentData) {
            String filePath = row.getFilePath();
            if (filePath != null && !filePath.isEmpty()) {
                String fileNameOnly = new File(filePath).getName();
                fileNames.add(fileNameOnly);
            }
        }
        return fileNames;
    }
    
    // Method to get experiment names for metadata panel updates
    public List<String> getExperimentNames() {
        List<String> expNames = new ArrayList<>();
        for (ExperimentRow row : experimentData) {
            String expName = row.getExperimentName();
            if (expName != null && !expName.isEmpty()) {
                expNames.add(expName);
            }
        }
        return expNames;
    }
    
    // Method to get fraction values
    public List<String> getFractionValues() {
        List<String> fractions = new ArrayList<>();
        for (ExperimentRow row : experimentData) {
            String fraction = row.getFraction();
            if (fraction != null && !fraction.isEmpty()) {
                fractions.add(fraction);
            }
        }
        return fractions;
    }

    public static class ExperimentRow {
        private final StringProperty filePath;
        private final StringProperty experimentName;
        private final StringProperty fraction;

        public ExperimentRow(String filePath, String experimentName) {
            this.filePath = new SimpleStringProperty(filePath == null ? "" : filePath);
            this.experimentName = new SimpleStringProperty(experimentName == null ? "" : experimentName);
            this.fraction = new SimpleStringProperty("1"); // Default fraction value
        }

        public String getFilePath() { return filePath.get(); }
        public void setFilePath(String filePath) { this.filePath.set(filePath == null ? "" : filePath); }
        public StringProperty filePathProperty() { return filePath; }

        public String getExperimentName() { return experimentName.get(); }
        public void setExperimentName(String experimentName) { 
            this.experimentName.set(experimentName == null ? "" : experimentName); 
        }
        public StringProperty experimentNameProperty() { return experimentName; }

        public String getFraction() { return fraction.get(); }
        public void setFraction(String fraction) { 
            this.fraction.set(fraction == null ? "1" : fraction); 
        }
        public StringProperty fractionProperty() { return fraction; }
        
        /**
         * Convert this ExperimentRow to a String array for parameter storage.
         * Order: [filePath, experimentName, fraction]
         */
        public String[] toArray() {
            return new String[]{getFilePath(), getExperimentName(), getFraction()};
        }
        
        /**
         * Create an ExperimentRow from a String array.
         * Expected order: [filePath, experimentName, fraction]
         */
        public static ExperimentRow fromArray(String[] data) {
            if (data == null || data.length < 2) {
                throw new IllegalArgumentException("Array must contain at least filePath and experimentName");
            }
            ExperimentRow row = new ExperimentRow(data[0], data[1]);
            if (data.length >= 3 && data[2] != null) {
                row.setFraction(data[2]);
            }
            return row;
        }
    }
    
    public static class PatternRow {
        private final BooleanProperty include;
        private final StringProperty content;

        public PatternRow(boolean include, String content) {
            this.include = new SimpleBooleanProperty(include);
            this.content = new SimpleStringProperty(content);
        }

        public boolean isInclude() { return include.get(); }
        public void setInclude(boolean include) { this.include.set(include); }
        public BooleanProperty includeProperty() { return include; }

        public String getContent() { return content.get(); }
        public void setContent(String content) { this.content.set(content); }
        public StringProperty contentProperty() { return content; }
    }
} 