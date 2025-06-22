package uk.ac.quadram.metapilot.task.mag.gui;

import javafx.application.Platform;
import javafx.beans.property.SimpleBooleanProperty;
import javafx.beans.property.SimpleIntegerProperty;
import javafx.beans.property.SimpleStringProperty;
import javafx.concurrent.Task;
import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.control.*;
import javafx.scene.control.cell.TextFieldTableCell;
import javafx.scene.layout.*;
import javafx.stage.DirectoryChooser;
import javafx.stage.FileChooser;
import javafx.scene.text.Font;
import javafx.scene.text.FontWeight;

import java.io.File;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import uk.ac.quadram.metapilot.core.MGYG.MGYGFtpTask;
import uk.ac.quadram.metapilot.core.MGYG.MGYGFtpTaskJavaFX;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;
import uk.ac.quadram.metapilot.task.mag.MagDbItem.MagVersion;

/**
 * JavaFX implementation of the MagDbConfigPanel
 * Provides the same functionality as the Swing version but with JavaFX components
 */
public class MagDbConfigPanelJavaFX extends BorderPane {

    private MagDbItem[][] magDbItems;
    private TableView<MagDbItem> mgnifyTable;
    private TableView<MagDbItem> localTable;
    
    // MGnify panel controls
    private ComboBox<MagDbItem> mgnifyComboBox;
    private TextField spCountTextField;
    private TextField latestVerTextField;
    private TextField mgnifyOldVerTextField;
    private TextField mgnifyRepoTextField;
    private Button downloadButton;
    private ProgressBar mgnifyDLBar;
    private TextArea mgnifyConsoleArea;
    
    // Local panel controls
    private TextField localNameTextField;
    private TextField identifierTextField;
    private TextField localRepoTextField;
    private TextField localDbTextField;
    private TextField localFuncTextField;
    private TextField localTaxaTextField;
    
    private SimpleDateFormat format = new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss");

    /**
     * Create the panel with the given database items
     */
    public MagDbConfigPanelJavaFX(MagDbItem[][] magDbItems) {
        this.magDbItems = magDbItems;
        setPadding(new Insets(10));
        
        // Create the main content panels
        TabPane tabPane = new TabPane();
        
        // MGnify Tab
        Tab mgnifyTab = new Tab("MGnify Database");
        mgnifyTab.setClosable(false);
        mgnifyTab.setContent(createMgnifyPanel());
        
        // Local Tab
        Tab localTab = new Tab("Local MAG Database");
        localTab.setClosable(false);
        localTab.setContent(createLocalPanel());
        
        tabPane.getTabs().addAll(mgnifyTab, localTab);
        
        // Set the tab pane as the center content
        setCenter(tabPane);
    }
    
    /**
     * Creates the MGnify database panel
     */
    private Pane createMgnifyPanel() {
        GridPane gridPane = new GridPane();
        gridPane.setHgap(20);
        gridPane.setVgap(10);
        gridPane.setPadding(new Insets(15));
        
        // Left side: MGnify database list
        VBox leftPanel = createMGnifyListPanel();
        
        // Right side: MGnify download panel
        VBox rightPanel = createMGnifyDownloadPanel();
        
        // Add panels to grid
        gridPane.add(leftPanel, 0, 0);
        gridPane.add(rightPanel, 1, 0);
        
        // Allow panels to grow to fill available space
        GridPane.setHgrow(leftPanel, Priority.ALWAYS);
        GridPane.setHgrow(rightPanel, Priority.ALWAYS);
        GridPane.setVgrow(leftPanel, Priority.ALWAYS);
        GridPane.setVgrow(rightPanel, Priority.ALWAYS);
        
        return gridPane;
    }
    
    /**
     * Creates the MGnify database list panel
     */
    private VBox createMGnifyListPanel() {
        VBox panel = new VBox(10);
        panel.setPadding(new Insets(10));
        panel.setStyle("-fx-border-color: #cccccc; -fx-border-radius: 5;");
        
        Label titleLabel = new Label("MGnify Database List");
        titleLabel.setFont(Font.font("System", FontWeight.BOLD, 14));
        
        // Create table view for MGnify databases
        mgnifyTable = new TableView<>();
        mgnifyTable.setEditable(false);
        mgnifyTable.setPlaceholder(new Label("No MGnify databases available"));
        
        // Create columns
        TableColumn<MagDbItem, Boolean> availableColumn = new TableColumn<>("Available");
        availableColumn.setCellValueFactory(data -> 
            new SimpleBooleanProperty(data.getValue().isAvailable()));
        availableColumn.setCellFactory(col -> new TableCell<MagDbItem, Boolean>() {
            private final CheckBox checkBox = new CheckBox();
            
            {
                // Make checkbox non-editable
                checkBox.setDisable(true);
                // Center the checkbox in the cell
                setAlignment(Pos.CENTER);
            }
            
            @Override
            protected void updateItem(Boolean item, boolean empty) {
                super.updateItem(item, empty);
                if (empty || item == null) {
                    setGraphic(null);
                } else {
                    checkBox.setSelected(item);
                    setGraphic(checkBox);
                }
            }
        });
        availableColumn.setPrefWidth(80);
        
        TableColumn<MagDbItem, String> idColumn = new TableColumn<>("Catalogue ID");
        idColumn.setCellValueFactory(data -> 
            new SimpleStringProperty(data.getValue().getCatalogueID()));
        idColumn.setPrefWidth(120);
        
        TableColumn<MagDbItem, String> nameColumn = new TableColumn<>("Identifier");
        nameColumn.setCellValueFactory(data -> 
            new SimpleStringProperty(data.getValue().getIdentifier()));
        nameColumn.setPrefWidth(150);
        
        TableColumn<MagDbItem, Number> speciesColumn = new TableColumn<>("Species Count");
        speciesColumn.setCellValueFactory(data -> 
            new SimpleIntegerProperty(data.getValue().getSpeciesCount()));
        speciesColumn.setPrefWidth(100);
        
        // Add columns to table
        mgnifyTable.getColumns().addAll(availableColumn, idColumn, nameColumn, speciesColumn);
        
        // Add data to table
        loadMgnifyTableData();
        
        // Set selection handler
        mgnifyTable.getSelectionModel().selectedItemProperty().addListener((obs, oldVal, newVal) -> {
            if (newVal != null) {
                updateMgnifySelection(newVal);
            }
        });
        
        // Add Check MGnify button and progress bar
        HBox buttonBox = new HBox(10);
        Button checkButton = new Button("Check MGnify");
        ProgressBar progressBar = new ProgressBar(0);
        progressBar.setPrefWidth(300);
        buttonBox.getChildren().addAll(checkButton, progressBar);
        
        // Set up the Check MGnify action
        checkButton.setOnAction(e -> {
            checkButton.setDisable(true);
            progressBar.setProgress(ProgressBar.INDETERMINATE_PROGRESS);
            
            // Create MGYGFtpTaskJavaFX for repository discovery
            MGYGFtpTaskJavaFX checkTask = new MGYGFtpTaskJavaFX();
            
            // Bind progress bar to task progress
            progressBar.progressProperty().bind(checkTask.progressProperty());
            
            // Handle task completion
            checkTask.setOnSucceeded(event -> {
                checkButton.setDisable(false);
                progressBar.progressProperty().unbind();
                progressBar.setProgress(0);
                
                // Update the table with discovered repositories
                MagDbItem[] repositories = checkTask.getRepositoryItems();
                if (repositories != null && repositories.length > 0) {
                    // Update magDbItems array
                    if (magDbItems == null) {
                        magDbItems = new MagDbItem[2][];
                    }
                    magDbItems[0] = repositories;
                    
                    // Refresh the table
                    Platform.runLater(() -> {
                        loadMgnifyTableData();
                        loadMgnifyComboBox();
                    });
                    
                    showAlert(Alert.AlertType.INFORMATION,
                            "Check Complete", 
                            "Found " + repositories.length + " MGnify databases");
                } else {
                    showAlert(Alert.AlertType.WARNING,
                            "Check Complete", 
                            "No MGnify databases found");
                }
            });
            
            checkTask.setOnFailed(event -> {
                checkButton.setDisable(false);
                progressBar.progressProperty().unbind();
                progressBar.setProgress(0);
                
                Throwable exception = checkTask.getException();
                showAlert(Alert.AlertType.ERROR,
                        "Check Failed", 
                        "Failed to check MGnify databases: " + 
                        (exception != null ? exception.getMessage() : "Unknown error"));
            });
            
            checkTask.setOnCancelled(event -> {
                checkButton.setDisable(false);
                progressBar.progressProperty().unbind();
                progressBar.setProgress(0);
                showAlert(Alert.AlertType.INFORMATION,
                        "Check Cancelled", 
                        "MGnify database check was cancelled");
            });
            
            // Start the task
            Thread thread = new Thread(checkTask);
            thread.setDaemon(true);
            thread.start();
        });
        
        // Make table take up all available space
        VBox.setVgrow(mgnifyTable, Priority.ALWAYS);
        
        // Add all components to panel
        panel.getChildren().addAll(titleLabel, mgnifyTable, buttonBox);
        
        return panel;
    }
    
    /**
     * Loads data into the MGnify table
     */
    private void loadMgnifyTableData() {
        mgnifyTable.getItems().clear();
        
        if (magDbItems != null && magDbItems[0] != null) {
            for (MagDbItem item : magDbItems[0]) {
                mgnifyTable.getItems().add(item);
            }
        }
        
        // If any items in the table, select the first one
        if (!mgnifyTable.getItems().isEmpty()) {
            mgnifyTable.getSelectionModel().select(0);
            updateMgnifySelection(mgnifyTable.getItems().get(0));
        }
    }
    
    /**
     * Creates the MGnify download panel
     */
    private VBox createMGnifyDownloadPanel() {
        VBox panel = new VBox(10);
        panel.setPadding(new Insets(10));
        panel.setStyle("-fx-border-color: #cccccc; -fx-border-radius: 5;");
        
        Label titleLabel = new Label("Download MGnify Database");
        titleLabel.setFont(Font.font("System", FontWeight.BOLD, 14));
        
        // Create the form grid
        GridPane formGrid = new GridPane();
        formGrid.setHgap(10);
        formGrid.setVgap(10);
        formGrid.setPadding(new Insets(5));
        
        // Database selection
        Label nameLabel = new Label("Name:");
        mgnifyComboBox = new ComboBox<>();
        mgnifyComboBox.setPrefWidth(200);
        loadMgnifyComboBox();
        
        // Species count
        Label speciesLabel = new Label("Species count:");
        spCountTextField = new TextField();
        spCountTextField.setEditable(false);
        
        // Latest version
        Label latestVerLabel = new Label("Latest version:");
        latestVerTextField = new TextField();
        latestVerTextField.setEditable(false);
        
        // Local version
        Label localVerLabel = new Label("Local version:");
        mgnifyOldVerTextField = new TextField();
        mgnifyOldVerTextField.setEditable(false);
        
        // Repository path
        Label repoLabel = new Label("Local repository:");
        mgnifyRepoTextField = new TextField();
        Button browseButton = new Button("Browse");
        
        // Set up browse button action
        browseButton.setOnAction(e -> {
            DirectoryChooser dirChooser = new DirectoryChooser();
            dirChooser.setTitle("Select MGnify Repository Directory");
            
            // If there's already a path, start from there
            if (!mgnifyRepoTextField.getText().isEmpty()) {
                File dir = new File(mgnifyRepoTextField.getText());
                if (dir.getParentFile() != null && dir.getParentFile().exists()) {
                    dirChooser.setInitialDirectory(dir.getParentFile());
                }
            }
            
            File selectedDir = dirChooser.showDialog(getScene().getWindow());
            if (selectedDir != null) {
                MagDbItem item = mgnifyComboBox.getValue();
                if (item != null) {
                    // Check if directory name matches catalogue ID
                    String fileName = selectedDir.getName();
                    if (!fileName.equals(item.getCatalogueID())) {
                        selectedDir = new File(selectedDir, item.getCatalogueID());
                    }
                    
                    mgnifyRepoTextField.setText(selectedDir.getAbsolutePath());
                    item.setRepositoryFile(selectedDir);
                    
                    // Update download button enablement
                    updateDownloadButton();
                }
            }
        });
        
        // Download button
        downloadButton = new Button("Download/Update");
        downloadButton.setDisable(true);
        mgnifyDLBar = new ProgressBar(0);
        mgnifyDLBar.setPrefWidth(350);
        
        // Console output area
        mgnifyConsoleArea = new TextArea();
        mgnifyConsoleArea.setEditable(false);
        mgnifyConsoleArea.setPrefHeight(200);
        
        // Set up download button action
        downloadButton.setOnAction(e -> {
            MagDbItem selectedItem = mgnifyComboBox.getValue();
            if (selectedItem == null) {
                showAlert(Alert.AlertType.WARNING, 
                         "Selection Required", 
                         "Please select a MGnify database to download");
                return;
            }
            
            File repoFile = selectedItem.getRepositoryFile();
            if (repoFile == null || !repoFile.getParentFile().exists()) {
                showAlert(Alert.AlertType.WARNING, 
                         "Repository Required", 
                         "Please specify a valid local repository path");
                return;
            }
            
            // Disable UI during download
            downloadButton.setDisable(true);
            browseButton.setDisable(true);
            mgnifyComboBox.setDisable(true);
            
            // Create download task using MGYGFtpTaskJavaFX
            MGYGFtpTaskJavaFX downloadTask = new MGYGFtpTaskJavaFX(selectedItem);
            
            // Bind progress and messages
            mgnifyDLBar.progressProperty().bind(downloadTask.progressProperty());
            mgnifyConsoleArea.textProperty().bind(downloadTask.messageProperty());
            
            // Handle task completion
            downloadTask.setOnSucceeded(event -> {
                // Re-enable UI
                downloadButton.setDisable(false);
                browseButton.setDisable(false);
                mgnifyComboBox.setDisable(false);
                
                // Unbind properties
                mgnifyDLBar.progressProperty().unbind();
                mgnifyConsoleArea.textProperty().unbind();
                
                // Set database as available
                selectedItem.setAvailable();
                selectedItem.setCurrentFile(repoFile);
                
                // Update table to show database is now available
                Platform.runLater(() -> {
                    mgnifyTable.refresh();
                });
                
                // Show success message
                showAlert(Alert.AlertType.INFORMATION, 
                         "Download Complete", 
                         "The MGnify database was downloaded successfully");
            });
            
            downloadTask.setOnFailed(event -> {
                // Re-enable UI
                downloadButton.setDisable(false);
                browseButton.setDisable(false);
                mgnifyComboBox.setDisable(false);
                
                // Unbind properties
                mgnifyDLBar.progressProperty().unbind();
                mgnifyConsoleArea.textProperty().unbind();
                
                // Show error message
                Throwable exception = downloadTask.getException();
                showAlert(Alert.AlertType.ERROR, 
                         "Download Failed", 
                         "Failed to download the MGnify database: " + 
                         (exception != null ? exception.getMessage() : "Unknown error"));
            });
            
            downloadTask.setOnCancelled(event -> {
                // Re-enable UI
                downloadButton.setDisable(false);
                browseButton.setDisable(false);
                mgnifyComboBox.setDisable(false);
                
                // Unbind properties
                mgnifyDLBar.progressProperty().unbind();
                mgnifyConsoleArea.textProperty().unbind();
                
                showAlert(Alert.AlertType.INFORMATION,
                        "Download Cancelled", 
                        "MGnify database download was cancelled");
            });
            
            // Start the task
            Thread thread = new Thread(downloadTask);
            thread.setDaemon(true);
            thread.start();
        });
        
        // Set up combobox selection handler
        mgnifyComboBox.setOnAction(e -> {
            MagDbItem item = mgnifyComboBox.getValue();
            if (item != null) {
                updateMgnifySelection(item);
            }
        });
        
        // Add form controls to grid
        int row = 0;
        formGrid.add(nameLabel, 0, row);
        formGrid.add(mgnifyComboBox, 1, row);
        formGrid.add(speciesLabel, 2, row);
        formGrid.add(spCountTextField, 3, row);
        
        row++;
        formGrid.add(latestVerLabel, 0, row);
        formGrid.add(latestVerTextField, 1, row);
        formGrid.add(localVerLabel, 2, row);
        formGrid.add(mgnifyOldVerTextField, 3, row);
        
        row++;
        formGrid.add(repoLabel, 0, row);
        formGrid.add(mgnifyRepoTextField, 1, row, 2, 1);
        formGrid.add(browseButton, 3, row);
        
        row++;
        formGrid.add(downloadButton, 0, row);
        formGrid.add(mgnifyDLBar, 1, row, 3, 1);
        
        // Make the console area expand to fill available space
        VBox.setVgrow(mgnifyConsoleArea, Priority.ALWAYS);
        
        // Add components to panel
        panel.getChildren().addAll(titleLabel, formGrid, mgnifyConsoleArea);
        
        return panel;
    }
    
    /**
     * Update UI when a MGnify database is selected
     */
    private void updateMgnifySelection(MagDbItem selectedItem) {
        if (selectedItem != null) {
            // Update combobox selection if available
            if (mgnifyComboBox != null && mgnifyComboBox.getValue() != selectedItem) {
                mgnifyComboBox.setValue(selectedItem);
            }
            
            // Update text fields
            if (spCountTextField != null) {
                spCountTextField.setText(String.valueOf(selectedItem.getSpeciesCount()));
            }
            
            if (latestVerTextField != null) {
                latestVerTextField.setText(selectedItem.getCatalogueVersion());
            }
            
            if (mgnifyOldVerTextField != null) {
                mgnifyOldVerTextField.setText(selectedItem.getOldVersion());
            }
            
            // Set repository path if available
            if (mgnifyRepoTextField != null) {
                File repoFile = selectedItem.getRepositoryFile();
                if (repoFile != null && repoFile.exists()) {
                    mgnifyRepoTextField.setText(repoFile.getAbsolutePath());
                } else {
                    mgnifyRepoTextField.clear();
                }
            }
            
            // Update download button state
            updateDownloadButton();
        }
    }
    
    /**
     * Update the state of the download button based on selection
     */
    private void updateDownloadButton() {
        if (downloadButton != null && mgnifyComboBox != null) {
            MagDbItem item = mgnifyComboBox.getValue();
            if (item != null) {
                downloadButton.setDisable(!item.updateAvailable());
            } else {
                downloadButton.setDisable(true);
            }
        }
    }
    
    /**
     * Load the MGnify combobox with database items
     */
    private void loadMgnifyComboBox() {
        mgnifyComboBox.getItems().clear();
        
        if (magDbItems != null && magDbItems[0] != null) {
            for (MagDbItem item : magDbItems[0]) {
                mgnifyComboBox.getItems().add(item);
            }
            
            if (!mgnifyComboBox.getItems().isEmpty()) {
                mgnifyComboBox.setValue(mgnifyComboBox.getItems().get(0));
                updateMgnifySelection(mgnifyComboBox.getValue());
            }
        }
    }
    
    /**
     * Creates the Local MAG database panel
     */
    private Pane createLocalPanel() {
        BorderPane panel = new BorderPane();
        panel.setPadding(new Insets(15));
        
        // Table view for local databases
        localTable = new TableView<>();
        localTable.setEditable(true);
        localTable.setPlaceholder(new Label("No local MAG databases configured"));
        
        // Create columns
        TableColumn<MagDbItem, String> idColumn = new TableColumn<>("Catalogue ID");
        idColumn.setCellValueFactory(data -> 
            new SimpleStringProperty(data.getValue().getCatalogueID()));
        idColumn.setCellFactory(TextFieldTableCell.forTableColumn());
        idColumn.setOnEditCommit(e -> {
            MagDbItem item = e.getRowValue();
            item.setCatalogueID(e.getNewValue());
        });
        idColumn.setPrefWidth(120);
        
        TableColumn<MagDbItem, String> nameColumn = new TableColumn<>("Identifier");
        nameColumn.setCellValueFactory(data -> 
            new SimpleStringProperty(data.getValue().getIdentifier()));
        nameColumn.setCellFactory(TextFieldTableCell.forTableColumn());
        nameColumn.setOnEditCommit(e -> {
            MagDbItem item = e.getRowValue();
            item.setIdentifier(e.getNewValue());
        });
        nameColumn.setPrefWidth(150);
        
        TableColumn<MagDbItem, String> pathColumn = new TableColumn<>("Path");
        pathColumn.setCellValueFactory(data -> {
            File file = data.getValue().getCurrentFile();
            return new SimpleStringProperty(file != null ? file.getAbsolutePath() : "");
        });
        pathColumn.setPrefWidth(300);
        
        // Actions column
        TableColumn<MagDbItem, Void> actionColumn = new TableColumn<>("Actions");
        actionColumn.setCellFactory(col -> new TableCell<MagDbItem, Void>() {
            private final Button browseButton = new Button("Browse");
            private final Button deleteButton = new Button("Delete");
            private final HBox buttonBox = new HBox(5, browseButton, deleteButton);
            
            {
                // Set up browse button
                browseButton.setOnAction(e -> {
                    MagDbItem item = getTableView().getItems().get(getIndex());
                    
                    FileChooser fileChooser = new FileChooser();
                    fileChooser.setTitle("Select Database File");
                    
                    // Set up extension filters
                    fileChooser.getExtensionFilters().addAll(
                        new FileChooser.ExtensionFilter("FASTA Files", "*.fasta", "*.fa"),
                        new FileChooser.ExtensionFilter("All Files", "*.*")
                    );
                    
                    // Set initial directory if available
                    File currentFile = item.getCurrentFile();
                    if (currentFile != null && currentFile.getParentFile() != null) {
                        fileChooser.setInitialDirectory(currentFile.getParentFile());
                    }
                    
                    File selectedFile = fileChooser.showOpenDialog(getScene().getWindow());
                    if (selectedFile != null) {
                        item.setCurrentFile(selectedFile);
                        getTableView().refresh();
                    }
                });
                
                // Set up delete button
                deleteButton.setOnAction(e -> {
                    MagDbItem item = getTableView().getItems().get(getIndex());
                    
                    Alert alert = new Alert(Alert.AlertType.CONFIRMATION);
                    alert.setTitle("Confirm Deletion");
                    alert.setHeaderText("Delete Database");
                    alert.setContentText("Are you sure you want to remove this database?");
                    
                    alert.showAndWait().ifPresent(result -> {
                        if (result == ButtonType.OK) {
                            getTableView().getItems().remove(item);
                        }
                    });
                });
            }
            
            @Override
            protected void updateItem(Void item, boolean empty) {
                super.updateItem(item, empty);
                setGraphic(empty ? null : buttonBox);
            }
        });
        actionColumn.setPrefWidth(150);
        
        // Add columns to table
        localTable.getColumns().addAll(idColumn, nameColumn, pathColumn, actionColumn);
        
        // Load data
        loadLocalTableData();
        
        // Add buttons at the bottom
        HBox buttonBox = new HBox(10);
        buttonBox.setPadding(new Insets(10, 0, 0, 0));
        
        Button addButton = new Button("Add New Database");
        addButton.setOnAction(e -> {
            // Add a new empty database
            MagDbItem newItem = new MagDbItem("", "New Database");
            
            if (magDbItems[1] == null) {
                magDbItems[1] = new MagDbItem[1];
                magDbItems[1][0] = newItem;
            } else {
                MagDbItem[] newArray = new MagDbItem[magDbItems[1].length + 1];
                System.arraycopy(magDbItems[1], 0, newArray, 0, magDbItems[1].length);
                newArray[magDbItems[1].length] = newItem;
                magDbItems[1] = newArray;
            }
            
            loadLocalTableData();
        });
        
        buttonBox.getChildren().add(addButton);
        
        // Set table as center, buttons at bottom
        panel.setCenter(localTable);
        panel.setBottom(buttonBox);
        
        return panel;
    }
    
    /**
     * Load data into the local MAG database table
     */
    private void loadLocalTableData() {
        localTable.getItems().clear();
        
        if (magDbItems != null && magDbItems[1] != null) {
            for (MagDbItem item : magDbItems[1]) {
                localTable.getItems().add(item);
            }
        }
    }
    
    /**
     * Helper method to show alerts
     */
    private void showAlert(Alert.AlertType type, String title, String message) {
        Alert alert = new Alert(type);
        alert.setTitle(title);
        alert.setHeaderText(null);
        alert.setContentText(message);
        alert.showAndWait();
    }
    
    /**
     * Get the current database configuration
     */
    public MagDbItem[][] getMagDbItems() {
        // Update local database array from table
        if (localTable != null && !localTable.getItems().isEmpty()) {
            magDbItems[1] = localTable.getItems().toArray(new MagDbItem[0]);
        }
        
        return magDbItems;
    }
    

} 