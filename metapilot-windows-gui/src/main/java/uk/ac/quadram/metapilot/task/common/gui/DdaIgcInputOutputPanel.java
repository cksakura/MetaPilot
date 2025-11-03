package uk.ac.quadram.metapilot.task.common.gui;

import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.control.*;
import javafx.scene.layout.*;
import javafx.stage.FileChooser;
import javafx.scene.text.Font;
import javafx.scene.text.FontWeight;
import javafx.scene.paint.Color;
import javafx.scene.control.Alert.AlertType;
import javafx.scene.control.Hyperlink;
import java.io.File;
import java.net.URI;
import java.io.InputStream;
import java.io.FileOutputStream;
import java.awt.Desktop;
import java.util.zip.GZIPInputStream;
import uk.ac.quadram.metapilot.task.par.MetaConstants;

/**
 * Unified Input/Output panel for DDA IGC workflows.
 * 
 * This component provides a consistent Input/Output interface across all four DDA IGC workflows:
 * - MaxQuant
 * - pFind
 * - Sage
 * - FragPipe
 * 
 * Features:
 * - File input selection via MetaPilotIOPanelFX
 * - Microbiome database selection
 * - Optional host database selection
 * - Consistent styling and behavior
 * 
 * @author MetaPilot Development Team
 * @version 1.0
 * @since MetaPilot 3.0
 */
public class DdaIgcInputOutputPanel extends VBox {
    
    // UI Components
    private MetaPilotIOPanelFX ioPanel;
    private TextField microbiomeDBField;
    private TextField hostDBField;
    private CheckBox appendHostDBCheckBox;
    private Button browseMicrobiomeDBButton;
    private Button browseHostDBButton;
    private Button downloadHumanGutButton;
    private Button downloadMouseGutButton;
    private Button downloadHumanSkinButton;
    private Button downloadPigGutButton;
    private Button downloadRatGutButton;
    private ProgressBar downloadProgressBar;
    private Label downloadStatusLabel;
    
    // Callback interface for parameter updates
    public interface ParameterUpdateCallback {
        void onParameterChanged();
    }
    
    private ParameterUpdateCallback parameterUpdateCallback;
    
    /**
     * Creates a new unified Input/Output panel for DDA IGC workflows.
     */
    public DdaIgcInputOutputPanel() {
        initializeComponents();
        setupLayout();
        setupEventHandlers();
    }
    
    /**
     * Creates a new unified Input/Output panel with parameter update callback.
     * 
     * @param callback Callback to invoke when parameters change
     */
    public DdaIgcInputOutputPanel(ParameterUpdateCallback callback) {
        this.parameterUpdateCallback = callback;
        initializeComponents();
        setupLayout();
        setupEventHandlers();
    }
    
    /**
     * Initialize all UI components.
     */
    private void initializeComponents() {
        // Initialize IO panel
        ioPanel = new MetaPilotIOPanelFX();
        
        // Initialize database fields
        microbiomeDBField = new TextField();
        microbiomeDBField.setPromptText("Select microbiome database file");
        microbiomeDBField.setPrefWidth(600);
        
        hostDBField = new TextField();
        hostDBField.setPromptText("Select host database file");
        hostDBField.setPrefWidth(600);
        hostDBField.setDisable(true);
        
        // Initialize checkboxes and buttons
        appendHostDBCheckBox = new CheckBox("Append Host Database");
        browseMicrobiomeDBButton = new Button("Browse");
        browseHostDBButton = new Button("Browse");
        browseHostDBButton.setDisable(true);
        
        // Initialize download buttons
        downloadHumanGutButton = new Button("Human Gut");
        downloadHumanGutButton.setStyle("-fx-background-color: #4CAF50; -fx-text-fill: white;");
        downloadMouseGutButton = new Button("Mouse Gut");
        downloadMouseGutButton.setStyle("-fx-background-color: #2196F3; -fx-text-fill: white;");
        downloadHumanSkinButton = new Button("Human Skin");
        downloadHumanSkinButton.setStyle("-fx-background-color: #FF9800; -fx-text-fill: white;");
        downloadPigGutButton = new Button("Pig Gut");
        downloadPigGutButton.setStyle("-fx-background-color: #9C27B0; -fx-text-fill: white;");
        downloadRatGutButton = new Button("Rat Gut");
        downloadRatGutButton.setStyle("-fx-background-color: #F44336; -fx-text-fill: white;");
        
        // Initialize download progress components
        downloadProgressBar = new ProgressBar(0);
        downloadProgressBar.setVisible(false);
        downloadStatusLabel = new Label("");
        downloadStatusLabel.setVisible(false);
    }
    
    /**
     * Setup the layout structure.
     */
    private void setupLayout() {
        setSpacing(10);
        setPadding(new Insets(10));
        
        // Create database section
        VBox databaseSection = createDatabaseSection();
        
        // Add all components to main layout
        getChildren().addAll(
            createSectionTitle("Input/Output Settings"),
            ioPanel,
            createSectionTitle("Database Settings"),
            databaseSection
        );
    }
    
    /**
     * Setup event handlers for all interactive components.
     */
    private void setupEventHandlers() {
        // Microbiome database field listener
        microbiomeDBField.textProperty().addListener((obs, oldVal, newVal) -> {
            if (parameterUpdateCallback != null) {
                parameterUpdateCallback.onParameterChanged();
            }
        });
        
        // Append host database checkbox
        appendHostDBCheckBox.setOnAction(e -> {
            updateHostDBVisibility();
            if (parameterUpdateCallback != null) {
                parameterUpdateCallback.onParameterChanged();
            }
        });
        
        // Browse buttons
        browseMicrobiomeDBButton.setOnAction(e -> browseMicrobiomeDatabase());
        browseHostDBButton.setOnAction(e -> browseHostDatabase());
        
        // Download buttons
        downloadHumanGutButton.setOnAction(e -> downloadHumanGutDatabase());
        downloadMouseGutButton.setOnAction(e -> downloadMouseGutDatabase());
        downloadHumanSkinButton.setOnAction(e -> downloadHumanSkinDatabase());
        downloadPigGutButton.setOnAction(e -> downloadPigGutDatabase());
        downloadRatGutButton.setOnAction(e -> downloadRatGutDatabase());
    }
    
    /**
     * Creates the database selection section.
     * 
     * @return VBox containing database selection components
     */
	private VBox createDatabaseSection() {
		VBox section = new VBox(10);
		section.setPadding(new Insets(10));
		section.setStyle(
				"-fx-background-color: #f8f8f8; -fx-border-color: #cccccc; -fx-border-radius: 5; -fx-border-width: 1;");

		// Microbiome Database
		HBox microbiomeDBBox = new HBox(10);
		microbiomeDBBox.setAlignment(Pos.CENTER_LEFT);
		microbiomeDBBox.getChildren().addAll(new Label("Microbiome Database:"), microbiomeDBField,
				browseMicrobiomeDBButton);

		// Popular Database Downloads
		VBox downloadSection = new VBox(8);
		downloadSection.setPadding(new Insets(10));
		downloadSection.setStyle(
				"-fx-background-color: #e8f5e8; -fx-border-color: #4CAF50; -fx-border-radius: 5; -fx-border-width: 1;");

		HBox dbTitleBax = new HBox(10);
		dbTitleBax.setAlignment(Pos.CENTER_LEFT);
		Label downloadTitle = new Label("Popular Microbiome Databases");
		downloadTitle.setStyle("-fx-font-weight: bold; -fx-text-fill: #2E7D32;");

		Hyperlink sourceLink = new Hyperlink("Source: China National GeneBank Database (CNGBdb)");
		sourceLink.setStyle("-fx-font-size: 14px; -fx-text-fill: #1976D2; -fx-font-style: italic;");
		sourceLink.setOnAction(e -> openCNGBWebsite());
		dbTitleBax.getChildren().addAll(downloadTitle, sourceLink);

		// First row of buttons
		HBox downloadButtonsRow = new HBox(10);
		downloadButtonsRow.setAlignment(Pos.CENTER_LEFT);
		downloadButtonsRow.getChildren().addAll(downloadHumanGutButton, downloadMouseGutButton, downloadHumanSkinButton,
				downloadPigGutButton, downloadRatGutButton);

		VBox downloadButtonsBox = new VBox(5);
		downloadButtonsBox.getChildren().addAll(downloadButtonsRow);

		// Database info labels
		Label humanGutInfo = new Label("Human Gut: ~1.6 GB • " + MetaConstants.humanGutProDb);
		humanGutInfo.setStyle("-fx-font-size: 12px; -fx-text-fill: #666666;");
		Label mouseGutInfo = new Label("Mouse Gut: ~420 MB • " + MetaConstants.mouseGutProDb);
		mouseGutInfo.setStyle("-fx-font-size: 12px; -fx-text-fill: #666666;");
		Label humanSkinInfo = new Label("Human Skin: ~1.0 GB • " + MetaConstants.humanSkinProDb);
		humanSkinInfo.setStyle("-fx-font-size: 12px; -fx-text-fill: #666666;");
		Label pigGutInfo = new Label("Pig Gut: ~750 MB • " + MetaConstants.pigGutProDb);
		pigGutInfo.setStyle("-fx-font-size: 12px; -fx-text-fill: #666666;");
		Label ratGutInfo = new Label("Rat Gut: ~400 MB • " + MetaConstants.ratGutProDb);
		ratGutInfo.setStyle("-fx-font-size: 12px; -fx-text-fill: #666666;");

		VBox infoBox = new VBox(1);
		infoBox.getChildren().addAll(humanGutInfo, mouseGutInfo, humanSkinInfo, pigGutInfo, ratGutInfo);

		downloadSection.getChildren().addAll(dbTitleBax, downloadButtonsBox, infoBox, downloadProgressBar,
				downloadStatusLabel);

		// Append Host Database
		HBox appendHostBox = new HBox(10);
		appendHostBox.getChildren().add(appendHostDBCheckBox);

		// Host Database
		HBox hostDBBox = new HBox(10);
		hostDBBox.setAlignment(Pos.CENTER_LEFT);
		hostDBBox.getChildren().addAll(new Label("Host Database:"), hostDBField, browseHostDBButton);

		section.getChildren().addAll(microbiomeDBBox, downloadSection, appendHostBox, hostDBBox);
		return section;
	}
    
    /**
     * Creates a section title with consistent styling.
     * 
     * @param title The title text
     * @return Label with section title styling
     */
    private Label createSectionTitle(String title) {
        Label titleLabel = new Label(title);
        titleLabel.setFont(Font.font("Arial", FontWeight.BOLD, 14));
        titleLabel.setTextFill(Color.DARKBLUE);
        titleLabel.setPadding(new Insets(10, 0, 5, 0));
        return titleLabel;
    }
    
    /**
     * Updates the visibility and enabled state of host database components.
     */
    private void updateHostDBVisibility() {
        boolean selected = appendHostDBCheckBox.isSelected();
        hostDBField.setDisable(!selected);
        browseHostDBButton.setDisable(!selected);
    }
    
    /**
     * Opens file chooser for microbiome database selection.
     */
    private void browseMicrobiomeDatabase() {
        FileChooser chooser = new FileChooser();
        chooser.setTitle("Select Microbiome Database File");
        chooser.getExtensionFilters().add(
            new FileChooser.ExtensionFilter("FASTA Files", "*.fasta", "*.fa")
        );
        
        File file = chooser.showOpenDialog(getScene().getWindow());
        if (file != null) {
            microbiomeDBField.setText(file.getAbsolutePath());
        }
    }
    
    /**
     * Opens file chooser for host database selection.
     */
    private void browseHostDatabase() {
        FileChooser chooser = new FileChooser();
        chooser.setTitle("Select Host Database File");
        chooser.getExtensionFilters().add(
            new FileChooser.ExtensionFilter("FASTA Files", "*.fasta", "*.fa")
        );
        
        File file = chooser.showOpenDialog(getScene().getWindow());
        if (file != null) {
            hostDBField.setText(file.getAbsolutePath());
        }
    }
    
    // Getters for accessing components and data
    
    /**
     * Gets the IO panel for file selection.
     * 
     * @return MetaPilotIOPanelFX instance
     */
    public MetaPilotIOPanelFX getIoPanel() {
        return ioPanel;
    }
    
    /**
     * Gets the microbiome database field.
     * 
     * @return TextField for microbiome database path
     */
    public TextField getMicrobiomeDBField() {
        return microbiomeDBField;
    }
    
    /**
     * Gets the host database field.
     * 
     * @return TextField for host database path
     */
    public TextField getHostDBField() {
        return hostDBField;
    }
    
    /**
     * Gets the append host database checkbox.
     * 
     * @return CheckBox for append host database option
     */
    public CheckBox getAppendHostDBCheckBox() {
        return appendHostDBCheckBox;
    }
    
    /**
     * Gets the browse microbiome database button.
     * 
     * @return Button for browsing microbiome database
     */
    public Button getBrowseMicrobiomeDBButton() {
        return browseMicrobiomeDBButton;
    }
    
    /**
     * Gets the browse host database button.
     * 
     * @return Button for browsing host database
     */
    public Button getBrowseHostDBButton() {
        return browseHostDBButton;
    }
    
    // Convenience methods for common operations
    
    /**
     * Gets the microbiome database path.
     * 
     * @return String path to microbiome database, or empty string if not set
     */
    public String getMicrobiomeDBPath() {
        return microbiomeDBField.getText();
    }
    
    /**
     * Sets the microbiome database path.
     * 
     * @param path Path to microbiome database
     */
    public void setMicrobiomeDBPath(String path) {
        microbiomeDBField.setText(path);
    }
    
    /**
     * Gets the host database path.
     * 
     * @return String path to host database, or empty string if not set
     */
    public String getHostDBPath() {
        return hostDBField.getText();
    }
    
    /**
     * Sets the host database path.
     * 
     * @param path Path to host database
     */
    public void setHostDBPath(String path) {
        hostDBField.setText(path);
    }
    
    /**
     * Checks if host database is enabled.
     * 
     * @return true if host database is enabled, false otherwise
     */
    public boolean isHostDBEnabled() {
        return appendHostDBCheckBox.isSelected();
    }
    
    /**
     * Sets the host database enabled state.
     * 
     * @param enabled true to enable host database, false to disable
     */
    public void setHostDBEnabled(boolean enabled) {
        appendHostDBCheckBox.setSelected(enabled);
        updateHostDBVisibility();
    }
    
    /**
     * Sets the parameter update callback.
     * 
     * @param callback Callback to invoke when parameters change
     */
    public void setParameterUpdateCallback(ParameterUpdateCallback callback) {
        this.parameterUpdateCallback = callback;
    }
    
    /**
     * Opens the CNGBdb Microbiome Database website in the default browser.
     */
    private void openCNGBWebsite() {
        try {
            Desktop.getDesktop().browse(new URI("https://db.cngb.org/microbiome/"));
        } catch (Exception e) {
            // Fallback: show alert with URL
            Alert alert = new Alert(AlertType.INFORMATION);
            alert.setTitle("CNGBdb Microbiome Database");
            alert.setHeaderText("China National GeneBank Database");
            alert.setContentText("Please visit: https://db.cngb.org/microbiome/\n\nThis website provides comprehensive microbiome databases including gene catalogs for various host-associated microbiomes.");
            alert.showAndWait();
        }
    }
    
    /**
     * Downloads the Human Gut Microbiome Database.
     */
    private void downloadHumanGutDatabase() {
        downloadDatabase(MetaConstants.humanGutProDbLink, "IGC.pep.gz", "Human Gut Microbiome Database");
    }
    
    /**
     * Downloads the Mouse Gut Microbiome Database.
     */
    private void downloadMouseGutDatabase() {
        downloadDatabase(MetaConstants.mouseGutProDbLink, "mouseGut_2.6M.GeneCatalog.pep.gz", "Mouse Gut Microbiome Database");
    }
    
    /**
     * Downloads the Human Skin Microbiome Database.
     */
    private void downloadHumanSkinDatabase() {
        downloadDatabase(MetaConstants.humanSkinProDbLink, "IHSMGC.pep.gz", "Human Skin Microbiome Database");
    }
    
    /**
     * Downloads the Pig Gut Microbiome Database.
     */
    private void downloadPigGutDatabase() {
        downloadDatabase(MetaConstants.pigGutProDbLink, "pigGut_7.7M.GeneSet.pep.gz", "Pig Gut Microbiome Database");
    }
    
    /**
     * Downloads the Rat Gut Microbiome Database.
     */
    private void downloadRatGutDatabase() {
        downloadDatabase(MetaConstants.ratGutProDbLink, "rat_geneset.pep.gz", "Rat Gut Microbiome Database");
    }
    
    /**
     * Downloads a database from the given URL.
     * 
     * @param urlString The URL to download from
     * @param fileName The local filename to save as
     * @param displayName The display name for progress messages
     */
    private void downloadDatabase(String urlString, String fileName, String displayName) {
        // Show file chooser to select download location
        FileChooser chooser = new FileChooser();
        chooser.setTitle("Select Download Location for " + displayName);
        // Change .pep.gz to .fasta for output
        String outputFileName = fileName.replace(".pep.gz", ".fasta");
        chooser.setInitialFileName(outputFileName);
        chooser.getExtensionFilters().add(
            new FileChooser.ExtensionFilter("FASTA Files", "*.fasta", "*.fa")
        );
        
        File downloadFile = chooser.showSaveDialog(getScene().getWindow());
        if (downloadFile == null) {
            return; // User cancelled
        }
        
        // Create temporary file for compressed download
        File compressedFile = new File(downloadFile.getParent(), downloadFile.getName() + ".gz");
        
        // Show progress UI
        downloadProgressBar.setVisible(true);
        downloadStatusLabel.setVisible(true);
        downloadProgressBar.setProgress(0);
        downloadStatusLabel.setText("Starting download...");
        
        // Disable download buttons during download
        downloadHumanGutButton.setDisable(true);
        downloadMouseGutButton.setDisable(true);
        downloadHumanSkinButton.setDisable(true);
        downloadPigGutButton.setDisable(true);
        downloadRatGutButton.setDisable(true);
        
        // Start download in background thread
        Thread downloadThread = new Thread(() -> {
            try {
                URI uri = URI.create(urlString);
                java.net.URL url = uri.toURL();
                
                // Get content length for progress calculation
                long contentLength = url.openConnection().getContentLengthLong();
                
                try (InputStream inputStream = url.openStream();
                     FileOutputStream outputStream = new FileOutputStream(compressedFile)) {
                    
                    byte[] buffer = new byte[8192];
                    int bytesRead;
                    long totalBytesRead = 0;
                    
                    while ((bytesRead = inputStream.read(buffer)) != -1) {
                        outputStream.write(buffer, 0, bytesRead);
                        totalBytesRead += bytesRead;
                        
                        // Update progress on JavaFX thread
                        final long currentBytes = totalBytesRead;
                        javafx.application.Platform.runLater(() -> {
                            if (contentLength > 0) {
                                double progress = (double) currentBytes / contentLength * 0.5; // 50% for download
                                downloadProgressBar.setProgress(progress);
                                downloadStatusLabel.setText(String.format("Downloading (compressed)... %.1f%% (%.1f MB / %.1f MB)", 
                                    progress * 100, 
                                    currentBytes / (1024.0 * 1024.0), 
                                    contentLength / (1024.0 * 1024.0)));
                            } else {
                                downloadStatusLabel.setText(String.format("Downloading (compressed)... %.1f MB", 
                                    currentBytes / (1024.0 * 1024.0)));
                            }
                        });
                    }
                }
                
                // Extract the downloaded file
                javafx.application.Platform.runLater(() -> {
                    downloadStatusLabel.setText("Extracting database...");
                });
                
                try (GZIPInputStream gzipInput = new GZIPInputStream(new java.io.FileInputStream(compressedFile));
                     FileOutputStream outputStream = new FileOutputStream(downloadFile)) {
                    
                    byte[] buffer = new byte[8192];
                    int bytesRead;
                    
                    while ((bytesRead = gzipInput.read(buffer)) != -1) {
                        outputStream.write(buffer, 0, bytesRead);
                    }
                }
                
                // Delete temporary compressed file
                compressedFile.delete();
                
                // Download completed successfully
                javafx.application.Platform.runLater(() -> {
                    downloadProgressBar.setProgress(1.0);
                    downloadStatusLabel.setText("Download completed successfully!");
                    
                    // Show success dialog
                    Alert alert = new Alert(AlertType.INFORMATION);
                    alert.setTitle("Download Complete");
                    alert.setHeaderText(displayName + " Downloaded Successfully");
                    alert.setContentText("The database has been saved to:\n" + downloadFile.getAbsolutePath() + 
                                       "\n\nWould you like to use this database for your analysis?");
                    
                    ButtonType useDatabaseButton = new ButtonType("Use This Database");
                    ButtonType cancelButton = new ButtonType("Cancel", ButtonBar.ButtonData.CANCEL_CLOSE);
                    alert.getButtonTypes().setAll(useDatabaseButton, cancelButton);
                    
                    alert.showAndWait().ifPresent(buttonType -> {
                        if (buttonType == useDatabaseButton) {
                            microbiomeDBField.setText(downloadFile.getAbsolutePath());
                            if (parameterUpdateCallback != null) {
                                parameterUpdateCallback.onParameterChanged();
                            }
                        }
                    });
                    
                    // Hide progress UI after a delay
                    new Thread(() -> {
                        try {
                            Thread.sleep(3000);
                            javafx.application.Platform.runLater(() -> {
                                downloadProgressBar.setVisible(false);
                                downloadStatusLabel.setVisible(false);
                                downloadHumanGutButton.setDisable(false);
                                downloadMouseGutButton.setDisable(false);
                                downloadHumanSkinButton.setDisable(false);
                                downloadPigGutButton.setDisable(false);
                                downloadRatGutButton.setDisable(false);
                            });
                        } catch (InterruptedException e) {
                            Thread.currentThread().interrupt();
                        }
                    }).start();
                });
                
            } catch (Exception e) {
                // Download failed
                javafx.application.Platform.runLater(() -> {
                    downloadProgressBar.setVisible(false);
                    downloadStatusLabel.setVisible(false);
                    downloadHumanGutButton.setDisable(false);
                    downloadMouseGutButton.setDisable(false);
                    downloadHumanSkinButton.setDisable(false);
                    downloadPigGutButton.setDisable(false);
                    downloadRatGutButton.setDisable(false);
                    
                    Alert alert = new Alert(AlertType.ERROR);
                    alert.setTitle("Download Failed");
                    alert.setHeaderText("Failed to download " + displayName);
                    alert.setContentText("Error: " + e.getMessage() + 
                                       "\n\nPlease check your internet connection and try again.");
                    alert.showAndWait();
                });
            }
        });
        
        downloadThread.setDaemon(true);
        downloadThread.start();
    }
}
