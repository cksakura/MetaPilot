package uk.ac.quadram.metapilot.task.common.gui;

import javafx.geometry.Insets;
import javafx.scene.control.*;
import javafx.scene.layout.*;
import javafx.stage.FileChooser;
import javafx.scene.text.Font;
import javafx.scene.text.FontWeight;
import javafx.scene.paint.Color;
import java.io.File;

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
        microbiomeDBField.setPrefWidth(400);
        
        hostDBField = new TextField();
        hostDBField.setPromptText("Select host database file");
        hostDBField.setPrefWidth(400);
        hostDBField.setDisable(true);
        
        // Initialize checkboxes and buttons
        appendHostDBCheckBox = new CheckBox("Append Host Database");
        browseMicrobiomeDBButton = new Button("Browse");
        browseHostDBButton = new Button("Browse");
        browseHostDBButton.setDisable(true);
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
    }
    
    /**
     * Creates the database selection section.
     * 
     * @return VBox containing database selection components
     */
    private VBox createDatabaseSection() {
        VBox section = new VBox(10);
        section.setPadding(new Insets(10));
        section.setStyle("-fx-background-color: #f8f8f8; -fx-border-color: #cccccc; -fx-border-radius: 5; -fx-border-width: 1;");
        
        // Microbiome Database
        HBox microbiomeDBBox = new HBox(10);
        microbiomeDBBox.getChildren().addAll(
            new Label("Microbiome Database:"), 
            microbiomeDBField, 
            browseMicrobiomeDBButton
        );
        
        // Append Host Database
        HBox appendHostBox = new HBox(10);
        appendHostBox.getChildren().add(appendHostDBCheckBox);
        
        // Host Database
        HBox hostDBBox = new HBox(10);
        hostDBBox.getChildren().addAll(
            new Label("Host Database:"), 
            hostDBField, 
            browseHostDBButton
        );
        
        section.getChildren().addAll(microbiomeDBBox, appendHostBox, hostDBBox);
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
}
