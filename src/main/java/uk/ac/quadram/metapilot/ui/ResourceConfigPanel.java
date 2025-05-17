package uk.ac.quadram.metapilot.ui;

import javafx.application.Platform;
import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.control.*;
import javafx.scene.layout.*;
import javafx.scene.paint.Color;
import javafx.scene.shape.Circle;
import javafx.stage.DirectoryChooser;
import javafx.stage.FileChooser;
import java.io.File;
import java.util.Map;

import uk.ac.quadram.metapilot.resources.MetaResourceIO;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters.ToolType;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters.ToolStatus;

/**
 * JavaFX panel for configuring resource paths
 */
public class ResourceConfigPanel extends BorderPane {
    
    private MetaResourceParameters resourceParams;
    private MetaResourceIO resourceIO;
    
    private TextField basePathField;
    private Map<ToolType, ToolPathControl> toolControls;
    
    /**
     * Create a new resource configuration panel
     */
    public ResourceConfigPanel() {
        // Initialize resource parameters and IO
        resourceIO = new MetaResourceIO();
        resourceParams = resourceIO.loadParameters();
        
        setPadding(new Insets(15));
        
        // Create content
        setTop(createHeaderPanel());
        setCenter(createToolsPanel());
        setBottom(createButtonPanel());
    }
    
    /**
     * Create the header panel with the base path control
     */
    private Pane createHeaderPanel() {
        VBox panel = new VBox(10);
        panel.setPadding(new Insets(0, 0, 15, 0));
        
        Label titleLabel = new Label("Resource Configuration");
        titleLabel.setStyle("-fx-font-size: 18px; -fx-font-weight: bold;");
        
        HBox basePathBox = new HBox(10);
        basePathBox.setAlignment(Pos.CENTER_LEFT);
        
        Label basePathLabel = new Label("Base Resources Path:");
        basePathField = new TextField(resourceParams.getBaseResourcePath());
        basePathField.setPrefWidth(400);
        
        Button browseButton = new Button("Browse...");
        browseButton.setOnAction(e -> {
            DirectoryChooser dirChooser = new DirectoryChooser();
            dirChooser.setTitle("Select Base Resources Directory");
            
            // Set initial directory if current path exists
            File currentDir = new File(basePathField.getText());
            if (currentDir.exists() && currentDir.isDirectory()) {
                dirChooser.setInitialDirectory(currentDir);
            }
            
            File selectedDir = dirChooser.showDialog(getScene().getWindow());
            if (selectedDir != null) {
                basePathField.setText(selectedDir.getAbsolutePath());
                updateBaseResourcePath(selectedDir.getAbsolutePath());
            }
        });
        
        basePathBox.getChildren().addAll(basePathLabel, basePathField, browseButton);
        
        Label instructionLabel = new Label(
            "Configure the paths to external tools used by MetaPilot. " +
            "Green indicators show that tools are available at the specified path."
        );
        instructionLabel.setWrapText(true);
        
        panel.getChildren().addAll(titleLabel, basePathBox, instructionLabel);
        
        return panel;
    }
    
    /**
     * Create the tools panel with all tool path controls
     */
    private ScrollPane createToolsPanel() {
        ScrollPane scrollPane = new ScrollPane();
        scrollPane.setFitToWidth(true);
        
        VBox toolsBox = new VBox(15);
        toolsBox.setPadding(new Insets(10));
        
        // Create DDA tools section
        TitledPane ddaSection = new TitledPane();
        ddaSection.setText("DDA Tools");
        ddaSection.setCollapsible(true);
        ddaSection.setExpanded(true);
        
        GridPane ddaGrid = new GridPane();
        ddaGrid.setHgap(20);
        ddaGrid.setVgap(15);
        
        int ddaRow = 0;
        for (ToolType toolType : ToolType.getDDATools()) {
            ToolPathControl control = new ToolPathControl(toolType);
            ddaGrid.add(control, 0, ddaRow++);
        }
        
        VBox ddaBox = new VBox(10);
        ddaBox.setPadding(new Insets(10));
        ddaBox.getChildren().add(ddaGrid);
        ddaSection.setContent(ddaBox);
        
        // Create DIA tools section
        TitledPane diaSection = new TitledPane();
        diaSection.setText("DIA Tools");
        diaSection.setCollapsible(true);
        diaSection.setExpanded(true);
        
        GridPane diaGrid = new GridPane();
        diaGrid.setHgap(20);
        diaGrid.setVgap(15);
        
        int diaRow = 0;
        for (ToolType toolType : ToolType.getDIATools()) {
            ToolPathControl control = new ToolPathControl(toolType);
            diaGrid.add(control, 0, diaRow++);
        }
        
        VBox diaBox = new VBox(10);
        diaBox.setPadding(new Insets(10));
        diaBox.getChildren().add(diaGrid);
        diaSection.setContent(diaBox);
        
        // Add both sections to the main box
        toolsBox.getChildren().addAll(ddaSection, diaSection);
        
        scrollPane.setContent(toolsBox);
        
        return scrollPane;
    }
    
    /**
     * Create the button panel with save and reset buttons
     */
    private Pane createButtonPanel() {
        HBox panel = new HBox(10);
        panel.setPadding(new Insets(15, 0, 0, 0));
        panel.setAlignment(Pos.CENTER_RIGHT);
        
        Button resetButton = new Button("Reset to Defaults");
        resetButton.setOnAction(e -> {
            Alert alert = new Alert(Alert.AlertType.CONFIRMATION);
            alert.setTitle("Reset Configuration");
            alert.setHeaderText("Reset to Default Paths");
            alert.setContentText("Are you sure you want to reset all tool paths to their defaults?");
            
            alert.showAndWait().ifPresent(result -> {
                if (result == ButtonType.OK) {
                    // Create new parameters with the current base path
                    String currentBasePath = basePathField.getText();
                    resourceParams = new MetaResourceParameters(currentBasePath);
                    
                    // Update UI
                    updateUI();
                }
            });
        });
        
        Button saveButton = new Button("Save Configuration");
        saveButton.setDefaultButton(true);
        saveButton.setOnAction(e -> {
            // Update base path from text field
            String basePath = basePathField.getText();
            if (!basePath.equals(resourceParams.getBaseResourcePath())) {
                updateBaseResourcePath(basePath);
            }
            
            // Save parameters
            boolean success = resourceIO.saveParameters(resourceParams);
            
            if (success) {
                showAlert(Alert.AlertType.INFORMATION, 
                        "Configuration Saved", 
                        "Resource configuration has been saved successfully.");
            } else {
                showAlert(Alert.AlertType.ERROR, 
                        "Save Failed", 
                        "Failed to save resource configuration.");
            }
        });
        
        panel.getChildren().addAll(resetButton, saveButton);
        
        return panel;
    }
    
    /**
     * Update the base resource path and refresh all tools
     */
    private void updateBaseResourcePath(String newPath) {
        resourceParams.setBaseResourcePath(newPath);
        updateUI();
    }
    
    /**
     * Update the UI to reflect current parameter values
     */
    private void updateUI() {
        basePathField.setText(resourceParams.getBaseResourcePath());
        
        // Force a refresh for all tool controls
        if (getScene() != null) {
            Platform.runLater(() -> {
                // This will rebuild the entire panel
                BorderPane parent = (BorderPane) getParent();
                if (parent != null) {
                    parent.setCenter(new ResourceConfigPanel());
                }
            });
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
     * Save the current configuration
     * This should be called when the window is closed
     */
    public void saveOnClose() {
        // Update base path from text field
        String basePath = basePathField.getText();
        if (!basePath.equals(resourceParams.getBaseResourcePath())) {
            updateBaseResourcePath(basePath);
        }
        
        // Save parameters
        resourceIO.saveParameters(resourceParams);
    }
    
    /**
     * Get the current resource parameters
     */
    public MetaResourceParameters getResourceParameters() {
        return resourceParams;
    }
    
    /**
     * Inner class representing a tool path control with status indicator
     */
    private class ToolPathControl extends HBox {
        private final ToolType toolType;
        private final Circle statusIndicator;
        private final TextField pathField;
        
        public ToolPathControl(ToolType toolType) {
            super(10);
            this.toolType = toolType;
            setAlignment(Pos.CENTER_LEFT);
            
            // Status indicator
            statusIndicator = new Circle(8);
            updateStatusIndicator();
            
            // Tool name label
            Label nameLabel = new Label(toolType.getDisplayName() + ":");
            nameLabel.setPrefWidth(120);
            
            // Path field
            pathField = new TextField(resourceParams.getToolPath(toolType));
            pathField.setPrefWidth(350);
            pathField.textProperty().addListener((obs, oldVal, newVal) -> {
                resourceParams.setToolPath(toolType, newVal);
                updateStatusIndicator();
            });
            
            // Browse button
            Button browseButton = new Button("Browse...");
            browseButton.setOnAction(e -> {
                // Determine if we need a file chooser or directory chooser
                if (toolType.toString().contains("DIR")) {
                    // For directory-based tools
                    DirectoryChooser dirChooser = new DirectoryChooser();
                    dirChooser.setTitle("Select " + toolType.getDisplayName() + " Directory");
                    
                    // Set initial directory if current path exists
                    File currentDir = new File(pathField.getText());
                    if (currentDir.exists()) {
                        if (currentDir.isDirectory()) {
                            dirChooser.setInitialDirectory(currentDir);
                        } else if (currentDir.getParentFile() != null) {
                            dirChooser.setInitialDirectory(currentDir.getParentFile());
                        }
                    }
                    
                    File selectedDir = dirChooser.showDialog(getScene().getWindow());
                    if (selectedDir != null) {
                        pathField.setText(selectedDir.getAbsolutePath());
                        resourceParams.setToolPath(toolType, selectedDir.getAbsolutePath());
                        updateStatusIndicator();
                    }
                } else {
                    // For executable tools
                    FileChooser fileChooser = new FileChooser();
                    fileChooser.setTitle("Select " + toolType.getDisplayName() + " Executable");
                    
                    // Add appropriate extension filters
                    MetaResourceParameters.ToolInfo toolInfo = resourceParams.getToolInfoMap().get(toolType);
                    String executableName = toolInfo != null ? toolInfo.getExecutableName() : "";
                    
                    if (executableName.endsWith(".exe")) {
                        fileChooser.getExtensionFilters().add(
                            new FileChooser.ExtensionFilter("Executable Files", "*.exe")
                        );
                    } else if (executableName.endsWith(".jar")) {
                        fileChooser.getExtensionFilters().add(
                            new FileChooser.ExtensionFilter("Java JAR Files", "*.jar")
                        );
                    } else if (executableName.endsWith(".bat")) {
                        fileChooser.getExtensionFilters().add(
                            new FileChooser.ExtensionFilter("Batch Files", "*.bat")
                        );
                    }
                    
                    fileChooser.getExtensionFilters().add(
                        new FileChooser.ExtensionFilter("All Files", "*.*")
                    );
                    
                    // Set initial directory if current path exists
                    File currentFile = new File(pathField.getText());
                    if (currentFile.exists()) {
                        if (currentFile.isDirectory()) {
                            fileChooser.setInitialDirectory(currentFile);
                        } else if (currentFile.getParentFile() != null) {
                            fileChooser.setInitialDirectory(currentFile.getParentFile());
                        }
                    }
                    
                    File selectedFile = fileChooser.showOpenDialog(getScene().getWindow());
                    if (selectedFile != null) {
                        pathField.setText(selectedFile.getAbsolutePath());
                        resourceParams.setToolPath(toolType, selectedFile.getAbsolutePath());
                        updateStatusIndicator();
                    }
                }
            });
            
            getChildren().addAll(statusIndicator, nameLabel, pathField, browseButton);
        }
        
        /**
         * Update the status indicator based on tool availability
         */
        private void updateStatusIndicator() {
            ToolStatus status = resourceParams.getToolStatus(toolType);
            
            if (status == ToolStatus.AVAILABLE) {
                statusIndicator.setFill(Color.GREEN);
                setTooltip(new Tooltip("Tool is available at the specified path"));
            } else {
                statusIndicator.setFill(Color.RED);
                setTooltip(new Tooltip("Tool is not available at the specified path"));
            }
        }
        
        /**
         * Set the tooltip for the control
         */
        private void setTooltip(Tooltip tooltip) {
            Tooltip.install(statusIndicator, tooltip);
        }
    }
} 