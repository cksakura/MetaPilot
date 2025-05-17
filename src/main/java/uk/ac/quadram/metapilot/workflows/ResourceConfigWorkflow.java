package uk.ac.quadram.metapilot.workflows;

import javafx.scene.Node;
import javafx.scene.layout.BorderPane;
import javafx.stage.Stage;
import javafx.stage.WindowEvent;

import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.resources.MetaResourceIO;
import uk.ac.quadram.metapilot.ui.ResourceConfigPanel;

/**
 * Workflow for configuring resources
 */
public class ResourceConfigWorkflow extends BaseWorkflow {
    
    private ResourceConfigPanel configPanel;
    private MetaResourceIO resourceIO;
    
    public ResourceConfigWorkflow() {
        super("Resources Configuration", "Configure paths to external tools");
        resourceIO = new MetaResourceIO();
    }
    
    @Override
    public void init() {
        // Nothing to initialize
    }
    
    @Override
    public Node createWorkflowUI() {
        BorderPane panel = new BorderPane();
        
        // Create configuration panel
        configPanel = new ResourceConfigPanel();
        panel.setCenter(configPanel);
        
        // Add window close event handler to save parameters
        if (getStage() != null) {
            getStage().setOnCloseRequest(this::handleClose);
        }
        
        return panel;
    }
    
    /**
     * Handle window close event
     */
    private void handleClose(WindowEvent event) {
        // Save parameters when the window is closed
        if (configPanel != null) {
            configPanel.saveOnClose();
        }
    }
    
    @Override
    public void start() {
        // Nothing to start
    }
    
    @Override
    public void showStage() {
        Stage stage = getStage();
        if (stage != null) {
            // Set window size
            stage.setWidth(800);
            stage.setHeight(600);
            
            // Show the stage
            stage.show();
            
            // Add close event handler
            stage.setOnCloseRequest(this::handleClose);
        }
    }
    
    @Override
    public void closeStage() {
        Stage stage = getStage();
        if (stage != null) {
            // Save parameters before closing
            if (configPanel != null) {
                configPanel.saveOnClose();
            }
            
            // Close the stage
            stage.close();
        }
    }
    
    /**
     * Get the resource parameters
     */
    public MetaResourceParameters getResourceParameters() {
        return configPanel != null ? configPanel.getResourceParameters() : null;
    }
} 