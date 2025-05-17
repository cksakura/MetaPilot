package uk.ac.quadram.metapilot.workflows;

import javafx.scene.Node;
import javafx.scene.Parent;
import javafx.scene.Scene;
import javafx.scene.layout.BorderPane;
import javafx.stage.Stage;

/**
 * Base class for workflow implementations
 */
public abstract class BaseWorkflow {
    
    private String name;
    private String description;
    private Stage stage;
    
    /**
     * Create a new workflow with the given name and description
     */
    public BaseWorkflow(String name, String description) {
        this.name = name;
        this.description = description;
        
        // Create a new stage for this workflow
        this.stage = new Stage();
        this.stage.setTitle(name);
    }
    
    /**
     * Get the name of this workflow
     */
    public String getName() {
        return name;
    }
    
    /**
     * Get the description of this workflow
     */
    public String getDescription() {
        return description;
    }
    
    /**
     * Get the stage for this workflow
     */
    public Stage getStage() {
        return stage;
    }
    
    /**
     * Initialize the workflow
     */
    public abstract void init();
    
    /**
     * Create the UI for this workflow
     */
    public abstract Node createWorkflowUI();
    
    /**
     * Start the workflow
     */
    public abstract void start();
    
    /**
     * Show the workflow stage
     */
    public void showStage() {
        if (stage != null) {
            // Create the UI and scene if not already created
            if (stage.getScene() == null) {
                Node ui = createWorkflowUI();
                if (ui != null) {
                    // Ensure UI is a Parent object
                    Parent root;
                    if (ui instanceof Parent) {
                        root = (Parent) ui;
                    } else {
                        // Wrap non-Parent nodes in a BorderPane
                        BorderPane wrapper = new BorderPane();
                        wrapper.setCenter(ui);
                        root = wrapper;
                    }
                    
                    Scene scene = new Scene(root);
                    stage.setScene(scene);
                }
            }
            
            // Show the stage
            stage.show();
        }
    }
    
    /**
     * Close the workflow stage
     */
    public void closeStage() {
        if (stage != null) {
            stage.close();
        }
    }
} 