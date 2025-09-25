package uk.ac.quadram.metapilot.task.common.gui;

import javafx.geometry.Insets;
import javafx.scene.control.Alert;
import javafx.scene.control.Alert.AlertType;
import javafx.scene.control.Tab;
import javafx.scene.control.TabPane;
import javafx.scene.layout.Priority;
import javafx.scene.layout.VBox;
import javafx.application.Platform;
import java.io.File;
import uk.ac.quadram.metapilot.task.par.MetaBaseParameter;
import uk.ac.quadram.metapilot.task.par.BaseParameterIO;

/**
 * Abstract base class for workflow panels that provides common functionality
 * across different workflow types (MAG, IGC, etc.).
 * 
 * @param <T> The type of parameters used by the workflow
 */
public abstract class AbstractWorkflowPanel<T extends MetaBaseParameter> extends VBox {
    
    // Common fields
    protected ConsolePanel consolePanel;
    protected BaseParameterIO parameterIO;
    
    // Current task reference for common operations
    protected Object currentTask;
    
    // Current analysis parameters for viewReport
    protected T currentAnalysisParams;
    
    /**
     * Constructor - subclasses should call this and then initialize their specific components
     */
    protected AbstractWorkflowPanel() {
        // Initialize common components
        initializeCommonComponents();
    }
    
    /**
     * Initialize common UI components that are shared across all workflow panels
     */
    protected void initializeCommonComponents() {
        setSpacing(10);
        setPadding(new Insets(10));
        setMinWidth(1200);
        setMinHeight(800);
        setPrefWidth(1600);
        setPrefHeight(900);
        
        // Initialize console panel
        consolePanel = new ConsolePanel();
        consolePanel.setStartCallback(this::startAnalysis);
        consolePanel.setStopCallback(this::stopAnalysis);
        consolePanel.setViewReportCallback(this::viewReport);
    }
    
    /**
     * Create the main tab pane structure that is common to all workflow panels
     */
    protected TabPane createMainTabPane() {
        TabPane tabPane = new TabPane();
        tabPane.setTabClosingPolicy(TabPane.TabClosingPolicy.UNAVAILABLE);
        tabPane.setMinWidth(980);
        tabPane.setMinHeight(780);
        VBox.setVgrow(tabPane, Priority.ALWAYS);
        
        // Create tabs
        Tab ioTab = new Tab("Input/Output");
        ioTab.setContent(createIOTab());
        
        Tab paramTab = new Tab("Parameters");
        paramTab.setContent(createParametersTab());
        
        Tab consoleTab = new Tab("Console");
        consoleTab.setContent(consolePanel);
        
        tabPane.getTabs().addAll(ioTab, paramTab, consoleTab);
        
        // Add tab change listener for common behavior
        setupTabChangeListener(tabPane, paramTab, consoleTab);
        
        return tabPane;
    }
    
    /**
     * Setup common tab change listener behavior
     */
    protected void setupTabChangeListener(TabPane tabPane, Tab paramTab, Tab consoleTab) {
        tabPane.getSelectionModel().selectedItemProperty().addListener((obs, oldTab, newTab) -> {
            if (newTab == paramTab) {
                onParametersTabSelected();
            } else if (newTab == consoleTab) {
                onConsoleTabSelected();
            }
        });
    }
    
    /**
     * Called when Parameters tab is selected - subclasses can override
     */
    protected void onParametersTabSelected() {
        // Default implementation - subclasses can override
    }
    
    /**
     * Called when Console tab is selected - common auto-save behavior
     */
    protected void onConsoleTabSelected() {
        try {
            consolePanel.appendConsoleText("Auto-saving parameters when switching to Console tab...\n");
            T currentParams = collectParametersFromUI();
            if (currentParams != null) {
                boolean saved = parameterIO.saveLastUsedParameters(currentParams);
                if (saved) {
                    consolePanel.appendConsoleText("✓ Parameters auto-saved to config folder\n");
                } else {
                    consolePanel.appendConsoleText("WARNING: Failed to auto-save parameters to config folder\n");
                }
                
                // Also save preferences
                boolean prefsSaved = parameterIO.savePreferences(currentParams);
                if (prefsSaved) {
                    consolePanel.appendConsoleText("✓ Preferences auto-saved\n");
                } else {
                    consolePanel.appendConsoleText("WARNING: Failed to auto-save preferences\n");
                }
            } else {
                consolePanel.appendConsoleText("WARNING: No parameters to auto-save\n");
            }
        } catch (Exception e) {
            String errorMsg = "ERROR: Failed to auto-save parameters: " + e.getMessage();
            consolePanel.appendConsoleText(errorMsg + "\n");
        }
    }
    
    /**
     * View report method - implemented exactly like MetaPilotDdaMagPfindPanel
     */
    protected void viewReport() {
        try {
            consolePanel.appendConsoleText("Opening report in default browser...\n");
            
            if (currentAnalysisParams == null) {
                consolePanel.appendConsoleText("ERROR: No analysis parameters available. Please run analysis first.\n");
                return;
            }
            
            File reportDir = currentAnalysisParams.getReportDir();
            if (reportDir == null || !reportDir.exists()) {
                consolePanel.appendConsoleText("ERROR: Report directory not available or does not exist.\n");
                consolePanel.appendConsoleText("Please ensure the analysis has completed successfully.\n");
                return;
            }
            
            File indexFile = new File(reportDir, "index.html");
            
            if (!indexFile.exists()) {
                consolePanel.appendConsoleText("ERROR: Main report index.html not found in " + reportDir.getAbsolutePath() + "\n");
                consolePanel.appendConsoleText("Please ensure the analysis has completed successfully.\n");
                return;
            }
            
            consolePanel.appendConsoleText("Opening main report: index.html\n");
            
            if (java.awt.Desktop.isDesktopSupported()) {
                java.awt.Desktop.getDesktop().browse(indexFile.toURI());
                consolePanel.appendConsoleText("✓ Report opened in default browser\n");
            } else {
                consolePanel.appendConsoleText("ERROR: Desktop not supported. Please open manually: " + indexFile.getAbsolutePath() + "\n");
            }
            
        } catch (Exception e) {
            String errorMsg = "ERROR: Failed to open report: " + e.getMessage();
            consolePanel.appendConsoleText(errorMsg + "\n");
            System.err.println("Failed to open report: " + e.getMessage());
            e.printStackTrace();
        }
    }
    
    /**
     * Show alert dialog - common utility method
     */
    protected void showAlert(String title, String message) {
        Platform.runLater(() -> {
            Alert alert = new Alert(AlertType.INFORMATION);
            alert.setTitle(title);
            alert.setHeaderText(null);
            alert.setContentText(message);
            alert.showAndWait();
        });
    }
    
    /**
     * Show error dialog - common utility method
     */
    protected void showError(String title, String content) {
        Platform.runLater(() -> {
            Alert alert = new Alert(AlertType.ERROR);
            alert.setTitle(title);
            alert.setHeaderText(null);
            alert.setContentText(content);
            alert.showAndWait();
        });
    }
    
    /**
     * Setup window close handler - common behavior
     */
    protected void setupWindowCloseHandler() {
        Platform.runLater(() -> {
            try {
                // Get the current scene and window
                if (getScene() != null) {
                    javafx.stage.Window window = getScene().getWindow();
                    if (window != null) {
                        window.setOnCloseRequest(event -> {
                            try {
                                consolePanel.appendConsoleText("Saving parameters on window close...\n");
                                T currentParams = collectParametersFromUI();
                                if (currentParams != null) {
                                    boolean saved = parameterIO.saveLastUsedParameters(currentParams);
                                    if (saved) {
                                        consolePanel.appendConsoleText("✓ Parameters saved to config folder\n");
                                    } else {
                                        consolePanel.appendConsoleText("WARNING: Failed to save parameters to config folder\n");
                                    }
                                    
                                    // Also save preferences
                                    boolean prefsSaved = parameterIO.savePreferences(currentParams);
                                    if (prefsSaved) {
                                        consolePanel.appendConsoleText("✓ Preferences saved\n");
                                    } else {
                                        consolePanel.appendConsoleText("WARNING: Failed to save preferences\n");
                                    }
                                } else {
                                    consolePanel.appendConsoleText("WARNING: No parameters to save on window close\n");
                                }
                            } catch (Exception e) {
                                String errorMsg = "ERROR: Failed to save parameters on window close: " + e.getMessage();
                                consolePanel.appendConsoleText(errorMsg + "\n");
                            }
                        });
                    }
                }
            } catch (Exception e) {
                // Handle exception silently to avoid disrupting UI
            }
        });
    }
    
    // Abstract methods that subclasses must implement
    
    /**
     * Create the Input/Output tab content
     */
    protected abstract VBox createIOTab();
    
    /**
     * Create the Parameters tab content
     */
    protected abstract VBox createParametersTab();
    
    /**
     * Collect all parameters from the UI into a parameters object
     */
    protected abstract T collectParametersFromUI();
    
    /**
     * Apply parameters to the UI components
     */
    protected abstract void applyParametersToUI(T params);
    
    /**
     * Start the analysis - implementation specific to each workflow
     */
    protected abstract void startAnalysis();
    
    /**
     * Stop the analysis - implementation specific to each workflow
     */
    protected abstract void stopAnalysis();
    
    /**
     * Load last used parameters - implementation specific to each workflow
     */
    protected abstract void loadLastUsedParameters();
    
    /**
     * Update parameter view - implementation specific to each workflow
     */
    protected abstract void updateParameterView();
    
    /**
     * Get the workflow name for logging purposes
     */
    protected abstract String getWorkflowName();
    
    // Getters for common components
    
    /**
     * Get the console panel
     */
    public ConsolePanel getConsolePanel() {
        return consolePanel;
    }
    
    /**
     * Get the current analysis parameters
     */
    public T getCurrentAnalysisParams() {
        return currentAnalysisParams;
    }
    
    /**
     * Set the current analysis parameters
     */
    public void setCurrentAnalysisParams(T params) {
        this.currentAnalysisParams = params;
    }
    
    /**
     * Get the current task
     */
    public Object getCurrentTask() {
        return currentTask;
    }
    
    /**
     * Set the current task
     */
    public void setCurrentTask(Object task) {
        this.currentTask = task;
    }
}
