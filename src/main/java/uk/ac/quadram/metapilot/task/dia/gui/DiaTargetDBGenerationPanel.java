package uk.ac.quadram.metapilot.task.dia.gui;

import javafx.application.Platform;
import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.control.*;
import javafx.scene.layout.*;
import javafx.scene.text.Font;
import javafx.scene.text.FontWeight;
import javafx.scene.paint.Color;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;
import java.io.File;
import javafx.animation.KeyFrame;
import javafx.animation.Timeline;
import javafx.util.Duration;
import uk.ac.quadram.metapilot.task.dia.MetaDiaTargetDBTask;
import uk.ac.quadram.metapilot.task.dia.par.MetaParameterDiaTg;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;

public class DiaTargetDBGenerationPanel extends VBox {
    
    // Database generation components
    private TextArea consoleOutput;
    
    // Current MAG database reference
    private MagDbItem currentMagDbItem;
    
    // Internal progress tracking
    private ProgressBar dbProgressBar;
    private Label progressLabel;
    
    // Step tracking for clear progress indication
    private int currentStep = 0;
    private int totalSteps = 4;
    
    // Callback interface for parameter updates
    public interface ParameterUpdateCallback {
        void onParametersChanged();
    }
    
    private ParameterUpdateCallback parameterCallback;
    
    // Required parameters for MetaDiaTargetDBTask
    private MetaParameterDiaTg metaParameter;
    private MetaResourceParameters resourceParameters;
    private MetaPilotDiaTargetPanel parentPanel;

    public DiaTargetDBGenerationPanel() {
        initializeUI();
    }
    
    public DiaTargetDBGenerationPanel(TextArea consoleOutput) {
        this.consoleOutput = consoleOutput;
        initializeUI();
    }
    
    public DiaTargetDBGenerationPanel(MetaParameterDiaTg metaParameter, MetaResourceParameters resourceParameters) {
        this.metaParameter = metaParameter;
        this.resourceParameters = resourceParameters;
        
        // CRITICAL: Initialize MagDbItem from metaParameter immediately
        if (metaParameter != null && metaParameter.getMagDbItem() != null) {
            this.currentMagDbItem = metaParameter.getMagDbItem();
            System.out.println("DEBUG: DiaTargetDBGenerationPanel constructor - initialized MagDbItem: " + currentMagDbItem.getCatalogueID());
        } else {
            System.err.println("ERROR: DiaTargetDBGenerationPanel constructor - MetaParameterDiaTg or MagDbItem is null!");
        }
        
        initializeUI();
    }
    
    public DiaTargetDBGenerationPanel(MetaParameterDiaTg metaParameter, MetaResourceParameters resourceParameters, TextArea consoleOutput) {
        this.metaParameter = metaParameter;
        this.resourceParameters = resourceParameters;
        this.consoleOutput = consoleOutput;
        
        // CRITICAL: Initialize MagDbItem from metaParameter immediately
        if (metaParameter != null && metaParameter.getMagDbItem() != null) {
            this.currentMagDbItem = metaParameter.getMagDbItem();
            System.out.println("DEBUG: DiaTargetDBGenerationPanel constructor - initialized MagDbItem: " + currentMagDbItem.getCatalogueID());
        } else {
            System.err.println("ERROR: DiaTargetDBGenerationPanel constructor - MetaParameterDiaTg or MagDbItem is null!");
        }
        
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
        System.out.println("DEBUG: DiaTargetDBGenerationPanel.setMagDbItem called with: " + 
            (magDbItem != null ? magDbItem.getCatalogueID() : "null"));
        if (magDbItem != null) {
            // Ensure currentFile is set - use repositoryFile if currentFile is null
            if (magDbItem.getCurrentFile() == null && magDbItem.getRepositoryFile() != null) {
                magDbItem.setCurrentFile(magDbItem.getRepositoryFile());
                System.out.println("DEBUG: Set currentFile from repositoryFile: " + magDbItem.getRepositoryFile().getAbsolutePath());
            }
            System.out.println("DEBUG: MagDbItem.getCurrentFile() = " + 
                (magDbItem.getCurrentFile() != null ? magDbItem.getCurrentFile().getAbsolutePath() : "null"));
        }
        
        // Note: Status check is only done once during panel initialization
        // to avoid duplicate status checks when database is updated
    }
    
    public void setMagDbItemWithVersion(MagDbItem magDbItem, MagDbItem.MagVersion selectedVersion) {
        this.currentMagDbItem = magDbItem;
        System.out.println("DEBUG: DiaTargetDBGenerationPanel.setMagDbItemWithVersion called with: " + 
            (magDbItem != null ? magDbItem.getCatalogueID() : "null") + 
            ", version: " + (selectedVersion != null ? selectedVersion.getCatalogueVersion() : "null"));
        
        if (magDbItem != null && selectedVersion != null) {
            // Set the version on the MagDbItem - this is critical for correct paths!
            magDbItem.setUsedVersion(selectedVersion);
            System.out.println("DEBUG: Set version " + selectedVersion.getCatalogueVersion() + " on MagDbItem " + magDbItem.getCatalogueID());
            
            // Ensure currentFile is set - use repositoryFile if currentFile is null
            if (magDbItem.getCurrentFile() == null && magDbItem.getRepositoryFile() != null) {
                magDbItem.setCurrentFile(magDbItem.getRepositoryFile());
                System.out.println("DEBUG: Set currentFile from repositoryFile: " + magDbItem.getRepositoryFile().getAbsolutePath());
            }
            System.out.println("DEBUG: MagDbItem.getCurrentFile() = " + 
                (magDbItem.getCurrentFile() != null ? magDbItem.getCurrentFile().getAbsolutePath() : "null"));
        }
        
        // Note: Status check is only done once during panel initialization
        // to avoid duplicate status checks when database is updated
    }
    
    public void setMetaParameter(MetaParameterDiaTg metaParameter) {
        this.metaParameter = metaParameter;
    }
    
    public void setResourceParameters(MetaResourceParameters resourceParameters) {
        this.resourceParameters = resourceParameters;
    }
    
    public void setParentPanel(MetaPilotDiaTargetPanel parentPanel) {
        this.parentPanel = parentPanel;
    }
    
    public MagDbItem getCurrentMagDbItem() {
        return currentMagDbItem;
    }
    
    public void updateMetaParameter(MetaParameterDiaTg metaParameter) {
        this.metaParameter = metaParameter;
    }

    private void initializeUI() {
        setSpacing(10);
        setPadding(new Insets(10));
        setStyle("-fx-background-color: #fff8dc; -fx-background-radius: 5;");
        setMinHeight(400);
        VBox.setVgrow(this, Priority.ALWAYS);

        Label titleLabel = createSectionTitle("Taxa-Function Database Generation");
        
        // Generation controls - only the Generate Database button
        HBox generationBox = createGenerationControls();
        
        // Progress and console section
        VBox progressSection = createProgressSection();
        
        getChildren().addAll(titleLabel, generationBox, progressSection);
        
        // Automatically check database status when panel is created
        Platform.runLater(() -> checkDatabaseStatusAutomatically());
    }
    
    private HBox createGenerationControls() {
        HBox generationBox = new HBox(10);
        generationBox.setAlignment(Pos.CENTER_LEFT);
        
        Button generateDbButton = new Button("Generate Database");
        generateDbButton.setPrefWidth(150);
        generateDbButton.setStyle("-fx-background-color: #4CAF50; -fx-text-fill: white; -fx-font-weight: bold;");
        generateDbButton.setOnAction(e -> generateDatabase());
        
        generationBox.getChildren().addAll(generateDbButton);
        return generationBox;
    }
    
    private VBox createProgressSection() {
        VBox progressSection = new VBox(8);
        progressSection.setPadding(new Insets(10));
        progressSection.setStyle("-fx-background-color: #f8f8f8; -fx-background-radius: 5;");
        progressSection.setMinHeight(200);
        progressSection.setPrefHeight(400);
        VBox.setVgrow(progressSection, Priority.ALWAYS);

        dbProgressBar = new ProgressBar(0);
        dbProgressBar.setMaxWidth(Double.MAX_VALUE);
        dbProgressBar.setPrefHeight(25);
        dbProgressBar.setStyle("-fx-accent: #4CAF50;");

        // Progress label to show current step
        progressLabel = new Label("Ready to generate database");
        progressLabel.setStyle("-fx-font-weight: bold; -fx-font-size: 12px; -fx-text-fill: #333333;");
        progressLabel.setMaxWidth(Double.MAX_VALUE);

        Label consoleTitle = new Label("Console Output:");
        consoleTitle.setStyle("-fx-font-weight: bold; -fx-font-size: 14px;");

        // Create internal console if external one not provided
        if (consoleOutput == null) {
            consoleOutput = new TextArea();
        }
        consoleOutput.setEditable(false);
        consoleOutput.setPrefHeight(300);
        consoleOutput.setMinHeight(200);
        consoleOutput.setStyle("-fx-font-family: 'Courier New', monospace; -fx-font-size: 12px;");
        VBox.setVgrow(consoleOutput, Priority.ALWAYS);

        progressSection.getChildren().addAll(
            progressLabel, dbProgressBar, consoleTitle, consoleOutput
        );
        
        return progressSection;
    }
    
    private void generateDatabase() {
        if (resourceParameters == null) {
            showError("Configuration Error", "MetaResourceParameters is not set. Please configure the panel properly.");
            return;
        }
        
        if (currentMagDbItem == null) {
            if (consoleOutput != null) {
                consoleOutput.appendText("ERROR: No MAG database selected\n");
                consoleOutput.appendText("Please select a MAG database from the Database Settings panel above.\n");
                consoleOutput.appendText("The database generation panel needs a valid MAG database to work.\n");
            }
            showError("Database Error", "MAG database item is not set. Please select a database first from the Database Settings panel above.");
            return;
        }
        
        // Create a simple MetaParameterDiaTg with default values for database generation
        // This is for internal use, so we don't need complex parameters
        try {
            metaParameter = new MetaParameterDiaTg(
                new uk.ac.quadram.metapilot.task.par.MetaBaseParameter(
                    uk.ac.quadram.metapilot.task.MetaLabWorkflowType.DiaNNMAGTG.name(),
                    new uk.ac.quadram.metapilot.task.par.MetaData(new String[0], new String[0], new int[0], new int[0]), // Empty metadata
                    "", // Empty result directory
                    currentMagDbItem.getCurrentFile() != null ? currentMagDbItem.getCurrentFile().getAbsolutePath() : "",
                    false, // appendHostDb - not used in DIA Target
                    "", // hostDB - not used in DIA Target
                    "DIA", // MS2 scan mode
                    1, // Core count
                    4, // Thread count
                    true // Meta workflow
                ),
                currentMagDbItem,
                false, false, new String[0], false, new String[0]
            );
            System.out.println("DEBUG: Created default MetaParameterDiaTg for database generation");
        } catch (Exception e) {
            System.err.println("ERROR: Failed to create MetaParameterDiaTg: " + e.getMessage());
            showError("Configuration Error", "Failed to create database generation parameters: " + e.getMessage());
            return;
        }
        
        // Create and run the MetaDiaTargetDBTask
        MetaDiaTargetDBTask dbTask = new MetaDiaTargetDBTask(
            metaParameter, 
            resourceParameters, 
            dbProgressBar
        );
        
        // Set up task event handlers
        dbTask.setOnRunning(e -> {
            Platform.runLater(() -> {
                // Initialize progress tracking
                currentStep = 0;
                updateProgressStep(0, "Starting database generation...");
                
                if (consoleOutput != null) {
                    consoleOutput.appendText("=== Taxa-Function Database Generation Started ===\n");
                    consoleOutput.appendText("Database name: pro_func\n");
                    consoleOutput.appendText("MAG database: " + currentMagDbItem.getCatalogueID() + "\n");
                    consoleOutput.appendText("Progress: Step 1/4 - Initializing...\n");
                }
            });
        });
        
        dbTask.setOnSucceeded(e -> {
            Platform.runLater(() -> {
                // Complete progress tracking
                updateProgressStep(totalSteps, "Database generation completed successfully!");
                dbProgressBar.setProgress(1.0);
                
                if (consoleOutput != null) {
                    consoleOutput.appendText("=== Database Generation Completed ===\n");
                    consoleOutput.appendText("✓ Database 'pro_func' generated successfully!\n");
                    consoleOutput.appendText("✓ Database has been added to the MAG database.\n");
                    consoleOutput.appendText("✓ All steps completed (4/4)\n");
                }
                
                notifyParameterUpdate();
                
                // Reset progress after a delay
                Timeline resetTimeline = new Timeline(new KeyFrame(Duration.seconds(3), ev -> {
                    dbProgressBar.setProgress(0);
                }));
                resetTimeline.play();
            });
        });
        
        dbTask.setOnFailed(e -> {
            Platform.runLater(() -> {
                // Reset progress tracking on failure
                updateProgressStep(0, "Database generation failed!");
                dbProgressBar.setProgress(0);
                
                Throwable exception = dbTask.getException();
                String errorMessage = exception != null ? exception.getMessage() : "Unknown error occurred";
                
                if (consoleOutput != null) {
                    consoleOutput.appendText("=== Database Generation Failed ===\n");
                    consoleOutput.appendText("✗ Error at step " + (currentStep + 1) + "/" + totalSteps + "\n");
                    consoleOutput.appendText("✗ Error message: " + errorMessage + "\n");
                    consoleOutput.appendText("✗ Please check the error details and try again.\n");
                }
                
                showError("Database Generation Failed", 
                    "Failed to generate database at step " + (currentStep + 1) + "/" + totalSteps + ": " + errorMessage);
            });
        });
        
        dbTask.setOnCancelled(e -> {
            Platform.runLater(() -> {
                // Reset progress tracking on cancellation
                updateProgressStep(0, "Database generation was cancelled");
                dbProgressBar.setProgress(0);
                
                if (consoleOutput != null) {
                    consoleOutput.appendText("=== Database Generation Cancelled ===\n");
                    consoleOutput.appendText("⚠ Database generation was cancelled by user.\n");
                    consoleOutput.appendText("⚠ Progress: " + currentStep + "/" + totalSteps + " steps completed.\n");
                }
            });
        });
        
        // Set up progress monitoring from the task
        dbTask.messageProperty().addListener((obs, oldMsg, newMsg) -> {
            if (newMsg != null && !newMsg.isEmpty()) {
                Platform.runLater(() -> {
                    // Parse progress messages and update step tracking
                    if (newMsg.contains("Initializing") || newMsg.contains("Starting")) {
                        updateProgressStep(1, "Step 1/4: " + newMsg);
                    } else if (newMsg.contains("Processing") || newMsg.contains("Reading")) {
                        updateProgressStep(2, "Step 2/4: " + newMsg);
                    } else if (newMsg.contains("Building") || newMsg.contains("Creating")) {
                        updateProgressStep(3, "Step 3/4: " + newMsg);
                    } else if (newMsg.contains("Finalizing") || newMsg.contains("Completing")) {
                        updateProgressStep(4, "Step 4/4: " + newMsg);
                    }
                    
                    if (consoleOutput != null) {
                        consoleOutput.appendText("Progress: " + newMsg + "\n");
                    }
                });
            }
        });
        
        // Start the task
        new Thread(dbTask).start();
    }
    
    /**
     * Automatically check database status when panel opens
     */
    private void checkDatabaseStatusAutomatically() {
        if (currentMagDbItem == null) {
            if (consoleOutput != null) {
                consoleOutput.appendText("=== Database Status Check ===\n");
                consoleOutput.appendText("INFO: No MAG database selected\n");
                consoleOutput.appendText("Please select a MAG database from the Database Settings panel above.\n");
                consoleOutput.appendText("The database generation panel needs a valid MAG database to work.\n");
            }
            return;
        }
        
        // Get the MAG database folder
        File magDbFolder = currentMagDbItem.getCurrentFile();
        if (magDbFolder == null || !magDbFolder.exists()) {
            if (consoleOutput != null) {
                consoleOutput.appendText("=== Database Status Check ===\n");
                consoleOutput.appendText("ERROR: MAG database folder not found: " + 
                    (magDbFolder != null ? magDbFolder.getAbsolutePath() : "null") + "\n");
                consoleOutput.appendText("Please check the database installation.\n");
            }
            return;
        }
        
        // Check for pro_func.db
        File proFuncDb = new File(magDbFolder, "pro_func.db");
        File eggnogDir = new File(magDbFolder, "eggNOG");
        File taxaFile = new File(magDbFolder, "genomes-all_metadata.tsv");
        
        if (consoleOutput != null) {
            consoleOutput.appendText("=== Database Status Check ===\n");
            consoleOutput.appendText("MAG Database: " + currentMagDbItem.getCatalogueID() + "\n");
            consoleOutput.appendText("Database Path: " + magDbFolder.getAbsolutePath() + "\n");
            consoleOutput.appendText("pro_func.db: " + (proFuncDb.exists() ? "✓ Found" : "✗ Missing") + "\n");
            consoleOutput.appendText("eggNOG directory: " + (eggnogDir.exists() ? "✓ Found" : "✗ Missing") + "\n");
            consoleOutput.appendText("genomes-all_metadata.tsv: " + (taxaFile.exists() ? "✓ Found" : "✗ Missing") + "\n");
            
            if (proFuncDb.exists()) {
                consoleOutput.appendText("STATUS: ✓ Taxa-function database is available\n");
                consoleOutput.appendText("ACTION: Ready for analysis\n");
            } else {
                consoleOutput.appendText("STATUS: ⚠ Taxa-function database not found\n");
                if (eggnogDir.exists() && taxaFile.exists()) {
                    consoleOutput.appendText("ACTION: Click 'Generate Database' to create pro_func.db\n");
                } else {
                    consoleOutput.appendText("ACTION: Required source files missing - cannot generate database\n");
                }
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
    
    /**
     * Update progress step with clear indication
     */
    private void updateProgressStep(int step, String message) {
        currentStep = step;
        double progress = (double) step / totalSteps;
        dbProgressBar.setProgress(progress);
        
        if (progressLabel != null) {
            if (step == 0) {
                progressLabel.setText("Ready to generate database");
                progressLabel.setStyle("-fx-font-weight: bold; -fx-font-size: 12px; -fx-text-fill: #333333;");
            } else if (step <= totalSteps) {
                progressLabel.setText("Step " + step + "/" + totalSteps + ": " + message);
                progressLabel.setStyle("-fx-font-weight: bold; -fx-font-size: 12px; -fx-text-fill: #0066CC;");
            } else {
                progressLabel.setText("✓ " + message);
                progressLabel.setStyle("-fx-font-weight: bold; -fx-font-size: 12px; -fx-text-fill: #009900;");
            }
        }
        
        System.out.println("DEBUG: Progress updated - Step " + step + "/" + totalSteps + ": " + message);
    }
}
