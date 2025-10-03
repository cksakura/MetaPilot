package uk.ac.quadram.metapilot.task.dia.gui;

import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.control.*;
import javafx.scene.layout.*;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.scene.text.Font;
import javafx.scene.text.FontWeight;
import javafx.scene.paint.Color;
import javafx.application.Platform;
import uk.ac.quadram.metapilot.task.common.gui.MetaPilotIOPanelFX;
import uk.ac.quadram.metapilot.task.common.gui.MagDatabaseSettingsPanel;
import uk.ac.quadram.metapilot.task.common.gui.MetadataConfigurationPanel;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.resources.MetaResourceIO;
import uk.ac.quadram.metapilot.task.dda.igc.gui.MetaPilotDdaIgcMaxQuantPanel;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;
import uk.ac.quadram.metapilot.task.mag.MagDbItem.MagVersion;
import uk.ac.quadram.metapilot.task.par.MetaBaseParameter;
import uk.ac.quadram.metapilot.task.par.MetaData;
import uk.ac.quadram.metapilot.task.dia.par.DiaParameterIO;
import uk.ac.quadram.metapilot.task.dia.par.MetaParameterDia;
import uk.ac.quadram.metapilot.task.dia.MetaDiaTask;
import uk.ac.quadram.metapilot.dbSearch.deepDetect.DeepDetectTask;
import uk.ac.quadram.metapilot.dbSearch.deepDetect.DeepDetectParameter;
import uk.ac.quadram.metapilot.task.common.gui.AbstractWorkflowPanel;
import uk.ac.quadram.metapilot.task.common.gui.ConsolePanel;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.File;
import java.time.LocalDateTime;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

public class MetaPilotDiaNormalPanel extends AbstractWorkflowPanel<MetaParameterDia> {
    private static final Logger LOGGER = LogManager.getLogger(MetaPilotDiaNormalPanel.class);
    
    // UI Components - Database specific for MAG (like DIA workflows)
    private MetaPilotIOPanelFX ioPanel;
    private MagDatabaseSettingsPanel magDatabasePanel;
    private VBox magDbStatusPanel;
    private TextArea predictionConsole;
    private ProgressBar predictionProgressBar;
    private Button predictButton;
    private TextField timeTextField;
    private Button updateButton;
    
    // DeepDetect task related fields
    private DeepDetectTask deepDetectTask;
    private ExecutorService taskExecutor;
    private MetaResourceParameters resourceParams;
    
    // MetaDiaTask for DIA analysis
    private MetaDiaTask currentTask;
    
    

    // DIA-specific mode selection
    private ComboBox<String> modeSelectionComboBox;
    private ComboBox<Integer> threadsComboBox;

    // Separate panels for library and model generation
    private DiaLibraryGenerationPanel libraryGenerationPanel;
    private DiaModelGenerationPanel modelGenerationPanel;
    
    // UI references for dynamic panel replacement
    private VBox leftColumn;
    private VBox libraryPlaceholder;

    // Metadata configuration components
    private MetadataConfigurationPanel metadataConfigurationPanel;


    // MAG database items (now handled by MagDatabaseSettingsPanel)
    
    // Add parameter IO instance

    public MetaPilotDiaNormalPanel(MagDbItem[] magDbItems, MetaResourceParameters resourceParams) {
        // Call parent constructor first
        super();
        
        // Initialize parameter IO instance
        parameterIO = new DiaParameterIO();
        this.resourceParams = resourceParams;
        initializeSpecificUI();
        // Set MAG databases after UI is initialized
        setAvailableMagDatabases(magDbItems);
        // Load last used parameters after UI and databases are initialized
        loadLastUsedParameters();
    }

    // Load available MAG databases
    public void setAvailableMagDatabases(MagDbItem[] magDbItems) {
        if (magDatabasePanel != null) {
            magDatabasePanel.setAvailableMagDatabases(magDbItems);
            // Update the panels with the selected MAG database
            if (magDbItems != null && magDbItems.length > 0) {
                updatePanelsWithMagDb(magDbItems[0]);
            }
            // Now that database panel is initialized, check the status
            updateMagDbStatus();
            
            // Create library generation panel now that MAG database is available
            createLibraryGenerationPanel();
        } else {
            System.out.println("DEBUG: magDatabasePanel is null - UI not initialized yet");
        }
    }

    private void updatePanelsWithMagDb(MagDbItem magDbItem) {
        if (libraryGenerationPanel != null) {
            libraryGenerationPanel.setMagDbItem(magDbItem);
        }
        if (modelGenerationPanel != null) {
            modelGenerationPanel.setMagDbItem(magDbItem);
        }
    }
    
    /**
     * Create the library generation panel after MAG database is available
     */
    private void createLibraryGenerationPanel() {
        if (libraryGenerationPanel != null) {
            // Already created, just update it
            return;
        }
        
        System.out.println("DEBUG: Creating library generation panel with MAG database available");
        
        // Try to load last used parameters first, then create with minimal parameters for library generation
        loadLastUsedParameters();
        MetaParameterDia initialMetaParameter = null;
        // Use minimal validation for library generation (doesn't require raw files)
        initialMetaParameter = generateMetaParameterDiaForLibraryGeneration();
        
        if (initialMetaParameter != null) {
            libraryGenerationPanel = new DiaLibraryGenerationPanel(initialMetaParameter, resourceParams, new File("MAGDB_1_0.json"));
            libraryGenerationPanel.setParameterUpdateCallback(() -> updateParameterView());
            libraryGenerationPanel.setParentPanel(this); // Set reference to this panel
            
            // Add the panel to the left column in the parameters tab
            // We need to find the left column and add the panel
            addLibraryPanelToUI();
            
            System.out.println("DEBUG: Library generation panel created successfully");
        } else {
            System.err.println("ERROR: Failed to create library generation panel - no valid parameters");
        }
    }
    
    /**
     * Ensure the library generation panel is created when switching to Parameters tab
     */
    private void ensureLibraryGenerationPanelCreated() {
        System.out.println("DEBUG: ensureLibraryGenerationPanelCreated() called");
        
        // Check if library generation panel already exists
        if (libraryGenerationPanel != null) {
            System.out.println("DEBUG: Library generation panel already exists");
            return;
        }
        
        // Check if MAG database is available
        if (magDatabasePanel == null || magDatabasePanel.getSelectedMagDatabase() == null) {
            System.out.println("DEBUG: No MAG database selected yet, cannot create library generation panel");
            return;
        }
        
        System.out.println("DEBUG: Creating library generation panel on tab switch");
        createLibraryGenerationPanel();
    }
    
    /**
     * Add the library generation panel to the UI after it's created
     */
    private void addLibraryPanelToUI() {
        if (libraryGenerationPanel == null || leftColumn == null || libraryPlaceholder == null) {
            return;
        }
        
        System.out.println("DEBUG: Replacing library placeholder with actual library generation panel");
        
        // Replace the placeholder with the actual library generation panel
        int placeholderIndex = leftColumn.getChildren().indexOf(libraryPlaceholder);
        if (placeholderIndex >= 0) {
            leftColumn.getChildren().set(placeholderIndex, libraryGenerationPanel);
            System.out.println("DEBUG: Successfully replaced library placeholder with actual panel");
        } else {
            System.err.println("ERROR: Could not find library placeholder in left column");
        }
    }

    private void initializeSpecificUI() {
        setSpacing(10);
        setPadding(new Insets(10));
        setMinWidth(1200);
        setMinHeight(900);
        setPrefWidth(1600);
        setPrefHeight(1000);
        setMaxHeight(Double.MAX_VALUE);

        // Initialize reusable components
        consolePanel = new ConsolePanel();
        magDatabasePanel = new MagDatabaseSettingsPanel(false);
        magDatabasePanel.setParameterUpdateCallback(() -> {
            updateParameterView();
            // Update status when database or version selection changes
            updateMagDbStatus();
        });
        
        // Set up console panel callbacks
        consolePanel.setStartCallback(this::startAnalysis);
        consolePanel.setStopCallback(this::stopAnalysis);
        consolePanel.setViewReportCallback(this::viewReport);
        
        // Initialize database status panel
        magDbStatusPanel = createMagDbStatusPanel();

        // Create tab pane with responsive sizing
        TabPane tabPane = new TabPane();
        tabPane.setTabClosingPolicy(TabPane.TabClosingPolicy.UNAVAILABLE);
        tabPane.setMinWidth(1080);
        tabPane.setMinHeight(900);
        VBox.setVgrow(tabPane, Priority.ALWAYS);

        // Create tabs
        Tab ioTab = new Tab("Input/Output");
        ioTab.setContent(createIOTab());

        Tab paramTab = new Tab("Parameters");
        paramTab.setContent(createParametersTab());

        Tab consoleTab = new Tab("Console");
        consoleTab.setContent(consolePanel);

        tabPane.getTabs().addAll(ioTab, paramTab, consoleTab);

        // Add tab change listener to update metadata when Parameters tab is selected
        // and save parameters when Console tab is selected
        tabPane.getSelectionModel().selectedItemProperty().addListener((obs, oldTab, newTab) -> {
            if (newTab == paramTab) {
                updateMetadataTableFromIOPanel();
                // Ensure library generation panel is created when switching to Parameters tab
                ensureLibraryGenerationPanelCreated();
                // Update parameter view when switching to Parameters tab
                updateParameterView();
            } else if (newTab == consoleTab) {
                // Update parameter view when switching to Console tab
                updateParameterView();
                
                // Save current parameters when Console tab is selected
                try {
                    consolePanel.appendConsoleText("Auto-saving parameters when switching to Console tab...\n");
                    MetaParameterDia currentParams = collectParametersFromUI();
                    if (currentParams == null) {
                        // If full validation fails (e.g., no raw files), use minimal validation for library generation
                        currentParams = generateMetaParameterDiaForLibraryGeneration();
                    }
                    if (currentParams != null) {
                        boolean saved = parameterIO.saveLastUsedParameters(currentParams);
                        if (saved) {
                            consolePanel.appendConsoleText("✓ Parameters auto-saved to config folder\n");
                        } else {
                            consolePanel.appendConsoleText("WARNING: Failed to auto-save parameters to config folder\n");
                        }
                    } else {
                        consolePanel.appendConsoleText("WARNING: No parameters to auto-save\n");
                    }
                } catch (Exception e) {
                    String errorMsg = "ERROR: Failed to auto-save parameters: " + e.getMessage();
                    consolePanel.appendConsoleText(errorMsg + "\n");
                }
            }
        });

        // Add tab pane to main panel
        getChildren().add(tabPane);
        
        // Initialize parameter view
        updateParameterView();
    }

    // Create IO tab (same as MAG workflows)
    @Override
    protected String getWorkflowName() {
        return "DIA Normal";
    }
    
    @Override
    public VBox createIOTab() {
        VBox tabContent = new VBox(10);
        tabContent.setPadding(new Insets(10));
        ioPanel = new MetaPilotIOPanelFX();
        tabContent.getChildren().addAll(
            createSectionTitle("Input/Output Settings"),
            ioPanel,
            createSectionTitle("Database Settings"),
            magDatabasePanel,
            createSectionTitle("MAG Library Status"),
            magDbStatusPanel
        );
        return tabContent;
    }

    // Create MAG Database Status Panel (Simplified for peptide prediction)
    private VBox createMagDbStatusPanel() {
        VBox statusPanel = new VBox(10);
        statusPanel.setPadding(new Insets(10));
        statusPanel.setStyle("-fx-background-color: #f8f8f8; -fx-border-color: #cccccc; -fx-border-radius: 5; -fx-border-width: 1;");
        statusPanel.setMinHeight(200);
        statusPanel.setPrefHeight(300);
        
        
        // Console output panel for prediction process
        VBox consolePanel = new VBox(5);
        consolePanel.setPadding(new Insets(5));
        consolePanel.setStyle("-fx-background-color: #f0f0f0; -fx-background-radius: 3; -fx-border-color: #cccccc; -fx-border-radius: 3; -fx-border-width: 1;");
        
        Label consoleTitle = new Label("MAG Database Status Console");
        consoleTitle.setStyle("-fx-font-weight: bold; -fx-font-size: 12px;");
        
        predictionConsole = new TextArea();
        predictionConsole.setEditable(false);
        predictionConsole.setWrapText(true);
        predictionConsole.setMinHeight(150);
        predictionConsole.setPrefHeight(180);
        predictionConsole.setMaxHeight(200);
        predictionConsole.setStyle("-fx-font-family: 'Courier New', monospace; -fx-font-size: 11px;");
        predictionConsole.setPromptText("Database status information will appear here...");
        predictionConsole.setId("predictionConsole");
        
        consolePanel.getChildren().addAll(consoleTitle, predictionConsole);
        
        // Prediction controls panel
        VBox predictionPanel = new VBox(5);
        predictionPanel.setPadding(new Insets(5));
        predictionPanel.setStyle("-fx-background-color: #f0f8ff; -fx-background-radius: 3; -fx-border-color: #b0d4f1; -fx-border-radius: 3; -fx-border-width: 1;");
        
        Label predictionLabel = new Label("Peptide Prediction Controls");
        predictionLabel.setStyle("-fx-font-weight: bold; -fx-font-size: 12px;");
        
        HBox predictionControls = new HBox(10);
        predictionControls.setAlignment(Pos.CENTER_LEFT);
        
        predictButton = new Button("Start Prediction");
        predictButton.setMinWidth(100);
        predictButton.setPrefWidth(120);
        predictButton.setMaxWidth(140);
        predictButton.setDisable(true); // Initially disabled until database is selected
        predictButton.setOnAction(e -> startPrediction());
        
        predictionProgressBar = new ProgressBar(0);
        predictionProgressBar.setMaxWidth(Double.MAX_VALUE);
        HBox.setHgrow(predictionProgressBar, Priority.ALWAYS);
        
        timeTextField = new TextField();
        timeTextField.setEditable(false);
        timeTextField.setMinWidth(100);
        timeTextField.setPrefWidth(120);
        timeTextField.setPromptText("Status");
        
        predictionControls.getChildren().addAll(predictButton, predictionProgressBar, timeTextField);
        
        // Update button
        updateButton = new Button("Check Status");
        updateButton.setMinWidth(80);
        updateButton.setPrefWidth(100);
        updateButton.setMaxWidth(120);
        updateButton.setOnAction(e -> updateMagDbStatus());
        
        HBox buttonPanel = new HBox(10);
        buttonPanel.setAlignment(Pos.CENTER_RIGHT);
        buttonPanel.getChildren().add(updateButton);
        
        predictionPanel.getChildren().addAll(predictionControls);
        
        statusPanel.getChildren().addAll(consolePanel, predictionPanel);
        
        // Initial status check will be done after database panel is initialized
        // This will be called from setAvailableMagDatabases() method
        
        return statusPanel;
    }

    @Override
    public VBox createParametersTab() {
        VBox tabContent = new VBox(10);
        tabContent.setPadding(new Insets(10));
        tabContent.setMinWidth(960);
        tabContent.setMinHeight(840);

        // Create main grid layout for parameters with 2-column responsive sizing
        GridPane mainGrid = new GridPane();
        mainGrid.setHgap(15);
        mainGrid.setVgap(15);
        mainGrid.setPadding(new Insets(10));
        mainGrid.setMinWidth(940);
        mainGrid.setMinHeight(740);
        VBox.setVgrow(mainGrid, Priority.ALWAYS);

        // Column constraints for 2-column layout
        ColumnConstraints col1 = new ColumnConstraints();
        col1.setPercentWidth(50);
        col1.setMinWidth(450);
        col1.setHgrow(Priority.ALWAYS);
        ColumnConstraints col2 = new ColumnConstraints();
        col2.setPercentWidth(50);
        col2.setMinWidth(450);
        col2.setHgrow(Priority.ALWAYS);
        mainGrid.getColumnConstraints().addAll(col1, col2);

        // Row constraints
        RowConstraints row1 = new RowConstraints();
        row1.setVgrow(Priority.ALWAYS);
        row1.setMinHeight(400);
        mainGrid.getRowConstraints().add(row1);

        // Left column: DIA Mode Selection, Thread Settings, Library Panel
        leftColumn = new VBox(15);
        leftColumn.setMinWidth(430);
        VBox diaSettingsSection = createDiaSettingsSection();
        
        // Create library generation panel placeholder - will be initialized after MAG database is set
        libraryPlaceholder = new VBox();
        libraryPlaceholder.setMinHeight(200);
        libraryPlaceholder.setStyle("-fx-background-color: #f0f0f0; -fx-border-color: #cccccc; -fx-border-radius: 5; -fx-border-width: 1;");
        Label placeholderLabel = new Label("Library Generation Panel\n(Will be available after MAG database is selected)");
        placeholderLabel.setStyle("-fx-text-alignment: center; -fx-font-style: italic; -fx-text-fill: #666666;");
        libraryPlaceholder.getChildren().add(placeholderLabel);
        libraryPlaceholder.setAlignment(Pos.CENTER);
        
        leftColumn.getChildren().addAll(diaSettingsSection, libraryPlaceholder);
        GridPane.setVgrow(leftColumn, Priority.ALWAYS);

        // Right column: Model Panel and Metadata Configuration
        VBox rightColumn = new VBox(15);
        rightColumn.setMinWidth(430);
        
        // Use the new reusable model generation panel
        modelGenerationPanel = new DiaModelGenerationPanel();
        modelGenerationPanel.setParameterUpdateCallback(() -> updateParameterView());
        modelGenerationPanel.setParentPanel(this); // Set reference to this panel
        // Ensure model panel has required resources and config file
        if (this.resourceParams == null) {
            try {
                MetaResourceIO resourceIO = new MetaResourceIO();
                this.resourceParams = resourceIO.loadParameters();
            } catch (Exception ignore) {
                // leave as null; panel will show an error if missing at run time
            }
        }
        modelGenerationPanel.setResourceParameters(this.resourceParams);
        modelGenerationPanel.setConfigFile(new File("MAGDB_1_0.json"));
        
        // Use the reusable metadata configuration panel
        metadataConfigurationPanel = new MetadataConfigurationPanel();
        metadataConfigurationPanel.setParameterUpdateCallback(() -> updateParameterView());
        
        rightColumn.getChildren().addAll(modelGenerationPanel, metadataConfigurationPanel);
        GridPane.setVgrow(rightColumn, Priority.ALWAYS);

        // Add columns to grid
        mainGrid.add(leftColumn, 0, 0);
        mainGrid.add(rightColumn, 1, 0);
        
        tabContent.getChildren().add(mainGrid);
        return tabContent;
    }

    private VBox createDiaSettingsSection() {
        VBox section = new VBox(10);
        section.setPadding(new Insets(10));
        section.setStyle("-fx-background-color: #f0f8ff; -fx-background-radius: 5;");

        Label titleLabel = createSectionTitle("DIA Analysis Settings");

        GridPane settingsGrid = new GridPane();
        settingsGrid.setHgap(15);
        settingsGrid.setVgap(8);

        ColumnConstraints col1 = new ColumnConstraints();
        col1.setMinWidth(150);
        col1.setHgrow(Priority.ALWAYS);
        ColumnConstraints col2 = new ColumnConstraints();
        col2.setMinWidth(150);
        col2.setHgrow(Priority.ALWAYS);
        settingsGrid.getColumnConstraints().addAll(col1, col2);

        // Mode Selection
        Label modeLabel = new Label("Analysis Mode:");
        modeLabel.setStyle("-fx-font-weight: bold;");
        modeSelectionComboBox = new ComboBox<>(FXCollections.observableArrayList(
            "Directly identification (no iterative search)",
            "Iterative search", 
            "Deep learning mode"));
        modeSelectionComboBox.setValue("Directly identification (no iterative search)");
        modeSelectionComboBox.setMinWidth(200);
        modeSelectionComboBox.setPrefWidth(300);
        modeSelectionComboBox.setMaxWidth(Double.MAX_VALUE);
        modeSelectionComboBox.setOnAction(e -> updateParameterView());

        VBox modeBox = new VBox(5);
        modeBox.getChildren().addAll(modeLabel, modeSelectionComboBox);

        // Thread Count
        Label threadsLabel = new Label("Thread Count:");
        threadsLabel.setStyle("-fx-font-weight: bold;");
        threadsComboBox = new ComboBox<>();
        for (int i = 1; i <= Runtime.getRuntime().availableProcessors(); i++) {
            threadsComboBox.getItems().add(i);
        }
        threadsComboBox.setValue(Math.min(4, Runtime.getRuntime().availableProcessors()));
        threadsComboBox.setMinWidth(100);
        threadsComboBox.setPrefWidth(120);
        threadsComboBox.setMaxWidth(Double.MAX_VALUE);
        threadsComboBox.setOnAction(e -> updateParameterView());

        VBox threadsBox = new VBox(5);
        threadsBox.getChildren().addAll(threadsLabel, threadsComboBox);

        settingsGrid.add(modeBox, 0, 0);
        settingsGrid.add(threadsBox, 1, 0);

        section.getChildren().addAll(titleLabel, settingsGrid);
        return section;
    }

    private Label createSectionTitle(String title) {
        Label label = new Label(title);
        label.setFont(Font.font("System", FontWeight.BOLD, 14));
        label.setTextFill(Color.DARKBLUE);
        return label;
    }


    private void updateMetadataTableFromIOPanel() {
        if (ioPanel == null || metadataConfigurationPanel == null) {
            return;
        }
        
        ObservableList<MetaPilotIOPanelFX.ExperimentRow> experiments = ioPanel.getExperimentData();
        metadataConfigurationPanel.updateMetadataFromIOPanel(experiments);
    }



    @Override
    protected void updateParameterView() {
        if (consolePanel == null) {
            return;
        }
        
        LOGGER.debug("updateParameterView() called");
        
        StringBuilder sb = new StringBuilder();
        sb.append("=== DIA Analysis Parameters ===\n\n");

        // Database settings
        sb.append("Database Settings:\n");
        if (magDatabasePanel != null) {
            MagDbItem selectedDb = magDatabasePanel.getSelectedMagDatabase();
            if (selectedDb != null) {
                sb.append("  MAG Database: ").append(selectedDb.getCatalogueID()).append("\n");
            }
            MagVersion selectedVersion = magDatabasePanel.getSelectedMagVersion();
            if (selectedVersion != null) {
                sb.append("  MAG Version: ").append(selectedVersion.getCatalogueVersion()).append("\n");
            }
            sb.append("  Append Host DB: ").append(magDatabasePanel.isAppendHostDb()).append("\n");
            if (magDatabasePanel.isAppendHostDb() && magDatabasePanel.getHostDatabase() != null && !magDatabasePanel.getHostDatabase().isEmpty()) {
                sb.append("  Host Database: ").append(magDatabasePanel.getHostDatabase()).append("\n");
            }
        }
        sb.append("\n");

        // DIA settings
        sb.append("DIA Settings:\n");
        if (modeSelectionComboBox != null) {
            sb.append("  Analysis Mode: ").append(modeSelectionComboBox.getValue()).append("\n");
        }
        if (threadsComboBox != null) {
            sb.append("  Thread Count: ").append(threadsComboBox.getValue()).append("\n");
        }
        sb.append("\n");

        // Library configuration
        sb.append("Spectral Libraries:\n");
        if (libraryGenerationPanel != null) {
            LOGGER.debug("updateParameterView - libraryGenerationPanel is not null");
            String[] selectedLibs = libraryGenerationPanel.getSelectedLibraryNames();
            LOGGER.debug("updateParameterView - selectedLibs.length = " + selectedLibs.length);
            if (selectedLibs.length > 0) {
                for (String libName : selectedLibs) {
                    sb.append("  ").append(libName).append("\n");
                }
            } else {
                sb.append("  No libraries selected\n");
            }
        } else {
            LOGGER.debug("updateParameterView - libraryGenerationPanel is null");
            sb.append("  Library panel not available\n");
        }
        sb.append("\n");

        // Model configuration
        sb.append("ML Models:\n");
        if (modelGenerationPanel != null) {
            String[] selectedModels = modelGenerationPanel.getSelectedModelNames();
            LOGGER.debug("updateParameterView - selectedModels.length = " + selectedModels.length);
            if (selectedModels.length > 0) {
                for (String modelName : selectedModels) {
                    sb.append("  ").append(modelName).append("\n");
                    LOGGER.debug("updateParameterView - selected model: " + modelName);
                }
            } else {
                sb.append("  No models selected (models are optional)\n");
                LOGGER.debug("updateParameterView - no models selected");
            }
        } else {
            LOGGER.debug("updateParameterView - modelGenerationPanel is null");
        }
        sb.append("\n");

        if (consolePanel != null) {
            consolePanel.setParameterText(sb.toString());
        }
        
        // Update the library generation panel with the current parameters
        if (libraryGenerationPanel != null) {
            // Try full validation first, fall back to minimal validation for library generation
            MetaParameterDia currentMetaParameter = collectParametersFromUI();
            if (currentMetaParameter == null) {
                // If full validation fails (e.g., no raw files), use minimal validation for library generation
                currentMetaParameter = generateMetaParameterDiaForLibraryGeneration();
            }
            if (currentMetaParameter != null) {
                libraryGenerationPanel.updateMetaParameter(currentMetaParameter);
            }
        } else {
            // Library generation panel not created yet, try to create it if MAG database is available
            ensureLibraryGenerationPanelCreated();
        }
        
        // Update the model generation panel with the current parameters
        if (modelGenerationPanel != null) {
            // Try full validation first, fall back to minimal validation for model generation
            MetaParameterDia currentMetaParameter = collectParametersFromUI();
            if (currentMetaParameter == null) {
                // If full validation fails (e.g., no raw files), use minimal validation for model generation
                currentMetaParameter = generateMetaParameterDiaForLibraryGeneration();
            }
            if (currentMetaParameter != null) {
                modelGenerationPanel.updateMetaParameter(currentMetaParameter);
            }
        }
    }


    @Override
    protected void startAnalysis() {
        try {
            consolePanel.appendConsoleText("Starting DIA Normal analysis...\n");
            LOGGER.info("startAnalysis called - DIA Normal workflow");
            
            // 1. Validate inputs (if applicable)
            if (!validateInputs()) {
                showError("Validation Error", "Please check all required inputs before starting analysis.");
                return;
            }
            
            // 2. Collect parameters from UI
            consolePanel.appendConsoleText("Collecting parameters from UI...\n");
            MetaParameterDia metaParameterDia = collectParametersFromUI();
            if (metaParameterDia == null) {
                consolePanel.appendConsoleText("ERROR: Failed to collect parameters from UI.\n");
                LOGGER.error("Failed to collect parameters from UI in startAnalysis");
                showError("Parameter Error", "Failed to generate analysis parameters. Please check your settings.");
                return;
            }
            
            consolePanel.appendConsoleText("✓ Successfully collected parameters from UI\n");
            LOGGER.info("Successfully collected parameters from UI for DIA Normal workflow");
            
            // 3. Store parameters for later use
            currentAnalysisParams = metaParameterDia;
            
            // 4. Load resource parameters from config file
            consolePanel.appendConsoleText("Loading resource parameters from config...\n");
            MetaResourceIO resourceIO = new MetaResourceIO();
            MetaResourceParameters mrp = resourceIO.loadParameters();
            if (mrp == null) {
                showError("Configuration Error", "Failed to load resource parameters from config file.");
                return;
            }
            consolePanel.appendConsoleText("✓ Resource parameters loaded from config file\n");
            
            // Clear console and start analysis
            consolePanel.clearConsole();
            consolePanel.appendConsoleText("Parameters validated successfully.\n");
            consolePanel.appendConsoleText("Result directory: " + metaParameterDia.getResult() + "\n");
        
            // 5. Create and configure task
            currentTask = new MetaDiaTask(
                metaParameterDia,
                mrp,
                consolePanel.getOverallProgressBar(),
                consolePanel.getCurrentTaskProgressBar()
            );
            
            // 6. Set up task event handlers
            currentTask.setOnSucceeded(e -> {
                Platform.runLater(() -> {
                    Boolean result = currentTask.getValue();
                    if (result != null && result) {
                        // Task completed successfully
                        consolePanel.appendConsoleText("✓ Analysis completed successfully!\n");
                        consolePanel.appendConsoleText("Results saved to: " + currentAnalysisParams.getResult() + "\n");
                        consolePanel.enableViewReportButton();
                        consolePanel.setButtonsEnabled(true, false, true);
                    } else {
                        // Task completed but returned false (failed)
                        consolePanel.appendConsoleText("✗ Analysis failed. Check console output for details.\n");
                        consolePanel.setButtonsEnabled(true, false, false);
                    }
                    consolePanel.stopStatusTimer(); // Stop the timer when task completes
                });
            });
            
            currentTask.setOnFailed(e -> {
                Platform.runLater(() -> {
                    Throwable exception = currentTask.getException();
                    String errorMessage = exception != null ? exception.getMessage() : "Unknown error occurred";
                    consolePanel.appendConsoleText("✗ Analysis failed: " + errorMessage + "\n");
                    consolePanel.setButtonsEnabled(true, false, false);
                    consolePanel.stopStatusTimer(); // Stop the timer when task fails
                    showError("Analysis Failed", "DIA Normal analysis failed: " + errorMessage);
                });
            });
            
            currentTask.setOnCancelled(e -> {
                Platform.runLater(() -> {
                    consolePanel.appendConsoleText("DIA Normal analysis was cancelled.\n");
                    consolePanel.setButtonsEnabled(true, false, false);
                    consolePanel.stopStatusTimer(); // Stop the timer when task is cancelled
                });
            });
            
            // 7. Run task in background thread
            Thread taskThread = new Thread(currentTask);
            taskThread.setDaemon(true);
            taskThread.start();
            
            consolePanel.appendConsoleText("✓ Analysis task started successfully.\n");
            LOGGER.info("DIA Normal analysis task started successfully");
            
        } catch (Exception e) {
            String errorMsg = "ERROR: Failed to start analysis: " + e.getMessage();
            consolePanel.appendConsoleText(errorMsg + "\n");
            LOGGER.error("Failed to start DIA Normal analysis", e);
            showError("Start Error", "Failed to start DIA Normal analysis: " + e.getMessage());
        }
    }

    @Override
    protected void stopAnalysis() {
        consolePanel.appendConsoleText("Stopping DIA analysis...\n");
        LOGGER.info("stopAnalysis called");
        
        if (currentTask != null && currentTask.isRunning()) {
            consolePanel.appendConsoleText("Force stopping DIA task...\n");
            currentTask.forceStop();
            consolePanel.appendConsoleText("DIA task stop requested.\n");
        } else {
            consolePanel.appendConsoleText("No running DIA task to stop.\n");
        }
    }
    
    
    /**
     * Update the MAG database status - check predicted peptide files and display in console
     */
    private void updateMagDbStatus() {
        if (predictionConsole == null) {
            return;
        }
        
        // Clear console and add timestamp
        predictionConsole.clear();
        predictionConsole.appendText("=== MAG Database Status Check ===\n");
        predictionConsole.appendText("Timestamp: " + LocalDateTime.now().toString() + "\n\n");
        
        if (magDatabasePanel == null) {
            predictionConsole.appendText("ERROR: Database panel not initialized\n");
            predictionConsole.appendText("Please restart the application.\n");
            updateButtonStates(false);
            return;
        }
        
        MagDbItem magDbItem = magDatabasePanel.getSelectedMagDatabase();
        if (magDbItem == null) {
            predictionConsole.appendText("INFO: No database selected\n");
            predictionConsole.appendText("Please select a MAG database from the dropdown above.\n");
            predictionConsole.appendText("If databases are available, they should appear in the dropdown.\n");
            updateButtonStates(false);
            return;
        }
        
        MagVersion usedVersion = magDatabasePanel.getSelectedMagVersion();
        if (usedVersion == null) {
            predictionConsole.appendText("INFO: No version selected for database: " + magDbItem.getCatalogueID() + "\n");
            predictionConsole.appendText("Please select a version from the version dropdown.\n");
            updateButtonStates(false);
            return;
        }
        
        // Display database information
        predictionConsole.appendText("Database Information:\n");
        predictionConsole.appendText("  Catalogue ID: " + magDbItem.getCatalogueID() + "\n");
        predictionConsole.appendText("  Version: " + usedVersion.getCatalogueVersion() + "\n");
        predictionConsole.appendText("  Species Count: " + usedVersion.getSpeciesCount() + "\n");
        predictionConsole.appendText("  Repository Path: " + (magDbItem.getRepositoryFile() != null ? magDbItem.getRepositoryFile().getAbsolutePath() : "Not set") + "\n\n");
        
        int speciesCount = usedVersion.getSpeciesCount();
        File predictedPepFile = magDbItem.getPredictedPepFile();
        
        // Check if predicted folder exists
        if (predictedPepFile == null || !predictedPepFile.exists()) {
            predictionConsole.appendText("STATUS: Predicted peptides folder does not exist\n");
            predictionConsole.appendText("Expected path: " + (predictedPepFile != null ? predictedPepFile.getAbsolutePath() : "null") + "\n");
            predictionConsole.appendText("ACTION: Click 'Start Prediction' to generate predicted peptides\n");
            updateButtonStates(false);
            return;
        }
        
        // Check if predicted folder contains the correct number of files
        File[] files = predictedPepFile.listFiles();
        if (files == null) {
            predictionConsole.appendText("ERROR: Cannot read predicted peptides folder\n");
            predictionConsole.appendText("Path: " + predictedPepFile.getAbsolutePath() + "\n");
            predictionConsole.appendText("ACTION: Check folder permissions\n");
            updateButtonStates(false);
            return;
        }
        
        int fileCount = files.length;
        predictionConsole.appendText("Predicted Peptides Status:\n");
        predictionConsole.appendText("  Expected files: " + speciesCount + "\n");
        predictionConsole.appendText("  Found files: " + fileCount + "\n");
        predictionConsole.appendText("  Completion: " + String.format("%.1f%%", (fileCount * 100.0 / speciesCount)) + "\n\n");
        
        if (fileCount == speciesCount) {
            predictionConsole.appendText("STATUS: ✓ All predicted peptides are available\n");
            predictionConsole.appendText("ACTION: Ready for analysis\n");
            updateButtonStates(true);
        } else if (fileCount > 0) {
            predictionConsole.appendText("STATUS: ⚠ Prediction incomplete (" + fileCount + "/" + speciesCount + " files)\n");
            predictionConsole.appendText("ACTION: Click 'Start Prediction' to continue\n");
            updateButtonStates(false);
        } else {
            predictionConsole.appendText("STATUS: ⚠ No predicted peptide files found\n");
            predictionConsole.appendText("ACTION: Click 'Start Prediction' to start prediction\n");
            updateButtonStates(false);
        }
    }
    
    /**
     * Update button states based on prediction status
     */
    private void updateButtonStates(boolean isComplete) {
        if (isComplete) {
            predictionProgressBar.setProgress(1.0);
            timeTextField.setText("Complete");
            predictButton.setDisable(true);
        } else {
            predictionProgressBar.setProgress(0);
            timeTextField.setText("Not Ready");
            predictButton.setDisable(false);
        }
    }
    
    // Getters for external access to status panel components
    public ProgressBar getPredictionProgressBar() {
        return predictionProgressBar;
    }
    
    public Button getPredictButton() {
        return predictButton;
    }
    
    public TextField getTimeTextField() {
        return timeTextField;
    }
    
    public Button getUpdateButton() {
        return updateButton;
    }
    
    public TextArea getPredictionConsole() {
        return predictionConsole;
    }
    
    /**
     * Start the peptide prediction process using DeepDetectTask
     */
    private void startPrediction() {
        if (magDatabasePanel == null || resourceParams == null) {
            predictionConsole.appendText("ERROR: Cannot start prediction - missing database or resource parameters\n");
            return;
        }
        
        MagDbItem magDbItem = magDatabasePanel.getSelectedMagDatabase();
        if (magDbItem == null) {
            predictionConsole.appendText("ERROR: No database selected for prediction\n");
            return;
        }

        // Check if prediction is already complete
        File predictedPepFile = magDbItem.getPredictedPepFile();
        if (predictedPepFile != null && predictedPepFile.exists()) {
            File[] files = predictedPepFile.listFiles();
            if (files != null && files.length == magDbItem.getUsedVersion().getSpeciesCount()) {
                predictionConsole.appendText("INFO: Prediction already complete. All files exist.\n");
                return;
            }
        }
        
        // Disable the predict button and update UI
        predictButton.setDisable(true);
        predictionProgressBar.setProgress(0);
        timeTextField.setText("Starting...");
        
        // Clear console and add prediction start message
        predictionConsole.clear();
        predictionConsole.appendText("=== Starting Peptide Prediction ===\n");
        predictionConsole.appendText("Timestamp: " + LocalDateTime.now().toString() + "\n");
        predictionConsole.appendText("Database: " + magDbItem.getCatalogueID() + " " + magDbItem.getUsedVersion().getCatalogueVersion() + "\n");
        predictionConsole.appendText("Species Count: " + magDbItem.getUsedVersion().getSpeciesCount() + "\n\n");
        
        MetaBaseParameter metaBaseParameter = new MetaBaseParameter();
        DeepDetectParameter deepDetectParam = new DeepDetectParameter();
        
        // Calculate thread count as 90% of available processors
        int availableProcessors = Runtime.getRuntime().availableProcessors();
        int threadCount = Math.max(1, (int) (availableProcessors * 0.6));
        
        predictionConsole.appendText("Thread Configuration:\n");
        predictionConsole.appendText("  Available Processors: " + availableProcessors + "\n");
        predictionConsole.appendText("  Using Threads: " + threadCount + " (60% of available)\n\n");
        
        // Create and configure DeepDetectTask with all required parameters
        deepDetectTask = new DeepDetectTask(
            metaBaseParameter,
            resourceParams,
            predictionProgressBar,
            deepDetectParam,
            magDbItem,
            threadCount, // Use 60% of available threads
            timeTextField // JavaFX TextField for time display
        );
        
        // Set up task completion callback
        deepDetectTask.setOnSucceeded(e -> {
            predictionConsole.appendText("\n=== Prediction Completed ===\n");
            predictionConsole.appendText("Timestamp: " + LocalDateTime.now().toString() + "\n");
            if (deepDetectTask.getValue()) {
                predictionConsole.appendText("STATUS: ✓ All predictions completed successfully\n");
                predictionConsole.appendText("ACTION: Ready for analysis\n");
                updateButtonStates(true);
            } else {
                predictionConsole.appendText("STATUS: ⚠ Prediction completed with some issues\n");
                predictionConsole.appendText("ACTION: Check individual files or restart prediction\n");
                updateButtonStates(false);
            }
            predictButton.setDisable(false);
        });
        
        deepDetectTask.setOnFailed(e -> {
            predictionConsole.appendText("\n=== Prediction Failed ===\n");
            predictionConsole.appendText("Timestamp: " + LocalDateTime.now().toString() + "\n");
            predictionConsole.appendText("ERROR: " + deepDetectTask.getException().getMessage() + "\n");
            predictionConsole.appendText("ACTION: Check error details and try again\n");
            updateButtonStates(false);
            predictButton.setDisable(false);
        });
        
        // Set up progress monitoring
        deepDetectTask.messageProperty().addListener((obs, oldMsg, newMsg) -> {
            if (newMsg != null && !newMsg.isEmpty()) {
                predictionConsole.appendText("Progress: " + newMsg + "\n");
            }
        });
        
        // Start the task in a background thread
        taskExecutor = Executors.newSingleThreadExecutor();
        taskExecutor.submit(deepDetectTask);
        
        predictionConsole.appendText("DeepDetect task started...\n");
        predictionConsole.appendText("Monitoring progress...\n\n");
    }
    
    /**
     * Generate MetaParameterDia for library generation (minimal validation)
     * This method creates a MetaParameterDia with only the essential parameters needed for library generation
     * 
     * @return MetaParameterDia object with minimal parameters, or null if essential validation fails
     */
    public MetaParameterDia generateMetaParameterDiaForLibraryGeneration() {
        System.out.println("DEBUG: generateMetaParameterDiaForLibraryGeneration() called");
        try {
            // Validate only essential inputs for library generation
            System.out.println("DEBUG: Starting minimal input validation for library generation...");
            if (!validateInputsForLibraryGeneration()) {
                System.err.println("Minimal input validation failed - cannot generate MetaParameterDia for library generation");
                return null;
            }
            System.out.println("DEBUG: Minimal input validation passed");
            
            // Create minimal MetaData (empty experiments for library generation)
            MetaData metaData = new MetaData(new String[0], new String[0], new int[0], new int[0]);
            
            // Get result directory (can be empty for library generation)
            String result = ioPanel != null ? ioPanel.getOutputDirectory() : "";
            
            // Get MAG database information
            MagDbItem selectedMagDb = magDatabasePanel.getSelectedMagDatabase();
            if (selectedMagDb == null) {
                System.err.println("No MAG database selected");
                return null;
            }
            
            MagVersion selectedVersion = magDatabasePanel.getSelectedMagVersion();
            if (selectedVersion == null) {
                System.err.println("No MAG database version selected");
                return null;
            }
            
            // Set the version on the MagDbItem - this is critical!
            selectedMagDb.setUsedVersion(selectedVersion);
            System.out.println("DEBUG: Set version " + selectedVersion.getCatalogueVersion() + " on MagDbItem " + selectedMagDb.getCatalogueID());
            
            // Get host database settings
            boolean appendHostDb = magDatabasePanel.isAppendHostDb();
            String hostDB = magDatabasePanel.getHostDatabase();
            
            // Get thread count
            int threadCount = threadsComboBox != null ? threadsComboBox.getValue() : 4;
            
            // Create MetaBaseParameter with minimal data
            MetaBaseParameter metaBasePar = new MetaBaseParameter(
                uk.ac.quadram.metapilot.task.MetaLabWorkflowType.DiaNNMAG.name(), // Workflow type
                metaData,
                result,
                selectedMagDb.getCurrentFile() != null ? selectedMagDb.getCurrentFile().getAbsolutePath() : "", // MicroDb path
                appendHostDb,
                hostDB,
                "DDA", // MS2 scan mode (default for DIA)
                1, // Core count
                threadCount,
                true // Meta workflow
            );
            
            // Create MetaParameterDia with minimal parameters for library generation
            MetaParameterDia metaParameterDia = new MetaParameterDia(
                metaBasePar,
                selectedMagDb,
                false, // librarySearch - will be set when libraries are actually selected
                false, // selfModelSearch - will be set when models are actually selected
                new String[0], // library - empty initially
                false, // useModel - will be set when models are actually selected
                new String[0] // models - empty initially
            );
            
            System.out.println("MetaParameterDia for library generation created successfully:");
            System.out.println("  - MAG Database: " + selectedMagDb.getCatalogueID());
            System.out.println("  - Threads: " + threadCount);
            
            return metaParameterDia;
            
        } catch (Exception e) {
            System.err.println("Error generating MetaParameterDia for library generation: " + e.getMessage());
            e.printStackTrace();
            return null;
        }
    }
    
    /**
     * Generate MetaParameterDia from the current UI parameters
     * This method collects all parameters from the UI and creates a MetaParameterDia object
     * 
     * @return MetaParameterDia object with all current UI parameters, or null if validation fails
     */
    public MetaParameterDia collectParametersFromUI() {
        try {
            consolePanel.appendConsoleText("Collecting parameters from UI...\n");
            LOGGER.info("Starting parameter collection from UI for DIA Normal workflow");
            
            // Validate inputs first
            consolePanel.appendConsoleText("Validating UI inputs...\n");
            LOGGER.info("Starting input validation...");
            if (!validateInputs()) {
                consolePanel.appendConsoleText("ERROR: Input validation failed\n");
                LOGGER.error("Input validation failed - cannot generate MetaParameterDia");
                return null;
            }
            consolePanel.appendConsoleText("✓ Input validation passed\n");
            LOGGER.info("Input validation passed");
            
            // Collect experiment data from ioPanel
            consolePanel.appendConsoleText("Collecting experiment data from IO panel...\n");
            ObservableList<MetaPilotIOPanelFX.ExperimentRow> fileRows = ioPanel.getExperimentData();
            consolePanel.appendConsoleText("Found " + fileRows.size() + " experiment rows\n");
            if (fileRows.isEmpty()) {
                consolePanel.appendConsoleText("ERROR: No experiment data found in IO panel\n");
                LOGGER.error("No experiment data found in IO panel");
                return null;
            }
            
            String[] rawFiles = new String[fileRows.size()];
            String[] expNames = new String[fileRows.size()];
            int[] fractions = new int[fileRows.size()];
            int[] replicates = new int[fileRows.size()];
            
            for (int i = 0; i < fileRows.size(); i++) {
                MetaPilotIOPanelFX.ExperimentRow row = fileRows.get(i);
                rawFiles[i] = row.getFilePath();
                expNames[i] = row.getExperimentName();
                // Get fraction from the IO panel, default to 1 if invalid
                try {
                    fractions[i] = Integer.parseInt(row.getFraction());
                    if (fractions[i] <= 0) {
                        fractions[i] = 1; // Ensure positive fraction
                    }
                } catch (NumberFormatException e) {
                    fractions[i] = 1; // Default to 1 if not a valid integer
                }
                replicates[i] = 1; // Default replicate
            }
            
            // Create MetaData object
            MetaData metaData = new MetaData(rawFiles, expNames, fractions, replicates);
            
        // Collect metadata information from the reusable panel
        String[][] metaInfo = new String[metadataConfigurationPanel.getMetadataTableData().size()][3];
        for (int i = 0; i < metadataConfigurationPanel.getMetadataTableData().size(); i++) {
            MetaPilotDdaIgcMaxQuantPanel.MetadataRow row = metadataConfigurationPanel.getMetadataTableData().get(i);
            metaInfo[i][0] = row.getAge() != null ? row.getAge() : "";
            metaInfo[i][1] = row.getGender() != null ? row.getGender() : "";
            metaInfo[i][2] = row.getCondition() != null ? row.getCondition() : "";
        }
        metaData.setMetaInfo(metaInfo);
            
            // Get result directory
            String result = ioPanel.getOutputDirectory();
            if (result == null || result.trim().isEmpty()) {
                System.err.println("Result directory not set");
                return null;
            }
            
            // Get MAG database information
            MagDbItem selectedMagDb = magDatabasePanel.getSelectedMagDatabase();
            if (selectedMagDb == null) {
                System.err.println("No MAG database selected");
                return null;
            }
            
            MagVersion selectedVersion = magDatabasePanel.getSelectedMagVersion();
            if (selectedVersion == null) {
                System.err.println("No MAG database version selected");
                return null;
            }
            
            // Set the version on the MagDbItem - this is critical!
            selectedMagDb.setUsedVersion(selectedVersion);
            System.out.println("DEBUG: Set version " + selectedVersion.getCatalogueVersion() + " on MagDbItem " + selectedMagDb.getCatalogueID());
            
            // Get host database settings
            boolean appendHostDb = magDatabasePanel.isAppendHostDb();
            String hostDB = magDatabasePanel.getHostDatabase();
            
            // Get thread count
            int threadCount = threadsComboBox != null ? threadsComboBox.getValue() : 4;
            
            // Create MetaBaseParameter
			MetaBaseParameter metaBasePar = new MetaBaseParameter(
					"DiaNNMAG", // Workflow type - using the correct DIA MAG workflow type
					metaData, result, selectedMagDb.getCurrentFile().getAbsolutePath(), // MicroDb path
					appendHostDb, hostDB, "DIA", // MS2 scan mode (default for DIA)
					1, // Core count
					threadCount, true // Meta workflow
			);
            
            // Get library information
            String[] library = new String[0];
            if (libraryGenerationPanel != null) {
                System.out.println("DEBUG: generateMetaParameterDia - libraryGenerationPanel is not null");
                String[] selectedLibs = libraryGenerationPanel.getSelectedLibraryNames();
                System.out.println("DEBUG: generateMetaParameterDia - selectedLibs.length = " + (selectedLibs != null ? selectedLibs.length : "null"));
                if (selectedLibs != null && selectedLibs.length > 0) {
                    library = selectedLibs;
                    System.out.println("DEBUG: generateMetaParameterDia - using " + library.length + " selected libraries");
                } else {
                    System.out.println("DEBUG: generateMetaParameterDia - no libraries selected");
                }
            } else {
                System.out.println("DEBUG: generateMetaParameterDia - libraryGenerationPanel is null");
            }
            
            // Get model information
            String[] models = new String[0];
            if (modelGenerationPanel != null) {
                String[] selectedModels = modelGenerationPanel.getSelectedModelNames();
                if (selectedModels != null && selectedModels.length > 0) {
                    models = selectedModels;
                }
            }
            
            // Determine analysis mode
            // librarySearch is true when the first item ("Directly identification (no iterative search)") is selected
            boolean librarySearch = modeSelectionComboBox.getSelectionModel().getSelectedIndex() == 0;
            boolean selfModelSearch = modeSelectionComboBox.getSelectionModel().getSelectedIndex() == 2;
            boolean useModel = models.length > 0;
            
            // Create MetaParameterDia
            MetaParameterDia metaParameterDia = new MetaParameterDia(
                metaBasePar,
                selectedMagDb,
                librarySearch,
                selfModelSearch,
                library,
                useModel,
                models
            );
            
            System.out.println("MetaParameterDia generated successfully:");
            System.out.println("  - Raw files: " + rawFiles.length);
            System.out.println("  - Libraries: " + library.length);
            System.out.println("  - Models: " + models.length);
            System.out.println("  - MAG Database: " + selectedMagDb.getCatalogueID());
            System.out.println("  - Threads: " + threadCount);
            
            consolePanel.appendConsoleText("✓ Successfully collected parameters from UI\n");
            LOGGER.info("Successfully collected parameters from UI for DIA Normal workflow");
            return metaParameterDia;
            
        } catch (Exception e) {
            String errorMsg = "ERROR: Failed to collect parameters from UI: " + e.getMessage();
            consolePanel.appendConsoleText(errorMsg + "\n");
            LOGGER.error("Error generating MetaParameterDia: " + e.getMessage(), e);
            return null;
        }
    }
    
    /**
     * Validate inputs for library generation (minimal validation)
     * Only validates essential parameters needed for library generation, not full analysis
     * 
     * @return true if essential inputs are valid, false otherwise
     */
    private boolean validateInputsForLibraryGeneration() {
        System.out.println("DEBUG: validateInputsForLibraryGeneration() called");
        
        // Check if MAG database is selected
        if (magDatabasePanel == null) {
            System.err.println("DEBUG: magDatabasePanel is null");
            return false;
        }
        if (magDatabasePanel.getSelectedMagDatabase() == null) {
            System.err.println("DEBUG: No MAG database selected");
            return false;
        }
        System.out.println("DEBUG: MAG database selected: " + magDatabasePanel.getSelectedMagDatabase().getCatalogueID());
        
        // Check if MAG database version is selected
        if (magDatabasePanel.getSelectedMagVersion() == null) {
            System.err.println("DEBUG: No MAG database version selected");
            return false;
        }
        System.out.println("DEBUG: MAG database version selected: " + magDatabasePanel.getSelectedMagVersion().getCatalogueVersion());
        
        System.out.println("DEBUG: All minimal validations passed for library generation");
        return true;
    }
    
    /**
     * Validate inputs before generating parameters
     * 
     * @return true if all required inputs are valid, false otherwise
     */
    private boolean validateInputs() {
        System.out.println("DEBUG: validateInputs() called");
        
        // Check if IO panel has experiment data
        if (ioPanel == null) {
            System.err.println("DEBUG: ioPanel is null");
            return false;
        }
        if (ioPanel.getExperimentData().isEmpty()) {
            System.err.println("DEBUG: No experiment data found in IO panel");
            return false;
        }
        System.out.println("DEBUG: IO panel has " + ioPanel.getExperimentData().size() + " experiments");
        
        // Check if result directory is set
        String outputDir = ioPanel.getOutputDirectory();
        if (outputDir == null || outputDir.trim().isEmpty()) {
            System.err.println("DEBUG: Result directory not set (value: " + outputDir + ")");
            return false;
        }
        System.out.println("DEBUG: Output directory is set: " + outputDir);
        
        // Check if MAG database is selected
        if (magDatabasePanel == null) {
            System.err.println("DEBUG: magDatabasePanel is null");
            return false;
        }
        if (magDatabasePanel.getSelectedMagDatabase() == null) {
            System.err.println("DEBUG: No MAG database selected");
            return false;
        }
        System.out.println("DEBUG: MAG database selected: " + magDatabasePanel.getSelectedMagDatabase().getCatalogueID());
        
        // Check if MAG database version is selected
        if (magDatabasePanel.getSelectedMagVersion() == null) {
            System.err.println("DEBUG: No MAG database version selected");
            return false;
        }
        System.out.println("DEBUG: MAG database version selected: " + magDatabasePanel.getSelectedMagVersion().getCatalogueVersion());
        
        System.out.println("DEBUG: All validations passed");
        return true;
    }
    
    /**
     * Load last used parameters and apply them to the UI
     */
    @Override
    protected void loadLastUsedParameters() {
        try {
            consolePanel.appendConsoleText("Loading last used parameters...\n");
            LOGGER.info("loadLastUsedParameters() called");
            MetaParameterDia params = (MetaParameterDia) parameterIO.loadLastUsedParameters();
            if (params != null) {
                consolePanel.appendConsoleText("✓ Successfully loaded last used parameters\n");
                LOGGER.info("Successfully loaded parameters from file");
                LOGGER.info("Parameters contain:");
                LOGGER.info("  - Output directory: {}", params.getResult() != null ? params.getResult() : "null");
                LOGGER.info("  - Raw files count: {}", params.getMetadata() != null && params.getMetadata().getRawFiles() != null ? params.getMetadata().getRawFiles().length : 0);
                LOGGER.info("  - MAG database: {}", params.getMagDbItem() != null ? params.getMagDbItem().getCatalogueID() : "null");
                LOGGER.info("  - Thread count: {}", params.getThreadCount());
                
                applyParametersToUI(params);
                LOGGER.info("Applied parameters to UI successfully");
            } else {
                consolePanel.appendConsoleText("No previous parameters found, using defaults\n");
                LOGGER.info("No parameters loaded from file, using defaults");
            }
        } catch (Exception e) {
            String errorMsg = "WARNING: Failed to load last used parameters: " + e.getMessage();
            consolePanel.appendConsoleText(errorMsg + "\n");
            LOGGER.warn("Failed to load last used DIA parameters: " + e.getMessage(), e);
            // Continue with default UI state
        }
    }
    
    /**
     * Apply parameters to the UI components
     */
    @Override
    protected void applyParametersToUI(MetaParameterDia params) {
        try {
            consolePanel.appendConsoleText("Applying parameters to UI...\n");
            if (params == null) {
                consolePanel.appendConsoleText("WARNING: No parameters to apply\n");
                LOGGER.warn("applyParametersToUI called with null params");
                return;
            }
            
            LOGGER.info("applyParametersToUI called with valid params");
        
        // Apply metadata to IO panel if available
        if (params.getMetadata() != null) {
            MetaData metadata = params.getMetadata();
            String[] rawFiles = metadata.getRawFiles();
            String[] expNames = metadata.getExpNames();
            
            System.out.println("DEBUG: Metadata found - rawFiles: " + (rawFiles != null ? rawFiles.length : "null") + 
                             ", expNames: " + (expNames != null ? expNames.length : "null"));
            
            if (rawFiles != null && expNames != null && rawFiles.length > 0) {
                ObservableList<MetaPilotIOPanelFX.ExperimentRow> experimentRows = FXCollections.observableArrayList();
                for (int i = 0; i < rawFiles.length; i++) {
                    if (rawFiles[i] != null && !rawFiles[i].trim().isEmpty()) {
                        experimentRows.add(new MetaPilotIOPanelFX.ExperimentRow(rawFiles[i], expNames[i]));
                        System.out.println("DEBUG: Added experiment row: " + rawFiles[i] + " -> " + expNames[i]);
                    }
                }
                if (ioPanel != null) {
                    ioPanel.setExperimentData(experimentRows);
                    System.out.println("DEBUG: Applied metadata with " + rawFiles.length + " raw files to IO panel");
                    
                    // Update the metadata configuration panel with the experiment data
                    if (metadataConfigurationPanel != null) {
                        metadataConfigurationPanel.updateMetadataFromIOPanel(experimentRows);
                    }
                } else {
                    System.err.println("ERROR: ioPanel is null, cannot apply experiment data");
                }
            } else {
                System.out.println("DEBUG: No valid raw files found in metadata");
            }
        } else {
            System.out.println("DEBUG: No metadata found in parameters");
        }
        
        // Apply output directory if available
        if (params.getResult() != null && !params.getResult().trim().isEmpty()) {
            if (ioPanel != null) {
                ioPanel.setOutputDirectory(params.getResult());
                System.out.println("DEBUG: Applied output directory: " + params.getResult());
            } else {
                System.err.println("ERROR: ioPanel is null, cannot apply output directory");
            }
        } else {
            System.out.println("DEBUG: No output directory found in parameters");
        }
        
        // Apply thread count
        if (threadsComboBox != null) {
            threadsComboBox.setValue(params.getThreadCount());
        }
        
        // Apply MAG database settings
        if (params.getMagDbItem() != null && magDatabasePanel != null) {
            MagDbItem magDbItem = params.getMagDbItem();
            // Find and select the matching MAG database
            for (MagDbItem item : magDatabasePanel.getAvailableMagDatabases()) {
                if (item.getCatalogueID().equals(magDbItem.getCatalogueID()) && 
                    item.getIdentifier().equals(magDbItem.getIdentifier())) {
                    magDatabasePanel.setSelectedMagDatabase(item);
                    
                    // Also apply version selection if available
                    if (magDbItem.getUsedVersion() != null) {
                        for (MagVersion version : item.getVersionList()) {
                            if (version.getCatalogueVersion().equals(magDbItem.getUsedVersion().getCatalogueVersion())) {
                                magDatabasePanel.setSelectedMagVersion(version);
                                System.out.println("DEBUG: Applied version selection: " + version.getCatalogueVersion());
                                break;
                            }
                        }
                    }
                    break;
                }
            }
            magDatabasePanel.setAppendHostDb(params.isAppendHostDb());
            if (params.getHostDB() != null) {
                magDatabasePanel.setHostDatabase(params.getHostDB());
            }
        }
        
        // Update parameter view after applying parameters
        Platform.runLater(() -> updateParameterView());
        
        consolePanel.appendConsoleText("✓ Successfully applied parameters to UI\n");
        LOGGER.info("Successfully applied parameters to UI for DIA Normal workflow");
        
        } catch (Exception e) {
            String errorMsg = "ERROR: Failed to apply parameters to UI: " + e.getMessage();
            consolePanel.appendConsoleText(errorMsg + "\n");
            LOGGER.error("Failed to apply parameters to UI for DIA Normal workflow", e);
        }
    }
    
    /**
     * Show error dialog
     */
    
    /**
     * Save current parameters to file
     */
    public void saveCurrentParameters() {
        try {
            // Try full validation first, fall back to minimal validation for library generation
            MetaParameterDia currentParams = collectParametersFromUI();
            if (currentParams == null) {
                // If full validation fails (e.g., no raw files), use minimal validation for library generation
                currentParams = generateMetaParameterDiaForLibraryGeneration();
            }
            
            if (currentParams != null) {
                boolean saved = parameterIO.saveLastUsedParameters(currentParams);
                if (saved) {
                    System.out.println("DEBUG: Successfully saved current parameters to last_used_dia_parameters.json");
                } else {
                    System.err.println("WARNING: Failed to save current parameters");
                }
            } else {
                System.err.println("WARNING: No parameters to save - both full and minimal validation failed");
            }
        } catch (Exception e) {
            System.err.println("ERROR: Failed to save current parameters: " + e.getMessage());
            e.printStackTrace();
        }
    }
    
    /**
     * Get the parameter class for this workflow
     */
    @Override
    protected Class<MetaParameterDia> getParameterClass() {
        return MetaParameterDia.class;
    }
} 
