package uk.ac.quadram.metapilot;

import javafx.application.Application;
import javafx.application.Platform;
import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.Scene;
import javafx.scene.control.*;
import javafx.scene.control.Alert;
import javafx.scene.control.Menu;
import javafx.scene.control.MenuBar;
import javafx.scene.control.MenuItem;
import javafx.scene.control.SeparatorMenuItem;
import javafx.scene.layout.*;
import javafx.stage.Stage;
import javafx.scene.image.Image;
import javafx.scene.image.ImageView;
import javafx.stage.FileChooser;
import java.io.File;
import java.util.ArrayList;
import java.util.prefs.Preferences;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.config.Configurator;
import org.apache.logging.log4j.core.config.builder.api.AppenderComponentBuilder;
import org.apache.logging.log4j.core.config.builder.api.ConfigurationBuilder;
import org.apache.logging.log4j.core.config.builder.api.ConfigurationBuilderFactory;
import org.apache.logging.log4j.core.config.builder.api.LayoutComponentBuilder;
import org.apache.logging.log4j.core.config.builder.api.RootLoggerComponentBuilder;
import org.apache.logging.log4j.core.config.builder.impl.BuiltConfiguration;

import uk.ac.quadram.metapilot.workflow.WorkflowType;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;
import uk.ac.quadram.metapilot.task.mag.MagDbItem.MagVersion;
import uk.ac.quadram.metapilot.task.mag.MagDbConfigIO;
import uk.ac.quadram.metapilot.task.mag.MagDbConfigPanelJavaFX;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.resources.MetaResourceIO;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters.ToolType;
import uk.ac.quadram.metapilot.ui.ResourceConfigPanel;
import uk.ac.quadram.metapilot.registration.RegistrationManager;
import uk.ac.quadram.metapilot.registration.UserRegistration;
import uk.ac.quadram.metapilot.ui.RegistrationDialog;
import uk.ac.quadram.metapilot.update.UpdateChecker;
import uk.ac.quadram.metapilot.task.dda.igc.gui.MetaPilotDdaIgcMaxQuantPanel;
import uk.ac.quadram.metapilot.task.dda.igc.gui.MetaPilotDdaIgcPFindPanel;
import uk.ac.quadram.metapilot.task.dda.igc.gui.MetaPilotDdaIgcFragPipePanel;
import uk.ac.quadram.metapilot.task.dda.igc.gui.MetaPilotDdaIgcSagePanel;
import uk.ac.quadram.metapilot.task.dda.mag.gui.MetaPilotDdaMagPfindPanel;
import uk.ac.quadram.metapilot.task.dda.mag.gui.MetaPilotDdaMagFragPipePanel;
import uk.ac.quadram.metapilot.task.dda.mag.gui.MetaPilotDdaMagSagePanel;
import uk.ac.quadram.metapilot.task.dia.gui.MetaPilotDiaNormalPanel;
import uk.ac.quadram.metapilot.task.dia.gui.MetaPilotDiaTargetPanel;
import uk.ac.quadram.metapilot.task.dia.gui.MetaPilotDiaRealTimePanel;
import uk.ac.quadram.metapilot.task.par.BaseParameterIO;
import uk.ac.quadram.metapilot.task.par.ParameterIORegistry;

/**
 * MetaPilot Main Application - A Comprehensive Metaproteomics Analysis Platform
 * 
 * <h2>Overview</h2>
 * MetaPilot is a sophisticated bioinformatics platform designed for metaproteomics
 * data analysis, supporting multiple mass spectrometry search engines and 
 * comprehensive taxonomic and functional annotation workflows.
 * 
 * <h2>Architecture</h2>
 * The application follows a modular JavaFX-based architecture with:
 * <ul>
 *   <li><strong>Task-based Processing:</strong> All analysis workflows are implemented 
 *       as background tasks extending {@link uk.ac.quadram.metapilot.task.MetaAbstractTask}</li>
 *   <li><strong>Resource Management:</strong> Centralized tool and database path management 
 *       through {@link uk.ac.quadram.metapilot.resources.MetaResourceParameters}</li>
 *   <li><strong>Workflow Integration:</strong> Support for multiple search engines including 
 *       MaxQuant, pFind, FragPipe, Sage, and DIA-NN</li>
 *   <li><strong>Database Configuration:</strong> MAG (Metagenome-Assembled Genomes) 
 *       database management for metagenomic analyses</li>
 * </ul>
 * 
 * <h2>Supported Workflows</h2>
 * <h3>Data-Dependent Acquisition (DDA) Workflows:</h3>
 * <ul>
 *   <li><strong>IGC (Integrated Gene Catalog) Workflows:</strong>
 *       <ul>
 *         <li>MaxQuant DDA IGC - Label-free quantification</li>
 *         <li>pFind DDA IGC - Chinese-developed search engine</li>
 *         <li>FragPipe DDA IGC - MSFragger-based pipeline</li>
 *         <li>Sage DDA IGC - Fast Rust-based search engine</li>
 *       </ul>
 *   </li>
 *   <li><strong>MAG (Metagenome-Assembled Genomes) Workflows:</strong>
 *       <ul>
 *         <li>pFind MAG - Targeted MAG database searches</li>
 *         <li>FragPipe MAG - MAG-specific FragPipe analysis</li>
 *         <li>Sage MAG - High-speed MAG proteome analysis</li>
 *       </ul>
 *   </li>
 * </ul>
 * 
 * <h3>Data-Independent Acquisition (DIA) Workflows:</h3>
 * <ul>
 *   <li><strong>DIA Normal:</strong> Standard DIA analysis with spectral libraries</li>
 *   <li><strong>DIA Target:</strong> Targeted DIA analysis for specific proteins</li>
 *   <li><strong>DIA Real-time:</strong> Real-time DIA analysis and monitoring</li>
 * </ul>
 * 
 * <h2>Core Features</h2>
 * <ul>
 *   <li><strong>Multi-engine Support:</strong> Integrates diverse search engines with 
 *       unified parameter management</li>
 *   <li><strong>Taxonomic Analysis:</strong> Comprehensive taxonomy assignment using 
 *       NCBI taxonomy and custom databases</li>
 *   <li><strong>Functional Annotation:</strong> Protein function prediction using 
 *       COG, KEGG, and custom functional databases</li>
 *   <li><strong>Quality Control:</strong> Built-in QC metrics and visualization</li>
 *   <li><strong>Database Management:</strong> Automated download and management of 
 *       reference databases</li>
 *   <li><strong>Registration System:</strong> User registration and license management</li>
 *   <li><strong>Update Management:</strong> Automatic software update checking</li>
 * </ul>
 * 
 * <h2>User Interface</h2>
 * The application provides a modern JavaFX interface with:
 * <ul>
 *   <li><strong>Main Workflow Panel:</strong> Easy access to all analysis workflows</li>
 *   <li><strong>Resource Configuration:</strong> Tool path management and validation</li>
 *   <li><strong>Database Configuration:</strong> MAG database setup and management</li>
 *   <li><strong>Progress Tracking:</strong> Real-time analysis progress monitoring</li>
 * </ul>
 * 
 * <h2>Technical Implementation</h2>
 * <ul>
 *   <li><strong>Threading:</strong> All CPU-intensive tasks run as background threads</li>
 *   <li><strong>Progress Monitoring:</strong> JavaFX Task-based progress tracking</li>
 *   <li><strong>Error Handling:</strong> Comprehensive exception handling with user feedback</li>
 *   <li><strong>Configuration Persistence:</strong> User preferences saved across sessions</li>
 * </ul>
 * 
 * <h2>Target Users</h2>
 * <ul>
 *   <li>Metaproteomics researchers analyzing microbiome samples</li>
 *   <li>Computational biologists working with mass spectrometry data</li>
 *   <li>Core facilities providing proteomics analysis services</li>
 *   <li>Students learning metaproteomics data analysis</li>
 * </ul>
 * 
 * @author Kai Cheng - Quadram Institute Bioscience
 * @version 1.0
 * @since 1.0
 * 
 * @see uk.ac.quadram.metapilot.task.MetaAbstractTask
 * @see uk.ac.quadram.metapilot.resources.MetaResourceParameters
 * @see uk.ac.quadram.metapilot.task.MetaLabWorkflowType
 */
public class MetaPilotApp extends Application {
    
    // Static block to configure application-level logging BEFORE Logger is initialized
    static {
        configureApplicationLogging();
    }
    
    private static final Logger LOGGER = LogManager.getLogger(MetaPilotApp.class);
    
    private MetaResourceParameters resourceParams;
    private MetaResourceIO resourceIO;
        
    // The MAG database items
    private MagDbItem[][] magDbItems;
    
    // Registration manager
    private RegistrationManager registrationManager;
    
    // Update checker
    private UpdateChecker updateChecker;
    
    // Registration email preference key
    private static final String PREF_REGISTERED_EMAIL = "registeredEmail";
    
    // Flag to control whether registration is required
    private static final boolean REGISTRATION_ENABLED = false; // Set to true when server is available
    
    /**
     * Configure application-level logging to logs/metapilot.log
     * This is called in static block before Logger instances are created.
     * Captures GUI startup, configuration loading, and errors before tasks start.
     */
    private static void configureApplicationLogging() {
        try {
            // Create logs directory if it doesn't exist
            File logsDir = new File("logs");
            if (!logsDir.exists()) {
                boolean created = logsDir.mkdirs();
                if (!created) {
                    System.err.println("Warning: Could not create logs directory: " + logsDir.getAbsolutePath());
                    return;
                }
            }
            
            // Build Log4j2 configuration programmatically
            ConfigurationBuilder<BuiltConfiguration> builder = 
                ConfigurationBuilderFactory.newConfigurationBuilder();
            
            builder.setStatusLevel(Level.INFO);
            builder.setConfigurationName("ApplicationLogger");
            
            // Pattern layout for log messages
            String pattern = "%d{yyyy-MM-dd HH:mm:ss} %-5level %logger{36} - %msg%n";
            LayoutComponentBuilder layoutBuilder = 
                builder.newLayout("PatternLayout")
                       .addAttribute("pattern", pattern);
            
            // File appender for application logs
            File logFile = new File(logsDir, "metapilot.log");
            AppenderComponentBuilder fileAppender = 
                builder.newAppender("ApplicationLogFile", "File")
                       .addAttribute("fileName", logFile.getAbsolutePath())
                       .addAttribute("append", true)
                       .add(layoutBuilder);
            
            // Console appender for immediate feedback
            AppenderComponentBuilder consoleAppender = 
                builder.newAppender("ApplicationConsole", "Console")
                       .addAttribute("target", "SYSTEM_OUT")
                       .add(layoutBuilder);
            
            builder.add(fileAppender);
            builder.add(consoleAppender);
            
            // Root logger - application level
            RootLoggerComponentBuilder rootLogger = builder.newRootLogger(Level.INFO);
            rootLogger.add(builder.newAppenderRef("ApplicationLogFile"));
            rootLogger.add(builder.newAppenderRef("ApplicationConsole"));
            builder.add(rootLogger);
            
            // Initialize configuration
            Configurator.initialize(builder.build());
            
            System.out.println("Application logging configured successfully to: " + logFile.getAbsolutePath());
            
        } catch (Exception e) {
            System.err.println("Failed to configure application logging: " + e.getMessage());
            e.printStackTrace();
        }
    }
    
    @Override
    public void start(Stage primaryStage) {
        LOGGER.info("MetaPilot GUI starting...");
        LOGGER.debug("Primary stage: {}", primaryStage.getTitle());
        
        // Initialize registration manager (but don't enforce registration yet)
        registrationManager = RegistrationManager.getInstance();
        
        // Initialize update checker
        updateChecker = new UpdateChecker(getHostServices());
        
        // Only check for updates, don't enforce registration at startup
        updateChecker.checkForUpdates(false);
        
        // Initialize the resource parameters
        resourceIO = new MetaResourceIO();
        resourceParams = resourceIO.loadParameters();

        // Load MAG database configuration
        try {
            // Try to load from the version file in the current directory
            File dbConfigFile = new File(new File("").getAbsolutePath(), MagDbConfigIO.versionFile);
            if (dbConfigFile.exists()) {
                magDbItems = MagDbConfigIO.parse(dbConfigFile);
            } else {
                // Initialize with empty configuration if file doesn't exist
                magDbItems = new MagDbItem[2][];
                magDbItems[0] = new MagDbItem[0]; // MGnify databases
                magDbItems[1] = new MagDbItem[0]; // Local databases
            }
        } catch (Exception e) {
            magDbItems = new MagDbItem[2][];
            magDbItems[0] = new MagDbItem[0];
            magDbItems[1] = new MagDbItem[0];
            System.err.println("Error loading MAG database configuration: " + e.getMessage());
        }
        
        // Create the main layout
        BorderPane root = new BorderPane();
        root.setPadding(new Insets(5));
        
        // Add menu bar
        MenuBar menuBar = createMenuBar(primaryStage);
        
        // Create a VBox to hold the menu bar and header
        VBox topContainer = new VBox(5);
        topContainer.getChildren().addAll(menuBar, createHeader());
        root.setTop(topContainer);
        
        // Create the main content panel with tabs
        TabPane tabPane = new TabPane();
        tabPane.setTabClosingPolicy(TabPane.TabClosingPolicy.UNAVAILABLE);
        
        // Create main tab content focused on workflows
        VBox mainContentPanel = new VBox(20);
        mainContentPanel.setPadding(new Insets(20, 10, 20, 10));
        
        // Create workflow selection
        TitledPane workflowPanel = createWorkflowPanel();
        
        // Add workflows to main content
        mainContentPanel.getChildren().add(workflowPanel);
        
        // Create Tab for Main menu
        Tab mainTab = new Tab("Main");
        mainTab.setContent(mainContentPanel);
        
        // Create Resources tab with the new ResourceConfigPanel
        Tab resourcesTab = new Tab("Resources");
        resourcesTab.setContent(createResourcesConfigPanel());
        
        // Create Database Configuration tab using SwingNode
        Tab databaseTab = new Tab("Database Configuration");
        databaseTab.setContent(createDatabaseConfigPanel());
        
        tabPane.getTabs().addAll(mainTab, resourcesTab, databaseTab);
        
        // Add a tab change listener to check for unsaved changes
        tabPane.getSelectionModel().selectedItemProperty().addListener((observable, oldTab, newTab) -> {
            if (oldTab == resourcesTab) {
                // When leaving resources tab, check for unsaved changes and update resourceParams
                Pane content = (Pane) oldTab.getContent();
                if (content instanceof ResourceConfigPanel) {
                    ResourceConfigPanel configPanel = (ResourceConfigPanel) content;
                    if (!configPanel.checkUnsavedChanges()) {
                        // User cancelled, go back to resources tab
                        Platform.runLater(() -> tabPane.getSelectionModel().select(oldTab));
                    } else {
                        // User didn't cancel, update resourceParams with current state
                        this.resourceParams = configPanel.getResourceParameters();
                        System.out.println("DEBUG: Resource parameters updated when leaving Resources tab");
                    }
                }
            }
        });
        
        root.setCenter(tabPane);
        
        // Create the scene
        Scene scene = new Scene(root, 1600, 1000);
        try {
            scene.getStylesheets().add(getClass().getResource("/styles/metapilot.css").toExternalForm());
        } catch (Exception e) {
            System.err.println("Could not load stylesheet: " + e.getMessage());
        }
        
        // Set up the stage
        primaryStage.setTitle("MetaPilot");
        primaryStage.setScene(scene);
        try {
            primaryStage.getIcons().add(new Image(getClass().getResourceAsStream("/MetaPilot.png")));
        } catch (Exception e) {
            System.err.println("Could not load icon: " + e.getMessage());
        }
        
        // Set window close handler to save configurations before exit
        primaryStage.setOnCloseRequest(event -> {
            try {
                // Check if there are unsaved changes in the resource panel
                Tab resourcesTabRef = tabPane.getTabs().get(1); // Resources tab is the second tab
                if (resourcesTabRef != null && resourcesTabRef.getContent() instanceof ResourceConfigPanel) {
                    ResourceConfigPanel configPanel = (ResourceConfigPanel) resourcesTabRef.getContent();
                    if (!configPanel.checkUnsavedChanges()) {
                        // User cancelled, prevent window from closing
                        event.consume();
                        return;
                    }
                }
                
                // Save the MAG database configuration
                saveMagDbConfiguration();
                
                // Save the resource parameters
                resourceIO.saveParameters(resourceParams);
            } catch (Exception e) {
                System.err.println("Failed to save configurations: " + e.getMessage());
            }
        });
        
        primaryStage.show();
    }
    
    private Pane createDatabaseConfigPanel() {
        // Create a JavaFX container to hold the database configuration panel
        BorderPane container = new BorderPane();
        container.setPadding(new Insets(20));
        container.getStyleClass().add("section-panel");
        
        // Add a JavaFX title label at the top
        Label titleLabel = new Label("Database Configuration");
        titleLabel.getStyleClass().add("section-title");
        container.setTop(titleLabel);
        BorderPane.setMargin(titleLabel, new Insets(0, 0, 20, 0));
        
        // Create our JavaFX implementation of the MagDbConfigPanel
        MagDbConfigPanelJavaFX dbPanel = new MagDbConfigPanelJavaFX(magDbItems);
        
        // Add the panel to the container
        container.setCenter(dbPanel);
        
        // Register a listener to save configuration when the app closes
        Platform.runLater(() -> {
            Stage stage = (Stage) container.getScene().getWindow();
            stage.setOnCloseRequest(event -> {
                magDbItems = dbPanel.getMagDbItems();
                saveMagDbConfiguration();
            });
        });
        
        return container;
    }
    
    private void saveMagDbConfiguration() {
        try {
            MagDbConfigIO.export(magDbItems, new File(new File("").getAbsolutePath(), MagDbConfigIO.versionFile));
        } catch (Exception e) {
            System.err.println("Error saving MAG database configuration: " + e.getMessage());
        }
    }
    
    private HBox createHeader() {
        HBox header = new HBox(20);
        header.setAlignment(Pos.CENTER_LEFT);
        header.setPadding(new Insets(0, 0, 20, 0));
        
        // Create title label
        Label titleLabel = new Label("MetaPilot");
        titleLabel.setStyle("-fx-font-size: 32px; -fx-font-weight: bold;");
        
        // Add a subtitle
        Label subtitleLabel = new Label("Next-generation Metaproteomics Analysis Platform");
        subtitleLabel.setStyle("-fx-font-size: 18px;");
        
        // Create a container for the title and subtitle
        VBox titleBox = new VBox(5, titleLabel, subtitleLabel);
        
        // Try to load and add the logo
        try {
            ImageView logoView = new ImageView(new Image(getClass().getResourceAsStream("/MetaPilot.png")));
            logoView.setFitHeight(80);
            logoView.setPreserveRatio(true);
            header.getChildren().add(logoView);
        } catch (Exception e) {
            // If logo loading fails, add a placeholder rectangle
            StackPane placeholder = new StackPane();
            placeholder.setPrefSize(80, 80);
            placeholder.setStyle("-fx-background-color: lightgray;");
            header.getChildren().add(placeholder);
        }
        
        header.getChildren().add(titleBox);
        
        return header;
    }
    
    private Pane createResourcesConfigPanel() {
        // Create the resource configuration panel
        ResourceConfigPanel configPanel = new ResourceConfigPanel();
        
        // Add a close listener to update our reference to the current resource parameters
        Platform.runLater(() -> {
            Stage stage = (Stage) configPanel.getScene().getWindow();
            stage.setOnCloseRequest(event -> {
                configPanel.saveOnClose();
                this.resourceParams = configPanel.getResourceParameters();
            });
        });
        
        return configPanel;
    }
    
    private TitledPane createWorkflowPanel() {
        TitledPane workflowPane = new TitledPane();
        workflowPane.setText("Workflows");
        workflowPane.setCollapsible(false);
        
        ScrollPane scrollPane = new ScrollPane();
        scrollPane.setFitToWidth(true);
        
        VBox workflowContainer = new VBox(20);
        workflowContainer.setPadding(new Insets(15));
        
        workflowContainer.getChildren().addAll(
            createDdaWorkflowPanel(),
            createDiaWorkflowPanel()
        );
        
        scrollPane.setContent(workflowContainer);
        workflowPane.setContent(scrollPane);
        
        return workflowPane;
    }
    
    private TitledPane createDdaWorkflowPanel() {
        TitledPane panel = new TitledPane();
        panel.setText("DDA Workflows");
        panel.setCollapsible(false);
        
        VBox workflowBox = new VBox(20);
        workflowBox.setPadding(new Insets(15));
        workflowBox.setAlignment(Pos.CENTER);
        
        // IGC Database Workflows
        TitledPane igcPanel = new TitledPane();
        igcPanel.setText("IGC Database");
        igcPanel.setCollapsible(false);
        
        HBox igcWorkflows = new HBox(20);
        igcWorkflows.setPadding(new Insets(10));
        igcWorkflows.setAlignment(Pos.CENTER);
        
        igcWorkflows.getChildren().addAll(
            createWorkflowBox("MaxQuant", "Run DDA workflow using MaxQuant", 
                e -> launchWorkflow(WorkflowType.DDA_IGC_MAXQUANT)),
            createWorkflowBox("pFind", "Run DDA workflow using pFind", 
                e -> launchWorkflow(WorkflowType.DDA_IGC_PFIND)),
            createWorkflowBox("FragPipe", "Run DDA workflow using FragPipe", 
                e -> launchWorkflow(WorkflowType.DDA_IGC_FRAGPIPE)),
            createWorkflowBox("Sage", "Run DDA workflow using Sage", 
                e -> launchWorkflow(WorkflowType.DDA_IGC_SAGE))
        );
        
        igcPanel.setContent(igcWorkflows);
        
        // MAG Database Workflows
        TitledPane magPanel = new TitledPane();
        magPanel.setText("MAG Database");
        magPanel.setCollapsible(false);
        
        HBox magWorkflows = new HBox(20);
        magWorkflows.setPadding(new Insets(10));
        magWorkflows.setAlignment(Pos.CENTER);
        
        magWorkflows.getChildren().addAll(
            createWorkflowBox("pFind", "Run DDA workflow using pFind", 
                e -> launchWorkflow(WorkflowType.DDA_MAG_PFIND)),
            createWorkflowBox("FragPipe", "Run DDA workflow using FragPipe", 
                e -> launchWorkflow(WorkflowType.DDA_MAG_FRAGPIPE)),
            createWorkflowBox("Sage", "Run DDA workflow using Sage", 
                e -> launchWorkflow(WorkflowType.DDA_MAG_SAGE))
        );
        
        magPanel.setContent(magWorkflows);
        
        workflowBox.getChildren().addAll(igcPanel, magPanel);
        panel.setContent(workflowBox);
        
        return panel;
    }
    
    private TitledPane createDiaWorkflowPanel() {
        TitledPane panel = new TitledPane();
        panel.setText("DIA Workflows");
        panel.setCollapsible(false);
        
        HBox workflowBox = new HBox(20);
        workflowBox.setPadding(new Insets(15));
        workflowBox.setAlignment(Pos.CENTER);
        
        workflowBox.getChildren().addAll(
            createWorkflowBox("Normal", "Run DIA workflow in normal mode", 
                e -> launchWorkflow(WorkflowType.DIA_NORMAL)),
            createWorkflowBox("Target", "Run DIA workflow in target mode", 
                e -> launchWorkflow(WorkflowType.DIA_TARGET)),
            createWorkflowBox("Real-time", "Run DIA workflow in real-time mode", 
                e -> launchWorkflow(WorkflowType.DIA_REALTIME))
        );
        
        panel.setContent(workflowBox);
        
        return panel;
    }
    
    private VBox createWorkflowBox(String title, String description, javafx.event.EventHandler<javafx.event.ActionEvent> action) {
        VBox box = new VBox(10);
        box.setAlignment(Pos.CENTER);
        box.setPadding(new Insets(15));
        box.setMinWidth(200);
        box.setMaxWidth(300);
        box.getStyleClass().add("workflow-box");
        
        Label titleLabel = new Label(title);
        titleLabel.setStyle("-fx-font-size: 16px; -fx-font-weight: bold;");
        
        Label descLabel = new Label(description);
        descLabel.setWrapText(true);
        descLabel.setAlignment(Pos.CENTER);
        
        Button launchButton = new Button("Launch");
        launchButton.setOnAction(action);
        
        box.getChildren().addAll(titleLabel, descLabel, launchButton);
        
        return box;
    }
    
    
    private boolean checkRequiredResources(WorkflowType type) {
        switch (type) {
            // DDA IGC Workflows
            case DDA_IGC_MAXQUANT:
                return resourceParams.isToolAvailable(ToolType.MAXQUANT);
            case DDA_IGC_PFIND:
                return resourceParams.isToolAvailable(ToolType.PFIND);
            case DDA_IGC_FRAGPIPE:
                return resourceParams.isToolAvailable(ToolType.FRAGPIPE);
            case DDA_IGC_SAGE:
                return resourceParams.isToolAvailable(ToolType.SAGE);
            
            // DDA MAG Workflows
            case DDA_MAG_PFIND:
                return resourceParams.isToolAvailable(ToolType.PFIND);
            case DDA_MAG_FRAGPIPE:
                return resourceParams.isToolAvailable(ToolType.FRAGPIPE);
            case DDA_MAG_SAGE:
                return resourceParams.isToolAvailable(ToolType.SAGE);
            
            // DIA Workflows
            case DIA_NORMAL:
            case DIA_TARGET:
            case DIA_REALTIME:
                return resourceParams.isToolAvailable(ToolType.DIANN);
            
            default:
                return false;
        }
    }
    
    private void launchWorkflow(WorkflowType type) {
        // First check if registration is required and enabled
        if (REGISTRATION_ENABLED && !checkRegistration()) {
            showWarning("Registration Required", 
                      "Registration is required to use this workflow. " +
                      "Please register the application first.");
            return;
        }
        
        // Then check if all required resources are available
        if (!checkRequiredResources(type)) {
            showWarning("Missing Resources", 
                      "Some required resources are missing for this workflow. " +
                      "Please configure all required resources first.");
            return;
        }
        
        ArrayList<MagDbItem> list = new ArrayList<MagDbItem>();
		if (magDbItems != null && magDbItems.length > 0) {
			for (int i = 0; i < magDbItems.length; i++) {
				for (int j = 0; j < magDbItems[i].length; j++) {
					boolean use = false;
					ArrayList<MagVersion> verlist = magDbItems[i][j].getVersionList();
					for (int k = 0; k < verlist.size(); k++) {
						if (verlist.get(k).isLocalAvailable()) {
							use = true;
							break;
						}
					}

					if (use) {
						list.add(magDbItems[i][j]);
					}
				}
			}
		}
		MagDbItem[] items = list.toArray(new MagDbItem[list.size()]);
        
        try {
            switch (type) {
                case DDA_IGC_MAXQUANT:
                    // Create a new stage for the MaxQuant workflow
                    Stage maxQuantStage = new Stage();
                    maxQuantStage.setTitle("DDA IGC MaxQuant Workflow");
                    
                    // Create the MaxQuant panel
                    MetaPilotDdaIgcMaxQuantPanel maxQuantPanel = new MetaPilotDdaIgcMaxQuantPanel();
                    
                    // Create a scene with the panel
                    Scene scene = new Scene(maxQuantPanel);
                    try {
                        scene.getStylesheets().add(getClass().getResource("/styles/metapilot.css").toExternalForm());
                    } catch (Exception e) {
                        System.err.println("Could not load stylesheet: " + e.getMessage());
                    }
                    
                    // Set up the stage
                    maxQuantStage.setScene(scene);
                    maxQuantStage.setMinWidth(800);
                    maxQuantStage.setMinHeight(600);
                    
                    // Show the stage
                    maxQuantStage.setMaximized(true);
                    maxQuantStage.show();
                    break;
                    
                case DDA_IGC_PFIND:
                    // Create a new stage for the pFind workflow
                    Stage pFindStage = new Stage();
                    pFindStage.setTitle("DDA IGC pFind Workflow");

                    // Create the pFind panel
                    MetaPilotDdaIgcPFindPanel pFindPanel = new MetaPilotDdaIgcPFindPanel();

                    // Create a scene with the panel
                    Scene pFindScene = new Scene(pFindPanel);
                    try {
                        pFindScene.getStylesheets().add(getClass().getResource("/styles/metapilot.css").toExternalForm());
                    } catch (Exception e) {
                        System.err.println("Could not load stylesheet: " + e.getMessage());
                    }

                    // Set up the stage
                    pFindStage.setScene(pFindScene);
                    pFindStage.setMinWidth(800);
                    pFindStage.setMinHeight(600);

                    // Show the stage
                    pFindStage.setMaximized(true);
                    pFindStage.show();
                    break;
                    
                case DDA_IGC_FRAGPIPE:
                    // Create a new stage for the FragPipe workflow
                    Stage fragPipeStage = new Stage();
                    fragPipeStage.setTitle("DDA IGC FragPipe Workflow");

                    // Create the FragPipe panel
                    MetaPilotDdaIgcFragPipePanel fragPipePanel = new MetaPilotDdaIgcFragPipePanel();

                    // Create a scene with the panel
                    Scene fragPipeScene = new Scene(fragPipePanel);
                    try {
                        fragPipeScene.getStylesheets().add(getClass().getResource("/styles/metapilot.css").toExternalForm());
                    } catch (Exception e) {
                        System.err.println("Could not load stylesheet: " + e.getMessage());
                    }

                    // Set up the stage
                    fragPipeStage.setScene(fragPipeScene);
                    fragPipeStage.setMinWidth(800);
                    fragPipeStage.setMinHeight(600);

                    // Show the stage
                    fragPipeStage.setMaximized(true);
                    fragPipeStage.show();
                    break;
                    
                case DDA_IGC_SAGE:
                    // Create a new stage for the IGC Sage workflow
                    Stage sageStage = new Stage();
                    sageStage.setTitle("DDA IGC Sage Workflow");

                    // Create the Sage panel
                    MetaPilotDdaIgcSagePanel sagePanel = new MetaPilotDdaIgcSagePanel();

                    // Create a scene with the panel
                    Scene sageScene = new Scene(sagePanel);
                    try {
                        sageScene.getStylesheets().add(getClass().getResource("/styles/metapilot.css").toExternalForm());
                    } catch (Exception e) {
                        System.err.println("Could not load stylesheet: " + e.getMessage());
                    }

                    // Set up the stage
                    sageStage.setScene(sageScene);
                    sageStage.setMinWidth(800);
                    sageStage.setMinHeight(600);

                    // Show the stage
                    sageStage.setMaximized(true);
                    sageStage.show();
                    break;
                    
                case DDA_MAG_PFIND:
                    // Launch the MAG main frame
                    // Note: Using default constructor if parameters are not supported
                	// Create a new stage for the pFind workflow
                    Stage pFindMagStage = new Stage();
                    pFindMagStage.setTitle("DDA MAG pFind Workflow");

                    // Create the pFind panel
                    MetaPilotDdaMagPfindPanel pFindMagPanel = new MetaPilotDdaMagPfindPanel(items);

                    // Create a scene with the panel
                    Scene pFindMagScene = new Scene(pFindMagPanel);
                    try {
                    	pFindMagScene.getStylesheets().add(getClass().getResource("/styles/metapilot.css").toExternalForm());
                    } catch (Exception e) {
                        System.err.println("Could not load stylesheet: " + e.getMessage());
                    }

                    // Set up the stage
                    pFindMagStage.setScene(pFindMagScene);
                    pFindMagStage.setMinWidth(800);
                    pFindMagStage.setMinHeight(600);

                    // Show the stage
                    pFindMagStage.setMaximized(true);
                    pFindMagStage.show();
                    break;
                    
                case DDA_MAG_FRAGPIPE:
                    // Create a new stage for the MAG FragPipe workflow
                    Stage fragPipeMagStage = new Stage();
                    fragPipeMagStage.setTitle("DDA MAG FragPipe Workflow");

                    // Create the MAG FragPipe panel
                    MetaPilotDdaMagFragPipePanel fragPipeMagPanel = new MetaPilotDdaMagFragPipePanel(items, resourceParams);

                    // Create a scene with the panel
                    Scene fragPipeMagScene = new Scene(fragPipeMagPanel);
                    try {
                        fragPipeMagScene.getStylesheets().add(getClass().getResource("/styles/metapilot.css").toExternalForm());
                    } catch (Exception e) {
                        System.err.println("Could not load stylesheet: " + e.getMessage());
                    }

                    // Set up the stage
                    fragPipeMagStage.setScene(fragPipeMagScene);
                    fragPipeMagStage.setMinWidth(800);
                    fragPipeMagStage.setMinHeight(600);

                    // Show the stage
                    fragPipeMagStage.setMaximized(true);
                    fragPipeMagStage.show();
                    break;
                    
                case DDA_MAG_SAGE:
                    Stage magSageStage = new Stage();
                    magSageStage.setTitle("DDA MAG Sage Workflow");
                    MetaPilotDdaMagSagePanel magSagePanel = new MetaPilotDdaMagSagePanel(items);
                    Scene magSageScene = new Scene(magSagePanel, 1600, 900);
                    try {
                        magSageScene.getStylesheets().add(getClass().getResource("/styles/metapilot.css").toExternalForm());
                    } catch (Exception e) {
                        System.err.println("Could not load stylesheet: " + e.getMessage());
                    }
                    magSageStage.setScene(magSageScene);
                    magSageStage.setMinWidth(800);
                    magSageStage.setMinHeight(600);
                    magSageStage.setMaximized(true);
                    magSageStage.show();
                    break;
                    
                case DIA_NORMAL:
                    Stage diaNormalStage = new Stage();
                    diaNormalStage.setTitle("DIA Normal Workflow");
                    MetaPilotDiaNormalPanel diaNormalPanel = new MetaPilotDiaNormalPanel(items, resourceParams);
                    Scene diaNormalScene = new Scene(diaNormalPanel, 1600, 900);
                    try {
                        diaNormalScene.getStylesheets().add(getClass().getResource("/styles/metapilot.css").toExternalForm());
                    } catch (Exception e) {
                        System.err.println("Could not load stylesheet: " + e.getMessage());
                    }
                    diaNormalStage.setScene(diaNormalScene);
                    diaNormalStage.setMinWidth(800);
                    diaNormalStage.setMinHeight(600);
                    diaNormalStage.setMaximized(true);
                    
                    // Add close handler to save parameters when window closes
                    diaNormalStage.setOnCloseRequest(e -> {
                        System.out.println("DEBUG: DIA Normal window closing - saving parameters");
                        diaNormalPanel.saveCurrentParameters();
                    });
                    
                    diaNormalStage.show();
                    break;
                    
                case DIA_TARGET:
                    // Launch the DIA Target main frame
                    Stage diaTargetStage = new Stage();
                    diaTargetStage.setTitle("DIA Target Workflow");
                    MetaPilotDiaTargetPanel diaTargetPanel = new MetaPilotDiaTargetPanel(items, resourceParams);
                    Scene diaTargetScene = new Scene(diaTargetPanel, 1600, 900);
                    try {
                        diaTargetScene.getStylesheets().add(getClass().getResource("/styles/metapilot.css").toExternalForm());
                    } catch (Exception e) {
                        System.err.println("Could not load stylesheet: " + e.getMessage());
                    }
                    diaTargetStage.setScene(diaTargetScene);
                    diaTargetStage.setMinWidth(800);
                    diaTargetStage.setMinHeight(600);
                    diaTargetStage.setMaximized(true);
                    
                    // Add close handler to save parameters
                    diaTargetStage.setOnCloseRequest(event -> {
                        try {
                            diaTargetPanel.saveCurrentParameters();
                        } catch (Exception e) {
                            System.err.println("Error saving DIA Target parameters on close: " + e.getMessage());
                        }
                    });
                    
                    diaTargetStage.show();
                    break;
                    
                case DIA_REALTIME:
                    // Launch the DIA Real-time main frame
                    Stage diaRealTimeStage = new Stage();
                    diaRealTimeStage.setTitle("DIA Real-time Workflow");
                    MetaPilotDiaRealTimePanel diaRealTimePanel = new MetaPilotDiaRealTimePanel(items, resourceParams);
                    Scene diaRealTimeScene = new Scene(diaRealTimePanel, 1600, 900);
                    try {
                        diaRealTimeScene.getStylesheets().add(getClass().getResource("/styles/metapilot.css").toExternalForm());
                    } catch (Exception e) {
                        System.err.println("Could not load stylesheet: " + e.getMessage());
                    }
                    diaRealTimeStage.setScene(diaRealTimeScene);
                    diaRealTimeStage.setMinWidth(800);
                    diaRealTimeStage.setMinHeight(600);
                    diaRealTimeStage.setMaximized(true);
                    
                    // Add close handler to save parameters when window closes
                    diaRealTimeStage.setOnCloseRequest(e -> {
                        System.out.println("DEBUG: DIA Real-time window closing - saving parameters");
                        diaRealTimePanel.saveCurrentParameters();
                    });
                    
                    diaRealTimeStage.show();
                    break;
            }
        } catch (Exception ex) {
            showError("Workflow Launch Error", 
                    "Failed to launch workflow", ex);
        }
    }
    
    /**
     * Load parameter file and launch the corresponding workflow panel
     */
    private void loadParameterFile(Stage primaryStage) {
        // Initialize parameter IO registry if not already done
        ParameterIORegistry.registerAllParameterIOs();
        
        // Create file chooser
        FileChooser fileChooser = new FileChooser();
        fileChooser.setTitle("Load Parameter File");
        fileChooser.getExtensionFilters().addAll(
            new FileChooser.ExtensionFilter("JSON Files", "*.json"),
            new FileChooser.ExtensionFilter("All Files", "*.*")
        );
        
        // Set initial directory to config folder if it exists
        File configDir = new File("config");
        if (configDir.exists() && configDir.isDirectory()) {
            fileChooser.setInitialDirectory(configDir);
        }
        
        // Show file chooser dialog
        File selectedFile = fileChooser.showOpenDialog(primaryStage);
        if (selectedFile == null) {
            return; // User cancelled
        }
        
        try {
            System.out.println("DEBUG: MetaPilotApp.loadParameterFile - Loading file: " + selectedFile.getAbsolutePath());
            
            // Load parameters from file
            Object parameters = BaseParameterIO.loadParametersFromFilePath(selectedFile.getAbsolutePath());
            System.out.println("DEBUG: MetaPilotApp.loadParameterFile - Loaded parameters: " + (parameters != null ? parameters.getClass().getSimpleName() : "NULL"));
            
            if (parameters == null) {
                System.out.println("DEBUG: MetaPilotApp.loadParameterFile - Parameters is null, showing warning");
                showWarning("Load Parameter File", 
                    "Failed to load parameters from file. The file may be corrupted or not a valid parameter file.");
                return;
            }
            
            // Get workflow type from parameters
            System.out.println("DEBUG: MetaPilotApp.loadParameterFile - Getting workflow type from parameters...");
            String workflowType = BaseParameterIO.getWorkflowTypeFromParameters(parameters);
            System.out.println("DEBUG: MetaPilotApp.loadParameterFile - Workflow type: '" + workflowType + "'");
            
            if (workflowType == null) {
                System.out.println("DEBUG: MetaPilotApp.loadParameterFile - Workflow type is null, showing warning");
                showWarning("Load Parameter File", 
                    "Could not determine workflow type from parameter file.");
                return;
            }
            
            // Map workflow type string to WorkflowType enum
            WorkflowType type = getWorkflowTypeFromString(workflowType);
            if (type == null) {
                showWarning("Load Parameter File", 
                    "Unsupported workflow type: " + workflowType);
                return;
            }
            
            // Launch the workflow with loaded parameters
            launchWorkflowWithParameters(type, parameters);
            
        } catch (Exception e) {
            showError("Load Parameter File Error", 
                "Failed to load parameter file: " + e.getMessage(), e);
        }
    }
    
    /**
     * Map workflow type string to WorkflowType enum
     */
    private WorkflowType getWorkflowTypeFromString(String workflowType) {
        System.out.println("DEBUG: getWorkflowTypeFromString called with: '" + workflowType + "'");
        
        // First try direct matching with analysis types
        for (WorkflowType type : WorkflowType.values()) {
            if (type.getAnalysisType().equals(workflowType)) {
                System.out.println("DEBUG: Found direct match: " + type);
                return type;
            }
        }
        
        // Handle DIA workflow type mapping from JSON to WorkflowType enum
        switch (workflowType) {
            case "DiaNNMAG":
                System.out.println("DEBUG: Mapping DiaNNMAG to DIA_NORMAL");
                return WorkflowType.DIA_NORMAL;
            case "DiaNNMAGRT":
                System.out.println("DEBUG: Mapping DiaNNMAGRT to DIA_REALTIME");
                return WorkflowType.DIA_REALTIME;
            case "DiaNNMAGTG":
                System.out.println("DEBUG: Mapping DiaNNMAGTG to DIA_TARGET");
                return WorkflowType.DIA_TARGET;
            default:
                System.out.println("DEBUG: No mapping found for workflow type: '" + workflowType + "'");
                return null;
        }
    }
    
    /**
     * Launch workflow with pre-loaded parameters
     */
    private void launchWorkflowWithParameters(WorkflowType type, Object parameters) {
        // First check if registration is required and enabled
        if (REGISTRATION_ENABLED && !checkRegistration()) {
            showWarning("Registration Required", 
                      "Registration is required to use this workflow. " +
                      "Please register the application first.");
            return;
        }
        
        // Then check if all required resources are available
        if (!checkRequiredResources(type)) {
            showWarning("Missing Resources", 
                      "Some required resources are missing for this workflow. " +
                      "Please configure all required resources first.");
            return;
        }
        
        ArrayList<MagDbItem> list = new ArrayList<MagDbItem>();
        if (magDbItems != null && magDbItems.length > 0) {
            for (int i = 0; i < magDbItems.length; i++) {
                for (int j = 0; j < magDbItems[i].length; j++) {
                    boolean use = false;
                    ArrayList<MagVersion> verlist = magDbItems[i][j].getVersionList();
                    for (int k = 0; k < verlist.size(); k++) {
                        if (verlist.get(k).isLocalAvailable()) {
                            use = true;
                            break;
                        }
                    }

                    if (use) {
                        list.add(magDbItems[i][j]);
                    }
                }
            }
        }
        MagDbItem[] items = list.toArray(new MagDbItem[list.size()]);
        
        try {
            switch (type) {
                case DDA_IGC_MAXQUANT:
                    launchDdaIgcMaxQuantWithParameters(parameters);
                    break;
                case DDA_IGC_PFIND:
                    launchDdaIgcPFindWithParameters(parameters);
                    break;
                case DDA_IGC_FRAGPIPE:
                    launchDdaIgcFragPipeWithParameters(parameters);
                    break;
                case DDA_IGC_SAGE:
                    launchDdaIgcSageWithParameters(parameters);
                    break;
                case DDA_MAG_PFIND:
                    launchDdaMagPFindWithParameters(parameters, items);
                    break;
                case DDA_MAG_FRAGPIPE:
                    launchDdaMagFragPipeWithParameters(parameters, items);
                    break;
                case DDA_MAG_SAGE:
                    launchDdaMagSageWithParameters(parameters, items);
                    break;
                case DIA_NORMAL:
                    launchDiaNormalWithParameters(parameters, items);
                    break;
                case DIA_TARGET:
                    launchDiaTargetWithParameters(parameters, items);
                    break;
                case DIA_REALTIME:
                    launchDiaRealTimeWithParameters(parameters, items);
                    break;
            }
        } catch (Exception ex) {
            showError("Workflow Launch Error", 
                    "Failed to launch workflow with parameters", ex);
        }
    }
    
    // Individual workflow launch methods with parameters
    
    private void launchDdaIgcMaxQuantWithParameters(Object parameters) {
        Stage maxQuantStage = new Stage();
        maxQuantStage.setTitle("DDA IGC MaxQuant Workflow");
        
        MetaPilotDdaIgcMaxQuantPanel maxQuantPanel = new MetaPilotDdaIgcMaxQuantPanel();
        
        Scene scene = new Scene(maxQuantPanel);
        try {
            scene.getStylesheets().add(getClass().getResource("/styles/metapilot.css").toExternalForm());
        } catch (Exception e) {
            System.err.println("Could not load stylesheet: " + e.getMessage());
        }
        
        maxQuantStage.setScene(scene);
        maxQuantStage.setMinWidth(800);
        maxQuantStage.setMinHeight(600);
        maxQuantStage.setMaximized(true);
        maxQuantStage.show();
        
        // Apply loaded parameters to the panel
        maxQuantPanel.loadParametersFromObject(parameters);
    }
    
    private void launchDdaIgcPFindWithParameters(Object parameters) {
        Stage pFindStage = new Stage();
        pFindStage.setTitle("DDA IGC pFind Workflow");
        
        MetaPilotDdaIgcPFindPanel pFindPanel = new MetaPilotDdaIgcPFindPanel();
        
        Scene pFindScene = new Scene(pFindPanel);
        try {
            pFindScene.getStylesheets().add(getClass().getResource("/styles/metapilot.css").toExternalForm());
        } catch (Exception e) {
            System.err.println("Could not load stylesheet: " + e.getMessage());
        }
        
        pFindStage.setScene(pFindScene);
        pFindStage.setMinWidth(800);
        pFindStage.setMinHeight(600);
        pFindStage.setMaximized(true);
        pFindStage.show();
        
        // Apply loaded parameters to the panel
        pFindPanel.loadParametersFromObject(parameters);
    }
    
    private void launchDdaIgcFragPipeWithParameters(Object parameters) {
        Stage fragPipeStage = new Stage();
        fragPipeStage.setTitle("DDA IGC FragPipe Workflow");
        
        MetaPilotDdaIgcFragPipePanel fragPipePanel = new MetaPilotDdaIgcFragPipePanel();
        
        Scene fragPipeScene = new Scene(fragPipePanel);
        try {
            fragPipeScene.getStylesheets().add(getClass().getResource("/styles/metapilot.css").toExternalForm());
        } catch (Exception e) {
            System.err.println("Could not load stylesheet: " + e.getMessage());
        }
        
        fragPipeStage.setScene(fragPipeScene);
        fragPipeStage.setMinWidth(800);
        fragPipeStage.setMinHeight(600);
        fragPipeStage.setMaximized(true);
        fragPipeStage.show();
        
        // Apply loaded parameters to the panel
        fragPipePanel.loadParametersFromObject(parameters);
    }
    
    private void launchDdaIgcSageWithParameters(Object parameters) {
        Stage sageStage = new Stage();
        sageStage.setTitle("DDA IGC Sage Workflow");
        
        MetaPilotDdaIgcSagePanel sagePanel = new MetaPilotDdaIgcSagePanel();
        
        Scene sageScene = new Scene(sagePanel);
        try {
            sageScene.getStylesheets().add(getClass().getResource("/styles/metapilot.css").toExternalForm());
        } catch (Exception e) {
            System.err.println("Could not load stylesheet: " + e.getMessage());
        }
        
        sageStage.setScene(sageScene);
        sageStage.setMinWidth(800);
        sageStage.setMinHeight(600);
        sageStage.setMaximized(true);
        sageStage.show();
        
        // Apply loaded parameters to the panel
        sagePanel.loadParametersFromObject(parameters);
    }
    
    private void launchDdaMagPFindWithParameters(Object parameters, MagDbItem[] items) {
        Stage pFindMagStage = new Stage();
        pFindMagStage.setTitle("DDA MAG pFind Workflow");
        
        MetaPilotDdaMagPfindPanel pFindMagPanel = new MetaPilotDdaMagPfindPanel(items);
        
        Scene pFindMagScene = new Scene(pFindMagPanel);
        try {
            pFindMagScene.getStylesheets().add(getClass().getResource("/styles/metapilot.css").toExternalForm());
        } catch (Exception e) {
            System.err.println("Could not load stylesheet: " + e.getMessage());
        }
        
        pFindMagStage.setScene(pFindMagScene);
        pFindMagStage.setMinWidth(800);
        pFindMagStage.setMinHeight(600);
        pFindMagStage.setMaximized(true);
        pFindMagStage.show();
        
        // Apply loaded parameters to the panel
        pFindMagPanel.loadParametersFromObject(parameters);
    }
    
    private void launchDdaMagFragPipeWithParameters(Object parameters, MagDbItem[] items) {
        Stage fragPipeMagStage = new Stage();
        fragPipeMagStage.setTitle("DDA MAG FragPipe Workflow");
        
        MetaPilotDdaMagFragPipePanel fragPipeMagPanel = new MetaPilotDdaMagFragPipePanel(items, resourceParams);
        
        Scene fragPipeMagScene = new Scene(fragPipeMagPanel);
        try {
            fragPipeMagScene.getStylesheets().add(getClass().getResource("/styles/metapilot.css").toExternalForm());
        } catch (Exception e) {
            System.err.println("Could not load stylesheet: " + e.getMessage());
        }
        
        fragPipeMagStage.setScene(fragPipeMagScene);
        fragPipeMagStage.setMinWidth(800);
        fragPipeMagStage.setMinHeight(600);
        fragPipeMagStage.setMaximized(true);
        fragPipeMagStage.show();
        
        // Apply loaded parameters to the panel
        fragPipeMagPanel.loadParametersFromObject(parameters);
    }
    
    private void launchDdaMagSageWithParameters(Object parameters, MagDbItem[] items) {
        Stage magSageStage = new Stage();
        magSageStage.setTitle("DDA MAG Sage Workflow");
        
        MetaPilotDdaMagSagePanel magSagePanel = new MetaPilotDdaMagSagePanel(items);
        
        Scene magSageScene = new Scene(magSagePanel, 1600, 900);
        try {
            magSageScene.getStylesheets().add(getClass().getResource("/styles/metapilot.css").toExternalForm());
        } catch (Exception e) {
            System.err.println("Could not load stylesheet: " + e.getMessage());
        }
        
        magSageStage.setScene(magSageScene);
        magSageStage.setMinWidth(800);
        magSageStage.setMinHeight(600);
        magSageStage.setMaximized(true);
        magSageStage.show();
        
        // Apply loaded parameters to the panel
        magSagePanel.loadParametersFromObject(parameters);
    }
    
    private void launchDiaNormalWithParameters(Object parameters, MagDbItem[] items) {
        Stage diaNormalStage = new Stage();
        diaNormalStage.setTitle("DIA Normal Workflow");
        
        MetaPilotDiaNormalPanel diaNormalPanel = new MetaPilotDiaNormalPanel(items, resourceParams);
        
        Scene diaNormalScene = new Scene(diaNormalPanel, 1600, 900);
        try {
            diaNormalScene.getStylesheets().add(getClass().getResource("/styles/metapilot.css").toExternalForm());
        } catch (Exception e) {
            System.err.println("Could not load stylesheet: " + e.getMessage());
        }
        
        diaNormalStage.setScene(diaNormalScene);
        diaNormalStage.setMinWidth(800);
        diaNormalStage.setMinHeight(600);
        diaNormalStage.setMaximized(true);
        
        // Add close handler to save parameters when window closes
        diaNormalStage.setOnCloseRequest(e -> {
            System.out.println("DEBUG: DIA Normal window closing - saving parameters");
            diaNormalPanel.saveCurrentParameters();
        });
        
        diaNormalStage.show();
        
        // Apply loaded parameters to the panel
        diaNormalPanel.loadParametersFromObject(parameters);
    }
    
    private void launchDiaTargetWithParameters(Object parameters, MagDbItem[] items) {
        Stage diaTargetStage = new Stage();
        diaTargetStage.setTitle("DIA Target Workflow");
        
        MetaPilotDiaTargetPanel diaTargetPanel = new MetaPilotDiaTargetPanel(items, resourceParams);
        
        Scene diaTargetScene = new Scene(diaTargetPanel, 1600, 900);
        try {
            diaTargetScene.getStylesheets().add(getClass().getResource("/styles/metapilot.css").toExternalForm());
        } catch (Exception e) {
            System.err.println("Could not load stylesheet: " + e.getMessage());
        }
        
        diaTargetStage.setScene(diaTargetScene);
        diaTargetStage.setMinWidth(800);
        diaTargetStage.setMinHeight(600);
        diaTargetStage.setMaximized(true);
        
        // Add close handler to save parameters
        diaTargetStage.setOnCloseRequest(event -> {
            try {
                diaTargetPanel.saveCurrentParameters();
            } catch (Exception e) {
                System.err.println("Error saving DIA Target parameters on close: " + e.getMessage());
            }
        });
        
        diaTargetStage.show();
        
        // Apply loaded parameters to the panel
        diaTargetPanel.loadParametersFromObject(parameters);
    }
    
    private void launchDiaRealTimeWithParameters(Object parameters, MagDbItem[] items) {
        Stage diaRealTimeStage = new Stage();
        diaRealTimeStage.setTitle("DIA Real-time Workflow");
        
        MetaPilotDiaRealTimePanel diaRealTimePanel = new MetaPilotDiaRealTimePanel(items, resourceParams);
        
        Scene diaRealTimeScene = new Scene(diaRealTimePanel, 1600, 900);
        try {
            diaRealTimeScene.getStylesheets().add(getClass().getResource("/styles/metapilot.css").toExternalForm());
        } catch (Exception e) {
            System.err.println("Could not load stylesheet: " + e.getMessage());
        }
        
        diaRealTimeStage.setScene(diaRealTimeScene);
        diaRealTimeStage.setMinWidth(800);
        diaRealTimeStage.setMinHeight(600);
        diaRealTimeStage.setMaximized(true);
        
        // Add close handler to save parameters when window closes
        diaRealTimeStage.setOnCloseRequest(e -> {
            System.out.println("DEBUG: DIA Real-time window closing - saving parameters");
            diaRealTimePanel.saveCurrentParameters();
        });
        
        diaRealTimeStage.show();
        
        // Apply loaded parameters to the panel (DIA Real-time has its own parameter loading mechanism)
        if (parameters instanceof uk.ac.quadram.metapilot.task.dia.par.MetaParameterDiaRT) {
            diaRealTimePanel.applyParametersToUI((uk.ac.quadram.metapilot.task.dia.par.MetaParameterDiaRT) parameters);
        }
    }
    
    private void showWarning(String title, String content) {
        Alert alert = new Alert(Alert.AlertType.WARNING);
        alert.setTitle(title);
        alert.setHeaderText(null);
        alert.setContentText(content);
        alert.showAndWait();
    }
    
    private void showError(String title, String content, Exception ex) {
        // Print detailed error information to console for debugging
        System.err.println("=== ERROR DETAILS ===");
        System.err.println("Title: " + title);
        System.err.println("Content: " + content);
        System.err.println("Exception: " + ex.getClass().getSimpleName());
        System.err.println("Message: " + ex.getMessage());
        System.err.println("Stack trace:");
        ex.printStackTrace();
        System.err.println("=====================");
        
        Alert alert = new Alert(Alert.AlertType.ERROR);
        alert.setTitle(title);
        alert.setHeaderText(null);
        alert.setContentText(content + "\n\nError: " + ex.getMessage());
        alert.showAndWait();
    }
    
    /**
     * Add a menu to the app
     */
    private MenuBar createMenuBar(Stage primaryStage) {
        MenuBar menuBar = new MenuBar();
        
        // File menu
        Menu fileMenu = new Menu("File");
        MenuItem loadParamItem = new MenuItem("Load Parameter File...");
        loadParamItem.setOnAction(e -> loadParameterFile(primaryStage));
        MenuItem exitItem = new MenuItem("Exit");
        exitItem.setOnAction(e -> Platform.exit());
        fileMenu.getItems().addAll(loadParamItem, new SeparatorMenuItem(), exitItem);
        
        // Tools menu
        Menu toolsMenu = new Menu("Tools");
        MenuItem resourcesItem = new MenuItem("Configure Resources...");
        resourcesItem.setOnAction(e -> {
            Scene scene = primaryStage.getScene();
            TabPane tabPane = (TabPane)((BorderPane)scene.getRoot()).getCenter();
            Tab resourcesTab = tabPane.getTabs().get(1);
            tabPane.getSelectionModel().select(resourcesTab);
        });
        
        MenuItem databaseItem = new MenuItem("Configure Database...");
        databaseItem.setOnAction(e -> {
            Scene scene = primaryStage.getScene();
            TabPane tabPane = (TabPane)((BorderPane)scene.getRoot()).getCenter();
            Tab databaseTab = tabPane.getTabs().get(2);
            tabPane.getSelectionModel().select(databaseTab);
        });
        
        toolsMenu.getItems().addAll(resourcesItem, databaseItem);
        
        // Help menu
        Menu helpMenu = new Menu("Help");
        MenuItem checkUpdateItem = new MenuItem("Check for Updates...");
        checkUpdateItem.setOnAction(e -> updateChecker.checkForUpdates(true));
        
        MenuItem aboutItem = new MenuItem("About");
        aboutItem.setOnAction(e -> showAboutDialog());
        
        helpMenu.getItems().addAll(checkUpdateItem, new SeparatorMenuItem(), aboutItem);
        
        // Add all menus
        menuBar.getMenus().addAll(fileMenu, toolsMenu, helpMenu);
        
        return menuBar;
    }
    
    /**
     * Show the about dialog
     */
    private void showAboutDialog() {
        Alert alert = new Alert(Alert.AlertType.INFORMATION);
        alert.setTitle("About MetaPilot");
        alert.setHeaderText("MetaPilot v" + UpdateChecker.getCurrentVersion());
        alert.setContentText(
            "Comprehensive Workflow Manager for Metaproteomics\n\n" +
            "© 2025 Quadram Institute Bioscience\n" +
            "All rights reserved.\n\n" +
            "Registered to: " + Preferences.userRoot().node(getClass().getPackage().getName()).get(PREF_REGISTERED_EMAIL, "")
        );
        alert.showAndWait();
    }
    
    /**
     * Check if the user is registered
     * @return true if registration is complete or successful
     */
    private boolean checkRegistration() {
        // Get the preferences node
        Preferences prefs = Preferences.userRoot().node(getClass().getPackage().getName());
        
        // Check if email is stored in preferences
        String email = prefs.get(PREF_REGISTERED_EMAIL, null);
        
        // If email is found, check if it's activated
        if (email != null && !email.isEmpty()) {
            if (registrationManager.isRegisteredAndActivated(email)) {
                // User is registered and activated
                return true;
            }
        }
        
        // User is not registered or not activated, show registration dialog
        RegistrationDialog dialog = new RegistrationDialog();
        
        // If we know the email, populate it
        if (email != null && !email.isEmpty()) {
            dialog.showExistingRegistration(email, "");
        }
        
        // Show the dialog and wait for result
        UserRegistration result = dialog.showAndWait().orElse(null);
        
        // If user completed registration
        if (result != null) {
            // Store email in preferences
            prefs.put(PREF_REGISTERED_EMAIL, result.getEmail());
            return true;
        }
        
        // Registration was cancelled
        return false;
    }
    
    @Override
    public void stop() {
        try {
            // Close registration database connection
            if (registrationManager != null) {
                registrationManager.close();
            }
            
            // ... any other cleanup code ...
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public static void main(String[] args) {
        // Application logging already configured by static block
        LOGGER.info("MetaPilot application starting...");
        LOGGER.debug("Command line arguments: {}", java.util.Arrays.toString(args));
        
        // Check if running in CLI mode
        if (args.length > 0 && args[0].equals("--cli")) {
            LOGGER.info("Starting in CLI mode");
            // Remove --cli from args and pass to CLI class
            String[] cliArgs = new String[args.length - 1];
            System.arraycopy(args, 1, cliArgs, 0, cliArgs.length);
            MetaPilotCLI.main(cliArgs);
        } else {
            LOGGER.info("Starting in GUI mode");
            launch(args);
        }
    }

} 