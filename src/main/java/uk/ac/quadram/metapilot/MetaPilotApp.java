package uk.ac.quadram.metapilot;

import javafx.application.Application;
import javafx.application.Platform;
import javafx.embed.swing.SwingNode;
import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.Scene;
import javafx.scene.control.*;
import javafx.scene.layout.*;
import javafx.stage.DirectoryChooser;
import javafx.stage.FileChooser;
import javafx.stage.Stage;
import javafx.scene.image.Image;
import javafx.scene.image.ImageView;
import javafx.scene.paint.Color;
import javafx.scene.shape.Circle;
import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.util.HashMap;
import java.util.Map;
import java.util.prefs.Preferences;
import javax.swing.SwingUtilities;
import javax.swing.JPanel;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import javafx.concurrent.Task;

import uk.ac.quadram.metapilot.task.pfind.gui.MetaPilotPFindMainFrame;
import uk.ac.quadram.metapilot.task.mag.gui.MetaPilotMainFrameMag;
import uk.ac.quadram.metapilot.task.dia.gui.MetaPilotMainFrameDia;
import uk.ac.quadram.metapilot.task.dia.target.gui.MetaPilotMainFrameDiaTg;
import uk.ac.quadram.metapilot.task.dia.monitor.gui.MetaPilotMainFrameDiaRT;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;
import uk.ac.quadram.metapilot.task.mag.MagDbConfigIO;
import uk.ac.quadram.metapilot.task.mag.gui.MagDbConfigPanel;
import uk.ac.quadram.metapilot.task.mag.gui.MagDbConfigPanelJavaFX;
import uk.ac.quadram.metapilot.core.MGYG.MGYGFtpTask;
import javafx.scene.control.cell.TextFieldTableCell;

import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.resources.MetaResourceIO;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters.ToolType;
import uk.ac.quadram.metapilot.ui.ResourceConfigPanel;
import uk.ac.quadram.metapilot.workflows.ResourceConfigWorkflow;

public class MetaPilotApp extends Application {
    
    private MetaResourceParameters resourceParams;
    private MetaResourceIO resourceIO;
    
    // The MAG database items
    private MagDbItem[][] magDbItems;
    
    @Override
    public void start(Stage primaryStage) {
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
        root.setPadding(new Insets(20));
        
        // Create a header with title and logo
        HBox header = createHeader();
        root.setTop(header);
        
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
                // When leaving resources tab, check for unsaved changes
                Pane content = (Pane) oldTab.getContent();
                if (content instanceof ResourceConfigPanel) {
                    ResourceConfigPanel configPanel = (ResourceConfigPanel) content;
                    if (!configPanel.checkUnsavedChanges()) {
                        // User cancelled, go back to resources tab
                        Platform.runLater(() -> tabPane.getSelectionModel().select(oldTab));
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
            primaryStage.getIcons().add(new Image(getClass().getResourceAsStream("/bmi/med/uOttawa/metalab/icon/METALAB_Icon_48.png")));
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
        Label subtitleLabel = new Label("Comprehensive Workflow Manager for Proteomics");
        subtitleLabel.setStyle("-fx-font-size: 18px;");
        
        // Create a container for the title and subtitle
        VBox titleBox = new VBox(5, titleLabel, subtitleLabel);
        
        // Try to load and add the logo
        try {
            ImageView logoView = new ImageView(new Image(getClass().getResourceAsStream("/bmi/med/uOttawa/metalab/icon/METALAB_Logo_Small.png")));
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
        
        HBox workflowBox = new HBox(20);
        workflowBox.setPadding(new Insets(15));
        workflowBox.setAlignment(Pos.CENTER);
        
        workflowBox.getChildren().addAll(
            createWorkflowBox("DDA IGC", "Run DDA workflow for IGC database", 
                e -> launchWorkflow(WorkflowType.DDA_IGC)),
            createWorkflowBox("DDA MAG", "Run DDA workflow for MAG database", 
                e -> launchWorkflow(WorkflowType.DDA_MAG))
        );
        
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
            createWorkflowBox("DIA Normal", "Run DIA workflow in normal mode", 
                e -> launchWorkflow(WorkflowType.DIA_NORMAL)),
            createWorkflowBox("DIA Target", "Run DIA workflow in target mode", 
                e -> launchWorkflow(WorkflowType.DIA_TARGET)),
            createWorkflowBox("DIA Real-time", "Run DIA workflow in real-time mode", 
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
    
    private enum WorkflowType {
        DDA_IGC, DDA_MAG, DIA_NORMAL, DIA_TARGET, DIA_REALTIME
    }
    
    private boolean checkRequiredResources(WorkflowType type) {
        switch (type) {
            case DDA_IGC:
                return resourceParams.isToolAvailable(ToolType.PFIND) && 
                       resourceParams.isToolAvailable(ToolType.DBREDUCER);
                       
            case DDA_MAG:
                return resourceParams.isToolAvailable(ToolType.MAXQUANT) && 
                       resourceParams.isToolAvailable(ToolType.DBREDUCER);
                       
            case DIA_NORMAL:
            case DIA_TARGET:
            case DIA_REALTIME:
                return resourceParams.isToolAvailable(ToolType.FRAGPIPE) && 
                       resourceParams.isToolAvailable(ToolType.DIANN);
                
            default:
                return false;
        }
    }
    
    private void launchWorkflow(WorkflowType type) {
        // First check if all required resources are available
        if (!checkRequiredResources(type)) {
            showWarning("Missing Resources", 
                      "Some required resources are missing for this workflow. " +
                      "Please configure all required resources first.");
            return;
        }
        
        try {
            switch (type) {
                case DDA_IGC:
                    // Launch the pFind main frame
                    MetaPilotPFindMainFrame pFindFrame = new MetaPilotPFindMainFrame();
                    pFindFrame.setVisible(true);
                    break;
                    
                case DDA_MAG:
                    // Launch the MAG main frame
                    // Note: Using default constructor if parameters are not supported
                    MetaPilotMainFrameMag magFrame = new MetaPilotMainFrameMag();
                    magFrame.setVisible(true);
                    break;
                    
                case DIA_NORMAL:
                    // Launch the DIA Normal main frame
                    // Note: Using default constructor if parameters are not supported
                    MetaPilotMainFrameDia diaFrame = new MetaPilotMainFrameDia();
                    diaFrame.setVisible(true);
                    break;
                    
                case DIA_TARGET:
                    // Launch the DIA Target main frame
                    // Note: Using default constructor if parameters are not supported
                    MetaPilotMainFrameDiaTg diaTgFrame = new MetaPilotMainFrameDiaTg();
                    diaTgFrame.setVisible(true);
                    break;
                    
                case DIA_REALTIME:
                    // Launch the DIA Real-time main frame
                    // Note: Using default constructor if parameters are not supported
                    MetaPilotMainFrameDiaRT diaRtFrame = new MetaPilotMainFrameDiaRT();
                    diaRtFrame.setVisible(true);
                    break;
            }
        } catch (Exception ex) {
            showError("Workflow Launch Error", 
                    "Failed to launch workflow", ex);
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
        Alert alert = new Alert(Alert.AlertType.ERROR);
        alert.setTitle(title);
        alert.setHeaderText(null);
        alert.setContentText(content + "\n\nError: " + ex.getMessage());
        alert.showAndWait();
    }
    
    public static void main(String[] args) {
        launch(args);
    }
} 