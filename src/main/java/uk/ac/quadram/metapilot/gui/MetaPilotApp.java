package uk.ac.quadram.metapilot.gui;

import javafx.application.Application;
import javafx.geometry.Insets;
import javafx.scene.Scene;
import javafx.scene.control.*;
import javafx.scene.layout.*;
import javafx.stage.Stage;
import javafx.scene.text.Font;
import javafx.scene.text.FontWeight;
import uk.ac.quadram.metapilot.gui.dda.igc.MetaPilotDdaIgcMaxQuantPanel;
import uk.ac.quadram.metapilot.gui.dda.igc.MetaPilotDdaIgcFragPipePanel;
import uk.ac.quadram.metapilot.gui.dda.igc.MetaPilotDdaIgcPFindPanel;
import uk.ac.quadram.metapilot.gui.dda.igc.MetaPilotDdaIgcSagePanel;
import uk.ac.quadram.metapilot.gui.dda.mag.MetaPilotDdaMagPFindPanel;
import uk.ac.quadram.metapilot.gui.dda.mag.MetaPilotDdaMagFragPipePanel;
import uk.ac.quadram.metapilot.gui.dda.mag.MetaPilotDdaMagSagePanel;
import uk.ac.quadram.metapilot.gui.dia.MetaPilotDiaNormalPanel;
import uk.ac.quadram.metapilot.gui.dia.MetaPilotDiaTargetPanel;
import uk.ac.quadram.metapilot.gui.dia.MetaPilotDiaRealTimePanel;

public class MetaPilotApp extends Application {
    private VBox mainLayout;
    private TabPane tabPane;
    private Tab resourcesTab;
    private Tab databaseTab;
    private Tab workflowsTab;

    @Override
    public void start(Stage primaryStage) {
        primaryStage.setTitle("MetaPilot");
        
        // Create main layout
        mainLayout = new VBox(10);
        mainLayout.setPadding(new Insets(20));
        mainLayout.setStyle("-fx-background-color: #f5f5f5;");
        
        // Create tab pane
        tabPane = new TabPane();
        tabPane.setStyle("-fx-background-color: white; -fx-border-color: #e0e0e0;");
        
        // Create tabs
        workflowsTab = new Tab("Workflows");
        resourcesTab = new Tab("Resources");
        databaseTab = new Tab("Database Configuration");
        
        // Set up workflows tab first
        workflowsTab.setContent(createWorkflowsPanel());
        
        // Set up resources tab
        resourcesTab.setContent(createResourcesPanel());
        
        // Set up database tab
        databaseTab.setContent(createDatabasePanel());
        
        // Add tabs to tab pane in the desired order
        tabPane.getTabs().addAll(workflowsTab, resourcesTab, databaseTab);
        
        // Add tab pane to main layout
        mainLayout.getChildren().add(tabPane);
        
        // Create scene with larger size
        Scene scene = new Scene(mainLayout, 1600, 900);
        scene.getStylesheets().add(getClass().getResource("/styles.css").toExternalForm());
        primaryStage.setScene(scene);
        primaryStage.show();
    }

    private VBox createResourcesPanel() {
        VBox panel = new VBox(10);
        panel.setPadding(new Insets(10));
        
        Label titleLabel = new Label("Resource Configuration");
        titleLabel.setFont(Font.font("System", FontWeight.BOLD, 14));
        
        // Add resource path configuration
        HBox resourcePathBox = new HBox(10);
        Label resourcePathLabel = new Label("Resource Path:");
        TextField resourcePathField = new TextField();
        Button browseButton = new Button("Browse");
        resourcePathBox.getChildren().addAll(resourcePathLabel, resourcePathField, browseButton);
        
        // Add tool paths configuration
        VBox toolsBox = new VBox(5);
        String[] tools = {
            "pFind", "DBReducer", "FragPipe", "Sage", "FlashLFQ", 
            "Function Definition", "MSConvert", "Dia-NN", "DeepDetect", "MaxQuant"
        };
        
        for (String tool : tools) {
            HBox toolBox = new HBox(10);
            Label toolLabel = new Label(tool + ":");
            TextField toolPathField = new TextField();
            Button toolBrowseButton = new Button("Browse");
            toolBox.getChildren().addAll(toolLabel, toolPathField, toolBrowseButton);
            toolsBox.getChildren().add(toolBox);
        }
        
        panel.getChildren().addAll(titleLabel, resourcePathBox, toolsBox);
        return panel;
    }

    private VBox createDatabasePanel() {
        VBox panel = new VBox(20);
        panel.setPadding(new Insets(20));
        panel.setStyle("-fx-background-color: white; -fx-border-color: #e0e0e0; -fx-border-radius: 5;");
        
        Label titleLabel = new Label("Database Configuration");
        titleLabel.setFont(Font.font("System", FontWeight.BOLD, 16));
        
        // Database Type Selection
        TitledPane dbTypePane = new TitledPane("Database Type", createDatabaseTypePanel());
        dbTypePane.setExpanded(true);
        
        // Database Connection Settings
        TitledPane connectionPane = new TitledPane("Connection Settings", createConnectionPanel());
        connectionPane.setExpanded(true);
        
        // Database Schema Settings
        TitledPane schemaPane = new TitledPane("Schema Settings", createSchemaPanel());
        schemaPane.setExpanded(true);
        
        panel.getChildren().addAll(titleLabel, dbTypePane, connectionPane, schemaPane);
        return panel;
    }

    private VBox createDatabaseTypePanel() {
        VBox panel = new VBox(10);
        panel.setPadding(new Insets(10));
        
        ToggleGroup dbTypeGroup = new ToggleGroup();
        RadioButton mysqlRadio = new RadioButton("MySQL");
        RadioButton postgresRadio = new RadioButton("PostgreSQL");
        RadioButton sqliteRadio = new RadioButton("SQLite");
        
        mysqlRadio.setToggleGroup(dbTypeGroup);
        postgresRadio.setToggleGroup(dbTypeGroup);
        sqliteRadio.setToggleGroup(dbTypeGroup);
        
        mysqlRadio.setSelected(true);
        
        panel.getChildren().addAll(mysqlRadio, postgresRadio, sqliteRadio);
        return panel;
    }

    private VBox createConnectionPanel() {
        VBox panel = new VBox(10);
        panel.setPadding(new Insets(10));
        
        GridPane grid = new GridPane();
        grid.setHgap(10);
        grid.setVgap(10);
        
        // Host
        grid.add(new Label("Host:"), 0, 0);
        TextField hostField = new TextField();
        hostField.setPromptText("localhost");
        grid.add(hostField, 1, 0);
        
        // Port
        grid.add(new Label("Port:"), 0, 1);
        TextField portField = new TextField();
        portField.setPromptText("3306");
        grid.add(portField, 1, 1);
        
        // Database
        grid.add(new Label("Database:"), 0, 2);
        TextField dbField = new TextField();
        grid.add(dbField, 1, 2);
        
        // Username
        grid.add(new Label("Username:"), 0, 3);
        TextField userField = new TextField();
        grid.add(userField, 1, 3);
        
        // Password
        grid.add(new Label("Password:"), 0, 4);
        PasswordField passField = new PasswordField();
        grid.add(passField, 1, 4);
        
        // Test Connection Button
        Button testButton = new Button("Test Connection");
        testButton.setStyle("-fx-background-color: #4CAF50; -fx-text-fill: white;");
        grid.add(testButton, 1, 5);
        
        panel.getChildren().add(grid);
        return panel;
    }

    private VBox createSchemaPanel() {
        VBox panel = new VBox(10);
        panel.setPadding(new Insets(10));
        
        // Schema Selection
        HBox schemaBox = new HBox(10);
        Label schemaLabel = new Label("Schema:");
        ComboBox<String> schemaCombo = new ComboBox<>();
        schemaCombo.setPromptText("Select schema");
        schemaBox.getChildren().addAll(schemaLabel, schemaCombo);
        
        // Table Selection
        HBox tableBox = new HBox(10);
        Label tableLabel = new Label("Table:");
        ComboBox<String> tableCombo = new ComboBox<>();
        tableCombo.setPromptText("Select table");
        tableBox.getChildren().addAll(tableLabel, tableCombo);
        
        // Schema Preview
        TextArea previewArea = new TextArea();
        previewArea.setEditable(false);
        previewArea.setPrefRowCount(10);
        previewArea.setStyle("-fx-font-family: monospace;");
        
        panel.getChildren().addAll(schemaBox, tableBox, previewArea);
        return panel;
    }

    private VBox createWorkflowsPanel() {
        VBox panel = new VBox(20);
        panel.setPadding(new Insets(20));
        panel.setStyle("-fx-background-color: white;");
        panel.setMaxHeight(600); // Set maximum height
        
        Label titleLabel = new Label("Workflow Selection");
        titleLabel.setFont(Font.font("System", FontWeight.BOLD, 16));
        
        // DDA IGC Workflows
        TitledPane ddaIgcPane = new TitledPane("DDA IGC", createDDAIGCWorkflows());
        ddaIgcPane.setExpanded(true);
        
        // DDA MAG Workflows
        TitledPane ddaMagPane = new TitledPane("DDA MAG", createDDAMAGWorkflows());
        ddaMagPane.setExpanded(true);
        
        // DIA Workflows
        TitledPane diaPane = new TitledPane("DIA", createDIAWorkflows());
        diaPane.setExpanded(true);
        
        // Create a ScrollPane to contain the workflow panels
        ScrollPane scrollPane = new ScrollPane();
        scrollPane.setFitToWidth(true);
        scrollPane.setStyle("-fx-background: white; -fx-border-color: white;");
        
        VBox content = new VBox(20);
        content.setPadding(new Insets(10));
        content.getChildren().addAll(ddaIgcPane, ddaMagPane, diaPane);
        
        scrollPane.setContent(content);
        
        panel.getChildren().addAll(titleLabel, scrollPane);
        return panel;
    }

    private VBox createDDAIGCWorkflows() {
        VBox panel = new VBox(10);
        panel.setPadding(new Insets(10));
        
        FlowPane flowPane = new FlowPane();
        flowPane.setHgap(10);
        flowPane.setVgap(10);
        flowPane.setPrefWrapLength(800); // Set preferred width for wrapping
        
        Button maxQuantButton = createWorkflowButton("MaxQuant", () -> {
            Stage stage = new Stage();
            stage.setTitle("DDA IGC MaxQuant");
            stage.setScene(new Scene(new MetaPilotDdaIgcMaxQuantPanel(), 1600, 900));
            stage.show();
        });
        
        Button pFindButton = createWorkflowButton("pFind", () -> {
            Stage stage = new Stage();
            stage.setTitle("DDA IGC PFind");
            stage.setScene(new Scene(new MetaPilotDdaIgcPFindPanel(), 1600, 900));
            stage.show();
        });
        
        Button fragPipeButton = createWorkflowButton("FragPipe", () -> {
            Stage stage = new Stage();
            stage.setTitle("DDA IGC FragPipe");
            stage.setScene(new Scene(new MetaPilotDdaIgcFragPipePanel(), 1600, 900));
            stage.show();
        });
        
        Button sageButton = createWorkflowButton("Sage", () -> {
            Stage stage = new Stage();
            stage.setTitle("DDA IGC Sage");
            stage.setScene(new Scene(new MetaPilotDdaIgcSagePanel(), 1600, 900));
            stage.show();
        });
        
        flowPane.getChildren().addAll(maxQuantButton, pFindButton, fragPipeButton, sageButton);
        panel.getChildren().add(flowPane);
        return panel;
    }

    private VBox createDDAMAGWorkflows() {
        VBox panel = new VBox(10);
        panel.setPadding(new Insets(10));
        
        FlowPane flowPane = new FlowPane();
        flowPane.setHgap(10);
        flowPane.setVgap(10);
        flowPane.setPrefWrapLength(800);
        
        Button pFindButton = createWorkflowButton("PFind", () -> {
            Stage stage = new Stage();
            stage.setTitle("DDA MAG PFind");
            stage.setScene(new Scene(new MetaPilotDdaMagPFindPanel(), 1600, 900));
            stage.show();
        });
        
        Button fragPipeButton = createWorkflowButton("FragPipe", () -> {
            Stage stage = new Stage();
            stage.setTitle("DDA MAG FragPipe");
            stage.setScene(new Scene(new MetaPilotDdaMagFragPipePanel(), 1600, 900));
            stage.show();
        });
        
        Button sageButton = createWorkflowButton("Sage", () -> {
            Stage stage = new Stage();
            stage.setTitle("DDA MAG Sage");
            stage.setScene(new Scene(new MetaPilotDdaMagSagePanel(), 1600, 900));
            stage.show();
        });
        
        flowPane.getChildren().addAll(pFindButton, fragPipeButton, sageButton);
        panel.getChildren().add(flowPane);
        return panel;
    }

    private VBox createDIAWorkflows() {
        VBox panel = new VBox(10);
        panel.setPadding(new Insets(10));
        
        FlowPane flowPane = new FlowPane();
        flowPane.setHgap(10);
        flowPane.setVgap(10);
        flowPane.setPrefWrapLength(800);
        
        Button normalButton = createWorkflowButton("Normal", () -> {
            Stage stage = new Stage();
            stage.setTitle("DIA Normal Workflow");
            stage.setScene(new Scene(new MetaPilotDiaNormalPanel(), 1600, 900));
            stage.show();
        });
        
        Button targetButton = createWorkflowButton("Target", () -> {
            Stage stage = new Stage();
            stage.setTitle("DIA Target Workflow");
            stage.setScene(new Scene(new MetaPilotDiaTargetPanel(), 1600, 900));
            stage.show();
        });
        
        Button realTimeButton = createWorkflowButton("Real-time", () -> {
            Stage stage = new Stage();
            stage.setTitle("DIA Real-time Workflow");
            stage.setScene(new Scene(new MetaPilotDiaRealTimePanel(), 1600, 900));
            stage.show();
        });
        
        flowPane.getChildren().addAll(normalButton, targetButton, realTimeButton);
        panel.getChildren().add(flowPane);
        return panel;
    }

    private Button createWorkflowButton(String text, Runnable action) {
        Button button = new Button(text);
        button.setPrefSize(200, 100);
        button.setStyle("-fx-background-color: #f8f9fa; " +
                       "-fx-border-color: #e0e0e0; " +
                       "-fx-border-width: 1; " +
                       "-fx-font-size: 14px; " +
                       "-fx-text-fill: black; " +
                       "-fx-alignment: center;");
        button.setOnAction(e -> action.run());
        return button;
    }

    public static void main(String[] args) {
        launch(args);
    }
} 