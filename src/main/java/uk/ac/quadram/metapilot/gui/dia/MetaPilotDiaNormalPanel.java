package uk.ac.quadram.metapilot.gui.dia;

import javafx.geometry.Insets;
import javafx.scene.control.*;
import javafx.scene.layout.*;
import javafx.stage.Stage;
import javafx.scene.Scene;
import uk.ac.quadram.metapilot.gui.common.IOPanel;
import uk.ac.quadram.metapilot.gui.common.ConsolePanel;

public class MetaPilotDiaNormalPanel extends VBox {
    private IOPanel ioPanel;
    private ConsolePanel consolePanel;
    private TabPane tabPane;
    private Tab parametersTab;
    private Tab runningTab;

    public MetaPilotDiaNormalPanel() {
        setPadding(new Insets(20));
        setSpacing(20);
        setStyle("-fx-background-color: white;");

        // Create main layout with three sections
        HBox mainLayout = new HBox(20);
        
        // Left panel - IO Panel
        VBox leftPanel = new VBox(10);
        leftPanel.setPrefWidth(400);
        leftPanel.setStyle("-fx-background-color: #f8f9fa; -fx-padding: 15; -fx-border-color: #e0e0e0; -fx-border-radius: 5;");
        
        Label ioTitle = new Label("Input/Output");
        ioTitle.setStyle("-fx-font-size: 16px; -fx-font-weight: bold;");
        ioPanel = new IOPanel();
        leftPanel.getChildren().addAll(ioTitle, ioPanel);

        // Center panel - Parameters
        VBox centerPanel = new VBox(10);
        centerPanel.setStyle("-fx-background-color: #f8f9fa; -fx-padding: 15; -fx-border-color: #e0e0e0; -fx-border-radius: 5;");
        
        Label paramsTitle = new Label("Parameters");
        paramsTitle.setStyle("-fx-font-size: 16px; -fx-font-weight: bold;");
        
        // Create parameters content
        VBox paramsContent = createParametersPanel();
        centerPanel.getChildren().addAll(paramsTitle, paramsContent);

        // Right panel - Running Status
        VBox rightPanel = new VBox(10);
        rightPanel.setStyle("-fx-background-color: #f8f9fa; -fx-padding: 15; -fx-border-color: #e0e0e0; -fx-border-radius: 5;");
        
        Label runningTitle = new Label("Running Status");
        runningTitle.setStyle("-fx-font-size: 16px; -fx-font-weight: bold;");
        
        // Create running status content
        VBox runningContent = createRunningPanel();
        rightPanel.getChildren().addAll(runningTitle, runningContent);

        // Add all panels to main layout
        mainLayout.getChildren().addAll(leftPanel, centerPanel, rightPanel);
        
        // Add main layout to this panel
        getChildren().add(mainLayout);
    }

    private VBox createParametersPanel() {
        VBox panel = new VBox(15);
        panel.setPadding(new Insets(10));

        // Thread count selection
        HBox threadBox = new HBox(10);
        Label threadLabel = new Label("Thread Count:");
        ComboBox<String> threadCombo = new ComboBox<>();
        threadCombo.getItems().addAll("1", "2", "4", "8", "16", "32", "64");
        threadCombo.setValue("4");
        threadBox.getChildren().addAll(threadLabel, threadCombo);

        // Mode selection
        VBox modeBox = new VBox(10);
        Label modeLabel = new Label("Mode:");
        ToggleGroup modeGroup = new ToggleGroup();
        RadioButton directButton = new RadioButton("Direct Identification");
        RadioButton iterativeButton = new RadioButton("Iterative Search");
        RadioButton deepLearningButton = new RadioButton("Deep Learning");
        directButton.setToggleGroup(modeGroup);
        iterativeButton.setToggleGroup(modeGroup);
        deepLearningButton.setToggleGroup(modeGroup);
        directButton.setSelected(true);
        modeBox.getChildren().addAll(modeLabel, directButton, iterativeButton, deepLearningButton);

        // Advanced parameters
        TitledPane advancedPane = new TitledPane("Advanced Parameters", createAdvancedParameters());
        advancedPane.setExpanded(false);

        panel.getChildren().addAll(threadBox, modeBox, advancedPane);
        return panel;
    }

    private VBox createAdvancedParameters() {
        VBox panel = new VBox(10);
        panel.setPadding(new Insets(10));

        // FDR threshold
        HBox fdrBox = new HBox(10);
        Label fdrLabel = new Label("FDR Threshold:");
        TextField fdrField = new TextField("0.01");
        fdrBox.getChildren().addAll(fdrLabel, fdrField);

        // Minimum peptide length
        HBox peptideBox = new HBox(10);
        Label peptideLabel = new Label("Min Peptide Length:");
        TextField peptideField = new TextField("7");
        peptideBox.getChildren().addAll(peptideLabel, peptideField);

        // Maximum peptide length
        HBox maxPeptideBox = new HBox(10);
        Label maxPeptideLabel = new Label("Max Peptide Length:");
        TextField maxPeptideField = new TextField("30");
        maxPeptideBox.getChildren().addAll(maxPeptideLabel, maxPeptideField);

        panel.getChildren().addAll(fdrBox, peptideBox, maxPeptideBox);
        return panel;
    }

    private VBox createRunningPanel() {
        VBox panel = new VBox(15);
        panel.setPadding(new Insets(10));

        // Progress bar
        ProgressBar progressBar = new ProgressBar(0);
        progressBar.setPrefWidth(Double.MAX_VALUE);

        // Status label
        Label statusLabel = new Label("Ready");
        statusLabel.setStyle("-fx-font-weight: bold;");

        // Console output
        TextArea consoleArea = new TextArea();
        consoleArea.setEditable(false);
        consoleArea.setPrefRowCount(10);
        consoleArea.setStyle("-fx-font-family: monospace;");

        // Control buttons
        HBox buttonBox = new HBox(10);
        Button startButton = new Button("Start");
        Button stopButton = new Button("Stop");
        Button resetButton = new Button("Reset");
        buttonBox.getChildren().addAll(startButton, stopButton, resetButton);

        panel.getChildren().addAll(progressBar, statusLabel, consoleArea, buttonBox);
        return panel;
    }
} 