package uk.ac.quadram.metapilot.gui.dia;

import javafx.geometry.Insets;
import javafx.scene.control.*;
import javafx.scene.layout.*;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.beans.property.*;
import javafx.stage.Stage;
import javafx.scene.Scene;
import uk.ac.quadram.metapilot.gui.common.IOPanel;
import uk.ac.quadram.metapilot.gui.common.ConsolePanel;
import uk.ac.quadram.metapilot.gui.dia.library.LibraryGenerationPanel;
import uk.ac.quadram.metapilot.gui.dia.model.ModelGenerationPanel;

public class MetaPilotDiaRealTimePanel extends VBox {
    private IOPanel ioPanel;
    private ConsolePanel consolePanel;
    private LibraryGenerationPanel libraryPanel;
    private ModelGenerationPanel modelPanel;
    private ComboBox<String> threadCountCombo;
    private ToggleGroup modeGroup;
    private TextArea libraryConsole;
    private TextField updateIntervalField;
    private CheckBox autoUpdateCheck;

    public MetaPilotDiaRealTimePanel() {
        setPadding(new Insets(10));
        setSpacing(10);

        // Create tab pane
        TabPane tabPane = new TabPane();

        // IO Tab
        Tab ioTab = new Tab("Input/Output");
        ioPanel = new IOPanel();
        ioTab.setContent(ioPanel);
        ioTab.setClosable(false);

        // Parameters Tab
        Tab paramsTab = new Tab("Parameters");
        paramsTab.setContent(createParametersPanel());
        paramsTab.setClosable(false);

        // Console Tab
        Tab consoleTab = new Tab("Console");
        consolePanel = new ConsolePanel();
        consoleTab.setContent(consolePanel);
        consoleTab.setClosable(false);

        tabPane.getTabs().addAll(ioTab, paramsTab, consoleTab);
        getChildren().add(tabPane);
    }

    private VBox createParametersPanel() {
        VBox panel = new VBox(10);
        panel.setPadding(new Insets(10));

        // Thread count selection
        HBox threadBox = new HBox(10);
        Label threadLabel = new Label("Thread Count:");
        ObservableList<String> threadOptions = FXCollections.observableArrayList(
            "1", "2", "4", "8", "16", "32", "64"
        );
        threadCountCombo = new ComboBox<>(threadOptions);
        threadCountCombo.setValue("4");
        threadBox.getChildren().addAll(threadLabel, threadCountCombo);

        // Mode selection
        VBox modeBox = new VBox(5);
        Label modeLabel = new Label("Mode:");
        modeGroup = new ToggleGroup();
        RadioButton directButton = new RadioButton("Direct Identification (No Iterative Search)");
        RadioButton iterativeButton = new RadioButton("Iterative Search");
        RadioButton deepLearningButton = new RadioButton("Deep Learning Mode");
        directButton.setToggleGroup(modeGroup);
        iterativeButton.setToggleGroup(modeGroup);
        deepLearningButton.setToggleGroup(modeGroup);
        directButton.setSelected(true);
        modeBox.getChildren().addAll(modeLabel, directButton, iterativeButton, deepLearningButton);

        // Real-time specific settings
        TitledPane realTimePane = new TitledPane("Real-time Settings", createRealTimeSettings());
        realTimePane.setExpanded(true);

        // Library Generation Panel
        TitledPane libraryPane = new TitledPane("Library Generation", createLibraryPanel());
        libraryPane.setExpanded(true);

        // Model Generation Panel
        TitledPane modelPane = new TitledPane("Model Generation", createModelPanel());
        modelPane.setExpanded(true);

        panel.getChildren().addAll(threadBox, modeBox, realTimePane, libraryPane, modelPane);
        return panel;
    }

    private VBox createRealTimeSettings() {
        VBox panel = new VBox(10);
        panel.setPadding(new Insets(5));

        // Update interval
        HBox intervalBox = new HBox(10);
        Label intervalLabel = new Label("Update Interval (seconds):");
        updateIntervalField = new TextField("60");
        updateIntervalField.setPrefWidth(100);
        intervalBox.getChildren().addAll(intervalLabel, updateIntervalField);

        // Auto update option
        autoUpdateCheck = new CheckBox("Enable Auto Update");
        autoUpdateCheck.setSelected(true);

        // Real-time console
        TitledPane consolePane = new TitledPane("Real-time Console", createLibraryConsole());
        consolePane.setExpanded(true);

        panel.getChildren().addAll(intervalBox, autoUpdateCheck, consolePane);
        return panel;
    }

    private VBox createLibraryPanel() {
        VBox panel = new VBox(10);
        panel.setPadding(new Insets(5));

        // Taxa selection
        HBox taxaBox = new HBox(10);
        Label taxaLabel = new Label("Taxa:");
        ComboBox<String> taxaCombo = new ComboBox<>();
        taxaCombo.setPromptText("Select taxa");
        taxaBox.getChildren().addAll(taxaLabel, taxaCombo);

        // Function selection
        HBox functionBox = new HBox(10);
        Label functionLabel = new Label("Function:");
        ComboBox<String> functionCombo = new ComboBox<>();
        functionCombo.setPromptText("Select function");
        functionBox.getChildren().addAll(functionLabel, functionCombo);

        // Library generation console
        TitledPane consolePane = new TitledPane("Library Generation Console", createLibraryConsole());
        consolePane.setExpanded(true);

        panel.getChildren().addAll(taxaBox, functionBox, consolePane);
        return panel;
    }

    private TextArea createLibraryConsole() {
        libraryConsole = new TextArea();
        libraryConsole.setEditable(false);
        libraryConsole.setWrapText(true);
        libraryConsole.setPrefRowCount(10);
        libraryConsole.setStyle("-fx-font-family: monospace;");
        return libraryConsole;
    }

    private VBox createModelPanel() {
        VBox panel = new VBox(10);
        panel.setPadding(new Insets(5));

        // Model generation options
        CheckBox useModelCheck = new CheckBox("Use Model in Data Processing");
        useModelCheck.setSelected(false);

        // Model selection
        HBox modelBox = new HBox(10);
        Label modelLabel = new Label("Model:");
        ComboBox<String> modelCombo = new ComboBox<>();
        modelCombo.setPromptText("Select model");
        modelBox.getChildren().addAll(modelLabel, modelCombo);

        panel.getChildren().addAll(useModelCheck, modelBox);
        return panel;
    }
} 