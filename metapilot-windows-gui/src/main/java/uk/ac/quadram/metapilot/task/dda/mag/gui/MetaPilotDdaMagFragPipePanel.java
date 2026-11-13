package uk.ac.quadram.metapilot.task.dda.mag.gui;

import uk.ac.quadram.metapilot.core.fragpipe.FragpipeModifications;
import uk.ac.quadram.metapilot.core.msfragger.MSFraggerHeadlessRunner;
import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.control.*;
import javafx.scene.layout.*;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.scene.text.Font;
import javafx.scene.text.FontWeight;
import javafx.scene.paint.Color;
import javafx.scene.control.Alert.AlertType;
import java.util.Optional;
import uk.ac.quadram.metapilot.task.common.gui.MetaPilotIOPanelFX;
import uk.ac.quadram.metapilot.task.common.gui.MetaPilotIOPanelFX.ExperimentRow;
import uk.ac.quadram.metapilot.task.common.gui.MagDatabaseSettingsPanel;
import uk.ac.quadram.metapilot.task.common.gui.ConsolePanel;
import uk.ac.quadram.metapilot.task.common.gui.MetadataConfigurationPanel;
import uk.ac.quadram.metapilot.task.common.gui.ModificationsPanel;
import uk.ac.quadram.metapilot.task.common.gui.ParameterSection;
import uk.ac.quadram.metapilot.task.dda.mag.par.DDAMagFragPipeParameters;
import uk.ac.quadram.metapilot.task.dda.mag.task.DDAMagFragPipeTask;
import uk.ac.quadram.metapilot.task.dda.mag.par.DDAMagFragPipeParameterIO;
import uk.ac.quadram.metapilot.task.mag.GenomeSelectionCut;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;
import uk.ac.quadram.metapilot.task.mag.MagDbItem.MagVersion;
import uk.ac.quadram.metapilot.task.par.MetaData;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.task.common.gui.AbstractWorkflowPanel;
import javafx.application.Platform;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.IOException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.List;
import java.util.ArrayList;
import java.io.File;
import javafx.stage.FileChooser;

public class MetaPilotDdaMagFragPipePanel extends AbstractWorkflowPanel<DDAMagFragPipeParameters> {
    // UI Components - Database specific for MAG
    private MetaPilotIOPanelFX ioPanel;
    private MagDatabaseSettingsPanel magDatabasePanel;

    // Reusable UI Components
    private MetadataConfigurationPanel metadataPanel;

    // Parameter Panel Components
    private ComboBox<String> enzymeComboBox;
    private Spinner<Integer> missedCleavagesSpinner;
    private ComboBox<String> digestionModeComboBox;

    // Reusable modifications panel
    private ModificationsPanel modificationsPanel;

    // Additional parameters
    private ComboBox<Integer> threadsComboBox;
    // Power-user escape hatch: explicit .workflow file overrides, blank by default.
    // Open search belongs to the dedicated DDA MAG open-search workflow; this panel
    // always runs FragPipe closed (LFQ-MBR) unless a custom .workflow is supplied.
    private TextField customBasicWorkflowField;
    private Button browseCustomBasicWorkflowButton;
    private TextField customCombinedWorkflowField;
    private Button browseCustomCombinedWorkflowButton;
    private ComboBox<String> msDataTypeComboBox;
    private ComboBox<String> genomeSelectorComboBox;
    private GenomeStopControls genomeStopControls;

    // Data
    private ObservableList<String> enzymes = FXCollections.observableArrayList();

    private MetaResourceParameters resourceParams;
    
    // Store FragpipeWorkflow to pass to task (not stored in parameters)
    private uk.ac.quadram.metapilot.dbSearch.fragpipe.FragpipeWorkflow currentFragpipeWorkflow;
    
    
    
    // Logger
    private static final Logger LOGGER = LogManager.getLogger(MetaPilotDdaMagFragPipePanel.class);

	public MetaPilotDdaMagFragPipePanel(MagDbItem[] magDbItems, MetaResourceParameters resourceParams) {
		// Call parent constructor first
		super();
		
		// Initialize parameter IO instance
		parameterIO = new DDAMagFragPipeParameterIO();
		this.resourceParams = resourceParams;
		
		loadEnzymeChoices();
		initializeSpecificUI();
		// Set MAG databases after UI is initialized
		setAvailableMagDatabases(magDbItems);
	}
	
    // Load available MAG databases
    public void setAvailableMagDatabases(MagDbItem[] magDbItems) {
        if (magDatabasePanel != null) {
            magDatabasePanel.setAvailableMagDatabases(magDbItems);
        }
    }

    private void initializeSpecificUI() {
        setSpacing(10);
        setPadding(new Insets(10));
        setMinWidth(1200);
        setMinHeight(800);
        setPrefWidth(1600);
        setPrefHeight(900);

        // Initialize reusable components
        consolePanel = new ConsolePanel();
        metadataPanel = new MetadataConfigurationPanel(true);
        magDatabasePanel = new MagDatabaseSettingsPanel();
        
        // Set up callbacks for parameter updates
        metadataPanel.setParameterUpdateCallback(this::updateParameterView);
        magDatabasePanel.setParameterUpdateCallback(this::updateParameterView);
        configureConsolePanel(consolePanel);

        // Create tab pane with responsive sizing
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

        // Add tab change listener to update metadata when Parameters tab is selected
        // and save parameters when Console tab is selected
        tabPane.getSelectionModel().selectedItemProperty().addListener((obs, oldTab, newTab) -> {
            if (newTab == paramTab) {
                updateMetadataTableFromIOPanel();
            } else if (newTab == consoleTab) {
                // Update parameter view first to show current modifications
                updateParameterView();
                
                // Save current parameters when Console tab is selected
                try {
                    consolePanel.appendConsoleText("Auto-saving parameters when switching to Console tab...\n");
                    DDAMagFragPipeParameters currentParams = collectParametersFromUI();
                    if (currentParams != null) {
                        boolean saved = parameterIO.saveLastUsedParameters(currentParams);
                        if (saved) {
                            consolePanel.appendConsoleText("✓ Parameters auto-saved to config folder\n");
                            LOGGER.info("Successfully auto-saved parameters to config folder when Console tab was selected");
                        } else {
                            consolePanel.appendConsoleText("WARNING: Failed to auto-save parameters to config folder\n");
                            LOGGER.warn("Failed to auto-save parameters to config folder when Console tab was selected");
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
                        LOGGER.warn("No parameters to auto-save when Console tab was selected");
                    }
                } catch (Exception e) {
                    String errorMsg = "ERROR: Failed to auto-save parameters: " + e.getMessage();
                    consolePanel.appendConsoleText(errorMsg + "\n");
                    LOGGER.warn("Failed to auto-save parameters when Console tab was selected: " + e.getMessage(), e);
                }
            }
        });

        // Add tab pane to main panel
        getChildren().add(tabPane);
        
        // Load parameters after UI is fully initialized
        Platform.runLater(() -> {
            loadLastUsedParameters();
            
            // Add window close event handler to save parameters
            setupWindowCloseHandler();
        });
    }

    // Create IO tab (similar to pFind)
    @Override
    protected String getWorkflowName() {
        return "DDA MAG FragPipe";
    }

    @Override
    protected boolean isHpcExportSupported() {
        return true;
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
            magDatabasePanel
        );
        return tabContent;
    }


    @Override
    public VBox createParametersTab() {
        VBox tabContent = new VBox(10);
        tabContent.setPadding(new Insets(10));
        tabContent.setMinWidth(960);

        // Create main grid layout for parameters with 2-column responsive sizing
        GridPane mainGrid = new GridPane();
        mainGrid.setHgap(15);
        mainGrid.setVgap(15);
        mainGrid.setPadding(new Insets(10));
        mainGrid.setMinWidth(940);

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
        mainGrid.getRowConstraints().add(row1);

        // Left column: Enzyme Settings, Additional Parameters, and Modifications Panel (NO taxonomic analysis for MAG)
        VBox leftColumn = new VBox(15);
        leftColumn.setMinWidth(430);
        TitledPane enzymeSection = createEnzymeSection();
        VBox additionalParamsSection = createAdditionalParametersSection();
        
        // Create horizontal modifications panel with fixed mods (narrow) and variable mods (wider)
        modificationsPanel = new ModificationsPanel(this::updateParameterView);
        modificationsPanel.setModifications(FragpipeModifications.DEFAULT_FIXED,
                FragpipeModifications.DEFAULT_VARIABLE);
        
        leftColumn.getChildren().addAll(enzymeSection, additionalParamsSection, modificationsPanel);
        GridPane.setVgrow(leftColumn, Priority.ALWAYS);
        
        // Right column: Metadata Configuration only
        VBox rightColumn = new VBox(15);
        rightColumn.setMinWidth(430);
        
        rightColumn.getChildren().addAll(metadataPanel);
        GridPane.setVgrow(rightColumn, Priority.ALWAYS);

        // Add columns to grid
        mainGrid.add(leftColumn, 0, 0);
        mainGrid.add(rightColumn, 1, 0);

        tabContent.getChildren().add(ParameterSection.createViewport(mainGrid));
        return tabContent;
    }

    // Add missing utility methods
    /**
     * Enzyme choices come from {@link MSFraggerHeadlessRunner#ENZYME_CHOICES}, the
     * same list the enzyme mapping switches on, rather than from
     * {@code enzyme_default.ini}.
     * <p>
     * That file is pFind's configuration vocabulary, in which {@code Trypsin} is
     * declared with no proline restriction and {@code Trypsin_P} with one - the
     * opposite of what those names mean to MSFragger, to MetaPilot's own
     * {@code core.enzyme.Enzyme} enum, and to its DIA-NN mapping. Populating this
     * dropdown from there is how that inversion reached the MSFragger workflows.
     * The pFind panels keep pFind's own names ({@code PFindEnzymes}), where those
     * semantics are the engine's.
     */
    private void loadEnzymeChoices() {
        enzymes.setAll(MSFraggerHeadlessRunner.ENZYME_CHOICES);
    }

    private Label createSectionTitle(String title) {
        Label label = new Label(title);
        label.setFont(Font.font("System", FontWeight.BOLD, 14));
        label.setTextFill(Color.DARKBLUE);
        return label;
    }

    private TitledPane createEnzymeSection() {

        HBox enzymeParamsBox = new HBox(15);
        enzymeParamsBox.setAlignment(Pos.CENTER_LEFT);

        VBox enzymeBox = new VBox(5);
        enzymeBox.setMinWidth(140);
        HBox.setHgrow(enzymeBox, Priority.ALWAYS);
        Label enzymeLabel = new Label("Enzyme:");
        enzymeLabel.setStyle("-fx-font-weight: bold;");
        enzymeBox.getChildren().add(enzymeLabel);

        enzymeComboBox = new ComboBox<>(enzymes);
        enzymeComboBox.setMinWidth(140);
        enzymeComboBox.setPrefWidth(160);
        enzymeComboBox.setMaxWidth(Double.MAX_VALUE);
        if (enzymes.size() > 0) {
            enzymeComboBox.setValue(enzymes.get(0));
        }
        enzymeComboBox.setOnAction(e -> updateParameterView());
        enzymeBox.getChildren().add(enzymeComboBox);

        VBox missedCleavagesBox = new VBox(5);
        missedCleavagesBox.setMinWidth(100);
        Label missedLabel = new Label("Missed Cleavages:");
        missedLabel.setStyle("-fx-font-weight: bold;");
        missedCleavagesSpinner = new Spinner<>(0, 5, 2);
        missedCleavagesSpinner.setMinWidth(100);
        missedCleavagesSpinner.setPrefWidth(100);
        missedCleavagesSpinner.setMaxWidth(Double.MAX_VALUE);
        missedCleavagesSpinner.valueProperty().addListener((obs, oldVal, newVal) -> updateParameterView());
        missedCleavagesBox.getChildren().add(missedLabel);
        missedCleavagesBox.getChildren().add(missedCleavagesSpinner);

        VBox digestionModeBox = new VBox(5);
        digestionModeBox.setMinWidth(140);
        HBox.setHgrow(digestionModeBox, Priority.ALWAYS);
        Label digestionLabel = new Label("Digestion Mode:");
        digestionLabel.setStyle("-fx-font-weight: bold;");
        digestionModeBox.getChildren().add(digestionLabel);

        digestionModeComboBox = new ComboBox<>(
                FXCollections.observableArrayList(MSFraggerHeadlessRunner.DIGEST_CHOICES));
        digestionModeComboBox.setMinWidth(140);
        digestionModeComboBox.setPrefWidth(160);
        digestionModeComboBox.setMaxWidth(Double.MAX_VALUE);
        digestionModeComboBox.setValue(MSFraggerHeadlessRunner.DIGEST_SPECIFIC);
        digestionModeComboBox.setOnAction(e -> updateParameterView());
        digestionModeBox.getChildren().add(digestionModeComboBox);

        enzymeParamsBox.getChildren().addAll(enzymeBox, missedCleavagesBox, digestionModeBox);

        return ParameterSection.create("Enzyme settings", enzymeParamsBox);
    }

    private VBox createAdditionalParametersSection() {
        GridPane paramGrid = new GridPane();
        paramGrid.setHgap(15);
        paramGrid.setVgap(8);
        
        ColumnConstraints col1 = new ColumnConstraints();
        col1.setPercentWidth(50);
        ColumnConstraints col2 = new ColumnConstraints();
        col2.setPercentWidth(50);
        paramGrid.getColumnConstraints().addAll(col1, col2);
        
        Label threadsLabel = new Label("Threads:");
        threadsLabel.setStyle("-fx-font-weight: bold;");
        threadsComboBox = new ComboBox<>();
        threadsComboBox.setMinWidth(100);
        threadsComboBox.setPrefWidth(100);
        threadsComboBox.setMaxWidth(Double.MAX_VALUE);
        
        // Populate with thread count options (1 to available processors)
        ObservableList<Integer> threadOptions = FXCollections.observableArrayList();
        for (int i = 1; i <= Runtime.getRuntime().availableProcessors(); i++) {
            threadOptions.add(i);
        }
        threadsComboBox.setItems(threadOptions);
        threadsComboBox.setValue(4); // Default value
        threadsComboBox.valueProperty().addListener((obs, oldVal, newVal) -> updateParameterView());
        
        VBox threadsBox = new VBox(5);
        threadsBox.getChildren().addAll(threadsLabel, threadsComboBox);
        
        FileChooser.ExtensionFilter workflowFilter =
                new FileChooser.ExtensionFilter("FragPipe workflow files", "*.workflow");
        FileChooser.ExtensionFilter allFilesFilter = new FileChooser.ExtensionFilter("All Files", "*.*");

        customBasicWorkflowField = new TextField();
        customBasicWorkflowField.setPromptText("Select custom marker/sample-search .workflow file");
        customBasicWorkflowField.setPrefWidth(400);
        browseCustomBasicWorkflowButton = new Button("Browse");
        browseCustomBasicWorkflowButton.setOnAction(e -> {
            FileChooser chooser = new FileChooser();
            chooser.setTitle("Select Custom Marker/Sample-Search Workflow File");
            chooser.getExtensionFilters().addAll(workflowFilter, allFilesFilter);
            File file = chooser.showOpenDialog(getScene().getWindow());
            if (file != null) {
                customBasicWorkflowField.setText(file.getAbsolutePath());
                updateParameterView();
            }
        });
        HBox customBasicRow = new HBox(5, new Label("Marker/sample-search:"), customBasicWorkflowField,
                browseCustomBasicWorkflowButton);
        customBasicRow.setAlignment(Pos.CENTER_LEFT);
        HBox.setHgrow(customBasicWorkflowField, Priority.ALWAYS);

        customCombinedWorkflowField = new TextField();
        customCombinedWorkflowField.setPromptText("Select custom combined/quant search .workflow file");
        customCombinedWorkflowField.setPrefWidth(400);
        browseCustomCombinedWorkflowButton = new Button("Browse");
        browseCustomCombinedWorkflowButton.setOnAction(e -> {
            FileChooser chooser = new FileChooser();
            chooser.setTitle("Select Custom Combined/Quant Search Workflow File");
            chooser.getExtensionFilters().addAll(workflowFilter, allFilesFilter);
            File file = chooser.showOpenDialog(getScene().getWindow());
            if (file != null) {
                customCombinedWorkflowField.setText(file.getAbsolutePath());
                updateParameterView();
            }
        });
        HBox customCombinedRow = new HBox(5, new Label("Combined/quant search:"), customCombinedWorkflowField,
                browseCustomCombinedWorkflowButton);
        customCombinedRow.setAlignment(Pos.CENTER_LEFT);
        HBox.setHgrow(customCombinedWorkflowField, Priority.ALWAYS);

        // Add MS data type option
        Label msDataTypeLabel = new Label("MS data type:");
        msDataTypeLabel.setStyle("-fx-font-weight: bold;");
        msDataTypeComboBox = new ComboBox<>(FXCollections.observableArrayList("Regular MS", "IM-MS (ion mobility, timsTOF only)"));
        msDataTypeComboBox.setValue("Regular MS");
        msDataTypeComboBox.setMinWidth(100);
        msDataTypeComboBox.setPrefWidth(200);
        msDataTypeComboBox.setMaxWidth(Double.MAX_VALUE);
        msDataTypeComboBox.setOnAction(e -> updateParameterView());
        
        VBox msDataTypeBox = new VBox(5);
        msDataTypeBox.getChildren().addAll(msDataTypeLabel, msDataTypeComboBox);

        Label genomeSelectorLabel = new Label("Marker set");
        genomeSelectorLabel.setStyle("-fx-font-weight: bold;");
        genomeSelectorComboBox = new ComboBox<>(FXCollections.observableArrayList(MagDbItem.MARKER2, MagDbItem.MARKER4));
        genomeSelectorComboBox.setValue(MagDbItem.MARKER2);
        genomeSelectorComboBox.setMinWidth(140);
        genomeSelectorComboBox.setPrefWidth(160);
        genomeSelectorComboBox.setMaxWidth(Double.MAX_VALUE);
        genomeSelectorComboBox.setPromptText("Genome selection workflow");
        genomeSelectorComboBox.setOnAction(e -> updateParameterView());
        VBox genomeSelectorBox = new VBox(5);
        genomeSelectorBox.getChildren().addAll(genomeSelectorLabel, genomeSelectorComboBox);

        genomeStopControls = new GenomeStopControls(this::updateParameterView);

        paramGrid.add(msDataTypeBox, 0, 0);
        paramGrid.add(threadsBox, 1, 0);
        TitledPane searchSection = ParameterSection.create("Search", paramGrid);

        VBox genomeSettings = new VBox(8, genomeSelectorBox, genomeStopControls);
        TitledPane genomeSection = ParameterSection.create("Genome selection", genomeSettings);

        Label workflowHint = new Label("Optional .workflow files");
        VBox workflowSettings = new VBox(6, workflowHint, customBasicRow, customCombinedRow);
        TitledPane workflowSection = ParameterSection.create("Workflow overrides", workflowSettings);
        return new VBox(10, searchSection, genomeSection, workflowSection);
    }
    
    // Combined modifications panel with fixed mods (narrow) and variable mods (wider) side by side

    private void updateMetadataTableFromIOPanel() {
        if (ioPanel == null || metadataPanel == null) {
            return;
        }
        
        ObservableList<MetaPilotIOPanelFX.ExperimentRow> experiments = ioPanel.getExperimentData();
        metadataPanel.updateMetadataFromIOPanel(experiments);
    }

    @Override
    protected void updateParameterView() {
        if (consolePanel == null) {
            return;
        }
        
        StringBuilder sb = new StringBuilder();
        sb.append("=== FragPipe Analysis Parameters ===\n\n");

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

        // Enzyme settings
        sb.append("Enzyme Settings:\n");
        if (enzymeComboBox != null) {
            sb.append("  Enzyme: ").append(enzymeComboBox.getValue()).append("\n");
        }
        if (missedCleavagesSpinner != null) {
            sb.append("  Missed Cleavages: ").append(missedCleavagesSpinner.getValue()).append("\n");
        }
        if (digestionModeComboBox != null) {
            sb.append("  Digestion Mode: ").append(digestionModeComboBox.getValue()).append("\n");
        }
        sb.append("\n");

        // Fixed modifications
        sb.append("Fixed Modifications:\n");
        if (modificationsPanel != null) {
            for (ModificationsPanel.FixedModificationRow row : modificationsPanel.getFixedModData()) {
                if (row.getDeltaMass() != 0.0) {
                    sb.append("  ").append(row.getSite()).append(": ").append(row.getDeltaMass()).append("\n");
                }
            }
        }
        sb.append("\n");

        // Variable modifications
        sb.append("Variable Modifications:\n");
        if (modificationsPanel != null) {
            for (ModificationsPanel.VariableModificationRow row : modificationsPanel.getVariableModData()) {
                if (row.getUse()) {
                    sb.append("  ").append(row.getModName()).append(" (").append(row.getSites()).append("): ")
                      .append(row.getDeltaMass()).append("\n");
                }
            }
        }
        sb.append("\n");

        // Additional parameters
        sb.append("Additional Parameters:\n");
        if (threadsComboBox != null) {
            sb.append("  Threads: ").append(threadsComboBox.getValue()).append("\n");
        }
        sb.append("  Search mode: Closed (open search uses the dedicated open-search workflow)\n");
        if (msDataTypeComboBox != null) {
            sb.append("  MS data type: ").append(msDataTypeComboBox.getValue()).append("\n");
        }
        if (genomeSelectorComboBox != null) {
            sb.append("  Marker set: ").append(genomeSelectorComboBox.getValue()).append("\n");
        }
        if (genomeStopControls != null) {
            sb.append("  When to stop adding genomes: ").append(genomeStopControls.summary()).append("\n");
        }
        sb.append("\n");

        // Note: Quantitative mode is determined by the selected workflow
        sb.append("Quantitative Mode: Determined by workflow selection\n\n");

        consolePanel.setParameterText(sb.toString());
    }


    private void updateMetadataTable() {
        if (metadataPanel == null || ioPanel == null) {
            return;
        }
        
        ObservableList<MetaPilotIOPanelFX.ExperimentRow> fileRows = ioPanel.getExperimentData();
        if (fileRows.isEmpty()) return;
        
        // Note: Quantitative mode is determined by workflow, use default label-free mode
        ObservableList<uk.ac.quadram.metapilot.task.common.gui.IsobaricLabelingPanel.IsobaricLabelRow> emptyIsobaricList = FXCollections.observableArrayList();
        metadataPanel.updateMetadataForQuantMode("Label-free", false, fileRows, emptyIsobaricList);
    }


    @Override
    protected void startAnalysis() {
        try {
            consolePanel.appendConsoleText("Starting MAG FragPipe analysis...\n");
            System.out.println("DEBUG: startAnalysis called - MAG FragPipe workflow");
            
            // 1. Collect parameters from UI
            consolePanel.appendConsoleText("Collecting parameters from UI...\n");
            DDAMagFragPipeParameters analysisParams = collectParametersFromUI();
            if (analysisParams == null) {
                consolePanel.appendConsoleText("ERROR: Failed to collect parameters from UI.\n");
                LOGGER.error("Failed to collect parameters from UI in startAnalysis");
                return;
            }
            
            consolePanel.appendConsoleText("✓ Successfully collected parameters from UI\n");
            LOGGER.info("Successfully collected parameters from UI for MAG FragPipe workflow");
            
            // Store analysis parameters for viewReport
            this.currentAnalysisParams = analysisParams;
            
            // 2. Save parameters to file
            String outputDir = analysisParams.getResult();
            if (outputDir != null && !outputDir.trim().isEmpty()) {
                consolePanel.appendConsoleText("Saving parameters to output directory...\n");
                
                // Save analysis-specific parameters
                boolean analysisParamsSaved = parameterIO.saveAnalysisParameters(outputDir, analysisParams);
                if (!analysisParamsSaved) {
                    consolePanel.appendConsoleText("WARNING: Failed to save analysis parameters to result directory\n");
                    LOGGER.warn("Failed to save analysis parameters to {}/analysis_parameters.json", outputDir);
                } else {
                    consolePanel.appendConsoleText("✓ Analysis parameters saved to result directory\n");
                    LOGGER.info("Saved analysis parameters to {}/analysis_parameters.json", outputDir);
                }
                
                // Save as last used parameters
                boolean lastParamsSaved = parameterIO.saveLastUsedParameters(analysisParams);
                if (!lastParamsSaved) {
                    consolePanel.appendConsoleText("WARNING: Failed to save last used parameters\n");
                    LOGGER.warn("Failed to save last used parameters");
                } else {
                    consolePanel.appendConsoleText("✓ Last used parameters saved\n");
                }
                
                // Also save to preferences for quick defaults
                boolean prefsSaved = parameterIO.savePreferences(analysisParams);
                if (prefsSaved) {
                    consolePanel.appendConsoleText("✓ Preferences saved\n");
                } else {
                    consolePanel.appendConsoleText("WARNING: Failed to save preferences\n");
                }
            } else {
                consolePanel.appendConsoleText("WARNING: No output directory specified, skipping parameter save\n");
                LOGGER.warn("No output directory specified, skipping parameter save");
            }
            
            // 3. Create MetaResourceParameters (load from config file)
            consolePanel.appendConsoleText("Loading resource parameters from config...\n");
            consolePanel.appendConsoleText("✓ Resource parameters loaded from config file\n");
            
            // 4. Create progress bars and status label
            consolePanel.appendConsoleText("Setting up progress tracking...\n");
            Label statusLabel = new Label("Ready to start analysis...");
            consolePanel.appendConsoleText("✓ Progress bars and status label initialized\n");
            
            // 5. Create and run DDAMagFragPipeTask - pass statusLabel and workflow directly to task (not through parameters)
            consolePanel.appendConsoleText("Creating MAG FragPipe task...\n");
			currentTask = new DDAMagFragPipeTask(analysisParams, resourceParams, consolePanel.getOverallProgressBar(),
					consolePanel.getCurrentTaskProgressBar(), statusLabel, this.currentFragpipeWorkflow);
            consolePanel.appendConsoleText("✓ Task created successfully\n");
            
            // 7. Set up task callbacks
            consolePanel.appendConsoleText("Setting up task callbacks...\n");
            currentTask.messageProperty().addListener((obs, oldMsg, newMsg) -> {
                Platform.runLater(() -> {
                    if (newMsg != null) {
                        consolePanel.appendConsoleText("Current Task: " + newMsg + "\n");
                        consolePanel.updateCurrentTaskFromMessage(newMsg);
                        statusLabel.setText(newMsg);
                    }
                });
            });
            
            // Set up task completion callback to enable View Report button
            currentTask.setOnSucceeded(e -> {
                Platform.runLater(() -> {
                    // Check if the task actually succeeded by examining the result
                    try {
                        Boolean result = currentTask.get();
                        if (result != null && result) {
                            consolePanel.appendConsoleText("✓ Analysis completed successfully!\n");
                            consolePanel.enableViewReportButton();
                        } else {
                            consolePanel.appendConsoleText("⚠ Analysis completed but may not have been fully successful.\n");
                        }
                    } catch (Exception ex) {
                        consolePanel.appendConsoleText("✗ Analysis failed: " + ex.getMessage() + "\n");
                    }
                    consolePanel.setButtonsEnabled(true, false, true);
                    consolePanel.stopStatusTimer(); // Stop the timer when task completes
                });
            });
            
            currentTask.setOnFailed(e -> {
                Platform.runLater(() -> {
                    consolePanel.appendConsoleText("✗ Analysis failed. Check console output for details.\n");
                    consolePanel.setButtonsEnabled(true, false, false);
                    consolePanel.stopStatusTimer(); // Stop the timer when task fails
                });
            });
            
            consolePanel.appendConsoleText("✓ Task callbacks configured\n");
            
            // 8. Run task in background thread
            consolePanel.appendConsoleText("Starting task execution...\n");
            Thread taskThread = new Thread(currentTask);
            taskThread.setDaemon(true);
            taskThread.start();
            
            consolePanel.appendConsoleText("✓ Analysis task started successfully.\n");
            LOGGER.info("MAG FragPipe analysis task started successfully");
            
        } catch (Exception e) {
            String errorMsg = "ERROR: Failed to start analysis: " + e.getMessage();
            consolePanel.appendConsoleText(errorMsg + "\n");
            LOGGER.error("Failed to start MAG FragPipe analysis", e);
            e.printStackTrace();
        }
    }

    @Override
    protected void stopAnalysis() {
        Alert alert = new Alert(AlertType.CONFIRMATION);
        alert.setTitle("Stop Analysis");
        alert.setHeaderText("Are you sure you want to stop the analysis?");
        alert.setContentText("This action cannot be undone.");
        
        Optional<ButtonType> result = alert.showAndWait();
        if (result.isPresent() && result.get() == ButtonType.OK) {
            consolePanel.appendConsoleText("Stopping FragPipe analysis...\n");
            LOGGER.info("stopAnalysis called - user confirmed stop");
            
            // Force stop the running task if it exists
            if (currentTask != null && currentTask.isRunning()) {
                consolePanel.appendConsoleText("Force stopping MAG FragPipe task...\n");
                currentTask.forceStop();
                consolePanel.appendConsoleText("MAG FragPipe task stop requested.\n");
            } else {
                consolePanel.appendConsoleText("No running MAG FragPipe task to stop.\n");
            }
            
            // Unbind progress bars to prevent memory leaks
            consolePanel.getOverallProgressBar().progressProperty().unbind();
            consolePanel.getCurrentTaskProgressBar().progressProperty().unbind();
            
            consolePanel.setButtonsEnabled(true, false, false);
            consolePanel.stopStatusTimer();
            consolePanel.getOverallProgressBar().setProgress(0);
            consolePanel.getCurrentTaskProgressBar().setProgress(0);
            consolePanel.setCurrentTaskLabel("Current Task: None");
        } else {
            LOGGER.info("stopAnalysis cancelled by user");
        }
    }
    
    
    /**
     * Collect all parameters from the UI into a parameters object
     */
    @Override
    protected DDAMagFragPipeParameters collectParametersFromUI() {
        try {
            consolePanel.appendConsoleText("Collecting parameters from UI...\n");
            LOGGER.info("Starting parameter collection from UI for MAG FragPipe workflow");
            
            // Validate inputs first
            consolePanel.appendConsoleText("Validating UI inputs...\n");
            if (!validateInputs()) {
                consolePanel.appendConsoleText("ERROR: Input validation failed\n");
                LOGGER.error("Input validation failed");
                return null;
            }
            consolePanel.appendConsoleText("✓ Input validation passed\n");
            LOGGER.info("Input validation passed");
            
            // Collect experiment data from ioPanel
            consolePanel.appendConsoleText("Collecting experiment data from IO panel...\n");
            ObservableList<MetaPilotIOPanelFX.ExperimentRow> fileRows = ioPanel.getExperimentData();
            consolePanel.appendConsoleText("Found " + fileRows.size() + " experiment rows\n");
            LOGGER.info("Found {} experiment rows", fileRows.size());
            
            if (fileRows.isEmpty()) {
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
                LOGGER.debug("Experiment {}: file={}, name={}, fraction={}", i, rawFiles[i], expNames[i], fractions[i]);
            }
            
            // Create MetaData object
            MetaData metaData = new MetaData(rawFiles, expNames, fractions, replicates);
            LOGGER.info("MetaData created successfully with {} raw files", rawFiles.length);
            
            // Collect basic parameters
            String result = ioPanel.getOutputDirectory();
            LOGGER.info("Output directory: {}", result);
            
            MagDbItem selectedMagDb = magDatabasePanel.getSelectedMagDatabase();
            String microDb = selectedMagDb != null ? selectedMagDb.getCatalogueID() : "";
            LOGGER.info("Selected MAG database: {}", microDb);

            MagVersion selectedVersion = magDatabasePanel.getSelectedMagVersion();
            
            boolean appendHostDb = magDatabasePanel.isAppendHostDb();
            String hostDB = appendHostDb ? magDatabasePanel.getHostDatabase() : "";
            
            String enzyme = enzymeComboBox.getValue() != null ? enzymeComboBox.getValue() : "Trypsin";
            int missCleavages = missedCleavagesSpinner.getValue() != null ? missedCleavagesSpinner.getValue() : 2;
            String digestMode = digestionModeComboBox.getValue() != null ? digestionModeComboBox.getValue() : "Specific";
            int threadCount = threadsComboBox.getValue() != null ? threadsComboBox.getValue() : 4;
            String customBasicPath = customBasicWorkflowField != null ? customBasicWorkflowField.getText() : null;
            String customCombinedPath = customCombinedWorkflowField != null ? customCombinedWorkflowField.getText() : null;
            // FragPipe MAG is closed-search only here. Open search lives in the dedicated
            // open-search workflow. Without a custom override, use the historical LFQ-MBR
            // template (marker-stage quant flags are unused).
            String workflow = (customBasicPath != null && !customBasicPath.isBlank()) ? null : "LFQ-MBR";
            String msDataType = msDataTypeComboBox.getValue() != null ? msDataTypeComboBox.getValue() : "Regular MS";
            
            // Create DDAMagFragPipeParameters with required constructor parameters
			DDAMagFragPipeParameters analysisParams = new DDAMagFragPipeParameters(
					uk.ac.quadram.metapilot.workflow.WorkflowType.DDA_MAG_FRAGPIPE.getAnalysisType(), metaData, result, microDb,
					appendHostDb, hostDB, msDataType, // Default MS2 scan mode for FragPipe
					threadCount);
            LOGGER.info("DDAMagFragPipeParameters object created successfully");
            
            // Set additional parameters
            analysisParams.setEnzyme(enzyme);
            analysisParams.setMissCleavages(missCleavages);
            analysisParams.setDigestMode(digestMode);
            
            // Set MAG-specific fields
            analysisParams.setMagDbCatalogueId(microDb);
            if (selectedMagDb != null && selectedMagDb.getIdentifier() != null
                    && !selectedMagDb.getIdentifier().isBlank()) {
                analysisParams.setMagDbIdentifier(selectedMagDb.getIdentifier());
            }
            if (selectedVersion != null) {
                analysisParams.setMagDbVersion(selectedVersion.getCatalogueVersion());
            } else if (selectedMagDb != null && selectedMagDb.getVersionList() != null && !selectedMagDb.getVersionList().isEmpty()) {
                // Fallback: keep previous behavior if version combobox isn't ready
                analysisParams.setMagDbVersion(selectedMagDb.getVersionList().get(0).getCatalogueVersion());
            }
            
            // Set the used version on the MagDbItem
            if (selectedVersion != null && selectedMagDb != null) {
                selectedMagDb.setUsedVersion(selectedVersion);
                LOGGER.info("Set used version on MagDbItem: {}", selectedVersion.getCatalogueVersion());
            }
            
            analysisParams.setMagDbItem(selectedMagDb);

            String genomeSelectorMode = genomeSelectorComboBox != null && genomeSelectorComboBox.getValue() != null
                    ? genomeSelectorComboBox.getValue() : MagDbItem.MARKER2;
            analysisParams.setGenomeSelectorMode(genomeSelectorMode);
            if (genomeStopControls != null) {
                analysisParams.setGenomeSelectionStrategy(genomeStopControls.strategy());
                analysisParams.setGenomeSelectionPsmCoverage(genomeStopControls.psmCoverage());
                analysisParams.setGenomeSelectionStepGuard(genomeStopControls.stepGuard());
                analysisParams.setGenomeSelectionMaxGenomes(genomeStopControls.maxGenomes());
            }

            analysisParams.setSearchMode("Closed");
            analysisParams.setOpenSearchAtMarkerStage(false);
            if (customBasicPath != null && !customBasicPath.isBlank()) {
                analysisParams.setCustomBasicWorkflowPath(customBasicPath.trim());
            }
            if (customCombinedPath != null && !customCombinedPath.isBlank()) {
                analysisParams.setCustomCombinedWorkflowPath(customCombinedPath.trim());
            }

            // Set quantitative mode based on workflow (default to label-free)
            analysisParams.setQuantitativeMode("Label-Free");
            
            analysisParams.setFixedModifications(modificationsPanel.getFixedModifications());
            analysisParams.setVariableModifications(modificationsPanel.getVariableModifications());
            
            // Set experiment data - convert JavaFX rows to String arrays
            List<String[]> experimentArrays = new ArrayList<>();
            for (MetaPilotIOPanelFX.ExperimentRow row : fileRows) {
                experimentArrays.add(row.toArray());
            }
            analysisParams.setExperimentData(experimentArrays);
            
            // Set metadata table data using reusable component - convert JavaFX rows to String arrays
            if (metadataPanel != null) {
                List<String[]> metadataArrays = new ArrayList<>();
                for (var row : metadataPanel.getMetadataTableData()) {
                    metadataArrays.add(row.toArray());
                }
                analysisParams.setMetadataTableData(metadataArrays);
                metadataPanel.applyTo(metaData);
            }
            
            // Create and configure FragPipe workflow (skipped entirely when the user
            // supplied an explicit custom basic-workflow override above - that file is
            // trusted as-is, per the Custom-mode escape hatch).
            try {
                // Get FragPipe path from MetaResourceParameters
				String fragpipePath = resourceParams.getToolPath(MetaResourceParameters.ToolType.FRAGPIPE);

                if (workflow == null) {
                    LOGGER.info("Custom basic workflow template supplied ({}); skipping GUI-side workflow build.",
                            customBasicPath);
                } else if (fragpipePath != null && !fragpipePath.trim().isEmpty()) {
                    File fragpipeFile = new File(fragpipePath);
                    if (fragpipeFile.exists()) {
                        File binDir = fragpipeFile.getParentFile();
                        if (binDir != null) {
                            File workflowsDir = new File(binDir.getParentFile(), "workflows");
                            if (workflowsDir.exists() && workflowsDir.isDirectory()) {
                                // Parse the selected workflow
                                uk.ac.quadram.metapilot.dbSearch.fragpipe.FragpipeWorkflow fragpipeWorkflow = 
                                    uk.ac.quadram.metapilot.dbSearch.fragpipe.FragpipeWorkflow.parse(workflowsDir, workflow);
                                
                                if (fragpipeWorkflow != null) {
                                    // Configure the workflow with UI parameters
                                    
                                    // The core runner writes enzyme, digestion and modifications
                                    // into every stage's workflow; set here only so this file agrees.
                                    String[] enzymeInfo = MSFraggerHeadlessRunner.enzymeInfo(enzyme);
                                    fragpipeWorkflow.setEnzyme(enzymeInfo[0], enzymeInfo[1], enzymeInfo[2],
                                            MSFraggerHeadlessRunner.enzymeTermini(digestMode, enzyme), missCleavages);

                                    // Precursor/fragment tolerance intentionally left untouched here -
                                    // it comes from the LFQ-MBR (or custom) template parsed above.

                                    // Set ion mobility based on MS data type
                                    boolean isIonMobility = msDataType.contains("IM-MS");
                                    fragpipeWorkflow.setIonMobility(isIonMobility);
                                    
                                    // Store the configured workflow as member variable (will be passed to task separately)
                                    this.currentFragpipeWorkflow = fragpipeWorkflow;
                                    
                                    LOGGER.info("Successfully configured FragPipe workflow: {}", workflow);
                                } else {
                                    LOGGER.warn("Failed to parse FragPipe workflow: {}", workflow);
                                }
                            } else {
                                LOGGER.warn("FragPipe workflows directory not found: {}", workflowsDir.getAbsolutePath());
                            }
                        }
                    } else {
                        LOGGER.warn("FragPipe executable not found: {}", fragpipePath);
                    }
                } else {
                    LOGGER.warn("FragPipe path not configured in MetaResourceParameters");
                }
            } catch (Exception e) {
                LOGGER.error("Failed to configure FragPipe workflow: " + e.getMessage(), e);
            }
            
            consolePanel.appendConsoleText("✓ Successfully collected parameters from UI\n");
            LOGGER.info("Successfully collected parameters from UI for MAG FragPipe workflow");
            return analysisParams;
            
        } catch (Exception e) {
            String errorMsg = "ERROR: Failed to collect parameters from UI: " + e.getMessage();
            consolePanel.appendConsoleText(errorMsg + "\n");
            LOGGER.error("Failed to collect parameters from UI: " + e.getMessage(), e);
            return null;
        }
    }
    
    /**
     * Validate inputs before starting analysis
     */
    private boolean validateInputs() {
        LOGGER.info("Starting input validation for MAG FragPipe workflow");
        
        // Check experiment data
        if (ioPanel == null) {
            LOGGER.error("IO panel is null");
            return false;
        }
        
        if (ioPanel.getExperimentData().isEmpty()) {
            LOGGER.error("No experiment data found in IO panel");
            return false;
        }
        LOGGER.info("Experiment data validation passed: {} experiments found", ioPanel.getExperimentData().size());
        
        // Check output directory
        String outputDir = ioPanel.getOutputDirectory();
        if (outputDir == null || outputDir.trim().isEmpty()) {
            LOGGER.error("Output directory is empty or null");
            return false;
        }
        LOGGER.info("Output directory validation passed: {}", outputDir);
        
        // Check MAG database selection using reusable component
        if (magDatabasePanel == null) {
            LOGGER.error("MAG database panel is null");
            return false;
        }
        
        if (!magDatabasePanel.validateSelection()) {
            String errorMsg = magDatabasePanel.getValidationErrorMessage();
            LOGGER.error("MAG database validation failed: {}", errorMsg);
            return false;
        }
        
        MagDbItem selectedMagDb = magDatabasePanel.getSelectedMagDatabase();
        MagVersion selectedVersion = magDatabasePanel.getSelectedMagVersion();
        LOGGER.info("MAG database validation passed: {} (version: {})", 
            selectedMagDb.getCatalogueID(), selectedVersion.getCatalogueVersion());
        
        LOGGER.info("All input validation checks passed successfully");
        return true;
    }
    
    /**
     * Load last used parameters and apply them to the UI
     */
    @Override
    protected void loadLastUsedParameters() {
        // Check if parameters have already been loaded from an external file
        // If so, skip loading last used parameters to avoid overwriting
        if (parametersLoadedFromFile) {
            System.out.println("DEBUG: Skipping loadLastUsedParameters() - parameters already loaded from file");
            consolePanel.appendConsoleText("Skipping last used parameters - using loaded file parameters\n");
            return;
        }
        
        try {
            consolePanel.appendConsoleText("Loading last used parameters...\n");
            LOGGER.info("Attempting to load last used parameters for DDA MAG FragPipe workflow");
            DDAMagFragPipeParameters params = (DDAMagFragPipeParameters) parameterIO.loadLastUsedParameters();
            if (params != null) {
                consolePanel.appendConsoleText("✓ Successfully loaded last used parameters\n");
                LOGGER.info("Successfully loaded parameters from file");
                applyParametersToUI(params);
                LOGGER.info("Applied parameters to UI for DDA MAG FragPipe workflow");
                
                // Also apply experiment data to the IO panel if available - convert String arrays back to JavaFX rows
                if (params.getExperimentData() != null && !params.getExperimentData().isEmpty()) {
                    List<MetaPilotIOPanelFX.ExperimentRow> experimentRows = new ArrayList<>();
                    for (String[] data : params.getExperimentData()) {
                        experimentRows.add(MetaPilotIOPanelFX.ExperimentRow.fromArray(data));
                    }
                    ioPanel.setExperimentData(experimentRows);
                    consolePanel.appendConsoleText("Applied " + params.getExperimentData().size() + " experiment rows to IO panel\n");
                    LOGGER.info("Applied {} experiment rows to IO panel", params.getExperimentData().size());
                }
                
                // Apply output directory if available
                if (params.getResult() != null && !params.getResult().trim().isEmpty()) {
                    ioPanel.setOutputDirectory(params.getResult());
                    consolePanel.appendConsoleText("Applied output directory: " + params.getResult() + "\n");
                    LOGGER.info("Applied output directory: {}", params.getResult());
                }
                
                // Update metadata table after loading parameters
                Platform.runLater(() -> {
                    updateMetadataTable();
                    updateParameterView();
                });
            } else {
                consolePanel.appendConsoleText("No previous parameters found, creating defaults\n");
                LOGGER.info("No parameters loaded from file, creating defaults");
                DDAMagFragPipeParameters defaultParams = new DDAMagFragPipeParameters();
                applyParametersToUI(defaultParams);
                consolePanel.appendConsoleText("✓ Default parameters loaded into UI\n");
                LOGGER.info("Applied default parameters to UI successfully");
            }
        } catch (Exception e) {
            String errorMsg = "WARNING: Failed to load last used parameters: " + e.getMessage();
            consolePanel.appendConsoleText(errorMsg + "\n");
            LOGGER.error("Failed to load last used MAG FragPipe parameters: " + e.getMessage(), e);
            // Continue with default UI state
        }
    }
    
    /**
     * Apply parameters to the UI components
     */
    @Override
    protected void applyParametersToUI(DDAMagFragPipeParameters params) {
        try {
            consolePanel.appendConsoleText("Applying parameters to UI...\n");
            if (params == null) {
                consolePanel.appendConsoleText("WARNING: No parameters to apply\n");
                LOGGER.warn("No parameters to apply to UI");
                return;
            }
        
        // Apply experiment data to IO panel if available - convert String arrays back to JavaFX rows
        if (params.getExperimentData() != null && !params.getExperimentData().isEmpty()) {
            List<ExperimentRow> experimentRows = new ArrayList<>();
            for (String[] data : params.getExperimentData()) {
                experimentRows.add(ExperimentRow.fromArray(data));
            }
            ioPanel.setExperimentData(experimentRows);
            consolePanel.appendConsoleText("Applied " + params.getExperimentData().size() + " experiment rows to IO panel\n");
            LOGGER.info("Applied {} experiment rows to IO panel", params.getExperimentData().size());
        } else if (params.getMetadata() != null) {
            // Fallback: Apply metadata to IO panel if experiment data is not available
            MetaData metadata = params.getMetadata();
            String[] rawFiles = metadata.getRawFiles();
            String[] expNames = metadata.getExpNames();
            
            if (rawFiles != null && expNames != null && rawFiles.length > 0) {
                ObservableList<MetaPilotIOPanelFX.ExperimentRow> experimentRows = FXCollections.observableArrayList();
                for (int i = 0; i < rawFiles.length; i++) {
                    if (rawFiles[i] != null && !rawFiles[i].trim().isEmpty()) {
                        experimentRows.add(new MetaPilotIOPanelFX.ExperimentRow(rawFiles[i], expNames[i]));
                    }
                }
                ioPanel.setExperimentData(experimentRows);
                LOGGER.info("Applied metadata with {} raw files to IO panel", rawFiles.length);
            }
        }
        
        // Apply output directory if available
        if (params.getResult() != null && !params.getResult().trim().isEmpty()) {
            ioPanel.setOutputDirectory(params.getResult());
            consolePanel.appendConsoleText("Applied output directory: " + params.getResult() + "\n");
            LOGGER.info("Applied output directory: {}", params.getResult());
        }
        
        // Apply only essential parameters from saved settings
        if (threadsComboBox != null) {
            threadsComboBox.setValue(params.getThreadCount());
        }
        if (params.getGenomeSelectorMode() != null && genomeSelectorComboBox != null) {
            genomeSelectorComboBox.setValue(params.getGenomeSelectorMode());
        }
        if (genomeStopControls != null) {
            genomeStopControls.apply(params.getGenomeSelectionStrategy(), params.getGenomeSelectionPsmCoverage(),
                    params.getGenomeSelectionStepGuard(), params.getGenomeSelectionMaxGenomes());
            String warning = GenomeSelectionCut.unusedCoverageKeysWarning(params.getGenomeSelectionStrategy(),
                    params.getGenomeSelectionPsmCoverage(), params.getGenomeSelectionStepGuard(),
                    params.getGenomeSelectionMaxGenomes());
            if (warning != null && consolePanel != null) {
                consolePanel.appendConsoleText("WARNING: " + warning + "\n");
            }
        }
        
        // Custom workflow path overrides are persisted; search mode is always Closed
        // in this panel (open search uses the dedicated open-search workflow).
        if (customBasicWorkflowField != null) {
            customBasicWorkflowField.setText(
                    params.getCustomBasicWorkflowPath() != null ? params.getCustomBasicWorkflowPath() : "");
        }
        if (customCombinedWorkflowField != null) {
            customCombinedWorkflowField.setText(
                    params.getCustomCombinedWorkflowPath() != null ? params.getCustomCombinedWorkflowPath() : "");
        }
        
        if (msDataTypeComboBox != null) {
            boolean ionMobility = params.getMs2ScanMode() != null
                    && params.getMs2ScanMode().toUpperCase().contains("IM-MS");
            msDataTypeComboBox.setValue(msDataTypeComboBox.getItems().get(ionMobility ? 1 : 0));
        }
        if (enzymeComboBox != null && params.getEnzyme() != null
                && MSFraggerHeadlessRunner.ENZYME_CHOICES.contains(params.getEnzyme())) {
            enzymeComboBox.setValue(params.getEnzyme());
        }
        if (missedCleavagesSpinner != null) {
            missedCleavagesSpinner.getValueFactory().setValue(params.getMissCleavages());
        }
        if (digestionModeComboBox != null) {
            int termini = 2;
            try {
                termini = MSFraggerHeadlessRunner.enzymeTermini(params.getDigestModeString(), null);
            } catch (IllegalArgumentException e) {
                LOGGER.warn("{}; showing {}", e.getMessage(), MSFraggerHeadlessRunner.DIGEST_SPECIFIC);
            }
            digestionModeComboBox.setValue(termini == 2 ? MSFraggerHeadlessRunner.DIGEST_SPECIFIC
                    : termini == 1 ? MSFraggerHeadlessRunner.DIGEST_SEMI_SPECIFIC
                    : MSFraggerHeadlessRunner.DIGEST_NONSPECIFIC);
        }
        
        // Apply database settings using reusable component
        if (magDatabasePanel != null) {
            String targetMagDbId = params.getMagDbCatalogueId();
            if (targetMagDbId == null || targetMagDbId.trim().isEmpty()) {
                // Backward compatibility: older saved parameters store MAG DB in MetaBaseParameter.microDb
                targetMagDbId = params.getMicroDb();
            }
            if (targetMagDbId != null && !targetMagDbId.trim().isEmpty()) {
                // Find and select the matching MAG database
                for (MagDbItem item : magDatabasePanel.getAvailableMagDatabases()) {
                    if (item.getCatalogueID().equals(targetMagDbId)) {
                        magDatabasePanel.setSelectedMagDatabase(item);
                        
                        // Also apply version selection if available
                        if (params.getMagDbVersion() != null) {
                            for (MagVersion version : magDatabasePanel.getSelectedMagDatabase().getVersionList()) {
                                if (version.getCatalogueVersion().equals(params.getMagDbVersion())) {
                                    magDatabasePanel.setSelectedMagVersion(version);
                                    LOGGER.info("Applied version selection: {}", version.getCatalogueVersion());
                                    break;
                                }
                            }
                        }
                        break;
                    }
                }
            }
            magDatabasePanel.setAppendHostDb(params.isAppendHostDb());
            if (params.getHostDB() != null) {
                magDatabasePanel.setHostDatabase(params.getHostDB());
            }
        }
        
        // A file saved before modifications reached FragPipe keeps the defaults shown,
        // which are the template modifications those runs searched.
        if (modificationsPanel != null) {
            modificationsPanel.setModifications(params.getFixedModifications(), params.getVariableModifications());
            LOGGER.info("Applied {} fixed and {} variable modifications", params.getFixedModifications().size(),
                    params.getVariableModifications().size());
        }
        
        // Note: Quantitative mode is determined by workflow selection, no UI settings needed

        // Apply metadata table data if available - convert String arrays back to JavaFX rows
        if (params.getMetadataTableData() != null && !params.getMetadataTableData().isEmpty() && metadataPanel != null) {
            List<uk.ac.quadram.metapilot.task.dda.igc.gui.MetaPilotDdaIgcMaxQuantPanel.MetadataRow> metadataRows = new ArrayList<>();
            for (String[] data : params.getMetadataTableData()) {
                metadataRows.add(uk.ac.quadram.metapilot.task.dda.igc.gui.MetaPilotDdaIgcMaxQuantPanel.MetadataRow.fromArray(data));
            }
            metadataPanel.setMetadataTableDataWithColumns(metadataRows);
            LOGGER.info("Applied {} metadata table rows with columns", params.getMetadataTableData().size());
        }
        
        // Update parameter view after applying parameters
        Platform.runLater(() -> updateParameterView());
        
        consolePanel.appendConsoleText("✓ Successfully applied parameters to UI\n");
        LOGGER.info("Successfully applied parameters to UI for MAG FragPipe workflow");
        
        } catch (Exception e) {
            String errorMsg = "ERROR: Failed to apply parameters to UI: " + e.getMessage();
            consolePanel.appendConsoleText(errorMsg + "\n");
            LOGGER.error("Failed to apply parameters to UI for MAG FragPipe workflow", e);
        }
    }
    
    /**
     * Load analysis parameters from a specific result directory
     */
    public void loadAnalysisParameters(String resultDirectory) {
        try {
            LOGGER.info("Loading MAG FragPipe analysis parameters from: {}", resultDirectory);
            DDAMagFragPipeParameters params = (DDAMagFragPipeParameters) parameterIO.loadAnalysisParameters(resultDirectory);
            if (params != null) {
                applyParametersToUI(params);
                LOGGER.info("Loaded MAG FragPipe analysis parameters from {}", resultDirectory);
                consolePanel.appendConsoleText("✓ Loaded MAG FragPipe analysis parameters from " + resultDirectory + "\n");
            } else {
                LOGGER.warn("No parameters found in {}", resultDirectory);
                consolePanel.appendConsoleText("WARNING: No parameters found in " + resultDirectory + "\n");
            }
        } catch (Exception e) {
            LOGGER.error("Failed to load MAG FragPipe analysis parameters from {}", resultDirectory, e);
            consolePanel.appendConsoleText("WARNING: Failed to load parameters from " + resultDirectory + ": " + e.getMessage() + "\n");
        }
    }
    
    /**
     * Get the parameter class for this workflow
     */
    @Override
    protected Class<DDAMagFragPipeParameters> getParameterClass() {
        return DDAMagFragPipeParameters.class;
    }
    
    /**
     * Show error dialog
     */
    
}
