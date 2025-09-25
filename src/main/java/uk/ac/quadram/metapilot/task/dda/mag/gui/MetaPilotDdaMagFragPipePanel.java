package uk.ac.quadram.metapilot.task.dda.mag.gui;

import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.control.*;
import javafx.scene.layout.*;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.scene.text.Font;
import javafx.scene.text.FontWeight;
import javafx.scene.paint.Color;
import uk.ac.quadram.metapilot.task.common.gui.MetaPilotIOPanelFX;
import uk.ac.quadram.metapilot.task.common.gui.MagDatabaseSettingsPanel;
import uk.ac.quadram.metapilot.task.common.gui.ConsolePanel;
import uk.ac.quadram.metapilot.task.common.gui.MetadataConfigurationPanel;
import uk.ac.quadram.metapilot.task.common.gui.ModificationsPanel;
import uk.ac.quadram.metapilot.task.dda.mag.par.DDAMagFragPipeParameters;
import uk.ac.quadram.metapilot.task.dda.mag.task.DDAMagFragPipeTask;
import uk.ac.quadram.metapilot.task.dda.mag.par.DDAMagFragPipeParameterIO;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;
import uk.ac.quadram.metapilot.task.mag.MagDbItem.MagVersion;
import uk.ac.quadram.metapilot.task.par.MetaData;
import uk.ac.quadram.metapilot.MetaPilotApp;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.task.common.gui.AbstractWorkflowPanel;
import javafx.application.Platform;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.IOException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.List;
import java.util.ArrayList;
import java.io.File;

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
    private Spinner<Integer> threadsSpinner;
    private ComboBox<String> workflowComboBox;
    private ComboBox<String> msDataTypeComboBox;

    // Data
    private ObservableList<String> enzymes = FXCollections.observableArrayList();
    private ObservableList<String> workflows = FXCollections.observableArrayList();

    // Add parameter IO instance
    private DDAMagFragPipeParameterIO parameterIO;
    private MetaResourceParameters resourceParams;
    
    // Store the last used analysis parameters for viewReport
    private DDAMagFragPipeParameters currentAnalysisParams;
    
    // Current task reference
    private DDAMagFragPipeTask currentTask;
    
    // Logger
    private static final Logger LOGGER = LoggerFactory.getLogger(MetaPilotDdaMagFragPipePanel.class);

	public MetaPilotDdaMagFragPipePanel(MagDbItem[] magDbItems, MetaResourceParameters resourceParams) {
		// Call parent constructor first
		super();
		
		// Initialize parameter IO instance
		parameterIO = new DDAMagFragPipeParameterIO();
		this.resourceParams = resourceParams;
		
		loadEnzymesFromIni();
		loadWorkflowsFromFragPipe();
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
        metadataPanel = new MetadataConfigurationPanel();
        magDatabasePanel = new MagDatabaseSettingsPanel();
        
        // Set up callbacks for parameter updates
        metadataPanel.setParameterUpdateCallback(this::updateParameterView);
        magDatabasePanel.setParameterUpdateCallback(this::updateParameterView);
        consolePanel.setStartCallback(this::startAnalysis);
        consolePanel.setStopCallback(this::stopAnalysis);
        consolePanel.setViewReportCallback(this::viewReport);

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
        tabContent.setMinHeight(760);

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

        // Left column: Enzyme Settings, Additional Parameters, and Modifications Panel (NO taxonomic analysis for MAG)
        VBox leftColumn = new VBox(15);
        leftColumn.setMinWidth(430);
        VBox enzymeSection = createEnzymeSection();
        VBox additionalParamsSection = createAdditionalParametersSection();
        
        // Create horizontal modifications panel with fixed mods (narrow) and variable mods (wider)
        modificationsPanel = new ModificationsPanel(this::updateParameterView);
        
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
        
        tabContent.getChildren().add(mainGrid);
        return tabContent;
    }

    // Add missing utility methods
    private void loadEnzymesFromIni() {
        try (InputStream in = getClass().getClassLoader().getResourceAsStream("enzyme_default.ini")) {
            if (in != null) {
                BufferedReader reader = new BufferedReader(new InputStreamReader(in));
                Pattern pattern = Pattern.compile("EN\\d+=(.+?) ");
                String line;
                while ((line = reader.readLine()) != null) {
                    Matcher m = pattern.matcher(line);
                    if (m.find()) {
                        String enzymeName = m.group(1).trim();
                        enzymes.add(enzymeName);
                    }
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

	private void loadWorkflowsFromFragPipe() {
		try {
			String fragpipePath = resourceParams.getToolPath(MetaResourceParameters.ToolType.FRAGPIPE);

			if (fragpipePath != null && !fragpipePath.trim().isEmpty()) {
				File fragpipeFile = new File(fragpipePath);
				if (fragpipeFile.exists()) {
					// Get parent directory (bin)
					File binDir = fragpipeFile.getParentFile();
					if (binDir != null) {
						// Get workflows directory (same level as bin)
						File workflowsDir = new File(binDir.getParentFile(), "workflows");
						if (workflowsDir.exists() && workflowsDir.isDirectory()) {
							// List all .workflow files
							File[] workflowFiles = workflowsDir
									.listFiles((dir, name) -> name.toLowerCase().endsWith(".workflow"));
							if (workflowFiles != null) {
								for (File workflowFile : workflowFiles) {
									String fileName = workflowFile.getName();
									// Remove .workflow extension
									String workflowName = fileName.substring(0, fileName.lastIndexOf('.'));
									workflows.add(workflowName);
								}
								LOGGER.info("Loaded {} workflows from FragPipe: {}", workflows.size(),
										workflowsDir.getAbsolutePath());
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
			LOGGER.error("Failed to load workflows from FragPipe: " + e.getMessage(), e);
		}

		// Add default workflow if none found
		if (workflows.isEmpty()) {
			workflows.add("Standard");
			LOGGER.warn("No workflows found, using default: Standard");
		}
	}

    private Label createSectionTitle(String title) {
        Label label = new Label(title);
        label.setFont(Font.font("System", FontWeight.BOLD, 14));
        label.setTextFill(Color.DARKBLUE);
        return label;
    }

    private VBox createEnzymeSection() {
        VBox section = new VBox(10);
        section.setPadding(new Insets(10));
        section.setStyle("-fx-background-color: #f5f5f5; -fx-background-radius: 5;");

        Label titleLabel = createSectionTitle("Enzyme Settings");

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

        digestionModeComboBox = new ComboBox<>(FXCollections.observableArrayList(
            "Specific", "Semi free N-term", "Semi free C-term", "Semi specific", "Unspecific", "None"));
        digestionModeComboBox.setMinWidth(140);
        digestionModeComboBox.setPrefWidth(160);
        digestionModeComboBox.setMaxWidth(Double.MAX_VALUE);
        digestionModeComboBox.setValue("Specific");
        digestionModeComboBox.setOnAction(e -> updateParameterView());
        digestionModeBox.getChildren().add(digestionModeComboBox);

        enzymeParamsBox.getChildren().addAll(enzymeBox, missedCleavagesBox, digestionModeBox);

        section.getChildren().addAll(titleLabel, enzymeParamsBox);
        return section;
    }

    private VBox createAdditionalParametersSection() {
        VBox section = new VBox(10);
        section.setPadding(new Insets(10));
        section.setStyle("-fx-background-color: #f0f8ff; -fx-background-radius: 5;");
        
        Label titleLabel = createSectionTitle("Additional Parameters");
        
        GridPane paramGrid = new GridPane();
        paramGrid.setHgap(15);
        paramGrid.setVgap(8);
        
        ColumnConstraints col1 = new ColumnConstraints();
        col1.setMinWidth(150);
        col1.setHgrow(Priority.ALWAYS);
        ColumnConstraints col2 = new ColumnConstraints();
        col2.setMinWidth(150);
        col2.setHgrow(Priority.ALWAYS);
        paramGrid.getColumnConstraints().addAll(col1, col2);
        
        Label threadsLabel = new Label("Threads:");
        threadsLabel.setStyle("-fx-font-weight: bold;");
        threadsSpinner = new Spinner<>(1, Runtime.getRuntime().availableProcessors(), 4);
        threadsSpinner.setMinWidth(100);
        threadsSpinner.setPrefWidth(100);
        threadsSpinner.setMaxWidth(Double.MAX_VALUE);
        threadsSpinner.valueProperty().addListener((obs, oldVal, newVal) -> updateParameterView());
        
        VBox threadsBox = new VBox(5);
        threadsBox.getChildren().addAll(threadsLabel, threadsSpinner);
        
        Label workflowLabel = new Label("Workflow:");
        workflowLabel.setStyle("-fx-font-weight: bold;");
        workflowComboBox = new ComboBox<>(workflows);
        workflowComboBox.setMinWidth(100);
        workflowComboBox.setPrefWidth(120);
        workflowComboBox.setMaxWidth(Double.MAX_VALUE);
		if (workflows.size() > 0) {
			if (workflows.contains("LFQ-MBR")) {
				workflowComboBox.setValue("LFQ-MBR");
				LOGGER.info("Set default workflow to: LFQ-MBR");
			} else if (!workflows.isEmpty()) {
				workflowComboBox.setValue(workflows.get(0));
				LOGGER.info("Set default workflow to: {}", workflows.get(0));
			}
		}
        workflowComboBox.setOnAction(e -> updateParameterView());
        
        VBox workflowBox = new VBox(5);
        workflowBox.getChildren().addAll(workflowLabel, workflowComboBox);
        
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
        
        paramGrid.add(threadsBox, 0, 0);
        paramGrid.add(workflowBox, 1, 0);
        paramGrid.add(msDataTypeBox, 0, 1, 2, 1); // Span both columns

        section.getChildren().addAll(titleLabel, paramGrid);
        return section;
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
        if (threadsSpinner != null) {
            sb.append("  Threads: ").append(threadsSpinner.getValue()).append("\n");
        }
        if (workflowComboBox != null) {
            sb.append("  Workflow: ").append(workflowComboBox.getValue()).append("\n");
        }
        if (msDataTypeComboBox != null) {
            sb.append("  MS data type: ").append(msDataTypeComboBox.getValue()).append("\n");
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
            ProgressBar overallBar = new ProgressBar();
            ProgressBar currentBar = new ProgressBar();
            Label statusLabel = new Label("Ready to start analysis...");
            consolePanel.appendConsoleText("✓ Progress bars and status label initialized\n");
            
            // 5. Set status label in parameters
            analysisParams.setStatusLabel(statusLabel);
            consolePanel.appendConsoleText("✓ Status label set in parameters\n");
            
            // 6. Create and run DDAMagFragPipeTask with DDAMagFragPipeParameters directly
            consolePanel.appendConsoleText("Creating MAG FragPipe task...\n");
			currentTask = new DDAMagFragPipeTask(analysisParams, resourceParams, consolePanel.getOverallProgressBar(),
					consolePanel.getCurrentTaskProgressBar());
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
                    consolePanel.appendConsoleText("✓ Analysis completed successfully!\n");
                    consolePanel.enableViewReportButton();
                    consolePanel.setButtonsEnabled(true, false, true);
                });
            });
            
            currentTask.setOnFailed(e -> {
                Platform.runLater(() -> {
                    consolePanel.appendConsoleText("✗ Analysis failed. Check console output for details.\n");
                    consolePanel.setButtonsEnabled(true, false, false);
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
        consolePanel.appendConsoleText("Stopping FragPipe analysis...\n");
        System.out.println("DEBUG: stopAnalysis called");
        
        if (currentTask != null && currentTask.isRunning()) {
            consolePanel.appendConsoleText("Force stopping MAG FragPipe task...\n");
            currentTask.forceStop();
            consolePanel.appendConsoleText("MAG FragPipe task stop requested.\n");
        } else {
            consolePanel.appendConsoleText("No running MAG FragPipe task to stop.\n");
        }
    }
    
    
    /**
     * Collect all parameters from the UI into a parameters object
     */
    @Override
    protected DDAMagFragPipeParameters collectParametersFromUI() {
        try {
            LOGGER.info("Starting parameter collection from UI for MAG FragPipe workflow");
            
            // Validate inputs first
            if (!validateInputs()) {
                LOGGER.error("Input validation failed");
                return null;
            }
            LOGGER.info("Input validation passed");
            
            // Collect experiment data from ioPanel
            ObservableList<MetaPilotIOPanelFX.ExperimentRow> fileRows = ioPanel.getExperimentData();
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
            
            boolean appendHostDb = magDatabasePanel.isAppendHostDb();
            String hostDB = appendHostDb ? magDatabasePanel.getHostDatabase() : "";
            
            String enzyme = enzymeComboBox.getValue() != null ? enzymeComboBox.getValue() : "Trypsin";
            int missCleavages = missedCleavagesSpinner.getValue() != null ? missedCleavagesSpinner.getValue() : 2;
            String digestMode = digestionModeComboBox.getValue() != null ? digestionModeComboBox.getValue() : "Specific";
            int threadCount = threadsSpinner.getValue() != null ? threadsSpinner.getValue() : 4;
            String workflow = workflowComboBox.getValue() != null ? workflowComboBox.getValue() : "Standard";
            String msDataType = msDataTypeComboBox.getValue() != null ? msDataTypeComboBox.getValue() : "Regular MS";
            
            // Create DDAMagFragPipeParameters with required constructor parameters
            DDAMagFragPipeParameters analysisParams = new DDAMagFragPipeParameters(
                MetaPilotApp.WorkflowType.DDA_MAG_FRAGPIPE.getAnalysisType(),
                metaData,
                result,
                microDb,
                appendHostDb,
                hostDB,
                "HCD-FTMS", // Default MS2 scan mode for FragPipe
                threadCount
            );
            LOGGER.info("DDAMagFragPipeParameters object created successfully");
            
            // Set additional parameters
            analysisParams.setEnzyme(enzyme);
            analysisParams.setMissCleavages(missCleavages);
            analysisParams.setDigestMode(digestMode);
            
            // Set MAG-specific fields
            analysisParams.setMagDbCatalogueId(microDb);
            if (selectedMagDb != null && selectedMagDb.getVersionList() != null && !selectedMagDb.getVersionList().isEmpty()) {
                String magDbVersion = selectedMagDb.getVersionList().get(0).getCatalogueVersion();
                analysisParams.setMagDbVersion(magDbVersion);
            }
            
            // Set the used version on the MagDbItem
            MagVersion selectedVersion = magDatabasePanel.getSelectedMagVersion();
            if (selectedVersion != null && selectedMagDb != null) {
                selectedMagDb.setUsedVersion(selectedVersion);
                LOGGER.info("Set used version on MagDbItem: {}", selectedVersion.getCatalogueVersion());
            }
            
            analysisParams.setMagDbItem(selectedMagDb);
            
            // Set quantitative mode based on workflow (default to label-free)
            analysisParams.setQuantitativeMode("Label-Free");
            
            // Set modifications (convert from FragPipe format to standard format)
            // Note: Modifications are set later in the workflow configuration section
            
            // Set experiment data
            analysisParams.setExperimentData(new ArrayList<>(fileRows));
            
            // Set metadata table data using reusable component
            if (metadataPanel != null) {
                analysisParams.setMetadataTableData(new ArrayList<>(metadataPanel.getMetadataTableData()));
            }
            
            // Create and configure FragPipe workflow
            try {
                // Get FragPipe path from MetaResourceParameters
				String fragpipePath = resourceParams.getToolPath(MetaResourceParameters.ToolType.FRAGPIPE);
                
                if (fragpipePath != null && !fragpipePath.trim().isEmpty()) {
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
                                    
                                    // Set enzyme configuration
                                    String[] enzymeInfo = getEnzymeInfo(enzyme);
                                    fragpipeWorkflow.setEnzyme(enzymeInfo[0], enzymeInfo[1], enzymeInfo[2], 2, missCleavages);
                                    
                                    // Set mass tolerances
                                    fragpipeWorkflow.setPrecursorMassTol(20.0, 1, 20.0, 1); // 20 ppm
                                    fragpipeWorkflow.setFragmentMassTol(0.5, 0); // 0.5 Da
                                    
                                    // Set ion mobility based on MS data type
                                    boolean isIonMobility = msDataType.contains("IM-MS");
                                    fragpipeWorkflow.setIonMobility(isIonMobility);
                                    
                                    // Set modifications from UI - convert from ModificationsPanel format to FragpipeWorkflow format
                                    List<String> fixedMods = new ArrayList<>();
                                    List<String> variableMods = new ArrayList<>();
                                    
                                    // Convert fixed modifications
                                    for (ModificationsPanel.FixedModificationRow row : modificationsPanel.getFixedModData()) {
                                        if (row.getDeltaMass() != 0.0) {
                                            fixedMods.add(row.getSite() + ":" + row.getDeltaMass());
                                        }
                                    }
                                    
                                    // Convert variable modifications
                                    for (ModificationsPanel.VariableModificationRow row : modificationsPanel.getVariableModData()) {
                                        if (row.getUse()) {
                                            variableMods.add(row.getModName() + ":" + row.getSites() + ":" + row.getDeltaMass());
                                        }
                                    }
                                    
                                    // Set modifications using the string format
                                    analysisParams.setFixedModifications(fixedMods);
                                    analysisParams.setVariableModifications(variableMods);
                                    
                                    // Set the configured workflow
                                    analysisParams.setFragpipeWorkflow(fragpipeWorkflow);
                                    
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
            
            LOGGER.info("Successfully collected parameters from UI for MAG FragPipe workflow");
            return analysisParams;
            
        } catch (Exception e) {
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
        try {
            LOGGER.info("Attempting to load last used parameters for DDA MAG FragPipe workflow");
            DDAMagFragPipeParameters params = (DDAMagFragPipeParameters) parameterIO.loadLastUsedParameters();
            if (params != null) {
                LOGGER.info("Successfully loaded parameters from file");
                applyParametersToUI(params);
                LOGGER.info("Applied parameters to UI for DDA MAG FragPipe workflow");
                
                // Also apply experiment data to the IO panel if available
                if (params.getExperimentData() != null && !params.getExperimentData().isEmpty()) {
                    ioPanel.setExperimentData(params.getExperimentData());
                    LOGGER.info("Applied {} experiment rows to IO panel", params.getExperimentData().size());
                }
                
                // Apply output directory if available
                if (params.getResult() != null && !params.getResult().trim().isEmpty()) {
                    ioPanel.setOutputDirectory(params.getResult());
                    LOGGER.info("Applied output directory: {}", params.getResult());
                }
                
                // Update metadata table after loading parameters
                Platform.runLater(() -> {
                    updateMetadataTable();
                    updateParameterView();
                });
            } else {
                LOGGER.warn("No parameters loaded from file, using defaults");
            }
        } catch (Exception e) {
            LOGGER.error("Failed to load last used MAG FragPipe parameters: " + e.getMessage(), e);
            // Continue with default UI state
        }
    }
    
    /**
     * Apply parameters to the UI components
     */
    @Override
    protected void applyParametersToUI(DDAMagFragPipeParameters params) {
        if (params == null) return;
        
        // Apply metadata to IO panel if available
        if (params.getMetadata() != null) {
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
        
        // Apply only essential parameters from saved settings
        // Thread count is the only parameter that should be loaded from saved parameters
        if (threadsSpinner != null) {
            threadsSpinner.getValueFactory().setValue(params.getThreadCount());
        }
        
        // Workflow selection - don't update since workflow is not saved in parameters
        // The workflow is only used to configure FragpipeWorkflow but not persisted
        // So we keep the current UI selection (which defaults to LFQ-MBR if available)
        
        // MS data type - use default
        if (msDataTypeComboBox != null) {
            msDataTypeComboBox.setValue("Regular MS");
        }
        
        // Note: Enzyme, modifications, and other UI parameters are not loaded from saved parameters
        // They remain as user-configurable settings with their default values
        
        // Apply database settings using reusable component
        if (magDatabasePanel != null) {
            if (params.getMagDbCatalogueId() != null) {
                // Find and select the matching MAG database
                for (MagDbItem item : magDatabasePanel.getAvailableMagDatabases()) {
                    if (item.getCatalogueID().equals(params.getMagDbCatalogueId())) {
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
        
        // Apply modifications
        if (modificationsPanel != null) {
            if (params.getFixedModifications() != null) {
                // Apply saved fixed modifications
                for (String mod : params.getFixedModifications()) {
                    if (mod.contains(":")) {
                        String[] parts = mod.split(":", 2);
                        if (parts.length == 2) {
                            String site = parts[0];
                            try {
                                double mass = Double.parseDouble(parts[1]);
                                // Find and update the corresponding row
                                for (ModificationsPanel.FixedModificationRow row : modificationsPanel.getFixedModData()) {
                                    if (row.getSite().equals(site)) {
                                        row.setDeltaMass(mass);
                                        break;
                                    }
                                }
                            } catch (NumberFormatException e) {
                                LOGGER.warn("Invalid mass format in fixed modification: {}", mod);
                            }
                        }
                    }
                }
                LOGGER.info("Applied {} fixed modifications", params.getFixedModifications().size());
            }
            
            if (params.getVariableModifications() != null) {
                // Apply saved variable modifications
                for (String mod : params.getVariableModifications()) {
                    if (mod.contains(":")) {
                        String[] parts = mod.split(":", 3);
                        if (parts.length == 3) {
                            String modName = parts[0];
                            String sites = parts[1];
                            try {
                                double mass = Double.parseDouble(parts[2]);
                                // Find and update the corresponding row
                                for (ModificationsPanel.VariableModificationRow row : modificationsPanel.getVariableModData()) {
                                    if (row.getModName().equals(modName)) {
                                        row.setUse(true);
                                        row.setSites(sites);
                                        row.setDeltaMass(mass);
                                        break;
                                    }
                                }
                            } catch (NumberFormatException e) {
                                LOGGER.warn("Invalid mass format in variable modification: {}", mod);
                            }
                        }
                    }
                }
                LOGGER.info("Applied {} variable modifications", params.getVariableModifications().size());
            }
        }
        
        // Note: Quantitative mode is determined by workflow selection, no UI settings needed

        // Apply metadata table data if available
        if (params.getMetadataTableData() != null && !params.getMetadataTableData().isEmpty() && metadataPanel != null) {
            metadataPanel.setMetadataTableData(params.getMetadataTableData());
            LOGGER.info("Applied {} metadata table rows", params.getMetadataTableData().size());
        }
        
        // Update parameter view after applying parameters
        Platform.runLater(() -> updateParameterView());
    }
    
    /**
     * Get enzyme information for FragPipe workflow configuration
     * Returns [enzymeName, cutAfter, notCut] array
     */
    private String[] getEnzymeInfo(String enzyme) {
        // Map common enzyme names to FragPipe format
        switch (enzyme.toLowerCase()) {
            case "trypsin":
                return new String[]{"stricttrypsin", "KR", ""};
            case "trypsin/p":
                return new String[]{"trypsin", "KR", "P"};
            case "chymotrypsin":
                return new String[]{"chymotrypsin", "FWY", ""};
            case "gluc":
                return new String[]{"gluc", "DE", ""};
            case "gluc-bicarb":
                return new String[]{"gluc-bicarb", "DE", ""};
            case "lysc":
                return new String[]{"lysc", "K", ""};
            case "lysc-p":
                return new String[]{"lysc-p", "K", "P"};
            case "argc":
                return new String[]{"argc", "R", ""};
            case "aspn":
                return new String[]{"aspn", "D", ""};
            case "elastase":
                return new String[]{"elastase", "AGILV", ""};
            case "proteinasek":
                return new String[]{"proteinasek", "AEFILMVWY", ""};
            case "thermolysin":
                return new String[]{"thermolysin", "AEFILMV", ""};
            default:
                // Default to trypsin
                return new String[]{"stricttrypsin", "KR", ""};
        }
    }
    
    /**
     * Show error dialog
     */
    
} 