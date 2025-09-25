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
import uk.ac.quadram.metapilot.task.common.gui.ConsolePanel;
import uk.ac.quadram.metapilot.task.common.gui.MetadataConfigurationPanel;
import uk.ac.quadram.metapilot.task.common.gui.IsobaricLabelingPanel;
import uk.ac.quadram.metapilot.task.common.gui.MagDatabaseSettingsPanel;
import uk.ac.quadram.metapilot.task.common.gui.PFindModificationPanel;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;
import uk.ac.quadram.metapilot.task.dda.igc.gui.MetaPilotDdaIgcMaxQuantPanel;
import uk.ac.quadram.metapilot.task.par.MetaData;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.resources.MetaResourceIO;
import uk.ac.quadram.metapilot.task.common.gui.MetaPilotIOPanelFX.ExperimentRow;
import uk.ac.quadram.metapilot.MetaPilotApp;
import javafx.application.Platform;
import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.io.File;
import java.io.FileWriter;
import org.json.JSONArray;
import org.json.JSONObject;
import java.text.SimpleDateFormat;
import java.util.Date;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import uk.ac.quadram.metapilot.task.dda.mag.par.DDAMagPfindParameters;
import uk.ac.quadram.metapilot.task.dda.mag.task.DDAMagPfindTask;
import uk.ac.quadram.metapilot.task.dda.mag.par.DDAMagPfindParameterIO;
import uk.ac.quadram.metapilot.task.mag.MagDbItem.MagVersion;

public class MetaPilotDdaMagPfindPanel extends VBox {
    // UI Components - Database specific
    private MetaPilotIOPanelFX ioPanel;
    private MagDatabaseSettingsPanel magDatabasePanel;
    private ToggleGroup quantModeGroup;

    // Reusable UI Components
    private ConsolePanel consolePanel;
    private MetadataConfigurationPanel metadataPanel;
    private IsobaricLabelingPanel isobaricPanel;
    // pFind uses a different modification system - no modifications panel needed

    // Parameter Panel Components
    private ComboBox<String> enzymeComboBox;
    private Spinner<Integer> missedCleavagesSpinner;
    private ComboBox<String> digestionModeComboBox;

    // Additional parameters
    private Spinner<Integer> threadsSpinner;
    private ComboBox<String> ms2ScanModeComboBox;
    private CheckBox generateSampleSpecificDBCheckBox;
    private CheckBox spectraClusteringCheckBox;

    // Reusable pFind modifications panel
    private PFindModificationPanel pfindModificationsPanel;

    // Data
    private ObservableList<String> enzymes = FXCollections.observableArrayList();
    
    // Add parameter IO instance
    private DDAMagPfindParameterIO parameterIO;
    
    // Store the last used analysis parameters for viewReport
    private DDAMagPfindParameters currentAnalysisParams;
    
    // Current task reference
    private DDAMagPfindTask currentTask;
    
    public MetaPilotDdaMagPfindPanel(MagDbItem[] magDbItems) {
        // Initialize parameter IO instance
        parameterIO = new DDAMagPfindParameterIO();
        
        System.out.println("DEBUG: MetaPilotDdaMagPfindPanel() constructor started");
        try {
            loadEnzymesFromIni();
            // Initialize MAG databases from parameter
            System.out.println("DEBUG: Loading MAG databases...");
            System.out.println("DEBUG: Found " + (magDbItems != null ? magDbItems.length : 0) + " MAG database items");
            initializeUI();
            // Now set MAG databases after UI is initialized
            setAvailableMagDatabases(magDbItems);
            System.out.println("DEBUG: MetaPilotDdaMagPfindPanel(MetaParameterMag) constructor completed successfully");
        } catch (Exception e) {
            System.err.println("ERROR: Exception in MetaPilotDdaMagPfindPanel constructor: " + e.getMessage());
            e.printStackTrace();
            throw e;
        }
    }

    private void initializeUI() {
        System.out.println("DEBUG: initializeUI() started");
        try {
            System.out.println("DEBUG: Setting basic panel properties...");
            setSpacing(10);
            setPadding(new Insets(10));
            setMinWidth(1200);
            setMinHeight(800);
            setPrefWidth(1600);
            setPrefHeight(900);

            System.out.println("DEBUG: Creating reusable UI components...");
            // Initialize reusable components
            consolePanel = new ConsolePanel();
            metadataPanel = new MetadataConfigurationPanel();
            isobaricPanel = new IsobaricLabelingPanel();
            // pFind doesn't use the ModificationsPanel - it has its own modification system
            magDatabasePanel = new MagDatabaseSettingsPanel();
            
            // Set up callbacks for parameter updates
            metadataPanel.setParameterUpdateCallback(this::updateParameterView);
            isobaricPanel.setParameterUpdateCallback(() -> {
                updateParameterView();
                updateMetadataTable(); // Also update metadata when isobaric settings change
            });
            // pFind doesn't use modifications panel callback
            magDatabasePanel.setParameterUpdateCallback(this::updateParameterView);
            consolePanel.setStartCallback(this::startAnalysis);
            consolePanel.setStopCallback(this::stopAnalysis);
            consolePanel.setViewReportCallback(this::viewReport);

            System.out.println("DEBUG: Creating tab pane...");
            // Create tab pane with responsive sizing
            TabPane tabPane = new TabPane();
            tabPane.setTabClosingPolicy(TabPane.TabClosingPolicy.UNAVAILABLE);
            tabPane.setMinWidth(980);
            tabPane.setMinHeight(780);
            VBox.setVgrow(tabPane, Priority.ALWAYS);

            System.out.println("DEBUG: Creating IO tab...");
            // Create tabs
            Tab ioTab = new Tab("Input/Output");
            ioTab.setContent(createIOTab());

            System.out.println("DEBUG: Creating Parameters tab...");
            Tab paramTab = new Tab("Parameters");
            paramTab.setContent(createParametersTab());

            System.out.println("DEBUG: Creating Console tab...");
            Tab consoleTab = new Tab("Console");
            consoleTab.setContent(consolePanel);

            System.out.println("DEBUG: Adding tabs to tab pane...");
            tabPane.getTabs().addAll(ioTab, paramTab, consoleTab);

            System.out.println("DEBUG: Adding tab change listener...");
            // Add tab change listener to update metadata when Parameters tab is selected
            // and save parameters when Console tab is selected
            tabPane.getSelectionModel().selectedItemProperty().addListener((obs, oldTab, newTab) -> {
                if (newTab == paramTab) {
                    updateMetadataTableFromIOPanel();
                } else if (newTab == consoleTab) {
                    // Save current parameters when Console tab is selected
                    try {
                        consolePanel.appendConsoleText("Auto-saving parameters when switching to Console tab...\n");
                        DDAMagPfindParameters currentParams = collectParametersFromUI();
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

            System.out.println("DEBUG: Adding tab pane to main panel...");
            // Add tab pane to main panel
            getChildren().add(tabPane);
            
            System.out.println("DEBUG: initializeUI() completed successfully");
            
            // Load parameters after UI is fully initialized
            Platform.runLater(() -> {
                loadLastUsedParameters();
                
                // Add window close event handler to save parameters
                setupWindowCloseHandler();
            });
        } catch (Exception e) {
            System.err.println("ERROR: Exception in initializeUI(): " + e.getMessage());
            e.printStackTrace();
            throw e;
        }
    }

    private VBox createIOTab() {
        System.out.println("DEBUG: createIOTab() started");
        try {
            VBox tabContent = new VBox(10);
            tabContent.setPadding(new Insets(10));
            System.out.println("DEBUG: Creating MetaLabIOPanelFX...");
            ioPanel = new MetaPilotIOPanelFX();
            System.out.println("DEBUG: Adding components to IO tab...");
            tabContent.getChildren().addAll(
                createSectionTitle("Input/Output Settings"),
                ioPanel,
                createSectionTitle("Database Settings"),
                magDatabasePanel
            );
            System.out.println("DEBUG: createIOTab() completed successfully");
            return tabContent;
        } catch (Exception e) {
            System.err.println("ERROR: Exception in createIOTab(): " + e.getMessage());
            e.printStackTrace();
            throw e;
        }
    }


    // Load available MAG databases
    public void setAvailableMagDatabases(MagDbItem[] magDbItems) {
        System.out.println("DEBUG: setAvailableMagDatabases called with " + (magDbItems != null ? magDbItems.length : 0) + " items");
        if (magDatabasePanel != null) {
            magDatabasePanel.setAvailableMagDatabases(magDbItems);
            System.out.println("DEBUG: Set MAG databases in reusable panel");
        } else {
            System.out.println("DEBUG: magDatabasePanel is null - UI not initialized yet");
        }
    }
    

    private VBox createParametersTab() {
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
        
        // Column constraints for 2-column layout (50% + 50%) with minimum widths
        ColumnConstraints col1 = new ColumnConstraints();
        col1.setPercentWidth(50);
        col1.setMinWidth(450);
        col1.setHgrow(Priority.ALWAYS);
        ColumnConstraints col2 = new ColumnConstraints();
        col2.setPercentWidth(50);
        col2.setMinWidth(450);
        col2.setHgrow(Priority.ALWAYS);
        mainGrid.getColumnConstraints().addAll(col1, col2);
        
        // Row constraints for proper height distribution
        RowConstraints row1 = new RowConstraints();
        row1.setVgrow(Priority.ALWAYS);
        row1.setMinHeight(400);
        mainGrid.getRowConstraints().add(row1);
        
        // Left column: Enzyme Settings, Additional Parameters, Modifications (stacked) - NO taxonomic analysis for MAG
        VBox leftColumn = new VBox(15);
        leftColumn.setMinWidth(430);
        VBox enzymeSection = createEnzymeSection();
        VBox additionalParamsSection = createAdditionalParametersSection();
        
        pfindModificationsPanel = new PFindModificationPanel(this::updateParameterView);
        leftColumn.getChildren().addAll(enzymeSection, additionalParamsSection, pfindModificationsPanel);
        GridPane.setVgrow(leftColumn, Priority.ALWAYS);
        
        // Right column: Quantitative Mode and Metadata Configuration (stacked)
        VBox rightColumn = new VBox(15);
        rightColumn.setMinWidth(430);
        VBox quantModeSection = createQuantitativeModeSection();
        
        rightColumn.getChildren().addAll(quantModeSection, metadataPanel);
        GridPane.setVgrow(rightColumn, Priority.ALWAYS);
        
        // Add columns to grid
        mainGrid.add(leftColumn, 0, 0);
        mainGrid.add(rightColumn, 1, 0);
        
        tabContent.getChildren().add(mainGrid);
        return tabContent;
    }
    
    private VBox createQuantitativeModeSection() {
        VBox section = new VBox(10);
        section.setPadding(new Insets(10));
        section.setStyle("-fx-background-color: #f5f5f5; -fx-background-radius: 5;");
        section.setMinHeight(300);
        VBox.setVgrow(section, Priority.ALWAYS);
        
        Label titleLabel = createSectionTitle("Quantitative Mode");
        quantModeGroup = new ToggleGroup();
        
        // Radio buttons for quantitative modes
        RadioButton labelFreeRadio = new RadioButton("Label-free");
        labelFreeRadio.setToggleGroup(quantModeGroup);
        labelFreeRadio.setSelected(true);
        RadioButton isobaricRadio = new RadioButton("Isobaric");
        isobaricRadio.setToggleGroup(quantModeGroup);
        
        HBox radioBox = new HBox(15);
        radioBox.setAlignment(Pos.CENTER_LEFT);
        radioBox.getChildren().addAll(labelFreeRadio, isobaricRadio);
        
        // StackPane for overlapping panels
        StackPane labelingStackPane = new StackPane();
        labelingStackPane.setMinHeight(200);
        labelingStackPane.setMaxHeight(Double.MAX_VALUE);
        VBox.setVgrow(labelingStackPane, Priority.ALWAYS);
        
        isobaricPanel.setMinHeight(200);
        isobaricPanel.setPrefHeight(220);
        isobaricPanel.setMaxHeight(Double.MAX_VALUE);
        isobaricPanel.setVisible(false);
        isobaricPanel.setManaged(false);
        labelingStackPane.getChildren().add(isobaricPanel);
        
        labelFreeRadio.setOnAction(e -> {
            isobaricPanel.setVisible(false);
            isobaricPanel.setManaged(false);
            updateParameterView();
            updateMetadataTable();
        });
        
        isobaricRadio.setOnAction(e -> {
            isobaricPanel.setVisible(true);
            isobaricPanel.setManaged(true);
            updateParameterView();
            updateMetadataTable();
        });
        
        section.getChildren().addAll(titleLabel, radioBox, labelingStackPane);
        return section;
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
        
        Label ms2Label = new Label("MS2 Scan Mode:");
        ms2Label.setStyle("-fx-font-weight: bold;");
        ms2ScanModeComboBox = new ComboBox<>(FXCollections.observableArrayList("HCD-FTMS", "HCD-ITMS", "CID-FTMS", "CID-ITMS"));
        ms2ScanModeComboBox.setValue("HCD-FTMS");
        ms2ScanModeComboBox.setMinWidth(100);
        ms2ScanModeComboBox.setPrefWidth(120);
        ms2ScanModeComboBox.setMaxWidth(Double.MAX_VALUE);
        ms2ScanModeComboBox.setOnAction(e -> updateParameterView());
        
        VBox ms2Box = new VBox(5);
        ms2Box.getChildren().addAll(ms2Label, ms2ScanModeComboBox);
        
        paramGrid.add(threadsBox, 0, 0);
        paramGrid.add(ms2Box, 1, 0);

        section.getChildren().addAll(titleLabel, paramGrid);
        return section;
    }
    
    
    private void startAnalysis() {
        try {
            consolePanel.appendConsoleText("Starting MAG pFind analysis...\n");
            System.out.println("DEBUG: startAnalysis called - MAG pFind workflow");
            
            // 1. Collect parameters from UI
            consolePanel.appendConsoleText("Collecting parameters from UI...\n");
            DDAMagPfindParameters analysisParams = collectParametersFromUI();
            if (analysisParams == null) {
                consolePanel.appendConsoleText("ERROR: Failed to collect parameters from UI.\n");
                LOGGER.error("Failed to collect parameters from UI in startAnalysis");
                return;
            }
            
            consolePanel.appendConsoleText("✓ Successfully collected parameters from UI\n");
            LOGGER.info("Successfully collected parameters from UI for MAG pFind workflow");
            
            // Store analysis parameters for viewReport
            this.currentAnalysisParams = analysisParams;
            
            // 2. Save parameters to file
            String outputDir = analysisParams.getResult();
            if (outputDir != null && !outputDir.trim().isEmpty()) {
                consolePanel.appendConsoleText("Saving parameters to file...\n");
                boolean paramsSaved = parameterIO.saveAnalysisParameters(outputDir, analysisParams);
                if (!paramsSaved) {
                    consolePanel.appendConsoleText("WARNING: Failed to save parameters to file. Analysis will continue.\n");
                    LOGGER.warn("Failed to save analysis parameters to file: {}", outputDir);
                } else {
                    consolePanel.appendConsoleText("✓ Parameters saved to: " + outputDir + "/mag_pfind_analysis_parameters.json\n");
                    LOGGER.info("Successfully saved analysis parameters to: {}", outputDir);
                }
                
                // Also save as last used parameters
                boolean lastParamsSaved = parameterIO.saveLastUsedParameters(analysisParams);
                if (lastParamsSaved) {
                    consolePanel.appendConsoleText("✓ Last used parameters saved\n");
                    LOGGER.info("Successfully saved last used parameters");
                } else {
                    consolePanel.appendConsoleText("WARNING: Failed to save last used parameters\n");
                    LOGGER.warn("Failed to save last used parameters");
                }
                
                // Save to preferences
                boolean prefsSaved = parameterIO.savePreferences(analysisParams);
                if (prefsSaved) {
                    consolePanel.appendConsoleText("✓ Preferences saved\n");
                    LOGGER.info("Successfully saved preferences");
                } else {
                    consolePanel.appendConsoleText("WARNING: Failed to save preferences\n");
                    LOGGER.warn("Failed to save preferences");
                }
            } else {
                consolePanel.appendConsoleText("WARNING: No output directory specified, skipping parameter save\n");
                LOGGER.warn("No output directory specified, skipping parameter save");
            }
            
            // 3. Create MetaResourceParameters (load from config file)
            consolePanel.appendConsoleText("Loading resource parameters from config...\n");
            MetaResourceIO resourceIO = new MetaResourceIO();
            MetaResourceParameters mrp = resourceIO.loadParameters();
            consolePanel.appendConsoleText("✓ Resource parameters loaded from config file\n");
            
            // 4. Create progress bars and status label
            consolePanel.appendConsoleText("Setting up progress tracking...\n");
            Label statusLabel = new Label("Ready to start analysis...");
            consolePanel.appendConsoleText("✓ Progress bars and status label initialized\n");
            
            // 5. Set status label in parameters
            analysisParams.setStatusLabel(statusLabel);
            consolePanel.appendConsoleText("✓ Status label set in parameters\n");
            
            // 6. Create and run DDAMagPfindTask with DDAMagPfindParameters directly
            consolePanel.appendConsoleText("Creating MAG pFind task...\n");
			currentTask = new DDAMagPfindTask(analysisParams, mrp, consolePanel.getOverallProgressBar(),
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
            LOGGER.info("MAG pFind analysis task started successfully");
            
        } catch (Exception e) {
            String errorMsg = "ERROR: Failed to start analysis: " + e.getMessage();
            consolePanel.appendConsoleText(errorMsg + "\n");
            LOGGER.error("Failed to start MAG pFind analysis", e);
            e.printStackTrace();
        }
    }

    private void stopAnalysis() {
        consolePanel.appendConsoleText("Stopping MAG pFind analysis...\n");
        System.out.println("DEBUG: stopAnalysis called");
        
        if (currentTask != null && currentTask.isRunning()) {
            consolePanel.appendConsoleText("Force stopping MAG pFind task...\n");
            currentTask.forceStop();
            consolePanel.appendConsoleText("MAG pFind task stop requested.\n");
        } else {
            consolePanel.appendConsoleText("No running MAG pFind task to stop.\n");
        }
    }

    private void viewReport() {
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

    private void updateMetadataTableFromIOPanel() {
        if (ioPanel == null || metadataPanel == null) {
            System.out.println("DEBUG: updateMetadataTableFromIOPanel called but some components are null - UI not initialized yet");
            return;
        }
        
        ObservableList<MetaPilotIOPanelFX.ExperimentRow> experiments = ioPanel.getExperimentData();
        metadataPanel.updateMetadataFromIOPanel(experiments);
    }

    private void updateMetadataTable() {
        if (metadataPanel == null || ioPanel == null) {
            System.out.println("DEBUG: updateMetadataTable called but metadata components are null - UI not initialized yet");
            return;
        }
        
        ObservableList<MetaPilotIOPanelFX.ExperimentRow> fileRows = ioPanel.getExperimentData();
        if (fileRows.isEmpty()) return;
        
        if (quantModeGroup != null && quantModeGroup.getSelectedToggle() != null) {
            RadioButton selected = (RadioButton) quantModeGroup.getSelectedToggle();
            String mode = selected.getText();
            boolean hasLabels = isobaricPanel.hasLabels();
            
            metadataPanel.updateMetadataForQuantMode(mode, hasLabels, fileRows, isobaricPanel.getIsobaricTableData());
        }
    }

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
    

    private void updateParameterView() {
        if (consolePanel == null) {
            System.out.println("DEBUG: updateParameterView called but consolePanel is null - UI not initialized yet");
            return;
        }
        
        StringBuilder sb = new StringBuilder();
        sb.append("=== Analysis Parameters ===\n\n");
        
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
        
        // Additional parameters
        sb.append("Additional Parameters:\n");
        if (threadsSpinner != null) {
            sb.append("  Threads: ").append(threadsSpinner.getValue()).append("\n");
        }
        if (ms2ScanModeComboBox != null) {
            sb.append("  MS2 Scan Mode: ").append(ms2ScanModeComboBox.getValue()).append("\n");
        }
        if (generateSampleSpecificDBCheckBox != null) {
            sb.append("  Generate Sample-specific DB: ").append(generateSampleSpecificDBCheckBox.isSelected()).append("\n");
        }
        if (spectraClusteringCheckBox != null) {
            sb.append("  Spectra Clustering: ").append(spectraClusteringCheckBox.isSelected()).append("\n");
        }
        sb.append("\n");
        
        // Quantitative mode
        sb.append("Quantitative Mode:\n");
        if (quantModeGroup != null && quantModeGroup.getSelectedToggle() != null) {
            RadioButton selected = (RadioButton) quantModeGroup.getSelectedToggle();
            sb.append("  Mode: ").append(selected.getText()).append("\n");
            
            if (selected.getText().contains("Isobaric") && isobaricPanel != null) {
                sb.append("  Isobaric Type: ").append(isobaricPanel.getSelectedIsobaricType()).append("\n");
            }
        }
        sb.append("\n");
        
        // Modifications
        sb.append("Modifications:\n");
        sb.append("  Fixed: ").append(pfindModificationsPanel.getFixedModifications().toString()).append("\n");
        sb.append("  Variable: ").append(pfindModificationsPanel.getVariableModifications().toString()).append("\n");
        sb.append("\n");
        
        consolePanel.setParameterText(sb.toString());
    }


    
    private static final Logger LOGGER = LoggerFactory.getLogger(MetaPilotDdaMagPfindPanel.class);
    
    /**
     * Collect all parameters from the UI into a parameters object
     */
    private DDAMagPfindParameters collectParametersFromUI() {
        try {
            LOGGER.info("Starting parameter collection from UI for MAG pFind workflow");
            
            // Validate inputs first
            LOGGER.info("Validating UI inputs...");
            if (!validateInputs()) {
                LOGGER.error("Input validation failed");
                return null;
            }
            LOGGER.info("Input validation passed");
            
            // Collect experiment data from ioPanel
            LOGGER.info("Collecting experiment data from IO panel...");
            ObservableList<ExperimentRow> fileRows = ioPanel.getExperimentData();
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
                ExperimentRow row = fileRows.get(i);
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
            LOGGER.info("Creating MetaData object...");
            MetaData metaData = new MetaData(rawFiles, expNames, fractions, replicates);
            LOGGER.info("MetaData created successfully with {} raw files", rawFiles.length);
            
            // Collect basic parameters
            LOGGER.info("Collecting basic parameters...");
            String result = ioPanel.getOutputDirectory();
            LOGGER.info("Output directory: {}", result);
            
            MagDbItem selectedMagDb = magDatabasePanel.getSelectedMagDatabase();
            String microDb = selectedMagDb != null ? selectedMagDb.getCatalogueID() : "";
            LOGGER.info("Selected MAG database: {}", microDb);
            
            boolean appendHostDb = magDatabasePanel.isAppendHostDb();
            LOGGER.info("Append host DB: {}", appendHostDb);
            
            String hostDB = appendHostDb ? magDatabasePanel.getHostDatabase() : "";
            LOGGER.info("Host database: {}", hostDB);
            
            String enzyme = enzymeComboBox.getValue() != null ? enzymeComboBox.getValue() : "Trypsin";
            LOGGER.info("Enzyme: {}", enzyme);
            
            int missCleavages = missedCleavagesSpinner.getValue() != null ? missedCleavagesSpinner.getValue() : 2;
            LOGGER.info("Missed cleavages: {}", missCleavages);
            
            String digestMode = digestionModeComboBox.getValue() != null ? digestionModeComboBox.getValue() : "Specific";
            LOGGER.info("Digestion mode: {}", digestMode);
            
            int threadCount = threadsSpinner.getValue() != null ? threadsSpinner.getValue() : 4;
            LOGGER.info("Thread count: {}", threadCount);
            
            String ms2ScanMode = ms2ScanModeComboBox.getValue() != null ? ms2ScanModeComboBox.getValue() : "HCD-FTMS";
            LOGGER.info("MS2 scan mode: {}", ms2ScanMode);
            
            boolean reduceDb = generateSampleSpecificDBCheckBox != null ? generateSampleSpecificDBCheckBox.isSelected() : false;
            LOGGER.info("Generate sample-specific DB: {}", reduceDb);
            
            boolean spectraCluster = spectraClusteringCheckBox != null ? spectraClusteringCheckBox.isSelected() : false;
            LOGGER.info("Spectra clustering: {}", spectraCluster);
            
            // Collect isobaric labeling information
            LOGGER.info("Collecting isobaric labeling information...");
            String[] labelTitles = new String[0];
            if (quantModeGroup != null && quantModeGroup.getSelectedToggle() != null) {
                RadioButton selected = (RadioButton) quantModeGroup.getSelectedToggle();
                LOGGER.info("Quantitative mode: {}", selected.getText());
                
                if (selected.getText().contains("Isobaric") && isobaricPanel != null && isobaricPanel.hasLabels()) {
                    ObservableList<IsobaricLabelingPanel.IsobaricLabelRow> isoLabels = isobaricPanel.getIsobaricTableData();
                    labelTitles = new String[isoLabels.size()];
                    for (int i = 0; i < isoLabels.size(); i++) {
                        labelTitles[i] = isoLabels.get(i).getLabelName();
                    }
                    LOGGER.info("Found {} isobaric labels", labelTitles.length);
                }
            }
            
            // Set label titles in metadata
            metaData.setLabelTitle(labelTitles);
            
            // Collect metadata information
            LOGGER.info("Collecting metadata information...");
            ObservableList<MetaPilotDdaIgcMaxQuantPanel.MetadataRow> metadataRows = metadataPanel != null ? metadataPanel.getMetadataTableData() : null;
            int metaTypeCount = metadataPanel != null ? metadataPanel.getMetadataColumnCount() : 0;
            LOGGER.info("Metadata rows: {}, columns: {}", metadataRows != null ? metadataRows.size() : 0, metaTypeCount);
            
            String[][] metaInfo = new String[metadataRows != null ? metadataRows.size() : 0][Math.max(3, metaTypeCount)];
            if (metadataRows != null) {
                for (int i = 0; i < metadataRows.size(); i++) {
                    MetaPilotDdaIgcMaxQuantPanel.MetadataRow row = metadataRows.get(i);
                    metaInfo[i][0] = row.getAge() != null ? row.getAge() : "";
                    metaInfo[i][1] = row.getGender() != null ? row.getGender() : "";
                    metaInfo[i][2] = row.getCondition() != null ? row.getCondition() : "";
                    // Add more metadata fields as needed
                }
            }
            metaData.setMetaInfo(metaInfo);
            
            // Create DDAMagPfindParameters with required constructor parameters
            LOGGER.info("Creating DDAMagPfindParameters object...");
            DDAMagPfindParameters analysisParams = new DDAMagPfindParameters(
                MetaPilotApp.WorkflowType.DDA_MAG_PFIND.getAnalysisType(),
                metaData,
                result,
                microDb,
                appendHostDb,
                hostDB,
                ms2ScanMode,
                threadCount
            );
            LOGGER.info("DDAMagPfindParameters object created successfully");
            
            // Set additional parameters
            LOGGER.info("Setting additional parameters...");
            analysisParams.setEnzyme(enzyme);
            analysisParams.setMissCleavages(missCleavages);
            analysisParams.setDigestMode(digestMode);
            analysisParams.setReduceDb(reduceDb);
            analysisParams.setSpectraCluster(spectraCluster);
            
            // Get MAG database version
            String magDbVersion = "";
            if (selectedMagDb != null && selectedMagDb.getVersionList() != null && !selectedMagDb.getVersionList().isEmpty()) {
                magDbVersion = selectedMagDb.getVersionList().get(0).getCatalogueVersion();
            }
            LOGGER.info("MAG database version: {}", magDbVersion);
            
            // Set MAG-specific fields
            analysisParams.setMagDbCatalogueId(microDb);
            analysisParams.setMagDbVersion(magDbVersion);
            
            // Set the used version on the MagDbItem
            MagVersion selectedVersion = magDatabasePanel.getSelectedMagVersion();
            if (selectedVersion != null && selectedMagDb != null) {
                selectedMagDb.setUsedVersion(selectedVersion);
                LOGGER.info("Set used version on MagDbItem: {}", selectedVersion.getCatalogueVersion());
            }
            
            analysisParams.setMagDbItem(selectedMagDb); // Store the complete MagDbItem
            
            // Get quantitative mode
            String quanMode = "Label-Free";
            if (quantModeGroup != null && quantModeGroup.getSelectedToggle() != null) {
                RadioButton selected = (RadioButton) quantModeGroup.getSelectedToggle();
                quanMode = selected.getText();
            }
            LOGGER.info("Quantitative mode: {}", quanMode);
            
            // Set modifications
            analysisParams.setFixedModifications(new ArrayList<>(pfindModificationsPanel.getFixedModifications()));
            analysisParams.setVariableModifications(new ArrayList<>(pfindModificationsPanel.getVariableModifications()));
            LOGGER.info("Modifications set - Fixed: {}, Variable: {}", 
                pfindModificationsPanel.getFixedModifications().size(), pfindModificationsPanel.getVariableModifications().size());
            
            // Set quantitative mode
            analysisParams.setQuantitativeMode(quanMode);
            
            // Set isobaric labeling
            analysisParams.setIsobaricType(quanMode.contains("Isobaric") ? quanMode : "");
            
            // Set isotopic labeling
            analysisParams.setLightLabels(new ArrayList<>()); // TODO: Get from isotopic labeling panel
            analysisParams.setMediumLabels(new ArrayList<>()); // TODO: Get from isotopic labeling panel
            analysisParams.setHeavyLabels(new ArrayList<>()); // TODO: Get from isotopic labeling panel
            
            // Set experiment data
            analysisParams.setExperimentData(new ArrayList<>(fileRows));
            
            // Set metadata table data
            if (metadataRows != null) {
                analysisParams.setMetadataTableData(new ArrayList<>(metadataRows));
            }
            
            LOGGER.info("Successfully collected parameters from UI for MAG pFind workflow");
            return analysisParams;
            
        } catch (Exception e) {
            LOGGER.error("Failed to collect parameters from UI: " + e.getMessage(), e);
            return null;
        }
    }
    
    /**
     * Convert digestion mode string to integer representation
     */
    private int getDigestModeInt(String digestMode) {
        switch (digestMode) {
            case "Specific": return 1;
            case "Semi free N-term": return 2;
            case "Semi free C-term": return 3;
            case "Semi specific": return 4;
            case "Unspecific": return 5;
            case "None": return 6;
            default: return 1;
        }
    }
    
    /**
     * Validate inputs before starting analysis
     */
    private boolean validateInputs() {
        LOGGER.info("Starting input validation for MAG pFind workflow");
        
        // Check experiment data
        if (ioPanel == null) {
            LOGGER.error("IO panel is null");
            showError("Validation Error", "IO panel is not initialized.");
            return false;
        }
        
        if (ioPanel.getExperimentData().isEmpty()) {
            LOGGER.error("No experiment data found in IO panel");
            showError("Validation Error", "Please add at least one experiment.");
            return false;
        }
        LOGGER.info("Experiment data validation passed: {} experiments found", ioPanel.getExperimentData().size());
        
        // Check output directory
        String outputDir = ioPanel.getOutputDirectory();
        if (outputDir == null || outputDir.trim().isEmpty()) {
            LOGGER.error("Output directory is empty or null");
            showError("Validation Error", "Please select an output directory.");
            return false;
        }
        LOGGER.info("Output directory validation passed: {}", outputDir);
        
        // Check MAG database selection using reusable component
        if (magDatabasePanel == null) {
            LOGGER.error("MAG database panel is null");
            showError("Validation Error", "MAG database panel is not initialized.");
            return false;
        }
        
        if (!magDatabasePanel.validateSelection()) {
            String errorMsg = magDatabasePanel.getValidationErrorMessage();
            LOGGER.error("MAG database validation failed: {}", errorMsg);
            showError("Validation Error", errorMsg);
            return false;
        }
        
        MagDbItem selectedMagDb = magDatabasePanel.getSelectedMagDatabase();
        MagVersion selectedVersion = magDatabasePanel.getSelectedMagVersion();
        LOGGER.info("MAG database validation passed: {} (version: {})", 
            selectedMagDb.getCatalogueID(), selectedVersion.getCatalogueVersion());
        
        // Check enzyme selection
        if (enzymeComboBox == null) {
            LOGGER.error("Enzyme combo box is null");
            showError("Validation Error", "Enzyme selection is not initialized.");
            return false;
        }
        
        String selectedEnzyme = enzymeComboBox.getValue();
        if (selectedEnzyme == null) {
            LOGGER.error("No enzyme selected");
            showError("Validation Error", "Please select an enzyme.");
            return false;
        }
        LOGGER.info("Enzyme validation passed: {}", selectedEnzyme);
        
        LOGGER.info("All input validation checks passed successfully");
        return true;
    }
    
    /**
     * Show error dialog
     */
    private void showError(String title, String content) {
        Platform.runLater(() -> {
            Alert alert = new Alert(Alert.AlertType.ERROR);
            alert.setTitle(title);
            alert.setHeaderText(null);
            alert.setContentText(content);
            alert.showAndWait();
        });
    }
    
    /**
     * Save parameters to a JSON file
     */
    private boolean saveParametersToFile(DDAMagPfindParameters parameters, String outputDir) {
        try {
            if (parameters == null || outputDir == null || outputDir.trim().isEmpty()) {
                LOGGER.warn("Cannot save parameters: parameters or output directory is null/empty");
                return false;
            }
            
            File outputFile = new File(outputDir, "mag_pfind_parameters.json");
            JSONObject json = new JSONObject();
            
            // Header information
            json.put("workflow", MetaPilotApp.WorkflowType.DDA_MAG_PFIND.getAnalysisType());
            json.put("version", "1.0");
            json.put("generated_at", new SimpleDateFormat("yyyy-MM-dd HH:mm:ss").format(new Date()));
            
            // Basic parameters
            json.put("output_directory", outputDir);
            json.put("mag_database", parameters.getMicroDb());
            json.put("mag_database_version", parameters.getMagDbVersion());
            json.put("append_host_db", parameters.isAppendHostDb());
            json.put("host_database", parameters.getHostDB());
            
            // Enzyme settings
            json.put("enzyme", parameters.getEnzyme());
            json.put("missed_cleavages", parameters.getMissCleavages());
            json.put("digestion_mode", parameters.getDigestMode());
            
            // Additional parameters
            json.put("threads", parameters.getThreadCount());
            json.put("ms2_scan_mode", parameters.getMs2ScanMode());
            json.put("generate_sample_specific_db", parameters.isReduceDb());
            
            // Quantitative mode
            json.put("quantitative_mode", parameters.getQuantitativeMode());
            
            // Modifications
            JSONArray fixedMods = new JSONArray();
            if (parameters.getFixMods() != null) {
                for (String mod : parameters.getFixMods()) {
                    fixedMods.put(mod);
                }
            }
            json.put("fixed_modifications", fixedMods);
            
            JSONArray varMods = new JSONArray();
            if (parameters.getVariMods() != null) {
                for (String mod : parameters.getVariMods()) {
                    varMods.put(mod);
                }
            }
            json.put("variable_modifications", varMods);
            
            // Experiment data
            if (parameters.getMetadata() != null) {
                JSONArray experiments = new JSONArray();
                String[] rawFiles = parameters.getMetadata().getRawFiles();
                String[] expNames = parameters.getMetadata().getExpNames();
                if (rawFiles != null && expNames != null) {
                    for (int i = 0; i < rawFiles.length; i++) {
                        JSONObject exp = new JSONObject();
                        exp.put("raw_file", rawFiles[i]);
                        exp.put("experiment_name", i < expNames.length ? expNames[i] : "");
                        experiments.put(exp);
                    }
                }
                json.put("experiments", experiments);
            }
            
            // Write to file
            try (FileWriter writer = new FileWriter(outputFile)) {
                writer.write(json.toString(2)); // Pretty print with 2 spaces
            }
            
            LOGGER.info("Successfully saved MAG pFind parameters to: " + outputFile.getAbsolutePath());
            return true;
            
        } catch (Exception e) {
            LOGGER.error("Failed to save parameters to file: " + e.getMessage(), e);
            return false;
        }
    }
    
    /**
     * Load last used parameters and apply them to the UI
     */
    private void loadLastUsedParameters() {
        try {
            DDAMagPfindParameters params = (DDAMagPfindParameters) parameterIO.loadLastUsedParameters();
            if (params != null) {
                applyParametersToUI(params);
                LOGGER.info("Loaded last used parameters for DDA MAG pFind workflow");
                
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
                LOGGER.warn("No parameters loaded, using defaults");
            }
        } catch (Exception e) {
            LOGGER.warn("Failed to load last used MAG pFind parameters: " + e.getMessage(), e);
            // Continue with default UI state
        }
    }
    
    /**
     * Apply parameters to the UI components
     */
    private void applyParametersToUI(DDAMagPfindParameters params) {
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
        
        // Apply basic settings that don't depend on file paths
        if (params.getEnzyme() != null && enzymeComboBox != null) {
            enzymeComboBox.setValue(params.getEnzyme());
        }
        if (missedCleavagesSpinner != null) {
            missedCleavagesSpinner.getValueFactory().setValue(params.getMissCleavages());
        }
        if (params.getDigestModeString() != null && digestionModeComboBox != null) {
            digestionModeComboBox.setValue(params.getDigestModeString());
        }
        if (threadsSpinner != null) {
            threadsSpinner.getValueFactory().setValue(params.getThreadCount());
        }
        if (params.getMs2ScanMode() != null && ms2ScanModeComboBox != null) {
            ms2ScanModeComboBox.setValue(params.getMs2ScanMode());
        }
        if (generateSampleSpecificDBCheckBox != null) {
            generateSampleSpecificDBCheckBox.setSelected(params.isReduceDb());
        }
        if (spectraClusteringCheckBox != null) {
            spectraClusteringCheckBox.setSelected(params.isSpectraCluster());
        }
        
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
        if (pfindModificationsPanel != null) {
            pfindModificationsPanel.setFixedModifications(params.getFixedModifications());
            pfindModificationsPanel.setVariableModifications(params.getVariableModifications());
            LOGGER.info("Applied modifications to pFind modifications panel");
        }
        
        // Apply quantitative mode settings
        if (params.getQuantitativeMode() != null && quantModeGroup != null) {
            String mode = params.getQuantitativeMode();
            if (mode.contains("Isobaric")) {
                // Find and select the isobaric radio button
                for (Toggle toggle : quantModeGroup.getToggles()) {
                    if (toggle instanceof RadioButton && ((RadioButton) toggle).getText().contains("Isobaric")) {
                        quantModeGroup.selectToggle(toggle);
                        break;
                    }
                }
            } else {
                // Find and select the label-free radio button
                for (Toggle toggle : quantModeGroup.getToggles()) {
                    if (toggle instanceof RadioButton && ((RadioButton) toggle).getText().contains("Label-free")) {
                        quantModeGroup.selectToggle(toggle);
                        break;
                    }
                }
            }
        }
        
        // Apply isobaric type if available
        if (params.getIsobaricType() != null && isobaricPanel != null) {
            isobaricPanel.setSelectedIsobaricType(params.getIsobaricType());
            LOGGER.info("Applied isobaric type: {}", params.getIsobaricType());
        }
        
        // Apply isobaric table data if available
        if (params.getIsobaricTableData() != null && !params.getIsobaricTableData().isEmpty() && isobaricPanel != null) {
            isobaricPanel.setIsobaricTableData(params.getIsobaricTableData());
            LOGGER.info("Applied {} isobaric table rows", params.getIsobaricTableData().size());
        }
        
        // Apply metadata table data if available
        if (params.getMetadataTableData() != null && !params.getMetadataTableData().isEmpty() && metadataPanel != null) {
            metadataPanel.setMetadataTableData(params.getMetadataTableData());
            LOGGER.info("Applied {} metadata table rows", params.getMetadataTableData().size());
        }
        
        // Update parameter view after applying parameters
        Platform.runLater(() -> updateParameterView());
    }
    
    /**
     * Load parameters from a specific analysis result directory
     */
    public void loadAnalysisParameters(String resultDirectory) {
        try {
            DDAMagPfindParameters params = (DDAMagPfindParameters) parameterIO.loadAnalysisParameters(resultDirectory);
            if (params != null) {
                applyParametersToUI(params);
                showError("Success", "Loaded MAG pFind analysis parameters from " + resultDirectory);
                LOGGER.info("Loaded MAG pFind analysis parameters from {}", resultDirectory);
            } else {
                showError("Warning", "No parameters found in " + resultDirectory);
                LOGGER.warn("No parameters found in {}", resultDirectory);
            }
        } catch (Exception e) {
            showError("Error", "Failed to load MAG pFind analysis parameters: " + e.getMessage());
            LOGGER.error("Failed to load MAG pFind analysis parameters from " + resultDirectory, e);
        }
    }

    private void setupWindowCloseHandler() {
        Platform.runLater(() -> {
            try {
                // Get the current scene and window
                if (getScene() != null) {
                    javafx.stage.Window window = getScene().getWindow();
                    if (window != null) {
                        window.setOnCloseRequest(event -> {
                            try {
                                consolePanel.appendConsoleText("Saving parameters on window close...\n");
                                DDAMagPfindParameters currentParams = collectParametersFromUI();
                                if (currentParams != null) {
                                    boolean saved = parameterIO.saveLastUsedParameters(currentParams);
                                    if (saved) {
                                        consolePanel.appendConsoleText("✓ Parameters saved to config folder: last_used_mag_pfind_parameters.json\n");
                                        LOGGER.info("Successfully saved parameters on window close to config folder");
                                    } else {
                                        consolePanel.appendConsoleText("WARNING: Failed to save parameters to config folder\n");
                                        LOGGER.warn("Failed to save parameters to config folder on window close");
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
                                    LOGGER.warn("No parameters to save on window close");
                                }
                            } catch (Exception e) {
                                String errorMsg = "ERROR: Failed to save parameters on window close: " + e.getMessage();
                                consolePanel.appendConsoleText(errorMsg + "\n");
                                LOGGER.error("Failed to save parameters on window close", e);
                            }
                        });
                        LOGGER.info("Window close handler set up successfully");
                    } else {
                        LOGGER.warn("Window is null, cannot set close handler");
                    }
                } else {
                    LOGGER.warn("Scene is null, cannot set close handler");
                }
            } catch (Exception e) {
                LOGGER.error("Failed to set up window close handler", e);
            }
        });
    }
} 