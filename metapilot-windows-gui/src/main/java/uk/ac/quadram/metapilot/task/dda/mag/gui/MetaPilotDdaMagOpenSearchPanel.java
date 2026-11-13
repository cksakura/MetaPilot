package uk.ac.quadram.metapilot.task.dda.mag.gui;

import uk.ac.quadram.metapilot.core.msfragger.MSFraggerHeadlessRunner;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Optional;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import javafx.application.Platform;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.geometry.Insets;
import javafx.scene.control.Alert;
import javafx.scene.control.Alert.AlertType;
import javafx.scene.control.ButtonType;
import javafx.scene.control.CheckBox;
import javafx.scene.control.ComboBox;
import javafx.scene.control.Label;
import javafx.scene.control.Spinner;
import javafx.scene.control.Tab;
import javafx.scene.control.TabPane;
import javafx.scene.control.TitledPane;
import javafx.scene.control.Tooltip;
import javafx.scene.layout.ColumnConstraints;
import javafx.scene.layout.GridPane;
import javafx.scene.layout.Priority;
import javafx.scene.layout.Region;
import javafx.scene.layout.VBox;
import javafx.scene.paint.Color;
import javafx.scene.text.Font;
import javafx.scene.text.FontWeight;

import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.task.common.gui.AbstractWorkflowPanel;
import uk.ac.quadram.metapilot.task.common.gui.ConsolePanel;
import uk.ac.quadram.metapilot.task.common.gui.MagDatabaseSettingsPanel;
import uk.ac.quadram.metapilot.task.common.gui.MetaPilotIOPanelFX;
import uk.ac.quadram.metapilot.task.common.gui.MetadataConfigurationPanel;
import uk.ac.quadram.metapilot.task.common.gui.ParameterSection;
import uk.ac.quadram.metapilot.task.common.gui.ModificationsPanel;
import uk.ac.quadram.metapilot.task.dda.igc.gui.MetaPilotDdaIgcMaxQuantPanel;
import uk.ac.quadram.metapilot.task.dda.mag.par.DDAMagOpenSearchParameterIO;
import uk.ac.quadram.metapilot.task.dda.mag.par.DDAMagOpenSearchParameters;
import uk.ac.quadram.metapilot.task.mag.GenomeSelectionCut;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;
import uk.ac.quadram.metapilot.task.mag.MagDbItem.MagVersion;
import uk.ac.quadram.metapilot.task.par.MetaData;

/**
 * Configuration panel for the DDA MAG open-search workflow.
 * <p>
 * Follows the same three-tab structure as the other MAG workflow panels:
 * <b>Input/Output</b> and <b>Console</b> are the shared components used
 * unchanged, and only the middle <b>Parameters</b> tab is specific to this
 * workflow. Extending {@link AbstractWorkflowPanel} keeps parameter
 * persistence, console wiring and the auto-save-on-Console-tab behaviour
 * consistent with pFind, FragPipe and Sage.
 * <p>
 * Note the Parameters tab has no search-engine selector, Closed/Open selector
 * and no rank-depth, peak-fitting or local-FDR knobs: those are locked to the
 * tested defaults. The engine is MSFragger. Round one is a closed marker
 * search; identification is joint closed+open.
 */
public class MetaPilotDdaMagOpenSearchPanel extends AbstractWorkflowPanel<DDAMagOpenSearchParameters> {

	private static final Logger LOGGER = LogManager.getLogger(MetaPilotDdaMagOpenSearchPanel.class);

	// Shared components (Input/Output tab)
	private MetaPilotIOPanelFX ioPanel;
	private MagDatabaseSettingsPanel magDatabasePanel;
	private MetadataConfigurationPanel metadataPanel;
	private ModificationsPanel modificationsPanel;
	private final ObservableList<String> enzymes = FXCollections.observableArrayList();

	private final MetaResourceParameters resourceParams;

	// Workflow-specific controls (Parameters tab). Rank depth, the open
	// precursor window, round-1 mode, the closed stage, and peak-fitting /
	// local-FDR internals are locked to the tested defaults and are not shown.
	private ComboBox<String> enzymeComboBox;
	private Spinner<Integer> missedCleavagesSpinner;
	private ComboBox<String> genomeSelectorComboBox;
	private GenomeStopControls genomeStopControls;
	private ComboBox<Integer> threadsComboBox;
	private ComboBox<String> heapComboBox;
	private CheckBox entrapmentCheckBox;

	public MetaPilotDdaMagOpenSearchPanel(MagDbItem[] magDbItems, MetaResourceParameters resourceParams) {
		super();
		this.resourceParams = resourceParams;
		try {
			parameterIO = new DDAMagOpenSearchParameterIO();
			loadEnzymeChoices();
			initializeSpecificUI();
			setAvailableMagDatabases(magDbItems);
			LOGGER.info("[DDA MAG Open Search] Panel constructed");
		} catch (Exception e) {
			LOGGER.error("[DDA MAG Open Search] Constructor failed", e);
			throw e;
		}
	}

	private void initializeSpecificUI() {
		setSpacing(10);
		setPadding(new Insets(10));
		setMinWidth(1200);
		setMinHeight(800);
		setPrefWidth(1600);
		setPrefHeight(900);

		consolePanel = new ConsolePanel();
		metadataPanel = new MetadataConfigurationPanel(true);
		magDatabasePanel = new MagDatabaseSettingsPanel();

		metadataPanel.setParameterUpdateCallback(this::updateParameterView);
		magDatabasePanel.setParameterUpdateCallback(this::updateParameterView);
		configureConsolePanel(consolePanel);

		TabPane tabPane = new TabPane();
		tabPane.setTabClosingPolicy(TabPane.TabClosingPolicy.UNAVAILABLE);
		tabPane.setMinWidth(980);
		tabPane.setMinHeight(780);
		VBox.setVgrow(tabPane, Priority.ALWAYS);

		Tab ioTab = new Tab("Input/Output");
		ioTab.setContent(createIOTab());

		Tab paramTab = new Tab("Parameters");
		paramTab.setContent(createParametersTab());

		Tab consoleTab = new Tab("Console");
		consoleTab.setContent(consolePanel);

		tabPane.getTabs().addAll(ioTab, paramTab, consoleTab);

		// Same auto-save-on-Console behaviour as the other MAG panels, so a run
		// started from any workflow leaves an equivalent parameter record.
		tabPane.getSelectionModel().selectedItemProperty().addListener((obs, oldTab, newTab) -> {
			if (newTab == paramTab) {
				updateMetadataTableFromIOPanel();
			} else if (newTab == consoleTab) {
				updateParameterView();
				try {
					DDAMagOpenSearchParameters current = collectParametersFromUI();
					if (current != null) {
						if (parameterIO.saveLastUsedParameters(current)) {
							consolePanel.appendConsoleText("Parameters auto-saved to config folder\n");
						} else {
							consolePanel.appendConsoleText("WARNING: failed to auto-save parameters\n");
						}
						parameterIO.savePreferences(current);
					}
				} catch (Exception e) {
					LOGGER.warn("Auto-save on Console tab failed", e);
					consolePanel.appendConsoleText("WARNING: auto-save failed: " + e.getMessage() + "\n");
				}
			}
		});

		getChildren().add(tabPane);

		// Deferred to runLater like the other MAG panels: the controls exist by
		// now, but applying parameters touches the scene, and this also lets
		// loadParametersFromObject(...) set parametersLoadedFromFile first when
		// the panel was opened from a saved parameter file.
		Platform.runLater(() -> {
			loadLastUsedParameters();
			setupWindowCloseHandler();
		});
	}

	// ------------------------------------------------------------- shared tabs

	/**
	 * Identical to the other MAG workflows: the shared I/O panel plus the shared
	 * MAG database settings panel.
	 */
	@Override
	public VBox createIOTab() {
		VBox tabContent = new VBox(10);
		tabContent.setPadding(new Insets(10));
		ioPanel = new MetaPilotIOPanelFX();
		tabContent.getChildren().addAll(createSectionTitle("Input/Output Settings"), ioPanel,
				createSectionTitle("Database Settings"), magDatabasePanel);
		return tabContent;
	}

	/**
	 * Enabled now that {@code mag-opensearch} (and {@code run} for
	 * {@code DDA_MAG_OPENSEARCH}) exist in the headless CLI.
	 */
	@Override
	protected boolean isHpcExportSupported() {
		return true;
	}

	public void setAvailableMagDatabases(MagDbItem[] magDbItems) {
		if (magDatabasePanel != null) {
			magDatabasePanel.setAvailableMagDatabases(magDbItems);
		}
	}

	// ------------------------------------------------ workflow-specific tab

	/**
	 * Same two-column layout the other MAG panels use: workflow settings and
	 * modifications on the left, metadata configuration on the right.
	 */
	@Override
	public VBox createParametersTab() {
		VBox tabContent = new VBox(10);
		tabContent.setPadding(new Insets(10));
		tabContent.setMinWidth(960);

		GridPane mainGrid = new GridPane();
		mainGrid.setHgap(15);
		mainGrid.setVgap(15);
		mainGrid.setPadding(new Insets(10));
		mainGrid.setMinWidth(940);

		ColumnConstraints col1 = new ColumnConstraints();
		col1.setPercentWidth(50);
		col1.setMinWidth(450);
		col1.setHgrow(Priority.ALWAYS);
		ColumnConstraints col2 = new ColumnConstraints();
		col2.setPercentWidth(50);
		col2.setMinWidth(450);
		col2.setHgrow(Priority.ALWAYS);
		mainGrid.getColumnConstraints().addAll(col1, col2);

		VBox leftColumn = new VBox(15);
		leftColumn.setMinWidth(430);
		modificationsPanel = new ModificationsPanel(this::updateParameterView);
		leftColumn.getChildren().addAll(createSearchSection(), createGenomeSection(), modificationsPanel);
		GridPane.setVgrow(leftColumn, Priority.ALWAYS);

		VBox rightColumn = new VBox(15);
		rightColumn.setMinWidth(430);
		rightColumn.getChildren().addAll(metadataPanel);
		GridPane.setVgrow(rightColumn, Priority.ALWAYS);

		mainGrid.add(leftColumn, 0, 0);
		mainGrid.add(rightColumn, 1, 0);

		tabContent.getChildren().add(ParameterSection.createViewport(mainGrid));
		return tabContent;
	}

	/**
	 * Enzyme names come from the bundled {@code enzyme_default.ini}, the same
	 * source the other panels read, rather than a hardcoded list that would
	 * drift from it.
	 */
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
	 * The pFind panel still reads the file, correctly: there its semantics are the
	 * engine's own.
	 */
	private void loadEnzymeChoices() {
		enzymes.setAll(MSFraggerHeadlessRunner.ENZYME_CHOICES);
	}

	private void updateMetadataTableFromIOPanel() {
		if (ioPanel == null || metadataPanel == null) {
			return;
		}
		metadataPanel.updateMetadataFromIOPanel(ioPanel.getExperimentData());
	}

	private TitledPane createSearchSection() {
		GridPane grid = paramGrid();
		int row = 0;

		grid.add(new Label("Enzyme:"), 0, row);
		enzymeComboBox = new ComboBox<>(enzymes);
		if (!enzymes.isEmpty()) {
			enzymeComboBox.setValue(enzymes.get(0));
		}
		enzymeComboBox.setMaxWidth(Double.MAX_VALUE);
		grid.add(enzymeComboBox, 1, row++);

		grid.add(new Label("Missed cleavages:"), 0, row);
		missedCleavagesSpinner = new Spinner<>(0, 5, 2);
		missedCleavagesSpinner.setMaxWidth(Double.MAX_VALUE);
		grid.add(missedCleavagesSpinner, 1, row++);

		grid.add(new Label("Threads:"), 0, row);
		threadsComboBox = new ComboBox<>();
		for (int i = 1; i <= Runtime.getRuntime().availableProcessors(); i++) {
			threadsComboBox.getItems().add(i);
		}
		threadsComboBox.getSelectionModel().select(Integer.valueOf(Runtime.getRuntime().availableProcessors()));
		threadsComboBox.setMaxWidth(Double.MAX_VALUE);
		grid.add(threadsComboBox, 1, row++);

		grid.add(new Label("MSFragger memory:"), 0, row);
		heapComboBox = new ComboBox<>();
		heapComboBox.getItems().addAll("8g", "16g", "32g", "64g", "128g");
		heapComboBox.getSelectionModel().select("32g");
		heapComboBox.setEditable(true);
		heapComboBox.setMaxWidth(Double.MAX_VALUE);
		heapComboBox.setTooltip(new Tooltip("Java -Xmx value. Open searches are memory hungry."));
		grid.add(heapComboBox, 1, row++);

		Label methodHint = new Label("Search method is fixed to the tested defaults: closed marker search "
				+ "(one candidate per spectrum), then joint closed+open identification of the selected genomes "
				+ "(three candidates per spectrum). Rank depth, precursor window, peak fitting and FDR internals "
				+ "are not shown; they can still be set in the parameter JSON for CLI and HPC.");
		methodHint.setWrapText(true);
		methodHint.setMaxWidth(Double.MAX_VALUE);
		methodHint.setMinHeight(Region.USE_PREF_SIZE);
		methodHint.setStyle("-fx-font-size: 11px; -fx-text-fill: #666;");
		grid.add(methodHint, 0, row++, 2, 1);

		entrapmentCheckBox = new CheckBox("Paired-shuffle entrapment (evaluation)");
		entrapmentCheckBox.setWrapText(true);
		entrapmentCheckBox.setTooltip(new Tooltip(
				"Same processing as mag-opensearch-entrapment: after genome selection, each sample database also contains a shuffled copy of the selected proteins.\n"
						+ "Results go to opensearch_entrapment_MAG_search/ and stop after identification.\n"
						+ "Requires Trypsin or Trypsin/P. Off by default."));
		entrapmentCheckBox.setOnAction(e -> updateParameterView());
		grid.add(entrapmentCheckBox, 0, row++, 2, 1);

		return ParameterSection.create("Search", grid);
	}

	private TitledPane createGenomeSection() {
		GridPane grid = paramGrid();
		int row = 0;

		grid.add(new Label("Marker set"), 0, row);
		genomeSelectorComboBox = new ComboBox<>();
		genomeSelectorComboBox.getItems().addAll(MagDbItem.MARKER2, MagDbItem.MARKER4);
		genomeSelectorComboBox.getSelectionModel().select(MagDbItem.MARKER4);
		genomeSelectorComboBox.setMaxWidth(Double.MAX_VALUE);
		genomeSelectorComboBox.setOnAction(e -> updateParameterView());
		grid.add(genomeSelectorComboBox, 1, row++);

		genomeStopControls = new GenomeStopControls(this::updateParameterView);
		grid.add(genomeStopControls, 0, row++, 2, 1);

		return ParameterSection.create("Genome discovery", grid);
	}

	private static GridPane paramGrid() {
		GridPane grid = new GridPane();
		grid.setHgap(10);
		grid.setVgap(8);
		return grid;
	}

	private Label createSectionTitle(String title) {
		Label label = new Label(title);
		label.setFont(Font.font("System", FontWeight.BOLD, 14));
		label.setTextFill(Color.DARKBLUE);
		return label;
	}

	// ------------------------------------------------------ parameter binding

	@Override
	protected DDAMagOpenSearchParameters collectParametersFromUI() {
		DDAMagOpenSearchParameters params = currentAnalysisParams == null ? new DDAMagOpenSearchParameters()
				: currentAnalysisParams;

		if (ioPanel != null) {
			params.setResult(ioPanel.getOutputDirectory());

			List<MetaPilotIOPanelFX.ExperimentRow> fileRows = new ArrayList<>();
			for (MetaPilotIOPanelFX.ExperimentRow row : ioPanel.getExperimentData()) {
				if (row.getFilePath() != null && !row.getFilePath().isBlank()) {
					fileRows.add(row);
				}
			}

			String[] rawFiles = new String[fileRows.size()];
			String[] expNames = new String[fileRows.size()];
			int[] fractions = new int[fileRows.size()];
			int[] replicates = new int[fileRows.size()];
			List<String[]> experiments = new ArrayList<>();

			for (int i = 0; i < fileRows.size(); i++) {
				MetaPilotIOPanelFX.ExperimentRow row = fileRows.get(i);
				rawFiles[i] = row.getFilePath();
				expNames[i] = row.getExperimentName();
				try {
					fractions[i] = Integer.parseInt(row.getFraction());
					if (fractions[i] <= 0) {
						fractions[i] = 1;
					}
				} catch (NumberFormatException e) {
					fractions[i] = 1;
				}
				replicates[i] = 1;
				experiments.add(new String[] { rawFiles[i], expNames[i], String.valueOf(fractions[i]) });
			}

			// MetaData is what the task hierarchy actually reads -
			// MetaPilotMagTask.createSeparateFolders() throws "Raw files are
			// required" without it. Setting only experimentData, as an earlier
			// version did, left it empty and the task failed in its constructor.
			params.setMetadata(new MetaData(rawFiles, expNames, fractions, replicates));
			params.setExperimentData(experiments);
			LOGGER.info("Collected {} raw files into MetaData", rawFiles.length);
		}

		if (magDatabasePanel != null) {
			MagDbItem selected = magDatabasePanel.getSelectedMagDatabase();
			MagVersion selectedVersion = magDatabasePanel.getSelectedMagVersion();

			String microDb = selected != null ? selected.getCatalogueID() : "";
			params.setMicroDb(microDb);
			params.setMagDbCatalogueId(microDb);

			if (selected != null && selected.getIdentifier() != null && !selected.getIdentifier().isBlank()) {
				params.setMagDbIdentifier(selected.getIdentifier());
			}

			// Version matters as much as the catalogue id: without it the reload
			// cannot restore which version was chosen.
			if (selectedVersion != null) {
				params.setMagDbVersion(selectedVersion.getCatalogueVersion());
			} else if (selected != null && selected.getVersionList() != null
					&& !selected.getVersionList().isEmpty()) {
				// Fallback for when the version combo has not initialised yet.
				params.setMagDbVersion(selected.getVersionList().get(0).getCatalogueVersion());
			}

			if (selectedVersion != null && selected != null) {
				selected.setUsedVersion(selectedVersion);
			}

			params.setMagDbItem(selected);
			params.setAppendHostDb(magDatabasePanel.isAppendHostDb());
			params.setHostDB(magDatabasePanel.getHostDatabase());
		}

		// Modifications, in the same string format the other MAG panels persist:
		// fixed as "site:mass", variable as "name:sites:mass".
		if (modificationsPanel != null) {
			List<String> fixedMods = new ArrayList<>();
			for (ModificationsPanel.FixedModificationRow row : modificationsPanel.getFixedModData()) {
				if (row.getDeltaMass() != 0.0) {
					fixedMods.add(row.getSite() + ":" + row.getDeltaMass());
				}
			}
			List<String> variableMods = new ArrayList<>();
			for (ModificationsPanel.VariableModificationRow row : modificationsPanel.getVariableModData()) {
				if (row.getUse()) {
					variableMods.add(row.getModName() + ":" + row.getSites() + ":" + row.getDeltaMass());
				}
			}
			params.setFixMods(fixedMods.toArray(new String[0]));
			params.setVariMods(variableMods.toArray(new String[0]));
		}

		if (metadataPanel != null) {
			List<String[]> metadataArrays = new ArrayList<>();
			for (MetaPilotDdaIgcMaxQuantPanel.MetadataRow row : metadataPanel.getMetadataTableData()) {
				metadataArrays.add(row.toArray());
			}
			params.setMetadataTableData(metadataArrays);
			metadataPanel.applyTo(params.getMetadata());
		}

		params.setEnzyme(enzymeComboBox.getValue());
		params.setMissCleavages(missedCleavagesSpinner.getValue());
		params.setGenomeSelectorMode(genomeSelectorComboBox.getValue());
		if (genomeStopControls != null) {
			params.setGenomeSelectionStrategy(genomeStopControls.strategy());
			params.setGenomeSelectionPsmCoverage(genomeStopControls.psmCoverage());
			params.setGenomeSelectionStepGuard(genomeStopControls.stepGuard());
			params.setGenomeSelectionMaxGenomes(genomeStopControls.maxGenomes());
		}
		params.setThreadCount(threadsComboBox.getValue() == null ? 1 : threadsComboBox.getValue());
		params.setMsfraggerMaxHeap(heapValue());
		// Rank depth, precursor window, round-1 mode, closed stage and
		// peak-fitting internals are not GUI settings.
		params.applyGuiLockedDefaults();
		params.setEntrapmentEnabled(entrapmentCheckBox != null && entrapmentCheckBox.isSelected());

		currentAnalysisParams = params;
		return params;
	}

	@Override
	protected void applyParametersToUI(DDAMagOpenSearchParameters params) {
		if (params == null) {
			return;
		}
		currentAnalysisParams = params;

		if (ioPanel != null) {
			if (params.getResult() != null) {
				ioPanel.setOutputDirectory(params.getResult());
			}
			if (params.getExperimentData() != null && !params.getExperimentData().isEmpty()) {
				List<MetaPilotIOPanelFX.ExperimentRow> rows = new ArrayList<>();
				for (String[] experiment : params.getExperimentData()) {
					rows.add(MetaPilotIOPanelFX.ExperimentRow.fromArray(experiment));
				}
				ioPanel.setExperimentData(rows);
			}
		}

		if (magDatabasePanel != null) {
			// Restore the saved MAG database selection. Without this the combo
			// keeps whatever setAvailableMagDatabases(...) selected first, which
			// is why a saved choice appeared to be ignored on every launch.
			// magDbIdentifier has no UI control of its own - it rides on the
			// params object (retained via currentAnalysisParams) and is refreshed
			// from the selected item on the next collect.
			String targetMagDbId = params.getMagDbCatalogueId();
			if (targetMagDbId == null || targetMagDbId.trim().isEmpty()) {
				// Older saved parameters kept the MAG DB in MetaBaseParameter.microDb.
				targetMagDbId = params.getMicroDb();
			}
			if (targetMagDbId != null && !targetMagDbId.trim().isEmpty()) {
				boolean matched = false;
				for (MagDbItem item : magDatabasePanel.getAvailableMagDatabases()) {
					if (targetMagDbId.equals(item.getCatalogueID())) {
						magDatabasePanel.setSelectedMagDatabase(item);
						matched = true;
						if (params.getMagDbVersion() != null && !params.getMagDbVersion().isBlank()) {
							for (MagVersion version : magDatabasePanel.getSelectedMagDatabase().getVersionList()) {
								if (params.getMagDbVersion().equals(version.getCatalogueVersion())) {
									magDatabasePanel.setSelectedMagVersion(version);
									LOGGER.info("Restored MAG database version: {}", version.getCatalogueVersion());
									break;
								}
							}
						}
						LOGGER.info("Restored MAG database selection: {}", targetMagDbId);
						break;
					}
				}
				if (!matched) {
					LOGGER.warn("Saved MAG database '{}' is not among the available databases; "
							+ "leaving the current selection", targetMagDbId);
					consolePanel.appendConsoleText("WARNING: saved MAG database '" + targetMagDbId
							+ "' is not available locally - please re-select it\n");
				}
			}
			magDatabasePanel.setAppendHostDb(params.isAppendHostDb());
			if (params.getHostDB() != null) {
				magDatabasePanel.setHostDatabase(params.getHostDB());
			}
		}

		if (params.getEnzyme() != null) {
			enzymeComboBox.getSelectionModel().select(params.getEnzyme());
		}
		missedCleavagesSpinner.getValueFactory().setValue(params.getMissCleavages());
		if (params.getGenomeSelectorMode() != null) {
			genomeSelectorComboBox.getSelectionModel().select(params.getGenomeSelectorMode());
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
		if (params.getThreadCount() > 0 && threadsComboBox.getItems().contains(params.getThreadCount())) {
			threadsComboBox.getSelectionModel().select(Integer.valueOf(params.getThreadCount()));
		}
		heapComboBox.getSelectionModel().select(params.getMsfraggerMaxHeap());
		if (entrapmentCheckBox != null) {
			entrapmentCheckBox.setSelected(params.isEntrapmentEnabled());
		}

		// Metadata rows: restored so a reopened panel shows the same table it
		// was saved with, rather than silently regenerating it from the IO panel.
		if (metadataPanel != null && params.getMetadataTableData() != null
				&& !params.getMetadataTableData().isEmpty()) {
			List<MetaPilotDdaIgcMaxQuantPanel.MetadataRow> metadataRows = new ArrayList<>();
			for (String[] data : params.getMetadataTableData()) {
				metadataRows.add(MetaPilotDdaIgcMaxQuantPanel.MetadataRow.fromArray(data));
			}
			metadataPanel.setMetadataTableDataWithColumns(metadataRows);
		}

		updateParameterView();
	}

	/** Applies a parameter object loaded from an external file. */
	public void loadParametersFromObject(Object loaded) {
		if (loaded instanceof DDAMagOpenSearchParameters) {
			parametersLoadedFromFile = true;
			applyParametersToUI((DDAMagOpenSearchParameters) loaded);
		} else if (loaded != null) {
			LOGGER.warn("Ignoring parameters of unexpected type: {}", loaded.getClass().getName());
		}
	}

	public void saveCurrentParameters() {
		try {
			DDAMagOpenSearchParameters params = collectParametersFromUI();
			if (params != null) {
				parameterIO.saveLastUsedParameters(params);
			}
		} catch (Exception e) {
			LOGGER.error("Failed to save open-search parameters", e);
		}
	}

	@Override
	protected void loadLastUsedParameters() {
		// Opened from an explicit parameter file: those win, don't overwrite them.
		if (parametersLoadedFromFile) {
			consolePanel.appendConsoleText("Skipping last used parameters - using loaded file parameters\n");
			LOGGER.info("Skipping loadLastUsedParameters(): parameters already loaded from file");
			return;
		}

		try {
			consolePanel.appendConsoleText("Loading last used parameters...\n");
			Object loaded = parameterIO.loadLastUsedParameters();

			if (loaded instanceof DDAMagOpenSearchParameters) {
				DDAMagOpenSearchParameters params = (DDAMagOpenSearchParameters) loaded;
				consolePanel.appendConsoleText("Successfully loaded last used parameters\n");
				applyParametersToUI(params);

				// The IO panel is not covered by applyParametersToUI's control
				// bindings, so restore its contents explicitly, as the other MAG
				// panels do.
				if (params.getExperimentData() != null && !params.getExperimentData().isEmpty()) {
					List<MetaPilotIOPanelFX.ExperimentRow> rows = new ArrayList<>();
					for (String[] data : params.getExperimentData()) {
						rows.add(MetaPilotIOPanelFX.ExperimentRow.fromArray(data));
					}
					ioPanel.setExperimentData(rows);
					consolePanel.appendConsoleText(
							"Applied " + params.getExperimentData().size() + " experiment rows to IO panel\n");
				}
				if (params.getResult() != null && !params.getResult().trim().isEmpty()) {
					ioPanel.setOutputDirectory(params.getResult());
					consolePanel.appendConsoleText("Applied output directory: " + params.getResult() + "\n");
				}

				Platform.runLater(this::updateParameterView);
				LOGGER.info("Applied last used parameters to UI");

			} else {
				// Explicit defaults rather than whatever the controls happen to
				// hold, so a first run starts from a known state.
				consolePanel.appendConsoleText("No previous parameters found, creating defaults\n");
				applyParametersToUI(new DDAMagOpenSearchParameters());
				consolePanel.appendConsoleText("Default parameters loaded into UI\n");
				LOGGER.info("No last used parameters found; applied defaults");
			}

		} catch (Exception e) {
			consolePanel.appendConsoleText("WARNING: Failed to load last used parameters: " + e.getMessage() + "\n");
			LOGGER.error("Failed to load last used open-search parameters", e);
			// Continue with the current UI state rather than failing the panel.
		}
	}

	// ---------------------------------------------------------------- run

	/**
	 * Pre-run validation, mirroring the other MAG panels: checks the shared IO
	 * and MAG database components through their own validators rather than
	 * re-implementing those checks here.
	 */
	private boolean validateInputs() {
		LOGGER.info("Starting input validation for DDA MAG open-search workflow");

		if (ioPanel == null || ioPanel.getExperimentData().isEmpty()) {
			warn("No spectra", "Add spectra files on the Input/Output tab before running.");
			return false;
		}
		String outputDir = ioPanel.getOutputDirectory();
		if (outputDir == null || outputDir.trim().isEmpty()) {
			warn("No result folder", "Choose a result folder on the Input/Output tab.");
			return false;
		}
		if (magDatabasePanel == null || !magDatabasePanel.validateSelection()) {
			String message = magDatabasePanel == null ? "MAG database panel is unavailable."
					: magDatabasePanel.getValidationErrorMessage();
			warn("MAG database", message + "\n\nGenome discovery needs it to build the refined per-sample "
					+ "database that round two searches against.");
			return false;
		}
		LOGGER.info("Input validation passed: {} experiments, output {}", ioPanel.getExperimentData().size(),
				outputDir);
		return true;
	}

	/**
	 * {@code ConsolePanel} flips the UI into its running state — buttons
	 * disabled, status timer started — <em>before</em> invoking this callback.
	 * So every path that does not actually start a task has to put that state
	 * back, or the window sits there looking busy forever. That covers both
	 * thrown exceptions (the task constructor validates raw files and throws)
	 * and ordinary early returns from validation.
	 */
	@Override
	protected void startAnalysis() {
		boolean started = false;
		try {
			started = startAnalysisInternal();
		} catch (Exception e) {
			LOGGER.error("Failed to start open-search analysis", e);
			// Full trace to the Console tab as well as the application log: a
			// failure before the task exists cannot reach the run's logging.log,
			// so the console is the user's only in-app view of the cause.
			consolePanel.appendConsoleText("ERROR: Failed to start analysis: "
					+ e.getClass().getSimpleName() + ": " + e.getMessage() + "\n");
			StringWriter trace = new StringWriter();
			e.printStackTrace(new PrintWriter(trace));
			consolePanel.appendConsoleText(trace.toString() + "\n");
			warn("Could not start analysis", e.getMessage() == null ? e.toString() : e.getMessage());
		} finally {
			if (!started) {
				abortStart();
			}
		}
	}

	/** Returns the UI to idle after a start that never began. */
	private void abortStart() {
		consolePanel.appendConsoleText("Analysis not started.\n");
		consolePanel.setButtonsEnabled(true, false, false);
		consolePanel.stopStatusTimer();
		consolePanel.getOverallProgressBar().progressProperty().unbind();
		consolePanel.getCurrentTaskProgressBar().progressProperty().unbind();
		consolePanel.getOverallProgressBar().setProgress(0);
		consolePanel.getCurrentTaskProgressBar().setProgress(0);
	}

	/** @return true only once a task thread has actually been started */
	private boolean startAnalysisInternal() {
		if (!validateInputs()) {
			return false;
		}

		DDAMagOpenSearchParameters params = collectParametersFromUI();
		if (params == null) {
			consolePanel.appendConsoleText("ERROR: Failed to collect parameters from UI.\n");
			return false;
		}
		// Kept for viewReport(), which resolves the report directory from it.
		this.currentAnalysisParams = params;

		// Every other MAG/IGC/DIA panel writes its analysis-specific parameters
		// JSON into the result folder here (see e.g. MetaPilotDdaMagSagePanel);
		// this panel never did, even though OpenSearchLayout.PARAMETER_FILE /
		// DDAMagOpenSearchParameterIO.getAnalysisParameterFileName() already
		// define "mag_opensearch_analysis_parameters.json" for it - so the file
		// handleTaskFailure() tells users to attach when reporting a bug never
		// actually existed for this workflow. Found and fixed 2026-09-18.
		if (params.getResult() != null && !params.getResult().trim().isEmpty()) {
			boolean analysisParamsSaved = parameterIO.saveAnalysisParameters(params.getResult(), params);
			if (!analysisParamsSaved) {
				consolePanel.appendConsoleText("WARNING: Failed to save analysis parameters to result directory\n");
				LOGGER.warn("Failed to save analysis parameters to {}", params.getResult());
			} else {
				consolePanel.appendConsoleText("\u2713 Analysis parameters saved to result directory\n");
			}
		}

		if (params.isEntrapmentEnabled()) {
			String enzyme = params.getEnzyme();
			if (!"Trypsin".equalsIgnoreCase(enzyme) && !"Trypsin/P".equalsIgnoreCase(enzyme)) {
				warn("Entrapment needs trypsin",
						"Paired-shuffle entrapment currently supports Trypsin or Trypsin/P. "
								+ "Change the enzyme, or turn the evaluation switch off.");
				return false;
			}
			consolePanel.appendConsoleText("Evaluation: paired-shuffle entrapment -> "
					+ new File(params.getResult(),
							uk.ac.quadram.metapilot.core.opensearch.OpenSearchEntrapment.RESULT_FOLDER)
									.getAbsolutePath()
					+ "\n");
		}
		if (params.getMagDbItem() == null) {
			warn("No MAG database", "Select a MAG database on the Input/Output tab.");
			return false;
		}
		if (params.getMetadata() == null || params.getMetadata().getRawFiles() == null
				|| params.getMetadata().getRawFiles().length == 0) {
			// The task constructor would otherwise throw "Raw files are required"
			// from MetaPilotMagTask.createSeparateFolders().
			warn("No spectra files", "No raw files were collected from the Input/Output tab. Add spectra files "
					+ "there before running.");
			return false;
		}
		File msfraggerJar = resourceParams == null ? null : resourceParams.resolveMsfraggerJar();
		if (msfraggerJar == null) {
			warn("MSFragger not found", "Set the FragPipe entry under Resources to its bin/fragpipe launcher. "
					+ "MetaPilot detects the bundled tools/MSFragger-*/MSFragger-*.jar automatically.");
			return false;
		}

		consolePanel.appendConsoleText("Starting DDA MAG open-search workflow ...\n");

		// Progress bars come from the console panel, as in the other MAG panels:
		// the overall bar tracks the task's progressProperty, the second bar
		// shows step activity.
		Label statusLabel = new Label("Ready to start analysis...");
		// currentTask is the base class's field, typed MetaPilotTask - this task
		// now extends MetaPilotMagTask, so it fits and inherits the per-run
		// logging.log, progress wiring and shared MAG post-analysis helpers.
		currentTask = new DDAMagOpenSearchTask(params, resourceParams, consolePanel.getOverallProgressBar(),
				consolePanel.getCurrentTaskProgressBar(), statusLabel);

		currentTask.messageProperty().addListener((obs, oldMsg, newMsg) -> Platform.runLater(() -> {
			if (newMsg != null) {
				consolePanel.appendConsoleText("Current Task: " + newMsg + "\n");
				consolePanel.updateCurrentTaskFromMessage(newMsg);
				statusLabel.setText(newMsg);
			}
		}));

		currentTask.setOnSucceeded(e -> Platform.runLater(() -> {
			try {
				Boolean ok = currentTask.get();
				if (ok != null && ok) {
					consolePanel.appendConsoleText("Analysis completed successfully.\n");
					consolePanel.enableViewReportButton();
					consolePanel.setButtonsEnabled(true, false, true);
				} else {
					consolePanel.appendConsoleText("Analysis finished but reported failure.\n");
					consolePanel.setButtonsEnabled(true, false, false);
				}
			} catch (Exception ex) {
				consolePanel.appendConsoleText("Analysis failed: " + ex.getMessage() + "\n");
				consolePanel.setButtonsEnabled(true, false, false);
			}
			consolePanel.stopStatusTimer();
		}));

		currentTask.setOnFailed(e -> Platform.runLater(() -> {
			consolePanel.appendConsoleText("Analysis failed. Check console output for details.\n");
			consolePanel.setButtonsEnabled(true, false, false);
			consolePanel.stopStatusTimer();
		}));

		consolePanel.setButtonsEnabled(false, true, false);
		consolePanel.startStatusTimer();

		Thread thread = new Thread(currentTask, "DDA-MAG-OpenSearch");
		thread.setDaemon(true);
		thread.start();
		return true;
	}

	@Override
	protected void stopAnalysis() {
		Alert alert = new Alert(AlertType.CONFIRMATION);
		alert.setTitle("Stop Analysis");
		alert.setHeaderText("Are you sure you want to stop the analysis?");
		alert.setContentText("Completed steps are kept, so a rerun resumes from where it stopped rather than "
				+ "starting over.");
		Optional<ButtonType> result = alert.showAndWait();
		if (result.isPresent() && result.get() == ButtonType.OK) {
			consolePanel.appendConsoleText("Stopping open-search analysis ...\n");
			if (currentTask != null && currentTask.isRunning()) {
				currentTask.forceStop();
			}
			// Same teardown as the other MAG panels: release any binding on the
			// bars before resetting, or setProgress throws on a bound property.
			consolePanel.getOverallProgressBar().progressProperty().unbind();
			consolePanel.getCurrentTaskProgressBar().progressProperty().unbind();
			consolePanel.getOverallProgressBar().setProgress(0);
			consolePanel.getCurrentTaskProgressBar().setProgress(0);
			consolePanel.setButtonsEnabled(true, false, false);
			consolePanel.stopStatusTimer();
		}
	}

	@Override
	protected void updateParameterView() {
		if (consolePanel != null) {
			consolePanel.setParameterText(getParameterDisplay());
		}
	}

	private String heapValue() {
		if (heapComboBox == null) {
			return "32g";
		}
		if (heapComboBox.isEditable() && heapComboBox.getEditor() != null) {
			String typed = heapComboBox.getEditor().getText();
			if (typed != null && !typed.isBlank()) {
				return typed.trim();
			}
		}
		return heapComboBox.getValue() == null || heapComboBox.getValue().isBlank() ? "32g" : heapComboBox.getValue();
	}

	private String getParameterDisplay() {
		StringBuilder sb = new StringBuilder();
		sb.append("=== DDA MAG Open Search ===\n");
		if (ioPanel != null) {
			sb.append("Result folder: ").append(ioPanel.getOutputDirectory()).append("\n");
			sb.append("Spectra files: ").append(ioPanel.getExperimentData().size()).append("\n");
		}
		if (magDatabasePanel != null && magDatabasePanel.getSelectedMagDatabase() != null) {
			sb.append("MAG database: ").append(magDatabasePanel.getSelectedMagDatabase().getCatalogueID())
					.append("\n");
		}
		if (enzymeComboBox != null) {
			sb.append("Search engine: MSFragger (closed marker, then joint closed+open)\n");
			sb.append("Enzyme: ").append(enzymeComboBox.getValue()).append(", missed cleavages ")
					.append(missedCleavagesSpinner.getValue()).append("\n");
			sb.append("Marker set: ").append(genomeSelectorComboBox.getValue()).append("\n");
			if (genomeStopControls != null) {
				sb.append("When to stop adding genomes: ").append(genomeStopControls.summary()).append("\n");
			}
			sb.append("Search method: tested defaults (closed marker search, then joint closed+open identification)\n");
			sb.append("Threads: ").append(threadsComboBox.getValue()).append(", MSFragger memory ")
					.append(heapValue()).append("\n");
			if (entrapmentCheckBox != null) {
				sb.append("Paired-shuffle entrapment: ").append(entrapmentCheckBox.isSelected()).append("\n");
			}
		}
		return sb.toString();
	}

	@Override
	protected Class<DDAMagOpenSearchParameters> getParameterClass() {
		return DDAMagOpenSearchParameters.class;
	}

	@Override
	protected String getWorkflowName() {
		return "DDA MAG Open Search";
	}

	private void warn(String header, String message) {
		Alert alert = new Alert(AlertType.WARNING);
		alert.setTitle("Open Search");
		alert.setHeaderText(header);
		alert.setContentText(message);
		alert.showAndWait();
	}
}
