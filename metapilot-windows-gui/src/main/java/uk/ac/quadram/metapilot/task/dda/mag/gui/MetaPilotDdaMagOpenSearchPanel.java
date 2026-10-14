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
import javafx.scene.control.ListCell;
import javafx.scene.control.Spinner;
import javafx.scene.control.SpinnerValueFactory;
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
import uk.ac.quadram.metapilot.task.common.gui.ModificationsPanel;
import uk.ac.quadram.metapilot.task.dda.igc.gui.MetaPilotDdaIgcMaxQuantPanel;
import uk.ac.quadram.metapilot.task.dda.mag.par.DDAMagOpenSearchParameterIO;
import uk.ac.quadram.metapilot.task.dda.mag.par.DDAMagOpenSearchParameters;
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
 * Note the Parameters tab has no Closed/Open selector: both search rounds are
 * open searches by construction, with round one consuming only its unmodified
 * PSMs for genome discovery. That is the design, not a setting.
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

	// Workflow-specific controls (Parameters tab)
	private ComboBox<String> engineComboBox;
	private Label engineHintLabel;
	private ComboBox<String> enzymeComboBox;
	private Spinner<Integer> missedCleavagesSpinner;
	private ComboBox<String> genomeSelectorComboBox;
	private ComboBox<String> genomeSelectionStrategyComboBox;
	private ComboBox<Integer> threadsComboBox;
	private Spinner<Double> precursorLowerSpinner;
	private Spinner<Double> precursorUpperSpinner;
	private Spinner<Integer> topNSpinner;
	private Spinner<Integer> round1TopNSpinner;
	private Spinner<Integer> round2TopNSpinner;
	private ComboBox<String> round1SearchModeComboBox;
	private CheckBox closedStageCheckBox;
	private ComboBox<String> heapComboBox;

	// Validation (advanced)
	private Spinner<Double> localFdrSpinner;
	private Spinner<Integer> classifierFoldsSpinner;
	private Spinner<Double> zeroPeakSigmasSpinner;
	private Spinner<Double> peakMinR2Spinner;
	private Spinner<Double> membershipFloorSpinner;
	private Spinner<Integer> peakMinAreaSpinner;
	private Spinner<Integer> maxPeaksToFitSpinner;
	private Spinner<Integer> maxFitIterationsSpinner;
	private Spinner<Integer> peakMinFitBinsSpinner;
	private Spinner<Double> peakMinSigmaBinsSpinner;
	private Spinner<Double> peakMinSnrSpinner;
	private Spinner<Double> peakMinSeparationBinsSpinner;
	private Spinner<Double> peakMinSeparationSigmasSpinner;

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
		metadataPanel = new MetadataConfigurationPanel();
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
	 * False until a {@code mag-opensearch} subcommand exists in
	 * {@code MetaPilotCliApp}. The workflow logic already lives in
	 * {@code MagOpenSearchHeadlessRunner} so it is CLI-ready, but exporting HPC
	 * parameters for a command that cannot yet be invoked would hand the user a
	 * file nothing consumes. The base class shows an explanatory message in
	 * place of the button.
	 */
	@Override
	protected boolean isHpcExportSupported() {
		return false;
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
		tabContent.setMinHeight(760);

		GridPane mainGrid = new GridPane();
		mainGrid.setHgap(15);
		mainGrid.setVgap(15);
		mainGrid.setPadding(new Insets(10));
		mainGrid.setMinWidth(940);
		mainGrid.setMinHeight(740);
		VBox.setVgrow(mainGrid, Priority.ALWAYS);

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
		leftColumn.getChildren().addAll(createSearchSection(), createGenomeSection(), createValidationSection(),
				modificationsPanel);
		GridPane.setVgrow(leftColumn, Priority.ALWAYS);

		VBox rightColumn = new VBox(15);
		rightColumn.setMinWidth(430);
		rightColumn.getChildren().addAll(metadataPanel);
		GridPane.setVgrow(rightColumn, Priority.ALWAYS);

		mainGrid.add(leftColumn, 0, 0);
		mainGrid.add(rightColumn, 1, 0);

		tabContent.getChildren().add(mainGrid);
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

		grid.add(new Label("Search engine(s):"), 0, row);
		engineComboBox = new ComboBox<>();
		engineComboBox.getItems().addAll(DDAMagOpenSearchParameters.ENGINES_MSFRAGGER,
				DDAMagOpenSearchParameters.ENGINES_MSFRAGGER_AND_SAGE);
		engineComboBox.getSelectionModel().select(DDAMagOpenSearchParameters.ENGINES_MSFRAGGER);
		engineComboBox.setMaxWidth(Double.MAX_VALUE);
		engineComboBox.setTooltip(new Tooltip("Whether a second engine improves coverage enough to justify its "
				+ "cost is an open question against a refined genome-specific database. The option exists so it "
				+ "can be measured rather than assumed."));

		// Unimplemented modes are shown greyed and cannot be picked, rather than
		// being selectable and rejected after the user has configured a whole run.
		// The item stays visible because its absence would hide the roadmap; it is
		// simply inert.
		engineComboBox.setCellFactory(list -> new ListCell<String>() {
			@Override
			protected void updateItem(String item, boolean empty) {
				super.updateItem(item, empty);
				setText(empty || item == null ? null : item);
				boolean unavailable = item != null && !DDAMagOpenSearchParameters.isEngineModeAvailable(item);
				setDisable(unavailable);
				setStyle(unavailable ? "-fx-text-fill: #9a9a9a;" : "");
			}
		});
		// Defence in depth: a disabled cell cannot be clicked, but the value can
		// still be set programmatically (loading a saved parameter file), so revert
		// rather than trust the cell factory alone.
		engineComboBox.valueProperty().addListener((obs, old, value) -> {
			if (value != null && !DDAMagOpenSearchParameters.isEngineModeAvailable(value)) {
				Platform.runLater(() -> engineComboBox.getSelectionModel()
						.select(DDAMagOpenSearchParameters.ENGINES_MSFRAGGER));
				return;
			}
			updateEngineHint(value);
		});
		grid.add(engineComboBox, 1, row++);

		engineHintLabel = new Label();
		engineHintLabel.setWrapText(true);
		engineHintLabel.setMaxWidth(Double.MAX_VALUE);
		// GridPane sizes a wrapping, column-spanning Label's row using its
		// pre-wrap preferred height unless min height is pinned to the
		// (post-wrap) preferred height - without this the row stays one line
		// tall and the rest of the text is clipped with an ellipsis.
		engineHintLabel.setMinHeight(Region.USE_PREF_SIZE);
		grid.add(engineHintLabel, 0, row++, 2, 1);
		updateEngineHint(engineComboBox.getValue());

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

		grid.add(new Label("Precursor window lower (Da):"), 0, row);
		precursorLowerSpinner = doubleSpinner(-1000, 0, -150, 10);
		grid.add(precursorLowerSpinner, 1, row++);

		grid.add(new Label("Precursor window upper (Da):"), 0, row);
		precursorUpperSpinner = doubleSpinner(0, 2000, 500, 10);
		grid.add(precursorUpperSpinner, 1, row++);

		Label windowHint = new Label("The wide window is what diverts modified spectra out of the zero peak, "
				+ "which is what makes the unmodified set usable for genome discovery.");
		windowHint.setWrapText(true);
		windowHint.setMaxWidth(Double.MAX_VALUE);
		windowHint.setMinHeight(Region.USE_PREF_SIZE);
		windowHint.setStyle("-fx-font-size: 11px; -fx-text-fill: #666;");
		grid.add(windowHint, 0, row++, 2, 1);

		grid.add(new Label("Candidates per spectrum:"), 0, row);
		topNSpinner = new Spinner<>(1, 20, 5);
		topNSpinner.setMaxWidth(Double.MAX_VALUE);
		topNSpinner.setTooltip(new Tooltip("MSFragger reports only the top hit by default; reporting more lets "
				+ "the classifier choose between candidates rather than the engine."));
		grid.add(topNSpinner, 1, row++);

		grid.add(new Label("Round 2 candidates per spectrum:"), 0, row);
		round2TopNSpinner = new Spinner<>(1, 20, 5);
		round2TopNSpinner.setMaxWidth(Double.MAX_VALUE);
		round2TopNSpinner.setTooltip(new Tooltip("Per-sample refined-database search (round two). Measured: "
				+ "raising this from 1 to 5 gains +15.3% PSMs and +5.4% peptides at the same FDR, since the extra "
				+ "ranks feed the peak-discovery votes and the classifier's negatives. Default 5."));
		grid.add(round2TopNSpinner, 1, row++);

		closedStageCheckBox = new CheckBox("Closed search first, open search on the remaining spectra");
		closedStageCheckBox.setSelected(true);
		closedStageCheckBox.setTooltip(new Tooltip("Each sample's refined database is also searched closed "
				+ "(+-20 ppm, isotope errors 0/1/2, M oxidation, protein N-terminal acetyl). Closed PSMs "
				+ "accepted at 1% FDR are identifications; the open search then identifies only the spectra "
				+ "they leave, as unmodified or allow-list modified peptides. Unexplained mass shifts are "
				+ "reported in unexplained_psms.tsv but are not identifications. Off: open search only."));
		grid.add(closedStageCheckBox, 0, row++, 2, 1);

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

		TitledPane pane = new TitledPane("Search", grid);
		pane.setCollapsible(false);
		return pane;
	}

	private TitledPane createGenomeSection() {
		GridPane grid = paramGrid();
		int row = 0;

		grid.add(new Label("Genome selector:"), 0, row);
		genomeSelectorComboBox = new ComboBox<>();
		genomeSelectorComboBox.getItems().addAll(MagDbItem.MARKER2, MagDbItem.MARKER4);
		genomeSelectorComboBox.getSelectionModel().select(MagDbItem.MARKER2);
		genomeSelectorComboBox.setMaxWidth(Double.MAX_VALUE);
		grid.add(genomeSelectorComboBox, 1, row++);

		grid.add(new Label("Selection strategy:"), 0, row);
		genomeSelectionStrategyComboBox = new ComboBox<>();
		genomeSelectionStrategyComboBox.getItems().addAll("Focus", "Comprehensive");
		genomeSelectionStrategyComboBox.getSelectionModel().select("Focus");
		genomeSelectionStrategyComboBox.setMaxWidth(Double.MAX_VALUE);
		grid.add(genomeSelectionStrategyComboBox, 1, row++);

		grid.add(new Label("Round 1 candidates per spectrum:"), 0, row);
		round1TopNSpinner = new Spinner<>(1, 20, 1);
		round1TopNSpinner.setMaxWidth(Double.MAX_VALUE);
		round1TopNSpinner.setTooltip(new Tooltip("Round one has no classifier, so nothing can promote a "
				+ "lower-ranked candidate: rank 1 is argmax(hyperscore) for 99.999% of spectra, and even a "
				+ "perfect oracle rescorer could add at most 4.95% true-genome hits. Stays at 1 by default."));
		grid.add(round1TopNSpinner, 1, row++);

		grid.add(new Label("Round 1 marker search:"), 0, row);
		round1SearchModeComboBox = new ComboBox<>();
		round1SearchModeComboBox.getItems().addAll(DDAMagOpenSearchParameters.ROUND1_SEARCH_CLOSED,
				DDAMagOpenSearchParameters.ROUND1_SEARCH_OPEN);
		round1SearchModeComboBox.getSelectionModel().select(DDAMagOpenSearchParameters.ROUND1_SEARCH_CLOSED);
		round1SearchModeComboBox.setMaxWidth(Double.MAX_VALUE);
		round1SearchModeComboBox.setTooltip(new Tooltip("closed: +-20 ppm with isotope errors 0/1/2 and "
				+ "variable M oxidation / protein N-terminal acetyl; about five times faster and selects the "
				+ "same genomes as open. open: the round-two precursor window, with Crystal-C."));
		grid.add(round1SearchModeComboBox, 1, row++);

		grid.add(new Label("Gaussian membership (sigmas):"), 0, row);
		zeroPeakSigmasSpinner = doubleSpinner(1, 6, 3, 0.5);
		zeroPeakSigmasSpinner.setTooltip(new Tooltip(
				"Membership width used only by the legacy Gaussian discovery option. "
						+ "The default KDE strategy uses empirical density-valley boundaries."));
		grid.add(zeroPeakSigmasSpinner, 1, row++);

		Label stageHint = new Label("Only round one's unmodified zero family is used to map genomes. PSMs "
				+ "near +/-1.00335 Da multiples are normalized onto the zero mode rather than treated as "
				+ "separate modifications.");
		stageHint.setWrapText(true);
		stageHint.setMaxWidth(Double.MAX_VALUE);
		stageHint.setMinHeight(Region.USE_PREF_SIZE);
		stageHint.setStyle("-fx-font-size: 11px; -fx-text-fill: #666;");
		grid.add(stageHint, 0, row++, 2, 1);

		TitledPane pane = new TitledPane("Genome discovery (round one)", grid);
		pane.setCollapsible(false);
		return pane;
	}

	private TitledPane createValidationSection() {
		GridPane grid = paramGrid();
		int row = 0;

		grid.add(new Label("Local FDR:"), 0, row);
		localFdrSpinner = doubleSpinner(0.001, 0.2, 0.01, 0.005);
		grid.add(localFdrSpinner, 1, row);

		grid.add(new Label("Classifier folds:"), 2, row);
		classifierFoldsSpinner = new Spinner<>(2, 10, 3);
		classifierFoldsSpinner.setMaxWidth(Double.MAX_VALUE);
		grid.add(classifierFoldsSpinner, 3, row++);

		grid.add(new Label("Peak min R-squared:"), 0, row);
		peakMinR2Spinner = doubleSpinner(0.5, 0.99, 0.9, 0.05);
		grid.add(peakMinR2Spinner, 1, row);

		grid.add(new Label("Membership floor (Da):"), 2, row);
		membershipFloorSpinner = doubleSpinner(0.0, 0.1, 0.010, 0.005);
		membershipFloorSpinner.setTooltip(new Tooltip(
				"Absolute lower bound on peak membership half-width.\n\n"
						+ "A PSM joins a peak when it is within max(3 sigma, this). Delta-mass peaks have a sharp "
						+ "core (~3 mDa FWHM) and heavy tails, so a pure 3-sigma rule cuts at ~3 mDa and discards "
						+ "real members measured with slightly worse mass accuracy.\n\n"
						+ "Set this to your instrument mass accuracy, not the fitted core width. 0 disables it."));
		grid.add(membershipFloorSpinner, 3, row++);

		grid.add(new Label("Peak min PSMs:"), 0, row);
		peakMinAreaSpinner = new Spinner<>(1, 1000, 15);
		peakMinAreaSpinner.setMaxWidth(Double.MAX_VALUE);
		peakMinAreaSpinner.setTooltip(new Tooltip(
				"Minimum PSMs in a delta-mass peak before a Gaussian is fitted to it.\n\n"
						+ "With 1:1 target:decoy a peak of n PSMs has a best-case binomial p of 0.5^n - the p-value "
						+ "it gets if EVERY PSM in it is a target. At n=3 that is 0.125, at n=7 it is 0.0078, so "
						+ "smaller peaks cannot reach significance however clean they are.\n\n"
						+ "Fitting them is not just slow: each one raises the multiple-testing burden applied to "
						+ "the real peaks. Lower this only if you expect very low-abundance modifications and "
						+ "accept the loss of power elsewhere."));
		grid.add(peakMinAreaSpinner, 1, row++);

		grid.add(new Label("Max peaks to fit:"), 0, row);
		maxPeaksToFitSpinner = new Spinner<>(0, 100000, 500);
		maxPeaksToFitSpinner.setMaxWidth(Double.MAX_VALUE);
		maxPeaksToFitSpinner.setTooltip(new Tooltip(
				"Hard cap on Gaussian fits per run; the largest peaks are kept. 0 means no cap.\n\n"
						+ "Bounds the worst case so this step cannot run for an unknown length of time on an "
						+ "unusual histogram."));
		grid.add(maxPeaksToFitSpinner, 1, row);

		grid.add(new Label("Max fit iterations:"), 2, row);
		maxFitIterationsSpinner = new Spinner<>(20, 5000, 200);
		maxFitIterationsSpinner.setMaxWidth(Double.MAX_VALUE);
		maxFitIterationsSpinner.setTooltip(new Tooltip(
				"Iteration ceiling for a single Gaussian fit.\n\n"
						+ "A genuinely peak-shaped segment converges quickly; one that is not runs to the limit "
						+ "and is discarded anyway."));
		grid.add(maxFitIterationsSpinner, 3, row++);

		grid.add(new Label("Peak min fit bins:"), 0, row);
		peakMinFitBinsSpinner = new Spinner<>(4, 101, 7);
		peakMinFitBinsSpinner.setMaxWidth(Double.MAX_VALUE);
		peakMinFitBinsSpinner.setTooltip(new Tooltip(
				"Minimum raw histogram bins used by a peak fit. Seven prevents a four-parameter fit from being "
						+ "accepted with no meaningful residual degrees of freedom."));
		grid.add(peakMinFitBinsSpinner, 1, row);

		grid.add(new Label("Peak min sigma (bins):"), 2, row);
		peakMinSigmaBinsSpinner = doubleSpinner(0.0, 10.0, 0.5, 0.1);
		peakMinSigmaBinsSpinner.setTooltip(new Tooltip(
				"Reject fitted widths below the histogram''s resolving power. A 0.5-bin default means a 0.5 mDa "
						+ "minimum sigma with the standard 1 mDa bins."));
		grid.add(peakMinSigmaBinsSpinner, 3, row++);

		grid.add(new Label("Peak min SNR:"), 0, row);
		peakMinSnrSpinner = doubleSpinner(0.0, 1000.0, 5.0, 0.5);
		peakMinSnrSpinner.setTooltip(new Tooltip(
				"Minimum fitted height divided by sqrt(local background + 1). Rejects fluctuations that fit well "
						+ "geometrically but do not rise meaningfully above the continuum."));
		grid.add(peakMinSnrSpinner, 1, row);

		grid.add(new Label("Min separation (bins):"), 2, row);
		peakMinSeparationBinsSpinner = doubleSpinner(0.0, 20.0, 2.0, 0.5);
		grid.add(peakMinSeparationBinsSpinner, 3, row++);

		grid.add(new Label("Min separation (sigma):"), 0, row);
		peakMinSeparationSigmasSpinner = doubleSpinner(0.0, 20.0, 3.0, 0.5);
		peakMinSeparationSigmasSpinner.setTooltip(new Tooltip(
				"Alternative fits closer than max(the bin floor, this many combined fitted sigmas) are unresolved; "
						+ "only the stronger fit is retained."));
		grid.add(peakMinSeparationSigmasSpinner, 1, row++);

		TitledPane pane = new TitledPane("Validation (round two)", grid);
		pane.setCollapsible(true);
		pane.setExpanded(false);
		return pane;
	}

	private void updateEngineHint(String value) {
		if (engineHintLabel == null) {
			return;
		}
		if (DDAMagOpenSearchParameters.ENGINES_MSFRAGGER_AND_SAGE.equals(value)) {
			engineHintLabel.setText("Sage joins round two only, against the same refined database, so genome "
					+ "selection stays identical and any difference is attributable to the second engine. Each "
					+ "engine is scored by its own classifier - Sage cannot supply the delta-mass localization "
					+ "features - then calibrated to a common scale before pooling. Results go to "
					+ "opensearch_MAG_search/msfragger_sage/, leaving a MSFragger-only result in "
					+ "msfragger_only/ intact; both reuse the searches in search_results/.");
			engineHintLabel.setStyle("-fx-font-size: 11px; -fx-text-fill: #666;");
		} else {
			engineHintLabel.setText("MSFragger invoked directly - no FragPipe, no Philosopher.");
			engineHintLabel.setStyle("-fx-font-size: 11px; -fx-text-fill: #666;");
		}
	}

	private static GridPane paramGrid() {
		GridPane grid = new GridPane();
		grid.setHgap(10);
		grid.setVgap(8);
		grid.setPadding(new Insets(10));
		return grid;
	}

	private static Spinner<Double> doubleSpinner(double min, double max, double initial, double step) {
		Spinner<Double> spinner = new Spinner<>();
		spinner.setValueFactory(new SpinnerValueFactory.DoubleSpinnerValueFactory(min, max, initial, step));
		spinner.setEditable(true);
		spinner.setMaxWidth(Double.MAX_VALUE);
		return spinner;
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
		}

		params.setSearchEngines(engineComboBox.getValue());
		params.setEnzyme(enzymeComboBox.getValue());
		params.setMissCleavages(missedCleavagesSpinner.getValue());
		params.setGenomeSelectorMode(genomeSelectorComboBox.getValue());
		params.setGenomeSelectionStrategy(genomeSelectionStrategyComboBox.getValue());
		params.setThreadCount(threadsComboBox.getValue() == null ? 1 : threadsComboBox.getValue());
		params.setPrecursorMassLower(precursorLowerSpinner.getValue());
		params.setPrecursorMassUpper(precursorUpperSpinner.getValue());
		params.setOutputReportTopN(topNSpinner.getValue());
		params.setRound1OutputReportTopN(round1TopNSpinner.getValue());
		params.setRound2OutputReportTopN(round2TopNSpinner.getValue());
		params.setRound1SearchMode(round1SearchModeComboBox.getValue());
		params.setClosedStageEnabled(closedStageCheckBox.isSelected());
		params.setMsfraggerMaxHeap(heapComboBox.getValue());

		params.setLocalFdr(localFdrSpinner.getValue());
		params.setClassifierFolds(classifierFoldsSpinner.getValue());
		params.setZeroPeakSigmas(zeroPeakSigmasSpinner.getValue());
		params.setPeakMinR2(peakMinR2Spinner.getValue());
		params.setMembershipFloorDa(membershipFloorSpinner.getValue());
		params.setPeakMinArea(peakMinAreaSpinner.getValue());
		params.setMaxPeaksToFit(maxPeaksToFitSpinner.getValue());
		params.setMaxFitIterations(maxFitIterationsSpinner.getValue());
		params.setPeakMinFitBins(peakMinFitBinsSpinner.getValue());
		params.setPeakMinSigmaBins(peakMinSigmaBinsSpinner.getValue());
		params.setPeakMinSnr(peakMinSnrSpinner.getValue());
		params.setPeakMinSeparationBins(peakMinSeparationBinsSpinner.getValue());
		params.setPeakMinSeparationSigmas(peakMinSeparationSigmasSpinner.getValue());

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

		// Saved files may name a mode that was selectable when written, or was
		// hand-edited. Coerce instead of restoring an unrunnable selection.
		String loadedEngine = DDAMagOpenSearchParameters.coerceEngineMode(params.getSearchEngines());
		if (!loadedEngine.equals(params.getSearchEngines())) {
			LOGGER.warn("Saved engine mode '{}' is not implemented; falling back to {}",
					params.getSearchEngines(), loadedEngine);
			params.setSearchEngines(loadedEngine);
		}
		engineComboBox.getSelectionModel().select(loadedEngine);
		updateEngineHint(loadedEngine);
		if (params.getEnzyme() != null) {
			enzymeComboBox.getSelectionModel().select(params.getEnzyme());
		}
		missedCleavagesSpinner.getValueFactory().setValue(params.getMissCleavages());
		if (params.getGenomeSelectorMode() != null) {
			genomeSelectorComboBox.getSelectionModel().select(params.getGenomeSelectorMode());
		}
		if (params.getGenomeSelectionStrategy() != null) {
			genomeSelectionStrategyComboBox.getSelectionModel().select(params.getGenomeSelectionStrategy());
		}
		if (params.getThreadCount() > 0 && threadsComboBox.getItems().contains(params.getThreadCount())) {
			threadsComboBox.getSelectionModel().select(Integer.valueOf(params.getThreadCount()));
		}
		precursorLowerSpinner.getValueFactory().setValue(params.getPrecursorMassLower());
		precursorUpperSpinner.getValueFactory().setValue(params.getPrecursorMassUpper());
		topNSpinner.getValueFactory().setValue(params.getOutputReportTopN());
		round1TopNSpinner.getValueFactory().setValue(params.getRound1OutputReportTopN());
		round2TopNSpinner.getValueFactory().setValue(params.getRound2OutputReportTopN());
		round1SearchModeComboBox.getSelectionModel().select(params.getRound1SearchMode());
		closedStageCheckBox.setSelected(params.isClosedStageEnabled());
		heapComboBox.getSelectionModel().select(params.getMsfraggerMaxHeap());

		localFdrSpinner.getValueFactory().setValue(params.getLocalFdr());
		classifierFoldsSpinner.getValueFactory().setValue(params.getClassifierFolds());
		zeroPeakSigmasSpinner.getValueFactory().setValue(params.getZeroPeakSigmas());
		peakMinR2Spinner.getValueFactory().setValue(params.getPeakMinR2());
		membershipFloorSpinner.getValueFactory().setValue(params.getMembershipFloorDa());
		peakMinAreaSpinner.getValueFactory().setValue(params.getPeakMinArea());
		maxPeaksToFitSpinner.getValueFactory().setValue(params.getMaxPeaksToFit());
		maxFitIterationsSpinner.getValueFactory().setValue(params.getMaxFitIterations());
		peakMinFitBinsSpinner.getValueFactory().setValue(params.getPeakMinFitBins());
		peakMinSigmaBinsSpinner.getValueFactory().setValue(params.getPeakMinSigmaBins());
		peakMinSnrSpinner.getValueFactory().setValue(params.getPeakMinSnr());
		peakMinSeparationBinsSpinner.getValueFactory().setValue(params.getPeakMinSeparationBins());
		peakMinSeparationSigmasSpinner.getValueFactory().setValue(params.getPeakMinSeparationSigmas());

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
		// Should be unreachable now the combo refuses the selection, but kept as
		// the last line of defence for any path that sets parameters directly.
		String engineProblem = DDAMagOpenSearchParameters.engineModeUnavailableReason(params.getSearchEngines());
		if (engineProblem != null) {
			warn("Engine mode not implemented", engineProblem);
			return false;
		}
		if (params.usesSage()) {
			String sagePath = resourceParams == null ? null
					: resourceParams.getToolPath(MetaResourceParameters.ToolType.SAGE);
			if (sagePath == null || sagePath.isBlank() || !new File(sagePath).exists()) {
				warn("Sage not found", "This mode runs Sage as a second engine in round two, but no Sage "
						+ "executable is configured. Set Resources > Sage, or choose \""
						+ DDAMagOpenSearchParameters.ENGINES_MSFRAGGER + "\".");
				return false;
			}
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
		if (engineComboBox != null) {
			sb.append("Engine(s): ").append(engineComboBox.getValue()).append("\n");
			sb.append("Enzyme: ").append(enzymeComboBox.getValue()).append(", missed cleavages ")
					.append(missedCleavagesSpinner.getValue()).append("\n");
			sb.append("Precursor window: ").append(precursorLowerSpinner.getValue()).append(" to ")
					.append(precursorUpperSpinner.getValue()).append(" Da\n");
			sb.append("Candidates per spectrum: ").append(topNSpinner.getValue()).append("\n");
			sb.append("Round 1 candidates per spectrum: ").append(round1TopNSpinner.getValue()).append("\n");
			sb.append("Round 2 candidates per spectrum: ").append(round2TopNSpinner.getValue()).append("\n");
			sb.append("Genome selector: ").append(genomeSelectorComboBox.getValue()).append(" (")
					.append(genomeSelectionStrategyComboBox.getValue()).append(")\n");
			sb.append("Round 1 marker search: ").append(round1SearchModeComboBox.getValue()).append("\n");
			sb.append("Round-one unmodified selection: constrained isotope-folded KDE\n");
			if (closedStageCheckBox.isSelected()) {
				sb.append("Identification: closed search first (1% FDR), open search on the remaining spectra; "
						+ "unmodified and allow-list shifts identified, other masses reported only\n");
			} else {
				sb.append("Identification: open search only\n");
			}
			sb.append("Peak screen: zero family accepted; allow-list shifts accepted when target-enriched; "
					+ "other masses " + (closedStageCheckBox.isSelected() ? "not identified" : "identification-only")
					+ "\n");
			sb.append("Local FDR: ").append(localFdrSpinner.getValue()).append("\n");
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
