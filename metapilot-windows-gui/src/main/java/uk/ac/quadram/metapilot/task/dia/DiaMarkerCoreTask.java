package uk.ac.quadram.metapilot.task.dia;

import static uk.ac.quadram.metapilot.core.dia.DiaMagScoringService.IDX_P_EMP;
import static uk.ac.quadram.metapilot.core.dia.DiaMagScoringService.IDX_Q_VALUE;
import static uk.ac.quadram.metapilot.core.dia.DiaMagScoringService.IDX_S;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.math.BigInteger;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutionException;

import org.apache.commons.lang3.math.NumberUtils;
import org.apache.logging.log4j.Logger;

import javafx.application.Platform;
import javafx.scene.control.ProgressBar;

import uk.ac.quadram.metapilot.core.dia.DiaHostLibrarySupport;
import uk.ac.quadram.metapilot.core.dia.DiaMagScoringService;
import uk.ac.quadram.metapilot.core.dia.Marker4LibrarySupport;
import uk.ac.quadram.metapilot.core.mag.MagGenomeSelectionService;
import uk.ac.quadram.metapilot.core.mag.Marker4ResourcePreparer;
import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNParquetReader;
import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNPrecursor;
import uk.ac.quadram.metapilot.dbSearch.diann.DiannPrMatrixReader;
import uk.ac.quadram.metapilot.core.runner.DiaMarkerHeadlessRunner;
import uk.ac.quadram.metapilot.core.runner.DiaNNNoQuantLibrarySearchService;
import uk.ac.quadram.metapilot.core.runner.MagFunctionalAnnotationService;
import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;
import uk.ac.quadram.metapilot.quant.flashLFQ.FlashLfqQuanPeptide;
import uk.ac.quadram.metapilot.core.tools.FormatTool;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters.ToolType;
import uk.ac.quadram.metapilot.task.dia.par.MetaParameterDia;
import uk.ac.quadram.metapilot.task.dia.par.MetaParameterDiaMarker;
import uk.ac.quadram.metapilot.task.io.ReportManifestGenerator;
import uk.ac.quadram.metapilot.task.io.pep.MetaPeptide;
import uk.ac.quadram.metapilot.task.io.pro.MetaProtein;
import uk.ac.quadram.metapilot.task.mag.MagDbItemResolver;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;
import uk.ac.quadram.metapilot.task.mag.MarkerGenomeSelector;
import uk.ac.quadram.metapilot.task.mag.MetaProteinAnnoMag;
import uk.ac.quadram.metapilot.task.MetaPythonReportTask;
import uk.ac.quadram.metapilot.task.MetaPythonReportCopyTask;

/**
 * DiaMarkerCoreTask - A specialized DIA task for marker-based DIA analysis
 * 
 * This class extends MetaDiaTask to provide specific functionality for DIA
 * Marker workflow (Focus or Comprehensive genome selection), which uses
 * pre-generated focus libraries and MarkerGenomeSelector with user-selected mode.
 */
public class DiaMarkerCoreTask extends MetaDiaTask {

	// Logger is inherited from MetaDiaTask -> MetaPilotMagTask -> MetaPilotTask as
	// an instance field
	private static final String taskName = "DIA Marker analysis";

	// Combined file folder (set during iterativeSearch and used in iden)
	protected File combinedFile;

	/**
	 * Constructor for DiaMarkerCoreTask
	 * 
	 * @param markerParams The DIA Marker (Focus/Comprehensive) specific parameters
	 * @param mrp         Resource parameters
	 * @param bar1        Main progress bar
	 * @param bar2        Secondary progress bar
	 */
	public DiaMarkerCoreTask(MetaParameterDiaMarker markerParams, MetaResourceParameters mrp, ProgressBar bar1,
			ProgressBar bar2) {
		this(withMode(markerParams, DiaRefinementMode.CORE), mrp, bar1, bar2, true);
	}

	protected DiaMarkerCoreTask(MetaParameterDiaMarker markerParams, MetaResourceParameters mrp, ProgressBar bar1,
			ProgressBar bar2, DiaRefinementMode refinementMode) {
		this(withMode(markerParams, refinementMode), mrp, bar1, bar2, true);
	}

	private DiaMarkerCoreTask(MetaParameterDiaMarker markerParams, MetaResourceParameters mrp, ProgressBar bar1,
			ProgressBar bar2, boolean modeAlreadyApplied) {
		super(markerParams, mrp, bar1, bar2);

		// Initialize result folder for marker workflows if not already set:
		// marker_focus / marker_comprehensive, based on genome selection strategy.
		initResultFolderIfNeeded();

		// Apply modifications from parameters to DiannParameter
		if (diannPar != null) {
			diannPar.setM_excision(markerParams.isMExcision());
			diannPar.setCarbami_C(markerParams.isCarbamiC());
			diannPar.setOxidation_M(markerParams.isOxidationM());
			diannPar.setAcety_proN(markerParams.isAcetyProN());
			diannPar.setPhospho_STY(markerParams.isPhosphoSTY());
			diannPar.setK_GG(markerParams.isKGG());
			logger.info(
					"Applied modifications from parameters: M_excision={}, Carbami_C={}, Oxidation_M={}, Acety_proN={}, Phospho_STY={}, K_GG={}",
					markerParams.isMExcision(), markerParams.isCarbamiC(), markerParams.isOxidationM(),
					markerParams.isAcetyProN(), markerParams.isPhosphoSTY(), markerParams.isKGG());
		}

		logger.info("Initializing DiaMarkerCoreTask...");
		System.out.println(new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss").format(new Date()) + "\t" + getTaskName()
				+ ": initializing task");

		logger.info("DIA Focus workflow initialized with library search mode");
	}

	private static MetaParameterDiaMarker withMode(MetaParameterDiaMarker markerParams, DiaRefinementMode refinementMode) {
		if (markerParams != null && refinementMode != null) {
			markerParams.setMarkerRefinementMode(refinementMode);
		}
		return markerParams;
	}

	/**
	 * Apply modification settings to DiannParameter Called from the UI panel after
	 * task creation
	 */
	public void applyModifications(boolean mExcision, boolean carbamiC, boolean oxidationM, boolean acetyProN,
			boolean phosphoSTY, boolean kGG) {
		if (diannPar != null) {
			diannPar.setM_excision(mExcision);
			diannPar.setCarbami_C(carbamiC);
			diannPar.setOxidation_M(oxidationM);
			diannPar.setAcety_proN(acetyProN);
			diannPar.setPhospho_STY(phosphoSTY);
			diannPar.setK_GG(kGG);
			logger.info(
					"Applied modifications to DiannParameter: M_excision={}, Carbami_C={}, Oxidation_M={}, Acety_proN={}, Phospho_STY={}, K_GG={}",
					mExcision, carbamiC, oxidationM, acetyProN, phosphoSTY, kGG);
		}
	}

	@Override
	protected String getTaskName() {
		return taskName;
	}

	/**
	 * Builds {@code libraries/Marker4/Marker4.speclib}, parquet, and {@code Marker4_info.txt} using the same logic as
	 * the Marker DIA library step, so Seed users can create the library without running a full Marker workflow.
	 *
	 * @param parameterSource thread count and modification flags (e.g. current Seed UI parameters)
	 * @return false if DiaNN is unavailable, inputs are null, or generation failed
	 */
	public static boolean generateMarker4LibraryOnly(MetaParameterDia parameterSource, MetaResourceParameters mrp,
			MagDbItem magDb, javafx.scene.control.ProgressBar bar1, javafx.scene.control.ProgressBar bar2) {
		if (parameterSource == null || mrp == null || magDb == null) {
			return false;
		}
		String diaNNPath = mrp.getToolPath(ToolType.DIANN);
		if (diaNNPath == null || diaNNPath.trim().isEmpty() || !new File(diaNNPath).exists()) {
			return false;
		}
		MetaParameterDiaMarker markerPar = buildMarkerParameterStubForLibraryGeneration(parameterSource, magDb);
		try {
			DiaMarkerCoreTask task = new DiaMarkerCoreTask(markerPar, mrp, bar1, bar2);
			if (task.diaNNTask == null) {
				return false;
			}
			return task.createSSLibrary();
		} catch (Exception e) {
			org.apache.logging.log4j.LogManager.getLogger(DiaMarkerCoreTask.class).error("generateMarker4LibraryOnly failed",
					e);
			return false;
		}
	}

	private static MetaParameterDiaMarker buildMarkerParameterStubForLibraryGeneration(MetaParameterDia src,
			MagDbItem magDb) {
		MetaParameterDiaMarker m = new MetaParameterDiaMarker();
		m.setMagDbItem(magDb);
		m.setThreadCount(src.getThreadCount());
		m.setMExcision(src.isMExcision());
		m.setCarbamiC(src.isCarbamiC());
		m.setOxidationM(src.isOxidationM());
		m.setAcetyProN(src.isAcetyProN());
		m.setPhosphoSTY(src.isPhosphoSTY());
		m.setKGG(src.isKGG());
		String res = src.getResult();
		if (res == null || res.trim().isEmpty()) {
			res = System.getProperty("java.io.tmpdir");
		}
		m.setResult(res);
		if (src.getMetadata() != null) {
			m.setMetadata(src.getMetadata());
		}
		return m;
	}

	@Override
	protected Logger getLogger() {
		return logger;
	}

	@Override
	protected void initResultFolderIfNeeded() {
		if (metaPar == null) {
			return;
		}
		if (metaPar.getDbSearchResultFile() != null) {
			this.resultFolderFile = metaPar.getDbSearchResultFile();
			return;
		}
		if (!(metaPar instanceof MetaParameterDiaMarker)) {
			// Fallback to base behavior for non-focus parameters
			super.initResultFolderIfNeeded();
			return;
		}
		MetaParameterDiaMarker markerParams = (MetaParameterDiaMarker) metaPar;
		// Use the canonical mode selector; persisted values are normalized by the parameter class.
		DiaRefinementMode markerMode = markerParams.getMarkerRefinementMode();
		String folderName;
		if (markerMode.usesComprehensiveSelection()) {
			folderName = "marker_comprehensive";
		} else {
			folderName = "marker_focus";
		}
		File base = new File(markerParams.getResult() != null ? markerParams.getResult() : ".");
		File dir = new File(base, folderName);
		if (!dir.exists()) {
			dir.mkdirs();
			logger.info("Created result directory for DIA Marker workflow ({}): {}", folderName, dir.getAbsolutePath());
		}
		markerParams.setDbSearchResultFile(dir);
		this.resultFolderFile = dir;
	}

	@Override
	protected Boolean call() throws Exception {
		updateMessage(getTaskName());
		updateProgress(0, 100);

		logger.info(getTaskName() + ": started");
		System.out.println(
				new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss").format(new Date()) + "\t" + getTaskName() + ": started");

		// Log detailed parameter information
		logger.info("DIA Focus workflow parameters:");
		logger.info("  Result directory: " + metaPar.getResult());
		logger.info("  Thread count: " + metaPar.getThreadCount());

		if (metadata != null) {
			String[] rawFiles = metadata.getRawFiles();
			logger.info("  Raw files count: " + (rawFiles != null ? rawFiles.length : 0));
			if (rawFiles != null) {
				for (int i = 0; i < rawFiles.length; i++) {
					logger.info("    Raw file " + (i + 1) + ": " + rawFiles[i]);
				}
			}
		}

		if (getClass().equals(DiaMarkerCoreTask.class)) {
			return callViaCoreRunner();
		}

		// Step 0: Copy report templates in background
		updateProgress(5, 100);
		updateMessage(getTaskName() + ": Initializing analysis...");
		if (Platform.isFxApplicationThread()) {
			secondBar.setProgress(ProgressBar.INDETERMINATE_PROGRESS);
		} else {
			Platform.runLater(() -> secondBar.setProgress(ProgressBar.INDETERMINATE_PROGRESS));
		}

		if (!copyReportTemplatesInBackground()) {
			logger.error(getTaskName() + ": Failed to copy report templates");
			return false;
		}

		// Step 1: Spectra Conversion (5% - 25%)
		logger.info("Step 1: Starting spectra format conversion...");
		updateMessage(getTaskName() + ": Converting spectra format...");
		updateProgress(10, 100);
		if (Platform.isFxApplicationThread()) {
			secondBar.setProgress(ProgressBar.INDETERMINATE_PROGRESS);
		} else {
			Platform.runLater(() -> secondBar.setProgress(ProgressBar.INDETERMINATE_PROGRESS));
		}

		boolean convert = spectraConvert();
		if (!convert) {
			logger.error(getTaskName() + ": error in spectra format conversion");
			System.err.println(new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss").format(new Date()) + "\t" + getTaskName()
					+ ": error in spectra format conversion");
			updateMessage(getTaskName() + ": Failed - Spectra conversion error");
			if (Platform.isFxApplicationThread()) {
				secondBar.setProgress(0.0);
			} else {
				Platform.runLater(() -> secondBar.setProgress(0.0));
			}
			return false;
		}
		updateProgress(25, 100);
		logger.info("Step 1: Spectra format conversion completed successfully");

		// Step 2: Focus library setup (25% - 35%)
		logger.info("Step 2: Using pre-generated focus library...");
		updateMessage(getTaskName() + ": Using focus library...");
		updateProgress(30, 100);
		if (Platform.isFxApplicationThread()) {
			secondBar.setProgress(ProgressBar.INDETERMINATE_PROGRESS);
		} else {
			Platform.runLater(() -> secondBar.setProgress(ProgressBar.INDETERMINATE_PROGRESS));
		}

		boolean useFocusLibrary = createSSLibrary();
		if (!useFocusLibrary) {
			logger.error(getTaskName() + ": error in focus library setup");
			System.err.println(new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss").format(new Date()) + "\t" + getTaskName()
					+ ": error in focus library setup");
			updateMessage(getTaskName() + ": Failed - Focus library setup error");
			if (Platform.isFxApplicationThread()) {
				secondBar.setProgress(0.0);
			} else {
				Platform.runLater(() -> secondBar.setProgress(0.0));
			}
			return false;
		}
		updateProgress(35, 100);
		logger.info("Step 2: Focus library setup completed successfully");

		// Step 3: Iterative search - genome selection and combined library generation
		// (35% - 50%)
		logger.info("Step 3: Starting iterative search for genome selection...");
		updateMessage(getTaskName() + ": Performing iterative search...");
		updateProgress(37, 100);
		if (Platform.isFxApplicationThread()) {
			secondBar.setProgress(ProgressBar.INDETERMINATE_PROGRESS);
		} else {
			Platform.runLater(() -> secondBar.setProgress(ProgressBar.INDETERMINATE_PROGRESS));
		}

		boolean iterative = iterativeSearch();
		if (!iterative) {
			logger.error(getTaskName() + ": error in iterative search");
			System.err.println(new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss").format(new Date()) + "\t" + getTaskName()
					+ ": error in iterative search");
			updateMessage(getTaskName() + ": Failed - Iterative search error");
			if (Platform.isFxApplicationThread()) {
				secondBar.setProgress(0.0);
			} else {
				Platform.runLater(() -> secondBar.setProgress(0.0));
			}
			return false;
		}
		updateProgress(50, 100);
		logger.info("Step 3: Iterative search completed successfully");

		// Step 4: Identification using combined library (50% - 75%)
		logger.info("Step 4: Starting identification using combined library...");
		updateMessage(getTaskName() + ": Running identification with combined library...");
		updateProgress(55, 100);
		if (Platform.isFxApplicationThread()) {
			secondBar.setProgress(ProgressBar.INDETERMINATE_PROGRESS);
		} else {
			Platform.runLater(() -> secondBar.setProgress(ProgressBar.INDETERMINATE_PROGRESS));
		}

		boolean iden = iden();
		if (!iden) {
			logger.error(getTaskName() + ": error in identification");
			System.err.println(new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss").format(new Date()) + "\t" + getTaskName()
					+ ": error in identification");
			updateMessage(getTaskName() + ": Failed - Identification error");
			if (Platform.isFxApplicationThread()) {
				secondBar.setProgress(0.0);
			} else {
				Platform.runLater(() -> secondBar.setProgress(0.0));
			}
			return false;
		}
		updateProgress(75, 100);
		logger.info("Step 4: Identification completed successfully");

		// Step 5: Export Report (75% - 100%)
		logger.info("Step 5: Starting report generation...");
		updateMessage(getTaskName() + ": Generating reports...");
		updateProgress(80, 100);
		if (Platform.isFxApplicationThread()) {
			secondBar.setProgress(ProgressBar.INDETERMINATE_PROGRESS);
		} else {
			Platform.runLater(() -> secondBar.setProgress(ProgressBar.INDETERMINATE_PROGRESS));
		}

		boolean exportReport = exportReport();
		if (!exportReport) {
			logger.error(getTaskName() + ": error in report generation");
			System.err.println(new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss").format(new Date()) + "\t" + getTaskName()
					+ ": error in report generation");
			updateMessage(getTaskName() + ": Failed - Report generation error");
			if (Platform.isFxApplicationThread()) {
				secondBar.setProgress(0.0);
			} else {
				Platform.runLater(() -> secondBar.setProgress(0.0));
			}
			return false;
		}
		updateProgress(95, 100);
		logger.info("Step 4: Report generation completed successfully");

		// Update progress for completion
		updateProgress(100, 100);
		updateMessage(getTaskName() + ": Analysis completed successfully!");
		if (Platform.isFxApplicationThread()) {
			secondBar.setProgress(1.0);
		} else {
			Platform.runLater(() -> secondBar.setProgress(1.0));
		}

		logger.info(getTaskName() + ": finished - analysis completed successfully");
		System.out.println(new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss").format(new Date()) + "\t" + getTaskName()
				+ ": finished - analysis completed successfully");
		return true;
	}

	private Boolean callViaCoreRunner() throws Exception {
		updateProgress(5, 100);
		updateMessage(getTaskName() + ": Initializing analysis...");
		updateSecondaryProgress(ProgressBar.INDETERMINATE_PROGRESS);

		if (!copyReportTemplatesInBackground()) {
			logger.error(getTaskName() + ": Failed to copy report templates");
			updateSecondaryProgress(0.0);
			return false;
		}

		updateMessage(getTaskName() + ": Converting spectra format...");
		updateProgress(10, 100);
		updateSecondaryProgress(ProgressBar.INDETERMINATE_PROGRESS);
		if (!spectraConvert()) {
			logger.error(getTaskName() + ": error in spectra format conversion");
			updateMessage(getTaskName() + ": Failed - Spectra conversion error");
			updateSecondaryProgress(0.0);
			return false;
		}
		updateProgress(25, 100);

		File diannExe = resolveGuiToolFile(ToolType.DIANN);
		if (diannExe == null || !diannExe.exists()) {
			logger.error(getTaskName() + ": DIA-NN executable not found");
			updateMessage(getTaskName() + ": Failed - DIA-NN executable not found");
			updateSecondaryProgress(0.0);
			return false;
		}

		try {
			DiaMarkerHeadlessRunner.Result result = DiaMarkerHeadlessRunner.runWithoutPostAnalysis(
					buildCoreRunnerParams(),
					diannExe,
					MagDbItemResolver.resolveConfigFile(null),
					createGuiProgressCallback(),
					inferCondaPrefixFromPython(resolveGuiToolFile(ToolType.PYTHON)),
					resolveGuiToolFile(ToolType.FUNCDB));

			this.resultFolderFile = result.resultDir;
			this.metaPar.setDbSearchResultFile(result.resultDir);
			this.combinedFile = result.combinedDir;
			this.quan_pep_file = result.prMatrixTsv;
			this.quan_pro_file = result.pgMatrixTsv;

			updateMessage(getTaskName() + ": Generating reports...");
			updateProgress(80, 100);
			updateSecondaryProgress(ProgressBar.INDETERMINATE_PROGRESS);
			if (!exportReport()) {
				logger.error(getTaskName() + ": error in report generation");
				updateMessage(getTaskName() + ": Failed - Report generation error");
				updateSecondaryProgress(0.0);
				return false;
			}

			updateProgress(100, 100);
			updateMessage(getTaskName() + ": Analysis completed successfully!");
			updateSecondaryProgress(1.0);
			return true;
		} catch (Exception e) {
			logger.error(getTaskName() + ": core runner execution failed", e);
			updateMessage(getTaskName() + ": Failed - " + e.getMessage());
			updateSecondaryProgress(0.0);
			throw e;
		}
	}

	protected ProgressCallback createGuiProgressCallback() {
		return MetaDiaMarkerCoreRunnerBridge.createGuiProgressCallback(this::updateProgress, this::updateMessage,
				logger, getTaskName());
	}

	protected MetaParameterDiaMarker buildCoreRunnerParams() {
		return MetaDiaMarkerCoreRunnerBridge.buildCoreRunnerParams((MetaParameterDiaMarker) metaPar, metadata,
				diaFiles);
	}

	protected File resolveGuiToolFile(ToolType toolType) {
		return MetaDiaMarkerCoreRunnerBridge.resolveGuiToolFile(mrp, toolType);
	}

	protected File inferCondaPrefixFromPython(File pythonExe) {
		return MetaDiaMarkerCoreRunnerBridge.inferCondaPrefixFromPython(pythonExe);
	}

	/**
	 * Libraries of the final combined search: {@code microbialLibrary}, plus the host .parquet when a host is
	 * configured. Host and microbial libraries are searched together only in this final search, in every DIA
	 * workflow ({@link DiaHostLibrarySupport}); the export then removes the host peptides. Returns null on error.
	 */
	protected String[] finalSearchLibraries(File microbialLibrary, File combinedDir) {
		try {
			File diannExe = resolveGuiToolFile(ToolType.DIANN);
			File condaPrefix = inferCondaPrefixFromPython(resolveGuiToolFile(ToolType.PYTHON));
			ProgressCallback cb = createGuiProgressCallback();
			MetaParameterDiaMarker runnerParams = buildCoreRunnerParams();
			List<File> hosts = resolveHostLibrariesForFinalSearch(combinedDir.getParentFile(), diannExe, cb,
					condaPrefix);
			if (!hosts.isEmpty() && (diannExe == null || !diannExe.exists())) {
				logger.error(getTaskName() + ": DIA-NN executable not found for the combined search with the host library");
				return null;
			}
			return DiaHostLibrarySupport.finalSearchLibraries(microbialLibrary, hosts, combinedDir,
					DiaMarkerHeadlessRunner.libraryBuildPar(runnerParams), diannExe, cb, condaPrefix);
		} catch (Exception e) {
			logger.error(getTaskName() + ": preparing the combined search libraries failed", e);
			updateMessage(getTaskName() + ": Failed - host library error: " + e.getMessage());
			return null;
		}
	}

	/**
	 * Host libraries for the final combined search: {@code host_database} with {@code append_host_db}; empty without.
	 * DIA Seed adds its Standalone libraries.
	 */
	protected List<File> resolveHostLibrariesForFinalSearch(File resultDir, File diannExe, ProgressCallback cb,
			File condaPrefix) throws Exception {
		return DiaMarkerHeadlessRunner.resolveHostParquets(buildCoreRunnerParams(),
				new File(resultDir, "host_library"), diannExe, cb, condaPrefix);
	}

	@Override
	protected boolean createSSLibrary() {
		logger.info("Starting focus library setup for DIA Focus workflow...");

		// Ensure Marker4 files (Marker4.fasta or Marker4_highpep.fasta, genomes-taxa.tsv)
		updateMessage(getTaskName() + ": Checking Marker4 files...");
		if (!ensureMarker4FilesIfNeeded()) {
			logger.error(getTaskName() + ": Failed to ensure Marker4 files");
			return false;
		}

		File fastaFile = magDb.getMarker4FastaForLibrary();
		File libFile = magDb.getMarker4LibFile();
		if (fastaFile == null || libFile == null) {
			logger.error(getTaskName() + ": Marker4 FASTA or library path is null");
			return false;
		}
		if (!fastaFile.exists()) {
			logger.error(getTaskName() + ": Marker4 FASTA not found: {}", fastaFile.getAbsolutePath());
			return false;
		}

		magDb.moveLegacyMarker4SpeclibToCanonical();

		if (!magDb.ensureMarker4LibraryFolderExists()) {
			logger.error(getTaskName() + ": Failed to create libraries/Marker4 folder");
			return false;
		}

		File diannExe = null;
		try {
			String diannPath = mrp != null ? mrp.getToolPath(ToolType.DIANN) : null;
			if (diannPath != null && !diannPath.isBlank()) {
				diannExe = new File(diannPath);
			}
		} catch (Exception ignored) {
			// version probe is best-effort
		}

		if (!(metaPar instanceof MetaParameterDia)) {
			logger.error(getTaskName() + ": metaPar is not MetaParameterDia; cannot finalize Marker4 library metadata");
			return false;
		}
		MetaParameterDia diaMeta = (MetaParameterDia) metaPar;

		Marker4LibrarySupport.logLibraryInUse(fastaFile, createGuiProgressCallback(), logger);
		try {
			Marker4LibrarySupport.invalidateSpeclibIfIdentityMismatch(magDb, fastaFile, diaMeta, diannExe,
					createGuiProgressCallback(), logger);
		} catch (IOException e) {
			logger.error(getTaskName() + ": Failed while checking Marker4.speclib identity", e);
			return false;
		}

		// If library file doesn't exist, generate it from FASTA
		if (!libFile.exists()) {
			logger.info("Focus library file does not exist, generating from {}...", fastaFile.getName());
			System.out.println(new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss").format(new Date()) + "\t" + getTaskName()
					+ ": generating focus library from FASTA");

			logger.info("Using FASTA file: {}", fastaFile.getAbsolutePath());
			logger.info("Generating library file: {}", libFile.getAbsolutePath());

			// Generate library using DiaNNTask
			diaNNTask.addTask(diannPar, fastaFile.getAbsolutePath(), libFile.getAbsolutePath(), false);
			diaNNTask.run(2);

			// Check if a temp file was created (DiaNN sometimes creates .predicted.speclib)
			File libDir = libFile.getParentFile();
			String libName = libFile.getName();
			if (libName.endsWith(".speclib")) {
				libName = libName.substring(0, libName.length() - ".speclib".length());
			}
			File tempLibFile = new File(libDir, libName + ".predicted.speclib");

			if (tempLibFile.exists() && !libFile.exists()) {
				logger.info("Renaming temp library file to final library file");
				tempLibFile.renameTo(libFile);
			}

			// Verify the library file was created
			if (!libFile.exists()) {
				logger.error(getTaskName() + ": Failed to generate focus library file: " + libFile.getAbsolutePath());
				System.err.println(new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss").format(new Date()) + "\t" + getTaskName()
						+ ": Failed to generate focus library file");
				return false;
			}
		}

		try {
			boolean metadataReady = Marker4LibrarySupport.ensureMarker4LibraryMetadata(magDb, diaMeta, fastaFile,
					diannExe, createGuiProgressCallback(), logger, false,
					(metadataPar, speclib, parquet, libFolder, cb, log) -> {
						diaNNTask.addTask(metadataPar, speclib.getAbsolutePath(), parquet.getAbsolutePath());
						diaNNTask.run(1);
						return parquet.exists();
					});
			if (!metadataReady) {
				logger.error(getTaskName() + ": Failed to build Marker4 parquet/info (library metadata)");
				return false;
			}
		} catch (Exception e) {
			logger.error(getTaskName() + ": Failed to build Marker4 parquet/info (library metadata)", e);
			return false;
		}

		logger.info("Focus library validated successfully: {}", libFile.getAbsolutePath());
		logger.info("Focus library file size: {} bytes", libFile.length());

		logger.info("DIA Focus library setup completed successfully");
		System.out.println(new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss").format(new Date()) + "\t" + getTaskName()
				+ ": focus library setup finished");

		return true;
	}

	/**
	 * Ensures Marker4.fasta, genomes-taxa.tsv, and (when required) Marker4_highpep.fasta via the
	 * shared headless preparer. Fails hard when the catalogue is above the highpep threshold and
	 * DeepDetect is not configured.
	 */
	private boolean ensureMarker4FilesIfNeeded() {
		if (magDb == null) {
			logger.error("MAG database is null");
			return false;
		}
		try {
			String exePath = mrp != null ? mrp.getToolPath(ToolType.DEEPDETECT) : null;
			File deepDetectTool = exePath != null && !exePath.isEmpty() ? new File(exePath) : null;
			String pythonPath = mrp != null ? mrp.getToolPath(ToolType.PYTHON) : null;
			File pythonExe = pythonPath != null && !pythonPath.isEmpty() ? new File(pythonPath) : null;
			Marker4ResourcePreparer.DeepDetectPredictor predictor =
					Marker4ResourcePreparer.predictorFromTool(deepDetectTool, pythonExe);
			ProgressCallback cb = createGuiProgressCallback();
			Marker4ResourcePreparer.ensureLibraryResources(magDb, cb, predictor);
			File effectiveFasta = magDb.getMarker4FastaForLibrary();
			logger.info("Marker4 files ready: FASTA={}",
					effectiveFasta != null ? effectiveFasta.getName() : "null");
			return true;
		} catch (IOException e) {
			logger.error("Failed to ensure Marker4 library resources: {}", e.getMessage());
			updateMessage(getTaskName() + ": " + e.getMessage());
			return false;
		}
	}

	@Override
	protected boolean iden() {
		logger.info("Starting identification process using combined library...");
		updateMessage(getTaskName() + ": peptide and protein identification with combined library");
		System.out.println(new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss").format(new Date()) + "\t" + getTaskName()
				+ ": peptide and protein identification started");

		combinedFile = new File(metaPar.getDbSearchResultFile(), "combined");
		if (!combinedFile.exists())
			combinedFile.mkdir();
		this.quan_pep_file = new File(combinedFile, "combined.pr_matrix.tsv");
		this.quan_pro_file = new File(combinedFile, "combined.pg_matrix.tsv");
		File combineParquetFile = new File(combinedFile, "combined.parquet");

		File markerFile = new File(metaPar.getDbSearchResultFile(), "marker");
		File combinedLibFile = new File(markerFile, "combined.speclib");
		if (!combinedLibFile.exists()) {
			logger.error(getTaskName() + ": Combined library file not found: " + combinedLibFile.getAbsolutePath());
			return false;
		}
		diaLibSearchPar.setMbr(true);
		diaLibSearchPar.setDiaFiles(diaFiles);
		diaLibSearchPar.setThreadCount(threadCount);
		diaLibSearchPar.setqValue(0.01);
		diannPar.setThreads(threadCount);
		diannPar.setMatrices(true);
		diannPar.setRelaxed_prot_inf(true);
		diannPar.setReanalyse(true);
		String[] finalLibs = finalSearchLibraries(combinedLibFile, combinedFile);
		if (finalLibs == null) {
			return false;
		}
		if (!quan_pep_file.exists() || !combineParquetFile.exists()) {
			diannPar.setExportDecoy(true);
			diaNNTask.addTask(diannPar, diaLibSearchPar, combineParquetFile.getAbsolutePath(), finalLibs);
			diaNNTask.run(diaFiles.length * 12);
		}

		if (!quan_pep_file.exists() || !combineParquetFile.exists()) {
			logCombinedSearchFailed(combineParquetFile, "library search");
			return false;
		}
		logger.info("Identification completed successfully");
		System.out.println(new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss").format(new Date()) + "\t" + getTaskName()
				+ ": identification finished");
		return true;
	}

	/** Returns a filesystem-safe folder name from a DIA file path. */
	protected static String getSafeFolderNameForDiaFile(String diaFilePath) {
		if (diaFilePath == null || diaFilePath.isEmpty()) return "file_unknown";
		String name = new File(diaFilePath).getName();
		int lastDot = name.lastIndexOf('.');
		if (lastDot > 0) name = name.substring(0, lastDot);
		name = name.replaceAll("[\\\\/:*?\"<>|]", "_");
		return name.isEmpty() ? "file_unknown" : name;
	}

	/** Reads genome IDs from a genome_selected.tsv file and adds them to the given set. */
	protected void readSelectedGenomesIntoSet(File selectedFile, HashSet<String> genomeSet) {
		if (selectedFile == null || !selectedFile.exists() || genomeSet == null) return;
		try (BufferedReader reader = new BufferedReader(new FileReader(selectedFile))) {
			reader.readLine(); // header
			String line;
			while ((line = reader.readLine()) != null) {
				line = line.trim();
				if (line.isEmpty() || line.startsWith("#")) break;
				String[] parts = line.split("\t");
				if (parts.length > 0 && !parts[0].trim().isEmpty()) genomeSet.add(parts[0].trim());
			}
		} catch (IOException e) {
			logger.warn("Error reading selected genomes from {}", selectedFile.getAbsolutePath(), e);
		}
	}

	@Override
	protected boolean iterativeSearch() {
		logger.info("Starting iterative search for DIA Focus workflow...");
		System.out.println(new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss").format(new Date()) + "\t" + getTaskName()
				+ ": iterative search started");

		// Step 1: Search the marker library to get initial PSMs
		File markerLibFile = magDb.resolveMarker4LibraryFile();
		if (markerLibFile == null || !markerLibFile.exists()) {
			logger.error(getTaskName() + ": Marker library file not found for iterative search");
			System.err.println(new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss").format(new Date()) + "\t" + getTaskName()
					+ ": Marker library file not found");
			return false;
		}

		// Create combined folder if it doesn't exist
		File markerFile = new File(metaPar.getDbSearchResultFile(), "marker");
		if (!markerFile.exists()) {
			markerFile.mkdir();
		}

		MetaParameterDiaMarker markerParams = (metaPar instanceof MetaParameterDiaMarker) ? (MetaParameterDiaMarker) metaPar : null;
		// Map markerMode to MarkerGenomeSelector knee mode (FOCUS vs COMPREHENSIVE)
		DiaRefinementMode markerMode = markerParams != null ? markerParams.getMarkerRefinementMode() : DiaRefinementMode.CORE;
		String mode;
		if (markerMode.usesComprehensiveSelection()) {
			mode = MarkerGenomeSelector.KNEE_MODE_COMPREHENSIVE;
		} else {
			mode = MarkerGenomeSelector.KNEE_MODE_FOCUS;
		}

		File currentFile = magDb.getCurrentFile();
		File markerFasta = new File(currentFile, "Marker4.fasta");
		File genomesTaxaTsv = new File(currentFile, "genomes-taxa.tsv");
		if (!markerFasta.exists() || !genomesTaxaTsv.exists()) {
			logger.error(getTaskName() + ": Marker FASTA or genomes-taxa.tsv not found for genome selection");
			return false;
		}

		File combinedFastaFile = new File(markerFile, "combined.fasta");
		HashSet<String> genomeSet = new HashSet<>();

		if (!combinedFastaFile.exists()) {
				File tempParquetFile = new File(markerFile, "marker.parquet");
				if (!tempParquetFile.exists()) {
					DiaNNNoQuantLibrarySearchService.configureNoQuantSearch(
							diannPar, diaLibSearchPar, diaFiles, threadCount, true, true);
					Marker4LibrarySupport.logLibraryInUse(magDb.getMarker4FastaForLibrary(),
							createGuiProgressCallback(), logger);
					logger.info("Searching marker library against DIA files...");
					diaNNTask.addNoQuantLibTask(diannPar, diaLibSearchPar, tempParquetFile.getAbsolutePath(),
							markerLibFile.getAbsolutePath());
					diaNNTask.run(diaFiles.length * 12);
				}
				DiaNNPrecursor[] precursors;
				if (tempParquetFile.exists()) {
					DiaNNParquetReader parquetReader = new DiaNNParquetReader(
							tempParquetFile.getAbsolutePath(), diannPar.getCut());
					precursors = parquetReader.getDiaNNPrecursors(0.01);
				} else {
					logCombinedSearchFailed(tempParquetFile, "marker library search");
					return false;
				}
				if (precursors == null || precursors.length == 0) {
					logger.warn(getTaskName() + ": No precursors found in iterative search, cannot select genomes");
					return false;
				}
				logger.info("Selecting genomes using MarkerGenomeSelector ({} mode)...", mode);
				File evidenceFile = new File(markerFile, "genome_evidence.tsv");
				File rankedFile = new File(markerFile, "genome_ranked.tsv");
				File selectedFile = new File(markerFile, "genome_selected.tsv");
				MagGenomeSelectionService.Request selectionRequest = new MagGenomeSelectionService.Request();
				selectionRequest.magDb = magDb;
				selectionRequest.genomeSelectorMode = MagDbItem.MARKER4;
				selectionRequest.genomeSelectionStrategy = mode;
				selectionRequest.diaPrecursors = precursors;
				selectionRequest.markerFasta = markerFasta;
				selectionRequest.genomesTaxaTsv = genomesTaxaTsv;
				selectionRequest.evidenceFile = evidenceFile;
				selectionRequest.rankedFile = rankedFile;
				selectionRequest.selectedFile = selectedFile;
				selectionRequest.maximumSelected = 300;
				try {
					MagGenomeSelectionService.selectGenomes(selectionRequest, createGuiProgressCallback());
				} catch (Exception e) {
					logger.error(getTaskName() + ": genome selection failed", e);
					return false;
				}
				readSelectedGenomesIntoSet(selectedFile, genomeSet);

			if (genomeSet.isEmpty()) {
				logger.warn(getTaskName() + ": No genomes were selected, cannot create combined database");
				return false;
			}
			logger.info("Selected {} genomes for combined database", genomeSet.size());
			logger.info("Creating combined FASTA file from {} genomes...", genomeSet.size());

			try (PrintWriter writer = new PrintWriter(combinedFastaFile)) {
				File dbFolder = magDb.getOriginalDbFolder();
				for (String genome : genomeSet) {
					File genomeDbFile = new File(dbFolder, genome + ".faa");
					if (!genomeDbFile.exists()) {
						genomeDbFile = new File(dbFolder, genome + ".fasta");
					}
					if (genomeDbFile.exists()) {
						try (BufferedReader reader = new BufferedReader(new FileReader(genomeDbFile))) {
							String line;
							while ((line = reader.readLine()) != null) {
								writer.println(line);
							}
						} catch (IOException e) {
							logger.error("Error reading genome file: {}", genomeDbFile.getAbsolutePath(), e);
						}
					} else {
						logger.warn("Genome file not found: {}", genome);
					}
				}
			} catch (IOException e) {
				logger.error("Error writing combined FASTA file: {}", combinedFastaFile.getAbsolutePath(), e);
				return false;
			}
		} else {
			logger.info("Combined fasta file already exists in {}", combinedFastaFile.getAbsolutePath());
		}

		if (!combinedFastaFile.exists() || combinedFastaFile.length() == 0) {
			logger.error(getTaskName() + ": Failed to create combined FASTA file");
			return false;
		}

		logger.info("Combined FASTA file created: {} ({} bytes)", combinedFastaFile.getAbsolutePath(),
				combinedFastaFile.length());

		// Step 6: Generate library from combined FASTA
		File combinedLibFile = new File(markerFile, "combined.speclib");
		logger.info("Generating combined library from FASTA...");
		
		if (!combinedLibFile.exists()) {
			File altLibFile = new File(markerFile, "combined.predicted.speclib");
			if (altLibFile.exists()) {
				altLibFile.renameTo(combinedLibFile);
			} else {
				diannPar.setThreads(threadCount);
				diaNNTask.addTask(diannPar, combinedFastaFile.getAbsolutePath(), combinedLibFile.getAbsolutePath(),
						false);
				diaNNTask.run(12);

				// Check if library was created (DiaNN might create it with .speclib extension)
				if (!combinedLibFile.exists()) {
					// Try without .predicted
					if (altLibFile.exists()) {
						altLibFile.renameTo(combinedLibFile);
					}
				}
			}
		} else {
			logger.info("Combined library file already exists in {}", combinedLibFile.getAbsolutePath());
		}

		if (!combinedLibFile.exists()) {
			logger.error(getTaskName() + ": Failed to generate combined library file");
			System.err.println(new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss").format(new Date()) + "\t" + getTaskName()
					+ ": Failed to generate combined library file");
			return false;
		}

		logger.info("Combined library generated successfully: {} ({} bytes)", combinedLibFile.getAbsolutePath(),
				combinedLibFile.length());
		logger.info("Iterative search completed successfully");
		System.out.println(new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss").format(new Date()) + "\t" + getTaskName()
				+ ": iterative search finished");

		return true;
	}

	/**
	 * Folder containing iterative search results (combined.fasta, etc.). Override in subclasses to use a different folder (e.g. seed_iterative).
	 */
	protected File getIterativeResultFolder() {
		return new File(metaPar.getDbSearchResultFile(), "marker");
	}

	protected void enrichPrecursorsWithPeptideProteinMapping(DiaNNPrecursor[] precursors, File parquetFolder) {
		// Default marker workflow: no additional peptide->protein remapping is required.
	}

	@Override
	protected File getExportFastaFolder() {
		return getIterativeResultFolder();
	}

	/** Folder for combined.parquet / combined.pr_matrix.tsv; must be set in iden() before exportReport(). */
	@Override
	protected File getExportParquetFolder() {
		return combinedFile;
	}

	@Override
	protected boolean exportReport() {
		return super.exportReport();
	}

	/**
	 * Run the full downstream report for any MetaDiaTask. Used by DiaSeedDirectTask.
	 * For Marker, runs the full implementation; for others (e.g. SeedDirect), returns true so that
	 * at least the new genome/protein score TSVs are considered success (full HTML report is Marker-only for now).
	 */
	protected static boolean runExportReportDownstream(MetaDiaTask task,
			java.util.HashMap<String, double[]> genomeScoreMap, java.util.HashMap<String, double[]> proteinScoreMap) {
		if (task instanceof DiaMarkerCoreTask)
			return ((DiaMarkerCoreTask) task).exportReportDownstreamFromNewScoring(genomeScoreMap, proteinScoreMap);
		if (task instanceof DiaSeedDirectTask)
			return runFullDownstreamFromNewScoring(task, genomeScoreMap, proteinScoreMap);
		return true;
	}

	/**
	 * Full downstream report (peptide report, intensity maps, taxonomy, HTML) for any MetaDiaTask.
	 * Used by DiaMarkerCoreTask and DiaSeedDirectTask.
	 */
	protected static boolean runFullDownstreamFromNewScoring(MetaDiaTask diaTask,
			java.util.HashMap<String, double[]> genomeScoreMap, java.util.HashMap<String, double[]> proteinScoreMap) {
        return diaTask.exportPreparedReport(diaTask.diaPp);
    }

	@Override
	protected boolean exportReportDownstreamFromNewScoring(java.util.HashMap<String, double[]> genomeScoreMap, java.util.HashMap<String, double[]> proteinScoreMap) {
		return runFullDownstreamFromNewScoring(this, genomeScoreMap, proteinScoreMap);
	}

}
