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
import java.util.Map;
import java.util.concurrent.ExecutionException;

import org.apache.commons.lang3.math.NumberUtils;
import org.apache.logging.log4j.Logger;

import javafx.application.Platform;
import javafx.scene.control.ProgressBar;

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
import uk.ac.quadram.metapilot.dbSearch.deepDetect.DeepDetectParameter;
import uk.ac.quadram.metapilot.dbSearch.deepDetect.DeepDetectTask;
import uk.ac.quadram.metapilot.core.MGYG.HapExtractor;
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
	private static final long LOW_MEMORY_EXPORT_PARQUET_THRESHOLD_BYTES = 1L * 1024L * 1024L * 1024L;

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

		if (!runOptionalHostSearch()) {
			return false;
		}

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

	protected boolean runOptionalHostSearch() {
		if (metaPar == null || !metaPar.isAppendHostDb()) {
			return true;
		}
		try {
			updateMessage(getTaskName() + ": Running optional host DB search...");
			File diannExe = resolveGuiToolFile(ToolType.DIANN);
			if (diannExe == null || !diannExe.exists()) {
				logger.error(getTaskName() + ": DIA-NN executable not found for host search");
				return false;
			}
			DiaMarkerHeadlessRunner.runHostSearchIfRequested(
					metaPar.getDbSearchResultFile(),
					diannExe,
					buildCoreRunnerParams(),
					diaFiles,
					createGuiProgressCallback(),
					inferCondaPrefixFromPython(resolveGuiToolFile(ToolType.PYTHON)));
			return true;
		} catch (Exception e) {
			logger.error(getTaskName() + ": host DB search failed", e);
			updateMessage(getTaskName() + ": Failed - Host DB search error: " + e.getMessage());
			return false;
		}
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

		if (!(metaPar instanceof MetaParameterDia)) {
			logger.error(getTaskName() + ": metaPar is not MetaParameterDia; cannot finalize Marker4 library metadata");
			return false;
		}
		try {
			boolean metadataReady = Marker4LibrarySupport.ensureMarker4LibraryMetadata(magDb,
					(MetaParameterDia) metaPar, createGuiProgressCallback(), logger, false,
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
			logger.error(getTaskName() + ": Failed to build Marker4 parquet/info (library metadata)");
			return false;
		}

		logger.info("Focus library validated successfully: {}", libFile.getAbsolutePath());
		logger.info("Focus library file size: {} bytes", libFile.length());

		logger.info("DIA Focus library setup completed successfully");
		System.out.println(new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss").format(new Date()) + "\t" + getTaskName()
				+ ": focus library setup finished");

		return true;
	}

	private static final double DEFAULT_DETECTABILITY_THRESHOLD = 0.8;
	private static final int MARKER4_HIGHPEP_PROTEIN_THRESHOLD = 100_000;

	/**
	 * Ensures Marker4.fasta and genomes-taxa.tsv exist. Marker4.fasta is auto-generated when a MAG
	 * database is downloaded; if missing, we extract it from the original DB. If Marker4.fasta has
	 * more than 100,000 proteins, we generate Marker4_highpep.fasta (high-detectability peptides only)
	 * to keep the library manageable; otherwise we use Marker4.fasta directly.
	 *
	 * @return true if Marker4.fasta (or Marker4_highpep.fasta) and genomes-taxa.tsv exist
	 */
	private boolean ensureMarker4FilesIfNeeded() {
		if (magDb == null) {
			logger.error("MAG database is null");
			return false;
		}
		File currentFile = magDb.getCurrentFile();
		if (currentFile == null || !currentFile.exists()) {
			logger.error("MAG database current file is null or does not exist");
			return false;
		}

		// Marker4.fasta is auto-generated on download; extract if missing (e.g. legacy DB)
		File marker4Fasta = new File(currentFile, "Marker4.fasta");
		if (!marker4Fasta.exists()) {
			updateMessage(getTaskName() + ": Extracting Marker4 proteins...");
			if (!ensureMarker4FastaExists(marker4Fasta)) {
				return false;
			}
		}

		// If Marker4.fasta has >100k proteins, generate Marker4_highpep.fasta for a smaller library
		int proteinCount = countFastaProteins(marker4Fasta);
		if (proteinCount < 0) {
			logger.error("Failed to read protein count from Marker4.fasta");
			return false;
		}
		if (proteinCount > MARKER4_HIGHPEP_PROTEIN_THRESHOLD) {
			logger.info("Marker4.fasta has {} proteins (>{}), generating Marker4_highpep.fasta...",
					proteinCount, MARKER4_HIGHPEP_PROTEIN_THRESHOLD);
			File markersHighpepFile = new File(currentFile, "Marker4_highpep.fasta");
			if (!markersHighpepFile.exists()) {
				ensureMarker4HighpepFasta(markersHighpepFile, marker4Fasta);
			}
		} else {
			logger.info("Marker4.fasta has {} proteins (≤{}), using directly for library", proteinCount,
					MARKER4_HIGHPEP_PROTEIN_THRESHOLD);
		}

		// Ensure genomes-taxa.tsv for genome selection
		File genomeTaxaFile = new File(currentFile, "genomes-taxa.tsv");
		if (!genomeTaxaFile.exists()) {
			ensureGenomesTaxaFile(genomeTaxaFile);
		}

		File effectiveFasta = magDb.getMarker4FastaForLibrary();
		if (effectiveFasta == null || !effectiveFasta.exists()) {
			logger.error("Marker4 FASTA not available: {}", effectiveFasta);
			return false;
		}
		if (!genomeTaxaFile.exists()) {
			logger.error("genomes-taxa.tsv not available");
			return false;
		}
		logger.info("Marker4 files ready: FASTA={}, genomes-taxa.tsv exists", effectiveFasta.getName());
		return true;
	}

	private boolean ensureMarker4FastaExists(File marker4Fasta) {
		File originalDbFolder = magDb.getOriginalDbFolder();
		if (originalDbFolder == null || !originalDbFolder.exists()) {
			logger.error("Cannot extract Marker4.fasta: original DB folder not available");
			return false;
		}
		HapExtractor extractor = new HapExtractor();
		try {
			extractor.extractFourCategories(originalDbFolder, marker4Fasta);
		} catch (IOException e) {
			logger.error("Failed to extract Marker4.fasta", e);
			return false;
		}
		if (!marker4Fasta.exists()) {
			logger.error("Marker4.fasta was not created");
			return false;
		}
		return true;
	}

	private void ensureMarker4HighpepFasta(File markersHighpepFile, File marker4Fasta) {
		File currentDir = markersHighpepFile.getParentFile();
		File predictedTsv = new File(currentDir, "Marker4.predicted.tsv");
		if (!predictedTsv.exists()) {
			updateMessage(getTaskName() + ": Running DeepDetect for peptide prediction...");
			String exePath = mrp != null ? mrp.getToolPath(ToolType.DEEPDETECT) : null;
			File deepDetectExe = exePath != null && !exePath.isEmpty() ? new File(exePath) : null;
			if (deepDetectExe == null || !deepDetectExe.exists()) {
				logger.error("DeepDetect executable not found; cannot generate Marker4.predicted.tsv");
				return;
			}
			DeepDetectParameter par = new DeepDetectParameter();
			par.setInput(marker4Fasta);
			par.setOutput(predictedTsv);
			int timeoutHours = 4;
			boolean ok = DeepDetectTask.testDeepDetect(deepDetectExe, par, timeoutHours);
			if (!ok || !predictedTsv.exists()) {
				logger.error("DeepDetect did not produce Marker4.predicted.tsv");
				return;
			}
		}
		if (!markersHighpepFile.exists()) {
			MarkerGenomeSelector.exportPeptideFasta(predictedTsv.getAbsolutePath(),
					markersHighpepFile.getAbsolutePath(), DEFAULT_DETECTABILITY_THRESHOLD);
		}
		if (markersHighpepFile.exists()) {
			logger.info("Generated Marker4_highpep.fasta: {}", markersHighpepFile.getAbsolutePath());
		} else {
			logger.warn("Marker4_highpep.fasta was not written (no peptides above threshold?)");
		}
	}

	private static int countFastaProteins(File fastaFile) {
		int count = 0;
		try (BufferedReader r = new BufferedReader(new FileReader(fastaFile))) {
			String line;
			while ((line = r.readLine()) != null) {
				if (line.startsWith(">")) {
					count++;
				}
			}
		} catch (IOException e) {
			return -1;
		}
		return count;
	}

	private void ensureGenomesTaxaFile(File genomeTaxaFile) {
		try {
			Marker4ResourcePreparer.ensureGenomesTaxaTsv(magDb, genomeTaxaFile, createGuiProgressCallback());
			logger.info("Generated genomes-taxa.tsv: {}", genomeTaxaFile.getAbsolutePath());
		} catch (IOException e) {
			logger.error("Could not generate genomes-taxa.tsv", e);
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
		if (!quan_pep_file.exists() || !combineParquetFile.exists()) {
			diannPar.setExportDecoy(true);
			diaNNTask.addTask(diannPar, diaLibSearchPar, combineParquetFile.getAbsolutePath(),
					new String[] { combinedLibFile.getAbsolutePath() });
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
		java.text.SimpleDateFormat fmt = new java.text.SimpleDateFormat("MM_dd_yyyy_HH_mm_ss");
		File quanPepFile = diaTask.getExportQuanPepFile();
		if (quanPepFile == null || !quanPepFile.exists() || !quanPepFile.canRead()) {
			String path = quanPepFile != null ? quanPepFile.getAbsolutePath() : "null";
			diaTask.getExportLogger().error(diaTask.getTaskName() + ": peptide matrix file missing or unreadable: " + path);
			return false;
		}
		File parquetFile = new File(quanPepFile.getParentFile(), "combined.parquet");
		boolean largeParquetExport = parquetFile.exists()
				&& parquetFile.length() >= LOW_MEMORY_EXPORT_PARQUET_THRESHOLD_BYTES;
		if (largeParquetExport) {
			diaTask.getExportLogger().info(diaTask.getTaskName()
					+ ": using compact low-memory DIA export path because parquet is "
					+ parquetFile.length() + " bytes (threshold "
					+ LOW_MEMORY_EXPORT_PARQUET_THRESHOLD_BYTES + " bytes)");
		}
		DiannPrMatrixReader prReader = new DiannPrMatrixReader(quanPepFile);
		FlashLfqQuanPeptide[] peps = prReader.getMetaPeptides(diaTask.diaPp);
		diaTask.setExportQuanExpNames(prReader.getIntensityTitle());
		// Capture the per-run summary counts while the precursor array is still
		// populated: exportSummary() further down only needs these small maps, so
		// diaPp can be released afterwards without a NullPointerException. Done on
		// every run, not only the low-memory path, so both paths behave identically.
		diaTask.captureSummaryStats();
		if (largeParquetExport) {
			diaTask.diaPp = null;
		}

		HashMap<String, FlashLfqQuanPeptide> modseqPepMap = new HashMap<String, FlashLfqQuanPeptide>();
		for (FlashLfqQuanPeptide pep : peps) {
			modseqPepMap.put(pep.getModSeq(), pep);
		}

		diaTask.getExportLogger().info(diaTask.getTaskName() + ": exporting peptide report started");
		System.out.println(fmt.format(new Date()) + "\t" + diaTask.getTaskName()
				+ ": exporting peptide report started");

		HashMap<String, double[]> genomeIntensityMap = new HashMap<String, double[]>();
		HashMap<String, double[]> proteinIntensityMap = new HashMap<String, double[]>();
		HashMap<String, Double> proteinTotalIntenMap = new HashMap<String, Double>();

		HashSet<String> pepSequenceSet = new HashSet<String>();
		HashMap<String, HashSet<String>> totalPepCountMap = new HashMap<String, HashSet<String>>();
		HashMap<String, HashSet<String>> razorPepCountMap = new HashMap<String, HashSet<String>>();
		ArrayList<FlashLfqQuanPeptide> peptideList = new ArrayList<FlashLfqQuanPeptide>();

		for (int i = 0; i < peps.length; i++) {
			String stripeSequence = peps[i].getSequence();
			String razorProtein = peps[i].getRazorProtein();

			pepSequenceSet.add(stripeSequence);
			peptideList.add(peps[i]);

			if (razorPepCountMap.containsKey(razorProtein)) {
				razorPepCountMap.get(razorProtein).add(stripeSequence);
			} else {
				HashSet<String> set = new HashSet<String>();
				set.add(stripeSequence);
				razorPepCountMap.put(razorProtein, set);
			}

			if (razorProtein.startsWith(diaTask.getExportMagDb().getIdentifier())) {
				String genome = razorProtein.split("_")[0];
				if (razorPepCountMap.containsKey(genome)) {
					razorPepCountMap.get(genome).add(stripeSequence);
				} else {
					HashSet<String> set = new HashSet<String>();
					set.add(stripeSequence);
					razorPepCountMap.put(genome, set);
				}
			}
		}

		HashMap<String, Integer> psmCountMap = new HashMap<String, Integer>();

		int[] filePsmCount = new int[diaTask.getExportQuanExpNames().length];
		HashSet<String>[] fileSequenceSets = new HashSet[diaTask.getExportQuanExpNames().length];
		for (int i = 0; i < fileSequenceSets.length; i++) {
			fileSequenceSets[i] = new HashSet<String>();
		}

		for (int i = 0; i < peps.length; i++) {
			String[] proteins = peps[i].getProteins();
			int count = 0;
			for (int j = 0; j < proteins.length; j++) {
				if (razorPepCountMap.containsKey(proteins[j])) {
					count++;
				}
			}

			if (count > 0) {
				String modSequence = peps[i].getModSeq();
				int[] idenType = peps[i].getIdenType();
				for (int j = 0; j < fileSequenceSets.length; j++) {
					if (idenType[j] == 0) {
						fileSequenceSets[j].add(modSequence);
					}
				}

				String stripeSequence = peps[i].getSequence();
				double[] pepIntensity = new double[peps[i].getIntensity().length];
				for (int j = 0; j < pepIntensity.length; j++) {
					pepIntensity[j] = peps[i].getIntensity()[j] / (double) count;
				}

				int pepMs2Count = 0;
				int[] ms2Count = peps[i].getMs2Counts();
				for (int j = 0; j < ms2Count.length; j++) {
					pepMs2Count += ms2Count[j];
					filePsmCount[j] += ms2Count[j];
				}

				for (int j = 0; j < proteins.length; j++) {
					if (razorPepCountMap.containsKey(proteins[j])) {

						if (psmCountMap.containsKey(proteins[j])) {
							psmCountMap.put(proteins[j], psmCountMap.get(proteins[j]) + pepMs2Count);
						} else {
							psmCountMap.put(proteins[j], pepMs2Count);
						}

						double proteinTotalInten;
						double[] proteinIntensity;
						if (proteinIntensityMap.containsKey(proteins[j])) {
							proteinIntensity = proteinIntensityMap.get(proteins[j]);
							proteinTotalInten = proteinTotalIntenMap.get(proteins[j]);
						} else {
							proteinIntensity = new double[pepIntensity.length];
							proteinIntensityMap.put(proteins[j], proteinIntensity);
							proteinTotalInten = 0.0;
						}
						for (int k = 0; k < pepIntensity.length; k++) {
							proteinIntensity[k] += pepIntensity[k];
							proteinTotalInten += pepIntensity[k];
						}
						proteinTotalIntenMap.put(proteins[j], proteinTotalInten);

						if (totalPepCountMap.containsKey(proteins[j])) {
							totalPepCountMap.get(proteins[j]).add(stripeSequence);
						} else {
							HashSet<String> set = new HashSet<String>();
							set.add(stripeSequence);
							totalPepCountMap.put(proteins[j], set);
						}

						if (proteins[j].startsWith(diaTask.getExportMagDb().getIdentifier())) {
							String genome = proteins[j].split("_")[0];
							if (MetaDiaTask.isConfidentDiaGenomeScore(genomeScoreMap, genome)) {
								double[] genomeIntensity;
								if (genomeIntensityMap.containsKey(genome)) {
									genomeIntensity = genomeIntensityMap.get(genome);
								} else {
									genomeIntensity = new double[pepIntensity.length];
									genomeIntensityMap.put(genome, genomeIntensity);
								}
								for (int k = 0; k < pepIntensity.length; k++) {
									genomeIntensity[k] += pepIntensity[k];
								}

								if (psmCountMap.containsKey(genome)) {
									psmCountMap.put(genome, psmCountMap.get(genome) + pepMs2Count);
								} else {
									psmCountMap.put(genome, pepMs2Count);
								}

								if (totalPepCountMap.containsKey(genome)) {
									totalPepCountMap.get(genome).add(stripeSequence);
								} else {
									HashSet<String> set = new HashSet<String>();
									set.add(stripeSequence);
									totalPepCountMap.put(genome, set);
								}
							}
						}
					}
				}
			}
		}

		// --- Pure razor: only razor protein list, each peptide contributes full intensity to its razor protein only ---
		HashMap<String, double[]> proteinIntensityMapRazor = new HashMap<>();
		HashMap<String, Double> proteinTotalIntenMapRazor = new HashMap<>();
		HashMap<String, Integer> psmCountMapRazor = new HashMap<>();
		for (int i = 0; i < peps.length; i++) {
			String razorProtein = peps[i].getRazorProtein();
			double[] fullIntensity = peps[i].getIntensity();
			int[] ms2Count = peps[i].getMs2Counts();
			int pepMs2Count = 0;
			for (int j = 0; j < ms2Count.length; j++) pepMs2Count += ms2Count[j];
			if (proteinIntensityMapRazor.containsKey(razorProtein)) {
				double[] intensity = proteinIntensityMapRazor.get(razorProtein);
				double total = proteinTotalIntenMapRazor.get(razorProtein);
				for (int k = 0; k < fullIntensity.length; k++) {
					intensity[k] += fullIntensity[k];
					total += fullIntensity[k];
				}
				proteinTotalIntenMapRazor.put(razorProtein, total);
				psmCountMapRazor.put(razorProtein, psmCountMapRazor.get(razorProtein) + pepMs2Count);
			} else {
				double[] intensity = new double[fullIntensity.length];
				double total = 0;
				for (int k = 0; k < fullIntensity.length; k++) {
					intensity[k] = fullIntensity[k];
					total += fullIntensity[k];
				}
				proteinIntensityMapRazor.put(razorProtein, intensity);
				proteinTotalIntenMapRazor.put(razorProtein, total);
				psmCountMapRazor.put(razorProtein, pepMs2Count);
			}
		}

		// --- Pure shared: all proteins, each peptide's intensity divided among all its proteins ---
		HashMap<String, double[]> proteinIntensityMapAll = new HashMap<>();
		HashMap<String, Double> proteinTotalIntenMapAll = new HashMap<>();
		HashMap<String, HashSet<String>> totalPepCountMapAll = new HashMap<>();
		HashMap<String, Integer> psmCountMapAll = new HashMap<>();
		for (int i = 0; i < peps.length; i++) {
			String[] prots = peps[i].getProteins();
			if (prots == null || prots.length == 0) continue;
			String stripeSequence = peps[i].getSequence();
			double[] fullIntensity = peps[i].getIntensity();
			int nProts = prots.length;
			double[] divided = new double[fullIntensity.length];
			for (int j = 0; j < fullIntensity.length; j++) divided[j] = fullIntensity[j] / (double) nProts;
			int[] ms2Count = peps[i].getMs2Counts();
			int pepMs2Count = 0;
			for (int j = 0; j < ms2Count.length; j++) pepMs2Count += ms2Count[j];
			for (String pro : prots) {
				if (proteinIntensityMapAll.containsKey(pro)) {
					double[] intensity = proteinIntensityMapAll.get(pro);
					double total = proteinTotalIntenMapAll.get(pro);
					for (int k = 0; k < divided.length; k++) {
						intensity[k] += divided[k];
						total += divided[k];
					}
					proteinTotalIntenMapAll.put(pro, total);
					totalPepCountMapAll.get(pro).add(stripeSequence);
					psmCountMapAll.put(pro, psmCountMapAll.get(pro) + pepMs2Count);
				} else {
					double[] intensity = new double[divided.length];
					double total = 0;
					for (int k = 0; k < divided.length; k++) {
						intensity[k] = divided[k];
						total += divided[k];
					}
					proteinIntensityMapAll.put(pro, intensity);
					proteinTotalIntenMapAll.put(pro, total);
					HashSet<String> set = new HashSet<>();
					set.add(stripeSequence);
					totalPepCountMapAll.put(pro, set);
					psmCountMapAll.put(pro, pepMs2Count);
				}
			}
		}

		// Ensure every reported protein has a score row (default 0, 1, 1 if missing)
		for (String pro : proteinTotalIntenMap.keySet()) {
			if (!proteinScoreMap.containsKey(pro)) {
				proteinScoreMap.put(pro, new double[] { 0.0, 1.0, 1.0 });
			}
		}
		for (String pro : proteinTotalIntenMapAll.keySet()) {
			if (!proteinScoreMap.containsKey(pro)) {
				proteinScoreMap.put(pro, new double[] { 0.0, 1.0, 1.0 });
			}
		}

		String[] proteins = proteinTotalIntenMap.keySet().toArray(String[]::new);
		Arrays.sort(proteins, new Comparator<String>() {

			@Override
			public int compare(String o1, String o2) {
				// sort by S descending (higher S first), then by total intensity descending
				double[] s1 = proteinScoreMap.get(o1);
				double[] s2 = proteinScoreMap.get(o2);
				if (s1 == null && s2 == null) return 0;
				if (s1 == null) return 1;
				if (s2 == null) return -1;
				if (s1[IDX_S] > s2[IDX_S]) return -1;
				if (s1[IDX_S] < s2[IDX_S]) return 1;
				if (proteinTotalIntenMap.get(o1) > proteinTotalIntenMap.get(o2)) return -1;
				if (proteinTotalIntenMap.get(o1) < proteinTotalIntenMap.get(o2)) return 1;
				return 0;
			}
		});

		HashMap<String, Integer> proteinIdMap = new HashMap<String, Integer>();
		for (int i = 0; i < proteins.length; i++) {
			proteinIdMap.put(proteins[i], i);
		}

		try (PrintWriter idWriter = new PrintWriter(new File(diaTask.getExportResultFolderFile(), "id.tsv"))) {
			idWriter.println("Raw file\tUnique peptide count\tPSM count");
			for (int i = 0; i < diaTask.getExportQuanExpNames().length; i++) {
				idWriter.println(diaTask.getExportQuanExpNames()[i] + "\t" + fileSequenceSets[i].size() + "\t" + filePsmCount[i]);
			}
			idWriter.close();
		} catch (IOException e) {
			System.out.println(fmt.format(new Date()) + "\t" + diaTask.getTaskName()
					+ ": error in writing identification information to " + new File(diaTask.getExportResultFolderFile(), "id.tsv"));
			diaTask.getExportLogger().error(diaTask.getTaskName() + ": error in writing identification information to "
					+ new File(diaTask.getExportResultFolderFile(), "id.tsv"), e);
		}

		diaTask.exportUpdateProgress(85, 100);

		HashMap<String, ArrayList<Integer>> proPepIdMap = new HashMap<String, ArrayList<Integer>>();
		HashMap<String, ArrayList<Integer>> proPepRazorIdMap = new HashMap<String, ArrayList<Integer>>();
		try {
			diaTask.setExportFinalPepTxt(new File(diaTask.getExportResultFolderFile(), "final_peptides.tsv"));

			PrintWriter writer = new PrintWriter(diaTask.getExportFinalPepTxt());

			StringBuilder titlesb = new StringBuilder();
			titlesb.append("Sequence").append("\t");
			titlesb.append("Base sequence").append("\t");
			titlesb.append("Length").append("\t");
			titlesb.append("Missed cleavages").append("\t");
			titlesb.append("Mass").append("\t");
			titlesb.append("Razor protein").append("\t");
			titlesb.append("Mapped proteins").append("\t");
			titlesb.append("Charges").append("\t");
			titlesb.append("Qvalue").append("\t");
			titlesb.append("PEP").append("\t");

			for (int i = 0; i < diaTask.getExportQuanExpNames().length; i++) {
				titlesb.append("Identification type " + diaTask.getExportQuanExpNames()[i]).append("\t");
			}

			for (int i = 0; i < diaTask.getExportQuanExpNames().length; i++) {
				titlesb.append("MS/MS count " + diaTask.getExportQuanExpNames()[i]).append("\t");
			}

			titlesb.append("Intensity").append("\t");

			for (int i = 0; i < diaTask.getExportQuanExpNames().length; i++) {
				titlesb.append("Intensity " + diaTask.getExportQuanExpNames()[i]).append("\t");
			}

			titlesb.append("Reverse").append("\t");
			titlesb.append("Potential contaminant").append("\t");
			titlesb.append("id").append("\t");
			titlesb.append("Protein group IDs").append("\t");
			titlesb.append("MS/MS Count");

			writer.println(titlesb);

			int id = 0;
			for (FlashLfqQuanPeptide peptide : peptideList) {

				String fullSeq = peptide.getModSeq();
				String sequence = peptide.getSequence();

				StringBuilder sb = new StringBuilder();
				sb.append(fullSeq).append("\t");
				sb.append(sequence).append("\t");
				sb.append(sequence.length()).append("\t");
				sb.append(peptide.getMissCleave()).append("\t");
				sb.append(peptide.getMass()).append("\t");
				String razorProtein = peptide.getRazorProtein();
				sb.append(razorProtein != null ? razorProtein : "").append("\t");

				String[] pros = peptide.getProteins();

				for (String pro : pros) {
					if (pro != null) {
						sb.append(pro).append(";");
					}
				}
				sb.deleteCharAt(sb.length() - 1);
				sb.append("\t");

				int[] charge = peptide.getCharges();
				StringBuilder chargesb = new StringBuilder();
				for (int j = 0; j < charge.length; j++) {
					if (charge[j] == 1) {
						chargesb.append(j).append(";");
					}
				}
				if (chargesb.length() > 1) {
					sb.append(chargesb.subSequence(0, chargesb.length() - 1)).append("\t");
				} else {
					sb.append("2").append("\t");
				}

				sb.append(peptide.getScore()).append("\t");
				sb.append(peptide.getPEP()).append("\t");

				int[] idenType = peptide.getIdenType();
				for (int j = 0; j < idenType.length; j++) {
					sb.append(MetaPeptide.idenTypeStrings[idenType[j]]).append("\t");
				}

				int totalMsCount = 0;
				int[] ms2Count = peptide.getMs2Counts();
				for (int j = 0; j < ms2Count.length; j++) {
					totalMsCount += ms2Count[j];
					sb.append(ms2Count[j]).append("\t");
				}

				double[] intensity = peptide.getIntensity();
				double totalIntensity = 0;
				for (int j = 0; j < intensity.length; j++) {
					totalIntensity += intensity[j];
				}

				sb.append(BigInteger.valueOf((long) Math.floor(totalIntensity))).append("\t");

				for (int j = 0; j < intensity.length; j++) {
					sb.append(BigInteger.valueOf((long) Math.floor(intensity[j]))).append("\t");
				}

				sb.append("").append("\t");
				sb.append("").append("\t");

				sb.append(id).append("\t");

				for (String pro : pros) {
					if (proteinIdMap.containsKey(pro)) {
						sb.append(proteinIdMap.get(pro)).append(";");

						if (proPepIdMap.containsKey(pro)) {
							proPepIdMap.get(pro).add(id);
						} else {
							ArrayList<Integer> list = new ArrayList<Integer>();
							list.add(id);
							proPepIdMap.put(pro, list);
						}

						HashSet<String> rSet = razorPepCountMap.get(pro);
						if (proPepRazorIdMap.containsKey(pro)) {
							if (rSet.contains(sequence)) {
								proPepRazorIdMap.get(pro).add(1);
							} else {
								proPepRazorIdMap.get(pro).add(0);
							}
						} else {
							ArrayList<Integer> list = new ArrayList<Integer>();
							if (rSet.contains(sequence)) {
								list.add(1);
							} else {
								list.add(0);
							}
							proPepRazorIdMap.put(pro, list);
						}
					}
				}
				sb.deleteCharAt(sb.length() - 1);
				sb.append("\t");

				sb.append(totalMsCount);

				writer.println(sb);
				id++;
			}

			writer.close();

		} catch (Exception e) {
			// TODO Auto-generated catch block
			diaTask.getExportLogger().error(diaTask.getTaskName() + ": error in writing the peptides to  " + diaTask.getExportFinalPepTxt(), e);
			System.err.println(fmt.format(new Date()) + "\t" + diaTask.getTaskName() + ": error in writing the peptides to  "
					+ diaTask.getExportFinalPepTxt());
		}

		diaTask.getExportLogger().info(diaTask.getTaskName() + ": peptide report has been exported to " + diaTask.getExportFinalPepTxt());
		System.out.println(fmt.format(new Date()) + "\t" + diaTask.getTaskName() + ": peptide report has been exported to "
				+ diaTask.getExportFinalPepTxt());

		if (!diaTask.getExportFinalPepTxt().exists() || diaTask.getExportFinalPepTxt().length() == 0) {
			return false;
		}

		diaTask.exportUpdateProgress(87, 100);

		diaTask.getExportLogger().info(diaTask.getTaskName() + ": exporting protein report started");
		System.out.println(fmt.format(new Date()) + "\t" + diaTask.getTaskName() + ": exporting protein report started");

		ArrayList<MetaProtein> proteinList = new ArrayList<MetaProtein>();
		try {

			diaTask.setExportFinalProTxt(new File(diaTask.getExportResultFolderFile(), "razor_proteins.tsv"));

			PrintWriter writer = new PrintWriter(diaTask.getExportFinalProTxt());

			StringBuilder titlesb = new StringBuilder();
			titlesb.append("Protein IDs").append("\t");
			titlesb.append("Peptide counts (all)").append("\t");
			titlesb.append("Peptide counts (razor)").append("\t");
			titlesb.append("Number of PSMs").append("\t");
			titlesb.append("protein_score_S").append("\t");
			titlesb.append("protein_p_emp").append("\t");
			titlesb.append("protein_q_value").append("\t");
			titlesb.append("Intensity").append("\t");
			for (int i = 0; i < diaTask.getExportQuanExpNames().length; i++) {
				titlesb.append("Intensity " + diaTask.getExportQuanExpNames()[i]).append("\t");
			}

			titlesb.append("Reverse").append("\t");
			titlesb.append("Potential contaminant").append("\t");
			titlesb.append("id").append("\t");
			titlesb.append("Peptide IDs").append("\t");
			titlesb.append("Peptide is razor");

			writer.println(titlesb);

			DecimalFormat df = FormatTool.getDFE4();
			DecimalFormat df4 = FormatTool.getDF4();
			for (int i = 0; i < proteins.length; i++) {

				StringBuilder sb = new StringBuilder();
				sb.append(proteins[i]).append("\t");
				sb.append(totalPepCountMap.get(proteins[i]).size()).append("\t");
				sb.append(razorPepCountMap.get(proteins[i]).size()).append("\t");
				sb.append(psmCountMap.get(proteins[i])).append("\t");

				double[] proScores = proteinScoreMap.get(proteins[i]);
				if (proScores != null) {
					sb.append(df.format(proScores[IDX_S])).append("\t");
					sb.append(df4.format(proScores[IDX_P_EMP])).append("\t");
					sb.append(df4.format(proScores[IDX_Q_VALUE])).append("\t");
				} else {
					sb.append("\t\t\t");
				}

				double totalIntensity = proteinTotalIntenMap.get(proteins[i]);
				sb.append(df.format(totalIntensity)).append("\t");

				double[] intensity = proteinIntensityMap.get(proteins[i]);
				for (int j = 0; j < intensity.length; j++) {
					sb.append(df.format(intensity[j])).append("\t");
				}
				sb.append("\t");
				sb.append("\t");
				sb.append(i).append("\t");

				ArrayList<Integer> pepIdList = proPepIdMap.get(proteins[i]);
				ArrayList<Integer> pepRazorList = proPepRazorIdMap.get(proteins[i]);

				for (int j = 0; j < pepIdList.size(); j++) {
					sb.append(pepIdList.get(j)).append(";");
				}
				sb.deleteCharAt(sb.length() - 1);
				sb.append("\t");

				for (int j = 0; j < pepRazorList.size(); j++) {
					if (pepRazorList.get(j) == 1) {
						sb.append("true;");
					} else {
						sb.append("false;");
					}
				}
				sb.deleteCharAt(sb.length() - 1);
				writer.println(sb);

				double scoreForMeta = (proScores != null) ? proScores[IDX_S] : 0.0;
				MetaProtein mp = new MetaProtein(i, 0, proteins[i], totalPepCountMap.get(proteins[i]).size(),
						psmCountMap.get(proteins[i]), scoreForMeta, intensity);
				proteinList.add(mp);
			}

			writer.close();

		} catch (IOException e) {
			// TODO Auto-generated catch block
			diaTask.getExportLogger().error(diaTask.getTaskName() + ": error in writing final protein result to " + diaTask.getExportFinalProTxt().getName(),
					e);
			System.err.println(fmt.format(new Date()) + "\t" + diaTask.getTaskName()
					+ ": error in writing final protein result to " + diaTask.getExportFinalProTxt().getName());
		}

		diaTask.getExportLogger().info(diaTask.getTaskName() + ": protein report has been exported to " + diaTask.getExportFinalProTxt().getName());
		System.out.println(fmt.format(new Date()) + "\t" + diaTask.getTaskName() + ": protein report has been exported to "
				+ diaTask.getExportFinalProTxt().getName());

		// --- all_proteins.tsv: all proteins, divided intensity per peptide ---
		String[] proteinsAll = proteinTotalIntenMapAll.keySet().toArray(new String[0]);
		Arrays.sort(proteinsAll, Comparator.<String>comparingDouble(p -> proteinScoreMap.containsKey(p) ? -proteinScoreMap.get(p)[IDX_S] : 0)
				.thenComparing((p1, p2) -> Double.compare(proteinTotalIntenMapAll.getOrDefault(p2, 0.0), proteinTotalIntenMapAll.getOrDefault(p1, 0.0))));
		HashMap<String, ArrayList<Integer>> proPepIdMapAll = new HashMap<>();
		HashMap<String, ArrayList<Integer>> proPepRazorIdMapAll = new HashMap<>();
		for (int id = 0; id < peptideList.size(); id++) {
			FlashLfqQuanPeptide pep = peptideList.get(id);
			String razorPro = pep.getRazorProtein();
			for (String pro : pep.getProteins()) {
				proPepIdMapAll.computeIfAbsent(pro, k -> new ArrayList<>()).add(id);
				proPepRazorIdMapAll.computeIfAbsent(pro, k -> new ArrayList<>()).add(razorPro.equals(pro) ? 1 : 0);
			}
		}
		try (PrintWriter writerAll = new PrintWriter(new File(diaTask.getExportResultFolderFile(), "all_proteins.tsv"))) {
			StringBuilder allHeader = new StringBuilder();
			allHeader.append("Protein IDs\tPeptide counts (all)\tPeptide counts (razor)\tNumber of PSMs\tprotein_score_S\tprotein_p_emp\tprotein_q_value\tIntensity\t");
			for (int h = 0; h < diaTask.getExportQuanExpNames().length; h++) allHeader.append("Intensity ").append(diaTask.getExportQuanExpNames()[h]).append("\t");
			allHeader.append("Reverse\tPotential contaminant\tid\tPeptide IDs\tPeptide is razor");
			writerAll.println(allHeader);
			DecimalFormat df = FormatTool.getDFE4();
			DecimalFormat df4 = FormatTool.getDF4();
			for (int i = 0; i < proteinsAll.length; i++) {
				String pro = proteinsAll[i];
				int pepAll = totalPepCountMapAll.containsKey(pro) ? totalPepCountMapAll.get(pro).size() : 0;
				int pepRazor = razorPepCountMap.containsKey(pro) ? razorPepCountMap.get(pro).size() : 0;
				StringBuilder sb = new StringBuilder();
				sb.append(pro).append("\t").append(pepAll).append("\t").append(pepRazor).append("\t");
				sb.append(psmCountMapAll.getOrDefault(pro, 0)).append("\t");
				double[] sc = proteinScoreMap.get(pro);
				if (sc != null) {
					sb.append(df.format(sc[IDX_S])).append("\t").append(df4.format(sc[IDX_P_EMP])).append("\t").append(df4.format(sc[IDX_Q_VALUE])).append("\t");
				} else sb.append("\t\t\t");
				sb.append(df.format(proteinTotalIntenMapAll.getOrDefault(pro, 0.0))).append("\t");
				double[] intensity = proteinIntensityMapAll.get(pro);
				if (intensity != null) for (double v : intensity) sb.append(df.format(v)).append("\t");
				sb.append("\t\t").append(i).append("\t");
				ArrayList<Integer> ids = proPepIdMapAll.get(pro);
				ArrayList<Integer> raz = proPepRazorIdMapAll.get(pro);
				if (ids != null) { for (Integer id : ids) sb.append(id).append(";"); if (!ids.isEmpty()) sb.setLength(sb.length() - 1); }
				sb.append("\t");
				if (raz != null) { for (Integer r : raz) sb.append(r == 1 ? "true;" : "false;"); if (!raz.isEmpty()) sb.setLength(sb.length() - 1); }
				writerAll.println(sb);
			}
		} catch (IOException e) {
			diaTask.getExportLogger().warn(diaTask.getTaskName() + ": could not write all_proteins.tsv", e);
		}

		diaTask.exportUpdateProgress(89, 100);

		diaTask.getExportLogger().info(diaTask.getTaskName() + ": functional annotation started");
		System.out.println(fmt.format(new Date()) + "\t" + diaTask.getTaskName() + ": functional annotation started");

		File funcFile = new File(diaTask.getExportResultFolderFile(), "functional_annotation");
		if (!funcFile.exists()) {
			funcFile.mkdir();
		}

		HashMap<String, Integer> proTaxIdMap = new HashMap<String, Integer>();
		MetaProtein[] metapros = proteinList.toArray(new MetaProtein[proteinList.size()]);
		MagFunctionalAnnotationService.Request funcRequest = new MagFunctionalAnnotationService.Request();
		funcRequest.resultDir = diaTask.getExportResultFolderFile();
		funcRequest.magDb = diaTask.getExportMagDb();
		funcRequest.microDbCatalogueId = diaTask.getExportMetaPar().getMicroDb();
		funcRequest.functionDefinitionDb = new File(diaTask.getExportMrp().getToolPath(ToolType.FUNCDB));
		funcRequest.quantitativeMode = diaTask.getExportMetaPar().getQuantitativeMode();
		funcRequest.intensityTitles = diaTask.getExportQuanExpNames();
		funcRequest.metaProteins = metapros;
		MagFunctionalAnnotationService.Result funcResult;
		try {
			funcResult = MagFunctionalAnnotationService.run(funcRequest);
		} catch (Exception e) {
			diaTask.getExportLogger().error(diaTask.getTaskName() + ": error in proteins functional annotations, task failed", e);
			System.err.println(fmt.format(new Date()) + "\t" + diaTask.getTaskName()
					+ ": error in proteins functional annotations, task failed");
			return false;
		}
		MetaProteinAnnoMag[] mpas = funcResult.annotatedProteins;
		for (int i = 0; i < mpas.length; i++) {
			String name = mpas[i].getPro().getName();
			if (NumberUtils.isCreatable(mpas[i].getTax_id())) {
				proTaxIdMap.put(name, Integer.parseInt(mpas[i].getTax_id()));
			} else {
				proTaxIdMap.put(name, 131567);
			}
		}

		File html = new File(funcFile, "functions.html");
		if (html.exists()) {
			try {
				File reportFunctionsHtml = new File(diaTask.getExportMetaPar().getReportDir(), "reports/functions.html");
				if (reportFunctionsHtml.getParentFile() != null && !reportFunctionsHtml.getParentFile().exists()) {
					reportFunctionsHtml.getParentFile().mkdirs();
				}
				Files.copy(html.toPath(), reportFunctionsHtml.toPath(), StandardCopyOption.REPLACE_EXISTING);
				diaTask.getExportLogger().info("Copied functions.html to reports directory: " + reportFunctionsHtml.getAbsolutePath());
			} catch (IOException e) {
				diaTask.getExportLogger().error(diaTask.getTaskName() + ": error copying functions.html to reports directory", e);
			}
		}

		// Copy functions.html to report directory for navigation access
		if (html.exists() && diaTask.getExportReportHtmlDir() != null) {
			File reportFunctionsHtml = new File(diaTask.getExportReportHtmlDir(), "functions.html");
			if (!reportFunctionsHtml.exists()) {
				try {
					java.nio.file.Files.copy(html.toPath(), reportFunctionsHtml.toPath(),
							java.nio.file.StandardCopyOption.REPLACE_EXISTING);
					diaTask.getExportLogger().info(diaTask.getTaskName() + ": copied functions.html to report directory");
					System.out.println(fmt.format(new Date()) + "\t" + diaTask.getTaskName()
							+ ": copied functions.html to report directory");
				} catch (IOException e) {
					diaTask.getExportLogger().error(diaTask.getTaskName() + ": error copying functions.html to report directory", e);
					System.err.println(fmt.format(new Date()) + "\t" + diaTask.getTaskName()
							+ ": error copying functions.html to report directory");
				}
			}
		}

		diaTask.exportUpdateProgress(91, 100);

		HashMap<String, HashSet<String>> uniquePepCountMapTax = new HashMap<String, HashSet<String>>();
		HashMap<String, HashSet<String>> totalPepWithILMap = new HashMap<String, HashSet<String>>();
		HashMap<String, HashSet<String>> razorPepWithILMap = new HashMap<String, HashSet<String>>();
		HashMap<String, HashSet<String>> uniquePepWithILMap = new HashMap<String, HashSet<String>>();
		DiaMagScoringService.fillDiaGenomeUniqueAndILMaps(peps, diaTask.getExportMagDb().getIdentifier(), razorPepCountMap,
				uniquePepCountMapTax, totalPepWithILMap, razorPepWithILMap, uniquePepWithILMap,
				genomeIntensityMap.keySet());

		HashMap<String, Double> genomeScoreForTaxonomy = new HashMap<>();
		for (String g : genomeIntensityMap.keySet()) {
			double[] st = genomeScoreMap.get(g);
			double p = st != null && st.length > DiaMagScoringService.IDX_P_EMP
					? Math.max(1e-300, st[DiaMagScoringService.IDX_P_EMP])
					: 1.0;
			genomeScoreForTaxonomy.put(g, Double.valueOf(p));
		}

		HashMap<String, double[]> prevDiaStats = diaTask.getLastDiaStyleGenomeStatsField();
		HashMap<String, double[]> statsCopy = new HashMap<>();
		for (java.util.Map.Entry<String, double[]> e : genomeScoreMap.entrySet()) {
			statsCopy.put(e.getKey(), e.getValue().clone());
		}
		diaTask.setLastDiaStyleGenomeStatsField(statsCopy);

		diaTask.exportUpdateProgress(94, 100);
		if (!diaTask.invokePerformTaxonomyAnalysisDdaFormat(genomeIntensityMap, genomeScoreForTaxonomy, psmCountMap,
				totalPepCountMap, razorPepCountMap, uniquePepCountMapTax, totalPepWithILMap, razorPepWithILMap,
				uniquePepWithILMap, diaTask.getExportExpNames(), diaTask.getTaskName())) {
			diaTask.setLastDiaStyleGenomeStatsField(prevDiaStats);
			return false;
		}
		diaTask.setLastDiaStyleGenomeStatsField(prevDiaStats);

		File taxFile = new File(diaTask.getExportResultFolderFile(), "taxonomy_analysis");
		File taxaReportFile = new File(taxFile, "Taxa_report.tsv");
		File genomeHtmlFile = new File(taxFile, "Genome.html");
		try {
			File reportGenomeHtml = new File(diaTask.getExportMetaPar().getReportDir(), "reports/Genome.html");
			if (reportGenomeHtml.getParentFile() != null && !reportGenomeHtml.getParentFile().exists()) {
				reportGenomeHtml.getParentFile().mkdirs();
			}
			if (genomeHtmlFile.exists()) {
				Files.copy(genomeHtmlFile.toPath(), reportGenomeHtml.toPath(), StandardCopyOption.REPLACE_EXISTING);
				diaTask.getExportLogger()
						.info("Copied Genome.html to reports directory: " + reportGenomeHtml.getAbsolutePath());
			}
		} catch (IOException e) {
			diaTask.getExportLogger().error(diaTask.getTaskName() + ": error copying Genome.html to reports directory",
					e);
		}
		if (genomeHtmlFile.exists() && diaTask.getExportReportHtmlDir() != null) {
			File reportGenomeHtmlLocal = new File(diaTask.getExportReportHtmlDir(), "Genome.html");
			if (!reportGenomeHtmlLocal.exists()) {
				try {
					java.nio.file.Files.copy(genomeHtmlFile.toPath(), reportGenomeHtmlLocal.toPath(),
							java.nio.file.StandardCopyOption.REPLACE_EXISTING);
					diaTask.getExportLogger().info(diaTask.getTaskName() + ": copied Genome.html to report directory");
				} catch (IOException e) {
					diaTask.getExportLogger().error(diaTask.getTaskName()
							+ ": error copying Genome.html to report directory", e);
				}
			}
		}

		diaTask.exportUpdateProgress(95, 100);

		diaTask.exportUpdateMessage(diaTask.getTaskName() + ": exporting report");

		diaTask.getExportLogger().info(diaTask.getTaskName() + ": exporting report started");
		System.out.println(fmt.format(new Date()) + "\t" + diaTask.getTaskName() + ": exporting report started");

		File reportDir = new File(diaTask.getExportResultFolderFile(), "report");
		if (!reportDir.exists()) {
			reportDir.mkdirs();
		}

		// Copy report templates using MetaPythonReportCopyTask (same as DDA workflow)
		MetaPythonReportCopyTask reportCopyTask = new MetaPythonReportCopyTask(diaTask.getExportResultFolderFile(), diaTask.getExportLogger());
		Thread copyThread = new Thread(reportCopyTask);
		copyThread.start();
		try {
			copyThread.join();
			if (!reportCopyTask.get()) {
				diaTask.getExportLogger().error(diaTask.getTaskName() + ": error in copying report template");
				System.err.println(
						fmt.format(new Date()) + "\t" + diaTask.getTaskName() + ": error in copying report template");
				return false;
			}
		} catch (Exception e) {
			diaTask.getExportLogger().error(diaTask.getTaskName() + ": error in copying report template", e);
			System.err.println(fmt.format(new Date()) + "\t" + diaTask.getTaskName() + ": error in copying report template: "
					+ e.getMessage());
			return false;
		}

		diaTask.getExportMetaPar().setReportDir(reportDir);

		diaTask.setExportReportHtmlDir(new File(reportDir + "\\reports"));

		if (!diaTask.getExportReportHtmlDir().exists()) {
			diaTask.getExportReportHtmlDir().mkdirs();
		}
		/*
		 * File jsOutput = new File(metaPar.getReportDir() + "\\reports\\tree.js"); if
		 * (!jsOutput.getParentFile().exists()) { jsOutput.getParentFile().mkdirs(); }
		 * 
		 * if (!jsOutput.exists() || jsOutput.length() == 0) {
		 * MetaTreeHandler.exportJS(taxRankMaps, taxIntensityMaps, cellularIntensity,
		 * expNames, jsOutput); }
		 * 
		 * // Generate bubble chart report from tree.js data (optional - don't fail if
		 * this doesn't work) File bubbleChartOutput = new File(metaPar.getReportDir() +
		 * "\\reports\\bubble_chart_report.html"); if (!bubbleChartOutput.exists() ||
		 * bubbleChartOutput.length() == 0) { try { generateBubbleChartReport(jsOutput,
		 * bubbleChartOutput); } catch (Exception e) { logger.warn(getTaskName() +
		 * ": Bubble chart generation failed, continuing with other reports: " +
		 * e.getMessage()); System.out.println(fmt.format(new Date()) + "\t" +
		 * getTaskName() +
		 * ": Bubble chart generation failed, continuing with other reports"); } }
		 */
		MetaPythonReportTask pyTask = new MetaPythonReportTask(diaTask.getExportMetaPar().getThreadCount(), diaTask.getExportMrp(), diaTask.getExportLogger());
		pyTask.setWorkflowType("DIA");
		pyTask.initializeManifestGenerator(new File(diaTask.getExportMetaPar().getReportDir(), "reports"));

		// Create initial manifest file to ensure index page can load
		ReportManifestGenerator manifestGen = new ReportManifestGenerator(new File(diaTask.getExportMetaPar().getReportDir(), "reports"));
		manifestGen.createInitialManifest();

		diaTask.setExportSummaryMetaFile(new File(diaTask.getExportResultFolderFile(), "metainfo.tsv"));
		try {
			diaTask.getExportMetadata().exportMetadata(diaTask.getExportSummaryMetaFile());
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		diaTask.exportSummary(pyTask);

		if (!diaTask.getExportFinalPepTxt().exists() || diaTask.getExportFinalPepTxt().length() == 0) {
			return false;
		} else {
			File report_peptides_summary = new File(diaTask.getExportReportHtmlDir(), "report_peptides_summary.html");
			if (!report_peptides_summary.exists() || report_peptides_summary.length() == 0) {
				if (diaTask.getExportSummaryMetaFile().exists()) {
					pyTask.addTask(MetaPythonReportTask.PEPTIDE, diaTask.getExportFinalPepTxt(), report_peptides_summary, diaTask.getExportSummaryMetaFile());
				} else {
					pyTask.addTask(MetaPythonReportTask.PEPTIDE, diaTask.getExportFinalPepTxt(), report_peptides_summary);
				}
			}
		}

		if (!diaTask.getExportFinalProTxt().exists() || diaTask.getExportFinalProTxt().length() == 0) {
			return false;
		} else {
			File report_proteinGroups_summary = new File(diaTask.getExportReportHtmlDir(), "report_proteinGroups_summary.html");
			if (!report_proteinGroups_summary.exists() || report_proteinGroups_summary.length() == 0) {
				if (diaTask.getExportSummaryMetaFile().exists()) {
					pyTask.addTask(MetaPythonReportTask.PROTEIN, diaTask.getExportFinalProTxt(), report_proteinGroups_summary,
							diaTask.getExportSummaryMetaFile());
				} else {
					pyTask.addTask(MetaPythonReportTask.PROTEIN, diaTask.getExportFinalProTxt(), report_proteinGroups_summary);
				}
			}
		}

		if (!taxaReportFile.exists() || taxaReportFile.length() == 0) {
			return false;
		} else {
			File report_taxonomy_summary = new File(diaTask.getExportReportHtmlDir(), "report_taxonomy_summary.html");
			if (!report_taxonomy_summary.exists() || report_taxonomy_summary.length() == 0) {

				if (diaTask.getExportSummaryMetaFile() != null && diaTask.getExportSummaryMetaFile().exists()) {
					pyTask.addTask(MetaPythonReportTask.TAXONOMY, taxaReportFile, report_taxonomy_summary,
							diaTask.getExportSummaryMetaFile());
				} else {
					pyTask.addTask(MetaPythonReportTask.TAXONOMY, taxaReportFile, report_taxonomy_summary);
				}
			}
		}

		if (!funcFile.exists() || funcFile.listFiles().length == 0) {
			return false;
		} else {
			File report_function_summary = new File(diaTask.getExportReportHtmlDir(), "report_function_summary.html");
			if (!report_function_summary.exists() || report_function_summary.length() == 0) {

				if (diaTask.getExportSummaryMetaFile() != null && diaTask.getExportSummaryMetaFile().exists()) {
					pyTask.addTask(MetaPythonReportTask.FUNCTION, funcFile, report_function_summary, diaTask.getExportSummaryMetaFile());
				} else {
					pyTask.addTask(MetaPythonReportTask.FUNCTION, funcFile, report_function_summary);
				}
			}
		}

		pyTask.run();

		try {
			boolean finish = pyTask.get();
			if (finish) {

				diaTask.getExportLogger().info(diaTask.getTaskName() + ": report generation successes");
				System.out.println(fmt.format(new Date()) + "\t" + diaTask.getTaskName() + ": report generation successes");

				diaTask.exportUpdateProgress(100, 100);
				return true;
			}

		} catch (InterruptedException | ExecutionException e) {
			// TODO Auto-generated catch block
			diaTask.getExportLogger().error(diaTask.getTaskName() + ": error in the report generation", e);
			System.out.println(fmt.format(new Date()) + "\t" + diaTask.getTaskName() + ": error in the report generation");
		}

		return true;
	}

	@Override
	protected boolean exportReportDownstreamFromNewScoring(java.util.HashMap<String, double[]> genomeScoreMap, java.util.HashMap<String, double[]> proteinScoreMap) {
		return runFullDownstreamFromNewScoring(this, genomeScoreMap, proteinScoreMap);
	}

}
