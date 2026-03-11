package uk.ac.quadram.metapilot.task.dia;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.math.BigInteger;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.sql.SQLException;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.concurrent.ExecutionException;

import org.apache.commons.lang3.math.NumberUtils;
import org.apache.logging.log4j.Logger;
import org.dom4j.DocumentException;

import javafx.application.Platform;
import javafx.scene.control.ProgressBar;

import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNParquetReader;
import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNPrecursor;
import uk.ac.quadram.metapilot.dbSearch.diann.DiannPrMatrixReader;
import uk.ac.quadram.metapilot.quant.flashLFQ.FlashLfqQuanPeptide;
import uk.ac.quadram.metapilot.core.prodb.FastaManager;
import uk.ac.quadram.metapilot.core.taxonomy.TaxonomyRanks;
import uk.ac.quadram.metapilot.core.tools.FormatTool;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters.ToolType;
import uk.ac.quadram.metapilot.task.dia.par.MetaParameterDiaMarker;
import uk.ac.quadram.metapilot.dbSearch.deepDetect.DeepDetectParameter;
import uk.ac.quadram.metapilot.dbSearch.deepDetect.DeepDetectTask;
import uk.ac.quadram.metapilot.core.MGYG.HapExtractor;
import uk.ac.quadram.metapilot.task.io.ReportManifestGenerator;
import uk.ac.quadram.metapilot.task.io.pep.MetaPeptide;
import uk.ac.quadram.metapilot.task.io.pro.MetaProtein;
import uk.ac.quadram.metapilot.task.io.pro.MetaProteinXMLReader2;
import uk.ac.quadram.metapilot.task.mag.HMGenomeHtmlWriter;
import uk.ac.quadram.metapilot.task.mag.MagFuncSearcher;
import uk.ac.quadram.metapilot.task.mag.MarkerGenomeSelector;
import uk.ac.quadram.metapilot.task.mag.MetaPilotMagTask;
import uk.ac.quadram.metapilot.task.mag.MetaProteinAnnoMag;
import uk.ac.quadram.metapilot.task.MetaPythonReportTask;
import uk.ac.quadram.metapilot.task.MetaPythonReportCopyTask;

/**
 * MetaDiaMarkerTask - A specialized DIA task for marker-based DIA analysis
 * 
 * This class extends MetaDiaTask to provide specific functionality for DIA
 * Marker workflow (Focus or Comprehensive genome selection), which uses
 * pre-generated focus libraries and MarkerGenomeSelector with user-selected mode.
 */
public class MetaDiaMarkerTask extends MetaDiaTask {

	// Logger is inherited from MetaDiaTask -> MetaPilotMagTask -> MetaPilotTask as
	// an instance field
	private static final String taskName = "DIA Marker analysis";

	// Combined file folder (set during iterativeSearch and used in iden)
	protected File combinedFile;

	/**
	 * Constructor for MetaDiaMarkerTask
	 * 
	 * @param markerParams The DIA Marker (Focus/Comprehensive) specific parameters
	 * @param mrp         Resource parameters
	 * @param bar1        Main progress bar
	 * @param bar2        Secondary progress bar
	 */
	public MetaDiaMarkerTask(MetaParameterDiaMarker markerParams, MetaResourceParameters mrp, ProgressBar bar1,
			ProgressBar bar2) {
		// Call parent constructor with the base DIA parameters
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

		logger.info("Initializing MetaDiaMarkerTask...");
		System.out.println(new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss").format(new Date()) + "\t" + getTaskName()
				+ ": initializing task");

		logger.info("DIA Focus workflow initialized with library search mode");
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
		// Use markerMode as the canonical selector; fall back to FOCUS.
		String markerMode = markerParams.getMarkerMode();
		String folderName;
		if ("COMPREHENSIVE".equalsIgnoreCase(markerMode) || "DEEP".equalsIgnoreCase(markerMode)) {
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
			secondBar.setProgress(0.0);
		} else {
			Platform.runLater(() -> secondBar.setProgress(0.0));
		}

		logger.info(getTaskName() + ": finished - analysis completed successfully");
		System.out.println(new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss").format(new Date()) + "\t" + getTaskName()
				+ ": finished - analysis completed successfully");
		return true;
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
		File genomeAllFile = magDb.getOriginalTaxaFile();
		if (genomeAllFile == null || !genomeAllFile.exists()) {
			logger.error("Cannot generate genomes-taxa.tsv: original taxa file not available");
			return;
		}
		try (BufferedReader reader = new BufferedReader(new FileReader(genomeAllFile))) {
			String headerLine = reader.readLine();
			if (headerLine == null) {
				return;
			}
			String[] title = headerLine.split("\t");
			int repId = -1;
			int lineageId = -1;
			for (int i = 0; i < title.length; i++) {
				if ("Species_rep".equals(title[i])) {
					repId = i;
				} else if ("Lineage".equals(title[i])) {
					lineageId = i;
				}
			}
			if (repId < 0 || lineageId < 0) {
				logger.error("Cannot generate genomes-taxa.tsv: Species_rep or Lineage column not found");
				return;
			}
			Map<String, String> map = new HashMap<>();
			String line;
			while ((line = reader.readLine()) != null) {
				String[] cs = line.split("\t");
				if (cs.length > Math.max(repId, lineageId) && !map.containsKey(cs[repId])) {
					map.put(cs[repId], cs[lineageId]);
				}
			}
			try (PrintWriter writer = new PrintWriter(genomeTaxaFile)) {
				writer.println("Genome\tLineage");
				for (String key : map.keySet()) {
					writer.println(key + "\t" + map.get(key));
				}
			}
			logger.info("Generated genomes-taxa.tsv: {}", genomeTaxaFile.getAbsolutePath());
		} catch (IOException e) {
			logger.error("Could not generate genomes-taxa.tsv", e);
		}
	}

	/**
	 * Ensures combined/ranked.tsv exists by running MarkerGenomeSelector on the combined precursors
	 * (same as DDA: ranking from combined result). No-op if combined/ranked.tsv already exists or
	 * precursors/magDb/Marker4.fasta/genomes-taxa are missing.
	 */
	static void ensureCombinedRankedFromPrecursors(MetaDiaTask diaTask, File resultFolder) {
		File combinedDir = new File(resultFolder, "combined");
		File rankedTsv = new File(combinedDir, "ranked.tsv");
		if (rankedTsv.exists()) return;
		if (diaTask.diaPp == null || diaTask.diaPp.length == 0) return;
		if (diaTask.getExportMagDb() == null) return;
		File currentFile = diaTask.getExportMagDb().getCurrentFile();
		if (currentFile == null || !currentFile.exists()) return;
		File markerFasta = new File(currentFile, "Marker4.fasta");
		File genomesTaxaTsv = new File(currentFile, "genomes-taxa.tsv");
		if (!markerFasta.exists() || !genomesTaxaTsv.exists()) return;
		String mode = MarkerGenomeSelector.KNEE_MODE_FOCUS;
		if (diaTask.getExportMetaPar() instanceof MetaParameterDiaMarker) {
			String markerMode = ((MetaParameterDiaMarker) diaTask.getExportMetaPar()).getMarkerMode();
			if ("COMPREHENSIVE".equalsIgnoreCase(markerMode) || "DEEP".equalsIgnoreCase(markerMode)) {
				mode = MarkerGenomeSelector.KNEE_MODE_COMPREHENSIVE;
			}
		}
		if (!combinedDir.exists()) combinedDir.mkdirs();
		File evidenceFile = new File(combinedDir, "evidence.tsv");
		File selectedFile = new File(combinedDir, "selected.tsv");
		try {
			MarkerGenomeSelector selector = new MarkerGenomeSelector();
			selector.setKneeSelectionMode(mode);
			selector.selectMarkerGenomesWithKnee(diaTask.diaPp, markerFasta, genomesTaxaTsv,
					evidenceFile.getAbsolutePath(), rankedTsv.getAbsolutePath(), selectedFile.getAbsolutePath(),
					500, true);
			diaTask.getExportLogger().info(diaTask.getTaskName() + ": generated combined/ranked.tsv from combined result for Selection_label");
		} catch (Exception e) {
			diaTask.getExportLogger().warn(diaTask.getTaskName() + ": could not generate combined/ranked.tsv: " + e.getMessage());
		}
	}

	/**
	 * Reads ranked.tsv/genome_ranked.tsv and fills export labels (core/mid/tail) for
	 * final Genome.tsv reporting.
	 *
	 * Export labeling rule (independent from selection stop strategy):
	 * 1) Rows before greedy starts (GreedyGain == 0): core.
	 * 2) Greedy starts at first row with GreedyGain != 0.
	 * 3) From greedy start, first focusKnee rows: core.
	 * 4) Next (comprehensiveKnee - focusKnee) rows: mid.
	 * 5) Remaining greedy rows: tail.
	 * 6) Stop at the first row where GreedyGain returns to 0; that row is not included.
	 */
	static void readCoreMidTailFromRankedTsv(File rankedFile,
			HashMap<String, String> genomeSelectionLabel, java.util.Set<String> genomeIntensityGenomes) {
		MetaPilotMagTask.RankedGenomeSelectionResult result =
				MetaPilotMagTask.buildExportSelectionFromRankedTsv(rankedFile, genomeIntensityGenomes);
		genomeSelectionLabel.putAll(result.genomeSelectionLabel);
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
			logger.error(getTaskName() + ": identification failed - result files not created");
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
		File markerLibFile = magDb.getMarker4LibFile();
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
		String markerMode = markerParams != null ? markerParams.getMarkerMode() : null;
		String mode;
		if ("COMPREHENSIVE".equalsIgnoreCase(markerMode) || "DEEP".equalsIgnoreCase(markerMode)) {
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
					diaLibSearchPar.setMbr(false);
					diaLibSearchPar.setDiaFiles(diaFiles);
					diaLibSearchPar.setqValue(0.01);
					diannPar.setThreads(threadCount);
					diannPar.setMatrices(false);
					diannPar.setRelaxed_prot_inf(true);
					diannPar.setReanalyse(true);
					diannPar.setExportDecoy(true);
					logger.info("Searching marker library against DIA files...");
					diaNNTask.addNoQuantLibTask(diannPar, diaLibSearchPar, tempParquetFile.getAbsolutePath(),
							markerLibFile.getAbsolutePath());
					diaNNTask.run(diaFiles.length * 12);
				}
				if (!tempParquetFile.exists() && diaFiles != null && diaFiles.length > 1) {
					List<String> damaged = runPerFileFallbackAndRetryCombinedNoQuant(tempParquetFile.getAbsolutePath(),
							new String[] { markerLibFile.getAbsolutePath() });
					if (!damaged.isEmpty()) {
						logger.warn(getTaskName() + ": The following raw files produced no results (may be damaged): {}", damaged);
						System.err.println(new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss").format(new Date()) + "\t" + getTaskName()
								+ ": The following raw files produced no results (may be damaged): " + String.join(", ", damaged));
					}
				}
				DiaNNPrecursor[] precursors;
				if (tempParquetFile.exists()) {
					DiaNNParquetReader parquetReader = new DiaNNParquetReader(
							tempParquetFile.getAbsolutePath(), diannPar.getCut());
					precursors = parquetReader.getDiaNNPrecursors(0.01);
				} else {
					// Combined marker.parquet missing after per-file fallback: use per_file_*.parquet directly
					logger.warn(getTaskName()
							+ ": Iterative marker search combined parquet missing; using per-file parquet results for genome selection");
					precursors = readPrecursorsFromPerFileParquets(markerFile, 0.01);
				}
				if (precursors == null || precursors.length == 0) {
					logger.warn(getTaskName() + ": No precursors found in iterative search, cannot select genomes");
					return false;
				}
				logger.info("Selecting genomes using MarkerGenomeSelector ({} mode)...", mode);
				MarkerGenomeSelector selector = new MarkerGenomeSelector();
				selector.setKneeSelectionMode(mode);
				File evidenceFile = new File(markerFile, "genome_evidence.tsv");
				File rankedFile = new File(markerFile, "genome_ranked.tsv");
				File selectedFile = new File(markerFile, "genome_selected.tsv");
				selector.selectMarkerGenomesWithKnee(precursors, markerFasta, genomesTaxaTsv,
						evidenceFile.getAbsolutePath(), rankedFile.getAbsolutePath(), selectedFile.getAbsolutePath(),
						300, true);
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
	 * Read DiaNN precursors from all non-empty per_file_*.parquet files in a directory.
	 * Used when the initial combined marker.parquet could not be created and per-file
	 * fallback was triggered.
	 *
	 * @param dir     directory containing per_file_*.parquet
	 * @param qValue  optional q-value cutoff; if null, use default in DiaNNParquetReader
	 * @return concatenated precursor array, or null if none
	 */
	private DiaNNPrecursor[] readPrecursorsFromPerFileParquets(File dir, Double qValue) {
		if (dir == null || !dir.isDirectory()) {
			return null;
		}
		String cut = diannPar != null && diannPar.getCut() != null ? diannPar.getCut() : "KR";
		List<DiaNNPrecursor> all = new ArrayList<>();
		for (int i = 0; i < (diaFiles != null ? diaFiles.length : Integer.MAX_VALUE); i++) {
			File pf = new File(dir, "per_file_" + i + ".parquet");
			if (!pf.exists()) {
				// per_file_* are written densely from 0..N-1; break at first missing index
				break;
			}
			if (pf.length() <= 0) {
				continue;
			}
			try {
				DiaNNParquetReader reader = new DiaNNParquetReader(pf.getAbsolutePath(), cut);
				DiaNNPrecursor[] prec = (qValue == null)
						? reader.getDiaNNPrecursors()
						: reader.getDiaNNPrecursors(qValue);
				if (prec != null) {
					for (DiaNNPrecursor p : prec) {
						if (p != null) {
							all.add(p);
						}
					}
				}
			} catch (Exception e) {
				logger.warn(getTaskName() + ": failed to read per-file parquet " + pf.getAbsolutePath(), e);
			}
		}
		return all.isEmpty() ? null : all.toArray(new DiaNNPrecursor[0]);
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
	 * Run the full downstream report for any MetaDiaTask. Used by MetaDiaSeedDirectTask.
	 * For Marker, runs the full implementation; for others (e.g. SeedDirect), returns true so that
	 * at least the new genome/protein score TSVs are considered success (full HTML report is Marker-only for now).
	 */
	protected static boolean runExportReportDownstream(MetaDiaTask task,
			java.util.HashMap<String, double[]> genomeScoreMap, java.util.HashMap<String, double[]> proteinScoreMap) {
		if (task instanceof MetaDiaMarkerTask)
			return ((MetaDiaMarkerTask) task).exportReportDownstreamFromNewScoring(genomeScoreMap, proteinScoreMap);
		if (task instanceof MetaDiaSeedDirectTask)
			return runFullDownstreamFromNewScoring(task, genomeScoreMap, proteinScoreMap);
		return true;
	}

	/**
	 * Full downstream report (peptide report, intensity maps, taxonomy, HTML) for any MetaDiaTask.
	 * Used by MetaDiaMarkerTask and MetaDiaSeedDirectTask.
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
		DiannPrMatrixReader prReader = new DiannPrMatrixReader(quanPepFile);
		FlashLfqQuanPeptide[] peps = prReader.getMetaPeptides(diaTask.diaPp);
		diaTask.setExportQuanExpNames(prReader.getIntensityTitle());

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
							if (genomeScoreMap.containsKey(genome)) {
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
			titlesb.append("Proteins").append("\t");
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

				String[] pros = peptide.getProteins();

				for (String pro : pros) {
					if (proteinIdMap.containsKey(pro)) {
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

			diaTask.setExportFinalProTxt(new File(diaTask.getExportResultFolderFile(), "final_proteins.tsv"));

			PrintWriter writer = new PrintWriter(diaTask.getExportFinalProTxt());

			StringBuilder titlesb = new StringBuilder();
			titlesb.append("Protein IDs").append("\t");
			titlesb.append("Peptide counts (all)").append("\t");
			titlesb.append("Peptide counts (razor)").append("\t");
			titlesb.append("Number of PSMs").append("\t");
			titlesb.append("S").append("\t");
			titlesb.append("p_emp").append("\t");
			titlesb.append("q_value").append("\t");
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
			allHeader.append("Protein IDs\tPeptide counts (all)\tPeptide counts (razor)\tNumber of PSMs\tS\tp_emp\tq_value\tIntensity\t");
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

		File funTsv = new File(funcFile, "functions.tsv");
		diaTask.setExportFinalProXmlFile(new File(funcFile, "functions.xml"));

		HashMap<String, Integer> proTaxIdMap = new HashMap<String, Integer>();
		MetaProtein[] metapros = proteinList.toArray(new MetaProtein[proteinList.size()]);
		MetaProteinAnnoMag[] mpas = diaTask.exportFuncAnnoSqlite(metapros, diaTask.getExportQuanExpNames());
		for (int i = 0; i < mpas.length; i++) {
			String name = mpas[i].getPro().getName();
			if (NumberUtils.isCreatable(mpas[i].getTax_id())) {
				proTaxIdMap.put(name, Integer.parseInt(mpas[i].getTax_id()));
			} else {
				proTaxIdMap.put(name, 131567);
			}
		}

		MetaProteinXMLReader2 funReader = null;
		try {
			funReader = new MetaProteinXMLReader2(diaTask.getExportFinalProXmlFile());
		} catch (DocumentException e) {
			// TODO Auto-generated catch block
			diaTask.getExportLogger().error(diaTask.getTaskName() + ": error in reading protein function information from "
					+ diaTask.getExportFinalProXmlFile().getName(), e);
			System.err.println(fmt.format(new Date()) + "\t" + diaTask.getTaskName()
					+ ": error in reading protein function information from " + diaTask.getExportFinalProXmlFile().getName());
			return false;
		}

		if (!funTsv.exists()) {
			System.out.println(fmt.format(new Date()) + "\t" + diaTask.getTaskName() + ": writing results to "
					+ funTsv.getName() + " started");
			try {
				funReader.exportTsv(funTsv);
				// Export separate functional category reports
				File goReportTsv = new File(funcFile, "Gene_Ontology_id_report.tsv");
				File ecReportTsv = new File(funcFile, "EC_id_report.tsv");
				File keggReportTsv = new File(funcFile, "KEGG_ko_report.tsv");
				File briteReportTsv = new File(funcFile, "BRITE_report.tsv");
				File cogReportTsv = new File(funcFile, "COG_accession_report.tsv");
				File nogReportTsv = new File(funcFile, "NOG_accession_report.tsv");
				File cogCatTsv = new File(funcFile, "COG_category_report.tsv");
				File nogCatTsv = new File(funcFile, "NOG_category_report.tsv");

				funReader.exportGOReport(goReportTsv);
				funReader.exportECReport(ecReportTsv);
				funReader.exportKeggReport(keggReportTsv);
				funReader.exportBriteReport(briteReportTsv);
				funReader.exportCOGReport(cogReportTsv);
				funReader.exportNOGReport(nogReportTsv);
				funReader.exportCOGCategoryReport(cogCatTsv);
				funReader.exportNOGCategoryReport(nogCatTsv);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				diaTask.getExportLogger().error(diaTask.getTaskName() + ": error in writing protein function information to " + funTsv.getName(),
						e);
				System.err.println(fmt.format(new Date()) + "\t" + diaTask.getTaskName()
						+ ": error in writing protein function information to " + funTsv.getName());
			}
			System.out.println(fmt.format(new Date()) + "\t" + diaTask.getTaskName() + ": writing results to "
					+ funTsv.getName() + " finished");
		}

		File html = new File(funcFile, "functions.html");
		if (!html.exists()) {
			System.out.println(fmt.format(new Date()) + "\t" + diaTask.getTaskName() + ": writing results to "
					+ html.getName() + " started");
			try {
				funReader.exportHtml(html);

				// Copy the generated functions.html to the reports directory
				File reportFunctionsHtml = new File(diaTask.getExportMetaPar().getReportDir(), "reports/functions.html");
				if (reportFunctionsHtml.getParentFile() != null && !reportFunctionsHtml.getParentFile().exists()) {
					reportFunctionsHtml.getParentFile().mkdirs();
				}
				Files.copy(html.toPath(), reportFunctionsHtml.toPath(), StandardCopyOption.REPLACE_EXISTING);
				diaTask.getExportLogger().info("Copied functions.html to reports directory: " + reportFunctionsHtml.getAbsolutePath());

			} catch (IOException e) {
				// TODO Auto-generated catch block
				diaTask.getExportLogger().error(diaTask.getTaskName() + ": error in writing protein function information to " + html.getName(), e);
				System.err.println(fmt.format(new Date()) + "\t" + diaTask.getTaskName()
						+ ": error in writing protein function information to " + html.getName());
			}
			System.out.println(fmt.format(new Date()) + "\t" + diaTask.getTaskName() + ": writing results to "
					+ html.getName() + " finished");
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

		// taxonomy analysis

		diaTask.getExportLogger().info(diaTask.getTaskName() + ": taxonomy analysis started");
		System.out.println(fmt.format(new Date()) + "\t" + diaTask.getTaskName() + ": taxonomy analysis started");

		File taxFile = new File(diaTask.getExportResultFolderFile(), "taxonomy_analysis");
		if (!taxFile.exists()) {
			taxFile.mkdir();
		}

		MagFuncSearcher searcher = null;
		try {
			searcher = new MagFuncSearcher(diaTask.getExportMetaPar().getMicroDb(), diaTask.getExportMrp().getToolPath(ToolType.FUNCDB), diaTask.getExportMagDb());

		} catch (NumberFormatException | SQLException | IOException e) {
			// TODO Auto-generated catch block
			diaTask.getExportLogger().error(diaTask.getTaskName() + ": error in genome analysis, task failed", e);
			System.err.println(
					fmt.format(new Date()) + "\t" + diaTask.getTaskName() + ": error in genome analysis, task failed");

			return false;
		}

		String[] genomes = genomeIntensityMap.keySet().toArray(new String[genomeIntensityMap.size()]);
		HashMap<String, String[]> genomeTaxaMap = null;
		try {
			genomeTaxaMap = searcher.matchGenome(genomes);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			diaTask.getExportLogger().error(diaTask.getTaskName() + ": Error in searching the genome2taxa database", e);
			System.err.println(
					fmt.format(new Date()) + "\t" + diaTask.getTaskName() + ": error in searching the genome2taxa database");
			return false;
		}

		File genomeHtmlFile = new File(taxFile, "Genome.html");
		try {
			HMGenomeHtmlWriter htmlWriter = new HMGenomeHtmlWriter(genomeHtmlFile);
			htmlWriter.write(genomeTaxaMap);
			htmlWriter.close();

			// Copy the generated Genome.html to the reports directory
			File reportGenomeHtml = new File(diaTask.getExportMetaPar().getReportDir(), "reports/Genome.html");
			if (reportGenomeHtml.getParentFile() != null && !reportGenomeHtml.getParentFile().exists()) {
				reportGenomeHtml.getParentFile().mkdirs();
			}
			Files.copy(genomeHtmlFile.toPath(), reportGenomeHtml.toPath(), StandardCopyOption.REPLACE_EXISTING);
			diaTask.getExportLogger().info("Copied Genome.html to reports directory: " + reportGenomeHtml.getAbsolutePath());
		} catch (IOException e) {
			// TODO Auto-generated catch block
			diaTask.getExportLogger().error("Error in writing genome information to " + genomeHtmlFile, e);
			System.err.println(
					fmt.format(new Date()) + "\t" + "Error in writing genome information to " + genomeHtmlFile);
			return false;
		}

		// Copy Genome.html to report directory for navigation access
		if (genomeHtmlFile.exists() && diaTask.getExportReportHtmlDir() != null) {
			File reportGenomeHtml = new File(diaTask.getExportReportHtmlDir(), "Genome.html");
			if (!reportGenomeHtml.exists()) {
				try {
					java.nio.file.Files.copy(genomeHtmlFile.toPath(), reportGenomeHtml.toPath(),
							java.nio.file.StandardCopyOption.REPLACE_EXISTING);
					diaTask.getExportLogger().info(diaTask.getTaskName() + ": copied Genome.html to report directory");
					System.out.println(fmt.format(new Date()) + "\t" + diaTask.getTaskName()
							+ ": copied Genome.html to report directory");
				} catch (IOException e) {
					diaTask.getExportLogger().error(diaTask.getTaskName() + ": error copying Genome.html to report directory", e);
					System.err.println(fmt.format(new Date()) + "\t" + diaTask.getTaskName()
							+ ": error copying Genome.html to report directory");
				}
			}
		}

		TaxonomyRanks[] mainRanks = TaxonomyRanks.getMainRanks();
		double[] cellularIntensity = new double[diaTask.getExportQuanExpNames().length];
		HashMap<String, String[]>[] taxRankMaps = new HashMap[mainRanks.length];
		HashMap<String, double[]>[] taxIntensityMaps = new HashMap[mainRanks.length];
		HashMap<String, Double>[] taxScoreMaps = new HashMap[mainRanks.length];
		for (int i = 0; i < taxRankMaps.length; i++) {
			taxRankMaps[i] = new HashMap<String, String[]>();
			taxIntensityMaps[i] = new HashMap<String, double[]>();
			taxScoreMaps[i] = new HashMap<String, Double>();
		}

		// Core/mid/tail classification from ranked.tsv (same as DDA: based on combined result).
		// Prefer combined/ranked.tsv (genome ranking from combined search); fallback to marker/ or seed/ ranked.
		HashMap<String, String> genomeSelectionLabel = new HashMap<>();
		File resultFolder = diaTask.getExportResultFolderFile();
		ensureCombinedRankedFromPrecursors(diaTask, resultFolder);
		File rankedFile = new File(resultFolder, "combined/ranked.tsv");
		if (!rankedFile.exists()) {
			rankedFile = new File(resultFolder, "marker/genome_ranked.tsv");
		}
		if (!rankedFile.exists()) {
			rankedFile = new File(resultFolder, "seed/genome_ranked.tsv");
		}
		if (rankedFile.exists()) {
			// Export labels for genomes that actually have intensity
			readCoreMidTailFromRankedTsv(rankedFile, genomeSelectionLabel, genomeIntensityMap.keySet());

			// Global statistics from ranked.tsv (no intensity filter)
			MetaPilotMagTask.RankedGenomeSelectionResult rankedAll =
					MetaPilotMagTask.buildExportSelectionFromRankedTsv(rankedFile, null);

			int contributingAll = rankedAll.orderedGenomes.size();
			int coreAll = 0, midAll = 0, tailAll = 0;
			for (String g : rankedAll.orderedGenomes) {
				String lbl = rankedAll.genomeSelectionLabel.get(g);
				if ("core".equals(lbl)) coreAll++;
				else if ("mid".equals(lbl)) midAll++;
				else if ("tail".equals(lbl)) tailAll++;
			}

			int labeledWithIntensity = genomeSelectionLabel.size();
			int coreWithIntensity = 0, midWithIntensity = 0, tailWithIntensity = 0;
			for (String lbl : genomeSelectionLabel.values()) {
				if ("core".equals(lbl)) coreWithIntensity++;
				else if ("mid".equals(lbl)) midWithIntensity++;
				else tailWithIntensity++;
			}

			java.util.HashSet<String> contributingSet = new java.util.HashSet<>(rankedAll.orderedGenomes);
			java.util.HashSet<String> intensitySet = new java.util.HashSet<>(genomeIntensityMap.keySet());
			java.util.HashSet<String> inWindowNoIntensity = new java.util.HashSet<>(contributingSet);
			inWindowNoIntensity.removeAll(intensitySet);
			java.util.HashSet<String> intensityBeyondWindow = new java.util.HashSet<>(intensitySet);
			intensityBeyondWindow.removeAll(contributingSet);

			// Detailed ranked.tsv summary
			diaTask.getExportLogger().info(
					"{}: Genome ranked.tsv summary – total_ranked={}, privilege_rows={}, greedy_rows={}, focusKnee={}, comprehensiveKnee={}, core_all={}, mid_all={}, tail_all={}",
					diaTask.getTaskName(),
					Integer.valueOf(rankedAll.totalRanked),
					Integer.valueOf(rankedAll.privilegeCount),
					Integer.valueOf(rankedAll.greedyCount),
					Integer.valueOf(rankedAll.focusKnee),
					Integer.valueOf(rankedAll.comprehensiveKnee),
					Integer.valueOf(coreAll),
					Integer.valueOf(midAll),
					Integer.valueOf(tailAll));

			// Export window vs intensity intersection summary
			diaTask.getExportLogger().info(
					"{}: Genome export summary – ranked_window={}, labeled_with_intensity={}, core_with_intensity={}, mid_with_intensity={}, tail_with_intensity={}, window_no_intensity={}, intensity_beyond_window={}",
					diaTask.getTaskName(),
					Integer.valueOf(contributingAll),
					Integer.valueOf(labeledWithIntensity),
					Integer.valueOf(coreWithIntensity),
					Integer.valueOf(midWithIntensity),
					Integer.valueOf(tailWithIntensity),
					Integer.valueOf(inWindowNoIntensity.size()),
					Integer.valueOf(intensityBeyondWindow.size()));
		}
		if (rankedFile.exists() && !genomeSelectionLabel.isEmpty()) {
			ArrayList<String> filteredGenomes = new ArrayList<>();
			for (String genome : genomes) {
				if (genomeSelectionLabel.containsKey(genome)) {
					filteredGenomes.add(genome);
				}
			}
			genomes = filteredGenomes.toArray(new String[0]);
		}

		// Sort genomes for Genome.tsv: core → mid → tail, then by q-value ascending.
		Arrays.sort(genomes, (g1, g2) -> {
			String label1 = genomeSelectionLabel.getOrDefault(g1, "tail");
			String label2 = genomeSelectionLabel.getOrDefault(g2, "tail");
			int group1 = "core".equals(label1) ? 0 : ("mid".equals(label1) ? 1 : 2);
			int group2 = "core".equals(label2) ? 0 : ("mid".equals(label2) ? 1 : 2);
			if (group1 != group2) {
				return Integer.compare(group1, group2);
			}
			double[] s1 = genomeScoreMap.get(g1);
			double[] s2 = genomeScoreMap.get(g2);
			double q1 = (s1 != null && s1.length > IDX_Q_VALUE) ? s1[IDX_Q_VALUE] : Double.POSITIVE_INFINITY;
			double q2 = (s2 != null && s2.length > IDX_Q_VALUE) ? s2[IDX_Q_VALUE] : Double.POSITIVE_INFINITY;
			int cmpQ = Double.compare(q1, q2);
			if (cmpQ != 0) {
				return cmpQ;
			}
			// Tie-breaker: higher S first, then lexicographic genome ID
			double sScore1 = (s1 != null && s1.length > IDX_S) ? s1[IDX_S] : 0.0;
			double sScore2 = (s2 != null && s2.length > IDX_S) ? s2[IDX_S] : 0.0;
			int cmpS = Double.compare(sScore2, sScore1);
			if (cmpS != 0) {
				return cmpS;
			}
			return g1.compareTo(g2);
		});

		File genomeFile = new File(taxFile, "Genome.tsv");
		PrintWriter genomeWriter = null;
		try {
			genomeWriter = new PrintWriter(genomeFile);
			StringBuilder genomeTitlesb = new StringBuilder();
			genomeTitlesb.append("Genome\t");
			genomeTitlesb.append("S\t");
			genomeTitlesb.append("p_emp\t");
			genomeTitlesb.append("q_value\t");
			genomeTitlesb.append("PSM count\t");
			genomeTitlesb.append("Peptide count\t");
			genomeTitlesb.append("Razor peptide count\t");
			genomeTitlesb.append("Selection_label\t");
			genomeTitlesb.append("GTDB_Name\t");
			genomeTitlesb.append("GTDB_Rank\t");

			for (int i = 0; i < mainRanks.length; i++) {
				genomeTitlesb.append("GTDB_").append(mainRanks[i].getName()).append("\t");
			}

			genomeTitlesb.append("NCBI_Name\t");
			genomeTitlesb.append("NCBI_Rank");

			for (int i = 0; i < mainRanks.length; i++) {
				genomeTitlesb.append("\t").append("NCBI_").append(mainRanks[i].getName());
			}

			for (String exp : diaTask.getExportExpNames()) {
				if (exp.startsWith("Intensity ")) {
					genomeTitlesb.append("\t").append(exp);
				} else {
					genomeTitlesb.append("\tIntensity ").append(exp);
				}
			}

			genomeWriter.println(genomeTitlesb);

			diaTask.exportUpdateProgress(94, 100);

			for (int i = 0; i < genomes.length; i++) {
				String[] taxaContent = genomeTaxaMap.get(genomes[i]);
				if (genomeIntensityMap.containsKey(genomes[i])) {
					double[] genomeIntensity = genomeIntensityMap.get(genomes[i]);
					for (int j = 0; j < genomeIntensity.length; j++) {
						cellularIntensity[j] += genomeIntensity[j];
					}

					StringBuilder genomeTaxaSb = new StringBuilder();
					genomeTaxaSb.append(genomes[i]).append("\t");

					double[] gScores = genomeScoreMap.get(genomes[i]);
					genomeTaxaSb.append(FormatTool.getDFE4().format(gScores[IDX_S])).append("\t");
					genomeTaxaSb.append(FormatTool.getDF4().format(gScores[IDX_P_EMP])).append("\t");
					genomeTaxaSb.append(FormatTool.getDF4().format(gScores[IDX_Q_VALUE])).append("\t");

					if (psmCountMap.containsKey(genomes[i])) {
						genomeTaxaSb.append(psmCountMap.get(genomes[i])).append("\t");
					} else {
						genomeTaxaSb.append("0").append("\t");
					}
					if (totalPepCountMap.containsKey(genomes[i])) {
						genomeTaxaSb.append(totalPepCountMap.get(genomes[i]).size()).append("\t");
					} else {
						genomeTaxaSb.append("0").append("\t");
					}
					if (razorPepCountMap.containsKey(genomes[i])) {
						genomeTaxaSb.append(razorPepCountMap.get(genomes[i]).size());
					} else {
						genomeTaxaSb.append("0");
					}
					genomeTaxaSb.append("\t").append(genomeSelectionLabel.getOrDefault(genomes[i], "tail"));

					for (int j = 1; j < taxaContent.length; j++) {
						genomeTaxaSb.append("\t").append(taxaContent[j]);
					}

					for (int j = 0; j < genomeIntensity.length; j++) {
						genomeTaxaSb.append("\t").append(genomeIntensity[j]);
					}
					genomeWriter.println(genomeTaxaSb);

					for (int j = 0; j < taxRankMaps.length; j++) {
						String taxon = taxaContent[j + 3];
						if (taxon == null || taxon.length() == 0) {
							continue;
						}
						if (!taxRankMaps[j].containsKey(taxon)) {
							String[] lineage = new String[j + 1];
							System.arraycopy(taxaContent, 3, lineage, 0, lineage.length);
							taxRankMaps[j].put(taxon, lineage);
						}

						double[] taxIntensity;
						if (taxIntensityMaps[j].containsKey(taxon)) {
							taxIntensity = taxIntensityMaps[j].get(taxon);
						} else {
							taxIntensity = new double[genomeIntensity.length];
						}

						for (int k = 0; k < taxIntensity.length; k++) {
							taxIntensity[k] += genomeIntensity[k];
						}
						taxIntensityMaps[j].put(taxon, taxIntensity);

						if (taxScoreMaps[j].containsKey(taxon)) {
							taxScoreMaps[j].put(taxon, taxScoreMaps[j].get(taxon) * gScores[IDX_S]);
						} else {
							taxScoreMaps[j].put(taxon, gScores[IDX_S]);
						}
					}
				}
			}
			genomeWriter.close();

		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			diaTask.getExportLogger().error("Error in writing genome information to " + genomeFile, e);
			System.err
					.println(fmt.format(new Date()) + "\t" + "Error in writing genome information to " + genomeFile);
			return false;
		}

		File taxaFile = new File(taxFile, "Taxa.tsv");
		File taxaReportFile = new File(taxFile, "Taxa_report.tsv");
		PrintWriter taxaWriter = null;
		PrintWriter taxaReportWriter = null;
		try {
			taxaWriter = new PrintWriter(taxaFile);
			taxaReportWriter = new PrintWriter(taxaReportFile);

			StringBuilder taxaTitlesb1 = new StringBuilder();
			taxaTitlesb1.append("Name\t");
			taxaTitlesb1.append("Score\t");
			taxaTitlesb1.append("Rank");

			for (int i = 0; i < mainRanks.length; i++) {
				taxaTitlesb1.append("\t").append(mainRanks[i].getName());
			}
			for (String exp : diaTask.getExportExpNames()) {
				if (exp.startsWith("Intensity ")) {
					taxaTitlesb1.append("\t").append(exp);
				} else {
					taxaTitlesb1.append("\tIntensity ").append(exp);
				}
			}
			StringBuilder taxaTitlesb2 = new StringBuilder();
			taxaTitlesb2.append("Name\t");
			taxaTitlesb2.append("Rank");

			for (int i = 0; i < mainRanks.length; i++) {
				taxaTitlesb2.append("\t").append(mainRanks[i].getName());
			}
			for (String exp : diaTask.getExportExpNames()) {
				if (exp.startsWith("Intensity ")) {
					taxaTitlesb2.append("\t").append(exp);
				} else {
					taxaTitlesb2.append("\tIntensity ").append(exp);
				}
			}
			taxaWriter.println(taxaTitlesb1);
			taxaReportWriter.println(taxaTitlesb2);

			for (int i = 0; i < taxRankMaps.length; i++) {
				for (String taxon : taxRankMaps[i].keySet()) {
					String[] taxaContent = taxRankMaps[i].get(taxon);
					double[] intensity = taxIntensityMaps[i].get(taxon);

					StringBuilder taxaSb = new StringBuilder();
					taxaSb.append(taxon).append("\t");

					double taxonScore = taxScoreMaps[i].get(taxon);
					if (taxonScore == 0) {
						taxaSb.append("1000").append("\t");
					} else {
						taxaSb.append(FormatTool.getDF4().format(-Math.log10(taxonScore))).append("\t");
					}

					taxaSb.append(mainRanks[i].getName());

					for (int j = 0; j < mainRanks.length; j++) {
						if (j < taxaContent.length) {
							taxaSb.append("\t").append(taxaContent[j]);
						} else {
							taxaSb.append("\t");
						}
					}

					for (int j = 0; j < intensity.length; j++) {
						taxaSb.append("\t").append(intensity[j]);
					}

					taxaWriter.println(taxaSb);
				}
			}

			taxaWriter.close();

			for (int i = 0; i < taxRankMaps.length; i++) {
				for (String taxon : taxRankMaps[i].keySet()) {
					String[] taxaContent = taxRankMaps[i].get(taxon);
					double[] intensity = taxIntensityMaps[i].get(taxon);

					StringBuilder taxaSb = new StringBuilder();
					taxaSb.append(taxon).append("\t");
					taxaSb.append(mainRanks[i].getName());

					for (int j = 0; j < mainRanks.length; j++) {
						if (j < taxaContent.length) {
							taxaSb.append("\t").append(taxaContent[j]);
						} else {
							taxaSb.append("\t");
						}
					}
					for (int j = 0; j < intensity.length; j++) {
						taxaSb.append("\t").append(intensity[j]);
					}

					taxaReportWriter.println(taxaSb);
				}
			}

			taxaReportWriter.close();

		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			diaTask.getExportLogger().error("Error in writing taxa information to " + taxaFile, e);
			System.err.println(fmt.format(new Date()) + "\t" + "Error in writing taxa information to " + taxaFile);
		}

		diaTask.getExportLogger().info(diaTask.getTaskName() + ": taxonomy analysis finished");
		System.out.println(fmt.format(new Date()) + "\t" + diaTask.getTaskName() + ": taxonomy analysis finished");

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

	/**
	 * Static test method: genome scoring via pooled decoy null.
	 * Uses UNIQUE peptides only (m=1); shared peptides are skipped.
	 * Combined PEP+intensity scoring: S = sum top-k (w * e * f(I1,I2)) with
	 * f(I1,I2) = 0.5 + 0.25*norm(log2(1+I1)) + 0.25*norm(log2(1+I2)), I1=Precursor.Quantity, I2=Ms1.Area.
	 * Reads parquet, collapses PSM→peptide, builds genome scores, computes empirical p-values
	 * by size-matched decoy resampling, BH q-values, exports genome rank and scores, and runinfo.
	 *
	 * @param parquetFile        parquet file (e.g. combined.parquet)
	 * @param uniquePeptidesOnly if true, use only genome-unique peptides (m=1); if false, use all peptides
	 *                           and assign shared peptides to each genome with weight w = 1/m
	 */
	/** Overload: unique peptides only (same as testReadPrecursorsFromParquet(parquetFile, true)). */
	public static void testReadPrecursorsFromParquet(File parquetFile) {
		testReadPrecursorsFromParquet(parquetFile, true);
	}

	public static void testReadPrecursorsFromParquet(File parquetFile, boolean uniquePeptidesOnly) {
		if (parquetFile == null || !parquetFile.exists()) {
			return;
		}
		try {
			String cut = "KR";
			DiaNNParquetReader reader = new DiaNNParquetReader(parquetFile.getAbsolutePath(), cut);
			DiaNNPrecursor[] precursors = reader.getDiaNNPrecursors();
			if (precursors == null) {
				return;
			}
			int nPrecursorsTotal = precursors.length;

			// --- 1. Collapse PSM → peptide ---
			// key: (seqString, charge, isDecoy); value: best CScore (max), best PEP (min), max Precursor.Quantity, max Ms1.Area
			HashMap<String, PeptideEvidence> peptideMap = new HashMap<>();
			for (DiaNNPrecursor pp : precursors) {
				String seq = pp.getSeqString();
				if (seq == null || seq.isEmpty()) {
					continue;
				}
				int charge = pp.getCharge();
				boolean isDecoy = pp.getIsDecoy();
				String key = seq + "|" + charge + "|" + isDecoy;
				double preQuan = pp.getPreQuan();
				double ms1Area = pp.getMs1Area();
				double cscore = pp.getCscore();
				double pep = pp.getPEP();

				PeptideEvidence ev = peptideMap.get(key);
				if (ev == null) {
					ev = new PeptideEvidence();
					ev.seq = seq;
					ev.bestCScore = cscore;
					ev.bestPEP = pep;
					ev.bestPreQuan = preQuan;
					ev.bestMs1Area = ms1Area;
					ev.isDecoy = isDecoy;
					ev.proteins = new HashSet<>();
					peptideMap.put(key, ev);
				} else {
					if (cscore > ev.bestCScore) {
						ev.bestCScore = cscore;
					}
					if (pep < ev.bestPEP) {
						ev.bestPEP = pep;
					}
					ev.bestPreQuan = Math.max(ev.bestPreQuan, preQuan);
					ev.bestMs1Area = Math.max(ev.bestMs1Area, ms1Area);
				}
				String[] proGroups = pp.getProGroups();
				if (proGroups != null) {
					for (String pro : proGroups) {
						if (pro != null && !pro.isEmpty()) {
							ev.proteins.add(pro);
						}
					}
				}
			}

			// --- 2a. Compute intensity factor: 0.5 + 0.25*norm(log2(1+I1)) + 0.25*norm(log2(1+I2)), I1=Precursor.Quantity, I2=Ms1.Area ---
			double cap = 6.0;
			double minLog1 = Double.MAX_VALUE, maxLog1 = -Double.MAX_VALUE;
			double minLog2 = Double.MAX_VALUE, maxLog2 = -Double.MAX_VALUE;
			for (PeptideEvidence ev : peptideMap.values()) {
				double logI1 = Math.log(1.0 + Math.max(0, ev.bestPreQuan)) / Math.log(2.0);
				double logI2 = Math.log(1.0 + Math.max(0, ev.bestMs1Area)) / Math.log(2.0);
				minLog1 = Math.min(minLog1, logI1);
				maxLog1 = Math.max(maxLog1, logI1);
				minLog2 = Math.min(minLog2, logI2);
				maxLog2 = Math.max(maxLog2, logI2);
			}
			double logRange1 = Math.max(1e-10, maxLog1 - minLog1);
			double logRange2 = Math.max(1e-10, maxLog2 - minLog2);

			// --- 2b. peptide → genome set, build target peptide list and decoy pool ---
			// Use only decoys with PEP < this for the null pool (avoids pool full of PEP=1, e=0)
			final double decoyPepMax = 0.99;
			ArrayList<Double> decoyEPool = new ArrayList<>();
			HashMap<String, GenomePeptideEntry> genomePeptides = new HashMap<>();
			int nTargetPeptides = 0;
			int nDecoyPeptides = 0;
			int nDecoyInPool = 0;  // decoys with PEP < decoyPepMax used in null
			int nUniquePeptides = 0;
			int nSharedPeptides = 0;
			HashMap<Integer, Integer> mDistribution = new HashMap<>();

			for (Map.Entry<String, PeptideEvidence> ent : peptideMap.entrySet()) {
				PeptideEvidence ev = ent.getValue();
				String seq = ev.seq;
				// Confidence: CScore if present (higher=better), else PEP (lower=better -> e = -log10(PEP))
				double e;
				if (ev.bestCScore > 0) {
					e = Math.min(cap, Math.max(0, ev.bestCScore));
				} else {
					e = Math.min(cap, Math.max(0, -Math.log10(Math.max(ev.bestPEP, 1e-10))));
				}
				double logI1 = Math.log(1.0 + Math.max(0, ev.bestPreQuan)) / Math.log(2.0);
				double logI2 = Math.log(1.0 + Math.max(0, ev.bestMs1Area)) / Math.log(2.0);
				double norm1 = (logI1 - minLog1) / logRange1;
				double norm2 = (logI2 - minLog2) / logRange2;
				double intensityFactor = 0.5 + 0.25 * norm1 + 0.25 * norm2;

				if (ev.isDecoy) {
					nDecoyPeptides++;
					// Only use decoys with PEP < threshold so the null pool is informative (not all e≈0)
					if (ev.bestPEP < decoyPepMax) {
						decoyEPool.add(e * intensityFactor);
						nDecoyInPool++;
					}
					continue;
				}

				// target: get genomes from proteins (skip decoy proteins)
				HashSet<String> genomes = new HashSet<>();
				for (String pro : ev.proteins) {
					if (pro.startsWith(FastaManager.shuffle)) {
						continue;
					}
					String[] parts = pro.split("_");
					if (parts.length > 0 && !parts[0].isEmpty()) {
						genomes.add(parts[0]);
					}
				}
				if (genomes.isEmpty()) {
					continue;
				}
				int m = genomes.size();
				if (m != 1) {
					nSharedPeptides++;
					mDistribution.merge(m, 1, Integer::sum);
					if (uniquePeptidesOnly) {
						continue;
					}
					// Use all peptides: add shared peptide to each genome with weight w = 1/m
					double w = 1.0 / m;
					double wei = w * e * intensityFactor;
					for (String g : genomes) {
						WeightedPeptideScore wps = new WeightedPeptideScore();
						wps.seq = seq;
						wps.e = e;
						wps.w = w;
						wps.we = wei;
						wps.intensityFactor = intensityFactor;
						GenomePeptideEntry ge = genomePeptides.get(g);
						if (ge == null) {
							ge = new GenomePeptideEntry();
							genomePeptides.put(g, ge);
						}
						ge.weightedScores.add(wps);
					}
					nTargetPeptides++;
					continue;
				}
				// Unique peptide (m=1)
				String g = genomes.iterator().next();
				double w = 1.0;
				double wei = w * e * intensityFactor;
				WeightedPeptideScore wps = new WeightedPeptideScore();
				wps.seq = seq;
				wps.e = e;
				wps.w = w;
				wps.we = wei;
				wps.intensityFactor = intensityFactor;

				GenomePeptideEntry ge = genomePeptides.get(g);
				if (ge == null) {
					ge = new GenomePeptideEntry();
					genomePeptides.put(g, ge);
				}
				ge.weightedScores.add(wps);

				nTargetPeptides++;
				nUniquePeptides++;
				mDistribution.merge(1, 1, Integer::sum);
			}

			// --- 3. Genome score S = sum of top-k (w*e*intensityFactor) ---
			int k = 5;
			Random rnd = new Random(42);
			int B = 10000;

			HashMap<String, Double> genomeScoreMap = new HashMap<>();
			HashMap<String, Integer> genomeNMap = new HashMap<>();
			HashMap<String, Integer> genomeTotalPeptidesMap = new HashMap<>();
			HashMap<String, List<WeightedPeptideScore>> genomeTopKMap = new HashMap<>();

			for (String g : genomePeptides.keySet()) {
				List<WeightedPeptideScore> scores = genomePeptides.get(g).weightedScores;
				genomeTotalPeptidesMap.put(g, scores.size());
				scores.sort((a, b) -> Double.compare(b.we, a.we));
				int n = Math.min(k, scores.size());
				double S = 0;
				List<WeightedPeptideScore> topK = new ArrayList<>();
				for (int i = 0; i < n; i++) {
					WeightedPeptideScore wps = scores.get(i);
					S += wps.we;
					topK.add(wps);
				}
				genomeScoreMap.put(g, S);
				genomeNMap.put(g, n);
				genomeTopKMap.put(g, topK);
			}

			// --- 4. Pooled decoy null: sample M, take top-k, sum ---
			if (decoyEPool.isEmpty()) {
				decoyEPool.add(0.0);
			}
			double[] decoyArr = new double[decoyEPool.size()];
			for (int i = 0; i < decoyArr.length; i++) {
				decoyArr[i] = decoyEPool.get(i);
			}
			int decoyLen = decoyArr.length;

			int maxM = 0;
			for (Integer m : genomeTotalPeptidesMap.values()) {
				maxM = Math.max(maxM, m);
			}
			double[] tmp = new double[maxM];

			HashMap<String, Double> genomePempMap = new HashMap<>();
			for (String g : genomeScoreMap.keySet()) {
				double Sg = genomeScoreMap.get(g);
				int M = genomeTotalPeptidesMap.get(g);
				int ng = genomeNMap.get(g);
				int count = 0;
				for (int b = 0; b < B; b++) {
					for (int i = 0; i < M; i++) {
						tmp[i] = decoyArr[rnd.nextInt(decoyLen)];
					}
					Arrays.sort(tmp, 0, M);
					double Sb = 0;
					for (int i = M - 1; i >= M - ng && i >= 0; i--) {
						Sb += tmp[i];
					}
					if (Sb >= Sg) {
						count++;
					}
				}
				double pemp = (1.0 + count) / (1.0 + B);
				genomePempMap.put(g, pemp);
			}

			// --- 5. BH correction for q-values ---
			String[] genomes = genomePempMap.keySet().toArray(new String[0]);
			Arrays.sort(genomes, (a, b) -> Double.compare(genomeScoreMap.get(a), genomeScoreMap.get(b)));
			List<String> sortedByScore = new ArrayList<>(Arrays.asList(genomes));
			Collections.reverse(sortedByScore);

			int N = sortedByScore.size();
			double[] pempArr = new double[N];
			for (int i = 0; i < N; i++) {
				pempArr[i] = genomePempMap.get(sortedByScore.get(i));
			}
			double[] qvals = bhCorrection(pempArr);

			// --- 6. Export genome rank and scores ---
			File outFile = new File(parquetFile.getParentFile(), "genome_rank_scores.tsv");
			try (PrintWriter pw = new PrintWriter(outFile)) {
				pw.println("Rank\tGenome\tS\tn_peptides\tp_emp\tq_value");
				for (int i = 0; i < N; i++) {
					String g = sortedByScore.get(i);
					pw.println((i + 1) + "\t" + g + "\t" + String.format("%.4f", genomeScoreMap.get(g)) + "\t"
							+ genomeNMap.get(g) + "\t" + String.format("%.6f", pempArr[i]) + "\t"
							+ String.format("%.6f", qvals[i]));
				}
			}

			// --- 6b. Filter precursors by p_emp=1 genomes: keep only precursors with at least one protein not from a genome with p_emp==1 ---
			HashSet<String> genomesPempOne = new HashSet<>();
			for (String g : genomePempMap.keySet()) {
				if (genomePempMap.get(g) >= 1.0 - 1e-9) {
					genomesPempOne.add(g);
				}
			}
			int nPrecursorsQ01 = 0;
			int nPrecursorsRetained = 0;
			HashSet<String> uniqueSeqQ01 = new HashSet<>();
			HashSet<String> uniqueSeqRetained = new HashSet<>();
			for (DiaNNPrecursor pp : precursors) {
				if (pp.getGlobalQvalue() >= 0.01) {
					continue;
				}
				nPrecursorsQ01++;
				String seq = pp.getSeqString();
				if (seq != null) {
					uniqueSeqQ01.add(seq);
				}
				String[] proGroups = pp.getProGroups();
				if (proGroups == null || proGroups.length == 0) {
					continue;
				}
				int kept = 0;
				for (String pro : proGroups) {
					if (pro == null || pro.isEmpty()) {
						continue;
					}
					if (pro.startsWith(FastaManager.shuffle)) {
						continue; // decoy protein, skip
					}
					String[] parts = pro.split("_");
					String genome = (parts.length > 0 && !parts[0].isEmpty()) ? parts[0] : null;
					if (genome != null && !genomesPempOne.contains(genome)) {
						kept++;
					}
				}
				if (kept > 0) {
					nPrecursorsRetained++;
					if (seq != null) {
						uniqueSeqRetained.add(seq);
					}
				}
			}
			double pctPrecursors = nPrecursorsQ01 > 0 ? 100.0 * nPrecursorsRetained / nPrecursorsQ01 : 0;
			double pctSequences = uniqueSeqQ01.size() > 0 ? 100.0 * uniqueSeqRetained.size() / uniqueSeqQ01.size() : 0;
			System.out.println("--- After removing proteins from genomes with p_emp=1 (precursors with Global.Q.Value < 0.01) ---");
			System.out.println("Precursors: " + nPrecursorsRetained + " / " + nPrecursorsQ01 + " retained (" + String.format("%.2f", pctPrecursors) + "%)");
			System.out.println("Unique sequences: " + uniqueSeqRetained.size() + " / " + uniqueSeqQ01.size() + " retained (" + String.format("%.2f", pctSequences) + "%)");

			// --- 7. Export run info (diagnostic numbers) ---
			int decoyPoolSize = decoyEPool.size();
			int nGenomes = genomePeptides.size();
			double minSid = Double.MAX_VALUE, maxSid = -Double.MAX_VALUE, sumSid = 0;
			HashMap<Double, Integer> sidCountMap = new HashMap<>();
			for (Double s : genomeScoreMap.values()) {
				minSid = Math.min(minSid, s);
				maxSid = Math.max(maxSid, s);
				sumSid += s;
				sidCountMap.merge(s, 1, Integer::sum);
			}
			double meanSid = nGenomes > 0 ? sumSid / nGenomes : 0;
			int nGenomesQ01 = 0;
			for (int i = 0; i < N; i++) {
				if (qvals[i] < 0.01) {
					nGenomesQ01++;
				}
			}
			double decoyEMin = Double.MAX_VALUE, decoyEMax = -Double.MAX_VALUE, decoyESum = 0;
			for (Double de : decoyEPool) {
				decoyEMin = Math.min(decoyEMin, de);
				decoyEMax = Math.max(decoyEMax, de);
				decoyESum += de;
			}
			double decoyEMean = decoyPoolSize > 0 ? decoyESum / decoyPoolSize : 0;

			// Intensity stats from precursors (for PEP+intensity scoring)
			double preNormMin = Double.MAX_VALUE, preNormMax = -Double.MAX_VALUE, preNormSum = 0;
			double preQuanMin = Double.MAX_VALUE, preQuanMax = -Double.MAX_VALUE, preQuanSum = 0;
			double ms1NormMin = Double.MAX_VALUE, ms1NormMax = -Double.MAX_VALUE, ms1NormSum = 0;
			double ms1AreaMin = Double.MAX_VALUE, ms1AreaMax = -Double.MAX_VALUE, ms1AreaSum = 0;
			for (DiaNNPrecursor pp : precursors) {
				double pn = pp.getPreTrans();
				double pq = pp.getPreQuan();
				double mn = pp.getMs1Normal();
				double ma = pp.getMs1Area();
				preNormMin = Math.min(preNormMin, pn);
				preNormMax = Math.max(preNormMax, pn);
				preNormSum += pn;
				preQuanMin = Math.min(preQuanMin, pq);
				preQuanMax = Math.max(preQuanMax, pq);
				preQuanSum += pq;
				ms1NormMin = Math.min(ms1NormMin, mn);
				ms1NormMax = Math.max(ms1NormMax, mn);
				ms1NormSum += mn;
				ms1AreaMin = Math.min(ms1AreaMin, ma);
				ms1AreaMax = Math.max(ms1AreaMax, ma);
				ms1AreaSum += ma;
			}
			int nP = precursors.length;
			double preNormMean = nP > 0 ? preNormSum / nP : 0;
			double preQuanMean = nP > 0 ? preQuanSum / nP : 0;
			double ms1NormMean = nP > 0 ? ms1NormSum / nP : 0;
			double ms1AreaMean = nP > 0 ? ms1AreaSum / nP : 0;

			// --- CScore and intensity distribution for diagnostics (target vs decoy) ---
			ArrayList<Double> targetCScore = new ArrayList<>();
			ArrayList<Double> decoyCScore = new ArrayList<>();
			ArrayList<Double> targetIntensity = new ArrayList<>();
			ArrayList<Double> decoyIntensity = new ArrayList<>();
			ArrayList<Double> targetPreQuan = new ArrayList<>();
			ArrayList<Double> decoyPreQuan = new ArrayList<>();
			ArrayList<Double> targetMs1Area = new ArrayList<>();
			ArrayList<Double> decoyMs1Area = new ArrayList<>();
			for (DiaNNPrecursor pp : precursors) {
				double cscore = pp.getCscore();
				double intensity = pp.getPreTrans();
				double preQuan = pp.getPreQuan();
				double ms1Area = pp.getMs1Area();
				if (pp.getIsDecoy()) {
					decoyCScore.add(cscore);
					decoyIntensity.add(intensity);
					decoyPreQuan.add(preQuan);
					decoyMs1Area.add(ms1Area);
				} else {
					targetCScore.add(cscore);
					targetIntensity.add(intensity);
					targetPreQuan.add(preQuan);
					targetMs1Area.add(ms1Area);
				}
			}
			// Export CScore and intensity distribution TSV
			File pepDistFile = new File(parquetFile.getParentFile(), "genome_rank_scores_cscore_distribution.tsv");
			try (PrintWriter pwPep = new PrintWriter(pepDistFile)) {
				pwPep.println("type\tCScore\tPrecursor.Normalised\tPrecursor.Quantity\tMs1.Area");
				for (int i = 0; i < targetCScore.size(); i++) {
					pwPep.println("target\t" + String.format("%.10f", targetCScore.get(i)) + "\t" + String.format("%.6f", targetIntensity.get(i))
							+ "\t" + String.format("%.6f", targetPreQuan.get(i)) + "\t" + String.format("%.6f", targetMs1Area.get(i)));
				}
				for (int i = 0; i < decoyCScore.size(); i++) {
					pwPep.println("decoy\t" + String.format("%.10f", decoyCScore.get(i)) + "\t" + String.format("%.6f", decoyIntensity.get(i))
							+ "\t" + String.format("%.6f", decoyPreQuan.get(i)) + "\t" + String.format("%.6f", decoyMs1Area.get(i)));
				}
			}
			// Summary stats for runinfo
			double[] targetCScoreStats = percentileStats(targetCScore);
			double[] decoyCScoreStats = percentileStats(decoyCScore);
			double[] targetIntStats = percentileStats(targetIntensity);
			double[] decoyIntStats = percentileStats(decoyIntensity);
			double[] targetPreQuanStats = percentileStats(targetPreQuan);
			double[] decoyPreQuanStats = percentileStats(decoyPreQuan);
			double[] targetMs1AreaStats = percentileStats(targetMs1Area);
			double[] decoyMs1AreaStats = percentileStats(decoyMs1Area);
			// Overlap: target [min,max] and decoy [min,max] intersect
			boolean cscoreOverlap = targetCScoreStats != null && decoyCScoreStats != null
					&& targetCScoreStats[0] <= decoyCScoreStats[1] && decoyCScoreStats[0] <= targetCScoreStats[1];

			ArrayList<String> parquetColumns = new ArrayList<>();
			try {
				parquetColumns = DiaNNParquetReader.getParquetColumnNames(parquetFile.getAbsolutePath());
			} catch (Exception ignore) {
			}
			List<String> intensityCols = new ArrayList<>();
			for (String col : parquetColumns) {
				String lc = col.toLowerCase();
				if (lc.contains("normalised") || lc.contains("quantity") || lc.contains("area")
						|| lc.contains("intensity") || lc.contains("translated") || lc.contains("lfq")) {
					intensityCols.add(col);
				}
			}

			File runInfoFile = new File(parquetFile.getParentFile(), "genome_rank_scores_runinfo.txt");
			try (PrintWriter pw = new PrintWriter(runInfoFile)) {
				pw.println("=================================================================================");
				pw.println("genome_rank_scores run info - diagnostic report");
				pw.println("=================================================================================");
				pw.println();
				pw.println("--- Input / Output ---");
				pw.println("input_parquet\t" + parquetFile.getAbsolutePath());
				pw.println("output_scores\t" + outFile.getAbsolutePath());
				pw.println();
				pw.println("--- Precursor / Peptide counts ---");
				pw.println("n_precursors_total\t" + nPrecursorsTotal + "\t(precursors with Q<0.01)");
				pw.println("n_peptides_target\t" + nTargetPeptides + "\t(after PSM collapse, mapping to >=1 genome)");
				pw.println("n_peptides_decoy\t" + nDecoyPeptides);
				pw.println();
				pw.println("--- Filter by p_emp=1 (precursors with Global.Q.Value < 0.01; exclude if all proteins from genomes with p_emp==1) ---");
				pw.println("n_genomes_p_emp_1\t" + genomesPempOne.size());
				pw.println("n_precursors_GlobalQ_lt_0.01\t" + nPrecursorsQ01);
				pw.println("n_precursors_retained\t" + nPrecursorsRetained + "\t(" + String.format("%.2f", pctPrecursors) + "%)");
				pw.println("n_unique_sequences_GlobalQ_lt_0.01\t" + uniqueSeqQ01.size());
				pw.println("n_unique_sequences_retained\t" + uniqueSeqRetained.size() + "\t(" + String.format("%.2f", pctSequences) + "%)");
				pw.println();
				pw.println("--- Peptide sharing (target) ---");
				pw.println("n_unique_peptides\t" + nUniquePeptides + "\t(m=1 genome)");
				pw.println("n_shared_peptides\t" + nSharedPeptides + "\t(m>1 genomes)");
				pw.print("m_distribution\t");
				List<Integer> mKeys = new ArrayList<>(mDistribution.keySet());
				Collections.sort(mKeys);
				for (int i = 0; i < mKeys.size(); i++) {
					if (i > 0) pw.print("; ");
					pw.print("m=" + mKeys.get(i) + ":" + mDistribution.get(mKeys.get(i)));
				}
				pw.println();
				pw.println();
				pw.println("--- CScore distribution (PSM-level, for diagnosing target/decoy overlap) ---");
				pw.println("cscore_distribution_file\t" + pepDistFile.getAbsolutePath());
				if (targetCScoreStats != null) {
					pw.println("target_CScore\tmin=" + String.format("%.8f", targetCScoreStats[0]) + "\tmax=" + String.format("%.8f", targetCScoreStats[1])
							+ "\tmean=" + String.format("%.8f", targetCScoreStats[2]) + "\tP1=" + String.format("%.8f", targetCScoreStats[3])
							+ "\tP50=" + String.format("%.8f", targetCScoreStats[6]) + "\tP99=" + String.format("%.8f", targetCScoreStats[9]));
				}
				if (decoyCScoreStats != null) {
					pw.println("decoy_CScore\tmin=" + String.format("%.8f", decoyCScoreStats[0]) + "\tmax=" + String.format("%.8f", decoyCScoreStats[1])
							+ "\tmean=" + String.format("%.8f", decoyCScoreStats[2]) + "\tP1=" + String.format("%.8f", decoyCScoreStats[3])
							+ "\tP50=" + String.format("%.8f", decoyCScoreStats[6]) + "\tP99=" + String.format("%.8f", decoyCScoreStats[9]));
				}
				pw.println("CScore_overlap\t" + (cscoreOverlap ? "YES (target and decoy CScore ranges intersect)" : "NO (no overlap - may cause p_emp floor)"));
				pw.println("e_score_note\tGenome scoring uses e = CScore (clamped to [0,6]). If decoy_e_max < target_e_min, decoy null cannot match target scores.");
				pw.println();
				pw.println("--- Intensity distribution (Precursor.Normalised) ---");
				if (targetIntStats != null) {
					pw.println("target_intensity\tmin=" + String.format("%.4f", targetIntStats[0]) + "\tmax=" + String.format("%.4f", targetIntStats[1])
							+ "\tmean=" + String.format("%.4f", targetIntStats[2]) + "\tP50=" + String.format("%.4f", targetIntStats[6]));
				}
				if (decoyIntStats != null) {
					pw.println("decoy_intensity\tmin=" + String.format("%.4f", decoyIntStats[0]) + "\tmax=" + String.format("%.4f", decoyIntStats[1])
							+ "\tmean=" + String.format("%.4f", decoyIntStats[2]) + "\tP50=" + String.format("%.4f", decoyIntStats[6]));
				}
				pw.println();
				pw.println("--- Precursor.Quantity distribution ---");
				if (targetPreQuanStats != null) {
					pw.println("target_Precursor.Quantity\tmin=" + String.format("%.4f", targetPreQuanStats[0]) + "\tmax=" + String.format("%.4f", targetPreQuanStats[1])
							+ "\tmean=" + String.format("%.4f", targetPreQuanStats[2]) + "\tP50=" + String.format("%.4f", targetPreQuanStats[6]));
				}
				if (decoyPreQuanStats != null) {
					pw.println("decoy_Precursor.Quantity\tmin=" + String.format("%.4f", decoyPreQuanStats[0]) + "\tmax=" + String.format("%.4f", decoyPreQuanStats[1])
							+ "\tmean=" + String.format("%.4f", decoyPreQuanStats[2]) + "\tP50=" + String.format("%.4f", decoyPreQuanStats[6]));
				}
				pw.println();
				pw.println("--- Ms1.Area distribution ---");
				if (targetMs1AreaStats != null) {
					pw.println("target_Ms1.Area\tmin=" + String.format("%.4f", targetMs1AreaStats[0]) + "\tmax=" + String.format("%.4f", targetMs1AreaStats[1])
							+ "\tmean=" + String.format("%.4f", targetMs1AreaStats[2]) + "\tP50=" + String.format("%.4f", targetMs1AreaStats[6]));
				}
				if (decoyMs1AreaStats != null) {
					pw.println("decoy_Ms1.Area\tmin=" + String.format("%.4f", decoyMs1AreaStats[0]) + "\tmax=" + String.format("%.4f", decoyMs1AreaStats[1])
							+ "\tmean=" + String.format("%.4f", decoyMs1AreaStats[2]) + "\tP50=" + String.format("%.4f", decoyMs1AreaStats[6]));
				}
				pw.println();
				pw.println("--- Decoy pool ---");
				pw.println("decoy_pep_max\t" + decoyPepMax + "\t(only decoys with PEP < this used in null pool)");
				pw.println("n_decoy_peptides_total\t" + nDecoyPeptides);
				pw.println("decoy_pool_size\t" + decoyPoolSize + "\t(n_decoy_in_pool=" + nDecoyInPool + ", PEP < " + decoyPepMax + ")");
				pw.println("decoy_e_min\t" + String.format("%.4f", decoyPoolSize > 0 ? decoyEMin : 0));
				pw.println("decoy_e_max\t" + String.format("%.4f", decoyPoolSize > 0 ? decoyEMax : 0));
				pw.println("decoy_e_mean\t" + String.format("%.4f", decoyEMean));
				pw.println();
				pw.println("--- Intensity (parquet columns, for CScore+intensity scoring) ---");
				pw.println("Precursor.Normalised\tmin=" + String.format("%.4f", nP > 0 ? preNormMin : 0)
						+ "\tmax=" + String.format("%.4f", nP > 0 ? preNormMax : 0)
						+ "\tmean=" + String.format("%.4f", preNormMean) + "\t(normalized precursor quantity)");
				pw.println("Precursor.Quantity\tmin=" + String.format("%.4f", nP > 0 ? preQuanMin : 0)
						+ "\tmax=" + String.format("%.4f", nP > 0 ? preQuanMax : 0)
						+ "\tmean=" + String.format("%.4f", preQuanMean) + "\t(MS2-based raw quantity)");
				pw.println("Ms1.Normalised\tmin=" + String.format("%.4f", nP > 0 ? ms1NormMin : 0)
						+ "\tmax=" + String.format("%.4f", nP > 0 ? ms1NormMax : 0)
						+ "\tmean=" + String.format("%.4f", ms1NormMean) + "\t(MS1 normalized)");
				pw.println("Ms1.Area\tmin=" + String.format("%.4f", nP > 0 ? ms1AreaMin : 0)
						+ "\tmax=" + String.format("%.4f", nP > 0 ? ms1AreaMax : 0)
						+ "\tmean=" + String.format("%.4f", ms1AreaMean) + "\t(MS1 raw peak area, if column present)");
				pw.print("parquet_intensity_related_columns\t");
				if (intensityCols.isEmpty()) {
					pw.println("(none)");
				} else {
					for (int i = 0; i < intensityCols.size(); i++) {
						if (i > 0) pw.print("; ");
						pw.print(intensityCols.get(i));
					}
					pw.println();
				}
				pw.println();
				pw.println("--- Scoring parameters ---");
				pw.println("peptide_mode\t" + (uniquePeptidesOnly ? "unique only (m=1, shared skipped)" : "all peptides (shared use weight w=1/m)"));
				pw.println("scoring_mode\tPEP/CScore+intensity (combined)");
				pw.println("intensity_factor\t0.5 + 0.25*norm(log2(1+I1)) + 0.25*norm(log2(1+I2))");
				pw.println("I1\tPrecursor.Quantity");
				pw.println("I2\tMs1.Area");
				pw.println("log2_I1_range\t" + String.format("%.2f", minLog1) + " to " + String.format("%.2f", maxLog1));
				pw.println("log2_I2_range\t" + String.format("%.2f", minLog2) + " to " + String.format("%.2f", maxLog2));
				pw.println("top_k\t" + k);
				pw.println("resamples_B\t" + B);
				pw.println("e_cap\t6.0");
				pw.println();
				pw.println("--- Genome scores (S) summary ---");
				pw.println("n_genomes\t" + nGenomes);
				pw.println("S_min\t" + String.format("%.4f", nGenomes > 0 ? minSid : 0));
				pw.println("S_max\t" + String.format("%.4f", nGenomes > 0 ? maxSid : 0));
				pw.println("S_mean\t" + String.format("%.4f", meanSid));
				pw.println("n_genomes_q_lt_0.01\t" + nGenomesQ01);
				pw.print("S_ties\t");
				List<Double> sidKeys = new ArrayList<>(sidCountMap.keySet());
				Collections.sort(sidKeys, Comparator.reverseOrder());
				for (int i = 0; i < sidKeys.size(); i++) {
					if (i > 0) pw.print("; ");
					pw.print(String.format("%.4f", sidKeys.get(i)) + ":" + sidCountMap.get(sidKeys.get(i)));
				}
				pw.println();
				pw.println();
				pw.println("--- Per-genome details (by rank) ---");
				pw.println("Format: Rank, Genome, S, n_peptides_total, n_top_k, p_emp, q_value");
				for (int i = 0; i < N; i++) {
					String g = sortedByScore.get(i);
					pw.println((i + 1) + "\t" + g + "\t" + String.format("%.4f", genomeScoreMap.get(g)) + "\t"
							+ genomeTotalPeptidesMap.get(g) + "\t" + genomeNMap.get(g)
							+ "\t" + String.format("%.6f", pempArr[i]) + "\t" + String.format("%.6f", qvals[i]));
				}
				pw.println();
				pw.println("--- Top-k peptides per genome ---");
				for (int i = 0; i < N; i++) {
					String g = sortedByScore.get(i);
					pw.println("## Rank " + (i + 1) + " " + g + " S=" + String.format("%.4f", genomeScoreMap.get(g)));
					List<WeightedPeptideScore> topK = genomeTopKMap.get(g);
					if (topK != null) {
						for (int j = 0; j < topK.size(); j++) {
							WeightedPeptideScore wps = topK.get(j);
							double mVal = (wps.w > 0) ? 1.0 / wps.w : 0;
							pw.println("  " + (j + 1) + "\t" + wps.seq + "\te=" + String.format("%.4f", wps.e)
									+ "\tf(I)=" + String.format("%.4f", wps.intensityFactor)
									+ "\tw=" + String.format("%.4f", wps.w) + "\tm=" + (int) mVal
									+ "\tw*e*f=" + String.format("%.4f", wps.we));
						}
					}
					pw.println();
				}
			}
		} catch (Exception e) {
			// propagate or log
		}
	}

	/** Returns [min, max, mean, p1, p5, p25, p50, p75, p95, p99] or null if empty. */
	private static double[] percentileStats(ArrayList<Double> list) {
		if (list == null || list.isEmpty()) {
			return null;
		}
		double[] arr = new double[list.size()];
		for (int i = 0; i < arr.length; i++) {
			arr[i] = list.get(i);
		}
		Arrays.sort(arr);
		double min = arr[0];
		double max = arr[arr.length - 1];
		double sum = 0;
		for (double v : arr) {
			sum += v;
		}
		double mean = sum / arr.length;
		double p = (double) (arr.length - 1);
		double p1 = arr[Math.min((int) Math.round(0.01 * p), arr.length - 1)];
		double p5 = arr[Math.min((int) Math.round(0.05 * p), arr.length - 1)];
		double p25 = arr[Math.min((int) Math.round(0.25 * p), arr.length - 1)];
		double p50 = arr[Math.min((int) Math.round(0.50 * p), arr.length - 1)];
		double p75 = arr[Math.min((int) Math.round(0.75 * p), arr.length - 1)];
		double p95 = arr[Math.min((int) Math.round(0.95 * p), arr.length - 1)];
		double p99 = arr[Math.min((int) Math.round(0.99 * p), arr.length - 1)];
		return new double[] { min, max, mean, p1, p5, p25, p50, p75, p95, p99 };
	}

	private static double[] bhCorrection(double[] p) {
		int n = p.length;
		Integer[] idx = new Integer[n];
		for (int i = 0; i < n; i++) {
			idx[i] = i;
		}
		Arrays.sort(idx, (a, b) -> Double.compare(p[a], p[b]));
		double[] q = new double[n];
		for (int i = 0; i < n; i++) {
			int ii = idx[i];
			q[ii] = Math.min(1.0, p[ii] * n / (i + 1));
		}
		for (int i = n - 2; i >= 0; i--) {
			q[idx[i]] = Math.min(q[idx[i]], q[idx[i + 1]]);
		}
		return q;
	}

	public static void main(String[] args) {
		testReadPrecursorsFromParquet(new File("Z:\\PRIDE_dataset\\PXD008738_Elo_DIA\\12mix\\MetaPilot\\focus\\combined"
				+ "\\combined.parquet"), false);
	}
}
