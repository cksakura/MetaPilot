package uk.ac.quadram.metapilot.task.dia;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.math.BigInteger;
import java.sql.SQLException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;

import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.math.NumberUtils;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.stat.inference.MannWhitneyUTest;
import org.dom4j.DocumentException;

import uk.ac.quadram.metapilot.core.math.MathTool;
import uk.ac.quadram.metapilot.core.prodb.FastaManager;
import uk.ac.quadram.metapilot.core.taxonomy.TaxonomyRanks;
import uk.ac.quadram.metapilot.core.tools.FormatTool;
import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNParquetReader;
import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNPrecursor;
import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNTask;
import uk.ac.quadram.metapilot.dbSearch.diann.DiannParameter;
import uk.ac.quadram.metapilot.dbSearch.diann.DiannPrMatrixReader;
import uk.ac.quadram.metapilot.quant.flashLFQ.FlashLfqQuanPeptide;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters.ToolType;
import uk.ac.quadram.metapilot.task.MetaPythonReportTask;
import uk.ac.quadram.metapilot.task.io.ReportManifestGenerator;
import uk.ac.quadram.metapilot.task.MetaPythonReportCopyTask;
import uk.ac.quadram.metapilot.task.dia.par.MetaParameterDia;
import uk.ac.quadram.metapilot.task.dia.par.MetaParameterDiaTg;
import uk.ac.quadram.metapilot.task.dia.par.MetaParameterDiaMarker;
import uk.ac.quadram.metapilot.task.par.MetaData;
import javafx.application.Platform;
import uk.ac.quadram.metapilot.task.io.pep.MetaPeptide;
import uk.ac.quadram.metapilot.task.io.pro.MetaProtein;
import uk.ac.quadram.metapilot.task.io.pro.MetaProteinAnnoEggNog;
import uk.ac.quadram.metapilot.task.io.pro.MetaProteinXMLReader2;
import uk.ac.quadram.metapilot.task.mag.HMGenomeHtmlWriter;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;
import uk.ac.quadram.metapilot.task.mag.MagFuncSearcher;
import uk.ac.quadram.metapilot.task.mag.MetaPilotMagTask;
import uk.ac.quadram.metapilot.task.mag.MetaProteinAnnoMag;
import javafx.scene.control.ProgressBar;

/**
 * Base class for DIA workflows. Abstract so that every task started from the UI
 * is a specific task (Marker, Seed Direct, Seed Iterative, Target, etc.).
 */
public abstract class MetaDiaTask extends MetaPilotMagTask {

	protected String[] diaFiles;
	protected DiaNNTask diaNNTask;
	protected DiaLibSearchPar diaLibSearchPar;
	protected DiannParameter diannPar;
	protected DiaNNPrecursor[] diaPp;
	protected MetaData metadata;
	protected String[] expNames;
	protected int[] fractions;

	protected File final_pro_xml_file;
	protected File reportHtmlDir;
	protected File final_summary_txt;
	protected File resultFolderFile;

	// Logger is inherited from MetaPilotMagTask -> MetaPilotTask as an instance
	// field
	private static final String taskName = "DIA data analysis";

	public MetaDiaTask(MetaParameterDia metaPar, MetaResourceParameters mrp, ProgressBar bar1, ProgressBar bar2) {
		super(metaPar, mrp, bar1, bar2);

		logger.info("Initializing MetaDiaTask...");
		logger.info("[{}] initializing task", getTaskName());

		this.metadata = metaPar.getMetadata();
		logger.info("Metadata loaded: {}", metadata != null ? "Yes" : "No");

		this.magDb = metaPar.getMagDbItem();
		if (magDb != null) {
			logger.info("MAG database loaded: {} (version: {})", magDb.getCatalogueID(),
					magDb.getUsedVersion() != null ? magDb.getUsedVersion().getCatalogueVersion() : "unknown");
		} else {
			logger.error("MAG database is null!");
		}

		this.threadCount = metaPar.getThreadCount();
		logger.info("Thread count set to: {}", threadCount);

		// Initialize DiaNN task
		String diaNNPath = mrp.getToolPath(ToolType.DIANN);
		logger.info("DiaNN path: {}", diaNNPath);

		// Validate DiaNN path
		if (diaNNPath == null || diaNNPath.trim().isEmpty()) {
			logger.error("DiaNN path is null or empty!");
			System.err.println("ERROR: DiaNN path is not configured in MetaResourceParameters");
			return;
		}

		File diaNNFile = new File(diaNNPath);
		if (!diaNNFile.exists()) {
			logger.error("DiaNN executable not found at path: {}", diaNNPath);
			System.err.println("ERROR: DiaNN executable not found at: " + diaNNPath);
			return;
		}

		this.diaNNTask = new DiaNNTask(diaNNPath);
		logger.info("DiaNN task initialized successfully");

		// Initialize DiaLibSearchPar
		this.diaLibSearchPar = new DiaLibSearchPar();
		this.diaLibSearchPar.setThreadCount(metaPar.getThreadCount());
		logger.info("DiaLibSearchPar initialized with thread count: {}", metaPar.getThreadCount());

		// Initialize DiannParameter
		this.diannPar = new DiannParameter();
		this.diannPar.setThreads(metaPar.getThreadCount());
		logger.info("DiannParameter initialized with thread count: {}", metaPar.getThreadCount());
		if (metaPar instanceof MetaParameterDiaMarker) {
			MetaParameterDiaMarker markerParams = (MetaParameterDiaMarker) metaPar;
			this.diannPar.setM_excision(markerParams.isMExcision());
			this.diannPar.setCarbami_C(markerParams.isCarbamiC());
			this.diannPar.setOxidation_M(markerParams.isOxidationM());
			this.diannPar.setAcety_proN(markerParams.isAcetyProN());
			this.diannPar.setPhospho_STY(markerParams.isPhosphoSTY());
			this.diannPar.setK_GG(markerParams.isKGG());
		} else if (metaPar instanceof MetaParameterDia) {
			MetaParameterDia diaParams = (MetaParameterDia) metaPar;
			this.diannPar.setM_excision(diaParams.isMExcision());
			this.diannPar.setCarbami_C(diaParams.isCarbamiC());
			this.diannPar.setOxidation_M(diaParams.isOxidationM());
			this.diannPar.setAcety_proN(diaParams.isAcetyProN());
			this.diannPar.setPhospho_STY(diaParams.isPhosphoSTY());
			this.diannPar.setK_GG(diaParams.isKGG());
		}

		// Set result directory if not already set (DIA Target and DIA Marker/Seed iterative set it here; Seed Direct sets it in initResultFolderIfNeeded())
		if (metaPar.getDbSearchResultFile() == null) {
			File diaResultDir = null;
			if (metaPar instanceof MetaParameterDiaTg) {
				diaResultDir = (new File(((MetaParameterDiaTg) metaPar).getLibraryPath())).getParentFile();
			} else if (metaPar instanceof MetaParameterDiaMarker) {
				if (this instanceof MetaDiaDeepComprehensiveTask) {
					diaResultDir = null;
				} else {
					MetaParameterDiaMarker markerParams = (MetaParameterDiaMarker) metaPar;
					String markerMode = markerParams.getMarkerMode();
					String folderName;
					if (this instanceof MetaDiaSeedIterativeTask) {
						if ("COMPREHENSIVE".equalsIgnoreCase(markerMode) || "DEEP".equalsIgnoreCase(markerMode)) {
							folderName = "seed_comprehensive";
						} else {
							folderName = "seed_focus";
						}
					} else if (this instanceof MetaDiaDeepMarkerTask) {
						folderName = "marker_deep";
					} else {
						if ("COMPREHENSIVE".equalsIgnoreCase(markerMode) || "DEEP".equalsIgnoreCase(markerMode)) {
							folderName = "marker_comprehensive";
						} else {
							folderName = "marker_focus";
						}
					}
					diaResultDir = new File(metaPar.getResult(), folderName);
					if (!diaResultDir.exists()) {
						diaResultDir.mkdirs();
						logger.info("Created result directory for DIA {} workflow ({}): {}",
								(this instanceof MetaDiaSeedIterativeTask) ? "Seed" : "Marker", folderName,
								diaResultDir.getAbsolutePath());
					}
				}
			}
			if (diaResultDir != null) {
				metaPar.setDbSearchResultFile(diaResultDir);
			} else {
				logger.info("Result directory will be set later by specific task implementation.");
			}
		} else {
			logger.info("Result directory already exist: {}", metaPar.getDbSearchResultFile().getAbsolutePath());
		}

		// Set the resultFolderFile field
		this.resultFolderFile = metaPar.getDbSearchResultFile();

		logger.info("Result folder file: {}", resultFolderFile != null ? resultFolderFile.getAbsolutePath() : "null");

		this.expNameMap = new HashMap<String, String>();
		logger.info("Experiment name map initialized");

		if (metadata != null) {
			this.rawFiles = metadata.getRawFiles();
			this.expNames = metadata.getExpNames();
			this.fractions = metadata.getFractions();

			logger.info("Metadata processing:");
			logger.info("  Raw files count: {}", rawFiles != null ? rawFiles.length : 0);
			logger.info("  Experiment names count: {}", expNames != null ? expNames.length : 0);
			logger.info("  Fractions count: {}", fractions != null ? fractions.length : 0);

			if (rawFiles != null && rawFiles.length > 0) {
				this.quanExpNames = new String[rawFiles.length];
				logger.info("Processing raw files:");
				for (int i = 0; i < rawFiles.length; i++) {
					this.quanExpNames[i] = rawFiles[i].substring(rawFiles[i].lastIndexOf("\\") + 1,
							rawFiles[i].lastIndexOf("."));
					this.expNameMap.put(quanExpNames[i], expNames[i]);
					logger.info("    Raw file {}: {} -> {}", i + 1, rawFiles[i], quanExpNames[i]);
				}
				logger.info("Experiment name mapping completed: {} entries", expNameMap.size());
			} else {
				logger.warn("No raw files found in metadata");
			}
		} else {
			logger.warn("Metadata is null - no raw files to process");
		}

		logger.info("MetaDiaTask initialization completed successfully");
		logger.info("[{}] initialization completed", getTaskName());
	}

	/**
	 * Initialize result folder if not already set. Subclasses override to create
	 * workflow-specific folders (e.g. seed_direct, marker_focus). Default: no-op.
	 */
	protected void initResultFolderIfNeeded() {
		// Base: no-op; subclasses set metaPar.setDbSearchResultFile() and resultFolderFile as needed.
	}

	protected boolean spectraConvert() {
		logger.info("Starting spectra conversion for {} files", rawFiles != null ? rawFiles.length : 0);
		logger.info("[{}] converting raw files started", getTaskName());

		ArrayList<String> list = new ArrayList<String>();
		this.diaFiles = new String[rawFiles.length];
		logger.info("Processing {} raw files for DIA format conversion", rawFiles.length);
		for (int i = 0; i < diaFiles.length; i++) {
			if (rawFiles[i].endsWith(".dia") || rawFiles[i].endsWith(".DIA")) {
				diaFiles[i] = rawFiles[i];
			} else if (rawFiles[i].endsWith(".d") || rawFiles[i].endsWith(".D")) {
				diaFiles[i] = rawFiles[i];
			} else if (rawFiles[i].endsWith(".mzml") || rawFiles[i].endsWith(".mzML")) {
				diaFiles[i] = rawFiles[i];
			} else if (rawFiles[i].endsWith(".raw") || rawFiles[i].endsWith(".RAW")) {
				diaFiles[i] = rawFiles[i] + ".dia";
				File diaFile = new File(diaFiles[i]);
				if (!diaFile.exists()) {
					list.add(rawFiles[i]);
				}
			} else {
				logger.error(getTaskName() + ": unknown spectra format");
				System.err.println(format.format(new Date()) + "\t" + getTaskName() + ": unknown spectra format");

				return false;
			}
		}

		if (list.size() > 0) {
			logger.info("Adding {} raw files to DiaNN conversion task", list.size());
			for (String file : list) {
				logger.info("  Converting: {} -> {}.dia", file, file);
			}
			diaNNTask.addTask(diannPar, list.toArray(String[]::new));
		}

		logger.info("Starting DiaNN conversion with {} threads", diaFiles.length);
		diaNNTask.run(1, diaFiles.length);
		logger.info("DiaNN conversion completed");

		for (int i = 0; i < diaFiles.length; i++) {
			File file = new File(diaFiles[i]);
			if (!file.exists()) {
				logger.error(getTaskName() + ": converting raw file to " + diaFiles[i] + " failed");
				logger.error("Expected output file does not exist: {}", file.getAbsolutePath());
				logger.error("Original raw file: {}", rawFiles[i]);
				System.err.println(format.format(new Date()) + "\t" + getTaskName() + ": converting raw file to "
						+ diaFiles[i] + " failed");
				System.err.println("Expected output file: " + file.getAbsolutePath());
				System.err.println("Original raw file: " + rawFiles[i]);
				System.err.println(
						"ERROR: cannot read .raw files, please download and install .NET Runtime .NET SDK 8.0.407 or later");
				System.err.println("https://dotnet.microsoft.com/en-us/download/dotnet/8.0");

				return false;
			} else {
				logger.info("Successfully converted: {} -> {} (size: {} bytes)", rawFiles[i], diaFiles[i],
						file.length());
			}
		}

		logger.info(getTaskName() + ": converting raw files finished");
		logger.info("[{}] converting raw files finished", getTaskName());

		return true;
	}

	protected boolean createSSLibrary() {
		logger.info("Starting sample-specific library creation...");
		updateProgress(0, 100);
		// Base implementation: library-search workflows (Target, Seed Direct) use provided libraries; no self-library.
		logger.info("Skipping sample-specific library (use provided libraries).");
		return true;
	}

	protected boolean iden() {
		logger.info("Starting identification process...");
		updateMessage(getTaskName() + ": peptide and protein identification");
		logger.info("[{}] peptide and protein identification started", getTaskName());

		logger.info("Configuring DiaLibSearchPar:");
		logger.info("  DIA files count: {}", diaFiles != null ? diaFiles.length : 0);
		logger.info("  Thread count: {}", threadCount);
		logger.info("  Q-value threshold: 0.01");
		logger.info("  MBR enabled: true");

		diaLibSearchPar.setMbr(true);
		diaLibSearchPar.setDiaFiles(diaFiles);
		diaLibSearchPar.setThreadCount(threadCount);
		diaLibSearchPar.setqValue(0.01);

		diannPar.setThreads(threadCount);
		diannPar.setMatrices(true);
		diannPar.setRelaxed_prot_inf(true);
		diannPar.setReanalyse(true);

		File combinedFile = new File(metaPar.getDbSearchResultFile(), "combined");
		if (!combinedFile.exists()) {
			combinedFile.mkdir();
		}
		this.quan_pep_file = new File(combinedFile, "combined.pr_matrix.tsv");
		this.quan_pro_file = new File(combinedFile, "combined.pg_matrix.tsv");

		File combineParquetFile = new File(combinedFile, "combined.parquet");

		if (!quan_pep_file.exists() || !combineParquetFile.exists()) {

			logger.info(getTaskName() + ": library search started");
			logger.info("[{}] library search started", getTaskName());

			updateProgress(10, 100);
			diannPar.setExportDecoy(true);

			String[] libPath = null;
			if (metaPar instanceof MetaParameterDiaTg) {
				String libraryPath = ((MetaParameterDiaTg) metaPar).getLibraryPath();
				libPath = new String[] { libraryPath };
				diaNNTask.addTask(diannPar, diaLibSearchPar, combineParquetFile.getAbsolutePath(), libPath);
			} else {
				File magDbFile = this.magDb.getCurrentFile();
				File libFolderFile = new File(magDbFile, "libraries");

				String[] usedLibs = ((MetaParameterDia) metaPar).getLibrary();
				libPath = new String[usedLibs.length];
				for (int i = 0; i < usedLibs.length; i++) {
					File directLibFile = new File(usedLibs[i]);
					if (directLibFile.exists()) {
						libPath[i] = directLibFile.getAbsolutePath();
					} else {
						File libFolder = new File(libFolderFile, usedLibs[i]);
						File speclibFile = new File(libFolder, usedLibs[i] + ".speclib");
						File tsvFile = new File(libFolder, usedLibs[i] + ".tsv");
						File parquetFile = new File(libFolder, usedLibs[i] + ".parquet");

						File chosen = null;
						if (speclibFile.exists()) {
							chosen = speclibFile;
						} else if (tsvFile.exists()) {
							chosen = tsvFile;
						} else if (parquetFile.exists()) {
							chosen = parquetFile;
						}

						if (chosen == null) {
							logger.error(getTaskName() + ": no library file (.speclib/.tsv/.parquet) found for library '{}'", usedLibs[i]);
							System.err.println(format.format(new Date()) + "\t" + getTaskName()
									+ ": no library file (.speclib/.tsv/.parquet) found for library '" + usedLibs[i] + "'");
							return false;
						}

						libPath[i] = chosen.getAbsolutePath();
					}
				}
				diaNNTask.addTask(diannPar, diaLibSearchPar, combineParquetFile.getAbsolutePath(), libPath);
			}
			diaNNTask.run(diaFiles.length * 12);

			// If combined result missing (e.g. one damaged raw file), retry by searching each file separately
			if (!combineParquetFile.exists() && diaFiles != null && diaFiles.length > 1) {
				String[] fallbackLibPath = metaPar instanceof MetaParameterDiaTg
						? new String[] { ((MetaParameterDiaTg) metaPar).getLibraryPath() }
						: libPath;
				if (fallbackLibPath != null) {
					List<String> damaged = runPerFileFallbackAndRetryCombined(combineParquetFile.getAbsolutePath(), fallbackLibPath);
					if (!damaged.isEmpty()) {
						logger.warn(getTaskName() + ": The following raw files produced no results (may be damaged): {}", damaged);
						System.err.println(format.format(new Date()) + "\t" + getTaskName()
								+ ": The following raw files produced no results (may be damaged): " + String.join(", ", damaged));
					}
				}
			}

			if (!combineParquetFile.exists()) {
				logger.error(getTaskName() + ": identification failed - result files not created");
				return false;
			}

			logger.info(getTaskName() + ": library search finished");
			logger.info("[{}] library search finished", getTaskName());
		}

		cleanTempFiles();

		// Stop the watch service thread
		stopWatchServiceThread();

		updateProgress(85, 100);

		return true;
	}

	/**
	 * When combined search produced no output (e.g. one damaged raw file), run each raw file
	 * separately, then re-run combined search with only successful files. Caller should log
	 * the returned damaged list to the user.
	 *
	 * @param combinedParquetPath path for the final combined.parquet
	 * @param libPath              library path(s) for the search
	 * @return list of raw file paths that produced no result (may be damaged)
	 */
	protected List<String> runPerFileFallbackAndRetryCombined(String combinedParquetPath, String[] libPath) {
		if (diaFiles == null || diaFiles.length == 0 || libPath == null || libPath.length == 0) {
			return new ArrayList<>();
		}
		File combinedDir = new File(combinedParquetPath).getParentFile();
		if (combinedDir == null || !combinedDir.exists()) {
			return new ArrayList<>();
		}
		logger.info(getTaskName() + ": combined search produced no output; retrying with per-file search to isolate damaged files");
		System.out.println(format.format(new Date()) + "\t" + getTaskName()
				+ ": retrying with per-file search (some raw files may be damaged)");

		// Run one search per raw file, but skip those that already have a non-empty per-file result
		int queued = 0;
		for (int i = 0; i < diaFiles.length; i++) {
			File pf = new File(combinedDir, "per_file_" + i + ".parquet");
			if (pf.exists() && pf.length() > 0) {
				continue; // already completed in a previous attempt
			}
			String perOut = pf.getAbsolutePath();
			diaLibSearchPar.setDiaFiles(new String[] { diaFiles[i] });
			diaNNTask.addTask(diannPar, diaLibSearchPar, perOut, libPath);
			queued++;
		}
		if (queued > 0) {
			diaNNTask.run(queued * 12);
		} else {
			logger.info(getTaskName() + ": all per-file parquet outputs already exist; skipping new per-file DIA-NN runs");
		}

		List<String> goodFiles = new ArrayList<>();
		List<String> damaged = new ArrayList<>();
		for (int i = 0; i < diaFiles.length; i++) {
			File pf = new File(combinedDir, "per_file_" + i + ".parquet");
			if (pf.exists() && pf.length() > 0) {
				goodFiles.add(diaFiles[i]);
			} else {
				damaged.add(diaFiles[i]);
			}
		}

		if (goodFiles.isEmpty()) {
			logger.error(getTaskName() + ": per-file fallback found no successful files");
			return damaged;
		}
		logger.info(getTaskName() + ": re-running combined search with {} successful file(s) ({} damaged or empty)",
				goodFiles.size(), damaged.size());
		diaLibSearchPar.setDiaFiles(goodFiles.toArray(new String[0]));
		diaNNTask.addTask(diannPar, diaLibSearchPar, combinedParquetPath, libPath);
		diaNNTask.run(goodFiles.size() * 12);
		return damaged;
	}

	/**
	 * Same as {@link #runPerFileFallbackAndRetryCombined(String, String[])} but uses addNoQuantLibTask
	 * (for first search: marker library or user library, no quantification).
	 * <p>
	 * NOTE: Unlike {@link #runPerFileFallbackAndRetryCombined(String, String[])}, this helper
	 * does <b>not</b> re-run a combined search. It only runs per-file searches and returns
	 * which raw files were damaged. Callers are expected to consume the per-file parquet
	 * outputs directly (e.g. by reading per_file_*.parquet).
	 */
	protected List<String> runPerFileFallbackAndRetryCombinedNoQuant(String outputParquetPath, String[] libPath) {
		if (diaFiles == null || diaFiles.length == 0 || libPath == null || libPath.length == 0) {
			return new ArrayList<>();
		}
		File outDir = new File(outputParquetPath).getParentFile();
		if (outDir == null || !outDir.exists()) {
			return new ArrayList<>();
		}
		logger.info(getTaskName() + ": first search produced no output; retrying with per-file search to isolate damaged files");
		System.out.println(format.format(new Date()) + "\t" + getTaskName()
				+ ": retrying with per-file search (some raw files may be damaged)");

		int queued = 0;
		for (int i = 0; i < diaFiles.length; i++) {
			File pf = new File(outDir, "per_file_" + i + ".parquet");
			if (pf.exists() && pf.length() > 0) {
				continue; // this raw file already has a per-file result
			}
			String perOut = pf.getAbsolutePath();
			diaLibSearchPar.setDiaFiles(new String[] { diaFiles[i] });
			diaNNTask.addNoQuantLibTask(diannPar, diaLibSearchPar, perOut, libPath);
			queued++;
		}
		if (queued > 0) {
			diaNNTask.run(queued * 12);
		} else {
			logger.info(getTaskName() + ": per-file fallback: all per_file_*.parquet already exist; no new DIA-NN runs scheduled");
		}

		List<String> goodFiles = new ArrayList<>();
		List<String> damaged = new ArrayList<>();
		for (int i = 0; i < diaFiles.length; i++) {
			File pf = new File(outDir, "per_file_" + i + ".parquet");
			if (pf.exists() && pf.length() > 0) {
				goodFiles.add(diaFiles[i]);
			} else {
				damaged.add(diaFiles[i]);
			}
		}

		if (goodFiles.isEmpty()) {
			logger.error(getTaskName() + ": per-file fallback found no successful files");
			return damaged;
		}
		logger.info(getTaskName() + ": per-file fallback completed with {} successful file(s) ({} damaged or empty); callers should consume per_file_*.parquet directly",
				goodFiles.size(), damaged.size());
		return damaged;
	}

	/**
	 * Same as {@link #runPerFileFallbackAndRetryCombinedNoQuant(String, String[])} but uses
	 * addTask with genLib (for first search with --out-lib, e.g. mixed microbiome+standalone).
	 * <p>
	 * NOTE: This helper also only runs per-file searches and does not re-run a combined
	 * first search. Callers are expected to read per_file_*.parquet (and any corresponding
	 * libraries) directly instead of relying on a single combined parquet.
	 */
	protected List<String> runPerFileFallbackAndRetryCombinedWithGenLib(String outputParquetPath, String[] libPath) {
		if (diaFiles == null || diaFiles.length == 0 || libPath == null || libPath.length == 0) {
			return new ArrayList<>();
		}
		File outDir = new File(outputParquetPath).getParentFile();
		if (outDir == null || !outDir.exists()) {
			return new ArrayList<>();
		}
		logger.info(getTaskName() + ": first search produced no output; retrying with per-file search to isolate damaged files");
		System.out.println(format.format(new Date()) + "\t" + getTaskName()
				+ ": retrying with per-file search (some raw files may be damaged)");

		int queued = 0;
		for (int i = 0; i < diaFiles.length; i++) {
			File pf = new File(outDir, "per_file_" + i + ".parquet");
			if (pf.exists() && pf.length() > 0) {
				continue; // already completed
			}
			String perOut = pf.getAbsolutePath();
			diaLibSearchPar.setDiaFiles(new String[] { diaFiles[i] });
			diaNNTask.addTask(diannPar, diaLibSearchPar, perOut, false, true, true, libPath);
			queued++;
		}
		if (queued > 0) {
			diaNNTask.run(queued * 12);
		} else {
			logger.info(getTaskName() + ": per-file fallback with genLib: all per_file_*.parquet already exist; no new DIA-NN runs scheduled");
		}

		List<String> goodFiles = new ArrayList<>();
		List<String> damaged = new ArrayList<>();
		for (int i = 0; i < diaFiles.length; i++) {
			File pf = new File(outDir, "per_file_" + i + ".parquet");
			if (pf.exists() && pf.length() > 0) {
				goodFiles.add(diaFiles[i]);
			} else {
				damaged.add(diaFiles[i]);
			}
		}

		if (goodFiles.isEmpty()) {
			logger.error(getTaskName() + ": per-file fallback found no successful files");
			return damaged;
		}
		logger.info(getTaskName() + ": per-file fallback (with genLib) completed with {} successful file(s) ({} damaged or empty); callers should consume per_file_*.parquet directly",
				goodFiles.size(), damaged.size());
		return damaged;
	}

	@SuppressWarnings("unused")
	private HashMap<String, HashSet<Integer>> exportPeptideTxt(FlashLfqQuanPeptide[] peps,
			HashMap<String, HashSet<String>> pepProMap, HashMap<String, DiaNNPrecursor> pepPpMap,
			HashMap<String, Integer> proIdMap) {

		logger.info(getTaskName() + ": exporting peptide report started");
		logger.info("[{}] exporting peptide report started", getTaskName());

		this.final_pep_txt = new File(resultFolderFile, "final_peptides.tsv");
		HashMap<String, HashSet<Integer>> proPepIdMap = new HashMap<String, HashSet<Integer>>();

		PrintWriter writer = null;
		try {
			writer = new PrintWriter(this.final_pep_txt);

			StringBuilder titlesb = new StringBuilder();
			titlesb.append("Sequence").append("\t");
			titlesb.append("Base sequence").append("\t");
			titlesb.append("Length").append("\t");
			titlesb.append("Missed cleavages").append("\t");
			titlesb.append("Mass").append("\t");
			titlesb.append("Proteins").append("\t");
			titlesb.append("Charges").append("\t");
			titlesb.append("Score").append("\t");
			titlesb.append("PEP").append("\t");
			titlesb.append("Q-value").append("\t");
			titlesb.append("Intensity").append("\t");

			for (int i = 0; i < quanExpNames.length; i++) {
				if (this.expNameMap.containsKey(quanExpNames[i])) {
					titlesb.append("Intensity " + expNameMap.get(quanExpNames[i])).append("\t");
				} else {
					titlesb.append("Intensity " + quanExpNames[i]).append("\t");
				}
			}

			titlesb.append("id").append("\t");
			titlesb.append("Protein group IDs").append("\t");

			writer.println(titlesb);

			int pepid = 0;
			for (int i = 0; i < peps.length; i++) {

				String sequence = peps[i].getSequence();
				String modSeq = peps[i].getModSeq();

				if (pepPpMap.containsKey(modSeq)) {
					DiaNNPrecursor diannPp = pepPpMap.get(modSeq);
					StringBuilder sb = new StringBuilder();
					sb.append(modSeq).append("\t");
					sb.append(sequence).append("\t");
					sb.append(sequence.length()).append("\t");
					sb.append(diannPp.getMissCleavage()).append("\t");
					sb.append(diannPp.getMass()).append("\t");

					HashSet<String> proSet = pepProMap.get(sequence);
					StringBuilder pidSb = new StringBuilder();
					if (proSet != null && proSet.size() > 0) {
						for (String pro : proSet) {
							sb.append(pro).append(";");
							if (proIdMap.containsKey(pro)) {
								pidSb.append(proIdMap.get(pro)).append(";");

								if (proPepIdMap.containsKey(pro)) {
									proPepIdMap.get(pro).add(i);
								} else {
									HashSet<Integer> pset = new HashSet<Integer>();
									pset.add(i);
									proPepIdMap.put(pro, pset);
								}
							}
						}
						sb.deleteCharAt(sb.length() - 1);
						if (pidSb.length() > 1) {
							pidSb.deleteCharAt(pidSb.length() - 1);
						}
					} else {
						continue;
					}

					sb.append("\t").append(diannPp.getCharge());
					sb.append("\t").append(diannPp.getCscore());
					sb.append("\t").append(diannPp.getPEP());
					sb.append("\t").append(diannPp.getQvalue());

					peps[i].setScore(diannPp.getCscore());

					double totalIntensity = 0;
					double[] intensity = peps[i].getIntensity();
					for (int j = 0; j < intensity.length; j++) {
						totalIntensity += intensity[j];
					}

					sb.append("\t").append(totalIntensity);
					for (int j = 0; j < intensity.length; j++) {
						sb.append("\t").append(intensity[j]);
					}

					sb.append("\t").append(pepid++);
					sb.append("\t").append(pidSb);
					writer.println(sb);
				}
			}

			writer.close();

		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			logger.error(getTaskName() + ": error in writing peptides to " + this.final_pep_txt.getName(), e);
			System.err.println(format.format(new Date()) + "\t" + getTaskName() + ": error in writing peptides to "
					+ this.final_pep_txt.getName());
		}

		logger.info(getTaskName() + ": peptide report has been exported to " + final_pep_txt);
		logger.info("[{}] peptide report has been exported to {}", getTaskName(), final_pep_txt);

		return proPepIdMap;
	}

	/**
	 * 
	 */
	@SuppressWarnings("unused")
	private void exportProteinsTxt(FlashLfqQuanPeptide[] peps, String[] allPros,
			HashMap<String, ProteinGroup> proGroupMap, HashMap<String, HashSet<String>> razorProPepMap,
			HashMap<String, HashSet<String>> allProPepMap, HashMap<String, double[]> proIntensityMap,
			HashMap<String, HashSet<Integer>> proPepIdMap) {

		logger.info(getTaskName() + ": exporting protein report started");
		logger.info("[{}] exporting protein report started", getTaskName());

		this.final_pro_txt = new File(resultFolderFile, "final_proteins.tsv");

		PrintWriter writer = null;
		try {
			writer = new PrintWriter(this.final_pro_txt);

			StringBuilder titlesb = new StringBuilder();
			titlesb.append("Protein IDs").append("\t");
			titlesb.append("Majority protein IDs").append("\t");
			titlesb.append("Peptide counts (all)").append("\t");
			titlesb.append("Peptide counts (razor)").append("\t");
			titlesb.append("Number of proteins").append("\t");
			titlesb.append("Peptides").append("\t");
			titlesb.append("Score").append("\t");
			titlesb.append("Intensity").append("\t");

			for (int i = 0; i < quanExpNames.length; i++) {
				if (this.expNameMap.containsKey(quanExpNames[i])) {
					titlesb.append("Intensity " + expNameMap.get(quanExpNames[i])).append("\t");
				} else {
					titlesb.append("Intensity " + quanExpNames[i]).append("\t");
				}
			}

			titlesb.append("id").append("\t");
			titlesb.append("Peptide IDs").append("\t");
			titlesb.append("Peptide is razor");

			writer.println(titlesb);

			DecimalFormat df = FormatTool.getDFE4();
			for (int i = 0; i < allPros.length; i++) {

				StringBuilder sb = new StringBuilder();
				double[] intensity = proIntensityMap.get(allPros[i]);

				if (intensity == null) {
					intensity = new double[quanExpNames.length];
					Arrays.fill(intensity, 0.0);
				}

				ProteinGroup pg = proGroupMap.get(allPros[i]);
				HashSet<String> sameSet = pg.getSameSet();
				HashSet<Integer> pepIdSet = proPepIdMap.get(allPros[i]);

				Integer[] pepIds = pepIdSet.toArray(new Integer[pepIdSet.size()]);
				Arrays.sort(pepIds);
				boolean[] isRazorPep = new boolean[pepIds.length];

				int proCount = 1 + sameSet.size();
				int[] pepCountAll = new int[proCount];
				int[] pepCountRazor = new int[proCount];
				double proScore = 0;

				for (int j = 0; j < pepIds.length; j++) {
					String sequence = peps[pepIds[j]].getSequence();
					pepCountAll[0]++;
					proScore += peps[pepIds[j]].getScore();
					if (razorProPepMap.get(allPros[i]).contains(sequence)) {
						pepCountRazor[0]++;
						isRazorPep[j] = true;
					}
				}

				int proNameId = 1;
				StringBuilder sameSb = new StringBuilder();
				for (String p : sameSet) {
					sameSb.append(";").append(p);
					for (int j = 0; j < pepIds.length; j++) {
						String sequence = peps[pepIds[j]].getSequence();
						if (allProPepMap.get(p).contains(sequence)) {
							pepCountAll[proNameId]++;
							if (isRazorPep[j]) {
								pepCountRazor[proNameId]++;
							}
						}
					}
					proNameId++;
				}

				if (sameSb.length() > 0) {
					sameSb.deleteCharAt(sameSb.length() - 1);
				}

				sb.append(allPros[i]).append(sameSb).append("\t");
				sb.append(allPros[i]).append(sameSb).append("\t");

				for (int j = 0; j < pepCountAll.length; j++) {
					sb.append(pepCountAll[j]).append(";");
				}
				sb.deleteCharAt(sb.length() - 1);
				sb.append("\t");

				for (int j = 0; j < pepCountRazor.length; j++) {
					sb.append(pepCountRazor[j]).append(";");
				}
				sb.deleteCharAt(sb.length() - 1);
				sb.append("\t");

				sb.append(proCount).append("\t");
				sb.append(pepIds.length).append("\t");

				sb.append(FormatTool.getDF2().format(proScore)).append("\t");

				double totalIntensity = 0;

				for (int j = 0; j < intensity.length; j++) {
					totalIntensity += intensity[j];
				}

				sb.append(df.format(totalIntensity)).append("\t");
				for (int j = 0; j < intensity.length; j++) {
					sb.append(df.format(intensity[j])).append("\t");
				}

				sb.append(i).append("\t");

				for (int j = 0; j < pepIds.length; j++) {
					sb.append(pepIds[j]).append(";");
				}
				sb.deleteCharAt(sb.length() - 1);
				sb.append("\t");

				for (int j = 0; j < isRazorPep.length; j++) {
					sb.append(isRazorPep[j]).append(";");
				}
				sb.deleteCharAt(sb.length() - 1);
				writer.println(sb);
			}

			writer.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			logger.error(getTaskName() + ": error in writing final protein result to " + this.final_pro_txt.getName(),
					e);
			System.err.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": error in writing final protein result to " + this.final_pro_txt.getName());
		}

		logger.info(getTaskName() + ": protein report has been exported to " + final_pro_txt.getName());
		logger.info("[{}] protein report has been exported to {}", getTaskName(), final_pro_txt.getName());
	}

	/**
	 * Default: no implementation. Subclasses (Marker, Seed, Target, etc.) override for their workflows.
	 */
	protected boolean iterativeSearch() {
		return false;
	}

	protected String getVersion() {
		return null;
	}

	/**
	 * Reverser proteins are also included to calculate the FDR of genomes
	 * 
	 * @param quanPeps
	 * @return
	 */
	@SuppressWarnings("unchecked")
	protected HashMap<String, double[]> refineGenomes(DiaNNPrecursor[] pps, String decoySymble, File genomeScoreFile) {

		String revIdentifier = decoySymble + this.magDb.getIdentifier();
		ArrayList<Double> globalDecoyList = new ArrayList<Double>();
		HashMap<String, ArrayList<Double>[]> genomeTDMap = new HashMap<String, ArrayList<Double>[]>();
		for (int i = 0; i < pps.length; i++) {
			double PEPi = pps[i].getPreQuan();
			if (PEPi == 0) {
				continue;
			}

			boolean isDecoy = false;
			HashSet<String> genomeSet = new HashSet<String>();
			String[] proteins = pps[i].getProGroups();
			for (String pro : proteins) {
				if (pro.startsWith(revIdentifier)) {
					String genome = (pro.substring(decoySymble.length()).split("_")[0]);
					genomeSet.add(genome);
					isDecoy = true;
				} else if (pro.startsWith(magDb.getIdentifier())) {
					String genome = pro.split("_")[0];
					genomeSet.add(genome);
				}
			}
			if (genomeSet.size() > 1) {
				continue;
			}

			if (isDecoy) {
				globalDecoyList.add(PEPi);
			}
			for (String genome : genomeSet) {
				if (isDecoy) {
					if (genomeTDMap.containsKey(genome)) {
						genomeTDMap.get(genome)[1].add(PEPi);
					} else {
						ArrayList<Double> targetList = new ArrayList<Double>();
						targetList.add(PEPi);
						ArrayList<Double> decoyList = new ArrayList<Double>();
						genomeTDMap.put(genome, new ArrayList[] { targetList, decoyList });
					}
				} else {
					if (genomeTDMap.containsKey(genome)) {
						genomeTDMap.get(genome)[0].add(PEPi);
					} else {
						ArrayList<Double> targetList = new ArrayList<Double>();
						targetList.add(PEPi);
						ArrayList<Double> decoyList = new ArrayList<Double>();
						genomeTDMap.put(genome, new ArrayList[] { targetList, decoyList });
					}
				}
			}
		}

		double[] globalDecoys = new double[globalDecoyList.size()];
		for (int i = 0; i < globalDecoys.length; i++) {
			globalDecoys[i] = Math.log10(globalDecoyList.get(i));
		}
		DescriptiveStatistics decoyStats = new DescriptiveStatistics(globalDecoys);
		double decoyMean = decoyStats.getMean();
		double decoySd = decoyStats.getStandardDeviation();

		MannWhitneyUTest mwTest = new MannWhitneyUTest();

		HashMap<String, double[]> genomeScoreMap = new HashMap<String, double[]>();
		DecimalFormat df4 = FormatTool.getDF4();
		DecimalFormat dfe4 = FormatTool.getDFE4();
		try (PrintWriter writer = new PrintWriter(genomeScoreFile)) {
			writer.println("Genome\tUnique_peptide_count\tDecoy_peptide_count\tTarget_mean\tDecoy_mean\tCohen's_d"
					+ "\tMann-Whitney_U\tMann-Whitney U p-value\tCliff's_delta");
			for (String genome : genomeTDMap.keySet()) {
				ArrayList<Double> targetList = genomeTDMap.get(genome)[0];
				ArrayList<Double> decoyList = genomeTDMap.get(genome)[1];

				double[] targetLogInten = new double[targetList.size()];
				double[] decoyLogInten = new double[decoyList.size()];

				for (int i = 0; i < targetLogInten.length; i++) {
					targetLogInten[i] = Math.log10(targetList.get(i));
				}

				for (int i = 0; i < decoyLogInten.length; i++) {
					decoyLogInten[i] = Math.log10(decoyList.get(i));
				}

				if (targetLogInten.length > 0) {
					DescriptiveStatistics stats1 = new DescriptiveStatistics(targetLogInten);
					double mean1 = stats1.getMean();
					double sd1 = stats1.getStandardDeviation();
					double pooledSD = Math
							.sqrt(((stats1.getN() - 1) * sd1 * sd1 + (decoyStats.getN() - 1) * decoySd * decoySd)
									/ (stats1.getN() + decoyStats.getN() - 2));
					double cohensD = (mean1 - decoyMean) / pooledSD;
					double mwU = mwTest.mannWhitneyU(targetLogInten, globalDecoys);
					double mwPValue = mwTest.mannWhitneyUTest(targetLogInten, globalDecoys); // p-value

					if (mean1 > decoyMean) {
						writer.println(genome + "\t" + targetLogInten.length + "\t" + globalDecoys.length + "\t"
								+ df4.format(mean1) + "\t" + df4.format(decoyMean) + "\t" + df4.format(cohensD) + "\t"
								+ df4.format(mwU) + "\t" + dfe4.format(mwPValue) + "\t"
								+ df4.format(MathTool.calculateCliffsDelta(targetLogInten, globalDecoys)));

						if (mwPValue < 0.01) {
							genomeScoreMap.put(genome, new double[] { mwPValue, 0.0 });
						}
					}
				}

				if (decoyLogInten.length > 0) {
					DescriptiveStatistics stats2 = new DescriptiveStatistics(decoyLogInten);
					double mean2 = stats2.getMean();
					double sd2 = stats2.getStandardDeviation();
					double pooledSD = Math
							.sqrt(((stats2.getN() - 1) * sd2 * sd2 + (decoyStats.getN() - 1) * decoySd * decoySd)
									/ (stats2.getN() + decoyStats.getN() - 2));
					double cohensD = (mean2 - decoyMean) / pooledSD;
					double mwU = mwTest.mannWhitneyU(decoyLogInten, globalDecoys);
					double mwPValue = mwTest.mannWhitneyUTest(decoyLogInten, globalDecoys); // p-value

					if (mean2 > decoyMean) {
						writer.println(FastaManager.shuffle + genome + "\t" + decoyLogInten.length + "\t"
								+ globalDecoys.length + "\t" + df4.format(mean2) + "\t" + df4.format(decoyMean) + "\t"
								+ df4.format(cohensD) + "\t" + df4.format(mwU) + "\t" + dfe4.format(mwPValue) + "\t"
								+ df4.format(MathTool.calculateCliffsDelta(decoyLogInten, globalDecoys)));

						if (mwPValue < 0.01) {
							genomeScoreMap.put(FastaManager.shuffle + genome, new double[] { mwPValue, 0.0 });
						}
					}
				}
			}
			writer.close();
		} catch (Exception e) {

		}

		String[] genomes = genomeScoreMap.keySet().toArray(String[]::new);
		Arrays.sort(genomes, Comparator.comparingDouble(genome -> genomeScoreMap.get(genome)[0]));
		int targetCount = 0;
		for (int i = 0; i < genomes.length; i++) {
			if (!genomes[i].startsWith(FastaManager.shuffle)) {
				targetCount++;
			}
			if (targetCount == 0) {
				genomeScoreMap.get(genomes[i])[1] = 1.0;
			} else {
				genomeScoreMap.get(genomes[i])[1] = (double) (i + 1 - targetCount) / (double) targetCount;
			}
		}

		return genomeScoreMap;
	}

	/**
	 * Reverser proteins are also included to calculate the FDR of genomes
	 * 
	 * @param quanPeps
	 * @return
	 */
	protected void calGenomeProteinScore(DiaNNPrecursor[] pps, HashMap<String, double[]> genomeScoreMap,
			HashMap<String, double[]> proteinScoreMap, File genomeScoreFile, File proteinScoreFile) {

		ArrayList<Double> globalDecoyList = new ArrayList<Double>();
		HashSet<String> fileSet = new HashSet<>();
		HashMap<String, ArrayList<Double>> genomeTDMap = new HashMap<String, ArrayList<Double>>();

		L: for (int i = 0; i < pps.length; i++) {
			double PEPi = pps[i].getPreQuan();
			if (PEPi == 0) {
				continue;
			}
			fileSet.add(pps[i].getRunString());

			boolean isTarget = true;
			String genome;
			String[] proteins = pps[i].getProGroups();
			if (proteins[0].startsWith(FastaManager.shuffle)) {
				genome = FastaManager.shuffle + (proteins[0].substring(FastaManager.shuffle.length()).split("_")[0]);
				isTarget = false;
			} else {
				genome = proteins[0].split("_")[0];
			}

			for (int j = 1; j < proteins.length; j++) {
				if (!proteins[j].startsWith(genome)) {
					continue L;
				}
			}

			if (!isTarget) {
				globalDecoyList.add(PEPi);
			}

			if (genomeTDMap.containsKey(genome)) {
				genomeTDMap.get(genome).add(PEPi);
			} else {
				ArrayList<Double> list = new ArrayList<Double>();
				list.add(PEPi);
				genomeTDMap.put(genome, list);
			}
		}

		double[] globalDecoys = new double[globalDecoyList.size()];
		for (int i = 0; i < globalDecoys.length; i++) {
			globalDecoys[i] = Math.log10(globalDecoyList.get(i));
		}
		DescriptiveStatistics decoyStats = new DescriptiveStatistics(globalDecoys);
		double decoyMean = decoyStats.getMean();
		double decoySd = decoyStats.getStandardDeviation();

		MannWhitneyUTest mwTest = new MannWhitneyUTest();

		HashSet<String> proset1 = new HashSet<String>();
		DecimalFormat df4 = FormatTool.getDF4();
		DecimalFormat dfe4 = FormatTool.getDFE4();

		int fileCount = fileSet.size();
		if (genomeScoreFile != null) {
		try (PrintWriter writer = new PrintWriter(genomeScoreFile)) {
			writer.println("Genome\tUnique_peptide_count\tDecoy_peptide_count\tTarget_mean\tDecoy_mean\tCohen's_d"
					+ "\tMann-Whitney_U\tMann-Whitney U p-value\tCliff's_delta");
			for (String genome : genomeTDMap.keySet()) {
				ArrayList<Double> targetList = genomeTDMap.get(genome);
				double[] targetLogInten = new double[targetList.size()];
				for (int i = 0; i < targetLogInten.length; i++) {
					targetLogInten[i] = Math.log10(targetList.get(i));
				}

				double[] usedLogInten;
				if (targetLogInten.length > globalDecoys.length) {
					Arrays.sort(targetLogInten);
					if (fileCount > globalDecoys.length) {
						usedLogInten = new double[fileCount];
					} else {
						usedLogInten = new double[globalDecoys.length];
					}
					for (int i = 0; i < globalDecoys.length; i++) {
						usedLogInten[i] = targetLogInten[targetLogInten.length - 1 - i];
					}
				} else {
					usedLogInten = targetLogInten;
				}

				if (usedLogInten.length > 0) {
					DescriptiveStatistics stats1 = new DescriptiveStatistics(usedLogInten);
					double mean1 = stats1.getMean();
					double sd1 = stats1.getStandardDeviation();
					double pooledSD = Math
							.sqrt(((stats1.getN() - 1) * sd1 * sd1 + (decoyStats.getN() - 1) * decoySd * decoySd)
									/ (stats1.getN() + decoyStats.getN() - 2));
					double cohensD = (mean1 - decoyMean) / pooledSD;
					double mwU = mwTest.mannWhitneyU(usedLogInten, globalDecoys);
					double mwPValue = mwTest.mannWhitneyUTest(usedLogInten, globalDecoys); // p-value

					if (mean1 > decoyMean) {
						writer.println(genome + "\t" + targetLogInten.length + "\t" + globalDecoys.length + "\t"
								+ df4.format(mean1) + "\t" + df4.format(decoyMean) + "\t" + df4.format(cohensD) + "\t"
								+ df4.format(mwU) + "\t" + dfe4.format(mwPValue) + "\t"
								+ df4.format(MathTool.calculateCliffsDelta(usedLogInten, globalDecoys)));

						if (mwPValue < 0.01) {
							genomeScoreMap.put(genome, new double[] { mwPValue, 0.0 });
						}
					}
				}
			}
			writer.close();
		} catch (Exception e) {
			logger.error(getTaskName() + ": error in calculating the genome score");
			System.out.println(
					format.format(new Date()) + "\t" + getTaskName() + ": error in calculating the genome score");
			return;
		}
		}

		String[] genomes = genomeScoreMap.keySet().toArray(String[]::new);
		Arrays.sort(genomes, Comparator.comparingDouble(genome -> genomeScoreMap.get(genome)[0]));
		int[] genomeTDCount = new int[2];
		for (String genome : genomes) {
			if (genome.startsWith(FastaManager.shuffle)) {
				genomeTDCount[0]++;
			} else {
				genomeTDCount[1]++;
			}
			double ratio = genomeTDCount[1] == 0 ? 1.0 : (double) genomeTDCount[0] / (double) genomeTDCount[1];
			if (ratio < 0.05) {
				genomeScoreMap.get(genome)[1] = ratio;
			} else {
				genomeScoreMap.remove(genome);
			}
		}

		HashMap<String, ArrayList<Double>> proteinTDMap = new HashMap<String, ArrayList<Double>>();
		for (int i = 0; i < pps.length; i++) {
			double PEPi = pps[i].getPreQuan();
			if (PEPi == 0) {
				continue;
			}

			String[] proteins = pps[i].getProGroups();
			for (int j = 0; j < proteins.length; j++) {
				String genome;
				if (proteins[j].startsWith(FastaManager.shuffle)) {
					genome = FastaManager.shuffle
							+ (proteins[j].substring(FastaManager.shuffle.length()).split("_")[0]);

				} else {
					genome = proteins[j].split("_")[0];
				}

				if (genomeScoreMap.containsKey(genome)) {
					if (proteinTDMap.containsKey(proteins[j])) {
						proteinTDMap.get(proteins[j]).add(PEPi);
					} else {
						ArrayList<Double> list = new ArrayList<Double>();
						list.add(PEPi);
						proteinTDMap.put(proteins[j], list);
					}
				}
			}
		}

		if (proteinScoreFile != null) {
		try (PrintWriter writer = new PrintWriter(proteinScoreFile)) {
			writer.println("Protein\tUnique_PSM_count\tDecoy_peptide_count\tTarget_mean\tDecoy_mean\tCohen's_d"
					+ "\tMann-Whitney_U\tMann-Whitney U p-value\tCliff's_delta");
			for (String protein : proteinTDMap.keySet()) {
				ArrayList<Double> targetList = proteinTDMap.get(protein);
				double[] targetLogInten = new double[targetList.size()];
				for (int i = 0; i < targetLogInten.length; i++) {
					targetLogInten[i] = Math.log10(targetList.get(i));
				}

				double[] usedLogInten;
				if (targetLogInten.length > fileCount) {
					Arrays.sort(targetLogInten);

					usedLogInten = new double[fileCount];
					for (int i = 0; i < fileCount; i++) {
						usedLogInten[i] = targetLogInten[targetLogInten.length - 1 - i];
					}
				} else {
					usedLogInten = targetLogInten;
				}

				if (usedLogInten.length > 0) {
					DescriptiveStatistics stats1 = new DescriptiveStatistics(targetLogInten);
					double mean1 = stats1.getMean();
					double sd1 = stats1.getStandardDeviation();
					double pooledSD = Math
							.sqrt(((stats1.getN() - 1) * sd1 * sd1 + (decoyStats.getN() - 1) * decoySd * decoySd)
									/ (stats1.getN() + decoyStats.getN() - 2));
					double cohensD = (mean1 - decoyMean) / pooledSD;
					double mwU = mwTest.mannWhitneyU(targetLogInten, globalDecoys);
					double mwPValue = mwTest.mannWhitneyUTest(targetLogInten, globalDecoys); // p-value
					proset1.add(protein);
					if (mean1 > decoyMean) {
						writer.println(protein + "\t" + targetLogInten.length + "\t" + globalDecoys.length + "\t"
								+ df4.format(mean1) + "\t" + df4.format(decoyMean) + "\t" + df4.format(cohensD) + "\t"
								+ df4.format(mwU) + "\t" + dfe4.format(mwPValue) + "\t"
								+ df4.format(MathTool.calculateCliffsDelta(targetLogInten, globalDecoys)));
						proteinScoreMap.put(protein, new double[] { mwPValue, 0.0 });
					}
				}
			}
			writer.close();
		} catch (Exception e) {
			logger.error(getTaskName() + ": error in calculating the protein score");
			System.out.println(
					format.format(new Date()) + "\t" + getTaskName() + ": error in calculating the protein score");
			return;
		}
		}

		String[] proteins = proteinScoreMap.keySet().toArray(String[]::new);
		Arrays.sort(proteins, Comparator.comparingDouble(pro -> proteinScoreMap.get(pro)[0]));
		int[] proTDCount = new int[2];
		for (String pro : proteins) {
			if (pro.startsWith(FastaManager.shuffle)) {
				proTDCount[0]++;
			} else {
				proTDCount[1]++;
			}
			double ratio = proTDCount[1] == 0 ? 1.0 : (double) proTDCount[0] / (double) proTDCount[1];
			proteinScoreMap.get(pro)[1] = ratio;
		}
	}

	protected void refineGenomeProteins(HashMap<String, HashSet<String>> genomePepMap,
			HashMap<String, HashSet<String>> proteinPepMap, int totalPeptideCount, HashSet<String> genomeSet,
			HashSet<String> proteinSet, File genomeRazorFile, File proteinRazorFile) {

		String[] genomes = genomePepMap.keySet().toArray(new String[genomePepMap.size()]);
		Arrays.sort(genomes, (g1, g2) -> {
			int s1 = genomePepMap.get(g1).size();
			int s2 = genomePepMap.get(g2).size();
			return s2 - s1;
		});

		int genomeCount = 0;
		for (int iteratorCount = 0; iteratorCount < 10; iteratorCount++) {
			HashSet<String> totalSet = new HashSet<String>();
			HashMap<String, Integer> orderMap = new HashMap<String, Integer>();
			int count = 0;
			for (int i = 0; i < genomes.length; i++) {
				int currentSize = totalSet.size();
				totalSet.addAll(genomePepMap.get(genomes[i]));
				int increaseSize = totalSet.size() - currentSize;
				orderMap.put(genomes[i], increaseSize);

				if (totalSet.size() == totalPeptideCount) {
					count = i;
				}
			}

			if (genomeCount == count) {
				break;
			} else {
				Arrays.sort(genomes, (g1, g2) -> {
					if (orderMap.get(g2) == orderMap.get(g1)) {
						return genomePepMap.get(g2).size() - genomePepMap.get(g1).size();
					} else {
						return orderMap.get(g2) - orderMap.get(g1);
					}
				});

				genomeCount = count;

				ArrayList<String> list = new ArrayList<String>();
				for (int i = 0; i < genomes.length; i++) {
					if (orderMap.get(genomes[i]) == 0) {
						break;
					} else {
						list.add(genomes[i]);
					}
				}
				genomes = list.toArray(String[]::new);
			}
		}

		HashSet<String> totalPepSet = new HashSet<String>();
		try (PrintWriter writer = new PrintWriter(genomeRazorFile)) {
			writer.println(
					"Rank\tGenome\tPeptide_count\tTotal_covered_peptide_count\tPercentage\tIncrement\tIncrement ratio");

			for (int i = 0; i < genomes.length; i++) {
				int currentSize = totalPepSet.size();
				totalPepSet.addAll(genomePepMap.get(genomes[i]));
				int add = totalPepSet.size() - currentSize;
				double additional = ((double) add / (double) totalPeptideCount);
				double percentage = ((double) totalPepSet.size() / (double) totalPeptideCount);

				writer.println(i + "\t" + genomes[i] + "\t" + genomePepMap.get(genomes[i]).size() + "\t"
						+ totalPepSet.size() + "\t" + percentage + "\t" + add + "\t" + additional);
				genomeSet.add(genomes[i]);
			}
			writer.close();
		} catch (IOException e) {
			logger.error(getTaskName() + ": error in writing genomes to " + genomeRazorFile, e);
			System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": error in writing genomes to "
					+ genomeRazorFile);
		}

		ArrayList<String> proList = new ArrayList<String>();
		for (String pro : proteinPepMap.keySet()) {
			String genome;
			if (pro.startsWith(FastaManager.shuffle)) {
				genome = FastaManager.shuffle + pro.substring(FastaManager.shuffle.length()).split("_")[0];
			} else {
				genome = pro.split("_")[0];
			}
			if (genomeSet.contains(genome)) {
				proList.add(pro);
			}
		}

		String[] proteins = proList.toArray(new String[proList.size()]);
		Arrays.sort(proteins, (g1, g2) -> {
			int s1 = proteinPepMap.get(g1).size();
			int s2 = proteinPepMap.get(g2).size();
			return s2 - s1;
		});

		int proteinCount = 0;
		for (int iteratorCount = 0; iteratorCount < 10; iteratorCount++) {
			HashSet<String> totalSet = new HashSet<String>();
			HashMap<String, Integer> orderMap = new HashMap<String, Integer>();
			int count = 0;
			for (int i = 0; i < proteins.length; i++) {
				int currentSize = totalSet.size();
				totalSet.addAll(proteinPepMap.get(proteins[i]));
				int increaseSize = totalSet.size() - currentSize;
				orderMap.put(proteins[i], increaseSize);

				if (totalSet.size() == totalPeptideCount) {
					count = i;
				}
			}

			if (proteinCount == count) {
				break;
			} else {
				Arrays.sort(proteins, (g1, g2) -> {
					if (orderMap.get(g2) == orderMap.get(g1)) {
						return proteinPepMap.get(g2).size() - proteinPepMap.get(g1).size();
					} else {
						return orderMap.get(g2) - orderMap.get(g1);
					}
				});
				proteinCount = count;

				ArrayList<String> list = new ArrayList<String>();
				for (int i = 0; i < proteins.length; i++) {
					if (orderMap.get(proteins[i]) == 0) {
						break;
					} else {
						list.add(proteins[i]);
					}
				}
				proteins = list.toArray(String[]::new);
			}
		}

		totalPepSet = new HashSet<String>();
		try (PrintWriter writer = new PrintWriter(proteinRazorFile)) {
			writer.println(
					"Rank\tProtein\tPeptide_count\tTotal_covered_peptide_count\tPercentage\tIncrement\tIncrement ratio");

			for (int i = 0; i < proteins.length; i++) {
				int currentSize = totalPepSet.size();
				totalPepSet.addAll(proteinPepMap.get(proteins[i]));
				int add = totalPepSet.size() - currentSize;
				double additional = ((double) add / (double) totalPeptideCount);
				double percentage = ((double) totalPepSet.size() / (double) totalPeptideCount);

				writer.println(i + "\t" + proteins[i] + "\t" + proteinPepMap.get(proteins[i]).size() + "\t"
						+ totalPepSet.size() + "\t" + percentage + "\t" + add + "\t" + additional);
				proteinSet.add(proteins[i]);
			}
			writer.close();
		} catch (IOException e) {
			logger.error(getTaskName() + ": error in writing genomes to " + proteinRazorFile, e);
			System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": error in writing genomes to "
					+ proteinRazorFile);
		}
	}

	@SuppressWarnings("unchecked")
	protected boolean exportReport() {
		return doExportReportInternal();
	}

	@SuppressWarnings("unchecked")
	protected boolean exportReportLegacy() {

		logger.info(getTaskName() + ": protein inference started");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": protein inference started");

		HashMap<String, HashSet<String>> hostPepProMap = new HashMap<String, HashSet<String>>();
		File hostFile = new File(resultFolderFile, "host");
		File hostTsvFile = new File(hostFile, "host.tsv");
		if (hostTsvFile.exists()) {
			File magDbFile = this.magDb.getCurrentFile();
			File libFolderFile = new File(magDbFile, "libraries");
			File hostParquetFile = new File(new File(libFolderFile, "host"), "host.parquet");

			if (hostParquetFile.exists()) {
				try {
					hostPepProMap = DiaNNParquetReader.getPepProMap(hostParquetFile.getAbsolutePath());
				} catch (IOException e) {
					// TODO Auto-generated catch block
					logger.error(getTaskName() + ": error in reading the host database from " + hostParquetFile);
					System.err.println(format.format(new Date()) + "\t" + getTaskName()
							+ ": error in reading the host database from " + hostParquetFile);
				}
			} else {
				logger.error(getTaskName() + ": host database " + hostParquetFile + " doesn't exist");
				System.err.println(format.format(new Date()) + "\t" + getTaskName() + ": host database "
						+ hostParquetFile + " doesn't exist");

				return false;
			}
		}

		MetaParameterDia diaPar = (MetaParameterDia) metaPar;
		String[] usedLibs = diaPar.getLibrary();

		File combinedFile = new File(metaPar.getDbSearchResultFile(), "combined");
		HashMap<String, HashMap<Integer, String>> fastaProPepMap = new HashMap<String, HashMap<Integer, String>>();

		File combinedFastaFile = new File(combinedFile, "combined.fasta");
		if (combinedFastaFile.exists()) {

			try (BufferedReader reader = new BufferedReader(new FileReader(combinedFastaFile))) {
				String line = null;
				while ((line = reader.readLine()) != null) {
					if (line.startsWith(">")) {
						String proName = line.split("\\s")[1];
						String sequence = reader.readLine();
						if (fastaProPepMap.containsKey(proName)) {
							fastaProPepMap.get(proName).put(sequence.length(), sequence);
						} else {
							HashMap<Integer, String> seqLengthMap = new HashMap<Integer, String>();
							seqLengthMap.put(sequence.length(), sequence);
							fastaProPepMap.put(proName, seqLengthMap);
						}
					}
				}
				reader.close();
			} catch (IOException e) {
				logger.error(getTaskName() + ": error in reading the library from " + combinedFastaFile, e);
				System.err.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": error in reading the library from " + combinedFastaFile);

				return false;
			}

		} else {

			File targetFastaFile = new File(metaPar.getDbSearchResultFile(),
					metaPar.getDbSearchResultFile().getName() + ".fasta");
			if (targetFastaFile.exists()) {
				try (BufferedReader reader = new BufferedReader(new FileReader(targetFastaFile))) {
					String line = null;
					while ((line = reader.readLine()) != null) {
						if (line.startsWith(">")) {
							String proName = line.split("\\s")[1];
							String sequence = reader.readLine();
							if (fastaProPepMap.containsKey(proName)) {
								fastaProPepMap.get(proName).put(sequence.length(), sequence);
							} else {
								HashMap<Integer, String> seqLengthMap = new HashMap<Integer, String>();
								seqLengthMap.put(sequence.length(), sequence);
								fastaProPepMap.put(proName, seqLengthMap);
							}
						}
					}
					reader.close();
				} catch (IOException e) {
					logger.error(getTaskName() + ": error in reading the library from " + targetFastaFile, e);
					System.err.println(format.format(new Date()) + "\t" + getTaskName()
							+ ": error in reading the library from " + targetFastaFile);

					return false;
				}
			} else {
				File magDbFile = this.magDb.getCurrentFile();
				File libFolderFile = new File(magDbFile, "libraries");
				for (int i = 0; i < usedLibs.length; i++) {
					File libFastaFile = new File(new File(libFolderFile, usedLibs[i]), usedLibs[i] + ".fasta");
					if (libFastaFile.exists()) {
						try (BufferedReader reader = new BufferedReader(new FileReader(libFastaFile))) {
							String line = null;
							while ((line = reader.readLine()) != null) {
								if (line.startsWith(">")) {
									String proName = line.split("\\s")[1];
									String sequence = reader.readLine();
									if (fastaProPepMap.containsKey(proName)) {
										fastaProPepMap.get(proName).put(sequence.length(), sequence);
									} else {
										HashMap<Integer, String> seqLengthMap = new HashMap<Integer, String>();
										seqLengthMap.put(sequence.length(), sequence);
										fastaProPepMap.put(proName, seqLengthMap);
									}
								}
							}
							reader.close();
						} catch (IOException e) {
							logger.error(getTaskName() + ": error in reading the library from " + libFastaFile, e);
							System.err.println(format.format(new Date()) + "\t" + getTaskName()
									+ ": error in reading the library from " + libFastaFile);

							return false;
						}
					} else {
						File libParquetFile = new File(new File(libFolderFile, usedLibs[i]), usedLibs[i] + ".parquet");
						if (libParquetFile.exists()) {
							try {
								HashMap<String, HashSet<String>> parquetProPepMap = DiaNNParquetReader
										.getProPepMap(libParquetFile.getAbsolutePath());

								for (String proName : parquetProPepMap.keySet()) {
									HashSet<String> seqSet = parquetProPepMap.get(proName);
									for (String sequence : seqSet) {
										if (fastaProPepMap.containsKey(proName)) {
											fastaProPepMap.get(proName).put(sequence.length(), sequence);
										} else {
											HashMap<Integer, String> seqLengthMap = new HashMap<Integer, String>();
											seqLengthMap.put(sequence.length(), sequence);
											fastaProPepMap.put(proName, seqLengthMap);
										}
									}
								}

							} catch (IOException e) {
								// TODO Auto-generated catch block
								logger.error(getTaskName() + ": error in reading the protein database from "
										+ libParquetFile);
								System.err.println(format.format(new Date()) + "\t" + getTaskName()
										+ ": error in reading the protein database from " + libParquetFile);
							}
						}

					}
				}
			}
		}

		File parquetFile = new File(combinedFile, "combined.parquet");
		DiaNNParquetReader reader = new DiaNNParquetReader(parquetFile.getAbsolutePath(), this.diannPar.getCut());
		this.diaPp = reader.getDiaNNPrecursors(0.01);

		logger.info(getTaskName() + ": " + this.diaPp.length + " PSMs were obtained with Global.Q.Value<0.01");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": " + this.diaPp.length
				+ " PSMs were obtained with Global.Q.Value<0.01");

		int originalDecoyCount = 0;
		HashSet<String> pepSet = new HashSet<String>();
		HashMap<String, String> decoyPepMap = new HashMap<String, String>();
		if (diaPp != null && diaPp.length > 0) {
			for (int i = 0; i < diaPp.length; i++) {
				String sequence = diaPp[i].getSeqString();
				pepSet.add(sequence);
				if (diaPp[i].getIsDecoy()) {
					String[] decoyPros = diaPp[i].getProGroups();
					String targetSeq = "";
					int lengthDiff = sequence.length();
					for (int j = 0; j < decoyPros.length; j++) {
						if (fastaProPepMap.containsKey(decoyPros[j])) {
							HashMap<Integer, String> seqLengthMap = fastaProPepMap.get(decoyPros[j]);
							if (seqLengthMap.containsKey(sequence.length())) {
								targetSeq = seqLengthMap.get(sequence.length());
							} else {
								for (Integer length : seqLengthMap.keySet()) {
									if (Math.abs(sequence.length() - length) < lengthDiff) {
										lengthDiff = Math.abs(sequence.length() - length);
										targetSeq = seqLengthMap.get(length);
									}
								}
							}
						}
					}
					if (targetSeq.length() > 0) {
						decoyPepMap.put(sequence, targetSeq);
						pepSet.add(targetSeq);
						originalDecoyCount++;
					}
				}
			}
		}

		logger.info(getTaskName() + ": total sequence count: " + pepSet.size() + "; decoy precursor count: "
				+ originalDecoyCount + "; decoy sequence count: " + decoyPepMap.size());
		System.out
				.println(format.format(new Date()) + "\t" + getTaskName() + ": total sequence count: " + pepSet.size());

		ConcurrentHashMap<String, HashMap<String, HashSet<String>>> matchedMap = getGenomeFromPep(pepSet, threadCount);
		if (matchedMap == null) {
			return false;
		}

		HashMap<String, HashSet<String>> pepProMap = new HashMap<String, HashSet<String>>();
		for (String genome : matchedMap.keySet()) {
			HashMap<String, HashSet<String>> proPepMap = matchedMap.get(genome);
			for (String pro : proPepMap.keySet()) {
				HashSet<String> pset = proPepMap.get(pro);
				for (String pep : pset) {
					if (pepProMap.containsKey(pep)) {
						pepProMap.get(pep).add(pro);
					} else {
						HashSet<String> proSet = new HashSet<String>();
						proSet.add(pro);
						pepProMap.put(pep, proSet);
					}
				}
			}
		}

		logger.info(getTaskName() + ": matched sequence count: " + pepProMap.size());
		System.out.println(
				format.format(new Date()) + "\t" + getTaskName() + ": matched sequence count: " + pepProMap.size());

		int decoyPPCount = 0;
		HashMap<String, HashSet<String>> genomeMap = new HashMap<String, HashSet<String>>();
		HashMap<String, HashSet<String>> proteinPepMap = new HashMap<String, HashSet<String>>();
		HashSet<String> tempPepSet = new HashSet<String>();
		HashSet<String> missSet = new HashSet<String>();
		for (DiaNNPrecursor pp : diaPp) {
			String sequence = pp.getSeqString();
			if (pepProMap.containsKey(sequence)) {
				for (String pro : pepProMap.get(sequence)) {
					String genome = pro.split("_")[0];
					if (genomeMap.containsKey(genome)) {
						genomeMap.get(genome).add(sequence);
					} else {
						HashSet<String> set = new HashSet<String>();
						set.add(sequence);
						genomeMap.put(genome, set);
					}
					if (proteinPepMap.containsKey(pro)) {
						proteinPepMap.get(pro).add(sequence);
					} else {
						HashSet<String> set = new HashSet<String>();
						set.add(sequence);
						proteinPepMap.put(pro, set);
					}
				}
				tempPepSet.add(sequence);
			} else if (hostPepProMap.containsKey(sequence)) {
				String[] pros = hostPepProMap.get(sequence).toArray(String[]::new);
				pp.setProGroups(pros);
			} else {
				if (decoyPepMap.containsKey(sequence)) {
					String targetSequence = decoyPepMap.get(sequence);
					if (pepProMap.containsKey(targetSequence)) {
						for (String pro : pepProMap.get(targetSequence)) {
							String genome = pro.split("_")[0];
							String decoyGenome = FastaManager.shuffle + genome;
							String decoyPro = FastaManager.shuffle + pro;

							if (genomeMap.containsKey(decoyGenome)) {
								genomeMap.get(decoyGenome).add(sequence);
							} else {
								HashSet<String> set = new HashSet<String>();
								set.add(sequence);
								genomeMap.put(decoyGenome, set);
							}

							if (proteinPepMap.containsKey(decoyPro)) {
								proteinPepMap.get(decoyPro).add(sequence);
							} else {
								HashSet<String> set = new HashSet<String>();
								set.add(sequence);
								proteinPepMap.put(decoyPro, set);
							}
						}
						tempPepSet.add(sequence);
						decoyPPCount++;
					}
				} else {
					missSet.add(sequence);
				}
			}
		}

		logger.info(getTaskName() + ": " + missSet.size()
				+ " peptide sequences were not found from the in-silico digested MAG database, "
				+ "maybe due to unspecific digestion, will be not be included in the following analysis");
		logger.info(getTaskName() + ": " + decoyPPCount
				+ " decoy precursor will be used for FDR control in the following steps.");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": " + missSet.size()
				+ " peptide sequences were not found from the in-silico digested MAG database, "
				+ "maybe due to unspecific digestion, will be not be included in the following analysis");

		File genomeRazorFile = new File(combinedFile, "genomes_razor.tsv");
		File proRazorFile = new File(combinedFile, "proteins_razor.tsv");

		HashSet<String> genomeSet = new HashSet<String>();
		HashSet<String> proteinSet = new HashSet<String>();

		refineGenomeProteins(genomeMap, proteinPepMap, tempPepSet.size(), genomeSet, proteinSet, genomeRazorFile,
				proRazorFile);

		ArrayList<DiaNNPrecursor> genomeFilteredList = new ArrayList<DiaNNPrecursor>();
		HashSet<String> temppepHashSet1 = new HashSet<String>();
		for (DiaNNPrecursor pp : diaPp) {
			String sequence = pp.getSeqString();
			if (hostPepProMap.containsKey(sequence)) {
				genomeFilteredList.add(pp);
			} else {
				if (tempPepSet.contains(sequence)) {
					ArrayList<String> proList = new ArrayList<String>();
					if (pepProMap.containsKey(sequence)) {
						HashSet<String> pros = pepProMap.get(sequence);
						for (String pro : pros) {
							String genome = pro.split("_")[0];
							if (genomeSet.contains(genome) && proteinSet.contains(pro)) {
								proList.add(pro);
							}
						}
					} else {
						if (decoyPepMap.containsKey(sequence)) {
							String targetSeq = decoyPepMap.get(sequence);
							if (pepProMap.containsKey(targetSeq)) {
								HashSet<String> pros = pepProMap.get(targetSeq);
								for (String pro : pros) {
									String genome = pro.split("_")[0];
									if (genomeSet.contains(FastaManager.shuffle + genome)
											&& proteinSet.contains(FastaManager.shuffle + pro)) {
										proList.add(FastaManager.shuffle + pro);
									}
								}
							}
						}
					}

					if (proList.size() > 0) {
						pp.setProGroups(proList.toArray(String[]::new));
						genomeFilteredList.add(pp);
						temppepHashSet1.add(sequence);
					}
				}
			}
		}

		DiaNNPrecursor[] filteredPPS = genomeFilteredList.toArray(DiaNNPrecursor[]::new);

		HashMap<String, double[]> genomeScoreMap = new HashMap<String, double[]>();
		HashMap<String, double[]> proteinScoreMap = new HashMap<String, double[]>();
		calGenomeProteinScore(filteredPPS, genomeScoreMap, proteinScoreMap, null, null);

		HashSet<String> tempPepHashSet = new HashSet<String>();
		ArrayList<DiaNNPrecursor> finalFilteredList = new ArrayList<DiaNNPrecursor>();
		for (DiaNNPrecursor pp : filteredPPS) {
			ArrayList<String> proList = new ArrayList<String>();
			String[] pros = pp.getProGroups();
			String razorProtein = "";
			double razorQvalue = 1.0;
			for (String pro : pros) {
				if (proteinScoreMap.containsKey(pro)) {
					proList.add(pro);
					if (proteinScoreMap.get(pro)[0] < razorQvalue) {
						razorQvalue = proteinScoreMap.get(pro)[0];
						razorProtein = pro;
					}
				}
			}

			if (proList.size() > 0) {
				tempPepHashSet.add(pp.getSeqString());
				pp.setProGroups(proList.toArray(String[]::new));
				pp.setRazorProtein(razorProtein);
				finalFilteredList.add(pp);
			}
		}

		this.diaPp = finalFilteredList.toArray(DiaNNPrecursor[]::new);
		DiannPrMatrixReader prReader = new DiannPrMatrixReader(this.quan_pep_file);
		FlashLfqQuanPeptide[] peps = prReader.getMetaPeptides(this.diaPp);
		this.quanExpNames = prReader.getIntensityTitle();

		HashMap<String, FlashLfqQuanPeptide> modseqPepMap = new HashMap<String, FlashLfqQuanPeptide>();
		for (FlashLfqQuanPeptide pep : peps) {
			modseqPepMap.put(pep.getModSeq(), pep);
		}

		logger.info(getTaskName() + ": protein inference finished, " + this.diaPp.length + " PSMs, "
				+ modseqPepMap.size() + " peptide, " + proteinScoreMap.size() + " proteins and " + genomeScoreMap.size()
				+ " genomes were kept");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": protein inference finished, "
				+ this.diaPp.length + " PSMs, " + modseqPepMap.size() + " peptide, " + proteinScoreMap.size()
				+ " proteins and " + genomeScoreMap.size() + " genomes were kept");

		logger.info(getTaskName() + ": exporting peptide report started");
		logger.info("[{}] exporting peptide report started", getTaskName());

		updateProgress(82, 100);

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

			if (razorProtein.startsWith(magDb.getIdentifier())) {
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

		int[] filePsmCount = new int[quanExpNames.length];
		HashSet<String>[] fileSequenceSets = new HashSet[quanExpNames.length];
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

						if (proteins[j].startsWith(magDb.getIdentifier())) {
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

		String[] proteins = proteinTotalIntenMap.keySet().toArray(String[]::new);
		Arrays.sort(proteins, new Comparator<String>() {

			@Override
			public int compare(String o1, String o2) {
				// TODO Auto-generated method stub
				if (proteinScoreMap.get(o1)[0] > proteinScoreMap.get(o2)[0]) {
					return 1;
				} else if (proteinScoreMap.get(o1)[0] < proteinScoreMap.get(o2)[0]) {
					return -1;
				} else {
					if (proteinTotalIntenMap.get(o1) > proteinTotalIntenMap.get(o2)) {
						return -1;
					} else if (proteinTotalIntenMap.get(o1) < proteinTotalIntenMap.get(o2)) {
						return 1;
					} else {
						return 0;
					}
				}
			}
		});

		HashMap<String, Integer> proteinIdMap = new HashMap<String, Integer>();
		for (int i = 0; i < proteins.length; i++) {
			proteinIdMap.put(proteins[i], i);
		}

		try (PrintWriter idWriter = new PrintWriter(new File(resultFolderFile, "id.tsv"))) {
			idWriter.println("Raw file\tUnique peptide count\tPSM count");
			for (int i = 0; i < quanExpNames.length; i++) {
				idWriter.println(quanExpNames[i] + "\t" + fileSequenceSets[i].size() + "\t" + filePsmCount[i]);
			}
			idWriter.close();
		} catch (IOException e) {
			System.out.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": error in writing identification information to " + new File(resultFolderFile, "id.tsv"));
			logger.error(getTaskName() + ": error in writing identification information to "
					+ new File(resultFolderFile, "id.tsv"), e);
		}

		updateProgress(85, 100);

		HashMap<String, ArrayList<Integer>> proPepIdMap = new HashMap<String, ArrayList<Integer>>();
		HashMap<String, ArrayList<Integer>> proPepRazorIdMap = new HashMap<String, ArrayList<Integer>>();
		try {
			this.final_pep_txt = new File(resultFolderFile, "final_peptides.tsv");

			PrintWriter writer = new PrintWriter(this.final_pep_txt);

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

			for (int i = 0; i < this.quanExpNames.length; i++) {
				titlesb.append("Identification type " + this.quanExpNames[i]).append("\t");
			}

			for (int i = 0; i < this.quanExpNames.length; i++) {
				titlesb.append("MS/MS count " + this.quanExpNames[i]).append("\t");
			}

			titlesb.append("Intensity").append("\t");

			for (int i = 0; i < quanExpNames.length; i++) {
				titlesb.append("Intensity " + this.quanExpNames[i]).append("\t");
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
			logger.error(getTaskName() + ": error in writing the peptides to  " + final_pep_txt, e);
			System.err.println(format.format(new Date()) + "\t" + getTaskName() + ": error in writing the peptides to  "
					+ final_pep_txt);
		}

		logger.info(getTaskName() + ": peptide report has been exported to " + final_pep_txt);
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": peptide report has been exported to "
				+ final_pep_txt);

		if (!final_pep_txt.exists() || final_pep_txt.length() == 0) {
			return false;
		}

		updateProgress(87, 100);

		logger.info(getTaskName() + ": exporting protein report started");
		logger.info("[{}] exporting protein report started", getTaskName());

		ArrayList<MetaProtein> proteinList = new ArrayList<MetaProtein>();
		try {

			this.final_pro_txt = new File(resultFolderFile, "final_proteins.tsv");

			PrintWriter writer = new PrintWriter(this.final_pro_txt);

			StringBuilder titlesb = new StringBuilder();
			titlesb.append("Protein IDs").append("\t");
			titlesb.append("Peptide counts (all)").append("\t");
			titlesb.append("Peptide counts (razor)").append("\t");
			titlesb.append("Number of PSMs").append("\t");
			titlesb.append("Qvalue").append("\t");
			titlesb.append("Local FDR").append("\t");
			titlesb.append("Intensity").append("\t");
			for (int i = 0; i < quanExpNames.length; i++) {
				titlesb.append("Intensity " + quanExpNames[i]).append("\t");
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

				sb.append(df.format(proteinScoreMap.get(proteins[i])[0])).append("\t");
				sb.append(df4.format(proteinScoreMap.get(proteins[i])[1])).append("\t");

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

				MetaProtein mp = new MetaProtein(i, 0, proteins[i], totalPepCountMap.get(proteins[i]).size(),
						psmCountMap.get(proteins[i]), proteinScoreMap.get(proteins[i])[0], intensity);
				proteinList.add(mp);
			}

			writer.close();

		} catch (IOException e) {
			// TODO Auto-generated catch block
			logger.error(getTaskName() + ": error in writing final protein result to " + this.final_pro_txt.getName(),
					e);
			System.err.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": error in writing final protein result to " + this.final_pro_txt.getName());
		}

		logger.info(getTaskName() + ": protein report has been exported to " + final_pro_txt.getName());
		logger.info("[{}] protein report has been exported to {}", getTaskName(), final_pro_txt.getName());

		updateProgress(89, 100);

		logger.info(getTaskName() + ": functional annotation started");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": functional annotation started");

		File funcFile = new File(resultFolderFile, "functional_annotation");
		if (!funcFile.exists()) {
			funcFile.mkdir();
		}

		File funTsv = new File(funcFile, "functions.tsv");
		final_pro_xml_file = new File(funcFile, "functions.xml");

		HashMap<String, Integer> proTaxIdMap = new HashMap<String, Integer>();
		MetaProtein[] metapros = proteinList.toArray(new MetaProtein[proteinList.size()]);
		MetaProteinAnnoMag[] mpas = funcAnnoSqlite(metapros, quanExpNames);
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
			funReader = new MetaProteinXMLReader2(final_pro_xml_file);
		} catch (DocumentException e) {
			// TODO Auto-generated catch block
			logger.error(getTaskName() + ": error in reading protein function information from "
					+ this.final_pro_xml_file.getName(), e);
			System.err.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": error in reading protein function information from " + this.final_pro_xml_file.getName());
			return false;
		}

		if (!funTsv.exists()) {
			System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": writing results to "
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
				logger.error(getTaskName() + ": error in writing protein function information to " + funTsv.getName(),
						e);
				System.err.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": error in writing protein function information to " + funTsv.getName());
			}
			System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": writing results to "
					+ funTsv.getName() + " finished");
		}

		File html = new File(funcFile, "functions.html");
		if (!html.exists()) {
			System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": writing results to "
					+ html.getName() + " started");
			try {
				funReader.exportHtml(html);

				// Copy the generated functions.html to the reports directory
				File reportFunctionsHtml = new File(metaPar.getReportDir(), "reports/functions.html");
				if (reportFunctionsHtml.getParentFile() != null && !reportFunctionsHtml.getParentFile().exists()) {
					reportFunctionsHtml.getParentFile().mkdirs();
				}
				Files.copy(html.toPath(), reportFunctionsHtml.toPath(), StandardCopyOption.REPLACE_EXISTING);
				logger.info("Copied functions.html to reports directory: " + reportFunctionsHtml.getAbsolutePath());

			} catch (IOException e) {
				// TODO Auto-generated catch block
				logger.error(getTaskName() + ": error in writing protein function information to " + html.getName(), e);
				System.err.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": error in writing protein function information to " + html.getName());
			}
			System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": writing results to "
					+ html.getName() + " finished");
		}

		// Copy functions.html to report directory for navigation access
		if (html.exists() && this.reportHtmlDir != null) {
			File reportFunctionsHtml = new File(this.reportHtmlDir, "functions.html");
			if (!reportFunctionsHtml.exists()) {
				try {
					java.nio.file.Files.copy(html.toPath(), reportFunctionsHtml.toPath(),
							java.nio.file.StandardCopyOption.REPLACE_EXISTING);
					logger.info(getTaskName() + ": copied functions.html to report directory");
					System.out.println(format.format(new Date()) + "\t" + getTaskName()
							+ ": copied functions.html to report directory");
				} catch (IOException e) {
					logger.error(getTaskName() + ": error copying functions.html to report directory", e);
					System.err.println(format.format(new Date()) + "\t" + getTaskName()
							+ ": error copying functions.html to report directory");
				}
			}
		}

		updateProgress(91, 100);

		// taxonomy analysis

		logger.info(getTaskName() + ": taxonomy analysis started");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": taxonomy analysis started");

		File taxFile = new File(resultFolderFile, "taxonomy_analysis");
		if (!taxFile.exists()) {
			taxFile.mkdir();
		}

		MagFuncSearcher searcher = null;
		try {
			searcher = new MagFuncSearcher(metaPar.getMicroDb(), mrp.getToolPath(ToolType.FUNCDB), magDb);

		} catch (NumberFormatException | SQLException | IOException e) {
			// TODO Auto-generated catch block
			logger.error(getTaskName() + ": error in genome analysis, task failed", e);
			System.err.println(
					format.format(new Date()) + "\t" + getTaskName() + ": error in genome analysis, task failed");

			return false;
		}

		String[] genomes = genomeIntensityMap.keySet().toArray(new String[genomeIntensityMap.size()]);
		HashMap<String, String[]> genomeTaxaMap = null;
		try {
			genomeTaxaMap = searcher.matchGenome(genomes);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			logger.error(getTaskName() + ": Error in searching the genome2taxa database", e);
			System.err.println(
					format.format(new Date()) + "\t" + getTaskName() + ": error in searching the genome2taxa database");
			return false;
		}

		File genomeHtmlFile = new File(taxFile, "Genome.html");
		try {
			HMGenomeHtmlWriter htmlWriter = new HMGenomeHtmlWriter(genomeHtmlFile);
			htmlWriter.write(genomeTaxaMap);
			htmlWriter.close();

			// Copy the generated Genome.html to the reports directory
			File reportGenomeHtml = new File(metaPar.getReportDir(), "reports/Genome.html");
			if (reportGenomeHtml.getParentFile() != null && !reportGenomeHtml.getParentFile().exists()) {
				reportGenomeHtml.getParentFile().mkdirs();
			}
			Files.copy(genomeHtmlFile.toPath(), reportGenomeHtml.toPath(), StandardCopyOption.REPLACE_EXISTING);
			logger.info("Copied Genome.html to reports directory: " + reportGenomeHtml.getAbsolutePath());
		} catch (IOException e) {
			// TODO Auto-generated catch block
			logger.error("Error in writing genome information to " + genomeHtmlFile, e);
			System.err.println(
					format.format(new Date()) + "\t" + "Error in writing genome information to " + genomeHtmlFile);
			return false;
		}

		// Copy Genome.html to report directory for navigation access
		if (genomeHtmlFile.exists() && this.reportHtmlDir != null) {
			File reportGenomeHtml = new File(this.reportHtmlDir, "Genome.html");
			if (!reportGenomeHtml.exists()) {
				try {
					java.nio.file.Files.copy(genomeHtmlFile.toPath(), reportGenomeHtml.toPath(),
							java.nio.file.StandardCopyOption.REPLACE_EXISTING);
					logger.info(getTaskName() + ": copied Genome.html to report directory");
					System.out.println(format.format(new Date()) + "\t" + getTaskName()
							+ ": copied Genome.html to report directory");
				} catch (IOException e) {
					logger.error(getTaskName() + ": error copying Genome.html to report directory", e);
					System.err.println(format.format(new Date()) + "\t" + getTaskName()
							+ ": error copying Genome.html to report directory");
				}
			}
		}

		TaxonomyRanks[] mainRanks = TaxonomyRanks.getMainRanks();
		double[] cellularIntensity = new double[quanExpNames.length];
		HashMap<String, String[]>[] taxRankMaps = new HashMap[mainRanks.length];
		HashMap<String, double[]>[] taxIntensityMaps = new HashMap[mainRanks.length];
		HashMap<String, Double>[] taxScoreMaps = new HashMap[mainRanks.length];
		for (int i = 0; i < taxRankMaps.length; i++) {
			taxRankMaps[i] = new HashMap<String, String[]>();
			taxIntensityMaps[i] = new HashMap<String, double[]>();
			taxScoreMaps[i] = new HashMap<String, Double>();
		}

		// Core/mid/tail from ranked.tsv (same as DDA: prefer combined result)
		HashMap<String, String> genomeSelectionLabel = new HashMap<>();
		MetaDiaMarkerTask.ensureCombinedRankedFromPrecursors(this, resultFolderFile);
		File rankedFile = new File(resultFolderFile, "combined/ranked.tsv");
		if (!rankedFile.exists()) rankedFile = new File(resultFolderFile, "marker/genome_ranked.tsv");
		if (!rankedFile.exists()) rankedFile = new File(resultFolderFile, "seed/genome_ranked.tsv");
		if (rankedFile.exists()) {
			MetaDiaMarkerTask.readCoreMidTailFromRankedTsv(rankedFile, genomeSelectionLabel, genomeIntensityMap.keySet());
			MetaPilotMagTask.RankedGenomeSelectionResult rankedAll =
					MetaPilotMagTask.buildExportSelectionFromRankedTsv(rankedFile, null);
			int contributingAll = rankedAll.orderedGenomes.size();
			int labeledWithIntensity = genomeSelectionLabel.size();
			int coreCount = 0, midCount = 0, tailCount = 0;
			for (String lbl : genomeSelectionLabel.values()) {
				if ("core".equals(lbl)) coreCount++;
				else if ("mid".equals(lbl)) midCount++;
				else tailCount++;
			}
			java.util.HashSet<String> contributingSet = new java.util.HashSet<>(rankedAll.orderedGenomes);
			java.util.HashSet<String> intensitySet = new java.util.HashSet<>(genomeIntensityMap.keySet());
			java.util.HashSet<String> inWindowNoIntensity = new java.util.HashSet<>(contributingSet);
			inWindowNoIntensity.removeAll(intensitySet);
			java.util.HashSet<String> intensityBeyondWindow = new java.util.HashSet<>(intensitySet);
			intensityBeyondWindow.removeAll(contributingSet);
			logger.info("{}: Genome export summary – ranked_window={}, labeled_with_intensity={}, core={}, mid={}, tail={}, window_no_intensity={}, intensity_beyond_window={}",
					getTaskName(), Integer.valueOf(contributingAll), Integer.valueOf(labeledWithIntensity),
					Integer.valueOf(coreCount), Integer.valueOf(midCount), Integer.valueOf(tailCount),
					Integer.valueOf(inWindowNoIntensity.size()), Integer.valueOf(intensityBeyondWindow.size()));
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

		File genomeFile = new File(taxFile, "Genome.tsv");
		PrintWriter genomeWriter = null;
		try {
			genomeWriter = new PrintWriter(genomeFile);
			StringBuilder genomeTitlesb = new StringBuilder();
			genomeTitlesb.append("Genome\t");
			genomeTitlesb.append("p-value\t");
			genomeTitlesb.append("Local FDR\t");
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

			for (String exp : expNames) {
				if (exp.startsWith("Intensity ")) {
					genomeTitlesb.append("\t").append(exp);
				} else {
					genomeTitlesb.append("\tIntensity ").append(exp);
				}
			}

			genomeWriter.println(genomeTitlesb);

			updateProgress(94, 100);

			for (int i = 0; i < genomes.length; i++) {
				String[] taxaContent = genomeTaxaMap.get(genomes[i]);
				if (genomeIntensityMap.containsKey(genomes[i])) {
					double[] genomeIntensity = genomeIntensityMap.get(genomes[i]);
					for (int j = 0; j < genomeIntensity.length; j++) {
						cellularIntensity[j] += genomeIntensity[j];
					}

					StringBuilder genomeTaxaSb = new StringBuilder();
					genomeTaxaSb.append(genomes[i]).append("\t");

					genomeTaxaSb.append(FormatTool.getDFE4().format(genomeScoreMap.get(genomes[i])[0])).append("\t");
					genomeTaxaSb.append(FormatTool.getDF4().format(genomeScoreMap.get(genomes[i])[1])).append("\t");

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
							taxScoreMaps[j].put(taxon, taxScoreMaps[j].get(taxon) * genomeScoreMap.get(genomes[i])[0]);
						} else {
							taxScoreMaps[j].put(taxon, genomeScoreMap.get(genomes[i])[0]);
						}
					}
				}
			}
			genomeWriter.close();

		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			logger.error("Error in writing genome information to " + genomeFile, e);
			System.err
					.println(format.format(new Date()) + "\t" + "Error in writing genome information to " + genomeFile);
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
			for (String exp : expNames) {
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
			for (String exp : expNames) {
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
			logger.error("Error in writing taxa information to " + taxaFile, e);
			System.err.println(format.format(new Date()) + "\t" + "Error in writing taxa information to " + taxaFile);
		}

		logger.info(getTaskName() + ": taxonomy analysis finished");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": taxonomy analysis finished");

		updateProgress(95, 100);

		updateMessage(getTaskName() + ": exporting report");

		logger.info(getTaskName() + ": exporting report started");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": exporting report started");

		File reportDir = new File(resultFolderFile, "report");
		if (!reportDir.exists()) {
			reportDir.mkdirs();
		}

		// Copy report templates using MetaPythonReportCopyTask (same as DDA workflow)
		MetaPythonReportCopyTask reportCopyTask = new MetaPythonReportCopyTask(resultFolderFile, logger);
		Thread copyThread = new Thread(reportCopyTask);
		copyThread.start();
		try {
			copyThread.join();
			if (!reportCopyTask.get()) {
				logger.error(getTaskName() + ": error in copying report template");
				System.err.println(
						format.format(new Date()) + "\t" + getTaskName() + ": error in copying report template");
				return false;
			}
		} catch (Exception e) {
			logger.error(getTaskName() + ": error in copying report template", e);
			System.err.println(format.format(new Date()) + "\t" + getTaskName() + ": error in copying report template: "
					+ e.getMessage());
			return false;
		}

		metaPar.setReportDir(reportDir);

		this.reportHtmlDir = new File(reportDir + "\\reports");

		if (!reportHtmlDir.exists()) {
			reportHtmlDir.mkdirs();
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
		 * this // doesn't work) File bubbleChartOutput = new
		 * File(metaPar.getReportDir() + "\\reports\\bubble_chart_report.html"); if
		 * (!bubbleChartOutput.exists() || bubbleChartOutput.length() == 0) { try {
		 * generateBubbleChartReport(jsOutput, bubbleChartOutput); } catch (Exception e)
		 * { logger.warn(getTaskName() +
		 * ": Bubble chart generation failed, continuing with other reports: " +
		 * e.getMessage()); System.out.println(format.format(new Date()) + "\t" +
		 * getTaskName() +
		 * ": Bubble chart generation failed, continuing with other reports"); } }
		 */
		MetaPythonReportTask task = new MetaPythonReportTask(metaPar.getThreadCount(), mrp, logger);
		task.setWorkflowType("DIA");
		task.initializeManifestGenerator(new File(metaPar.getReportDir(), "reports"));

		// Create initial manifest file to ensure index page can load
		ReportManifestGenerator manifestGen = new ReportManifestGenerator(new File(metaPar.getReportDir(), "reports"));
		manifestGen.createInitialManifest();

		this.summaryMetaFile = new File(this.resultFolderFile, "metainfo.tsv");
		try {
			metadata.exportMetadata(summaryMetaFile);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		this.exportSummary(task);

		if (!final_pep_txt.exists() || final_pep_txt.length() == 0) {
			return false;
		} else {
			File report_peptides_summary = new File(reportHtmlDir, "report_peptides_summary.html");
			if (!report_peptides_summary.exists() || report_peptides_summary.length() == 0) {
				if (summaryMetaFile.exists()) {
					task.addTask(MetaPythonReportTask.PEPTIDE, final_pep_txt, report_peptides_summary, summaryMetaFile);
				} else {
					task.addTask(MetaPythonReportTask.PEPTIDE, final_pep_txt, report_peptides_summary);
				}
			}
		}

		if (!final_pro_txt.exists() || final_pro_txt.length() == 0) {
			return false;
		} else {
			File report_proteinGroups_summary = new File(reportHtmlDir, "report_proteinGroups_summary.html");
			if (!report_proteinGroups_summary.exists() || report_proteinGroups_summary.length() == 0) {
				if (summaryMetaFile.exists()) {
					task.addTask(MetaPythonReportTask.PROTEIN, final_pro_txt, report_proteinGroups_summary,
							summaryMetaFile);
				} else {
					task.addTask(MetaPythonReportTask.PROTEIN, final_pro_txt, report_proteinGroups_summary);
				}
			}
		}

		if (!taxaReportFile.exists() || taxaReportFile.length() == 0) {
			return false;
		} else {
			File report_taxonomy_summary = new File(this.reportHtmlDir, "report_taxonomy_summary.html");
			if (!report_taxonomy_summary.exists() || report_taxonomy_summary.length() == 0) {

				if (summaryMetaFile != null && summaryMetaFile.exists()) {
					task.addTask(MetaPythonReportTask.TAXONOMY, taxaReportFile, report_taxonomy_summary,
							summaryMetaFile);
				} else {
					task.addTask(MetaPythonReportTask.TAXONOMY, taxaReportFile, report_taxonomy_summary);
				}
			}
		}

		if (!funcFile.exists() || funcFile.listFiles().length == 0) {
			return false;
		} else {
			File report_function_summary = new File(this.reportHtmlDir, "report_function_summary.html");
			if (!report_function_summary.exists() || report_function_summary.length() == 0) {

				if (summaryMetaFile != null && summaryMetaFile.exists()) {
					task.addTask(MetaPythonReportTask.FUNCTION, funcFile, report_function_summary, summaryMetaFile);
				} else {
					task.addTask(MetaPythonReportTask.FUNCTION, funcFile, report_function_summary);
				}
			}
		}

		task.run();

		try {
			boolean finish = task.get();
			if (finish) {

				logger.info(getTaskName() + ": report generation successes");
				System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": report generation successes");

				updateProgress(100, 100);
				return true;
			}

		} catch (InterruptedException | ExecutionException e) {
			// TODO Auto-generated catch block
			logger.error(getTaskName() + ": error in the report generation", e);
			System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": error in the report generation");
		}

		return true;
	}

	/**
	 * The version is different from the parent class
	 */
	protected boolean exportReportDiaNN() {

		updateMessage(getTaskName() + ": exporting report");

		logger.info(getTaskName() + ": exporting report started");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": exporting report started");

		File reportDir = new File(resultFolderFile, "report");
		if (!reportDir.exists()) {
			reportDir.mkdirs();
		}

		// Copy report templates using MetaPythonReportCopyTask (same as DDA workflow)
		MetaPythonReportCopyTask reportCopyTask = new MetaPythonReportCopyTask(resultFolderFile, logger);
		Thread copyThread = new Thread(reportCopyTask);
		copyThread.start();
		try {
			copyThread.join();
			if (!reportCopyTask.get()) {
				logger.error(getTaskName() + ": error in copying report template");
				System.err.println(
						format.format(new Date()) + "\t" + getTaskName() + ": error in copying report template");
				return false;
			}
		} catch (Exception e) {
			logger.error(getTaskName() + ": error in copying report template", e);
			System.err.println(format.format(new Date()) + "\t" + getTaskName() + ": error in copying report template: "
					+ e.getMessage());
			return false;
		}

		metaPar.setReportDir(reportDir);

		this.reportHtmlDir = new File(reportDir + "\\reports");

		if (!reportHtmlDir.exists()) {
			reportHtmlDir.mkdirs();
		}

		MetaPythonReportTask task = new MetaPythonReportTask(metaPar.getThreadCount(), mrp, logger);
		task.setWorkflowType("DIA");
		task.initializeManifestGenerator(new File(metaPar.getReportDir(), "reports"));

		// Create initial manifest file to ensure index page can load
		ReportManifestGenerator manifestGen = new ReportManifestGenerator(new File(metaPar.getReportDir(), "reports"));
		manifestGen.createInitialManifest();

		this.summaryMetaFile = new File(resultFolderFile, "metainfo.tsv");
		try {
			metadata.exportMetadata(summaryMetaFile);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		boolean summayrFinished = this.exportSummary(task);

		updateProgress(82, 100);

		this.final_pep_txt = new File(resultFolderFile, "final_peptides.tsv");

		boolean pepReportOk = false;
		if (!final_pep_txt.exists() || final_pep_txt.length() == 0) {
			return false;
		} else {
			File report_peptides_summary = new File(reportHtmlDir, "report_peptides_summary.html");
			if (!report_peptides_summary.exists() || report_peptides_summary.length() == 0) {
				if (summaryMetaFile.exists()) {
					pepReportOk = task.addTask(MetaPythonReportTask.PEPTIDE, final_pep_txt, report_peptides_summary,
							summaryMetaFile);
				} else {
					pepReportOk = task.addTask(MetaPythonReportTask.PEPTIDE, final_pep_txt, report_peptides_summary);
				}
			} else {
				pepReportOk = true;
			}
		}

		updateProgress(84, 100);

		this.final_pro_txt = new File(resultFolderFile, "final_proteins.tsv");

		boolean proReportOk = false;
		if (!final_pro_txt.exists() || final_pro_txt.length() == 0) {
			return false;
		} else {
			File report_proteinGroups_summary = new File(reportHtmlDir, "report_proteinGroups_summary.html");
			if (!report_proteinGroups_summary.exists() || report_proteinGroups_summary.length() == 0) {
				if (summaryMetaFile.exists()) {
					proReportOk = task.addTask(MetaPythonReportTask.PROTEIN, final_pro_txt,
							report_proteinGroups_summary, summaryMetaFile);
				} else {
					proReportOk = task.addTask(MetaPythonReportTask.PROTEIN, final_pro_txt,
							report_proteinGroups_summary);
				}
			} else {
				proReportOk = true;
			}
		}

		updateProgress(86, 100);

		boolean taxFuncReportOk = exportTaxaFunc(task);

		if (!summayrFinished || !pepReportOk || !proReportOk || !taxFuncReportOk) {
			updateProgress(99, 100);

			logger.info(getTaskName() + ": error in generating the reports");
			System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": error in generating the reports");

			return false;
		}

		task.run();

		try {
			boolean finish = task.get();
			if (finish) {

				logger.info(getTaskName() + ": report generation successes");
				System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": report generation successes");

				updateProgress(100, 100);
				return true;
			}

		} catch (InterruptedException | ExecutionException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		return false;
	}

	protected boolean exportSummary(MetaPythonReportTask task) {

		this.final_summary_txt = new File(resultFolderFile, "final_summary.tsv");
		try {

			PrintWriter countWriter = new PrintWriter(this.final_summary_txt);
			StringBuilder titlesb = new StringBuilder();
			titlesb.append("Raw file").append("\t");
			titlesb.append("Experiment").append("\t");
			titlesb.append("MS/MS").append("\t");
			titlesb.append("MS/MS Identified").append("\t");
			titlesb.append("MS/MS Identified [%]").append("\t");
			titlesb.append("Peptide Sequences Identified");

			countWriter.println(titlesb);

			int totalms2count = 0;

			HashMap<String, Integer> totalms2Map = new HashMap<String, Integer>();
			HashMap<String, HashSet<String>> totalUnipepMap = new HashMap<String, HashSet<String>>();
			HashSet<String> totalPepSet = new HashSet<String>();
			for (int i = 0; i < this.diaPp.length; i++) {
				String run = this.diaPp[i].getRunString();
				String sequence = this.diaPp[i].getSeqString();
				if (totalms2Map.containsKey(run)) {
					totalms2Map.put(run, totalms2Map.get(run) + 1);
					totalUnipepMap.get(run).add(sequence);
				} else {
					totalms2Map.put(run, 1);
					HashSet<String> set = new HashSet<String>();
					set.add(sequence);
					totalUnipepMap.put(run, set);
				}
				totalPepSet.add(sequence);
			}

			for (String file : totalUnipepMap.keySet()) {
				StringBuilder sb = new StringBuilder();
				sb.append(file).append("\t");

				if (this.expNameMap.containsKey(file)) {
					sb.append(expNameMap.get(file)).append("\t");
				} else {
					sb.append(file).append("\t");
				}
				int idCounts = totalms2Map.get(file);
				int pepCount = totalUnipepMap.get(file).size();
				sb.append(idCounts).append("\t");
				sb.append(idCounts).append("\t");
				sb.append("1.0").append("\t");
				sb.append(pepCount);

				countWriter.println(sb);

				totalms2count += idCounts;
			}

			StringBuilder sb = new StringBuilder();
			sb.append("Total").append("\t");
			sb.append("\t");
			sb.append(totalms2count).append("\t");
			sb.append(totalms2count).append("\t");
			sb.append("1.0").append("\t");
			sb.append(totalPepSet.size());

			countWriter.print(sb);

			countWriter.close();

		} catch (IOException e) {
			// TODO Auto-generated catch block
			logger.error(getTaskName() + ": error in writing PSM information to " + final_summary_txt.getName(), e);
			System.err.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": error in writing PSM information to " + final_summary_txt.getName());
		}

		if (!final_summary_txt.exists() || final_summary_txt.length() == 0) {
			return false;

		} else {
			boolean finish = false;
			File report_ID_summary = new File(reportHtmlDir, "report_ID_summary.html");
			if (summaryMetaFile != null && summaryMetaFile.exists()) {
				finish = task.addTask(MetaPythonReportTask.ID_SUMMARY, final_summary_txt, report_ID_summary,
						summaryMetaFile);
			} else {
				finish = task.addTask(MetaPythonReportTask.ID_SUMMARY, final_summary_txt, report_ID_summary);
			}

			return finish;
		}
	}

	protected boolean exportTaxaFunc(MetaPythonReportTask task) {

		logger.info(getTaskName() + ": taxonomy analysis and functional annotation started");
		System.out.println(format.format(new Date()) + "\t" + getTaskName()
				+ ": taxonomy analysis and functional annotation started");

		// functional annotation

		HashMap<String, Integer> proTaxIdMap = new HashMap<String, Integer>();

		logger.info(getTaskName() + ": functional annotation started");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": functional annotation started");

		BufferedReader quanPepReader = null;
		try {
			quanPepReader = new BufferedReader(new FileReader(this.final_pep_txt));
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			logger.error(getTaskName() + ": error in reading quantified peptides from " + this.quan_pep_file.getName(),
					e);
			System.err.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": error in reading quantified peptides from " + this.quan_pep_file.getName());
			return false;
		}

		HashMap<String, ArrayList<Double>> genomeQValueMap = new HashMap<String, ArrayList<Double>>();
		try {

			String line = quanPepReader.readLine();
			String[] title = line.split("\t");
			int proId = -1;
			int qvalueId = -1;

			for (int i = 0; i < title.length; i++) {
				if (title[i].equals("Proteins")) {
					proId = i;
				} else if (title[i].equals("PEP")) {
					qvalueId = i;
				}
			}

			while ((line = quanPepReader.readLine()) != null) {
				String[] cs = line.split("\t");
				double qvalue = Double.parseDouble(cs[qvalueId]);
				String[] pros = cs[proId].split(";");
				for (int i = 0; i < pros.length; i++) {
					if (pros[i].startsWith(magDb.getIdentifier())) {
						String genome = pros[i].split("_")[0];
						ArrayList<Double> qList;
						if (genomeQValueMap.containsKey(genome)) {
							qList = genomeQValueMap.get(genome);
						} else {
							qList = new ArrayList<Double>();
							genomeQValueMap.put(genome, qList);
						}
						qList.add(qvalue);
					}
				}
			}
			quanPepReader.close();

		} catch (IOException e) {
			// TODO Auto-generated catch block
			logger.error(getTaskName() + ": error in reading quantified peptides from " + this.quan_pep_file.getName(),
					e);
			System.err.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": error in reading quantified peptides from " + this.quan_pep_file.getName());
		}

		File funcFile = new File(resultFolderFile, "functional_annotation");
		if (!funcFile.exists()) {
			funcFile.mkdir();
		}

		final_pro_xml_file = new File(funcFile, "functions.xml");

		HashMap<String, double[]> proIntensityMap = new HashMap<String, double[]>();
		HashMap<String, double[]> genomeIntensityMap = new HashMap<String, double[]>();

		if (final_pro_xml_file.exists() && final_pro_xml_file.length() > 0) {

			try (BufferedReader quanProReader = new BufferedReader(new FileReader(this.final_pro_txt))) {
				String line = quanProReader.readLine();
				String[] title = line.split("\t");
				int proId = -1;

				ArrayList<String> tlist = new ArrayList<String>();
				for (int i = 0; i < title.length; i++) {
					if (title[i].startsWith("Intensity ")) {
						tlist.add(title[i]);
					} else if (title[i].equals("Majority protein IDs")) {
						proId = i;
					}
				}

				String[] fileNames = tlist.toArray(new String[tlist.size()]);

				int[] fileIds = new int[fileNames.length];
				for (int i = 0; i < fileNames.length; i++) {
					for (int j = 0; j < title.length; j++) {
						if (fileNames[i].equals(title[j])) {
							fileIds[i] = j;
						}
					}
				}

				while ((line = quanProReader.readLine()) != null) {

					String[] cs = line.split("\t");
					double[] intensity = new double[fileIds.length];
					for (int i = 0; i < fileIds.length; i++) {
						intensity[i] = Double.parseDouble(cs[fileIds[i]]);
					}

					String[] pros = cs[proId].split(";");
					for (int i = 0; i < pros.length; i++) {
						proIntensityMap.put(pros[i], intensity);

						if (pros[i].startsWith(magDb.getIdentifier())) {
							String genome = pros[i].split("_")[0];
							double[] genomeIntensity;
							if (genomeIntensityMap.containsKey(genome)) {
								genomeIntensity = genomeIntensityMap.get(genome);
							} else {
								genomeIntensity = new double[intensity.length];
							}

							for (int j = 0; j < genomeIntensity.length; j++) {
								genomeIntensity[j] += intensity[j] / (double) pros.length;
							}
							genomeIntensityMap.put(genome, genomeIntensity);
						}
					}
				}
				quanProReader.close();

			} catch (IOException e) {
				// TODO Auto-generated catch block
				logger.error(
						getTaskName() + ": error in reading quantified proteins from " + this.final_pro_txt.getName(),
						e);
				System.err.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": error in reading quantified proteins from " + this.final_pro_txt.getName());
			}

			logger.info(getTaskName() + ": proteins with functional annotations have been exported to "
					+ final_pro_xml_file.getName());
			System.out.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": proteins with functional annotations have been exported to " + final_pro_xml_file.getName());
		} else {

			ArrayList<MetaProtein> list = new ArrayList<MetaProtein>();
			String[] fileNames = null;
			try (BufferedReader quanProReader = new BufferedReader(new FileReader(this.final_pro_txt))) {
				String line = quanProReader.readLine();
				String[] title = line.split("\t");
				int id = -1;
				int proId = -1;
				int pepCountId = -1;
//				int rzPepCountId = -1;
				int scoreId = -1;

				ArrayList<String> tlist = new ArrayList<String>();
				for (int i = 0; i < title.length; i++) {
					if (title[i].startsWith("Intensity ")) {
						tlist.add(title[i]);
					} else if (title[i].equals("id")) {
						id = i;
					} else if (title[i].equals("Majority protein IDs")) {
						proId = i;
					} else if (title[i].equals("Peptide counts (all)")) {
						pepCountId = i;
					} else if (title[i].equals("Peptide counts (razor)")) {
//						rzPepCountId = i;
					} else if (title[i].equals("Score")) {
						scoreId = i;
					}
				}

				fileNames = tlist.toArray(new String[tlist.size()]);

				int[] fileIds = new int[fileNames.length];
				for (int i = 0; i < fileNames.length; i++) {
					for (int j = 0; j < title.length; j++) {
						if (fileNames[i].equals(title[j])) {
							fileIds[i] = j;
						}
					}
				}

				while ((line = quanProReader.readLine()) != null) {

					String[] cs = line.split("\t");
					String[] pros = cs[proId].split(";");

					if (pros[0].startsWith(magDb.getIdentifier())) {

						int groupId = Integer.parseInt(cs[id]);

						double[] intensity = new double[fileIds.length];
						for (int i = 0; i < fileIds.length; i++) {
							intensity[i] = Double.parseDouble(cs[fileIds[i]]);
						}

						String[] pepCountSs = cs[pepCountId].split(";");
//						String[] rzPepCountSs = cs[rzPepCountId].split(";");

						for (int i = 0; i < pros.length; i++) {

							int pepCount = i < pepCountSs.length ? Integer.parseInt(pepCountSs[i])
									: Integer.parseInt(pepCountSs[0]);
//							int rzPepCount = i < rzPepCountSs.length ? Integer.parseInt(rzPepCountSs[i])
//									: Integer.parseInt(rzPepCountSs[0]);

							MetaProtein mp = new MetaProtein(groupId, i + 1, pros[i], pepCount, pepCount,
									Double.parseDouble(cs[scoreId]), intensity);
							list.add(mp);

							if (pros[i].startsWith(magDb.getIdentifier())) {
								String genome = pros[i].split("_")[0];
								double[] genomeIntensity;
								if (genomeIntensityMap.containsKey(genome)) {
									genomeIntensity = genomeIntensityMap.get(genome);
								} else {
									genomeIntensity = new double[intensity.length];
								}

								for (int j = 0; j < genomeIntensity.length; j++) {
									genomeIntensity[j] += intensity[j] / (double) pros.length;
								}
								genomeIntensityMap.put(genome, genomeIntensity);
							}
						}
					}
				}
				quanProReader.close();

			} catch (IOException e) {
				// TODO Auto-generated catch block
				logger.error(
						getTaskName() + ": error in reading quantified proteins from " + this.final_pro_txt.getName(),
						e);
				System.err.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": error in reading quantified proteins from " + this.final_pro_txt.getName());
			}

			updateProgress(88, 100);

			MetaProtein[] metapros = list.toArray(new MetaProtein[list.size()]);
			MetaProteinAnnoMag[] mpas = funcAnnoSqlite(metapros, fileNames);
			for (int i = 0; i < mpas.length; i++) {
				String name = mpas[i].getPro().getName();
				if (NumberUtils.isCreatable(mpas[i].getTax_id())) {
					proTaxIdMap.put(name, Integer.parseInt(mpas[i].getTax_id()));
				} else {
					proTaxIdMap.put(name, 131567);
				}
			}
		}

		File funcResultFile = final_pro_xml_file.getParentFile();
		File funTsv = new File(funcResultFile, "functions.tsv");

		MetaProteinXMLReader2 funReader = null;
		try {
			funReader = new MetaProteinXMLReader2(final_pro_xml_file);
		} catch (DocumentException e) {
			// TODO Auto-generated catch block
			logger.error(getTaskName() + ": error in reading protein function information from "
					+ this.final_pro_xml_file.getName(), e);
			System.err.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": error in reading protein function information from " + this.final_pro_xml_file.getName());
			return false;
		}

		if (proTaxIdMap.size() == 0) {
			MetaProteinAnnoEggNog[] mpas = funReader.getProteins();
			for (int i = 0; i < mpas.length; i++) {
				String name = mpas[i].getPro().getName();
				if (NumberUtils.isCreatable(mpas[i].getTax_id())) {
					proTaxIdMap.put(name, Integer.parseInt(mpas[i].getTax_id()));
				} else {
					proTaxIdMap.put(name, 131567);
				}
			}
		}

		if (!funTsv.exists()) {
			System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": writing results to "
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
				logger.error(getTaskName() + ": error in writing protein function information to " + funTsv.getName(),
						e);
				System.err.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": error in writing protein function information to " + funTsv.getName());
			}
			System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": writing results to "
					+ funTsv.getName() + " finished");
		}

		File html = new File(funcResultFile, "functions.html");
		if (!html.exists()) {
			System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": writing results to "
					+ html.getName() + " started");
			try {
				funReader.exportHtml(html);
				// Copy the generated functions.html to the reports directory
				File reportFunctionsHtml = new File(metaPar.getReportDir(), "reports/functions.html");
				if (reportFunctionsHtml.getParentFile() != null && !reportFunctionsHtml.getParentFile().exists()) {
					reportFunctionsHtml.getParentFile().mkdirs();
				}
				Files.copy(html.toPath(), reportFunctionsHtml.toPath(), StandardCopyOption.REPLACE_EXISTING);
				logger.info("Copied functions.html to reports directory: " + reportFunctionsHtml.getAbsolutePath());
			} catch (IOException e) {
				// TODO Auto-generated catch block
				logger.error(getTaskName() + ": error in writing protein function information to " + html.getName(), e);
				System.err.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": error in writing protein function information to " + html.getName());
			}
			System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": writing results to "
					+ html.getName() + " finished");
		}

		File report_function_summary = new File(this.reportHtmlDir, "report_function_summary.html");
		if (!report_function_summary.exists() || report_function_summary.length() == 0) {

			if (summaryMetaFile != null && summaryMetaFile.exists()) {
				task.addTask(MetaPythonReportTask.FUNCTION, funcResultFile, report_function_summary, summaryMetaFile);
			} else {
				task.addTask(MetaPythonReportTask.FUNCTION, funcResultFile, report_function_summary);
			}
		}

		logger.info(getTaskName() + ": functional annotation finished");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": functional annotation finished");

		updateProgress(92, 100);

		// taxonomy analysis

		logger.info(getTaskName() + ": taxonomy analysis started");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": taxonomy analysis started");

		File taxFile = new File(resultFolderFile, "taxonomy_analysis");
		if (!taxFile.exists()) {
			taxFile.mkdir();
		}

		HashMap<String, Integer> psmCountMap = new HashMap<String, Integer>();
		HashMap<String, HashSet<String>> pepCountMap = new HashMap<String, HashSet<String>>();
		HashMap<String, HashSet<String>> uniquePepCountMap = new HashMap<String, HashSet<String>>();

		HashMap<String, HashSet<String>> genomeProMap = new HashMap<String, HashSet<String>>();
		for (int i = 0; i < diaPp.length; i++) {

			String baseSequence = diaPp[i].getSeqString();
			String[] pros = diaPp[i].getProGroups();

			for (int j = 0; j < pros.length; j++) {

				if (psmCountMap.containsKey(pros[j])) {
					psmCountMap.put(pros[j], psmCountMap.get(pros[j]) + 1);

					pepCountMap.get(pros[j]).add(baseSequence);
					if (pros.length == 1) {
						uniquePepCountMap.get(pros[j]).add(baseSequence);
					}
				} else {
					psmCountMap.put(pros[j], 1);

					HashSet<String> pepSet = new HashSet<String>();
					pepSet.add(baseSequence);
					pepCountMap.put(pros[j], pepSet);

					HashSet<String> unipepSet = new HashSet<String>();
					if (pros.length == 1) {
						unipepSet.add(baseSequence);
					}
					uniquePepCountMap.put(pros[j], unipepSet);
				}

				if (pros[j].startsWith(magDb.getIdentifier())) {
					String genome = pros[j].split("_")[0];

					if (psmCountMap.containsKey(genome)) {
						psmCountMap.put(genome, psmCountMap.get(genome) + 1);

						pepCountMap.get(genome).add(baseSequence);
						if (pros.length == 1) {
							uniquePepCountMap.get(genome).add(baseSequence);
						}
					} else {
						psmCountMap.put(genome, 1);

						HashSet<String> pepSet = new HashSet<String>();
						pepSet.add(baseSequence);
						pepCountMap.put(genome, pepSet);

						HashSet<String> unipepSet = new HashSet<String>();
						if (pros.length == 1) {
							unipepSet.add(baseSequence);
						}
						uniquePepCountMap.put(genome, unipepSet);
					}

					if (genomeProMap.containsKey(genome)) {
						genomeProMap.get(genome).add(pros[j]);
					} else {
						HashSet<String> proSet = new HashSet<String>();
						proSet.add(pros[j]);
						genomeProMap.put(genome, proSet);
					}
				}
			}
		}

		return exportGenome(genomeIntensityMap, psmCountMap, pepCountMap, uniquePepCountMap, genomeQValueMap, taxFile,
				expNames, task);
	}

	protected void cleanTempFiles() {
		if (this.separateResultFolders != null) {
			for (int i = 0; i < this.separateResultFolders.length; i++) {
				File secondFile = new File(separateResultFolders[i], "secondSearch");
				File secondLibFile = new File(secondFile, "secondSearch.speclib");
				if (secondLibFile.exists()) {
					FileUtils.deleteQuietly(secondLibFile);
				}

				File thirdFile = new File(separateResultFolders[i], "thirdSearch");
				File thirdLibFile = new File(thirdFile, "thirdSearch.speclib");
				if (thirdLibFile.exists()) {
					FileUtils.deleteQuietly(thirdLibFile);
				}

				File lashLibFile = new File(separateResultFolders[i], "lastSearch.speclib");
				if (lashLibFile.exists()) {
					FileUtils.deleteQuietly(lashLibFile);
				}
			}
		}

		File combineLibFile = new File(resultFolderFile, "combined.speclib");
		if (combineLibFile.exists()) {
			FileUtils.deleteQuietly(combineLibFile);
		}

		File shuffledLibFile = new File(resultFolderFile, "shuffled.speclib");
		if (shuffledLibFile.exists()) {
			FileUtils.deleteQuietly(shuffledLibFile);
		}
	}

	@SuppressWarnings("unused")
	private void dbFromPredictedPeptide(File[] predictedFiles, HashSet<String> pepSet, File outputDbFile,
			int upperLimit) {

		HashMap<String, HashMap<String, Double>> pepGenomeMap = new HashMap<String, HashMap<String, Double>>();

		for (int i = 0; i < predictedFiles.length; i++) {
			try (BufferedReader reader = new BufferedReader(new FileReader(predictedFiles[i]))) {
				String pline = reader.readLine();
				while ((pline = reader.readLine()) != null) {
					String[] cs = pline.split("\t");
					String genome = cs[0].substring(0, cs[0].indexOf("_"));
					double score = Double.parseDouble(cs[2]);
					if (pepGenomeMap.containsKey(cs[1])) {
						pepGenomeMap.get(cs[1]).put(genome, score);
					} else {
						HashMap<String, Double> gMap = new HashMap<String, Double>();
						gMap.put(genome, score);
						pepGenomeMap.put(cs[1], gMap);
					}
				}
				reader.close();
			} catch (IOException e) {

				logger.error(getTaskName() + ": error in reading " + predictedFiles[i], e);
				System.err.println(
						format.format(new Date()) + "\t" + getTaskName() + ": error in reading " + predictedFiles[i]);
				return;
			}
		}

		HashMap<String, Double> genomeCountMap = new HashMap<String, Double>();
		for (String pep : pepGenomeMap.keySet()) {
			HashMap<String, Double> map = pepGenomeMap.get(pep);
			double score = 1.0 / (double) map.size();
			for (String genome : map.keySet()) {
				if (genomeCountMap.containsKey(genome)) {
					genomeCountMap.put(genome, genomeCountMap.get(genome) + score);
				} else {
					genomeCountMap.put(genome, score);
				}
			}
		}

		HashMap<String, Double> pepScoreMap = new HashMap<String, Double>();
		for (String pep : pepGenomeMap.keySet()) {
			if (pepGenomeMap.containsKey(pep)) {
				HashMap<String, Double> map = pepGenomeMap.get(pep);
				double score = 0;

				for (String genome : map.keySet()) {
					score += genomeCountMap.get(genome) * map.get(genome);
				}

				score = score / (double) map.size();

				pepScoreMap.put(pep, score);
			}
		}

		String[] peps = pepScoreMap.keySet().toArray(String[]::new);
		Arrays.sort(peps, new Comparator<String>() {

			@Override
			public int compare(String o1, String o2) {
				// TODO Auto-generated method stub

				if (pepScoreMap.get(o1) < pepScoreMap.get(o2)) {
					return 1;
				}
				if (pepScoreMap.get(o1) > pepScoreMap.get(o2)) {
					return -1;
				}

				return 0;
			}
		});

		HashSet<String> totalPepSet = new HashSet<String>(pepSet);
		int end = peps.length > upperLimit ? upperLimit : peps.length;
		for (int j = 0; j < end; j++) {
			totalPepSet.add(peps[j]);
		}

		try (PrintWriter writer = new PrintWriter(outputDbFile)) {
			int id = 0;
			for (String pep : totalPepSet) {
				writer.println(">pro" + id + " pro" + id);
				writer.println(pep);
				id++;
			}
			writer.close();
		} catch (IOException e) {

			logger.error(getTaskName() + ": error in writing to " + outputDbFile, e);
			System.err.println(
					format.format(new Date()) + "\t" + getTaskName() + ": error in writing to " + outputDbFile);
		}
	}

	private void dbFromPredictedPeptideFunc(File predictedFolder, File funcFolder, File outputDbFile,
			HashSet<String> genomeSet, HashSet<String> pepSet, HashSet<String> excludeSet, int upperLimit,
			boolean identifiedNOG) {

		HashMap<String, HashSet<String>> genomeProMap = new HashMap<String, HashSet<String>>();
		for (String genome : genomeSet) {
			if (identifiedNOG) {
				File predictedFile = new File(predictedFolder, genome + ".predicted.tsv");
				HashSet<String> proSet = new HashSet<String>();
				if (predictedFile.exists()) {
					try (BufferedReader reader = new BufferedReader(new FileReader(predictedFile))) {
						String pline = reader.readLine();
						while ((pline = reader.readLine()) != null) {
							String[] cs = pline.split("\t");
							if (pepSet.contains(cs[1])) {
								proSet.add(cs[0]);
							}
						}
						reader.close();
					} catch (IOException e) {

						logger.error(getTaskName() + ": error in reading " + predictedFile, e);
						System.err.println(format.format(new Date()) + "\t" + getTaskName() + ": error in reading "
								+ predictedFile);
						return;
					}
				}

				HashSet<String> eggnogSet = new HashSet<String>();
				File funcFile = new File(funcFolder, genome + "_eggNOG.tsv");
				if (funcFile.exists()) {
					try (BufferedReader reader = new BufferedReader(new FileReader(funcFile))) {
						String pline = reader.readLine();
						String[] title = pline.split("\t");
						int eggnogId = -1;
						for (int i = 0; i < title.length; i++) {
							if (title[i].equals("eggNOG_OGs")) {
								eggnogId = i;
								break;
							}
						}
						while ((pline = reader.readLine()) != null) {
							String[] cs = pline.split("\t");
							if (proSet.contains(cs[0])) {
								String[] eggnogs = cs[eggnogId].split(",");
								for (int j = 0; j < eggnogs.length; j++) {
									eggnogSet.add(eggnogs[j]);
								}
							}
						}
						reader.close();
					} catch (IOException e) {

						logger.error(getTaskName() + ": error in reading " + funcFile, e);
						System.err.println(
								format.format(new Date()) + "\t" + getTaskName() + ": error in reading " + funcFile);
						return;
					}

					HashSet<String> proSet2 = new HashSet<String>();
					try (BufferedReader reader = new BufferedReader(new FileReader(funcFile))) {
						String pline = reader.readLine();
						String[] title = pline.split("\t");
						int eggnogId = -1;
						for (int i = 0; i < title.length; i++) {
							if (title[i].equals("eggNOG_OGs")) {
								eggnogId = i;
								break;
							}
						}
						while ((pline = reader.readLine()) != null) {
							String[] cs = pline.split("\t");
							boolean use = false;
							String[] eggnogs = cs[eggnogId].split(",");
							for (int j = 0; j < eggnogs.length; j++) {
								if (eggnogSet.contains(eggnogs[j])) {
									use = true;
									break;
								}
							}

							if (use) {
								proSet2.add(cs[0]);
							}
						}
						reader.close();
					} catch (IOException e) {

						logger.error(getTaskName() + ": error in reading " + funcFile, e);
						System.err.println(
								format.format(new Date()) + "\t" + getTaskName() + ": error in reading " + funcFile);
						return;
					}

					genomeProMap.put(genome, proSet2);
				}
			} else {
				File funcFile = new File(funcFolder, genome + "_eggNOG.tsv");
				HashSet<String> proSet2 = new HashSet<String>();
				try (BufferedReader reader = new BufferedReader(new FileReader(funcFile))) {
					String pline = reader.readLine();
					while ((pline = reader.readLine()) != null) {
						String[] cs = pline.split("\t");
						proSet2.add(cs[0]);
					}
					reader.close();
				} catch (IOException e) {

					logger.error(getTaskName() + ": error in reading " + funcFile, e);
					System.err.println(
							format.format(new Date()) + "\t" + getTaskName() + ": error in reading " + funcFile);
					return;
				}

				genomeProMap.put(genome, proSet2);
			}
		}

		HashMap<String, HashMap<String, Double>> pepGenomeMap = new HashMap<String, HashMap<String, Double>>();
		for (String genome : genomeSet) {
			File predictedFile = new File(predictedFolder, genome + ".predicted.tsv");
			HashSet<String> proSet = genomeProMap.get(genome);
			try (BufferedReader reader = new BufferedReader(new FileReader(predictedFile))) {
				String pline = reader.readLine();
				while ((pline = reader.readLine()) != null) {
					String[] cs = pline.split("\t");
					if (proSet.contains(cs[0])) {
						double score = Double.parseDouble(cs[2]);
						if (pepGenomeMap.containsKey(cs[1])) {
							pepGenomeMap.get(cs[1]).put(genome, score);
						} else {
							HashMap<String, Double> gMap = new HashMap<String, Double>();
							gMap.put(genome, score);
							pepGenomeMap.put(cs[1], gMap);
						}
					}
				}
				reader.close();
			} catch (IOException e) {

				logger.error(getTaskName() + ": error in reading " + predictedFile, e);
				System.err.println(
						format.format(new Date()) + "\t" + getTaskName() + ": error in reading " + predictedFile);
				return;
			}
		}

		HashMap<String, Double> genomeCountMap = new HashMap<String, Double>();
		for (String pep : pepGenomeMap.keySet()) {
			HashMap<String, Double> map = pepGenomeMap.get(pep);
			double score = 1.0 / (double) map.size();
			for (String genome : map.keySet()) {
				if (genomeCountMap.containsKey(genome)) {
					genomeCountMap.put(genome, genomeCountMap.get(genome) + score);
				} else {
					genomeCountMap.put(genome, score);
				}
			}
		}

		HashMap<String, Double> pepScoreMap = new HashMap<String, Double>();
		for (String pep : pepGenomeMap.keySet()) {
			if (pepGenomeMap.containsKey(pep)) {
				HashMap<String, Double> map = pepGenomeMap.get(pep);
				double score = 0;

				for (String genome : map.keySet()) {
					score += genomeCountMap.get(genome) * map.get(genome);
				}

				score = score / (double) map.size();

				pepScoreMap.put(pep, score);
			}
		}

		String[] peps = pepScoreMap.keySet().toArray(String[]::new);
		Arrays.sort(peps, new Comparator<String>() {

			@Override
			public int compare(String o1, String o2) {
				// TODO Auto-generated method stub

				if (pepScoreMap.get(o1) < pepScoreMap.get(o2)) {
					return 1;
				}
				if (pepScoreMap.get(o1) > pepScoreMap.get(o2)) {
					return -1;
				}

				return 0;
			}
		});

		try (PrintWriter writer = new PrintWriter(outputDbFile)) {
			int id = 0;
			for (String pep : pepSet) {
				writer.println(">pro" + id + " pro" + id);
				writer.println(pep);
				id++;
			}
			for (String pep : peps) {
				if (!pepSet.contains(pep) && !excludeSet.contains(pep)) {
					writer.println(">pro" + id + " pro" + id);
					writer.println(pep);
					id++;
				}
				if (id >= upperLimit) {
					break;
				}
			}

			writer.close();
		} catch (IOException e) {

			logger.error(getTaskName() + ": error in writing to " + outputDbFile, e);
			System.err.println(
					format.format(new Date()) + "\t" + getTaskName() + ": error in writing to " + outputDbFile);
		}
	}

	@SuppressWarnings("unused")
	private void dbFromPredictedPeptide(File[] predictedFiles, HashSet<String> pepSet, HashSet<String> notUsedSet,
			File outputDbFile, int upperLimit) {

		HashMap<String, HashMap<String, Double>> pepGenomeMap = new HashMap<String, HashMap<String, Double>>();

		for (int i = 0; i < predictedFiles.length; i++) {
			try (BufferedReader reader = new BufferedReader(new FileReader(predictedFiles[i]))) {
				String pline = reader.readLine();
				while ((pline = reader.readLine()) != null) {
					String[] cs = pline.split("\t");
					String genome = cs[0].substring(0, cs[0].indexOf("_"));
					double score = Double.parseDouble(cs[2]);
					if (pepGenomeMap.containsKey(cs[1])) {
						pepGenomeMap.get(cs[1]).put(genome, score);
					} else {
						HashMap<String, Double> gMap = new HashMap<String, Double>();
						gMap.put(genome, score);
						pepGenomeMap.put(cs[1], gMap);
					}
				}
				reader.close();
			} catch (IOException e) {

				logger.error(getTaskName() + ": error in reading " + predictedFiles[i], e);
				System.err.println(
						format.format(new Date()) + "\t" + getTaskName() + ": error in reading " + predictedFiles[i]);
				return;
			}
		}

		HashMap<String, Double> genomeCountMap = new HashMap<String, Double>();
		for (String pep : pepGenomeMap.keySet()) {
			HashMap<String, Double> map = pepGenomeMap.get(pep);
			double score = 1.0 / (double) map.size();
			for (String genome : map.keySet()) {
				if (genomeCountMap.containsKey(genome)) {
					genomeCountMap.put(genome, genomeCountMap.get(genome) + score);
				} else {
					genomeCountMap.put(genome, score);
				}
			}
		}

		HashMap<String, Double> pepScoreMap = new HashMap<String, Double>();
		for (String pep : pepGenomeMap.keySet()) {
			if (pepGenomeMap.containsKey(pep)) {
				HashMap<String, Double> map = pepGenomeMap.get(pep);
				double score = 0;

				for (String genome : map.keySet()) {
					score += genomeCountMap.get(genome) * map.get(genome);
				}

				score = score / (double) map.size();

				pepScoreMap.put(pep, score);
			}
		}

		String[] peps = pepScoreMap.keySet().toArray(String[]::new);
		Arrays.sort(peps, new Comparator<String>() {

			@Override
			public int compare(String o1, String o2) {
				// TODO Auto-generated method stub

				if (pepScoreMap.get(o1) < pepScoreMap.get(o2)) {
					return 1;
				}
				if (pepScoreMap.get(o1) > pepScoreMap.get(o2)) {
					return -1;
				}

				return 0;
			}
		});

		HashSet<String> totalPepSet = new HashSet<String>(pepSet);
		for (int j = 0; j < peps.length; j++) {
			if (!notUsedSet.contains(peps[j])) {
				totalPepSet.add(peps[j]);
				if (totalPepSet.size() > upperLimit) {
					break;
				}
			}
		}

		try (PrintWriter writer = new PrintWriter(outputDbFile)) {
			int id = 0;
			for (String pep : totalPepSet) {
				writer.println(">pro" + id + " pro" + id);
				writer.println(pep);
				id++;
			}
			writer.close();
		} catch (IOException e) {

			logger.error(getTaskName() + ": error in writing to " + outputDbFile, e);
			System.err.println(
					format.format(new Date()) + "\t" + getTaskName() + ": error in writing to " + outputDbFile);
		}
	}

	protected void dbFromModelPeptide(String[] modelRankFiles, HashSet<String> genomeSet, HashSet<String> pepSet,
			File outputDbFile, int upperLimit) {

		logger.info(getTaskName() + ": generating database to " + outputDbFile + " started");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": generating database to " + outputDbFile
				+ " started");

		HashSet<String> totalPepSet = new HashSet<String>();
		totalPepSet.addAll(pepSet);

		for (int i = 0; i < modelRankFiles.length; i++) {

			File modelRankFile = new File(modelRankFiles[i]);
			File modelIntensityFile = new File(modelRankFiles[i].replace("_pep_rank.tsv", "_pytorch_predictions.tsv"));

			HashMap<String, Double> pepIntensityMap = new HashMap<String, Double>();
			HashMap<String, HashSet<String>> pepGenomeMap = new HashMap<String, HashSet<String>>();

			try (BufferedReader reader = new BufferedReader(new FileReader(modelRankFile))) {
				String pline = reader.readLine();
				String genome = "";
				boolean have = false;
				while ((pline = reader.readLine()) != null) {
					if (pline.startsWith("pro")) {
						continue;
					} else if (genome.length() > 0 && pline.startsWith(genome)) {
						if (have) {
							String[] cs = pline.split("\t");

							pepIntensityMap.put(cs[1], 0.0);

							if (pepGenomeMap.containsKey(cs[1])) {
								pepGenomeMap.get(cs[1]).add(genome);
							} else {
								HashSet<String> gset = new HashSet<String>();
								gset.add(genome);
								pepGenomeMap.put(cs[1], gset);
							}
						}
					} else {
						String[] cs = pline.split("\t");
						genome = cs[0].substring(0, cs[0].indexOf("_"));
						if (genomeSet.contains(genome)) {
							have = true;
						} else {
							have = false;
						}
						if (have) {
							pepIntensityMap.put(cs[1], 0.0);

							if (pepGenomeMap.containsKey(cs[1])) {
								pepGenomeMap.get(cs[1]).add(genome);
							} else {
								HashSet<String> gset = new HashSet<String>();
								gset.add(genome);
								pepGenomeMap.put(cs[1], gset);
							}
						}
					}
				}
				reader.close();
			} catch (IOException e) {

				logger.error(getTaskName() + ": error in reading model information from " + modelRankFile, e);
				System.err.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": error in reading model information from  " + modelRankFile);
				return;
			}

			try (BufferedReader reader = new BufferedReader(new FileReader(modelIntensityFile))) {
				String line = reader.readLine();
				while ((line = reader.readLine()) != null) {
					String[] cs = line.split(",");
					if (pepIntensityMap.containsKey(cs[1])) {
						double intensity = Double.parseDouble(cs[2]);
						if (intensity > pepIntensityMap.get(cs[1])) {
							pepIntensityMap.put(cs[1], intensity);
						}
					}
				}
				reader.close();
			} catch (IOException e) {
				logger.error(getTaskName() + ": error in reading model information from  " + modelIntensityFile, e);
				System.err.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": error in reading model information from  " + modelIntensityFile);
				return;
			}

			HashMap<String, Double> genomeCountMap = new HashMap<String, Double>();
			for (String pep : pepGenomeMap.keySet()) {
				if (pepGenomeMap.containsKey(pep)) {
					HashSet<String> set = pepGenomeMap.get(pep);
					double score = 1.0 / (double) set.size();
					for (String genome : set) {
						if (genomeCountMap.containsKey(genome)) {
							genomeCountMap.put(genome, genomeCountMap.get(genome) + score);
						} else {
							genomeCountMap.put(genome, score);
						}
					}
				}
			}

			HashMap<String, Double> pepScoreMap = new HashMap<String, Double>();
			for (String pep : pepIntensityMap.keySet()) {
				if (pepGenomeMap.containsKey(pep)) {
					HashSet<String> set = pepGenomeMap.get(pep);
					double score = 0;

					for (String genome : set) {
						score += genomeCountMap.get(genome);
					}

					score = score / (double) set.size() * pepIntensityMap.get(pep);

					pepScoreMap.put(pep, score);
				}
			}

			String[] peps = pepScoreMap.keySet().toArray(String[]::new);
			Arrays.sort(peps, new Comparator<String>() {

				@Override
				public int compare(String o1, String o2) {
					// TODO Auto-generated method stub

					if (pepScoreMap.get(o1) < pepScoreMap.get(o2)) {
						return 1;
					}
					if (pepScoreMap.get(o1) > pepScoreMap.get(o2)) {
						return -1;
					}

					return 0;
				}
			});

			int end = peps.length > upperLimit ? upperLimit : peps.length;
			for (int j = 0; j < end; j++) {
				totalPepSet.add(peps[j]);
			}
		}

		try (PrintWriter writer = new PrintWriter(outputDbFile)) {
			int id = 0;
			for (String pep : totalPepSet) {
				writer.println(">pro" + id + " pro" + id);
				writer.println(pep);
				id++;
			}
			writer.close();
		} catch (IOException e) {

			logger.error(getTaskName() + ": error in writing to " + outputDbFile, e);
			System.err.println(
					format.format(new Date()) + "\t" + getTaskName() + ": error in writing to " + outputDbFile);
		}

		logger.info(getTaskName() + ": generating database to " + outputDbFile + " finished");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": generating database to " + outputDbFile
				+ " finished");
	}

	protected void dbFromModelPeptide(String[] modelRankFiles, HashSet<String> genomeSet, HashSet<String> pepSet,
			HashSet<String> notUsedSet, File outputDbFile, int upperLimit) {

		logger.info(getTaskName() + ": generating database to " + outputDbFile + " started");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": generating database to " + outputDbFile
				+ " started");

		HashSet<String> totalPepSet = new HashSet<String>();
		totalPepSet.addAll(pepSet);

		for (int i = 0; i < modelRankFiles.length; i++) {
			File modelRankFile = new File(modelRankFiles[i]);
			File modelIntensityFile = new File(modelRankFiles[i].replace("_pep_rank.tsv", "_pytorch_predictions.tsv"));

			HashMap<String, Double> pepIntensityMap = new HashMap<String, Double>();
			HashMap<String, HashSet<String>> pepGenomeMap = new HashMap<String, HashSet<String>>();

			try (BufferedReader reader = new BufferedReader(new FileReader(modelRankFile))) {
				String pline = reader.readLine();
				String genome = "";
				boolean have = false;
				while ((pline = reader.readLine()) != null) {
					if (pline.startsWith("pro")) {
						continue;
					} else if (genome.length() > 0 && pline.startsWith(genome)) {
						if (have) {
							String[] cs = pline.split("\t");

							pepIntensityMap.put(cs[1], 0.0);

							if (pepGenomeMap.containsKey(cs[1])) {
								pepGenomeMap.get(cs[1]).add(genome);
							} else {
								HashSet<String> gset = new HashSet<String>();
								gset.add(genome);
								pepGenomeMap.put(cs[1], gset);
							}
						}
					} else {
						String[] cs = pline.split("\t");
						genome = cs[0].substring(0, cs[0].indexOf("_"));
						if (genomeSet.contains(genome)) {
							have = true;
						} else {
							have = false;
						}
						if (have) {
							pepIntensityMap.put(cs[1], 0.0);

							if (pepGenomeMap.containsKey(cs[1])) {
								pepGenomeMap.get(cs[1]).add(genome);
							} else {
								HashSet<String> gset = new HashSet<String>();
								gset.add(genome);
								pepGenomeMap.put(cs[1], gset);
							}
						}
					}
				}
				reader.close();
			} catch (IOException e) {

				logger.error(getTaskName() + ": error in reading model information from  " + modelRankFile, e);
				System.err.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": error in reading model information from  " + modelRankFile);
				return;
			}

			try (BufferedReader reader = new BufferedReader(new FileReader(modelIntensityFile))) {
				String line = reader.readLine();
				while ((line = reader.readLine()) != null) {
					String[] cs = line.split("[\t,]");
					if (pepIntensityMap.containsKey(cs[1])) {
						double intensity = Double.parseDouble(cs[2]);
						if (intensity > pepIntensityMap.get(cs[1])) {
							pepIntensityMap.put(cs[1], intensity);
						}
					}
				}
				reader.close();
			} catch (IOException e) {
				logger.error(getTaskName() + ": error in reading model information from  " + modelIntensityFile, e);
				System.err.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": error in reading model information from  " + modelIntensityFile);
				return;
			}

			HashMap<String, Double> genomeCountMap = new HashMap<String, Double>();
			for (String pep : pepGenomeMap.keySet()) {
				if (pepGenomeMap.containsKey(pep)) {
					HashSet<String> set = pepGenomeMap.get(pep);
					double score = 1.0 / (double) set.size();
					for (String genome : set) {
						if (genomeCountMap.containsKey(genome)) {
							genomeCountMap.put(genome, genomeCountMap.get(genome) + score);
						} else {
							genomeCountMap.put(genome, score);
						}
					}
				}
			}

			HashMap<String, Double> pepScoreMap = new HashMap<String, Double>();
			for (String pep : pepIntensityMap.keySet()) {
				if (pepGenomeMap.containsKey(pep)) {
					HashSet<String> set = pepGenomeMap.get(pep);
					double score = 0;

					for (String genome : set) {
						score += genomeCountMap.get(genome);
					}

					score = score / (double) set.size() * pepIntensityMap.get(pep);

					pepScoreMap.put(pep, score);
				}
			}

			String[] peps = pepScoreMap.keySet().toArray(String[]::new);
			Arrays.sort(peps, new Comparator<String>() {

				@Override
				public int compare(String o1, String o2) {
					// TODO Auto-generated method stub

					if (pepScoreMap.get(o1) < pepScoreMap.get(o2)) {
						return 1;
					}
					if (pepScoreMap.get(o1) > pepScoreMap.get(o2)) {
						return -1;
					}

					return 0;
				}
			});

			for (int j = 0; j < peps.length; j++) {
				if (!notUsedSet.contains(peps[j])) {
					totalPepSet.add(peps[j]);
					if (totalPepSet.size() >= upperLimit) {
						break;
					}
				}
			}
		}

		try (PrintWriter writer = new PrintWriter(outputDbFile)) {
			int id = 0;
			for (String pep : totalPepSet) {
				writer.println(">pro" + id + " pro" + id);
				writer.println(pep);
				id++;
			}
			writer.close();
		} catch (IOException e) {

			logger.error(getTaskName() + ": error in writing to " + outputDbFile, e);
			System.err.println(
					format.format(new Date()) + "\t" + getTaskName() + ": error in writing to " + outputDbFile);
		}

		logger.info(getTaskName() + ": generating database to " + outputDbFile + " finished");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": generating database to " + outputDbFile
				+ " finished");
	}

	@SuppressWarnings({ "unused", "unchecked" })
	private void dbFromModelPeptideNotUsed(String[] modelPath, HashSet<String> genomeSet, HashSet<String> pepSet,
			HashSet<String> notUsedSet, File outputDbFile) {

		logger.info(getTaskName() + ": generating database to " + outputDbFile + " started");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": generating database to " + outputDbFile
				+ " started");

		HashSet<String>[] addSets = new HashSet[10];
		for (int i = 0; i < addSets.length; i++) {
			addSets[i] = new HashSet<String>();
		}

		for (int i = 0; i < modelPath.length; i++) {
			try (BufferedReader reader = new BufferedReader(new FileReader(modelPath[i]))) {
				String pline = reader.readLine();
				String genome = "";
				boolean have = false;
				while ((pline = reader.readLine()) != null) {
					if (genome.length() > 0 && pline.startsWith(genome)) {
						if (have) {
							String[] cs = pline.split("\t");
							if (!notUsedSet.contains(cs[1])) {
								double score = Double.parseDouble(cs[2]);

								int id = (int) (score * 10.0);

								if (id < addSets.length) {
									addSets[id].add(cs[1]);
								} else {
									addSets[addSets.length - 1].add(cs[1]);
								}
							}
						}
					} else {
						String[] cs = pline.split("\t");
						genome = cs[0].substring(0, cs[0].indexOf("_"));
						if (genomeSet.contains(genome)) {
							have = true;
						} else {
							have = false;
						}
						if (have) {
							if (!notUsedSet.contains(cs[1])) {
								double score = Double.parseDouble(cs[2]);
								int id = (int) (score * 10.0);

								if (id < addSets.length) {
									addSets[id].add(cs[1]);
								} else {
									addSets[addSets.length - 1].add(cs[1]);
								}
							}
						}
					}
				}
				reader.close();
			} catch (IOException e) {

				logger.error(getTaskName() + ": error in reading " + modelPath[i], e);
				System.err.println(
						format.format(new Date()) + "\t" + getTaskName() + ": error in reading " + modelPath[i]);
			}
		}

		HashSet<String> addSet = new HashSet<String>();
		for (int i = 1; i <= addSets.length; i++) {
			addSet.addAll(addSets[addSets.length - i]);
			if (addSet.size() > 1000000) {
				break;
			}
		}

		notUsedSet.addAll(addSet);
		try (PrintWriter writer = new PrintWriter(outputDbFile)) {
			int id = 0;
			for (String pep : pepSet) {
				writer.println(">pro" + id + " pro" + id);
				writer.println(pep);
				id++;
			}
			for (String pep : addSet) {
				writer.println(">pro" + id + " pro" + id);
				writer.println(pep);
				id++;
			}
			writer.close();
		} catch (IOException e) {

			logger.error(getTaskName() + ": error in writing to " + outputDbFile, e);
			System.err.println(
					format.format(new Date()) + "\t" + getTaskName() + ": error in writing to " + outputDbFile);
		}

		logger.info(getTaskName() + ": generating database to " + outputDbFile + " finished");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": generating database to " + outputDbFile
				+ " finished");
	}

	protected HashSet<String> getGenomeSet(File resultFile, double filterThreshold, HashSet<String> pepSet,
			int maxIterativeCount, boolean breakAdvance, int threadCount) {
		ConcurrentHashMap<String, HashMap<String, HashSet<String>>> matchedMap = getGenomeFromPep(pepSet, threadCount);
		if (matchedMap == null) {
			return null;
		} else {
			HashMap<String, HashSet<String>> genomePepMap = new HashMap<String, HashSet<String>>();
			for (String genome : matchedMap.keySet()) {
				HashSet<String> genomePepSet = new HashSet<String>();
				HashMap<String, HashSet<String>> proPepMap = matchedMap.get(genome);
				for (String pro : proPepMap.keySet()) {
					genomePepSet.addAll(proPepMap.get(pro));
				}
				genomePepMap.put(genome, genomePepSet);
			}

			File magsFile = new File(resultFile, "mags.tsv");
			HashSet<String> genomeSet = refineGenomes(genomePepMap, pepSet.size(), filterThreshold, maxIterativeCount,
					breakAdvance, magsFile);

			File usedMagsFile = new File(resultFile, "usedMags.tsv");
			try (PrintWriter writer = new PrintWriter(usedMagsFile)) {
				writer.println("Genome");
				for (String genome : genomeSet) {
					writer.println(genome);
				}
				writer.close();
			} catch (IOException e) {
				logger.error(getTaskName() + ": error in writing genomes to " + usedMagsFile, e);
				System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": error in writing genomes to "
						+ usedMagsFile);
			}

			return genomeSet;
		}
	}

	protected HashSet<String> refineGenomes(HashMap<String, HashSet<String>> genomePepMap, int totalPepCount,
			double threshold, int maxIteratorCount, boolean breakAdvance, File out) {

		String[] genomes = genomePepMap.keySet().toArray(new String[genomePepMap.size()]);
		Arrays.sort(genomes, (g1, g2) -> {
			int s1 = genomePepMap.get(g1).size();
			int s2 = genomePepMap.get(g2).size();
			return s2 - s1;
		});

		for (int iteratorCount = 0; iteratorCount < maxIteratorCount; iteratorCount++) {
			HashSet<String> totalSet = new HashSet<String>();
			HashMap<String, Integer> orderMap = new HashMap<String, Integer>();
			for (int i = 0; i < genomes.length; i++) {
				int currentSize = totalSet.size();
				totalSet.addAll(genomePepMap.get(genomes[i]));
				int increaseSize = totalSet.size() - currentSize;
				orderMap.put(genomes[i], increaseSize);
			}

			Arrays.sort(genomes, (g1, g2) -> {
				return orderMap.get(g2) - orderMap.get(g1);
			});
		}

		double[] percentage = new double[genomes.length];
		HashSet<String> genomeSet = new HashSet<String>();
		HashSet<String> totalSet = new HashSet<String>();
		int breakPoint1 = -1;
		int breakPoint2 = -1;

		try (PrintWriter writer = new PrintWriter(out)) {
			writer.println("Rank\tGenome\tTotal_covered_PSM_count\tPercentage\tIncrement\tIncrement ratio");
			for (int i = 0; i < genomes.length; i++) {

				int currentSize = totalSet.size();
				totalSet.addAll(genomePepMap.get(genomes[i]));
				int add = totalSet.size() - currentSize;
				double additional = ((double) add / (double) totalPepCount);

				percentage[i] = ((double) totalSet.size() / (double) totalPepCount);

				writer.println(i + "\t" + genomes[i] + "\t" + totalSet.size() + "\t" + percentage[i] + "\t" + add + "\t"
						+ additional);

				if (add > 1) {
					breakPoint1 = i;
				}
				if (additional > 0.0001) {
					breakPoint2 = i;
				}
			}
			writer.close();
		} catch (IOException e) {
			logger.error(getTaskName() + ": error in writing genomes to " + out, e);
			System.out
					.println(format.format(new Date()) + "\t" + getTaskName() + ": error in writing genomes to " + out);
		}

		for (int i = 0; i < genomes.length; i++) {
			genomeSet.add(genomes[i]);
			if (percentage[i] >= threshold) {
				break;
			}
			if (breakAdvance) {
				if (i == breakPoint1 || i == breakPoint2) {
					break;
				}
			}

			if (genomeSet.size() >= 300) {
				break;
			}
		}

		return genomeSet;
	}

	protected HashSet<String> refineProtein(HashMap<String, HashSet<String>> proPepMap, int totalPepCount, File out) {

		String[] proteins = proPepMap.keySet().toArray(new String[proPepMap.size()]);
		Arrays.sort(proteins, (g1, g2) -> {
			int s1 = proPepMap.get(g1).size();
			int s2 = proPepMap.get(g2).size();
			return s2 - s1;
		});

		for (int iteratorCount = 0; iteratorCount < 5; iteratorCount++) {
			HashSet<String> totalSet = new HashSet<String>();
			HashMap<String, Integer> orderMap = new HashMap<String, Integer>();
			for (int i = 0; i < proteins.length; i++) {
				int currentSize = totalSet.size();
				totalSet.addAll(proPepMap.get(proteins[i]));
				int increaseSize = totalSet.size() - currentSize;
				orderMap.put(proteins[i], increaseSize);
			}

			Arrays.sort(proteins, (g1, g2) -> {
				return orderMap.get(g2) - orderMap.get(g1);
			});
		}

		double[] percentage = new double[proteins.length];
		HashSet<String> totalSet = new HashSet<String>();
		int breakPoint1 = -1;
		HashSet<String> proSet = new HashSet<String>();
		try (PrintWriter writer = new PrintWriter(out)) {
			writer.println("Rank\tProtein\tTotal_covered_PSM_count\tPercentage\tIncrement\tIncrement ratio");
			for (int i = 0; i < proteins.length; i++) {

				int currentSize = totalSet.size();
				totalSet.addAll(proPepMap.get(proteins[i]));
				int add = totalSet.size() - currentSize;
				double additional = ((double) add / (double) totalPepCount);
				if (add > 0) {
					breakPoint1 = i;
				}

				percentage[i] = ((double) totalSet.size() / (double) totalPepCount);

				writer.println(i + "\t" + proteins[i] + "\t" + totalSet.size() + "\t" + percentage[i] + "\t" + add
						+ "\t" + additional);
			}
			writer.close();
		} catch (IOException e) {
			logger.error(getTaskName() + ": error in writing genomes to " + out, e);
			System.out
					.println(format.format(new Date()) + "\t" + getTaskName() + ": error in writing genomes to " + out);
		}

		for (int i = 0; i < proteins.length; i++) {
			proSet.add(proteins[i]);
			if (i == breakPoint1) {
				break;
			}
		}

		return proSet;
	}

	protected ConcurrentHashMap<String, HashMap<String, HashSet<String>>> getGenomeFromPep(HashSet<String> allPepSet,
			int threadCount) {
		ConcurrentHashMap<String, HashMap<String, HashSet<String>>> genomePepMap = new ConcurrentHashMap<String, HashMap<String, HashSet<String>>>();
		HashSet<String> matchedPepSet = new HashSet<String>();
		ExecutorService executor = Executors.newFixedThreadPool(threadCount);
		File currentFile = super.magDb.getCurrentFile();
		File predictPepFile = new File(currentFile, "predicted");
		File[] files = predictPepFile.listFiles();
		for (File file : files) {
			executor.execute(() -> {
				if (file.getName().endsWith(".predicted.tsv")) {
					try (BufferedReader reader = new BufferedReader(new FileReader(file))) {
						HashMap<String, HashSet<String>> proPepMap = new HashMap<String, HashSet<String>>();
						String pline = reader.readLine();
						while ((pline = reader.readLine()) != null) {
							String[] cs = pline.split("\t");
							if (allPepSet.contains(cs[1])) {
								matchedPepSet.add(cs[1]);
								if (proPepMap.containsKey(cs[0])) {
									proPepMap.get(cs[0]).add(cs[1]);
								} else {
									HashSet<String> pepSet = new HashSet<String>();
									pepSet.add(cs[1]);
									proPepMap.put(cs[0], pepSet);
								}
							}
						}
						if (proPepMap.size() > 0) {
							String preFileName = file.getName();
							String genome = preFileName.substring(0, preFileName.length() - ".predicted.tsv".length());
							genomePepMap.put(genome, proPepMap);
						}
					} catch (IOException e) {

						logger.error(getTaskName() + ": error in reading " + file, e);
						System.err.println(
								format.format(new Date()) + "\t" + getTaskName() + ": error in reading " + file);
					}
				}
			});
		}

		executor.shutdown();

		try {
			if (executor.awaitTermination(1, TimeUnit.HOURS)) {

				allPepSet.clear();
				allPepSet.addAll(matchedPepSet);

				logger.info(getTaskName() + ": matching peptides finished, " + allPepSet.size() + " peptides from "
						+ genomePepMap.size() + " genomes were matched");
				System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": matching peptides finished, "
						+ allPepSet.size() + " peptides from " + genomePepMap.size() + " genomes were matched");

			} else {
				logger.error(getTaskName()
						+ ": the program didn't responce in a long time, increase the thread count will be helpful");
				System.err.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": the program didn't responce in a long time, increase the thread count will be helpful");

				return null;
			}
		} catch (InterruptedException e) {

			logger.error(getTaskName() + ": error in reading predicted peptides", e);
			System.err.println(
					format.format(new Date()) + "\t" + getTaskName() + ": error in reading predicted peptides");

			return null;
		}

		return genomePepMap;
	}

	@SuppressWarnings("unused")
	private ConcurrentHashMap<String, HashSet<String>> getGenomeFromPep(HashSet<String> allPepSet,
			HashSet<String> genomeSet, int threadCount) {
		ConcurrentHashMap<String, HashSet<String>> genomePepMap = new ConcurrentHashMap<String, HashSet<String>>();
		ExecutorService executor = Executors.newFixedThreadPool(threadCount);
		File currentFile = super.magDb.getCurrentFile();
		File predictPepFile = new File(currentFile, "predicted");
		File[] files = predictPepFile.listFiles(f -> {
			String name = f.getName();
			name = name.substring(0, name.indexOf("."));
			return genomeSet.contains(name);
		});
		for (File file : files) {
			executor.execute(() -> {
				try (BufferedReader reader = new BufferedReader(new FileReader(file))) {
					HashSet<String> pepSet = new HashSet<String>();
					String pline = reader.readLine();
					while ((pline = reader.readLine()) != null) {
						String[] cs = pline.split("\t");
						if (allPepSet.contains(cs[1])) {
							pepSet.add(cs[1]);
						}
					}
					if (pepSet.size() > 0) {
						String preFileName = file.getName();
						String genome = preFileName.substring(0, preFileName.length() - ".predicted.tsv".length());
						genomePepMap.put(genome, pepSet);
					}
				} catch (IOException e) {

					logger.error(getTaskName() + ": error in reading " + file, e);
					System.err.println(format.format(new Date()) + "\t" + getTaskName() + ": error in reading " + file);
				}
			});
		}

		executor.shutdown();

		try {
			if (executor.awaitTermination(1, TimeUnit.HOURS)) {
				logger.info(getTaskName() + ": matching peptides finished, " + allPepSet.size() + " peptides from "
						+ genomePepMap.size() + " genomes were matched");
				System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": matching peptides finished, "
						+ allPepSet.size() + " peptides from " + genomePepMap.size() + " genomes were matched");

			} else {
				logger.error(getTaskName()
						+ ": the program didn't responce in a long time, increase the thread count will be helpful");
				System.err.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": the program didn't responce in a long time, increase the thread count will be helpful");

				return null;
			}
		} catch (InterruptedException e) {

			logger.error(getTaskName() + ": error in reading predicted peptides", e);
			System.err.println(
					format.format(new Date()) + "\t" + getTaskName() + ": error in reading predicted peptides");

			return null;
		}

		return genomePepMap;
	}

	@Override
	protected String getTaskName() {
		return taskName;
	}

	/**
	 * Thread-safe method to update the main progress bar
	 */
	private void updateMainProgress(double progress) {
		if (Platform.isFxApplicationThread()) {
			progressBar.setProgress(progress);
		} else {
			Platform.runLater(() -> progressBar.setProgress(progress));
		}
	}

	/**
	 * Thread-safe method to update the secondary progress bar
	 */
	protected void updateSecondaryProgress(double progress) {
		if (Platform.isFxApplicationThread()) {
			secondBar.setProgress(progress);
		} else {
			Platform.runLater(() -> secondBar.setProgress(progress));
		}
	}

	@Override
	protected Boolean call() throws Exception {
		updateMessage(getTaskName());
		updateProgress(0, 100);

		logger.info(getTaskName() + ": started");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": started");

		// Log detailed parameter information
		logger.info("Workflow parameters:");
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
		if (magDb != null) {
			logger.info("  MAG database: " + magDb.getCatalogueID());
			logger.info("  MAG database version: "
					+ (magDb.getUsedVersion() != null ? magDb.getUsedVersion().getCatalogueVersion() : "unknown"));
		}

		// Log resource information
		logger.info("Resource parameters:");
		logger.info("  DiaNN path: " + mrp.getToolPath(MetaResourceParameters.ToolType.DIANN));
		logger.info("  MetaResourceParameters object: {}", mrp);
		logger.info("  MetaResourceParameters class: {}", mrp.getClass().getName());

		// Update main progress and keep second bar active
		updateProgress(5, 100);
		updateSecondaryProgress(ProgressBar.INDETERMINATE_PROGRESS); // Indeterminate progress for continuous activity

		// Step 0: Copy report templates in background
		if (!copyReportTemplatesInBackground()) {
			logger.error(getTaskName() + ": Failed to copy report templates");
			return false;
		}

		// Step 1: Spectra Conversion
		logger.info("Step 1: Starting spectra format conversion...");
		updateMessage(getTaskName() + ": Converting spectra format...");
		updateProgress(10, 100);
		updateSecondaryProgress(ProgressBar.INDETERMINATE_PROGRESS);

		boolean convert = spectraConvert();
		if (!convert) {
			logger.error(getTaskName() + ": error in spectra format conversion");
			System.err
					.println(format.format(new Date()) + "\t" + getTaskName() + ": error in spectra format conversion");
			System.err.println("ERROR: cannot read .raw files, please download and install .NET SDK 8.0.407 or later");
			System.err.println("https://dotnet.microsoft.com/en-us/download/dotnet/8.0");
			// Update progress on error
			updateMessage(getTaskName() + ": Failed - Spectra conversion error");
			updateSecondaryProgress(0.0); // Stop the indeterminate progress
			return false;
		}
		logger.info("Step 1: Spectra format conversion completed successfully");

		// Step 2: sample-specific library generation
		logger.info("Step 2: Starting sample-specific library generation...");
		updateProgress(15, 100);
		updateMessage(getTaskName() + ": Creating library...");
		updateSecondaryProgress(ProgressBar.INDETERMINATE_PROGRESS);
		updateMessage("Creating sample-specific library...");

		System.out.println(
				format.format(new Date()) + "\t" + getTaskName() + ": creating sample-specific library started");

		boolean createSSDB = createSSLibrary();
		if (!createSSDB) {
			logger.error(getTaskName() + ": error in creating sample-specific library");
			System.err.println(
					format.format(new Date()) + "\t" + getTaskName() + ": error in creating sample-specific library");
			// Update progress on error
			updateMessage(getTaskName() + ": Failed - sample-specific library creation error");
			updateSecondaryProgress(0.0); // Stop the indeterminate progress
			return false;
		}
		logger.info("Step 2: sample-specific library creation completed successfully");
		System.out.println(
				format.format(new Date()) + "\t" + getTaskName() + ": creating sample-specific library finished");

		// Step 3: Identification
		logger.info("Step 3: Starting identification...");
		updateProgress(20, 100);
		updateMessage(getTaskName() + ": Running identification...");
		updateSecondaryProgress(ProgressBar.INDETERMINATE_PROGRESS);
		updateMessage("Performing protein identification...");

		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": identification started");

		boolean iden = iden();
		if (!iden) {
			logger.error(getTaskName() + ": error in identification");
			System.err.println(format.format(new Date()) + "\t" + getTaskName() + ": error in identification");
			// Update progress on error
			updateMessage(getTaskName() + ": Failed - Identification error");
			updateSecondaryProgress(0.0); // Stop the indeterminate progress
			return false;
		}
		logger.info("Step 3: Identification completed successfully");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": identification finished");

		// Step 4: Export Report
		logger.info("Step 4: Starting report generation...");
		updateProgress(85, 100);
		updateMessage(getTaskName() + ": Generating reports...");
		updateSecondaryProgress(ProgressBar.INDETERMINATE_PROGRESS);
		updateMessage("Preparing final reports...");

		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": report generation started");

		boolean exportReport = exportReport();
		if (!exportReport) {
			logger.error(getTaskName() + ": error in report generation");
			System.err.println(format.format(new Date()) + "\t" + getTaskName() + ": error in report generation");
			// Update progress on error
			updateMessage(getTaskName() + ": Failed - Report generation error");
			updateSecondaryProgress(0.0); // Stop the indeterminate progress
			return false;
		}
		logger.info("Step 4: Report generation completed successfully");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": report generation finished");

		// Update progress for completion
		updateProgress(100, 100);
		updateMessage(getTaskName() + ": Completed successfully!");
		updateSecondaryProgress(0.0); // Stop the indeterminate progress
		updateMessage("Workflow completed!");

		logger.info(getTaskName() + ": finished - all workflow steps completed successfully");
		System.out.println(format.format(new Date()) + "\t" + getTaskName()
				+ ": finished - all workflow steps completed successfully");
		return true;
	}

	@Override
	public void forceStop() {
		logger.info(getTaskName() + ": force stop requested");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": force stop requested");

		// Cancel DiaNN task if running
		if (diaNNTask != null && !diaNNTask.isDone()) {
			logger.info(getTaskName() + ": Cancelling DiaNN task");
			System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": Cancelling DiaNN task");
			diaNNTask.cancel();
		}

		// Cancel this task itself
		this.cancel(true);

		logger.info(getTaskName() + ": force stop completed");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": force stop completed");
	}

	protected class ProteinGroup {
		@SuppressWarnings("unused")
		String proName;
		HashSet<String> sameSet;
		HashSet<String> subSet;

		ProteinGroup(String proName) {
			this.proName = proName;
			this.sameSet = new HashSet<String>();
			this.subSet = new HashSet<String>();
		}

		public void addSameGroup(String pro) {
			this.sameSet.add(pro);
		}

		public HashSet<String> getSameSet() {
			return sameSet;
		}
	}

	/**
	 * Static method to copy first search results from input folder to output folder
	 * This is useful for testing different workflows where the initial search
	 * results are identical
	 * 
	 * @param inputFolder  The source folder containing the first search results
	 * @param outputFolder The destination folder where results will be copied
	 * @return true if copy was successful, false otherwise
	 */
	public static boolean copyFirstSearchResults(File inputFolder, File outputFolder) {
		if (inputFolder == null || outputFolder == null) {
			System.err.println("Error: Input or output folder is null");
			return false;
		}

		if (!inputFolder.exists() || !inputFolder.isDirectory()) {
			System.err.println(
					"Error: Input folder does not exist or is not a directory: " + inputFolder.getAbsolutePath());
			return false;
		}

		// Create output folder if it doesn't exist
		if (!outputFolder.exists()) {
			if (!outputFolder.mkdirs()) {
				System.err.println("Error: Failed to create output folder: " + outputFolder.getAbsolutePath());
				return false;
			}
		}

		System.out.println("Copying first search results from: " + inputFolder.getAbsolutePath());
		System.out.println("To: " + outputFolder.getAbsolutePath());

		boolean success = true;

		try {
			// Copy all separate result folders (one for each DIA file)
			File[] inputFiles = inputFolder.listFiles();
			if (inputFiles != null) {
				for (File inputFile : inputFiles) {
					if (inputFile.isDirectory()) {
						// Check if this folder contains a "firstSearch" subfolder
						File firstSearchInput = new File(inputFile, "firstSearch");
						if (firstSearchInput.exists() && firstSearchInput.isDirectory()) {
							// Only process folders that contain firstSearch
							File outputFile = new File(outputFolder, inputFile.getName());

							// Create the output directory
							if (!outputFile.exists()) {
								outputFile.mkdirs();
							}

							File firstSearchOutput = new File(outputFile, "firstSearch");
							if (!firstSearchOutput.exists()) {
								firstSearchOutput.mkdirs();
							}

							// Copy firstSearch.parquet
							File parquetInput = new File(firstSearchInput, "firstSearch.parquet");
							File parquetOutput = new File(firstSearchOutput, "firstSearch.parquet");
							if (parquetInput.exists()) {
								Files.copy(parquetInput.toPath(), parquetOutput.toPath(),
										StandardCopyOption.REPLACE_EXISTING);
								System.out.println("Copied: " + parquetInput.getAbsolutePath() + " -> "
										+ parquetOutput.getAbsolutePath());
							}

							// Copy all other files in firstSearch folder (excluding firstSearch.fasta)
							File[] firstSearchFiles = firstSearchInput.listFiles();
							if (firstSearchFiles != null) {
								for (File file : firstSearchFiles) {
									if (file.isFile() && !file.getName().equals("firstSearch.fasta")) {
										File outputFile1 = new File(firstSearchOutput, file.getName());
										Files.copy(file.toPath(), outputFile1.toPath(),
												StandardCopyOption.REPLACE_EXISTING);
										System.out.println("Copied: " + file.getAbsolutePath() + " -> "
												+ outputFile1.getAbsolutePath());
									}
								}
							}
						} else {
							// Skip folders that don't contain firstSearch
							System.out.println("Skipping folder (no firstSearch found): " + inputFile.getName());
						}
					}
				}
			}

			System.out.println("First search results copied successfully!");

		} catch (IOException e) {
			System.err.println("Error copying files: " + e.getMessage());
			e.printStackTrace();
			success = false;
		}

		return success;
	}

	/**
	 * Generate bubble chart report from tree.js data using Python.
	 */
	protected void generateBubbleChartReport(File treeJsFile, File outputFile) {
		try {
			logger.info(taskName + ": Generating bubble chart report from " + treeJsFile.getName());
			System.out.println(format.format(new Date()) + "\t" + taskName + ": Generating bubble chart report");

			// Use MetaPythonReportTask to generate the bubble chart
			MetaPythonReportTask reportTask = new MetaPythonReportTask(1, mrp, logger);
			boolean success = reportTask.addTask(MetaPythonReportTask.BUBBLE_CHART, treeJsFile, outputFile);

			if (success) {
				// Wait for completion
				reportTask.run();
				boolean result = reportTask.get();

				if (result) {
					logger.info(taskName + ": Bubble chart report generated successfully: " + outputFile.getName());
					System.out.println(format.format(new Date()) + "\t" + taskName
							+ ": Bubble chart report generated successfully");

					// Add to manifest
					File reportsDir = new File(metaPar.getReportDir(), "reports");
					ReportManifestGenerator manifestGen = new ReportManifestGenerator(reportsDir);
					manifestGen.addReportEntry("Bubble Chart Report", "bubble_chart_report.html",
							"Interactive hierarchical bubble chart showing taxonomic abundance with treemap, sunburst, and distribution visualizations.",
							new String[] { "bubble-chart", "taxonomy", "hierarchical", "treemap" });
				} else {
					logger.error(taskName + ": Failed to generate bubble chart report");
					System.err.println(
							format.format(new Date()) + "\t" + taskName + ": Failed to generate bubble chart report");
				}
			} else {
				logger.error(taskName + ": Failed to add bubble chart task");
				System.err.println(format.format(new Date()) + "\t" + taskName + ": Failed to add bubble chart task");
			}
		} catch (Exception e) {
			logger.error(taskName + ": Error generating bubble chart report", e);
			System.err.println(format.format(new Date()) + "\t" + taskName + ": Error generating bubble chart report: "
					+ e.getMessage());
		}
	}

	/**
	 * Index in score array: S, p_emp, q_value. Used by new genome/protein scoring
	 * scheme.
	 */
	protected static final int IDX_S = 0, IDX_P_EMP = 1, IDX_Q_VALUE = 2;

	/**
	 * Folder containing combined.fasta for new-scheme export. Default:
	 * result/combined. Override in subclasses (e.g. MetaDiaMarkerTask) to use
	 * iterative result folder.
	 */
	protected File getExportFastaFolder() {
		return new File(metaPar.getDbSearchResultFile(), "combined");
	}

	/**
	 * Folder containing combined.parquet and combined.pr_matrix.tsv for new-scheme
	 * export. Default: result/combined. Override in subclasses (e.g.
	 * MetaDiaMarkerTask) to use combinedFile.
	 */
	protected File getExportParquetFolder() {
		return new File(metaPar.getDbSearchResultFile(), "combined");
	}

	/**
	 * Used by computeGenomeScoreMap / computeProteinScoreMap (new scoring scheme).
	 */
	protected static class PeptideEvidence {
		String seq;
		double bestCScore = 0;
		double bestPEP = 1.0;
		double bestPreQuan = 0;
		double bestMs1Area = 0;
		boolean isDecoy;
		HashSet<String> proteins = new HashSet<>();
	}

	/** Used by computeGenomeScoreMap (new scoring scheme). */
	protected static class WeightedPeptideScore {
		String seq;
		double e;
		double w;
		double we;
		double intensityFactor = 1.0;
	}

	/** Used by computeGenomeScoreMap (new scoring scheme). */
	protected static class GenomePeptideEntry {
		ArrayList<WeightedPeptideScore> weightedScores = new ArrayList<>();
	}

	/**
	 * Compute per-genome scores (S, p_emp, q_value) from precursors using
	 * PEP/CScore+intensity scoring and decoy resampling. Uses unique peptides only
	 * when uniquePeptidesOnly is true. Value is double[]{S, p_emp, q_value};
	 * q_value is Benjamini–Hochberg corrected.
	 */
	protected HashMap<String, double[]> computeGenomeScoreMap(DiaNNPrecursor[] precursors, boolean uniquePeptidesOnly) {
		HashMap<String, double[]> genomeScoreMap = new HashMap<>();
		if (precursors == null || precursors.length == 0) {
			return genomeScoreMap;
		}
		HashMap<String, PeptideEvidence> peptideMap = new HashMap<>();
		for (DiaNNPrecursor pp : precursors) {
			String seq = pp.getSeqString();
			if (seq == null || seq.isEmpty())
				continue;
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
				ev.bestCScore = Math.max(ev.bestCScore, cscore);
				ev.bestPEP = Math.min(ev.bestPEP, pep);
				ev.bestPreQuan = Math.max(ev.bestPreQuan, preQuan);
				ev.bestMs1Area = Math.max(ev.bestMs1Area, ms1Area);
			}
			String[] proGroups = pp.getProGroups();
			if (proGroups != null) {
				for (String pro : proGroups) {
					if (pro != null && !pro.isEmpty())
						ev.proteins.add(pro);
				}
			}
		}
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
		final double decoyPepMax = 0.99;
		ArrayList<Double> decoyEPool = new ArrayList<>();
		HashMap<String, GenomePeptideEntry> genomePeptides = new HashMap<>();
		for (Map.Entry<String, PeptideEvidence> ent : peptideMap.entrySet()) {
			PeptideEvidence ev = ent.getValue();
			String seq = ev.seq;
			double e = ev.bestCScore > 0 ? Math.min(cap, Math.max(0, ev.bestCScore))
					: Math.min(cap, Math.max(0, -Math.log10(Math.max(ev.bestPEP, 1e-10))));
			double logI1 = Math.log(1.0 + Math.max(0, ev.bestPreQuan)) / Math.log(2.0);
			double logI2 = Math.log(1.0 + Math.max(0, ev.bestMs1Area)) / Math.log(2.0);
			double norm1 = (logI1 - minLog1) / logRange1;
			double norm2 = (logI2 - minLog2) / logRange2;
			double intensityFactor = 0.5 + 0.25 * norm1 + 0.25 * norm2;
			if (ev.isDecoy) {
				if (ev.bestPEP < decoyPepMax)
					decoyEPool.add(e * intensityFactor);
				continue;
			}
			HashSet<String> genomes = new HashSet<>();
			for (String pro : ev.proteins) {
				if (pro.startsWith(FastaManager.shuffle))
					continue;
				String[] parts = pro.split("_");
				if (parts.length > 0 && !parts[0].isEmpty())
					genomes.add(parts[0]);
			}
			if (genomes.isEmpty())
				continue;
			int m = genomes.size();
			if (m != 1 && uniquePeptidesOnly)
				continue;
			double w = (m == 1) ? 1.0 : (1.0 / m);
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
		}
		int k = 5;
		Random rnd = new Random(42);
		int B = 10000;
		HashMap<String, Double> SMap = new HashMap<>();
		HashMap<String, Integer> genomeNMap = new HashMap<>();
		HashMap<String, Integer> genomeTotalPeptidesMap = new HashMap<>();
		for (String g : genomePeptides.keySet()) {
			List<WeightedPeptideScore> scores = genomePeptides.get(g).weightedScores;
			genomeTotalPeptidesMap.put(g, scores.size());
			scores.sort((a, b) -> Double.compare(b.we, a.we));
			int n = Math.min(k, scores.size());
			double S = 0;
			for (int i = 0; i < n; i++)
				S += scores.get(i).we;
			SMap.put(g, S);
			genomeNMap.put(g, n);
		}
		if (decoyEPool.isEmpty())
			decoyEPool.add(0.0);
		double[] decoyArr = new double[decoyEPool.size()];
		for (int i = 0; i < decoyArr.length; i++)
			decoyArr[i] = decoyEPool.get(i);
		int decoyLen = decoyArr.length;
		int maxM = 0;
		for (Integer m : genomeTotalPeptidesMap.values())
			maxM = Math.max(maxM, m);
		double[] tmp = new double[Math.max(maxM, 1)];
		HashMap<String, Double> pempMap = new HashMap<>();
		for (String g : SMap.keySet()) {
			double Sg = SMap.get(g);
			int M = genomeTotalPeptidesMap.get(g);
			int ng = genomeNMap.get(g);
			int count = 0;
			for (int b = 0; b < B; b++) {
				for (int i = 0; i < M; i++)
					tmp[i] = decoyArr[rnd.nextInt(decoyLen)];
				Arrays.sort(tmp, 0, M);
				double Sb = 0;
				for (int i = M - 1; i >= M - ng && i >= 0; i--)
					Sb += tmp[i];
				if (Sb >= Sg)
					count++;
			}
			double pemp = (1.0 + count) / (1.0 + B);
			pempMap.put(g, pemp);
		}
		String[] genomes = pempMap.keySet().toArray(new String[0]);
		Arrays.sort(genomes, Comparator.comparingDouble(pempMap::get));
		int nGenomes = genomes.length;
		for (int i = 0; i < nGenomes; i++) {
			double pemp = pempMap.get(genomes[i]);
			double q = Math.min(1.0, pemp * nGenomes / (i + 1));
			genomeScoreMap.put(genomes[i], new double[] { SMap.get(genomes[i]), pemp, q });
		}
		return genomeScoreMap;
	}

	/**
	 * Compute per-protein scores (S, p_emp, q_value) from precursors using
	 * PEP/CScore+intensity scoring and decoy resampling. S = sum of ALL weighted
	 * scores per protein (no top-k).
	 */
	protected HashMap<String, double[]> computeProteinScoreMap(DiaNNPrecursor[] precursors) {
		HashMap<String, double[]> proteinScoreMap = new HashMap<>();
		if (precursors == null || precursors.length == 0) {
			return proteinScoreMap;
		}
		HashMap<String, PeptideEvidence> peptideMap = new HashMap<>();
		for (DiaNNPrecursor pp : precursors) {
			String seq = pp.getSeqString();
			if (seq == null || seq.isEmpty())
				continue;
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
				ev.bestCScore = Math.max(ev.bestCScore, cscore);
				ev.bestPEP = Math.min(ev.bestPEP, pep);
				ev.bestPreQuan = Math.max(ev.bestPreQuan, preQuan);
				ev.bestMs1Area = Math.max(ev.bestMs1Area, ms1Area);
			}
			String[] proGroups = pp.getProGroups();
			if (proGroups != null) {
				for (String pro : proGroups) {
					if (pro != null && !pro.isEmpty())
						ev.proteins.add(pro);
				}
			}
		}
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
		final double decoyPepMax = 0.99;
		ArrayList<Double> decoyEPool = new ArrayList<>();
		HashMap<String, ArrayList<Double>> proteinPeptides = new HashMap<>();
		for (Map.Entry<String, PeptideEvidence> ent : peptideMap.entrySet()) {
			PeptideEvidence ev = ent.getValue();
			String seq = ev.seq;
			double e = ev.bestCScore > 0 ? Math.min(cap, Math.max(0, ev.bestCScore))
					: Math.min(cap, Math.max(0, -Math.log10(Math.max(ev.bestPEP, 1e-10))));
			double logI1 = Math.log(1.0 + Math.max(0, ev.bestPreQuan)) / Math.log(2.0);
			double logI2 = Math.log(1.0 + Math.max(0, ev.bestMs1Area)) / Math.log(2.0);
			double norm1 = (logI1 - minLog1) / logRange1;
			double norm2 = (logI2 - minLog2) / logRange2;
			double intensityFactor = 0.5 + 0.25 * norm1 + 0.25 * norm2;
			double wei = e * intensityFactor;
			if (ev.isDecoy) {
				if (ev.bestPEP < decoyPepMax)
					decoyEPool.add(wei);
				for (String pro : ev.proteins) {
					if (!pro.startsWith(FastaManager.shuffle))
						continue;
					ArrayList<Double> weiList = proteinPeptides.get(pro);
					if (weiList == null) {
						weiList = new ArrayList<>();
						proteinPeptides.put(pro, weiList);
					}
					weiList.add(wei);
				}
				continue;
			}
			for (String pro : ev.proteins) {
				if (pro.startsWith(FastaManager.shuffle))
					continue;
				ArrayList<Double> weiList = proteinPeptides.get(pro);
				if (weiList == null) {
					weiList = new ArrayList<>();
					proteinPeptides.put(pro, weiList);
				}
				weiList.add(wei);
			}
		}
		Random rnd = new Random(42);
		int B = 10000;
		HashMap<String, Double> SMap = new HashMap<>();
		HashMap<String, Integer> proteinNMap = new HashMap<>();
		for (String pro : proteinPeptides.keySet()) {
			ArrayList<Double> weiList = proteinPeptides.get(pro);
			double S = 0;
			for (double w : weiList)
				S += w;
			SMap.put(pro, S);
			proteinNMap.put(pro, weiList.size());
		}
		if (decoyEPool.isEmpty())
			decoyEPool.add(0.0);
		double[] decoyArr = new double[decoyEPool.size()];
		for (int i = 0; i < decoyArr.length; i++)
			decoyArr[i] = decoyEPool.get(i);
		int decoyLen = decoyArr.length;
		int maxM = 0;
		for (Integer m : proteinNMap.values())
			maxM = Math.max(maxM, m);
		double[] tmp = new double[Math.max(maxM, 1)];
		HashMap<String, Double> pempMap = new HashMap<>();
		for (String pro : SMap.keySet()) {
			double Sp = SMap.get(pro);
			int M = proteinNMap.get(pro);
			int count = 0;
			for (int b = 0; b < B; b++) {
				double Sb = 0;
				for (int i = 0; i < M; i++) {
					tmp[i] = decoyArr[rnd.nextInt(decoyLen)];
					Sb += tmp[i];
				}
				if (Sb >= Sp)
					count++;
			}
			pempMap.put(pro, (1.0 + count) / (1.0 + B));
		}
		String[] proteins = pempMap.keySet().toArray(new String[0]);
		Arrays.sort(proteins, Comparator.comparingDouble(pempMap::get));
		int nProteins = proteins.length;
		for (int i = 0; i < nProteins; i++) {
			double pemp = pempMap.get(proteins[i]);
			double q = Math.min(1.0, pemp * nProteins / (i + 1));
			proteinScoreMap.put(proteins[i], new double[] { SMap.get(proteins[i]), pemp, q });
		}
		return proteinScoreMap;
	}

	/**
	 * Core DIA export implementation (new genome/protein scoring: S, p_emp,
	 * q_value from PEP/CScore+intensity, decoy resampling). Uses
	 * {@link #getExportFastaFolder()} and {@link #getExportParquetFolder()} for
	 * paths. Called by {@link #exportReport()} and shared across DIA workflows.
	 */
	protected boolean doExportReportInternal() {

		// Prefer workflow-set peptide matrix file when already set and valid
		if (this.quan_pep_file != null && this.quan_pep_file.exists() && this.quan_pep_file.canRead()) {
			// use existing quan_pep_file; parquet folder same parent for parquet path below
		} else {
			File parquetFolder = getExportParquetFolder();
			if (parquetFolder == null && metaPar != null) {
				parquetFolder = new File(metaPar.getDbSearchResultFile(), "combined");
				logger.warn(getTaskName() + ": getExportParquetFolder() was null, using fallback {}",
						parquetFolder.getAbsolutePath());
			}
			this.quan_pep_file = new File(parquetFolder != null ? parquetFolder : new File("combined"),
					"combined.pr_matrix.tsv");
		}
		if (!this.quan_pep_file.exists() || !this.quan_pep_file.canRead()) {
			logger.error(getTaskName() + ": peptide matrix file missing or unreadable: {}",
					this.quan_pep_file.getAbsolutePath());
			return false;
		}

		File parquetFolder = this.quan_pep_file.getParentFile();

		logger.info(getTaskName() + ": protein inference started (new genome scoring)");
		System.out.println(new java.text.SimpleDateFormat("MM_dd_yyyy_HH_mm_ss").format(new Date()) + "\t"
				+ getTaskName() + ": protein inference started");

		File parquetFile = new File(parquetFolder != null ? parquetFolder : new File("combined"), "combined.parquet");
		if (!parquetFile.exists() || !parquetFile.canRead()) {
			logger.error(getTaskName() + ": parquet file missing or unreadable: {}", parquetFile.getAbsolutePath());
			return false;
		}
		DiaNNParquetReader diannReader = new DiaNNParquetReader(parquetFile.getAbsolutePath(), this.diannPar.getCut());
		DiaNNPrecursor[] precursorsQ01 = diannReader.getDiaNNPrecursors();
		if (precursorsQ01 == null)
			precursorsQ01 = new DiaNNPrecursor[0];
		precursorsQ01 = mapPrecursorsToGenomes(precursorsQ01);
		if (precursorsQ01 == null)
			precursorsQ01 = new DiaNNPrecursor[0];
		// Allow subclasses (e.g. deep marker/seed) to enrich ProGroups (peptide→protein mapping)
		// from their own intermediate round outputs before scoring.
		precursorsQ01 = applyPostParquetTransforms(precursorsQ01, parquetFolder);
		if (precursorsQ01 == null)
			precursorsQ01 = new DiaNNPrecursor[0];
		HashMap<String, double[]> genomeScoreMapForFilter = computeGenomeScoreMap(precursorsQ01, true);
		HashSet<String> genomesPempOne = new HashSet<>();
		for (String g : genomeScoreMapForFilter.keySet()) {
			if (genomeScoreMapForFilter.get(g)[IDX_P_EMP] >= 1.0 - 1e-9)
				genomesPempOne.add(g);
		}
		ArrayList<DiaNNPrecursor> retainedList = new ArrayList<>();
		for (DiaNNPrecursor pp : precursorsQ01) {
			if (pp.getGlobalQvalue() >= 0.01)
				continue;
			String[] proGroups = pp.getProGroups();
			if (proGroups == null || proGroups.length == 0)
				continue;
			ArrayList<String> keptPro = new ArrayList<>();
			for (String pro : proGroups) {
				if (pro == null || pro.isEmpty())
					continue;
				if (pro.startsWith(FastaManager.shuffle)) {
					keptPro.add(pro);
					continue;
				}
				String genome = pro.split("_").length > 0 ? pro.split("_")[0] : null;
				if (genome != null && !genomesPempOne.contains(genome))
					keptPro.add(pro);
			}
			if (!keptPro.isEmpty()) {
				pp.setProGroups(keptPro.toArray(new String[0]));
				retainedList.add(pp);
			}
		}
		this.diaPp = retainedList.toArray(new DiaNNPrecursor[0]);
		HashMap<String, double[]> genomeScoreMap = computeGenomeScoreMap(this.diaPp, true);
		HashMap<String, double[]> proteinScoreMap = computeProteinScoreMap(this.diaPp);
		ArrayList<DiaNNPrecursor> finalFilteredList = new ArrayList<>();
		for (DiaNNPrecursor pp : this.diaPp) {
			ArrayList<String> proList = new ArrayList<>();
			String[] pros = pp.getProGroups();
			String razorProtein = "";
			double bestS = -Double.MAX_VALUE;
			for (String pro : pros) {
				if (proteinScoreMap.containsKey(pro)) {
					proList.add(pro);
					double S = proteinScoreMap.get(pro)[IDX_S];
					if (S > bestS) {
						bestS = S;
						razorProtein = pro;
					}
				}
			}
			if (!proList.isEmpty()) {
				pp.setProGroups(proList.toArray(new String[0]));
				pp.setRazorProtein(razorProtein);
				finalFilteredList.add(pp);
			}
		}
		this.diaPp = finalFilteredList.toArray(new DiaNNPrecursor[0]);
		// Delegate the rest of report generation (peptide report, HTML, etc.) to the
		// existing export flow in subclasses if they override it; otherwise run the
		// same downstream as Marker.
		return exportReportDownstreamFromNewScoring(genomeScoreMap, proteinScoreMap);
	}

	/**
	 * Getters/setters for use by static downstream report (e.g.
	 * runFullDownstreamFromNewScoring).
	 */
	protected org.apache.logging.log4j.Logger getExportLogger() {
		return logger;
	}

	/**
	 * Peptide matrix file for report export. DIA workflows must set this (e.g. in
	 * iden/combine) before report generation.
	 * 
	 * @return the combined peptide matrix file (e.g.
	 *         combined/combined.pr_matrix.tsv), or null if not set
	 */
	protected File getExportQuanPepFile() {
		return quan_pep_file;
	}

	protected String[] getExportQuanExpNames() {
		return quanExpNames;
	}

	protected void setExportQuanExpNames(String[] v) {
		this.quanExpNames = v;
	}

	protected File getExportResultFolderFile() {
		return resultFolderFile;
	}

	protected MagDbItem getExportMagDb() {
		return magDb;
	}

	protected File getExportFinalPepTxt() {
		return final_pep_txt;
	}

	protected void setExportFinalPepTxt(File f) {
		this.final_pep_txt = f;
	}

	protected File getExportFinalProTxt() {
		return final_pro_txt;
	}

	protected void setExportFinalProTxt(File f) {
		this.final_pro_txt = f;
	}

	protected File getExportReportHtmlDir() {
		return reportHtmlDir;
	}

	protected void setExportReportHtmlDir(File f) {
		this.reportHtmlDir = f;
	}

	protected File getExportFinalProXmlFile() {
		return final_pro_xml_file;
	}

	protected void setExportFinalProXmlFile(File f) {
		this.final_pro_xml_file = f;
	}

	protected File getExportSummaryMetaFile() {
		return summaryMetaFile;
	}

	protected void setExportSummaryMetaFile(File f) {
		this.summaryMetaFile = f;
	}

	protected MetaParameterDia getExportMetaPar() {
		return (MetaParameterDia) metaPar;
	}

	protected MetaResourceParameters getExportMrp() {
		return mrp;
	}

	protected String[] getExportExpNames() {
		return expNames;
	}

	protected MetaData getExportMetadata() {
		return metadata;
	}

	protected void exportUpdateProgress(long a, long b) {
		updateProgress(a, b);
	}

	protected void exportUpdateMessage(String msg) {
		updateMessage(msg);
	}

	protected MetaProteinAnnoMag[] exportFuncAnnoSqlite(MetaProtein[] pros, String[] fileNames) {
		return funcAnnoSqlite(pros, fileNames);
	}

	/**
	 * Optional mapping of precursors to genome proteins (e.g. from MAG predicted
	 * TSVs). Default: return precursors unchanged. Override in Seed Direct to fill
	 * proGroups from predicted folder.
	 */
	protected DiaNNPrecursor[] mapPrecursorsToGenomes(DiaNNPrecursor[] precursors) {
		return precursors;
	}

	/**
	 * Hook to apply additional transformations to precursors after reading the combined parquet
	 * and optional genome mapping, but before genome/protein scoring and report export.
	 * Default: invoke {@link #enrichPrecursorsWithPeptideProteinMapping(DiaNNPrecursor[], File)}
	 * and return the same array.
	 */
	protected DiaNNPrecursor[] applyPostParquetTransforms(DiaNNPrecursor[] precursors, File parquetFolder) {
		enrichPrecursorsWithPeptideProteinMapping(precursors, parquetFolder);
		return precursors;
	}

	/**
	 * Optional hook for subclasses (e.g. deep marker/seed) to rebuild peptide→protein mapping
	 * from their own iterative round outputs. Default: no-op.
	 */
	protected void enrichPrecursorsWithPeptideProteinMapping(DiaNNPrecursor[] precursors, File parquetFolder) {
		// default: do nothing
	}

	/**
	 * Downstream report generation after new genome/protein scores (peptide TSV,
	 * intensity maps, HTML). Default: run the same logic as MetaDiaMarkerTask
	 * (shared implementation). Override to customize.
	 */
	protected boolean exportReportDownstreamFromNewScoring(HashMap<String, double[]> genomeScoreMap,
			HashMap<String, double[]> proteinScoreMap) {
		// Default: call the base exportReport() which does the old scheme - not what we
		// want. Subclasses that use new scoring (Marker, SeedDirect) override this.
		return false;
	}
}
