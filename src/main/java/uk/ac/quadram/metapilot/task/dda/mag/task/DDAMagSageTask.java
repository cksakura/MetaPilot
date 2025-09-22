package uk.ac.quadram.metapilot.task.dda.mag.task;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileFilter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.sql.SQLException;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.Arrays;
import java.util.Comparator;
import java.util.concurrent.ExecutionException;
import java.math.BigInteger;
import java.text.DecimalFormat;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import uk.ac.quadram.metapilot.dbSearch.sage.SageTask;
import uk.ac.quadram.metapilot.quant.flashLFQ.FlashLfqTask;
import uk.ac.quadram.metapilot.core.prodb.FastaManager;
import uk.ac.quadram.metapilot.core.tools.FormatTool;
import uk.ac.quadram.metapilot.dbSearch.sage.SagePeptide;
import uk.ac.quadram.metapilot.dbSearch.sage.SageResultReader;
import uk.ac.quadram.metapilot.task.io.pep.MetaPeptide;
import uk.ac.quadram.metapilot.task.io.pro.MetaProtein;
import uk.ac.quadram.metapilot.task.io.pro.MetaProteinXMLReader2;
import uk.ac.quadram.metapilot.task.mag.HMGenomeHtmlWriter;
import uk.ac.quadram.metapilot.task.mag.MagFuncSearcher;
import uk.ac.quadram.metapilot.task.mag.MetaProteinAnnoMag;
import uk.ac.quadram.metapilot.task.MetaPythonReportTask;
import uk.ac.quadram.metapilot.task.io.MetaTreeHandler;
import uk.ac.quadram.metapilot.core.taxonomy.TaxonomyRanks;
import org.apache.commons.lang3.math.NumberUtils;
import org.dom4j.DocumentException;
import uk.ac.quadram.metapilot.task.par.MetaConstants;
import javafx.scene.control.ProgressBar;
import javafx.scene.control.Label;
import javafx.application.Platform;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters.ToolType;
import uk.ac.quadram.metapilot.task.mag.MetaPilotMagTask;
import uk.ac.quadram.metapilot.task.dda.mag.par.DDAMagSageParameters;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;
import uk.ac.quadram.metapilot.spectra.io.MSConvertor;

public class DDAMagSageTask extends MetaPilotMagTask {

	private String[] mzMLFiles;
	private File sageExe;
	private File msConvertExe;
	private SageTask sageTask;
	private FlashLfqTask flashLfqTask;
	private MetaPythonReportTask reportTask;
	private Thread reportCopyThread;
	private Thread reportThread;
	private ExecutorService executorService;

	// Progress tracking components
	private Label statusLabel;
	
	// Add direct reference to DDAMagSageParameters
	private DDAMagSageParameters magSageParams;

	private static final Logger LOGGER = LogManager.getLogger(DDAMagSageTask.class);
	private static final String taskName = "DDA MAG database search by Sage workflow";

	public DDAMagSageTask(DDAMagSageParameters magSageParams, MetaResourceParameters mrp, ProgressBar progressBar,
			ProgressBar secondBar) {
		super(magSageParams, mrp, progressBar, secondBar);
		
		// Store the direct reference to DDAMagSageParameters
		this.magSageParams = magSageParams;
		
		// Initialize status label from parameters
		this.statusLabel = magSageParams.getStatusLabel();
		
		// Initialize the magDb field from parameters
		MagDbItem selectedMagDb = magSageParams.getMagDbItem();
		if (selectedMagDb != null) {
			this.magDb = selectedMagDb;
			LOGGER.info("Initialized MAG database: {} (version: {})", selectedMagDb.getCatalogueID(),
					selectedMagDb.getVersionList() != null && !selectedMagDb.getVersionList().isEmpty()
							? selectedMagDb.getVersionList().get(0).getCatalogueVersion()
							: "unknown");

			// Validate that the MAG database has the required files
			File hapFastaFile = selectedMagDb.getHapFastaFile();
			if (hapFastaFile == null || !hapFastaFile.exists()) {
				LOGGER.error("MAG database HAP FASTA file is null or does not exist: {}",
						hapFastaFile != null ? hapFastaFile.getAbsolutePath() : "null");
				LOGGER.error("Repository file: {}", selectedMagDb.getRepositoryFile());
				LOGGER.error("Current file: {}", selectedMagDb.getCurrentFile());
				throw new IllegalArgumentException("MAG database HAP FASTA file is required and must exist");
			}
			LOGGER.info("MAG database HAP FASTA file validated: {}", hapFastaFile.getAbsolutePath());
			
			File currentDbFile = selectedMagDb.getCurrentFile();
			if (currentDbFile == null || !currentDbFile.exists()) {
				LOGGER.error("MAG database file is null or does not exist: {}",
						currentDbFile != null ? currentDbFile.getAbsolutePath() : "null");
				LOGGER.error("Repository file: {}", selectedMagDb.getRepositoryFile());
				LOGGER.error("Current file: {}", selectedMagDb.getCurrentFile());
				throw new IllegalArgumentException("MAG database file is required and must exist");
			}

			LOGGER.info("MAG database file validated: {}", currentDbFile.getAbsolutePath());
			this.metaPar.setMicroDb(selectedMagDb.getCurrentFile().getAbsolutePath());
		}
		
		// Initialize Sage task
		this.sageExe = new File(mrp.getToolPath(ToolType.SAGE));
		this.sageTask = new SageTask(sageExe.getAbsolutePath());
		
		// Validate that the MAG database has the required files
		if (selectedMagDb != null) {
			File hapFastaFile = selectedMagDb.getHapFastaFile();
			if (hapFastaFile == null || !hapFastaFile.exists()) {
				LOGGER.error("MAG database HAP FASTA file is null or does not exist: {}", 
					hapFastaFile != null ? hapFastaFile.getAbsolutePath() : "null");
				LOGGER.error("Repository file: {}", selectedMagDb.getRepositoryFile());
				LOGGER.error("Current file: {}", selectedMagDb.getCurrentFile());
				throw new IllegalArgumentException("MAG database HAP FASTA file is required and must exist");
			}
			LOGGER.info("MAG database HAP FASTA file validated: {}", hapFastaFile.getAbsolutePath());
		} else {
			LOGGER.error("MAG database item is null in parameters");
			throw new IllegalArgumentException("MAG database item is required");
		}
		
		// Initialize using the DDAMagSageParameters directly
		String[] rawFiles = magSageParams.getMetadata().getRawFiles();
		this.mzMLFiles = new String[rawFiles.length];
		this.expNameMap = new HashMap<String, String>();
		this.spectraDirFile = magSageParams.getSpectraFile();
		for (int i = 0; i < rawFiles.length; i++) {
			File rawi = new File(rawFiles[i]);
			this.rawDir = rawi.getParentFile();
			String fileName = rawi.getName();
			String baseName = fileName;
			
			// Handle different file formats
			if (fileName.toLowerCase().endsWith(".mzml")) {
				// File is already in mzML format, use original path directly
				mzMLFiles[i] = rawFiles[i];
			} else if (fileName.toLowerCase().endsWith(".raw")) {
				// File is in .raw format
				baseName = fileName.substring(0, fileName.length() - ".raw".length());
				File mzMLFile = new File(spectraDirFile, baseName + ".mzML");
				mzMLFiles[i] = mzMLFile.getAbsolutePath();
			} else if (fileName.toLowerCase().endsWith(".d")) {
				// File is in .d format
				baseName = fileName.substring(0, fileName.length() - ".d".length());
				File mzMLFile = new File(spectraDirFile, baseName + ".mzML");
				mzMLFiles[i] = mzMLFile.getAbsolutePath();
			} else {
				// Unknown format, assume it needs conversion and use the original name
				File mzMLFile = new File(spectraDirFile, baseName + ".mzML");
				mzMLFiles[i] = mzMLFile.getAbsolutePath();
			}
		}

		// Initialize the result directory structure in the constructor
		File sageMagSearchFile = new File(magSageParams.getResult(), "sage_MAG_search");
		if (!sageMagSearchFile.exists()) {
			sageMagSearchFile.mkdir();
			LOGGER.info("Created Sage MAG search directory: {}", sageMagSearchFile.getAbsolutePath());
		} else {
			LOGGER.info("Sage MAG search directory already exists: {}", sageMagSearchFile.getAbsolutePath());
		}
		magSageParams.setDbSearchResultFile(sageMagSearchFile);

		this.sageExe = new File(mrp.getToolPath(ToolType.SAGE));
		this.msConvertExe = new File(mrp.getToolPath(ToolType.MSCONVERT));
		
		LOGGER.info("Sage path from MetaResourceParameters: {}", mrp.getToolPath(ToolType.SAGE));
		LOGGER.info("Sage executable file: {}", sageExe.getAbsolutePath());
		LOGGER.info("Sage exists: {}", sageExe.exists());
		
		LOGGER.info("MSConvert path from MetaResourceParameters: {}", mrp.getToolPath(ToolType.MSCONVERT));
		LOGGER.info("MSConvert executable file: {}", msConvertExe.getAbsolutePath());
		LOGGER.info("MSConvert exists: {}", msConvertExe.exists());
	}

	@Override
	public void forceStop() {
		LOGGER.info(getTaskName() + ": Force stop requested");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": Force stop requested");
		
		// Cancel Sage task if running
		if (sageTask != null) {
			LOGGER.info(getTaskName() + ": Stopping Sage identification task");
			// SageTask doesn't have forceStop method, but we can interrupt any running threads
			// The task will be stopped when this main task is cancelled
		}
		
		// Cancel FlashLFQ task if running
		if (flashLfqTask != null && !flashLfqTask.isDone()) {
			LOGGER.info(getTaskName() + ": Cancelling FlashLFQ quantification task");
			flashLfqTask.cancel();
		}
		
		// Cancel report generation task if running
		if (reportTask != null && !reportTask.isDone()) {
			LOGGER.info(getTaskName() + ": Cancelling report generation task");
			reportTask.cancel(true);
		}
		
		// Cancel executor service if running
		if (executorService != null && !executorService.isShutdown()) {
			LOGGER.info(getTaskName() + ": Shutting down executor service");
			executorService.shutdownNow();
		}
		
		// Cancel report copy thread if running
		if (reportCopyThread != null && reportCopyThread.isAlive()) {
			LOGGER.info(getTaskName() + ": Interrupting report copy thread");
			reportCopyThread.interrupt();
		}
		
		// Cancel report generation thread if running
		if (reportThread != null && reportThread.isAlive()) {
			LOGGER.info(getTaskName() + ": Interrupting report generation thread");
			reportThread.interrupt();
		}
		
		// Cancel this task itself
		this.cancel(true);
		
		LOGGER.info(getTaskName() + ": Force stop completed");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": Force stop completed");
	}

	@Override
	protected String getTaskName() {
		return taskName;
	}

	@Override
	protected Logger getLogger() {
		return LOGGER;
	}

	@Override
	protected Boolean call() throws Exception {
		updateMessage(getTaskName());
		updateProgress(0, 100);

		// Copy report templates in background
		if (!copyReportTemplatesInBackground()) {
			LOGGER.error(getTaskName() + ": Failed to copy report templates");
			return false;
		}

		LOGGER.info(getTaskName() + ": started");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": started");
		
		// Log detailed parameter information
		LOGGER.info("Workflow parameters:");
		LOGGER.info("  Result directory: " + magSageParams.getResult());
		LOGGER.info("  Quantitative mode: " + magSageParams.getQuantitativeMode());
		LOGGER.info("  Thread count: " + magSageParams.getThreadCount());
		LOGGER.info("  MS2 scan mode: " + magSageParams.getMs2ScanMode());
		if (magSageParams.getMetadata() != null) {
			String[] rawFiles = magSageParams.getMetadata().getRawFiles();
			LOGGER.info("  Raw files count: " + (rawFiles != null ? rawFiles.length : 0));
			if (rawFiles != null) {
				for (int i = 0; i < rawFiles.length; i++) {
					LOGGER.info("    Raw file " + (i+1) + ": " + rawFiles[i]);
				}
			}
		}
		
		// Log resource information
		LOGGER.info("Resource parameters:");
		LOGGER.info("  Sage path: " + mrp.getToolPath(MetaResourceParameters.ToolType.SAGE));
		
		// Update main progress and keep second bar active
		updateProgress(5, 100);
		if (statusLabel != null) {
			Platform.runLater(() -> statusLabel.setText(getTaskName() + ": Initializing..."));
		}
		Platform.runLater(() -> secondBar.setProgress(ProgressBar.INDETERMINATE_PROGRESS));

		// Step 1: Spectra Conversion
		LOGGER.info("Step 1: Starting spectra format conversion to mzML...");
		if (statusLabel != null) {
			Platform.runLater(() -> statusLabel.setText(getTaskName() + ": Converting spectra format..."));
		}
		updateProgress(10, 100);
		Platform.runLater(() -> secondBar.setProgress(ProgressBar.INDETERMINATE_PROGRESS));
		
		boolean convert = convertSpectraToMzML();
		if (!convert) {
			LOGGER.error(getTaskName() + ": error in spectra format conversion");
			System.err.println(format.format(new Date()) + "\t" + getTaskName() + ": error in spectra format conversion");
			updateMessage(getTaskName() + ": Failed - Spectra conversion error");
			Platform.runLater(() -> secondBar.setProgress(0.0));
			return false;
		}
		LOGGER.info("Step 1: Spectra format conversion completed successfully");

		// Step 2: High Abundant Database Search
		LOGGER.info("Step 2: Starting high abundant database search...");
		updateProgress(10, 100);
		updateMessage(getTaskName() + ": Searching high abundant database...");
		Platform.runLater(() -> secondBar.setProgress(ProgressBar.INDETERMINATE_PROGRESS));
		
		System.out.println(
				format.format(new Date()) + "\t" + getTaskName() + ": searching high abundant database started");

		// Initialize separateResultFolders
		LOGGER.info("Initializing result directories in: {}", magSageParams.getDbSearchResultFile().getAbsolutePath());
		this.separateResultFolders = new File[this.mzMLFiles.length];
		for (int i = 0; i < this.mzMLFiles.length; i++) {
			String rawName = new File(this.mzMLFiles[i]).getName();
			rawName = rawName.substring(0, rawName.length() - ".mzML".length());
			this.separateResultFolders[i] = new File(magSageParams.getDbSearchResultFile(), rawName);
			if (!this.separateResultFolders[i].exists()) {
				this.separateResultFolders[i].mkdir();
				LOGGER.info("  Created result directory: {}", this.separateResultFolders[i].getAbsolutePath());
			} else {
				LOGGER.info("  Result directory already exists: {}", this.separateResultFolders[i].getAbsolutePath());
			}
		}
		
		boolean searchHap = this.separateSearchHap();
		if (!searchHap) {
			LOGGER.error(getTaskName() + ": error in searching high abundant database");
			return false;
		}

		// Step 3: Protein Extraction from HAP Results
		LOGGER.info("Step 3: Starting protein extraction from HAP results...");
		updateProgress(40, 100);
		boolean extract = this.fastaExtract4Iterator();
		if (!extract) {
			LOGGER.error(getTaskName() + ": error in protein extraction");
			return false;
		}

		// Step 4: Database Search
		LOGGER.info("Step 4: Starting database search through sample-specific database...");
		updateProgress(65, 100);
		boolean search = this.separateSearch();
		if (!search) {
			LOGGER.error(getTaskName() + ": error in database search");
			return false;
		}

		// Step 5: Combine Results
		LOGGER.info("Step 5: Starting combination of separate search results...");
		updateProgress(70, 100);
		boolean combine = this.combineResults();
		if (!combine) {
			LOGGER.error(getTaskName() + ": error in combining separate search results");
			return false;
		}

		// Step 6: Quantification
		LOGGER.info("Step 6: Starting quantification...");
		updateProgress(75, 100);
		boolean quant = false;

		if (magSageParams.getQuantitativeMode().equals(MetaConstants.isobaricLabel)) {
			LOGGER.info("  Using isobaric quantification mode");
			updateProgress(80, 100);
			quant = lfQuant();
			if (!quant) {
				LOGGER.error(getTaskName() + ": error in isobaric quantification");
				return false;
			}
		}

		// Step 7: Export Report
		LOGGER.info("Step 7: Starting report generation...");
		updateProgress(85, 100);
		boolean report = this.exportReport();
		if (!report) {
			LOGGER.error(getTaskName() + ": error in report generation");
			return false;
		}

		// Update progress for completion
		updateProgress(100, 100);
		updateMessage(getTaskName() + ": Completed successfully!");
		Platform.runLater(() -> secondBar.setProgress(0.0));
		
		LOGGER.info(getTaskName() + ": finished - all workflow steps completed successfully");
		return true;
	}

	@Override
	protected String getVersion() {
		return null;
	}
	
	/**
	 * Convert spectra files to mzML format for Sage processing using multiple threads
	 */
	private boolean convertSpectraToMzML() {
		LOGGER.info("Starting spectra conversion to mzML format using MSConvert");
		
		// Validate MSConvert executable
		if (msConvertExe == null || !msConvertExe.exists()) {
			LOGGER.error("MSConvert executable not found: {}", 
				msConvertExe != null ? msConvertExe.getAbsolutePath() : "null");
			return false;
		}
		
		String[] rawFiles = magSageParams.getMetadata().getRawFiles();
		if (rawFiles == null || rawFiles.length == 0) {
			LOGGER.error("No raw files found in metadata");
			return false;
		}
		
		// Ensure spectra directory exists
		File spectraDir = magSageParams.getSpectraFile();
		if (!spectraDir.exists()) {
			spectraDir.mkdirs();
			LOGGER.info("Created spectra directory: {}", spectraDir.getAbsolutePath());
		}
		
		// mzMLFiles array is already initialized in constructor
		
		// Collect files that need conversion
		List<Integer> filesToConvert = new ArrayList<>();
		for (int i = 0; i < rawFiles.length; i++) {
			String rawFile = rawFiles[i];
			File rawFileObj = new File(rawFile);
			
			if (!rawFileObj.exists()) {
				LOGGER.error("Raw file does not exist: {}", rawFile);
				return false;
			}
			
			// Check if the mzML file already exists (mzMLFiles[i] was set in constructor)
			File mzMLFile = new File(mzMLFiles[i]);
			
			if (mzMLFile.exists() && mzMLFile.length() > 0) {
				LOGGER.info("mzML file already exists, skipping conversion: {}", mzMLFile.getAbsolutePath());
			} else {
				filesToConvert.add(i);
			}
		}
		
		// If no files need conversion, return success
		if (filesToConvert.isEmpty()) {
			LOGGER.info("All mzML files already exist, no conversion needed");
			return true;
		}
		
		// Use multiple threads for conversion
		int threadCount = Math.min(filesToConvert.size(), magSageParams.getThreadCount());
		ExecutorService executor = Executors.newFixedThreadPool(threadCount);
		
		LOGGER.info("Processing {} files using {} threads", filesToConvert.size(), threadCount);
		
		// Submit conversion/copy tasks
		for (int fileIndex : filesToConvert) {
			executor.submit(() -> {
				String rawFile = rawFiles[fileIndex];
				File rawFileObj = new File(rawFile);
				String fileName = rawFileObj.getName().toLowerCase();
				
				try {
					String convertedPath = null;
					
					if (fileName.endsWith(".raw")) {
						LOGGER.info("Converting .raw file to mzML: {}", rawFile);
						convertedPath = MSConvertor.raw2mzML(msConvertExe.getAbsolutePath(), rawFile);
					} else if (fileName.endsWith(".d")) {
						LOGGER.info("Converting .d file to mzML: {}", rawFile);
						convertedPath = MSConvertor.d2mzML(msConvertExe.getAbsolutePath(), rawFile, spectraDir.getAbsolutePath());
					} else {
						LOGGER.error("Unsupported file format: {}. Only .raw and .d files are supported for conversion.", fileName);
						return;
					}
					
					if (convertedPath != null) {
						// Verify the converted file exists
						File convertedFile = new File(convertedPath);
						if (convertedFile.exists() && convertedFile.length() > 0) {
							// Move the converted file to the spectra directory
							File targetFile = new File(mzMLFiles[fileIndex]);
							try {
								Files.move(convertedFile.toPath(), targetFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
								LOGGER.info("Successfully moved converted file from {} to {}", convertedPath, mzMLFiles[fileIndex]);
							} catch (IOException e) {
								LOGGER.error("Failed to move converted file from {} to {}: {}", convertedPath, mzMLFiles[fileIndex], e.getMessage());
							}
						} else {
							LOGGER.error("Converted mzML file does not exist or is empty: {}", convertedPath);
						}
					} else {
						LOGGER.error("Failed to convert file to mzML: {}", rawFile);
					}
				} catch (Exception e) {
					LOGGER.error("Error processing file {}: {}", rawFile, e.getMessage());
				}
			});
		}
		
		// Shutdown executor and wait for completion
		executor.shutdown();
		
		try {
			// Wait for all conversions to complete (with timeout)
			boolean finished = executor.awaitTermination(filesToConvert.size() * 4, TimeUnit.HOURS);
			if (finished) {
				LOGGER.info("All file conversions completed successfully");
			} else {
				LOGGER.error("File conversion timeout - some files may not have been converted");
				executor.shutdownNow();
				return false;
			}
		} catch (InterruptedException e) {
			LOGGER.error("File conversion was interrupted: {}", e.getMessage());
			executor.shutdownNow();
			return false;
		}
		
		// Verify all files were converted successfully
		for (int i = 0; i < mzMLFiles.length; i++) {
			if (mzMLFiles[i] == null) {
				LOGGER.error("Failed to convert file: {}", rawFiles[i]);
				return false;
			}
		}
		
		LOGGER.info("Spectra conversion completed successfully. Converted {} files to mzML format.", rawFiles.length);
		return true;
	}

	/**
	 * Perform high abundant database search using Sage
	 */
	private boolean separateSearchHap() {
		LOGGER.info("Starting separate search for high abundant database using Sage");
		
		// Get thread configuration
		int threadCount = magSageParams.getThreadCount();
		LOGGER.info("Thread configuration: {}", threadCount);
		
		File hapFastaFile = this.magDb.getHapFastaFile();
		if (hapFastaFile == null || !hapFastaFile.exists()) {
			LOGGER.error("HAP FASTA file is null or does not exist for MAG database: {}",
					this.magDb.getCatalogueID());
			return false;
		}

		// Check if we need to perform new searches
		boolean needNewSearch = false;
		int filesNeedingSearch = 0;
		for (int i = 0; i < this.mzMLFiles.length; i++) {
			File hapFile = new File(separateResultFolders[i], "hap");
			File resultFile = new File(hapFile, "results.sage.tsv");
			if (!resultFile.exists() || resultFile.length() == 0) {
				needNewSearch = true;
				filesNeedingSearch++;
			}
		}
		
		if (needNewSearch) {
			LOGGER.info("Found {} files requiring new search out of {} total files", 
				filesNeedingSearch, this.mzMLFiles.length);
			
			for (int i = 0; i < this.mzMLFiles.length; i++) {
				File hapFile = new File(separateResultFolders[i], "hap");
				if (!hapFile.exists()) {
					hapFile.mkdirs();
				}
				File mzMLFile = new File(mzMLFiles[i]);
				if (mzMLFile.exists()) {
					LOGGER.info("  Searching file {} with Sage", i+1);
					
					// Use the enhanced searchMzML method with separate min/max precursor tolerance
					sageTask.searchMzML(new String[] { mzMLFile.getAbsolutePath() }, 
							hapFile.getAbsolutePath(),
							hapFastaFile.getAbsolutePath(),
							magSageParams.getPrecursorToleranceMin(),
							magSageParams.getPrecursorToleranceMax(),
							magSageParams.getFragmentTolerance(),
							magSageParams.getEnzyme(),
							magSageParams.getMissCleavages(),
							magSageParams.getFixMods(),
							magSageParams.getVariMods(),
							false); // isQuant = false for HAP search
				} else {
					LOGGER.warn("mzML file {} was not found, skipping", mzMLFiles[i]);
				}
		}

		// Set up WatchService to monitor results.sage.tsv file generation
		FileFilter filter = new FileFilter() {
			@Override
			public boolean accept(File file) {
				return file.getName().equals("results.sage.tsv");
			}
		};

		// Use the generic WatchServiceLoop from the superclass
		LOGGER.info("Starting WatchService monitoring for high abundant database search");
		startWatchServiceLoop(separateResultFolders, "hap", filter, filesNeedingSearch, 10, 40);
		
		LOGGER.info("Starting Sage execution");
		sageTask.run(magSageParams.getThreadCount(), separateResultFolders.length * 3);
		} else {
			LOGGER.info("No new searches required, all results already exist");
		}

		// Verify search results
		LOGGER.info("Verifying search results for {} folders", separateResultFolders.length);
		for (int i = 0; i < this.separateResultFolders.length; i++) {
			File hapFile = new File(separateResultFolders[i], "hap");
//			LOGGER.info("  Cleaning and verifying results for folder {}: {}", i+1, hapFile.getAbsolutePath());
//			SageTask.cleanSageResultFolder(hapFile, true);
			
			File resultFile = new File(hapFile, "results.sage.tsv");
			if (!resultFile.exists() || resultFile.length() == 0) {
				LOGGER.error("Sage search result was not found from {}", hapFile.getAbsolutePath());
				return false;
			}
			LOGGER.info("  Results verified successfully for folder {}", i+1);
		}

		LOGGER.info("High abundant database search completed successfully for all {} files", separateResultFolders.length);
		
		// Stop the watch service thread
		stopWatchServiceThread();
		
		return true;
	}

	/**
	 * Create Sage configuration file for database search
	 */
	private File createSageConfigFile(File mzMLFile, File fastaFile, File outputDir, boolean isQuant) {
		File configFile = new File(outputDir, "sage_config.json");
		
		try (PrintWriter writer = new PrintWriter(configFile)) {
			writer.println("{");
			writer.println("  \"database\": {");
			writer.println("    \"fasta\": \"" + fastaFile.getAbsolutePath().replace("\\", "\\\\") + "\",");
			writer.println("    \"decoy_tag\": \"" + magSageParams.getDecoysTag() + "\",");
			writer.println("    \"generate_decoys\": " + magSageParams.isGenerateDecoys());
			writer.println("  },");
			writer.println("  \"precursor_tol\": {");
			writer.println("    \"ppm\": [" + magSageParams.getPrecursorToleranceMin() + ", " + magSageParams.getPrecursorToleranceMax() + "]");
			writer.println("  },");
			writer.println("  \"fragment_tol\": {");
			writer.println("    \"ppm\": [0, " + magSageParams.getFragmentTolerance() + "]");
			writer.println("  },");
			writer.println("  \"enzyme\": {");
			writer.println("    \"missed_cleavages\": " + magSageParams.getMissCleavages());
			writer.println("  },");
			writer.println("  \"output_directory\": \"" + outputDir.getAbsolutePath().replace("\\", "\\\\") + "\",");
			writer.println("  \"mzml_paths\": [\"" + mzMLFile.getAbsolutePath().replace("\\", "\\\\") + "\"]");
			
			// Add quantification settings if isQuant is true
			if (isQuant) {
				writer.println(",");
				writer.println("  \"quant\": {");
				writer.println("    \"lfq\": true,");
				writer.println("    \"lfq_settings\": {");
				writer.println("      \"peak_scoring\": \"" + magSageParams.getPeakScoring() + "\",");
				writer.println("      \"integration\": \"" + magSageParams.getIntegration() + "\",");
				writer.println("      \"spectral_angle\": " + magSageParams.getSpectralAngle() + ",");
				writer.println("      \"ppm_tolerance\": " + magSageParams.getLfqPpmTolerance());
				writer.println("    }");
				writer.println("  }");
			}
			
			writer.println("}");
			
			LOGGER.info("Created Sage config file: {}", configFile.getAbsolutePath());
			return configFile;
			
		} catch (IOException e) {
			LOGGER.error("Failed to create Sage config file: {}", configFile.getAbsolutePath(), e);
			return null;
		}
	}

	/**
	 * Extract protein sequences from HAP search results
	 */
	protected boolean fastaExtract4Iterator() {
		LOGGER.info("Starting FASTA extraction for {} result folders", separateResultFolders.length);

		for (int i = 0; i < this.separateResultFolders.length; i++) {
			File hapFile = new File(separateResultFolders[i], "hap");
			File dbFile = new File(hapFile, separateResultFolders[i].getName() + ".fasta");

			if (dbFile.exists() && dbFile.length() > 0) {
				LOGGER.info("Sample-specific database for {} already exists", separateResultFolders[i].getName());
				continue;
			}

			File spectraFile = new File(hapFile, "results.sage.tsv");
			HashMap<String, HashSet<String>> genomeMap = new HashMap<String, HashSet<String>>();
			HashSet<String> pepSet = new HashSet<String>();

			try (BufferedReader hapReader = new BufferedReader(new FileReader(spectraFile))) {
				String headerLine = hapReader.readLine(); // Read header
				String[] header = headerLine.split("\t");

				// Find the index of "proteins" column
				int proteinsIndex = -1;
				for (int j = 0; j < header.length; j++) {
					if (header[j].toLowerCase().contains("proteins")) {
						proteinsIndex = j;
						break;
					}
				}

				if (proteinsIndex == -1) {
					LOGGER.error("Could not find 'proteins' column in header: {}", headerLine);
					continue;
				}

				LOGGER.info("Found proteins column at index {}: {}", proteinsIndex, header[proteinsIndex]);

				String line;
				while ((line = hapReader.readLine()) != null) {
					String[] cs = line.split("\t");
					if (cs.length > proteinsIndex) {
						String[] pros = cs[proteinsIndex].split(";");
						for (String pro : pros) {
							if (pro.startsWith(magDb.getIdentifier())) {
								String genome = pro.split("_")[0];
								genomeMap.computeIfAbsent(genome, k -> new HashSet<>()).add(cs[0]);
							}
						}
						pepSet.add(cs[0]);
					}
				}
			} catch (IOException e) {
				LOGGER.error("Error reading Sage search result file from {}", spectraFile.getAbsolutePath(), e);
				return false;
			}

			LOGGER.info("Sample-specific database for {} - Found {} peptides and {} genomes",
					separateResultFolders[i].getName(), pepSet.size(), genomeMap.size());
			// Create sample-specific database
			File hapGenomeFile = new File(hapFile, "hap_genomes.tsv");
			HashSet<String> genomeSet = refineGenomes(genomeMap, pepSet.size(), psmCoverThreshold, hapGenomeFile);

			try (PrintWriter writer = new PrintWriter(dbFile)) {
				File dbFolder = this.magDb.getOriginalDbFolder();
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
						}
					}
				}
			} catch (IOException e) {
				LOGGER.error("Error writing sample-specific database to {}", dbFile.getAbsolutePath(), e);
				return false;
			}
		}

		return true;
	}

	/**
	 * Perform sample-specific database search using Sage
	 */
	protected boolean separateSearch() {
		int taskCount = 0;
		for (int i = 0; i < this.separateResultFolders.length; i++) {
			File resultFile = new File(separateResultFolders[i], "results.sage.tsv");
			if (resultFile.exists() && resultFile.length() > 0) {
				System.out.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": Sage search result already exists in " + separateResultFolders[i].getAbsolutePath());

				LOGGER.info(getTaskName() + ": Sage search result already exists in "
						+ separateResultFolders[i].getAbsolutePath());
			} else {
				File hapFile = new File(separateResultFolders[i], "hap");
				File db = new File(hapFile, separateResultFolders[i].getName() + ".fasta");
				if (!db.exists()) {
					System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": sample-specific database"
							+ db + " was not found");

					LOGGER.error(getTaskName() + ": sample-specific database" + db + " was not found");
					return false;
				}
				taskCount++;
			}
		}

		if (taskCount > 0) {
			for (int i = 0; i < this.mzMLFiles.length; i++) {
				File mzMLFile = new File(mzMLFiles[i]);
				File hapFile = new File(separateResultFolders[i], "hap");
				File db = new File(hapFile, separateResultFolders[i].getName() + ".fasta");

				if (mzMLFile.exists()) {
					if (db.exists() && db.length() > 0) {
						// Use the enhanced searchMzML method for sample-specific database search
						sageTask.searchMzML(new String[] { mzMLFile.getAbsolutePath() },
								separateResultFolders[i].getAbsolutePath(), db.getAbsolutePath(),
								magSageParams.getPrecursorToleranceMin(), magSageParams.getPrecursorToleranceMax(),
								magSageParams.getFragmentTolerance(), magSageParams.getEnzyme(),
								magSageParams.getMissCleavages(), magSageParams.getFixMods(),
								magSageParams.getVariMods(), false);
					} else {
						System.out.println(
								format.format(new Date()) + "\t" + getTaskName() + ": no protein sequences were in "
										+ db + ", which means no peptides/proteins can be identifed from "
										+ separateResultFolders[i].getName());

						LOGGER.error(getTaskName() + ": no protein sequences were in " + db
								+ ", which means no peptides/proteins can be identifed from "
								+ separateResultFolders[i].getName());
					}
				}
		}

		// Set up WatchService to monitor results.sage.tsv file generation
		FileFilter filter = new FileFilter() {
			@Override
			public boolean accept(File file) {
				return file.getName().equals("results.sage.tsv");
			}
		};

		// Use the generic WatchServiceLoop from the superclass
		LOGGER.info("Starting WatchService monitoring for sample-specific database search");
		startWatchServiceLoop(separateResultFolders, "", filter, taskCount, 40, 65);

		sageTask.run(magSageParams.getThreadCount(), this.separateResultFolders.length * 10);
		}

		// Clean and verify results
		for (int i = 0; i < this.separateResultFolders.length; i++) {
//			SageTask.cleanSageResultFolder(separateResultFolders[i], true);
			
			File hapFile = new File(separateResultFolders[i], "hap");
			if (hapFile.exists()) {
//				SageTask.cleanSageResultFolder(hapFile, true);
			}
			
			File resultFile = new File(separateResultFolders[i], "results.sage.tsv");
			if (!resultFile.exists() || resultFile.length() == 0) {
				LOGGER.error("Sage search result was not found in {}", separateResultFolders[i].getAbsolutePath());
				return false;
			}
		}

		// Stop the watch service thread
		stopWatchServiceThread();

		return true;
	}

	/**
	 * Combine search results from all samples
	 */
	protected boolean combineResults() {
		HashSet<String> totalPsmSet = new HashSet<String>();
		HashMap<String, HashSet<String>> genomeMap = new HashMap<String, HashSet<String>>();
		
		// Parse individual Sage result files from each sample
		for (int i = 0; i < separateResultFolders.length; i++) {
			File pepResultFile = new File(separateResultFolders[i], "results.sage.tsv");
			try {
				BufferedReader combineTsvReader = new BufferedReader(new FileReader(pepResultFile));
				String line = combineTsvReader.readLine();
				String[] title = line.split("\t");
				int proid = -1;
				int pepid = -1;
				int qvalueid = -1;
				for (int j = 0; j < title.length; j++) {
					if (title[j].equals("proteins")) {
						proid = j;
					} else if (title[j].equals("peptide")) {
						pepid = j;
					} else if (title[j].equals("spectrum_q")) {
						qvalueid = j;
					}
				}
				while ((line = combineTsvReader.readLine()) != null) {
					String[] cs = line.split("\t");
					if (Double.parseDouble(cs[qvalueid]) > 0.01 || cs[proid].contains("REV_")) {
						continue;
					}
					String[] pros = cs[proid].split(";");

					for (int j = 0; j < pros.length; j++) {
						if (genomeMap.containsKey(pros[j])) {
							genomeMap.get(pros[j]).add(cs[pepid]);
						} else {
							HashSet<String> set = new HashSet<String>();
							set.add(cs[pepid]);
							genomeMap.put(pros[j], set);
						}
						totalPsmSet.add(cs[pepid]);
					}
				}
				combineTsvReader.close();
			} catch (IOException e) {
				LOGGER.error(getTaskName() + ": error in reading database search result file " + pepResultFile, e);
				System.err.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": error in reading database search result file " + pepResultFile);
			}
		}

		File combinedFile = new File(magSageParams.getDbSearchResultFile(), "combined");
		if(!combinedFile.exists()) {
			combinedFile.mkdirs();
		}
		// Create combined database
		File combinedFastaFile = new File(combinedFile, "combined.fasta");
		if (!combinedFastaFile.exists() || combinedFastaFile.length() == 0) {
			File hapGenomeFile = new File(combinedFile, "razorProteins.tsv");
			HashSet<String> proteinSet = refineGenomes(genomeMap, totalPsmSet.size(), 1.0, hapGenomeFile);

			try (PrintWriter writer = new PrintWriter(combinedFastaFile)) {
				Random random = new Random();
				for (int i = 0; i < separateResultFolders.length; i++) {
					File hapFile = new File(separateResultFolders[i], "hap");
					File dbFile = new File(hapFile, separateResultFolders[i].getName() + ".fasta");
					try (BufferedReader reader = new BufferedReader(new FileReader(dbFile))) {
						String line = null;
						boolean use = false;
						String shufName = "";
						List<Character> sequenceChars = new ArrayList<Character>();
						while ((line = reader.readLine()) != null) {
							if (line.startsWith(">")) {
								if (sequenceChars.size() > 0) {
									Collections.shuffle(sequenceChars, random);

									StringBuilder shuffledSequence = new StringBuilder(sequenceChars.size());
									for (Character c : sequenceChars) {
										shuffledSequence.append(c);
									}
									writer.println(shufName);
									writer.println(shuffledSequence);

									sequenceChars = new ArrayList<Character>();
								}

								String proName = line.substring(1, line.indexOf(" "));
								if (proteinSet.contains(proName)) {
									use = true;
									writer.println(line);
									proteinSet.remove(proName);

									shufName = ">" + FastaManager.shuffle + line.substring(1);
								} else {
									use = false;
								}
							} else {
								if (use) {
									writer.println(line);
									for (int j = 0; j < line.length(); j++) {
										sequenceChars.add(line.charAt(j));
									}
								}
							}
						}
						reader.close();

						if (sequenceChars.size() > 0) {
							Collections.shuffle(sequenceChars, random);

							StringBuilder shuffledSequence = new StringBuilder(sequenceChars.size());
							for (Character c : sequenceChars) {
								shuffledSequence.append(c);
							}
							writer.println(shufName);
							writer.println(shuffledSequence);
						}

					} catch (IOException e) {
						LOGGER.error(getTaskName() + ": error in reading database file " + dbFile, e);
						System.err.println(format.format(new Date()) + "\t" + getTaskName()
								+ ": error in reading database file " + dbFile);
					}
				}
				writer.close();
			} catch (IOException e) {
				LOGGER.error(getTaskName() + ": error in writing the combined database file to " + combinedFastaFile, e);
				System.err.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": error in writing the combined database file to " + combinedFastaFile);
			}
		}

		// Run final Sage search against combined database
		File pepResultFile = new File(combinedFile, "results.sage.tsv");
		File quanResultFile = new File(combinedFile, "lfq.tsv");

		if (!pepResultFile.exists() && !quanResultFile.exists()) {
			magSageParams.setCurrentDb(combinedFastaFile.getAbsolutePath());

			// Run Sage search against combined database
			try {
				sageTask.searchMzML(mzMLFiles, combinedFile.getAbsolutePath(), combinedFastaFile.getAbsolutePath(),
						magSageParams.getPrecursorToleranceMin(), magSageParams.getPrecursorToleranceMax(),
						magSageParams.getFragmentTolerance(), magSageParams.getEnzyme(), magSageParams.getMissCleavages(), 
						magSageParams.getFixMods(), magSageParams.getVariMods(), true);
			} catch (Exception e) {
				LOGGER.error(getTaskName() + ": error in running Sage search against combined database", e);
				System.out.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": error in running Sage search against combined database");
				return false;
			}

			LOGGER.info(getTaskName() + ": searching combined database started");
			System.out.println(
					format.format(new Date()) + "\t" + getTaskName() + ": searching combined database started");
			
			LOGGER.info("Starting Sage execution");
			sageTask.run(separateResultFolders.length * 3);
			
		} else {
			System.out.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": Sage search result already exists in " + magSageParams.getDbSearchResultFile());
			LOGGER.info(getTaskName() + ": Sage search result already exists in " + magSageParams.getDbSearchResultFile());
		}

		// Check if results exist
		if (pepResultFile.exists() && pepResultFile.length() > 0 && quanResultFile.exists()
				&& quanResultFile.length() > 0) {
			LOGGER.info(getTaskName() + ": combined database search finished successfully");
			System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": combined database search finished successfully");
			return true;
		} else {
			LOGGER.info(getTaskName() + ": combined database search failed");
			System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": combined database search failed");
			return false;
		}
	}

	/**
	 * Perform label-free quantification
	 */
	protected boolean lfQuant() {
		String flashLFQ = mrp.getToolPath(ToolType.FLASHLFQ);
		flashLfqTask = new FlashLfqTask(flashLFQ, magSageParams.getThreadCount());
		this.quan_pep_file = new File(magSageParams.getDbSearchResultFile(), "QuantifiedPeptides.tsv");
		
		if (!quan_pep_file.exists() || quan_pep_file.length() == 0) {
			File msms = new File(magSageParams.getDbSearchResultFile(), "psms.tsv");
			boolean lfqQuan = flashLfqTask.run(magSageParams.getMetadata(), msms.getAbsolutePath(),
					rawDir.getAbsolutePath(), true, true);
			return lfqQuan && quan_pep_file.exists();
		}
		return true;
	}

	/**
	 * Export final reports
	 */
	@SuppressWarnings("unchecked")
	protected boolean exportReport() {
		LOGGER.info(getTaskName() + ": exporting peptide report started");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": exporting peptide report started");

		ArrayList<String> fileList = new ArrayList<String>();
		HashMap<String, Integer> fileIdMap = new HashMap<String, Integer>();
		File expDesignFile = new File(rawDir, "ExperimentalDesign.tsv");

		// Build file mapping exactly like pFind - try experimental design first, then result folders
		if (expDesignFile.exists()) {
			try (BufferedReader reader = new BufferedReader(new FileReader(expDesignFile))) {
				String line = reader.readLine();
				while ((line = reader.readLine()) != null) {
					String[] cs = line.split("\t");
					fileIdMap.put(cs[0], fileIdMap.size());
					fileList.add(cs[0]);
				}
			} catch (IOException e) {
				LOGGER.error(getTaskName() + ": error in reading the experimental design from " + expDesignFile, e);
				System.err.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": error in reading the experimental design from " + expDesignFile);
			}
		} else {
			File[] files = magSageParams.getDbSearchResultFile().listFiles();
			for (int i = 0; i < files.length; i++) {
				if (files[i].isDirectory() && (new File(files[i], "results.sage.tsv")).exists()) {
					fileIdMap.put(files[i].getName(), fileIdMap.size());
					fileList.add(files[i].getName());
				}
			}
		}

		File combinedFile = new File(magSageParams.getDbSearchResultFile(), "combined");
		// PSM-level data processing - exactly like pFind but using Sage PSM format
		HashMap<String, SagePeptide> sagePsmMap = new HashMap<String, SagePeptide>();
		HashMap<String, int[]> pepPsmCountMap = new HashMap<String, int[]>();
		HashMap<String, int[]> pepChargeMap = new HashMap<String, int[]>();
		HashMap<String, double[]> pepExpectMap = new HashMap<String, double[]>();
		HashMap<String, Integer> psmCountMap = new HashMap<String, Integer>();

		int[] filePsmCount = new int[fileIdMap.size()];
		HashSet<String>[] fileSequenceSets = new HashSet[fileIdMap.size()];
		for (int i = 0; i < fileSequenceSets.length; i++) {
			fileSequenceSets[i] = new HashSet<String>();
		}

		// Process PSMs from Sage results (equivalent to FlashLfqPsmReader for pFind)
		File psmsFile = new File(combinedFile, "results.sage.tsv");
		if (psmsFile.exists()) {
			// Read Sage PSM data - this would need proper Sage PSM reader implementation
			// For now, we'll use SageResultReader but process it like PSMs
			SageResultReader sageReader = new SageResultReader(combinedFile);
			SagePeptide[] sagePsms = sageReader.getMetaPeptides();
			
			for (SagePeptide psm : sagePsms) {
				String sequence = psm.getSequence();
				String modseq = psm.getModSeq();
				// For Sage, we need to extract file information from the peptide data
				// This is format-specific and may need adjustment based on actual Sage output
				
				int[] psmCount;
				if (pepPsmCountMap.containsKey(modseq)) {
					psmCount = pepPsmCountMap.get(modseq);
				} else {
					psmCount = new int[fileIdMap.size()];
					pepPsmCountMap.put(modseq, psmCount);
				}
				
				// Update PSM counts and file mappings
				int[] ms2Counts = psm.getMs2Counts();
				for (int j = 0; j < ms2Counts.length && j < fileIdMap.size(); j++) {
					psmCount[j] += ms2Counts[j];
					filePsmCount[j] += ms2Counts[j];
					if (ms2Counts[j] > 0) {
						fileSequenceSets[j].add(sequence);
					}
				}

				int[] charges = psm.getCharges();
				if (pepChargeMap.containsKey(modseq)) {
					int[] existingCharges = pepChargeMap.get(modseq);
					for (int j = 0; j < charges.length && j < existingCharges.length; j++) {
						if (charges[j] == 1) existingCharges[j] = 1;
					}
				} else {
					pepChargeMap.put(modseq, charges.clone());
				}
				
				double[] expect;
				if (pepExpectMap.containsKey(modseq)) {
					expect = pepExpectMap.get(modseq);
				} else {
					expect = new double[fileIdMap.size()];
					Arrays.fill(expect, 1.0);
					pepExpectMap.put(modseq, expect);
				}
				
				double pep = psm.getPEP();
				for (int j = 0; j < expect.length; j++) {
					if (pep < expect[j]) {
						expect[j] = pep;
					}
				}

				// Keep best PSM for each peptide (lowest PEP)
				if (sagePsmMap.containsKey(modseq)) {
					if (pep < sagePsmMap.get(modseq).getPEP()) {
						sagePsmMap.put(modseq, psm);
					}
				} else {
					sagePsmMap.put(modseq, psm);
				}

				// Count PSMs per protein and genome
				String[] pros = psm.getProteins();
				for (String pro : pros) {
					if (psmCountMap.containsKey(pro)) {
						psmCountMap.put(pro, psmCountMap.get(pro) + 1);
					} else {
						psmCountMap.put(pro, 1);
					}

					if (pro.startsWith(this.magDb.getIdentifier())) {
						String genome = pro.split("_")[0];
						if (psmCountMap.containsKey(genome)) {
							psmCountMap.put(genome, psmCountMap.get(genome) + 1);
						} else {
							psmCountMap.put(genome, 1);
						}
					}
				}
			}
		}

		// Quantification data processing - use quan_pep_file like pFind
		SageResultReader quanPepReader;
		String quanMode = magSageParams.getQuantitativeMode();
		if (quanMode.equals(MetaConstants.labelFree)) {
			// For Sage, we use the quantified peptides from FlashLFQ or Sage's own quantification
			quanPepReader = new SageResultReader(combinedFile);
		} else {
			LOGGER.error(getTaskName() + ": Isobaric quantification not yet supported for Sage workflow");
			System.err.println(format.format(new Date()) + "\t" + getTaskName() + ": Isobaric quantification not yet supported");
			return false;
		}

		SagePeptide[] quanPeps = quanPepReader.getMetaPeptides();
		
		// Set individual PEP values from PSM data
		for (int i = 0; i < quanPeps.length; i++) {
			String seq = quanPeps[i].getModSeq();
			if (pepExpectMap.containsKey(seq)) {
				// For Sage, we need to set this information if the API supports it
				// This may require extending SagePeptide class
			}
		}

		// Use the exact same protein inference strategy as pFind
		HashMap<String, Double> genomeScoreMap = this.refineGenomes(quanPeps);
		String[] razorProteins = this.refineProteins(quanPeps, genomeScoreMap);

		updateProgress(82, 100);

		HashMap<String, double[]> genomeIntensityMap = new HashMap<String, double[]>();
		HashMap<String, double[]> proteinIntensityMap = new HashMap<String, double[]>();
		HashMap<String, Double> proteinTotalIntenMap = new HashMap<String, Double>();

		HashSet<String> pepSequenceSet = new HashSet<String>();
		HashMap<String, HashSet<String>> totalPepCountMap = new HashMap<String, HashSet<String>>();
		HashMap<String, HashSet<String>> razorPepCountMap = new HashMap<String, HashSet<String>>();

		ArrayList<SagePeptide> peptideList = new ArrayList<SagePeptide>();

		// First pass: build razor peptide maps (exactly like pFind)
		for (int i = 0; i < quanPeps.length; i++) {
			if (razorProteins[i] != null) {
				String stripeSequence = quanPeps[i].getSequence();

				pepSequenceSet.add(stripeSequence);
				peptideList.add(quanPeps[i]);

				if (razorPepCountMap.containsKey(razorProteins[i])) {
					razorPepCountMap.get(razorProteins[i]).add(stripeSequence);
				} else {
					HashSet<String> set = new HashSet<String>();
					set.add(stripeSequence);
					razorPepCountMap.put(razorProteins[i], set);
				}

				if (razorProteins[i].startsWith(magDb.getIdentifier())) {
					String genome = razorProteins[i].split("_")[0];
					if (razorPepCountMap.containsKey(genome)) {
						razorPepCountMap.get(genome).add(stripeSequence);
					} else {
						HashSet<String> set = new HashSet<String>();
						set.add(stripeSequence);
						razorPepCountMap.put(genome, set);
					}
				}
			}
		}

		// Second pass: intensity calculation and protein mapping (exactly like pFind)
		for (int i = 0; i < quanPeps.length; i++) {
			String[] proteins = quanPeps[i].getProteins();
			int count = 0;
			for (int j = 0; j < proteins.length; j++) {
				if (razorPepCountMap.containsKey(proteins[j])) {
					count++;
				}
			}

			if (count > 0) {
				String stripeSequence = quanPeps[i].getSequence();
				double[] pepIntensity = new double[quanPeps[i].getIntensity().length];
				for (int j = 0; j < pepIntensity.length; j++) {
					pepIntensity[j] = quanPeps[i].getIntensity()[j] / (double) count;
				}

				for (int j = 0; j < proteins.length; j++) {
					if (razorPepCountMap.containsKey(proteins[j])) {
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

								if (totalPepCountMap.containsKey(genome)) {
									totalPepCountMap.get(genome).add(stripeSequence);
								} else {
									HashSet<String> set = new HashSet<String>();
									set.add(stripeSequence);
									totalPepCountMap.put(genome, set);
								}

								if (totalPepCountMap.containsKey(proteins[j])) {
									totalPepCountMap.get(proteins[j]).add(stripeSequence);
								} else {
									HashSet<String> set = new HashSet<String>();
									set.add(stripeSequence);
									totalPepCountMap.put(proteins[j], set);
								}
							}
						} else {
							// Handle non-MAG proteins (host proteins)
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
						}
					}
				}
			}
		}

		String[] proteins = proteinTotalIntenMap.keySet().toArray(String[]::new);
		Arrays.sort(proteins, new Comparator<String>() {
			@Override
			public int compare(String o1, String o2) {
				if (proteinTotalIntenMap.get(o1) > proteinTotalIntenMap.get(o2))
					return -1;
				if (proteinTotalIntenMap.get(o1) < proteinTotalIntenMap.get(o2))
					return 1;
				return 0;
			}
		});

		HashMap<String, Integer> proteinIdMap = new HashMap<String, Integer>();
		for (int i = 0; i < proteins.length; i++) {
			proteinIdMap.put(proteins[i], i);
		}

		try (PrintWriter idWriter = new PrintWriter(new File(magSageParams.getResult(), "id.tsv"))) {
			idWriter.println("Raw file\tUnique peptide count\tPSM count");

			String[] rawNames = fileIdMap.keySet().toArray(new String[fileIdMap.size()]);
			Arrays.sort(rawNames, new Comparator<String>() {
				@Override
				public int compare(String o1, String o2) {
					return fileIdMap.get(o1) - fileIdMap.get(o2);
				}
			});

			for (int i = 0; i < rawNames.length; i++) {
				idWriter.println(rawNames[i] + "\t" + fileSequenceSets[i].size() + "\t" + filePsmCount[i]);
			}
		} catch (IOException e) {
			System.out.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": error in writing identification information to " + new File(magSageParams.getResult(), "id.tsv"));
			LOGGER.error(getTaskName() + ": error in writing identification information to "
					+ new File(magSageParams.getResult(), "id.tsv"), e);
		}

		updateProgress(85, 100);

		HashMap<String, ArrayList<Integer>> proPepIdMap = new HashMap<String, ArrayList<Integer>>();
		HashMap<String, ArrayList<Integer>> proPepRazorIdMap = new HashMap<String, ArrayList<Integer>>();
		HashMap<String, Double> proScoreMap = new HashMap<String, Double>();
		try {
			this.final_pep_txt = new File(magSageParams.getDbSearchResultFile(), "final_peptides.tsv");

			PrintWriter writer = new PrintWriter(this.final_pep_txt);

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

			String[] idenNames = fileIdMap.keySet().toArray(String[]::new);
			Arrays.sort(idenNames, new Comparator<String>() {
				@Override
				public int compare(String o1, String o2) {
					return fileIdMap.get(o1) - fileIdMap.get(o2);
				}
			});

			for (int i = 0; i < idenNames.length; i++) {
				titlesb.append("Identification type " + idenNames[i]).append("\t");
			}

			for (int i = 0; i < idenNames.length; i++) {
				titlesb.append("MS/MS count " + idenNames[i]).append("\t");
			}

			titlesb.append("Intensity").append("\t");

			for (int i = 0; i < idenNames.length; i++) {
				titlesb.append("Intensity " + idenNames[i]).append("\t");
			}

			titlesb.append("Reverse").append("\t");
			titlesb.append("Potential contaminant").append("\t");
			titlesb.append("id").append("\t");
			titlesb.append("Protein group IDs").append("\t");
			titlesb.append("MS/MS Count");

			writer.println(titlesb);

			int id = 0;
			for (SagePeptide peptide : peptideList) {

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

						if (proScoreMap.containsKey(pro)) {
							proScoreMap.put(pro, proScoreMap.get(pro) * peptide.getPEP());
						} else {
							proScoreMap.put(pro, peptide.getPEP());
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
			LOGGER.error(getTaskName() + ": error in writing the peptides to  " + final_pep_txt, e);
			System.err.println(format.format(new Date()) + "\t" + getTaskName() + ": error in writing the peptides to  "
					+ final_pep_txt);
		}

		LOGGER.info(getTaskName() + ": peptide report has been exported to " + final_pep_txt);
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": peptide report has been exported to "
				+ final_pep_txt);

		if (!final_pep_txt.exists() || final_pep_txt.length() == 0) {
			return false;
		}

		updateProgress(87, 100);

		LOGGER.info(getTaskName() + ": exporting protein report started");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": exporting protein report started");

		ArrayList<MetaProtein> proteinList = new ArrayList<MetaProtein>();
		try {

			this.final_pro_txt = new File(magSageParams.getDbSearchResultFile(), "final_proteins.tsv");

			PrintWriter writer = new PrintWriter(this.final_pro_txt);

			StringBuilder titlesb = new StringBuilder();
			titlesb.append("Protein IDs").append("\t");
			titlesb.append("Peptide counts (all)").append("\t");
			titlesb.append("Peptide counts (razor)").append("\t");
			titlesb.append("Number of PSMs").append("\t");
			titlesb.append("PEP").append("\t");
			titlesb.append("Intensity").append("\t");
			
			String[] idenNames = fileIdMap.keySet().toArray(String[]::new);
			Arrays.sort(idenNames, new Comparator<String>() {
				@Override
				public int compare(String o1, String o2) {
					return fileIdMap.get(o1) - fileIdMap.get(o2);
				}
			});
			
			for (int i = 0; i < idenNames.length; i++) {
				titlesb.append("Intensity " + idenNames[i]).append("\t");
			}

			titlesb.append("Reverse").append("\t");
			titlesb.append("Potential contaminant").append("\t");
			titlesb.append("id").append("\t");
			titlesb.append("Peptide IDs").append("\t");
			titlesb.append("Peptide is razor");

			writer.println(titlesb);

			DecimalFormat df = FormatTool.getDFE4();
			for (int i = 0; i < proteins.length; i++) {

				StringBuilder sb = new StringBuilder();
				sb.append(proteins[i]).append("\t");
				sb.append(totalPepCountMap.get(proteins[i]).size()).append("\t");
				sb.append(razorPepCountMap.get(proteins[i]).size()).append("\t");
				sb.append(psmCountMap.get(proteins[i])).append("\t");
				sb.append(proScoreMap.get(proteins[i])).append("\t");

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
						psmCountMap.get(proteins[i]), proScoreMap.get(proteins[i]), intensity);
				proteinList.add(mp);
			}

			writer.close();

		} catch (IOException e) {
			LOGGER.error(getTaskName() + ": error in writing final protein result to " + this.final_pro_txt.getName(),
					e);
			System.err.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": error in writing final protein result to " + this.final_pro_txt.getName());
		}

		LOGGER.info(getTaskName() + ": protein report has been exported to " + final_pro_txt.getName());
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": protein report has been exported to "
				+ final_pro_txt.getName());

		updateProgress(89, 100);

		// Functional annotation and taxonomy analysis (merged from exportTaxaFunc)
		LOGGER.info(getTaskName() + ": functional annotation started");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": functional annotation started");

		File funcFile = new File(magSageParams.getDbSearchResultFile(), "functional_annotation");
		if (!funcFile.exists()) {
			funcFile.mkdir();
		}

		File funTsv = new File(funcFile, "functions.tsv");
		File funReportTsv = new File(funcFile, "functions_report.tsv");
		File final_pro_xml_file = new File(funcFile, "functions.xml");

		HashMap<String, Integer> proTaxIdMap = new HashMap<String, Integer>();
		MetaProtein[] metapros = proteinList.toArray(new MetaProtein[proteinList.size()]);
		
		String[] intensityTitles = fileList.toArray(new String[fileList.size()]);
		MetaProteinAnnoMag[] mpas = funcAnnoSqlite(metapros, intensityTitles);
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
			LOGGER.error(getTaskName() + ": error in reading protein function information from "
					+ final_pro_xml_file.getName(), e);
			System.err.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": error in reading protein function information from " + final_pro_xml_file.getName());
			return false;
		}

		if (!funTsv.exists() || !funReportTsv.exists()) {
			System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": writing results to "
					+ funTsv.getName() + " started");
			try {
				funReader.exportTsv(funTsv);
				funReader.exportTsvReport(funReportTsv);
			} catch (IOException e) {
				LOGGER.error(getTaskName() + ": error in writing protein function information to " + funTsv.getName(),
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
			} catch (IOException e) {
				LOGGER.error(getTaskName() + ": error in writing protein function information to " + html.getName(), e);
				System.err.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": error in writing protein function information to " + html.getName());
			}
			System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": writing results to "
					+ html.getName() + " finished");
		}

		updateProgress(91, 100);

		// taxonomy analysis
		LOGGER.info(getTaskName() + ": taxonomy analysis started");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": taxonomy analysis started");

		File taxFile = new File(magSageParams.getDbSearchResultFile(), "taxonomy_analysis");
		if (!taxFile.exists()) {
			taxFile.mkdir();
		}

		MagFuncSearcher searcher = null;
		try {
			searcher = new MagFuncSearcher(magSageParams.getMicroDb(), mrp.getToolPath(ToolType.FUNCDB));
		} catch (NumberFormatException | SQLException | IOException e) {
			LOGGER.error(getTaskName() + ": error in genome analysis, task failed", e);
			System.err.println(
					format.format(new Date()) + "\t" + getTaskName() + ": error in genome analysis, task failed");
			return false;
		}

		// Use the genome intensity map from the exportReport method
		HashMap<String, String[]> genomeTaxaMap = null;
		try {
			String[] genomes = genomeIntensityMap.keySet().toArray(new String[genomeIntensityMap.size()]);
			genomeTaxaMap = searcher.matchGenome(genomes);
		} catch (IOException e) {
			LOGGER.error(getTaskName() + ": Error in searching the genome2taxa database", e);
			System.err.println(
					format.format(new Date()) + "\t" + getTaskName() + ": error in searching the genome2taxa database");
			return false;
		}

		File genomeHtmlFile = new File(taxFile, "Genome.html");
		try {
			HMGenomeHtmlWriter htmlWriter = new HMGenomeHtmlWriter(genomeHtmlFile);
			htmlWriter.write(genomeTaxaMap);
			htmlWriter.close();
		} catch (IOException e) {
			LOGGER.error("Error in writing genome information to " + genomeHtmlFile, e);
			System.err.println(
					format.format(new Date()) + "\t" + "Error in writing genome information to " + genomeHtmlFile);
			return false;
		}

		// Generate detailed taxonomy reports (Taxa.tsv and Taxa_report.tsv)
		TaxonomyRanks[] mainRanks = TaxonomyRanks.getMainRanks();
		double[] cellularIntensity = new double[intensityTitles.length];
		HashMap<String, String[]>[] taxRankMaps = new HashMap[mainRanks.length];
		HashMap<String, double[]>[] taxIntensityMaps = new HashMap[mainRanks.length];
		HashMap<String, Double>[] taxScoreMaps = new HashMap[mainRanks.length];
		for (int i = 0; i < taxRankMaps.length; i++) {
			taxRankMaps[i] = new HashMap<String, String[]>();
			taxIntensityMaps[i] = new HashMap<String, double[]>();
			taxScoreMaps[i] = new HashMap<String, Double>();
		}

		File genomeFile = new File(taxFile, "Genome.tsv");
		try (PrintWriter genomeWriter = new PrintWriter(genomeFile)) {
			StringBuilder genomeTitlesb = new StringBuilder();
			genomeTitlesb.append("Genome\t");
			genomeTitlesb.append("Minus_log(p-value)\t");
			genomeTitlesb.append("PSM count\t");
			genomeTitlesb.append("Peptide count\t");
			genomeTitlesb.append("Razor peptide count\t");
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

			for (String exp : intensityTitles) {
				if (exp.startsWith("Intensity ")) {
					genomeTitlesb.append("\t").append(exp);
				} else {
					genomeTitlesb.append("\tIntensity ").append(exp);
				}
			}

			genomeWriter.println(genomeTitlesb);

			String[] genomes = genomeIntensityMap.keySet().toArray(new String[genomeIntensityMap.size()]);
			for (int i = 0; i < genomes.length; i++) {
				String[] taxaContent = genomeTaxaMap.get(genomes[i]);
				if (genomeIntensityMap.containsKey(genomes[i])) {
					double[] genomeIntensity = genomeIntensityMap.get(genomes[i]);
					for (int j = 0; j < genomeIntensity.length; j++) {
						cellularIntensity[j] += genomeIntensity[j];
					}

					double score = genomeScoreMap.get(genomes[i]);
					StringBuilder genomeTaxaSb = new StringBuilder();
					genomeTaxaSb.append(genomes[i]).append("\t");

					genomeTaxaSb.append(FormatTool.getDF2().format(-Math.log10(score))).append("\t");

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

					for (int j = 1; j < taxaContent.length; j++) {
						genomeTaxaSb.append("\t").append(taxaContent[j]);
					}

					for (int j = 0; j < genomeIntensity.length; j++) {
						genomeTaxaSb.append("\t").append(genomeIntensity[j]);
					}
					genomeWriter.println(genomeTaxaSb);

					// Process taxonomy ranks for aggregated reports
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
							taxScoreMaps[j].put(taxon, taxScoreMaps[j].get(taxon) - Math.log10(score));
						} else {
							taxScoreMaps[j].put(taxon, -Math.log10(score));
						}
					}
				}
			}
		} catch (FileNotFoundException e) {
			LOGGER.error("Error in writing genome information to " + genomeFile, e);
			System.err.println(format.format(new Date()) + "\t" + "Error in writing genome information to " + genomeFile);
			return false;
		}

		// Generate Taxa.tsv and Taxa_report.tsv files
		File taxaFile = new File(taxFile, "Taxa.tsv");
		File taxaReportFile = new File(taxFile, "Taxa_report.tsv");
		try (PrintWriter taxaWriter = new PrintWriter(taxaFile);
			 PrintWriter taxaReportWriter = new PrintWriter(taxaReportFile)) {

			StringBuilder taxaTitlesb1 = new StringBuilder();
			taxaTitlesb1.append("Name\t");
			taxaTitlesb1.append("Minus_log(p-value)\t");
			taxaTitlesb1.append("Rank");

			for (int i = 0; i < mainRanks.length; i++) {
				taxaTitlesb1.append("\t").append(mainRanks[i].getName());
			}
			for (String exp : intensityTitles) {
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
			for (String exp : intensityTitles) {
				if (exp.startsWith("Intensity ")) {
					taxaTitlesb2.append("\t").append(exp);
				} else {
					taxaTitlesb2.append("\tIntensity ").append(exp);
				}
			}
			
			taxaWriter.println(taxaTitlesb1);
			taxaReportWriter.println(taxaTitlesb2);

			// Write taxonomy data for Taxa.tsv
			for (int i = 0; i < taxRankMaps.length; i++) {
				for (String taxon : taxRankMaps[i].keySet()) {
					String[] taxaContent = taxRankMaps[i].get(taxon);
					double[] intensity = taxIntensityMaps[i].get(taxon);

					StringBuilder taxaSb = new StringBuilder();
					taxaSb.append(taxon).append("\t");

					double taxonScore = taxScoreMaps[i].get(taxon);
					taxaSb.append(FormatTool.getDF2().format(taxonScore)).append("\t");

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

			// Write taxonomy data for Taxa_report.tsv
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

		} catch (FileNotFoundException e) {
			LOGGER.error("Error in writing taxa information to " + taxaFile, e);
			System.err.println(format.format(new Date()) + "\t" + "Error in writing taxa information to " + taxaFile);
			return false;
		}

		LOGGER.info(getTaskName() + ": taxonomy analysis finished");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": taxonomy analysis finished");

		updateProgress(95, 100);

		updateMessage(getTaskName() + ": exporting report");

		LOGGER.info(getTaskName() + ": exporting report started");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": exporting report started");

		// Generate JavaScript tree data for interactive taxonomy visualization
		File jsOutput = new File(magSageParams.getReportDir() + "\\reports\\metamap\\data\\tree.js");
		if (!jsOutput.getParentFile().exists()) {
			jsOutput.getParentFile().mkdirs();
		}

		if (!jsOutput.exists() || jsOutput.length() == 0) {
			MetaTreeHandler.exportJS(taxRankMaps, taxIntensityMaps, cellularIntensity, intensityTitles, jsOutput);
		}

		File summaryMetaFile = new File(magSageParams.getDbSearchResultFile(), "metainfo.tsv");
		try {
			magSageParams.getMetadata().exportMetadata(summaryMetaFile);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

		// Generate summary report first
		MetaPythonReportTask reportTask = new MetaPythonReportTask(magSageParams.getThreadCount(), mrp);
		this.exportSummary(reportTask);
		
		// Add individual report tasks
		boolean reportSuccess = true;
		
		// Add summary report with defensive check
		File summaryReportFile = new File(magSageParams.getDbSearchResultFile(), "final_summary.tsv");
		if (summaryReportFile.exists() && summaryReportFile.length() > 0) {
			File report_ID_summary = new File(reportHtmlDir, "report_ID_summary.html");
			try {
				if (summaryMetaFile != null && summaryMetaFile.exists()) {
					reportSuccess &= reportTask.addTask(MetaPythonReportTask.ID_SUMMARY, summaryReportFile, report_ID_summary, summaryMetaFile);
				} else {
					reportSuccess &= reportTask.addTask(MetaPythonReportTask.ID_SUMMARY, summaryReportFile, report_ID_summary);
				}
				LOGGER.info(getTaskName() + ": ID summary report added successfully");
			} catch (Exception e) {
				LOGGER.warn(getTaskName() + ": Failed to add ID summary report, skipping: " + e.getMessage());
				System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": Failed to add ID summary report, skipping");
			}
		} else {
			LOGGER.warn(getTaskName() + ": Summary report file not found or empty, skipping ID summary report");
			System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": Summary report file not found or empty, skipping ID summary report");
		}

		// Add peptide report with defensive check
		if (final_pep_txt.exists() && final_pep_txt.length() > 0) {
			File report_peptides_summary = new File(reportHtmlDir, "report_peptides_summary.html");
			try {
				if (summaryMetaFile != null && summaryMetaFile.exists()) {
					reportSuccess &= reportTask.addTask(MetaPythonReportTask.PEPTIDE, final_pep_txt, report_peptides_summary, summaryMetaFile);
				} else {
					reportSuccess &= reportTask.addTask(MetaPythonReportTask.PEPTIDE, final_pep_txt, report_peptides_summary);
				}
				LOGGER.info(getTaskName() + ": Peptide report added successfully");
			} catch (Exception e) {
				LOGGER.warn(getTaskName() + ": Failed to add peptide report, skipping: " + e.getMessage());
				System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": Failed to add peptide report, skipping");
			}
		} else {
			LOGGER.warn(getTaskName() + ": Peptide report file not found or empty, skipping peptide report");
			System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": Peptide report file not found or empty, skipping peptide report");
		}

		// Add protein report with defensive check
		if (final_pro_txt.exists() && final_pro_txt.length() > 0) {
			File report_proteinGroups_summary = new File(reportHtmlDir, "report_proteinGroups_summary.html");
			try {
				if (summaryMetaFile != null && summaryMetaFile.exists()) {
					reportSuccess &= reportTask.addTask(MetaPythonReportTask.PROTEIN, final_pro_txt, report_proteinGroups_summary, summaryMetaFile);
				} else {
					reportSuccess &= reportTask.addTask(MetaPythonReportTask.PROTEIN, final_pro_txt, report_proteinGroups_summary);
				}
				LOGGER.info(getTaskName() + ": Protein report added successfully");
			} catch (Exception e) {
				LOGGER.warn(getTaskName() + ": Failed to add protein report, skipping: " + e.getMessage());
				System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": Failed to add protein report, skipping");
			}
		} else {
			LOGGER.warn(getTaskName() + ": Protein report file not found or empty, skipping protein report");
			System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": Protein report file not found or empty, skipping protein report");
		}

		if (taxaReportFile.exists() && taxaReportFile.length() > 0) {
			File report_taxonomy_summary = new File(reportHtmlDir, "report_taxonomy_summary.html");
			try {
				if (summaryMetaFile != null && summaryMetaFile.exists()) {
					reportSuccess &= reportTask.addTask(MetaPythonReportTask.TAXONOMY, taxaReportFile, report_taxonomy_summary, summaryMetaFile);
				} else {
					reportSuccess &= reportTask.addTask(MetaPythonReportTask.TAXONOMY, taxaReportFile, report_taxonomy_summary);
				}
				LOGGER.info(getTaskName() + ": Taxonomy report added successfully");
			} catch (Exception e) {
				LOGGER.warn(getTaskName() + ": Failed to add taxonomy report, skipping: " + e.getMessage());
				System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": Failed to add taxonomy report, skipping");
			}
		} else {
			LOGGER.warn(getTaskName() + ": Taxa_report.tsv file not found or empty, skipping taxonomy report");
			System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": Taxa_report.tsv file not found or empty, skipping taxonomy report");
			System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": Note: Taxonomy report generation may need to be implemented");
		}

		// Add function report with defensive check
		if (funReportTsv.exists() && funReportTsv.length() > 0) {
			File report_function_summary = new File(reportHtmlDir, "report_function_summary.html");
			try {
				if (summaryMetaFile != null && summaryMetaFile.exists()) {
					reportSuccess &= reportTask.addTask(MetaPythonReportTask.FUNCTION, funReportTsv, report_function_summary, summaryMetaFile);
				} else {
					reportSuccess &= reportTask.addTask(MetaPythonReportTask.FUNCTION, funReportTsv, report_function_summary);
				}
				LOGGER.info(getTaskName() + ": Function report added successfully");
			} catch (Exception e) {
				LOGGER.warn(getTaskName() + ": Failed to add function report, skipping: " + e.getMessage());
				System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": Failed to add function report, skipping");
			}
		} else {
			LOGGER.warn(getTaskName() + ": Function report file not found or empty, skipping function report");
			System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": Function report file not found or empty, skipping function report");
		}

		if (!reportSuccess) {
			LOGGER.warn(getTaskName() + ": Some report tasks could not be added");
			System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": Some report tasks could not be added");
		}

		// Run the report generation task
		Thread reportThread = new Thread(reportTask);
		reportThread.start();
		try {
			reportThread.join();
		} catch (InterruptedException e) {
			LOGGER.error(getTaskName() + ": Report generation thread was interrupted", e);
			System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": Report generation thread was interrupted");
			Thread.currentThread().interrupt();
			return false;
		}

		try {
			boolean finish = reportTask.get();
			if (finish) {
				LOGGER.info(getTaskName() + ": report generation successes");
				System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": report generation successes");

				// Final progress update
				updateProgress(100, 100);
				return true;
			}

		} catch (InterruptedException | ExecutionException e) {
			LOGGER.error(getTaskName() + ": error in the report generation", e);
			System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": error in the report generation");
			
			// Update progress on error
			return false;
		}

		// Final progress update even if report generation fails
		updateProgress(100, 100);
		return true;
	}


}
