package uk.ac.quadram.metapilot.task.dda.mag.task;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileFilter;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
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

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import uk.ac.quadram.metapilot.dbSearch.sage.SageTask;
import uk.ac.quadram.metapilot.quant.flashLFQ.FlashLfqTask;
import uk.ac.quadram.metapilot.core.prodb.FastaManager;
import uk.ac.quadram.metapilot.dbSearch.sage.SagePeptide;
import uk.ac.quadram.metapilot.dbSearch.sage.SagePSM;
import uk.ac.quadram.metapilot.dbSearch.sage.SageResultReader;
import uk.ac.quadram.metapilot.task.MetaPythonReportTask;
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
	protected boolean exportReport() {
		LOGGER.info(getTaskName() + ": exporting peptide report started");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": exporting peptide report started");

		// Build file list and mapping
		ArrayList<String> fileList = new ArrayList<String>();
		HashMap<String, Integer> fileIdMap = new HashMap<String, Integer>();
		File expDesignFile = new File(rawDir, "ExperimentalDesign.tsv");

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

		// Read PSMs
		SagePSM[] psms = null;
		try {
			SageResultReader sageReader = new SageResultReader(combinedFile);
			psms = sageReader.getPSMs();
		} catch (Exception e) {
			LOGGER.error(getTaskName() + ": error in reading the PSMs from " + combinedFile, e);
			System.err.println(format.format(new Date()) + "\t" + getTaskName() + ": error in reading the PSMs from "
					+ combinedFile);
			return false;
		}

		// Read quantified peptides
		SageResultReader quanPepReader;
		String quanMode = magSageParams.getQuantitativeMode();
		if (quanMode.equals(MetaConstants.labelFree)) {
			quanPepReader = new SageResultReader(combinedFile);
		} else {
			LOGGER.error(getTaskName() + ": Isobaric quantification not yet supported for Sage workflow");
			System.err.println(format.format(new Date()) + "\t" + getTaskName() + ": Isobaric quantification not yet supported");
			return false;
		}

		SagePeptide[] quanPeps = quanPepReader.getMetaPeptides();

		// Build intensity title array
		String[] intensityTitle = fileList.toArray(String[]::new);

		// Call the common method to process peptides and generate reports
		return processPeptideListAndGenerateReports(quanPeps, psms, fileIdMap, fileList, quanMode, intensityTitle,
				getTaskName());
	}

	@Override
	protected String getVersion() {
		// TODO Auto-generated method stub
		return null;
	}
}
