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
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.concurrent.ExecutionException;

import org.apache.logging.log4j.Logger;

import uk.ac.quadram.metapilot.dbSearch.sage.SageTask;
import uk.ac.quadram.metapilot.quant.flashLFQ.FlashLfqTask;
import uk.ac.quadram.metapilot.core.prodb.DecoyShuffler;
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
import uk.ac.quadram.metapilot.task.mag.MarkerGenomeSelector;
import uk.ac.quadram.metapilot.task.mag.MarkerGenomeSelector.RepresentativeLevel;
import uk.ac.quadram.metapilot.task.dda.mag.par.DDAMagSageParameters;
import uk.ac.quadram.metapilot.task.io.pep.AbstractMetaPeptideReader;
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

	// Logger is inherited from MetaPilotMagTask -> MetaPilotTask as an instance field
	private static final String taskName = "DDA MAG database search by Sage workflow";

	public DDAMagSageTask(DDAMagSageParameters magSageParams, MetaResourceParameters mrp, ProgressBar progressBar,
			ProgressBar secondBar, Label statusLabel) {
		super(magSageParams, mrp, progressBar, secondBar);

		// Store the direct reference to DDAMagSageParameters
		this.magSageParams = magSageParams;

		// Initialize status label from constructor parameter (passed directly from GUI)
		this.statusLabel = statusLabel;

		// Initialize the magDb field from parameters
		MagDbItem selectedMagDb = magSageParams.getMagDbItem();
		if (selectedMagDb != null) {
			this.magDb = selectedMagDb;
			logger.info("Initialized MAG database: {} (version: {})", selectedMagDb.getCatalogueID(),
					selectedMagDb.getVersionList() != null && !selectedMagDb.getVersionList().isEmpty()
							? selectedMagDb.getVersionList().get(0).getCatalogueVersion()
							: "unknown");

			// Validate that the MAG database has the required files
			File hapFastaFile = selectedMagDb.getHapFastaFile();
			if (hapFastaFile == null || !hapFastaFile.exists()) {
				logger.error("MAG database HAP FASTA file is null or does not exist: {}",
						hapFastaFile != null ? hapFastaFile.getAbsolutePath() : "null");
				logger.error("Repository file: {}", selectedMagDb.getRepositoryFile());
				logger.error("Current file: {}", selectedMagDb.getCurrentFile());
				throw new IllegalArgumentException("MAG database HAP FASTA file is required and must exist");
			}
			logger.info("MAG database HAP FASTA file validated: {}", hapFastaFile.getAbsolutePath());

			File currentDbFile = selectedMagDb.getCurrentFile();
			if (currentDbFile == null || !currentDbFile.exists()) {
				logger.error("MAG database file is null or does not exist: {}",
						currentDbFile != null ? currentDbFile.getAbsolutePath() : "null");
				logger.error("Repository file: {}", selectedMagDb.getRepositoryFile());
				logger.error("Current file: {}", selectedMagDb.getCurrentFile());
				throw new IllegalArgumentException("MAG database file is required and must exist");
			}

			logger.info("MAG database file validated: {}", currentDbFile.getAbsolutePath());
			this.metaPar.setMicroDb(selectedMagDb.getCurrentFile().getAbsolutePath());
		}

		// Initialize Sage task
		this.sageExe = new File(mrp.getToolPath(ToolType.SAGE));
		this.sageTask = new SageTask(sageExe.getAbsolutePath());

		// Validate that the MAG database has the required files
		if (selectedMagDb != null) {
			File hapFastaFile = selectedMagDb.getHapFastaFile();
			if (hapFastaFile == null || !hapFastaFile.exists()) {
				logger.error("MAG database HAP FASTA file is null or does not exist: {}",
						hapFastaFile != null ? hapFastaFile.getAbsolutePath() : "null");
				logger.error("Repository file: {}", selectedMagDb.getRepositoryFile());
				logger.error("Current file: {}", selectedMagDb.getCurrentFile());
				throw new IllegalArgumentException("MAG database HAP FASTA file is required and must exist");
			}
			logger.info("MAG database HAP FASTA file validated: {}", hapFastaFile.getAbsolutePath());
		} else {
			logger.error("MAG database item is null in parameters");
			throw new IllegalArgumentException("MAG database item is required");
		}

		// Initialize using the DDAMagSageParameters directly
		String[] rawFiles = magSageParams.getMetadata().getRawFiles();
		this.mzMLFiles = new String[rawFiles.length];
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
				baseName = fileName.substring(0, fileName.length() - ".mzml".length());
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
			logger.info("Created Sage MAG search directory: {}", sageMagSearchFile.getAbsolutePath());
		} else {
			logger.info("Sage MAG search directory already exists: {}", sageMagSearchFile.getAbsolutePath());
		}
		magSageParams.setDbSearchResultFile(sageMagSearchFile);

		createSeparateFolders();

		this.sageExe = new File(mrp.getToolPath(ToolType.SAGE));
		this.msConvertExe = new File(mrp.getToolPath(ToolType.MSCONVERT));

		logger.info("Sage path from MetaResourceParameters: {}", mrp.getToolPath(ToolType.SAGE));
		logger.info("Sage executable file: {}", sageExe.getAbsolutePath());
		logger.info("Sage exists: {}", sageExe.exists());

		logger.info("MSConvert path from MetaResourceParameters: {}", mrp.getToolPath(ToolType.MSCONVERT));
		logger.info("MSConvert executable file: {}", msConvertExe.getAbsolutePath());
		logger.info("MSConvert exists: {}", msConvertExe.exists());
	}

	@Override
	public void forceStop() {
		logger.info(getTaskName() + ": Force stop requested");
		logger.info("[{}] Force stop requested", getTaskName());

		// Cancel Sage task if running
		if (sageTask != null) {
			logger.info(getTaskName() + ": Stopping Sage identification task");
			// SageTask doesn't have forceStop method, but we can interrupt any running
			// threads
			// The task will be stopped when this main task is cancelled
		}

		// Cancel FlashLFQ task if running
		if (flashLfqTask != null && !flashLfqTask.isDone()) {
			logger.info(getTaskName() + ": Cancelling FlashLFQ quantification task");
			flashLfqTask.cancel();
		}

		// Cancel report generation task if running
		if (reportTask != null && !reportTask.isDone()) {
			logger.info(getTaskName() + ": Cancelling report generation task");
			reportTask.cancel(true);
		}

		// Cancel executor service if running
		if (executorService != null && !executorService.isShutdown()) {
			logger.info(getTaskName() + ": Shutting down executor service");
			executorService.shutdownNow();
		}

		// Cancel report copy thread if running
		if (reportCopyThread != null && reportCopyThread.isAlive()) {
			logger.info(getTaskName() + ": Interrupting report copy thread");
			reportCopyThread.interrupt();
		}

		// Cancel report generation thread if running
		if (reportThread != null && reportThread.isAlive()) {
			logger.info(getTaskName() + ": Interrupting report generation thread");
			reportThread.interrupt();
		}

		// Cancel this task itself
		this.cancel(true);

		logger.info(getTaskName() + ": Force stop completed");
		logger.info("[{}] Force stop completed", getTaskName());
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
	protected Boolean call() throws Exception {
		updateMessage(getTaskName());
		updateProgress(0, 100);

		logger.info(getTaskName() + ": started");
		logger.info("[{}] Started", getTaskName());

		// Log detailed parameter information
		logger.info("Workflow parameters:");
		logger.info("  Result directory: " + magSageParams.getResult());
		logger.info("  Quantitative mode: " + magSageParams.getQuantitativeMode());
		logger.info("  Thread count: " + magSageParams.getThreadCount());
		logger.info("  MS2 scan mode: " + magSageParams.getMs2ScanMode());
		if (magSageParams.getMetadata() != null) {
			String[] rawFiles = magSageParams.getMetadata().getRawFiles();
			logger.info("  Raw files count: " + (rawFiles != null ? rawFiles.length : 0));
			if (rawFiles != null) {
				for (int i = 0; i < rawFiles.length; i++) {
					logger.info("    Raw file " + (i + 1) + ": " + rawFiles[i]);
				}
			}
		}

		// Log resource information
		logger.info("Resource parameters:");
		logger.info("  Sage path: " + mrp.getToolPath(MetaResourceParameters.ToolType.SAGE));

		// Update main progress and keep second bar active
		updateProgress(5, 100);
		if (statusLabel != null) {
			Platform.runLater(() -> statusLabel.setText(getTaskName() + ": Initializing..."));
		}
		Platform.runLater(() -> secondBar.setProgress(ProgressBar.INDETERMINATE_PROGRESS));

		// Step 1: Spectra Conversion
		logger.info("Step 1: Starting spectra format conversion to mzML...");
		if (statusLabel != null) {
			Platform.runLater(() -> statusLabel.setText(getTaskName() + ": Converting spectra format..."));
		}
		updateProgress(10, 100);
		Platform.runLater(() -> secondBar.setProgress(ProgressBar.INDETERMINATE_PROGRESS));

		boolean convert = convertSpectraToMzML();
		if (!convert) {
			logger.error(getTaskName() + ": error in spectra format conversion");
			System.err
					.println(format.format(new Date()) + "\t" + getTaskName() + ": error in spectra format conversion");
			updateMessage(getTaskName() + ": Failed - Spectra conversion error");
			Platform.runLater(() -> secondBar.setProgress(0.0));
			return false;
		}
		logger.info("Step 1: Spectra format conversion completed successfully");

		// Step 2: High Abundant Database Search
		logger.info("Step 2: Starting high abundant database search...");
		updateProgress(10, 100);
		updateMessage(getTaskName() + ": Searching high abundant database...");
		Platform.runLater(() -> secondBar.setProgress(ProgressBar.INDETERMINATE_PROGRESS));

		logger.info("[{}] Searching high abundant database started", getTaskName());

		// Initialize separateResultFolders
		logger.info("Initializing result directories in: {}", magSageParams.getDbSearchResultFile().getAbsolutePath());

		boolean searchHap = this.separateSearchHap();
		if (!searchHap) {
			logger.error(getTaskName() + ": error in searching high abundant database");
			return false;
		}

		// Step 3: Protein Extraction from HAP Results
		logger.info("Step 3: Starting protein extraction from HAP results...");
		updateProgress(40, 100);
		boolean extract = this.fastaExtract4MarkerSelector();
		if (!extract) {
			logger.error(getTaskName() + ": error in protein extraction");
			return false;
		}

		// Step 4: Database Search
		logger.info("Step 4: Starting database search through sample-specific database...");
		updateProgress(65, 100);
		boolean search = this.separateSearch();
		if (!search) {
			logger.error(getTaskName() + ": error in database search");
			return false;
		}

		// Step 5: Combine Results
		logger.info("Step 5: Starting combination of separate search results...");
		updateProgress(75, 100);
		boolean combine = this.combineResults();
		if (!combine) {
			logger.error(getTaskName() + ": error in combining separate search results");
			return false;
		}

		// Step 6: Export Report
		logger.info("Step 6: Starting report generation...");
		updateProgress(85, 100);
		boolean report = this.exportReport();
		if (!report) {
			logger.error(getTaskName() + ": error in report generation");
			return false;
		}

		// Update progress for completion
		updateProgress(100, 100);
		updateMessage(getTaskName() + ": Completed successfully!");
		Platform.runLater(() -> secondBar.setProgress(0.0));

		logger.info(getTaskName() + ": finished - all workflow steps completed successfully");
		return true;
	}

	/**
	 * Convert spectra files to mzML format for Sage processing using multiple
	 * threads
	 * @throws IOException 
	 */
	private boolean convertSpectraToMzML() throws IOException {
		logger.info("Starting spectra conversion to mzML format using MSConvert");

		// Validate MSConvert executable
		if (msConvertExe == null || !msConvertExe.exists()) {
			logger.error("MSConvert executable not found: {}",
					msConvertExe != null ? msConvertExe.getAbsolutePath() : "null");
			return false;
		}

		String[] rawFiles = magSageParams.getMetadata().getRawFiles();
		if (rawFiles == null || rawFiles.length == 0) {
			logger.error("No raw files found in metadata");
			return false;
		}

		// Ensure spectra directory exists
		File spectraDir = magSageParams.getSpectraFile();
		if (!spectraDir.exists()) {
			spectraDir.mkdirs();
			logger.info("Created spectra directory: {}", spectraDir.getAbsolutePath());
		}

		// mzMLFiles array is already initialized in constructor

		// Collect files that need conversion
		List<Integer> filesToConvert = new ArrayList<>();
		for (int i = 0; i < rawFiles.length; i++) {
			String rawFile = rawFiles[i];
			File rawFileObj = new File(rawFile);

			if (!rawFileObj.exists()) {
				logger.error("Raw file does not exist: {}", rawFile);
				return false;
			}

			// Check if the mzML file already exists (mzMLFiles[i] was set in constructor)
			File mzMLFile = new File(mzMLFiles[i]);

			if (mzMLFile.exists() && mzMLFile.length() > 0) {
				logger.info("mzML file already exists, skipping conversion: {}", mzMLFile.getAbsolutePath());
			} else {
				filesToConvert.add(i);
			}
		}

		// If no files need conversion, return success
		if (filesToConvert.isEmpty()) {
			logger.info("All mzML files already exist, no conversion needed");

		} else {

			// Use multiple threads for conversion
			int threadCount = Math.min(filesToConvert.size(), magSageParams.getThreadCount());
			ExecutorService executor = Executors.newFixedThreadPool(threadCount);

			logger.info("Processing {} files using {} threads", filesToConvert.size(), threadCount);

			// Submit conversion/copy tasks
			for (int fileIndex : filesToConvert) {
				executor.submit(() -> {
					String rawFile = rawFiles[fileIndex];
					File rawFileObj = new File(rawFile);
					String fileName = rawFileObj.getName().toLowerCase();

					try {
						String convertedPath = null;

						if (fileName.endsWith(".raw")) {
							logger.info("Converting .raw file to mzML: {}", rawFile);
							convertedPath = MSConvertor.raw2mzML(msConvertExe.getAbsolutePath(), rawFile);
						} else if (fileName.endsWith(".d")) {
							logger.info("Converting .d file to mzML: {}", rawFile);
							convertedPath = MSConvertor.d2mzML(msConvertExe.getAbsolutePath(), rawFile,
									spectraDir.getAbsolutePath());
						} else {
							logger.error(
									"Unsupported file format: {}. Only .raw and .d files are supported for conversion.",
									fileName);
							return;
						}

						if (convertedPath != null) {
							// Verify the converted file exists
							File convertedFile = new File(convertedPath);
							if (convertedFile.exists() && convertedFile.length() > 0) {
								// Move the converted file to the spectra directory
								File targetFile = new File(mzMLFiles[fileIndex]);
								try {
									Files.move(convertedFile.toPath(), targetFile.toPath(),
											StandardCopyOption.REPLACE_EXISTING);
									logger.info("Successfully moved converted file from {} to {}", convertedPath,
											mzMLFiles[fileIndex]);
								} catch (IOException e) {
									logger.error("Failed to move converted file from {} to {}: {}", convertedPath,
											mzMLFiles[fileIndex], e.getMessage());
								}
							} else {
								logger.error("Converted mzML file does not exist or is empty: {}", convertedPath);
							}
						} else {
							logger.error("Failed to convert file to mzML: {}", rawFile);
						}
					} catch (Exception e) {
						logger.error("Error processing file {}: {}", rawFile, e.getMessage());
					}
				});
			}

			// Shutdown executor and wait for completion
			executor.shutdown();

			try {
				// Wait for all conversions to complete (with timeout)
				boolean finished = executor.awaitTermination(filesToConvert.size() * 4, TimeUnit.HOURS);
				if (finished) {
					logger.info("All file conversions completed successfully");
				} else {
					logger.error("File conversion timeout - some files may not have been converted");
					executor.shutdownNow();
					return false;
				}
			} catch (InterruptedException e) {
				logger.error("File conversion was interrupted: {}", e.getMessage());
				executor.shutdownNow();
				return false;
			}
		}

		File rtFile = new File(spectraDirFile, "rt.txt");
		File ms2CountFile = new File(spectraDirFile, "ms2Count.txt");
		
		if (!rtFile.exists() || rtFile.length() == 0 || !ms2CountFile.exists() || ms2CountFile.length() == 0) {

			PrintWriter rtWriter = new PrintWriter(rtFile);
			PrintWriter ms2CountWriter = new PrintWriter(ms2CountFile);

			rtWriter.println("FileName\tScannum\tRt");
			ms2CountWriter.println("FileName\tMS2 count");

			Pattern scanPattern = Pattern.compile("scan=(\\d+)");
			Pattern rtPattern = Pattern.compile("name=\"scan start time\"\\s+value=\"([^\"]+)\"");

			// Verify all files were converted successfully
			for (int i = 0; i < mzMLFiles.length; i++) {
				if (mzMLFiles[i] == null) {
					logger.error("Failed to convert file: {}", rawFiles[i]);
					rtWriter.close();
					ms2CountWriter.close();
					
					return false;
				}

				String fileName = (new File(mzMLFiles[i])).getName();
				fileName = fileName.substring(0, fileName.length() - ".mzML".length());

				int count = 0;
				int scannum = 0;
				boolean isMS2 = false;
				BufferedReader reader = new BufferedReader(new FileReader(mzMLFiles[i]));
				String line = null;
				while ((line = reader.readLine()) != null) {
					line = line.trim();
					if (line.startsWith("<spectrum index=")) {
						Matcher matcher = scanPattern.matcher(line);
						if (matcher.find()) {
							scannum = Integer.parseInt(matcher.group(1));
						}
						isMS2 = false;
					} else if (line.endsWith("name=\"MSn spectrum\" value=\"\"/>")) {
						count++;
						isMS2 = true;
					} else if (line.startsWith("<cvParam cvRef=\"MS\" accession=\"MS:1000016\" name=\"scan start time\"")) {
						Matcher matcher = rtPattern.matcher(line);
						double rt = -1;
						if (matcher.find()) {
							rt = Double.parseDouble(matcher.group(1));
						}
						if (isMS2) {
							rtWriter.println(fileName + "\t" + scannum + "\t" + rt);
						}
					}
				}
				reader.close();
				ms2CountWriter.println(fileName + "\t" + count);
			}

			rtWriter.close();
			ms2CountWriter.close();
		}

		logger.info("Spectra conversion completed successfully. Converted {} files to mzML format.", rawFiles.length);
		return true;
	}

	/**
	 * Perform high abundant database search using Sage
	 */
	private boolean separateSearchHap() {
		logger.info("Starting separate search for high abundant database using Sage");

		// Get thread configuration
		int threadCount = magSageParams.getThreadCount();
		logger.info("Thread configuration: {}", threadCount);

		File hapFastaFile = this.magDb.getHapFastaFile();
		if (hapFastaFile == null || !hapFastaFile.exists()) {
			logger.error("HAP FASTA file is null or does not exist for MAG database: {}", this.magDb.getCatalogueID());
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
			logger.info("Found {} files requiring new search out of {} total files", filesNeedingSearch,
					this.mzMLFiles.length);

			for (int i = 0; i < this.mzMLFiles.length; i++) {
				File hapFile = new File(separateResultFolders[i], "hap");
				if (!hapFile.exists()) {
					hapFile.mkdirs();
				}

				File resultFile = new File(hapFile, "results.sage.tsv");
				if (!resultFile.exists() || resultFile.length() == 0) {
					File mzMLFile = new File(mzMLFiles[i]);
					if (mzMLFile.exists()) {
						logger.info("  Searching file {} with Sage", i + 1);

						// Use the new searchMzML method with DDAMagSageParameters
						sageTask.searchMzML(new String[] { mzMLFile.getAbsolutePath() }, hapFile.getAbsolutePath(),
								hapFastaFile.getAbsolutePath(), magSageParams, false); // isQuant = false for HAP search
					} else {
						logger.warn("mzML file {} was not found, skipping", mzMLFiles[i]);
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
			logger.info("Starting WatchService monitoring for high abundant database search");
			startWatchServiceLoop(separateResultFolders, "hap", filter, filesNeedingSearch, 10, 40);

			// Calculate optimal thread pool based on RAM (1 thread per 16GB)
			int threadPool = maxTaskCount > filesNeedingSearch ? filesNeedingSearch : maxTaskCount;
			logger.info("Starting Sage execution with {} thread pool for {} tasks", threadPool, this.mzMLFiles.length);
			sageTask.run(threadPool, separateResultFolders.length * 12);
		} else {
			logger.info("No new searches required, all results already exist");
		}

		// Verify search results
		logger.info("Verifying search results for {} folders", separateResultFolders.length);
		for (int i = 0; i < this.separateResultFolders.length; i++) {
			File hapFile = new File(separateResultFolders[i], "hap");
//			logger.info("  Cleaning and verifying results for folder {}: {}", i+1, hapFile.getAbsolutePath());
//			SageTask.cleanSageResultFolder(hapFile, true);

			File resultFile = new File(hapFile, "results.sage.tsv");
			if (!resultFile.exists() || resultFile.length() == 0) {
				logger.error("Sage search result was not found from {}", hapFile.getAbsolutePath());
				return false;
			}
			logger.info("  Results verified successfully for folder {}", i + 1);
		}

		logger.info("High abundant database search completed successfully for all {} files",
				separateResultFolders.length);

		// Stop the watch service thread
		stopWatchServiceThread();

		return true;
	}

	/**
	 * Extract protein sequences from HAP search results
	 */
	protected boolean fastaExtract4Iterator() {
		logger.info("Starting FASTA extraction for {} result folders", separateResultFolders.length);

		for (int i = 0; i < this.separateResultFolders.length; i++) {
			File hapFile = new File(separateResultFolders[i], "hap");
			File dbFile = new File(hapFile, separateResultFolders[i].getName() + ".fasta");

			if (dbFile.exists() && dbFile.length() > 0) {
				logger.info("Sample-specific database for {} already exists", separateResultFolders[i].getName());
				continue;
			}

			File spectraFile = new File(hapFile, "results.sage.tsv");
			HashMap<String, HashSet<String>> genomeMap = new HashMap<String, HashSet<String>>();
			HashSet<String> pepSet = new HashSet<String>();

			try (BufferedReader hapReader = new BufferedReader(new FileReader(spectraFile))) {
				String headerLine = hapReader.readLine(); // Read header
				String[] header = headerLine.split("\t");

				// Find the index of "proteins" column
				int peptidesIndex = -1;
				int proteinsIndex = -1;
				int labelIndex = -1;
				int peptideQIndex = -1;
				for (int j = 0; j < header.length; j++) {
					if (header[j].toLowerCase().equals("peptide")) {
						peptidesIndex = j;
					} else if (header[j].toLowerCase().equals("proteins")) {
						proteinsIndex = j;
					} else if (header[j].toLowerCase().equals("label")) {
						labelIndex = j;
					} else if (header[j].toLowerCase().equals("peptide_q")) {
						peptideQIndex = j;
					}
				}

				if (proteinsIndex == -1 || labelIndex == -1 || peptideQIndex == -1) {
					logger.error("Could not find 'proteins' 'label' or 'peptide_q' column in header: {}", headerLine);
					continue;
				}

				logger.info("Found peptide column at index {}: {}", peptidesIndex, header[peptidesIndex]);
				logger.info("Found proteins column at index {}: {}", proteinsIndex, header[proteinsIndex]);
				logger.info("Found label column at index {}: {}", labelIndex, header[labelIndex]);
				logger.info("Found peptide_q column at index {}: {}", peptideQIndex, header[peptideQIndex]);

				String line;
				while ((line = hapReader.readLine()) != null) {
					String[] cs = line.split("\t");
					if (cs.length > labelIndex) {
						if (!cs[labelIndex].equals("1")) {
							continue;
						}
					} else {
						continue;
					}

					if (cs.length > peptideQIndex) {
						if (Double.parseDouble(cs[peptideQIndex]) > 0.01) {
							continue;
						}
					} else {
						continue;
					}

					if (cs.length > proteinsIndex) {
						String[] pros = cs[proteinsIndex].split(";");
						for (String pro : pros) {
							if (pro.startsWith(magDb.getIdentifier())) {
								String genome = pro.split("_")[0];
								if (genomeMap.containsKey(genome)) {
									genomeMap.get(genome).add(cs[peptidesIndex]);
								} else {
									HashSet<String> set = new HashSet<String>();
									set.add(cs[peptidesIndex]);
									genomeMap.put(genome, set);
								}

							}
						}
						pepSet.add(cs[peptidesIndex]);
					}
				}
			} catch (IOException e) {
				logger.error("Error reading Sage search result file from {}", spectraFile.getAbsolutePath(), e);
				return false;
			}

			logger.info("Sample-specific database for {} - Found {} peptides and {} genomes",
					separateResultFolders[i].getName(), pepSet.size(), genomeMap.size());
			// Create sample-specific database
			File hapGenomeFile = new File(hapFile, "hap_genomes.tsv");
			HashSet<String> genomeSet = applyOccamsRazor(genomeMap, pepSet, psmCoverThreshold, hapGenomeFile);

			if (genomeSet.size() == 0) {
				logger.info("No genomes were found from {}", separateResultFolders[i].getName());
			}

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
				logger.error("Error writing sample-specific database to {}", dbFile.getAbsolutePath(), e);
				return false;
			}
		}

		return true;
	}

	/**
	 * Extract protein sequences from HAP search results using MarkerGenomeSelector
	 */
	protected boolean fastaExtract4MarkerSelector() {
		logger.info("Starting FASTA extraction (marker selector) for {} result folders", separateResultFolders.length);

		File currentFile = magDb.getCurrentFile();
		File markerFasta = new File(currentFile, "4markers.fasta");
		if (markerFasta == null || !markerFasta.exists()) {
			logger.error("Marker FASTA file was not found: {}", markerFasta != null ? markerFasta.getAbsolutePath() : "null");
			return false;
		}

		File genomesTaxaTsv = new File(currentFile, "genomes-taxa.tsv");
		if (genomesTaxaTsv == null || !genomesTaxaTsv.exists()) {
			logger.error("Genomes taxa TSV was not found: {}", genomesTaxaTsv != null ? genomesTaxaTsv.getAbsolutePath() : "null");
			return false;
		}

		MarkerGenomeSelector selector = new MarkerGenomeSelector();
		int stepThreshold = 1;
		int maximumSelected = 300;
		double coverage = 0.85;

		for (int i = 0; i < this.separateResultFolders.length; i++) {
			File hapFile = new File(separateResultFolders[i], "hap");
			File dbFile = new File(hapFile, separateResultFolders[i].getName() + ".fasta");

			if (dbFile.exists() && dbFile.length() > 0) {
				logger.info("Sample-specific database for {} already exists", separateResultFolders[i].getName());
				continue;
			}

			File spectraFile = new File(hapFile, "results.sage.tsv");
			if (!spectraFile.exists() || spectraFile.length() == 0) {
				logger.error("Sage search result was not found from {}", spectraFile.getAbsolutePath());
				return false;
			}

			File evidenceFile = new File(hapFile, "marker_genome_evidence.tsv");
			File rankedFile = new File(hapFile, "marker_genome_ranked.tsv");
			File selectedFile = new File(hapFile, "marker_genome_selected.tsv");

			selector.peptideTestWithStrategy(spectraFile.getAbsolutePath(), markerFasta, genomesTaxaTsv,
					evidenceFile.getAbsolutePath(), rankedFile.getAbsolutePath(), selectedFile.getAbsolutePath(), coverage,
					stepThreshold, maximumSelected);

			HashSet<String> genomeSet = new HashSet<>();
			try (BufferedReader reader = new BufferedReader(new FileReader(selectedFile))) {
				String line = reader.readLine(); // header
				while ((line = reader.readLine()) != null) {
					line = line.trim();
					if (line.isEmpty() || line.startsWith("#")) {
						break;
					}
					String[] parts = line.split("\t");
					if (parts.length > 0 && !parts[0].trim().isEmpty()) {
						genomeSet.add(parts[0].trim());
					}
				}
			} catch (IOException e) {
				logger.error("Error reading marker selected genomes from {}", selectedFile.getAbsolutePath(), e);
				return false;
			}

			logger.info("Marker selector for {} - Selected {} genomes", separateResultFolders[i].getName(), genomeSet.size());

			if (genomeSet.isEmpty()) {
				logger.info("No genomes were found from {}", separateResultFolders[i].getName());
			}

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
				logger.error("Error writing sample-specific database to {}", dbFile.getAbsolutePath(), e);
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
				logger.info("[{}] Sage search result already exists in {}", getTaskName(), separateResultFolders[i].getAbsolutePath());

				logger.info(getTaskName() + ": Sage search result already exists in "
						+ separateResultFolders[i].getAbsolutePath());
			} else {
				File hapFile = new File(separateResultFolders[i], "hap");
				File db = new File(hapFile, separateResultFolders[i].getName() + ".fasta");
				if (!db.exists()) {
					logger.error("[{}] sample-specific database {} was not found", getTaskName(), db);
					logger.error(getTaskName() + ": sample-specific database" + db + " was not found");
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
				File resultFile = new File(separateResultFolders[i], "results.sage.tsv");
				if (!resultFile.exists() || resultFile.length() == 0) {
					if (mzMLFile.exists()) {
						if (db.exists() && db.length() > 0) {
							// Use the new searchMzML method with DDAMagSageParameters for sample-specific
							// database search
							sageTask.searchMzML(new String[] { mzMLFile.getAbsolutePath() },
									separateResultFolders[i].getAbsolutePath(), db.getAbsolutePath(), magSageParams,
									false);
						} else {
							logger.warn(
									"[{}] no protein sequences were in {}, which means no peptides/proteins can be identified from {}",
									getTaskName(), db, separateResultFolders[i].getName());

							logger.error(getTaskName() + ": no protein sequences were in " + db
									+ ", which means no peptides/proteins can be identifed from "
									+ separateResultFolders[i].getName());
						}
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
			logger.info("Starting WatchService monitoring for sample-specific database search");
			startWatchServiceLoop(separateResultFolders, "", filter, taskCount, 40, 65);

			// Calculate optimal thread pool based on RAM (1 thread per 16GB)
			int threadPool = maxTaskCount > taskCount ? taskCount : maxTaskCount;
			logger.info("Starting Sage execution with {} thread pool for {} tasks", threadPool, taskCount);
			sageTask.run(threadPool, this.separateResultFolders.length * 12);
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
				
				File db = new File(hapFile, separateResultFolders[i].getName() + ".fasta");
				if(db.length()>0) {
					logger.error("Sage search result was not found in {}", separateResultFolders[i].getAbsolutePath());
					return false;
				}
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
		combinedFile = new File(magSageParams.getDbSearchResultFile(), "combined");
		if (!combinedFile.exists()) {
			combinedFile.mkdirs();
		}
		
		File combinedFastaFile = new File(combinedFile, "combined.fasta");
		if (!combinedFastaFile.exists() || combinedFastaFile.length() == 0) {
			
			File combinedResultTsv = new File(combinedFile, "combined_results.sage.tsv");
			HashSet<String> totalProteinSet = new HashSet<String>();
			if(!combinedResultTsv.exists()) {
				try (PrintWriter writer = new PrintWriter(combinedResultTsv)) {
					boolean wroteHeader = false;
					for (int i = 0; i < separateResultFolders.length; i++) {
						File pepResultFile = new File(separateResultFolders[i], "results.sage.tsv");
						if (!pepResultFile.exists() || pepResultFile.length() == 0) {
							File fastaFile = new File(new File(separateResultFolders[i], "hap"),
									separateResultFolders[i].getName() + ".fasta");
							if (fastaFile.length() == 0) {
								logger.info("No proteins found from {}, skip...", separateResultFolders[i].getName());
							}
							continue;
						}
						try (BufferedReader combineTsvReader = new BufferedReader(new FileReader(pepResultFile))) {
							String line = combineTsvReader.readLine();
							if (line == null) {
								continue;
							}

							String[] title = line.split("\t");
							int proid = -1;
							int qvalueid = -1;
							for (int j = 0; j < title.length; j++) {
								if (title[j].equals("proteins")) {
									proid = j;
								} else if (title[j].equals("spectrum_q")) {
									qvalueid = j;
								}
							}
							if (!wroteHeader) {
								writer.println(line);
								wroteHeader = true;
							}
							while ((line = combineTsvReader.readLine()) != null) {
								String[] cs = line.split("\t");
								if (Double.parseDouble(cs[qvalueid]) > 0.01 || cs[proid].contains("REV_")
										|| cs[proid].contains("rev_")) {
									continue;
								}
								writer.println(line);
								String[] pros = cs[proid].split(";");
								for (int j = 0; j < pros.length; j++) {
									totalProteinSet.add(pros[j]);
								}
							}
						}
					}
				} catch (IOException e) {
					logger.error(getTaskName() + ": error in writing combined results file " + combinedResultTsv, e);
					return false;
				}
			} else {
				try (BufferedReader combineTsvReader = new BufferedReader(new FileReader(combinedResultTsv))) {
					String line = combineTsvReader.readLine();
					String[] title = line.split("\t");
					int proid = -1;
					int qvalueid = -1;
					for (int j = 0; j < title.length; j++) {
						if (title[j].equals("proteins")) {
							proid = j;
						} else if (title[j].equals("spectrum_q")) {
							qvalueid = j;
						}
					}
					while ((line = combineTsvReader.readLine()) != null) {
						String[] cs = line.split("\t");
						if (Double.parseDouble(cs[qvalueid]) > 0.01 || cs[proid].contains("REV_")
								|| cs[proid].contains("rev_")) {
							continue;
						}
						String[] pros = cs[proid].split(";");
						for (int j = 0; j < pros.length; j++) {
							totalProteinSet.add(pros[j]);
						}
					}
				} catch (IOException e) {
					logger.error(getTaskName() + ": error in reading combined results file " + combinedResultTsv, e);
					return false;
				}
			}
			
			logger.info("Before filtering: {} proteins were found from the separate search result", totalProteinSet.size());
			
			// Create combined database
			File currentFile = magDb.getCurrentFile();
			File markerFasta = new File(currentFile, "4markers.fasta");
			if (markerFasta == null || !markerFasta.exists()) {
				logger.error("Marker FASTA file was not found: {}",
						markerFasta != null ? markerFasta.getAbsolutePath() : "null");
				return false;
			}

			File genomesTaxaTsv = new File(currentFile, "genomes-taxa.tsv");
			if (genomesTaxaTsv == null || !genomesTaxaTsv.exists()) {
				logger.error("Genomes taxa TSV was not found: {}",
						genomesTaxaTsv != null ? genomesTaxaTsv.getAbsolutePath() : "null");
				return false;
			}

			File evidenceFile = new File(combinedFile, "marker_genome_evidence.tsv");
			File rankedFile = new File(combinedFile, "marker_genome_ranked.tsv");
			File selectedFile = new File(combinedFile, "marker_genome_selected.tsv");

			MarkerGenomeSelector selector = new MarkerGenomeSelector();
			selector.setRepresentativeLevel("CLASS");
			int stepThreshold = 1;
			double coverage = 0.9;
			int maximumSelected = 1000;
			selector.peptideTestWithStrategy(combinedResultTsv.getAbsolutePath(), markerFasta, genomesTaxaTsv,
					evidenceFile.getAbsolutePath(), rankedFile.getAbsolutePath(), selectedFile.getAbsolutePath(),
					coverage, stepThreshold, maximumSelected);

			HashSet<String> genomeSet = new HashSet<>();
			try (BufferedReader reader = new BufferedReader(new FileReader(selectedFile))) {
				String line = reader.readLine(); // header
				while ((line = reader.readLine()) != null) {
					line = line.trim();
					if (line.isEmpty() || line.startsWith("#")) {
						break;
					}
					String[] parts = line.split("\t");
					if (parts.length > 0 && !parts[0].trim().isEmpty()) {
						genomeSet.add(parts[0].trim());
					}
				}
			} catch (IOException e) {
				logger.error("Error reading marker selected genomes from {}", selectedFile.getAbsolutePath(), e);
				return false;
			}
			
			logger.info("MarkerGenomeSelector select {} genomes", genomeSet.size());
			
			
			HashMap<String, HashSet<String>> genomePepMap = new HashMap<String, HashSet<String>>();
			HashMap<String, HashSet<String>> proPepMap = new HashMap<String, HashSet<String>>();
			HashSet<String> pepSet = new HashSet<String>();
			
			try (BufferedReader combineTsvReader = new BufferedReader(new FileReader(combinedResultTsv))) {
				String line = combineTsvReader.readLine();
				String[] title = line.split("\t");
				int proid = -1;
				int pepid = -1;
				for (int j = 0; j < title.length; j++) {
					if (title[j].equals("proteins")) {
						proid = j;
					} else if (title[j].equals("peptide")) {
						pepid = j;
					}
				}
				while ((line = combineTsvReader.readLine()) != null) {
					String[] cs = line.split("\t");
					String[] pros = cs[proid].split(";");
					for (int j = 0; j < pros.length; j++) {
						int loc = pros[j].indexOf("_");
						if (loc > 0) {
							String genome = pros[j].substring(0, pros[j].indexOf("_"));
							if (genomeSet.contains(genome)) {
								genomePepMap.computeIfAbsent(genome, k -> new HashSet<>()).add(cs[pepid]);
								proPepMap.computeIfAbsent(pros[j], k -> new HashSet<>()).add(cs[pepid]);
								pepSet.add(cs[pepid]);
							}
						}
					}
				}
			} catch (IOException e) {
				logger.error(getTaskName() + ": error in reading combined results file " + combinedResultTsv, e);
				return false;
			}			
			
			File razorGenomeFile = new File(combinedFile, "razorGenomes.tsv");
			HashSet<String> selecteGenomeSet = this.applyOccamsRazor(genomePepMap, pepSet, 0.99, razorGenomeFile);
			Iterator<String> it = proPepMap.keySet().iterator();
			while (it.hasNext()) {
				String pro = it.next();
				String genome = pro.substring(0, pro.indexOf("_"));
				if (!selecteGenomeSet.contains(genome)) {
					it.remove();
				}
			}

			File razorProFile = new File(combinedFile, "razorProteins.tsv");
			HashSet<String> selecteProSet = this.applyOccamsRazor(proPepMap, pepSet, 1.0, razorProFile);

			HashSet<String> containSet = new HashSet<String>();
			try (PrintWriter writer = new PrintWriter(combinedFastaFile)) {
				for (int i = 0; i < separateResultFolders.length; i++) {
					File hapFile = new File(separateResultFolders[i], "hap");
					File dbFile = new File(hapFile, separateResultFolders[i].getName() + ".fasta");
					try (BufferedReader reader = new BufferedReader(new FileReader(dbFile))) {
						String line = null;
						String proName = "";
						boolean use = false;
						StringBuilder sb = new StringBuilder();
						while ((line = reader.readLine()) != null) {
							if (line.startsWith(">")) {

								if (sb.length() > 0) {
									String shuffleProName = ">shuf_" + proName;
									String shuffle = DecoyShuffler.generateShuffledDecoy(sb.toString(),
											this.sageTask.getEnzymeCleavageSites(this.magSageParams.getEnzyme()), 42);
									writer.println(shuffleProName);
									writer.println(shuffle);
									sb = new StringBuilder();
								}

								proName = line.substring(1, line.indexOf(" "));
								if (selecteProSet.contains(proName) && !containSet.contains(proName)) {
									use = true;
									writer.println(line);
									containSet.add(proName);
								} else {
									use = false;
								}

							} else {
								if (use) {
									writer.println(line);
									sb.append(line);
								}
							}
						}

						if (sb.length() > 0) {
							String shuffleProName = ">shuf_" + proName;
							String shuffle = DecoyShuffler.generateShuffledDecoy(sb.toString(),
									this.sageTask.getEnzymeCleavageSites(this.magSageParams.getEnzyme()), 42);
							writer.println(shuffleProName);
							writer.println(shuffle);
						}

						reader.close();
					} catch (IOException e) {
						logger.error(getTaskName() + ": error in reading database file " + dbFile, e);
						System.err.println(format.format(new Date()) + "\t" + getTaskName()
								+ ": error in reading database file " + dbFile);
					}
				}

				if (metaPar.isAppendHostDb()) {
					BufferedReader hostReader = new BufferedReader(new FileReader(this.metaPar.getHostDB()));
					String line = null;
					while ((line = hostReader.readLine()) != null) {
						writer.println(line);
					}
					hostReader.close();
				}
				writer.close();

				magSageParams.setCurrentDb(combinedFastaFile.getAbsolutePath());

			} catch (IOException e) {
				logger.error(getTaskName() + ": error in writing the combined database file to " + combinedFastaFile,
						e);
				System.err.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": error in writing the combined database file to " + combinedFastaFile);
			}
		}

		// Run final Sage search against combined database
		File pepResultFile = new File(combinedFile, "results.sage.tsv");
		File quanResultFile = new File(combinedFile, "lfq.tsv");

		if (!pepResultFile.exists() && !quanResultFile.exists()) {
			// Run Sage search against combined database
			try {
				sageTask.searchMzML(mzMLFiles, combinedFile.getAbsolutePath(), combinedFastaFile.getAbsolutePath(),
						magSageParams, true);
			} catch (Exception e) {
				logger.error(getTaskName() + ": error in running Sage search against combined database", e);
				logger.error("[{}] error in running Sage search against combined database", getTaskName());
				return false;
			}

			logger.info(getTaskName() + ": searching combined database started");
			logger.info("[{}] searching combined database started", getTaskName());

			logger.info("Starting Sage execution");
			sageTask.run(separateResultFolders.length * 12);

		} else {
			logger.info("[{}] Sage search result already exists in {}", getTaskName(), magSageParams.getDbSearchResultFile());
			logger.info(
					getTaskName() + ": Sage search result already exists in " + magSageParams.getDbSearchResultFile());
		}

		// Check if results exist
		if (pepResultFile.exists() && pepResultFile.length() > 0 && quanResultFile.exists()
				&& quanResultFile.length() > 0) {
			logger.info(getTaskName() + ": combined database search finished successfully");
			logger.info("[{}] combined database search finished successfully", getTaskName());
			return true;
		} else {
			logger.info(getTaskName() + ": combined database search failed");
			logger.error("[{}] combined database search failed", getTaskName());
			return false;
		}
	}

	/**
	 * Export final reports
	 */
	protected boolean exportReport() {
		logger.info(getTaskName() + ": exporting peptide report started");
		logger.info("[{}] exporting peptide report started", getTaskName());

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
				logger.error(getTaskName() + ": error in reading the experimental design from " + expDesignFile, e);
				System.err.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": error in reading the experimental design from " + expDesignFile);
			}
		} else {
			File[] files = magSageParams.getDbSearchResultFile().listFiles();
			for (int i = 0; i < files.length; i++) {
				if (files[i].isDirectory() && !files[i].getName().equals("combined")
						&& (new File(files[i], "hap")).exists()) {
					fileIdMap.put(files[i].getName(), fileIdMap.size());
					fileList.add(files[i].getName());
				}
			}
		}

		// Read quantified peptides
		SageResultReader quanPepReader;
		String quanMode = magSageParams.getQuantitativeMode();
		if (quanMode.equals(MetaConstants.labelFree)) {
			quanPepReader = new SageResultReader(combinedFile);
		} else {
			logger.error(getTaskName() + ": Isobaric quantification not yet supported for Sage workflow");
			System.err.println(
					format.format(new Date()) + "\t" + getTaskName() + ": Isobaric quantification not yet supported");
			return false;
		}

		SagePeptide[] quanPeps = quanPepReader.getMetaPeptides();
		SagePSM[] psms = quanPepReader.getPSMs();

		logger.info(getTaskName() + ": " + quanPeps.length + " peptides and " + psms.length + " PSMs were obtained");

		// Build intensity title array
		String[] intensityTitle = fileList.toArray(String[]::new);

		// Create and execute the post-analysis task directly
		DDAMagPostAnalysisTask postAnalysisTask = new DDAMagPostAnalysisTask(metaPar, mrp, secondBar, quanPeps, psms,
				fileIdMap, fileList, quanMode, intensityTitle, getTaskName(), AbstractMetaPeptideReader.SHUF, this);

		// Execute the task and wait for completion
		Thread taskThread = new Thread(postAnalysisTask);
		taskThread.setDaemon(true);
		taskThread.start();

		try {
			taskThread.join();
			boolean success = postAnalysisTask.get();

			if (success) {
				logger.info(getTaskName() + ": Post-analysis completed successfully");
				return true;
			} else {
				logger.error(getTaskName() + ": Post-analysis task failed");
				return false;
			}
		} catch (InterruptedException e) {
			logger.error(getTaskName() + ": Post-analysis task was interrupted", e);
			Thread.currentThread().interrupt();
			return false;
		} catch (ExecutionException e) {
			logger.error(getTaskName() + ": Post-analysis task failed with exception", e);
			return false;
		}
	}

	@Override
	protected String getVersion() {
		// TODO Auto-generated method stub
		return null;
	}
	public static void getGenomeCountFromFasta(String in) {
		HashSet<String> genomeSet = new HashSet<String>();
		try (BufferedReader reader = new BufferedReader(new FileReader(in))) {
			String line = null;
			while ((line = reader.readLine()) != null) {
				if (line.startsWith(">MGYG")) {
					String genomeName = line.substring(1, line.indexOf("_"));
					genomeSet.add(genomeName);
				}
			}
			reader.close();
		} catch (IOException e) {
		}
		System.out.println(genomeSet.size());
	}
	
	public static void comparePep(String in1, String in2, String out) {
		HashSet<String> s2 = new HashSet<String>();
		try (BufferedReader reader = new BufferedReader(new FileReader(in2))) {
			String line = reader.readLine();
			while ((line = reader.readLine()) != null) {
				String[] cs = line.split("\t");
				s2.add(cs[0]);
			}
			reader.close();
		} catch (IOException e) {
		}
		
		HashMap<String, HashSet<String>> genomePepMap1 = new HashMap<String, HashSet<String>>();
		HashMap<String, HashSet<String>> genomePepMap2 = new HashMap<String, HashSet<String>>();
		HashSet<String> s1 = new HashSet<String>();
		try (BufferedReader reader = new BufferedReader(new FileReader(in1))) {
			PrintWriter writer = new PrintWriter(out);
			String line = reader.readLine();
			writer.println(line);
			while ((line = reader.readLine()) != null) {
				String[] cs = line.split("\t");
				s1.add(cs[0]);
				if(!s2.contains(cs[0])) {
					writer.println(line);
					String[] pros = cs[2].split(";");
					for(int i=0;i<pros.length;i++) {
						String genome = pros[i].substring(0, pros[i].indexOf("_"));
						genomePepMap1.computeIfAbsent(genome, k -> new HashSet<>()).add(cs[0]);
					}
				} else {
					String[] pros = cs[2].split(";");
					for(int i=0;i<pros.length;i++) {
						String genome = pros[i].substring(0, pros[i].indexOf("_"));
						genomePepMap2.computeIfAbsent(genome, k -> new HashSet<>()).add(cs[0]);
					}
				}
			}
			reader.close();
			writer.close();
		} catch (IOException e) {
		}
		
		HashSet<String> allset = new HashSet<String>();
		allset.addAll(s1);
		allset.addAll(s2);
		int common = s1.size() + s2.size() - allset.size();
		System.out.println(s1.size() + "\t" + s2.size() + "\t" + allset.size() + "\t" + (s1.size() - common) + "\t"
				+ (s2.size() - common));
		
		for (String genome : genomePepMap1.keySet()) {
			if (!genomePepMap2.containsKey(genome)) {
				System.out.println(genome + "\t" + genomePepMap1.get(genome).size());
			}
		}
	}
	
	public static void main(String[] args) {
//		getGenomeCountFromFasta("Z:\\PRIDE_dataset\\PXD006118_3mix_type\\Run5\\MetaPilot\\sage_MAG_search\\combined\\combined.fasta");
//		getGenomeCountFromFasta("Z:\\PRIDE_dataset\\PXD006118_3mix_type\\Run1\\MetaPilot\\sage_MAG_search_rib\\Run1_C3_1952ng\\hap\\Run1_C3_1952ng.fasta");
		comparePep("Z:\\PRIDE_dataset\\PXD006118_3mix_type\\Run1\\MetaPilot\\sage_MAG_search\\combined\\lfq.tsv", 
				"Z:\\PRIDE_dataset\\PXD006118_3mix_type\\Run1\\refdb\\Run1_lfq.tsv", 
				"Z:\\PRIDE_dataset\\PXD006118_3mix_type\\Run1\\MetaPilot\\sage_MAG_search\\combined\\wrong.tsv");
	}
}
