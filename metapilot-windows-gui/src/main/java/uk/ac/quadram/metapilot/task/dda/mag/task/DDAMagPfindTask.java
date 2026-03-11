package uk.ac.quadram.metapilot.task.dda.mag.task;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileFilter;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.lang.management.ManagementFactory;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.ExecutionException;

import org.apache.logging.log4j.Logger;

import com.sun.management.OperatingSystemMXBean;

import uk.ac.quadram.metapilot.core.mod.IsobaricTag;
import uk.ac.quadram.metapilot.dbSearch.dbReducer.DBReducerTask;
import uk.ac.quadram.metapilot.dbSearch.open.io.OpenPsmTsvWriter;
import uk.ac.quadram.metapilot.dbSearch.pfind.PFindIsobaricTask;
import uk.ac.quadram.metapilot.dbSearch.pfind.PFindPSM;
import uk.ac.quadram.metapilot.dbSearch.pfind.PFindPsmReader;
import uk.ac.quadram.metapilot.dbSearch.pfind.PFindResultReader;
import uk.ac.quadram.metapilot.dbSearch.pfind.PFindTask;
import javafx.scene.control.ProgressBar;
import javafx.scene.control.Label;
import javafx.application.Platform;
import uk.ac.quadram.metapilot.quant.flashLFQ.FlashLfqPsmReader;
import uk.ac.quadram.metapilot.quant.flashLFQ.FlashLfqQuanPepReader;
import uk.ac.quadram.metapilot.quant.flashLFQ.FlashLfqQuanPeptide;
import uk.ac.quadram.metapilot.quant.flashLFQ.FlashLfqTask;
import uk.ac.quadram.metapilot.quant.flashLFQ.FlashLfqPsm;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters.ToolType;
import uk.ac.quadram.metapilot.task.dda.PFindSpectraConvertTask;
import uk.ac.quadram.metapilot.task.io.par.MetaOpenPar;
import uk.ac.quadram.metapilot.task.mag.MetaPilotMagTask;
import uk.ac.quadram.metapilot.task.par.MetaConstants;
import uk.ac.quadram.metapilot.task.dda.mag.par.DDAMagPfindParameters;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;
import uk.ac.quadram.metapilot.task.mag.MarkerGenomeSelector;

public class DDAMagPfindTask extends MetaPilotMagTask {

	private String[] mgfFiles;
	private String[] pf2Files;
	private File dbReducerExe;
	private File pFindExe;
	private boolean useDBReducer = false;

	// Fields to track running subtasks
	private PFindSpectraConvertTask spectraConvertTask;
	private DBReducerTask dbReducerTask;
	private PFindTask pFindTask;
	private FlashLfqTask flashLfqTask;
	private Thread reportCopyThread;

	// Progress tracking components
	private Label statusLabel;
	
	// Add direct reference to DDAMagPfindParameters
	private DDAMagPfindParameters magPfindParams;

	// Logger is inherited from MetaPilotMagTask -> MetaPilotTask as an instance field
	private static final String taskName = "DDA MAG database search by pFind workflow";

	public DDAMagPfindTask(DDAMagPfindParameters magPfindParams, MetaResourceParameters mrp, ProgressBar progressBar,
			ProgressBar secondBar, Label statusLabel) {
		super(magPfindParams, mrp, progressBar, secondBar);

		// Store the direct reference to DDAMagPfindParameters
		this.magPfindParams = magPfindParams;

		// Initialize status label from constructor parameter (passed directly from GUI)
		this.statusLabel = statusLabel;
		
		// Initialize the magDb field from parameters
		MagDbItem selectedMagDb = magPfindParams.getMagDbItem();
		if (selectedMagDb != null) {
			this.magDb = selectedMagDb;
			logger.info("Initialized MAG database: {} (version: {})", selectedMagDb.getCatalogueID(),
					selectedMagDb.getVersionList() != null && !selectedMagDb.getVersionList().isEmpty()
							? selectedMagDb.getVersionList().get(0).getCatalogueVersion()
							: "unknown");

			// Validate that the MAG database has the required files
			File hapFastaFile = selectedMagDb.getHapFastaFile(magPfindParams.getGenomeSelectorMode());
			if (hapFastaFile == null || !hapFastaFile.exists()) {
				logger.error("MAG database HAP FASTA file is null or does not exist: {}",
						hapFastaFile != null ? hapFastaFile.getAbsolutePath() : "null");
				logger.error("Repository file: {}", selectedMagDb.getRepositoryFile());
				logger.error("Current file: {}", selectedMagDb.getCurrentFile());
				throw new IllegalArgumentException("MAG database HAP FASTA file is required and must exist");
			}
			logger.info("MAG database HAP FASTA file validated: {}", hapFastaFile.getAbsolutePath());
			logger.info("Genome selector mode: {}", magPfindParams.getGenomeSelectorMode());
			
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
		} else {
			logger.error("MAG database item is null in parameters");
			throw new IllegalArgumentException("MAG database item is required");
		}
		// Initialize the result directory structure in the constructor
		File pfindMagSearchFile = new File(magPfindParams.getResult(), "pfind_MAG_search");
		if (!pfindMagSearchFile.exists()) {
			pfindMagSearchFile.mkdir();
			logger.info("Created pFind MAG search directory: {}", pfindMagSearchFile.getAbsolutePath());
		} else {
			logger.info("pFind MAG search directory already exists: {}", pfindMagSearchFile.getAbsolutePath());
		}
		magPfindParams.setDbSearchResultFile(pfindMagSearchFile);

		createSeparateFolders();

		// Initialize using the DDAMagPfindParameters directly
		String ms2Type = this.magPfindParams.getMs2ScanMode();
		String suffix = "";
		if (ms2Type.equals(MetaConstants.HCD_FTMS)) {
			suffix = "_HCDFT";
		} else if (ms2Type.equals(MetaConstants.HCD_ITMS)) {
			suffix = "_HCDIT";
		} else if (ms2Type.equals(MetaConstants.CID_FTMS)) {
			suffix = "_CIDFT";
		} else if (ms2Type.equals(MetaConstants.CID_ITMS)) {
			suffix = "_CIDIT";
		}
		this.mgfFiles = new String[rawFiles.length];
		this.pf2Files = new String[rawFiles.length];

		for (int i = 0; i < rawFiles.length; i++) {
			File mgfFile = new File(spectraDirFile, rawFileNames[i] + ".mgf");
			File pf2 = new File(spectraDirFile, rawFileNames[i] + suffix + ".pf2");
			mgfFiles[i] = mgfFile.getAbsolutePath();
			pf2Files[i] = pf2.getAbsolutePath();
		}

		this.dbReducerExe = new File(mrp.getToolPath(ToolType.DBREDUCER));
		this.pFindExe = new File(mrp.getToolPath(ToolType.PFIND));
		this.pFindTask = new PFindTask(pFindExe.getParentFile());
		
		logger.info("DBReducer path from MetaResourceParameters: {}", mrp.getToolPath(ToolType.DBREDUCER));
		logger.info("pFind path from MetaResourceParameters: {}", mrp.getToolPath(ToolType.PFIND));
		logger.info("DBReducer executable file: {}", dbReducerExe.getAbsolutePath());
		logger.info("pFind executable file: {}", pFindExe.getAbsolutePath());
		logger.info("DBReducer exists: {}", dbReducerExe.exists());
		logger.info("pFind exists: {}", pFindExe.exists());

		if (dbReducerExe.getName().equals("DBReducer.exe") && dbReducerExe.exists()
				&& magPfindParams.getMs2ScanMode().equals(MetaConstants.HCD_FTMS)) {
//			useDBReducer = true;
		}
	}

	protected void setTaskCount() {
		long totalMemorySize = ((OperatingSystemMXBean) ManagementFactory.getOperatingSystemMXBean())
				.getTotalMemorySize() / 1024 / 1024 / 1024;
		this.maxTaskCount = (int) (totalMemorySize / 5);
		if (this.maxTaskCount == 0) {
			this.maxTaskCount = 1;
		}

		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": the number of max memory is "
				+ totalMemorySize + " GB, the max number of tasks is set as " + maxTaskCount);
	}
	
	@Override
	public void forceStop() {
		logger.info(getTaskName() + ": Force stop requested");
		logger.info("[{}] Force stop requested", getTaskName());
		
		// Cancel spectra conversion task if running
		if (spectraConvertTask != null && !spectraConvertTask.isDone()) {
			logger.info(getTaskName() + ": Cancelling spectra conversion task");
			spectraConvertTask.cancel(true);
		}
		
		// Cancel DBReducer task if running
		if (dbReducerTask != null) {
			logger.info(getTaskName() + ": Stopping DBReducer task");
			// DBReducerTask doesn't have forceStop method, but we can interrupt any running threads
			// The task will be stopped when this main task is cancelled
		}
		
		// Cancel pFind task if running
		if (pFindTask != null) {
			logger.info(getTaskName() + ": Stopping pFind identification task");
			// PFindTask doesn't have forceStop method, but we can interrupt any running threads
			// The task will be stopped when this main task is cancelled
		}
		
		// Cancel FlashLFQ task if running
		if (flashLfqTask != null && !flashLfqTask.isDone()) {
			logger.info(getTaskName() + ": Cancelling FlashLFQ quantification task");
			flashLfqTask.cancel();
		}
		
		// Cancel watch service thread if running
		if (watchServiceThread != null && watchServiceThread.isAlive()) {
			logger.info(getTaskName() + ": Interrupting watch service thread");
			watchServiceThread.interrupt();
		}
		
		// Cancel report copy thread if running
		if (reportCopyThread != null && reportCopyThread.isAlive()) {
			logger.info(getTaskName() + ": Interrupting report copy thread");
			reportCopyThread.interrupt();
		}
		
		// Cancel this task itself
		this.cancel(true);
		
		logger.info(getTaskName() + ": Force stop completed");
		logger.info("[{}] Force stop completed", getTaskName());
	}

	@Override
	protected String getTaskName() {
		// TODO Auto-generated method stub
		return taskName;
	}

	@Override
	protected Logger getLogger() {
		// TODO Auto-generated method stub
		return logger;
	}

	@Override
	protected Boolean call() throws Exception {
		// TODO Auto-generated method stub
		updateMessage(getTaskName());
		updateProgress(0, 100);

		logger.info(getTaskName() + ": started");
		logger.info("[{}] started", getTaskName());

		// Log detailed parameter information
		logger.info("Workflow parameters:");
		logger.info("  Result directory: " + magPfindParams.getResult());
		logger.info("  Quantitative mode: " + magPfindParams.getQuantitativeMode());
		logger.info("  Thread count: " + magPfindParams.getThreadCount());
		logger.info("  MS2 scan mode: " + magPfindParams.getMs2ScanMode());
		logger.info("  Use DBReducer: " + useDBReducer);
		if (magPfindParams.getMetadata() != null) {
			String[] rawFiles = magPfindParams.getMetadata().getRawFiles();
			logger.info("  Raw files count: " + (rawFiles != null ? rawFiles.length : 0));
			if (rawFiles != null) {
				for (int i = 0; i < rawFiles.length; i++) {
					logger.info("    Raw file " + (i + 1) + ": " + rawFiles[i]);
				}
			}
		}

		// Log resource information
		logger.info("Resource parameters:");
		logger.info("  pFind path: " + mrp.getToolPath(MetaResourceParameters.ToolType.PFIND));
		logger.info("  DBReducer path: " + mrp.getToolPath(MetaResourceParameters.ToolType.DBREDUCER));
		logger.info("  MetaResourceParameters object: {}", mrp);
		logger.info("  MetaResourceParameters class: {}", mrp.getClass().getName());

		// Update main progress and keep second bar active
		updateProgress(5, 100);
		if (statusLabel != null) {
			Platform.runLater(() -> statusLabel.setText(getTaskName() + ": Initializing..."));
		}
		Platform.runLater(() -> secondBar.setProgress(ProgressBar.INDETERMINATE_PROGRESS)); // Indeterminate progress
																							// for continuous activity

		// Step 1: Spectra Conversion
		logger.info("Step 1: Starting spectra format conversion...");
		if (statusLabel != null) {
			Platform.runLater(() -> statusLabel.setText(getTaskName() + ": Converting spectra format..."));
		}
		updateProgress(10, 100);
		Platform.runLater(() -> secondBar.setProgress(ProgressBar.INDETERMINATE_PROGRESS));

		spectraConvertTask = new PFindSpectraConvertTask(magPfindParams, mrp, secondBar, logger);
		spectraConvertTask.execute();
		try {
			boolean convert = spectraConvertTask.get(); // Wait for completion and get result
			if (!convert) {
				logger.error(getTaskName() + ": error in spectra format conversion");
				System.err.println(
						format.format(new Date()) + "\t" + getTaskName() + ": error in spectra format conversion");
				// Update progress on error
				updateMessage(getTaskName() + ": Failed - Spectra conversion error");
				Platform.runLater(() -> secondBar.setProgress(0.0)); // Stop the indeterminate progress
				return false;
			}
			logger.info("Step 1: Spectra format conversion completed successfully");
		} catch (InterruptedException | ExecutionException e) {
			logger.error(getTaskName() + ": exception in spectra format conversion", e);
			System.err.println(
					format.format(new Date()) + "\t" + getTaskName() + ": exception in spectra format conversion");
			// Update progress on error
			updateMessage(getTaskName() + ": Failed - Spectra conversion exception");
			Platform.runLater(() -> secondBar.setProgress(0.0)); // Stop the indeterminate progress
			return false;
		}

		// Step 2: High Abundant Database Search
		logger.info("Step 2: Starting high abundant database search...");
		updateProgress(10, 100);
		updateMessage(getTaskName() + ": Searching high abundant database...");
		Platform.runLater(() -> secondBar.setProgress(ProgressBar.INDETERMINATE_PROGRESS));
		updateMessage("Performing high abundant protein search...");

		logger.info("[{}] searching high abundant database started", getTaskName());

		boolean searchHap = this.separateSearchHap();

		if (searchHap) {
			logger.info("Step 2: High abundant database search completed successfully");
			logger.info("[{}] searching high abundant database finished", getTaskName());
		} else {
			logger.error(getTaskName() + ": error in searching high abundant database");
			System.err.println(
					format.format(new Date()) + "\t" + getTaskName() + ": error in searching high abundant database");

			return false;
		}

		// Step 3: Protein Extraction from HAP Results
		logger.info("Step 3: Starting protein extraction from HAP results...");
		updateProgress(20, 100);
		updateMessage(getTaskName() + ": Extracting proteins from HAP results...");
		Platform.runLater(() -> secondBar.setProgress(ProgressBar.INDETERMINATE_PROGRESS));
		updateMessage("Extracting protein sequences...");
		updateProgress(40, 100);
		updateMessage("Extracting protein sequences...");

		System.out.println(
				format.format(new Date()) + "\t" + getTaskName() + ": extracting proteins from the HAP result started");

		boolean extract = this.fastaExtractWithGenomeSelectorMode();

		if (extract) {
			logger.info("Step 3: Protein extraction from HAP results completed successfully");
			logger.info("[{}] extracting proteins from the HAP result finished", getTaskName());
		} else {
			logger.error(getTaskName() + ": error in protein extraction");
			System.err.println(format.format(new Date()) + "\t" + getTaskName() + ": error in protein extraction");
			// Update progress on error
			updateProgress(0, 100);
			updateMessage(getTaskName() + ": Failed - Protein extraction error");
			Platform.runLater(() -> secondBar.setProgress(0.0)); // Stop the indeterminate progress
			updateMessage("Workflow failed!");
			return false;
		}

		// Step 4: Database Search
		logger.info("Step 4: Starting database search through sample-specific database...");
		updateProgress(40, 100);
		updateMessage(getTaskName() + ": Searching sample-specific database...");
		Platform.runLater(() -> secondBar.setProgress(ProgressBar.INDETERMINATE_PROGRESS));
		updateMessage("Preparing for database search...");

		logger.info("[{}] searching sample-specific database started", getTaskName());

		boolean search = this.separateSearch();

		if (search) {
			logger.info("Step 4: Database search completed successfully");
			System.out.println(
					format.format(new Date()) + "\t" + getTaskName() + ": searching sample-specific database finished");

			// Update progress after successful search
			updateProgress(65, 100);
			updateMessage(getTaskName() + ": Database search completed");
			Platform.runLater(() -> secondBar.setProgress(ProgressBar.INDETERMINATE_PROGRESS));
			updateMessage("Database search completed successfully");
		} else {
			logger.error(getTaskName() + ": error in database search");
			System.err.println(format.format(new Date()) + "\t" + getTaskName() + ": error in database search");

			return false;
		}

		// Step 5: Combine Results
		logger.info("Step 5: Starting combination of separate search results...");
		updateProgress(70, 100);
		updateMessage(getTaskName() + ": Combining search results...");
		Platform.runLater(() -> secondBar.setProgress(ProgressBar.INDETERMINATE_PROGRESS));
		updateMessage("Merging identification results...");

		System.out.println(
				format.format(new Date()) + "\t" + getTaskName() + ": combining separate search results started");

		boolean combine = this.combineResults();

		if (combine) {
			logger.info("Step 5: Combining separate search results completed successfully");
			System.out.println(
					format.format(new Date()) + "\t" + getTaskName() + ": combining separate search results finished");
		} else {
			logger.error(getTaskName() + ": error in combining separate search results");
			System.err.println(
					format.format(new Date()) + "\t" + getTaskName() + ": error in combining separate search results");

			return false;
		}

		// Step 6: Quantification
		logger.info("Step 6: Starting quantification...");
		updateProgress(75, 100);
		updateMessage(getTaskName() + ": Quantification...");
		Platform.runLater(() -> secondBar.setProgress(ProgressBar.INDETERMINATE_PROGRESS));
		updateMessage("Starting quantification analysis...");

		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": quantification started");

		boolean quant = false;

		if (magPfindParams.getQuantitativeMode().equals(MetaConstants.labelFree)) {
			logger.info("  Using label-free quantification mode");
			updateProgress(80, 100);
			updateMessage(getTaskName() + ": Label-free quantification...");
			Platform.runLater(() -> secondBar.setProgress(ProgressBar.INDETERMINATE_PROGRESS));
			updateMessage("Performing label-free quantification...");
			quant = lfQuant();
			if (!quant) {
				logger.error(getTaskName() + ": error in label-free quantification");
				System.err.println(
						format.format(new Date()) + "\t" + getTaskName() + ": error in label-free quantification");

				// Update progress on error
				updateProgress(0, 100);
				updateMessage(getTaskName() + ": Failed - Label-free quantification error");
				Platform.runLater(() -> secondBar.setProgress(0.0)); // Stop the indeterminate progress
				updateMessage("Workflow failed!");
				return false;
			}
			logger.info("  Label-free quantification completed successfully");

			// Update progress after successful label-free quantification
			updateProgress(80, 100);
			updateMessage(getTaskName() + ": Label-free quantification completed");
			Platform.runLater(() -> secondBar.setProgress(ProgressBar.INDETERMINATE_PROGRESS));
			updateMessage("Label-free quantification completed successfully");
		} else if (magPfindParams.getQuantitativeMode().equals(MetaConstants.isobaricLabel)) {
			logger.info("  Using isobaric labeling quantification mode");
			updateProgress(80, 100);
			updateMessage(getTaskName() + ": Label-free quantification for isobaric mode...");
			Platform.runLater(() -> secondBar.setProgress(ProgressBar.INDETERMINATE_PROGRESS));
			updateMessage("Performing label-free quantification...");
			quant = lfQuant();
			if (!quant) {
				logger.error(getTaskName() + ": error in label-free quantification for isobaric mode");
				System.err.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": error in label-free quantification for isobaric mode");

				return false;
			}
			logger.info("  Label-free quantification for isobaric mode completed successfully");

			updateProgress(80, 100);
			updateMessage(getTaskName() + ": Isobaric quantification...");
			Platform.runLater(() -> secondBar.setProgress(ProgressBar.INDETERMINATE_PROGRESS));
			updateMessage("Performing isobaric quantification...");

			logger.info("  Starting isobaric quantification...");
			System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": isobaric quantification started");

			quant = isobaricQuant();

			if (!quant) {
				logger.error(getTaskName() + ": error in isobaric quantification");
				System.err.println(
						format.format(new Date()) + "\t" + getTaskName() + ": error in isobaric quantification");

				return false;
			}
			logger.info("  Isobaric quantification completed successfully");
			System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": isobaric quantification finished");
		}

		// Step 7: Export Report
		logger.info("Step 7: Starting report generation...");
		updateProgress(85, 100);
		updateMessage(getTaskName() + ": Generating reports...");
		Platform.runLater(() -> secondBar.setProgress(ProgressBar.INDETERMINATE_PROGRESS));
		updateMessage("Preparing final reports...");

		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": report generation started");

		boolean report = this.exportReport();

		if (report) {
			logger.info("Step 7: Report generation completed successfully");
			System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": report generation finished");
		} else {
			logger.error(getTaskName() + ": error in report generation");
			System.err.println(format.format(new Date()) + "\t" + getTaskName() + ": error in report generation");

			return false;
		}

		// Update progress for completion
		updateProgress(100, 100);
		updateMessage(getTaskName() + ": Completed successfully!");
		Platform.runLater(() -> secondBar.setProgress(0.0)); // Stop the indeterminate progress
		updateMessage("Workflow completed!");
		// Stop the indeterminate progress

		logger.info(getTaskName() + ": finished - all workflow steps completed successfully");
		System.out.println(format.format(new Date()) + "\t" + getTaskName()
				+ ": finished - all workflow steps completed successfully");
		return true;
	}

	private boolean separateSearchHap() {
		logger.info("Starting separate search for high abundant database");

		// Get thread configuration
		int threadCount = magPfindParams.getThreadCount();

		File hapFastaFile = this.magDb.getHapFastaFile(magPfindParams.getGenomeSelectorMode());
		if (hapFastaFile == null || !hapFastaFile.exists()) {
			logger.error("HAP FASTA file is null or does not exist for MAG database: {}", this.magDb.getCatalogueID());
			return false;
		}

		// Check if we need to perform new searches
		boolean needNewSearch = false;
		int filesNeedingSearch = 0;
		for (int i = 0; i < this.pf2Files.length; i++) {
			File hapFile = new File(separateResultFolders[i], "hap");

			File resultFile = new File(hapFile, "pFind-Filtered.spectra");
			if (!resultFile.exists() || resultFile.length() == 0) {
				needNewSearch = true;
				filesNeedingSearch++;
			}
		}

		if (needNewSearch) {
			logger.info("Found {} files requiring new search out of {} total files", filesNeedingSearch,
					this.pf2Files.length);

			int threadPoolCount = maxTaskCount > filesNeedingSearch ? filesNeedingSearch : maxTaskCount;
			logger.info("Thread configuration: total={}, per-task={}, pool-size={}", threadCount,
					threadCount / threadPoolCount + 1, threadPoolCount);

			// Validate pFind installation and configuration
			File pFindBinDir = pFindExe.getParentFile();
			logger.info("Using pFind for high abundant database search");

			// Check for required pFind files
			File searcherExe = new File(pFindBinDir, "Searcher.exe");
			File templateParamDir = new File(pFindBinDir, "template/param");
			File pFindCfgTemplate = new File(templateParamDir, "pFind.cfg");
			File pParseCfgTemplate = new File(templateParamDir, "pParse.cfg");
			File pQuantCfgTemplate = new File(templateParamDir, "pQuant.cfg");

			logger.info("Checking pFind directory structure:");
			logger.info("  Searcher.exe exists: {} ({})", searcherExe.exists(), searcherExe.getAbsolutePath());
			logger.info("  template/param directory exists: {} ({})", templateParamDir.exists(),
					templateParamDir.getAbsolutePath());
			logger.info("  pFind.cfg template exists: {} ({})", pFindCfgTemplate.exists(),
					pFindCfgTemplate.getAbsolutePath());
			logger.info("  pParse.cfg template exists: {} ({})", pParseCfgTemplate.exists(),
					pParseCfgTemplate.getAbsolutePath());
			logger.info("  pQuant.cfg template exists: {} ({})", pQuantCfgTemplate.exists(),
					pQuantCfgTemplate.exists());

			// Check if required template files exist
			if (!pFindCfgTemplate.exists()) {
				logger.error("Required pFind template file not found: {}", pFindCfgTemplate.getAbsolutePath());
				logger.error("Please ensure pFind is properly installed with template configuration files");
				return false;
			}

			if (!pParseCfgTemplate.exists()) {
				logger.error("Required pParse template file not found: {}", pParseCfgTemplate.getAbsolutePath());
				logger.error("Please ensure pFind is properly installed with template configuration files");
				return false;
			}

			if (!pQuantCfgTemplate.exists()) {
				logger.error("Required pQuant template file not found: {}", pQuantCfgTemplate.getAbsolutePath());
				logger.error("Please ensure pFind is properly installed with template configuration files");
				return false;
			}

			// Create PFindTask and perform searches
			pFindTask.setTotalThread(threadCount / threadPoolCount + 1);

			// Get required parameters
			MetaOpenPar mosp = magPfindParams.getMop();
			String ms2Type = magPfindParams.getMs2ScanMode();

			for (int i = 0; i < this.pf2Files.length; i++) {
				File hapFile = new File(separateResultFolders[i], "hap");
				if (!hapFile.exists()) {
					hapFile.mkdirs();
				}

				File resultFile = new File(hapFile, "pFind-Filtered.spectra");
				if (!resultFile.exists() || resultFile.length() == 0) {
					File pf2 = new File(pf2Files[i]);
					if (pf2.exists()) {
						// Safety check for HAP FASTA file
						if (hapFastaFile == null || !hapFastaFile.exists()) {
							logger.error("HAP FASTA file is null or does not exist for MAG database: {}",
									this.magDb.getCatalogueID());
							return false;
						}

						logger.info("  Searching file {} with pFind", i + 1);
						pFindTask.searchPF2(new String[] { pf2.getAbsolutePath() }, hapFile.getAbsolutePath(),
								hapFastaFile.getAbsolutePath(), ms2Type, magPfindParams.isOpenSearch(), true, mosp.getOpenPsmFDR(),
								magPfindParams.getFixedModifications().toArray(new String[0]),
								magPfindParams.getVariableModifications().toArray(new String[0]),
								magPfindParams.getEnzyme(), magPfindParams.getMissCleavages(),
								magPfindParams.getDigestMode(), magPfindParams.getIsobaricTag());
					} else {
						logger.warn("PF2 file {} was not found, skipping", pf2Files[i]);
					}
				}
			}
			logger.info("Starting pFind execution with {} threads", threadPoolCount);
			pFindTask.run(threadPoolCount, separateResultFolders.length * 12);
		} else {
			logger.info("No new searches required, all results already exist");
		}

		FileFilter filter = new FileFilter() {
			@Override
			public boolean accept(File file) {
				return file.getName().equals("pFind-Filtered.spectra");
			}
		};

		// Use the generic WatchServiceLoop from the superclass
		logger.info("Starting WatchService monitoring for high abundant database search");
		int taskCount = this.pf2Files.length; // Calculate taskCount for the WatchService
		startWatchServiceLoop(separateResultFolders, "hap", filter, taskCount, 10, 38);

		logger.info("Verifying search results for {} folders", separateResultFolders.length);
		for (int i = 0; i < this.separateResultFolders.length; i++) {
			File hapFile = new File(separateResultFolders[i], "hap");
			logger.info("  Cleaning and verifying results for folder {}: {}", i + 1, hapFile.getAbsolutePath());
			PFindTask.cleanPFindResultFolder(hapFile, true);

			PFindResultReader reader = new PFindResultReader(hapFile, spectraDirFile);
			if (!reader.hasResult()) {
				System.out.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": pFind search result was not found from " + hapFile.getAbsolutePath());

				logger.error("pFind search result was not found from {}", hapFile.getAbsolutePath());
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

	protected boolean fastaExtractWithGenomeSelectorMode() {
		boolean useMarker4 = MagDbItem.MARKER4.equalsIgnoreCase(magPfindParams.getGenomeSelectorMode());
		logger.info("Starting FASTA extraction (genome selector mode={}) for {} result folders",
				useMarker4 ? MagDbItem.MARKER4 : MagDbItem.MARKER2, separateResultFolders.length);

		File currentFile = magDb.getCurrentFile();
		File markerFasta = null;
		File genomesTaxaTsv = null;
		if (useMarker4) {
			markerFasta = new File(currentFile, "Marker4.fasta");
			if (markerFasta == null || !markerFasta.exists()) {
				logger.error("Marker FASTA file was not found: {}",
						markerFasta != null ? markerFasta.getAbsolutePath() : "null");
				return false;
			}

			genomesTaxaTsv = new File(currentFile, "genomes-taxa.tsv");
			if (genomesTaxaTsv == null || !genomesTaxaTsv.exists()) {
				logger.error("Genomes taxa TSV was not found: {}",
						genomesTaxaTsv != null ? genomesTaxaTsv.getAbsolutePath() : "null");
				return false;
			}
		}

		MarkerGenomeSelector selector = new MarkerGenomeSelector();
		selector.setKneeSelectionMode(magPfindParams.getGenomeSelectionStrategy());
		int maximumSelected = 300;
		boolean enableAutoRep = useMarker4;
		String prefix = useMarker4 ? "marker_genome" : "genome";

		for (int i = 0; i < this.separateResultFolders.length; i++) {
			File hapFile = new File(separateResultFolders[i], "hap");
			File dbFile = new File(hapFile, separateResultFolders[i].getName() + ".fasta");

			if (dbFile.exists() && dbFile.length() > 0) {
				logger.info("Sample-specific database for {} already exists", separateResultFolders[i].getName());
				continue;
			}

			File spectraFile = new File(hapFile, "pFind-Filtered.spectra");
			if (!spectraFile.exists() || spectraFile.length() == 0) {
				logger.error("pFind search result was not found from {}", spectraFile.getAbsolutePath());
				return false;
			}

			PFindPsmReader psmReader = new PFindPsmReader(spectraFile);
			PFindPSM[] psms = psmReader.getPsms();

			File evidenceFile = new File(hapFile, prefix + "_evidence.tsv");
			File rankedFile = new File(hapFile, prefix + "_ranked.tsv");
			File selectedFile = new File(hapFile, prefix + "_selected.tsv");

			if (useMarker4) {
				selector.selectMarkerGenomesWithKnee(psms, markerFasta, genomesTaxaTsv, evidenceFile.getAbsolutePath(),
						rankedFile.getAbsolutePath(), selectedFile.getAbsolutePath(), maximumSelected, enableAutoRep);
			} else {
				selector.selectMarkerGenomesWithKnee(psms, evidenceFile.getAbsolutePath(), rankedFile.getAbsolutePath(),
						selectedFile.getAbsolutePath(), maximumSelected, enableAutoRep);
			}

			HashSet<String> genomeSet = new HashSet<>();
			if (selectedFile.exists()) {
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
					logger.error("Error reading selected genomes from {}", selectedFile.getAbsolutePath(), e);
					return false;
				}
			}

			logger.info("Marker selector ({}) for {} - Selected {} genomes",
					useMarker4 ? MagDbItem.MARKER4 : MagDbItem.MARKER2, separateResultFolders[i].getName(),
					genomeSet.size());

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

				if (magPfindParams.isAppendHostDb()) {
					BufferedReader hostReader = new BufferedReader(new FileReader(this.magPfindParams.getHostDB()));
					String linei;
					while ((linei = hostReader.readLine()) != null) {
						writer.println(linei);
					}
					hostReader.close();
				}
			} catch (IOException e) {
				logger.error("Error writing sample-specific database to {}", dbFile.getAbsolutePath(), e);
				return false;
			}
		}

		return true;
	}

	protected boolean separateSearch() {

		int taskCount = 0;
		for (int i = 0; i < this.separateResultFolders.length; i++) {

			PFindResultReader reader = new PFindResultReader(separateResultFolders[i], spectraDirFile);
			if (reader.hasResult()) {
				System.out.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": pFind search result already exists in " + separateResultFolders[i].getAbsolutePath());

				logger.info(getTaskName() + ": pFind search result already exists in "
						+ separateResultFolders[i].getAbsolutePath());
			} else {
				File hapFile = new File(separateResultFolders[i], "hap");
				File db = new File(hapFile, separateResultFolders[i].getName() + ".fasta");
				if (!db.exists()) {
					System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": sample-specific database"
							+ db + " was not found");

					logger.error(getTaskName() + ": sample-specific database" + db + " was not found");

					return false;
				}
				taskCount++;
			}
		}
		MetaOpenPar mosp = magPfindParams.getMop();
		String ms2Type = magPfindParams.getMs2ScanMode();
		if (taskCount > 0) {
			int totalThreadCount = magPfindParams.getThreadCount();
			int threadPoolCount = maxTaskCount > taskCount ? taskCount : maxTaskCount;
			int threadCount = totalThreadCount / threadPoolCount + 1;

			if (this.useDBReducer) {
				dbReducerTask = new DBReducerTask(dbReducerExe);
				dbReducerTask.setTotalThread(threadCount);

				for (int i = 0; i < this.mgfFiles.length; i++) {
					File dbReducerFile = new File(separateResultFolders[i], "DBReducer");
					if (!dbReducerFile.exists()) {
						dbReducerFile.mkdir();
					}
					File reducedDb = new File(dbReducerFile, "DBReducer.fasta");
					if (!reducedDb.exists()) {
						File mgf = new File(mgfFiles[i]);
						if (mgf.exists()) {
							File hapFile = new File(separateResultFolders[i], "hap");
							File db = new File(hapFile, separateResultFolders[i].getName() + ".fasta");
							if (db.length() > 0) {
								dbReducerTask.searchMgf(new String[] { mgf.getAbsolutePath() },
										dbReducerFile.getAbsolutePath(), db.getAbsolutePath(), ms2Type,
										mosp.getOpenPsmFDR(),
										magPfindParams.getFixedModifications().toArray(new String[0]),
										magPfindParams.getVariableModifications().toArray(new String[0]));
							} else {
								System.out.println(format.format(new Date()) + "\t" + getTaskName()
										+ ": no protein sequences were in " + db
										+ ", which means no peptides/proteins can be identifed from "
										+ separateResultFolders[i].getName());

								logger.error(getTaskName() + ": no protein sequences were in " + db
										+ ", which means no peptides/proteins can be identifed from "
										+ separateResultFolders[i].getName());
							}
						} else {

							System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": spectra file "
									+ mgf + " was not found");

							logger.error(getTaskName() + ": spectra file " + mgf + " was not found");

							return false;
						}
					}
				}
				dbReducerTask.run(threadPoolCount, separateResultFolders.length * 12);

				pFindTask = new PFindTask(pFindExe.getParentFile());
				pFindTask.setTotalThread(threadCount);

				for (int i = 0; i < this.pf2Files.length; i++) {
					File pf2 = new File(pf2Files[i]);
					File dbReducerFile = new File(separateResultFolders[i], "DBReducer");
					File db = new File(dbReducerFile, "DBReducer.fasta");
					File resultFile = new File(separateResultFolders[i], "pFind-Filtered.spectra");

					if (!resultFile.exists() || resultFile.length() == 0) {
						if (pf2.exists()) {

							if (db.exists() && db.length() > 0) {
								pFindTask.searchPF2(new String[] { pf2.getAbsolutePath() },
										separateResultFolders[i].getAbsolutePath(), db.getAbsolutePath(), ms2Type,
										magPfindParams.isOpenSearch(), true, mosp.getOpenPsmFDR(),
										magPfindParams.getFixedModifications().toArray(new String[0]),
										magPfindParams.getVariableModifications().toArray(new String[0]),
										magPfindParams.getEnzyme(), magPfindParams.getMissCleavages(),
										magPfindParams.getDigestMode(), magPfindParams.getIsobaricTag());

							} else {

								System.out.println(format.format(new Date()) + "\t" + getTaskName()
										+ ": no protein sequences were in " + db
										+ ", which means no peptides/proteins can be identifed from "
										+ separateResultFolders[i].getName());

								logger.error(getTaskName() + ": no protein sequences were in " + db
										+ ", which means no peptides/proteins can be identifed from "
										+ separateResultFolders[i].getName());
							}
						}
					}
				}

				pFindTask.run(threadPoolCount, this.separateResultFolders.length * 12);

			} else {
				PFindTask pFindTask = new PFindTask(pFindExe.getParentFile());
				pFindTask.setTotalThread(threadCount);

				for (int i = 0; i < this.pf2Files.length; i++) {
					File pf2 = new File(pf2Files[i]);
					File hapFile = new File(separateResultFolders[i], "hap");
					File db = new File(hapFile, separateResultFolders[i].getName() + ".fasta");
					File resultFile = new File(separateResultFolders[i], "pFind-Filtered.spectra");
					if (!resultFile.exists() || resultFile.length() == 0) {
						if (pf2.exists()) {
							if (db.exists() && db.length() > 0) {
								pFindTask.searchPF2(new String[] { pf2.getAbsolutePath() },
										separateResultFolders[i].getAbsolutePath(), db.getAbsolutePath(), ms2Type,
										magPfindParams.isOpenSearch(), true, mosp.getOpenPsmFDR(),
										magPfindParams.getFixedModifications().toArray(new String[0]),
										magPfindParams.getVariableModifications().toArray(new String[0]),
										magPfindParams.getEnzyme(),
										magPfindParams.getMissCleavages(),
										magPfindParams.getDigestMode(),
										magPfindParams.getIsobaricTag());
							} else {

								System.out.println(
										format.format(new Date()) + "\t" + getTaskName() + ": no protein sequences were in "
												+ db + ", which means no peptides/proteins can be identifed from "
												+ separateResultFolders[i].getName());

								logger.error(getTaskName() + ": no protein sequences were in " + db
										+ ", which means no peptides/proteins can be identifed from "
										+ separateResultFolders[i].getName());
							}
						}
					}
				}

				pFindTask.run(threadPoolCount, this.separateResultFolders.length * 12);
			}
		}

		FileFilter filter = new FileFilter() {
			@Override
			public boolean accept(File file) {
				return file.getName().equals("pFind-Filtered.spectra");
			}
		};

		// Use the generic WatchServiceLoop from the superclass
		logger.info("Starting WatchService monitoring for sample-specific database search");
		startWatchServiceLoop(separateResultFolders, "", filter, taskCount, 40, 68);

		for (int i = 0; i < this.separateResultFolders.length; i++) {
			PFindTask.cleanPFindResultFolder(separateResultFolders[i], true);

			File hapFile = new File(separateResultFolders[i], "hap");
			if (hapFile.exists()) {
				PFindTask.cleanPFindResultFolder(hapFile, true);
			}

			File dbReduceFile = new File(separateResultFolders[i], "DBReducer");
			if (dbReduceFile.exists()) {
				PFindTask.cleanPFindResultFolder(dbReduceFile, true);
			}
		}

		// Stop the watch service thread
		stopWatchServiceThread();

		return true;
	}

	protected boolean searchHostDatabase() {
		File hostDb = new File(magPfindParams.getHostDB());
		if (hostDb == null || !hostDb.exists()) {
			logger.error(getTaskName() + ": host database was not found: {}",
					hostDb != null ? hostDb.getAbsolutePath() : "null");
			return false;
		}

		File hostFolder = new File(metaPar.getDbSearchResultFile(), "host");
		if (!hostFolder.exists()) {
			hostFolder.mkdir();
		}
		File resultFile = new File(hostFolder, "pFind-Filtered.spectra");
		if (resultFile.exists() && resultFile.length() > 0) {
			logger.info(getTaskName() + ": host search results already exist");
			return true;
		}

		MetaOpenPar mosp = magPfindParams.getMop();
		String ms2Type = magPfindParams.getMs2ScanMode();
		int totalThreadCount = magPfindParams.getThreadCount();

		PFindTask pFindTask = new PFindTask(pFindExe.getParentFile());
		pFindTask.setTotalThread(totalThreadCount);

		String[] pf2Paths = new String[this.pf2Files.length];
		for (int i = 0; i < this.pf2Files.length; i++) {
			File pf2 = new File(pf2Files[i]);
			if (!pf2.exists()) {
				logger.error(getTaskName() + ": spectra file {} was not found", pf2.getAbsolutePath());
				return false;
			}
			pf2Paths[i] = pf2.getAbsolutePath();
		}
		pFindTask.searchPF2(pf2Paths, hostFolder.getAbsolutePath(), hostDb.getAbsolutePath(), ms2Type,
				magPfindParams.isOpenSearch(), true, mosp.getOpenPsmFDR(),
				magPfindParams.getFixedModifications().toArray(new String[0]),
				magPfindParams.getVariableModifications().toArray(new String[0]), magPfindParams.getEnzyme(),
				magPfindParams.getMissCleavages(), magPfindParams.getDigestMode(), magPfindParams.getIsobaricTag());

		pFindTask.run(1, 12);

		FileFilter filter = new FileFilter() {
			@Override
			public boolean accept(File file) {
				return file.getName().equals("pFind-Filtered.spectra");
			}
		};

		logger.info("Starting WatchService monitoring for host database search");
		startWatchServiceLoop(new File[] { hostFolder }, "", filter, 1, 68, 75);

		PFindTask.cleanPFindResultFolder(hostFolder, true);

		stopWatchServiceThread();
		return true;
	}

	protected boolean lfQuant() {

		updateMessage(getTaskName() + ": peptide and protein quantification");

		logger.info(getTaskName() + ": peptide and protein label-free quantification started");
		System.out.println(format.format(new Date()) + "\t" + getTaskName()
				+ ": peptide and protein label-free quantification started");

		String flashLFQ = mrp.getToolPath(ToolType.FLASHLFQ);
		flashLfqTask = new FlashLfqTask(flashLFQ, magPfindParams.getThreadCount());
		this.quan_pep_file = new File(combinedFile, "QuantifiedPeptides.tsv");
		
		if (!quan_pep_file.exists() || quan_pep_file.length() == 0) {

			File msms = new File(combinedFile, "psms.tsv");
			boolean lfqQuan = flashLfqTask.run(magPfindParams.getMetadata(), msms.getAbsolutePath(),
					rawDir.getAbsolutePath(), true, true);
			if (!lfqQuan) {

				logger.info(getTaskName() + ": error in peptide and protein label-free quantification");
				System.out.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": error in peptide and protein label-free quantification");

				return false;

			} else {
				if (!quan_pep_file.exists()) {

					logger.info(getTaskName() + ": error in peptide and protein label-free quantification");
					System.out.println(format.format(new Date()) + "\t" + getTaskName()
							+ ": error in peptide and protein label-free quantification");

					return false;
				}
			}
		} else {
			logger.info(getTaskName() + ": peptide quantification result has been found in " + quan_pep_file
					+ ", go to next step");
			System.out.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": peptide quantification result has been found in " + quan_pep_file + ", go to next step");
		}

		logger.info(getTaskName() + ": peptide and protein quantification label-free finished");
		System.out.println(format.format(new Date()) + "\t" + getTaskName()
				+ ": peptide and protein label-free quantification finished");

		return true;
	}

	protected boolean isobaricQuant() {

		updateMessage(getTaskName() + ": peptide and protein quantification");

		logger.info(getTaskName() + ": peptide and protein isobaric quantification started");
		System.out.println(format.format(new Date()) + "\t" + getTaskName()
				+ ": peptide and protein isobaric quantification started");

		ArrayList<PFindPSM> psmList = new ArrayList<PFindPSM>();
		for (int i = 0; i < this.separateResultFolders.length; i++) {
			PFindPsmReader psmReader = new PFindPsmReader(new File(this.separateResultFolders[i], "pFind.spectra"));
			PFindPSM[] psms = psmReader.getPsms();
			for (int j = 0; j < psms.length; j++) {
				psmList.add(psms[j]);
			}
		}
		PFindPSM[] psms = psmList.toArray(new PFindPSM[psmList.size()]);
		PFindIsobaricTask pfindIsobaricTask = new PFindIsobaricTask(psms, magPfindParams.getSpectraFile(), combinedFile,
				magPfindParams.getIsobaricTag(), magPfindParams.getIsobaric(), magPfindParams.getIsoCorFactor(),
				magPfindParams.isCombineLabel(), magPfindParams.getMetadata());

		try {
			pfindIsobaricTask.runOnlyPeptide();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			logger.error(getTaskName() + ": error in quantifying PSMs", e);
			System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": error in quantifying PSMs");
		}

		File quan_pep_file = pfindIsobaricTask.getQuanPepFile();

		if (quan_pep_file.exists()) {
			logger.info(getTaskName() + ": quantification results have been fount in "
					+ magPfindParams.getDbSearchResultFile());
			System.out.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": quantification results have been fount in " + magPfindParams.getDbSearchResultFile());
			return true;
		}

		if (!quan_pep_file.exists()) {
			return false;
		}

		logger.info(getTaskName() + ": peptide and protein isobaric quantification finished");
		System.out.println(format.format(new Date()) + "\t" + getTaskName()
				+ ": peptide and protein isobaric quantification finished");

		return true;
	}

	protected boolean combineResults() {
		boolean useMarker4 = MagDbItem.MARKER4.equalsIgnoreCase(magPfindParams.getGenomeSelectorMode());
		String modeLabel = useMarker4 ? "marker" : "combined";
		System.out.println(
				format.format(new Date()) + "\t" + getTaskName() + ": " + modeLabel + " PSM combination started");
		logger.info(getTaskName() + ": {} PSM combination started", modeLabel);

		File combinedPsms = new File(combinedFile, "psms.tsv");
		if (combinedPsms.exists() && combinedPsms.length() > 200) {
			logger.info(getTaskName() + ": PSMs already combined in {}", combinedPsms.getAbsolutePath());
			return true;
		}

		File rtFile = new File(magPfindParams.getSpectraFile(), "rt.txt");
		HashMap<String, Double> rtMap = new HashMap<String, Double>();
		try (BufferedReader reader = new BufferedReader(new FileReader(rtFile))) {
			String line = reader.readLine();
			while ((line = reader.readLine()) != null) {
				String[] cs = line.split("\t");
				rtMap.put(cs[0] + "_" + cs[1], Double.parseDouble(cs[2]));
			}
		} catch (IOException e) {
			System.out.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": error in reading the retention time information from " + rtFile);
			logger.error(getTaskName() + ": error in reading the retention time information from " + rtFile);
		}

		ArrayList<PFindPSM> psmList = new ArrayList<>();
		ArrayList<PFindPSM> decoyPsmList = new ArrayList<>();
		ArrayList<PFindPSM> hostPsmList = new ArrayList<>();
		for (int i = 0; i < this.separateResultFolders.length; i++) {
			File spectraFile = new File(this.separateResultFolders[i], "pFind-Filtered.spectra");
			if (!spectraFile.exists() || spectraFile.length() == 0) {
				logger.warn("PSM source file not found: {}", spectraFile.getAbsolutePath());
				continue;
			}

			PFindPsmReader psmReader = new PFindPsmReader(spectraFile);
			PFindPSM[] psms = psmReader.getPsms();

			for (int j = 0; j < psms.length; j++) {
				String sequence = psms[j].getBaseSequence();
				if (sequence.contains("X")) {
					continue;
				}

				ArrayList<String> revList = new ArrayList<String>();
				ArrayList<String> hostList = new ArrayList<String>();
				String[] proteins = psms[j].getProteins();
				for (int k = 0; k < proteins.length; k++) {
					if (proteins[k].startsWith("REV_")) {
						revList.add(proteins[k]);
					} else {
						if (!proteins[k].startsWith(this.magDb.getIdentifier())) {
							hostList.add(proteins[k]);
						}
					}
				}

				if (revList.size() > 0) {
					String[] revProteins = revList.toArray(String[]::new);
					psms[j].setProteins(revProteins);
					decoyPsmList.add(psms[j]);
				} else {
					if (hostList.size() > 0) {
						String[] hostProteins = hostList.toArray(String[]::new);
						psms[j].setProteins(hostProteins);
						hostPsmList.add(psms[j]);
					} else {
						psmList.add(psms[j]);
					}
				}
			}

			logger.info("{} PSMs were found from {}, now the total number of the PSMs is {}", psms.length,
					this.separateResultFolders[i].getName(), psmList.size());
		}

		if (psmList.isEmpty()) {
			logger.warn(getTaskName() + ": no PSMs found to combine");
			return false;
		}

		PFindPSM[] psms = psmList.toArray(PFindPSM[]::new);

		// Generate evidence.tsv, ranked.tsv, selected.tsv in combined/ for taxonomy
		// analysis
		// (Marker_selection column and genome_ranked.tsv). We do NOT filter PSMs by
		// selected genomes here.
		File evidenceTsv = new File(combinedFile, "evidence.tsv");
		File rankedTsv = new File(combinedFile, "ranked.tsv");
		File selectedTsv = new File(combinedFile, "selected.tsv");
		File currentFile = magDb.getCurrentFile();
		File markerFasta = null;
		File genomesTaxaTsv = null;
		if (useMarker4) {
			markerFasta = new File(currentFile, "Marker4.fasta");
			if (markerFasta == null || !markerFasta.exists()) {
				logger.error("Marker FASTA file was not found: {}",
						markerFasta != null ? markerFasta.getAbsolutePath() : "null");
				return false;
			}
			genomesTaxaTsv = new File(currentFile, "genomes-taxa.tsv");
			if (genomesTaxaTsv == null || !genomesTaxaTsv.exists()) {
				logger.error("Genomes taxa TSV was not found: {}",
						genomesTaxaTsv != null ? genomesTaxaTsv.getAbsolutePath() : "null");
				return false;
			}
		}
		// Generate evidence.tsv, ranked.tsv, selected.tsv using user's strategy.
		// The ranked.tsv contains knee columns (KneeStep, KneeKneedle, KneeCurvature,
		// KneePiecewise, KneeSlope)
		// which allow taxonomy analysis to compute core/mid/tail labels without running
		// selector twice.
		if (!evidenceTsv.exists() || !rankedTsv.exists() || !selectedTsv.exists()) {
			MarkerGenomeSelector selector = new MarkerGenomeSelector();
			selector.setKneeSelectionMode(MarkerGenomeSelector.KNEE_MODE_COMPREHENSIVE);
			boolean enableAutoRep = useMarker4;
			if (useMarker4) {
				selector.selectMarkerGenomesWithKnee(psms, markerFasta, genomesTaxaTsv, evidenceTsv.getAbsolutePath(),
						rankedTsv.getAbsolutePath(), selectedTsv.getAbsolutePath(), 500, enableAutoRep);
			} else {
				selector.selectMarkerGenomesWithKnee(psms, evidenceTsv.getAbsolutePath(), rankedTsv.getAbsolutePath(),
						selectedTsv.getAbsolutePath(), 500, enableAutoRep);
			}
		}
		if (rankedTsv.exists()) {
			logger.info("Marker selector ({}) generated evidence/ranked/selected genome lists in combined/",
					useMarker4 ? MagDbItem.MARKER4 : MagDbItem.MARKER2);
		} else {
			logger.warn("Marker selector did not produce ranked.tsv; Selection_label in Genome.tsv will be incomplete");
		}

		// Write ALL target PSMs (no filtering by selected genomes); quantification uses
		// complete list.
		int[] targetPSMCount = new int[2];
		OpenPsmTsvWriter msmsWriter = new OpenPsmTsvWriter(combinedPsms, OpenPsmTsvWriter.pFindFormat);
		for (PFindPSM psm : psms) {
			targetPSMCount[0]++;
			String scanKey = psm.getFileName() + "_" + psm.getScan();
			if (rtMap.containsKey(scanKey)) {
				psm.setRt(rtMap.get(scanKey));
				msmsWriter.addPFindPSM(psm);
				targetPSMCount[1]++;
			}
		}

		for (PFindPSM psm : hostPsmList) {
			targetPSMCount[0]++;
			String scanKey = psm.getFileName() + "_" + psm.getScan();
			if (rtMap.containsKey(scanKey)) {
				psm.setRt(rtMap.get(scanKey));
				msmsWriter.addPFindPSM(psm);
				targetPSMCount[1]++;
			}
		}

		logger.info("{} PSMs written to combined file; {} PSMs have assigned retention time", targetPSMCount[0],
				targetPSMCount[1]);

		for (PFindPSM psm : decoyPsmList) {
			String scanKey = psm.getFileName() + "_" + psm.getScan();
			if (rtMap.containsKey(scanKey)) {
				psm.setRt(rtMap.get(scanKey));
				msmsWriter.addPFindPSM(psm);
			}
			msmsWriter.addPFindPSM(psm);
		}

		msmsWriter.close();

		if (combinedPsms.exists() && combinedPsms.length() > 200) {
			System.out.println(
					format.format(new Date()) + "\t" + getTaskName() + ": " + modeLabel + " PSM combination finished");
			logger.info(getTaskName() + ": {} PSM combination finished", modeLabel);
			return true;
		}

		System.out.println(
				format.format(new Date()) + "\t" + getTaskName() + ": " + modeLabel + " PSM combination failed");
		logger.info(getTaskName() + ": {} PSM combination failed", modeLabel);
		return false;
	}

	protected boolean exportReport() {

		logger.info(getTaskName() + ": exporting peptide report started");
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
				reader.close();
			} catch (IOException e) {
				logger.error(getTaskName() + ": error in reading the experimental design from " + expDesignFile, e);
				System.err.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": error in reading the experimental design from " + expDesignFile);
			}
		} else {
			File[] files = magPfindParams.getDbSearchResultFile().listFiles();
			for (int i = 0; i < files.length; i++) {
				if (files[i].isDirectory() && !files[i].getName().equals("combined") && (new File(files[i], "pFind-Filtered.spectra")).exists()) {
					fileIdMap.put(files[i].getName(), fileIdMap.size());
					fileList.add(files[i].getName());
				}
			}
		}

		// Read PSMs
		FlashLfqPsm[] psms = null;
		try {
			FlashLfqPsmReader flashLfqPsmReader = new FlashLfqPsmReader(new File(combinedFile, "psms.tsv"));
			psms = flashLfqPsmReader.getPSMs();
		} catch (IOException e) {
			logger.error(getTaskName() + ": error in reading the PSMs from " + (new File(combinedFile, "psms.tsv")), e);
			System.err.println(format.format(new Date()) + "\t" + getTaskName() + ": error in reading the PSMs from "
					+ (new File(combinedFile, "psms.tsv")));
			return false;
		}

		// Read quantified peptides
		FlashLfqQuanPepReader quanPepReader;
		String quanMode = magPfindParams.getQuantitativeMode();
		if (quanMode.equals(MetaConstants.labelFree)) {
			quanPepReader = new FlashLfqQuanPepReader(this.quan_pep_file, fileList);
		} else if (quanMode.equals(MetaConstants.isobaricLabel)) {
			IsobaricTag isoTag = magPfindParams.getIsobaricTag();
			quanPepReader = new FlashLfqQuanPepReader(this.quan_pep_file, fileList, isoTag.getReportTitle());
		} else {
			logger.error(getTaskName() + ": unknown quantitative type");
			System.err.println(format.format(new Date()) + "\t" + getTaskName() + ": unknown quantitative type");
			return false;
		}

		FlashLfqQuanPeptide[] quanPeps = quanPepReader.getMetaPeptides();

		// Build intensity title array
		String[] intensityTitle;
		if (quanMode.equals(MetaConstants.labelFree)) {
			intensityTitle = fileList.toArray(String[]::new);
		} else if (quanMode.equals(MetaConstants.isobaricLabel)) {
			IsobaricTag isoTag = magPfindParams.getIsobaricTag();
			String[] repoTitle = isoTag.getReportTitle();
			intensityTitle = new String[fileList.size() * repoTitle.length];
			for (int i = 0; i < fileList.size(); i++) {
				for (int j = 0; j < repoTitle.length; j++) {
					intensityTitle[i * repoTitle.length + j] = fileList.get(i) + "_" + repoTitle[j];
				}
			}
		} else {
			logger.error(getTaskName() + ": unknown quantitative type");
			return false;
		}
		
		logger.info(getTaskName() + ": " + quanPeps.length + " peptides and " + psms.length + " PSMs were obtained");

		// Create and execute the post-analysis task directly
		DDAMagPostAnalysisTask postAnalysisTask = new DDAMagPostAnalysisTask(metaPar, mrp, secondBar, quanPeps, psms,
				fileIdMap, fileList, quanMode, intensityTitle, getTaskName(), PFindResultReader.REV, this);
		
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

}
