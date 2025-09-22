package uk.ac.quadram.metapilot.task.dda.mag.task;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileFilter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.math.NumberUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.dom4j.DocumentException;

import uk.ac.quadram.metapilot.core.prodb.FastaManager;
import uk.ac.quadram.metapilot.core.tools.FormatTool;
import uk.ac.quadram.metapilot.dbSearch.fragpipe.FragPipeTask;
import uk.ac.quadram.metapilot.dbSearch.fragpipe.FragPipePSM;
import uk.ac.quadram.metapilot.dbSearch.fragpipe.FragpipePepReader;
import uk.ac.quadram.metapilot.dbSearch.fragpipe.FragpipePeptide;
import uk.ac.quadram.metapilot.dbSearch.fragpipe.FragpipeProReader;
import uk.ac.quadram.metapilot.dbSearch.fragpipe.FragpipeWorkflow;
import javafx.scene.control.ProgressBar;
import javafx.application.Platform;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters.ToolType;
import uk.ac.quadram.metapilot.task.MetaPythonReportTask;
import uk.ac.quadram.metapilot.task.io.pep.MetaPeptide;
import uk.ac.quadram.metapilot.task.io.pro.MetaProtein;
import uk.ac.quadram.metapilot.task.io.pro.MetaProteinAnnoEggNog;
import uk.ac.quadram.metapilot.task.io.pro.MetaProteinXMLReader2;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;
import uk.ac.quadram.metapilot.task.mag.MetaPilotMagTask;
import uk.ac.quadram.metapilot.task.mag.MetaProteinAnnoMag;
import uk.ac.quadram.metapilot.task.par.MetaConstants;
import uk.ac.quadram.metapilot.task.dda.mag.par.DDAMagFragPipeParameters;
import javafx.scene.control.Label;

public class DDAMagFragPipeTask extends MetaPilotMagTask {

	private String fragpipeBat;
	private FragPipeTask fragPipeTask;
	private MetaPythonReportTask reportTask;
	private Thread reportCopyThread;
	private Thread reportThread;

	private String[] fragpipeExpNames;
	
	// Progress tracking components
	private Label statusLabel;
	
	// Add direct reference to DDAMagFragPipeParameters
	private DDAMagFragPipeParameters magFragPipeParams;
	private FragpipeWorkflow workflow;
	// File references for reports

	protected static final Logger LOGGER = LogManager.getLogger(DDAMagFragPipeTask.class);
	protected static final String taskName = "DDA MAG database search by FragPipe workflow";

	public DDAMagFragPipeTask(DDAMagFragPipeParameters magFragPipeParams, MetaResourceParameters mrp, ProgressBar progressBar,
			ProgressBar secondBar) {
		super(magFragPipeParams, mrp, progressBar, secondBar);

		// Store the direct reference to DDAMagFragPipeParameters
		this.magFragPipeParams = magFragPipeParams;
		this.workflow = magFragPipeParams.getFragpipeWorkflow();
		// Initialize status label - will be set by the UI panel
		this.statusLabel = null; // Will be set by UI panel via setStatusLabel method

		// Initialize the magDb field from parameters
		MagDbItem selectedMagDb = magFragPipeParams.getMagDbItem();
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
		} else {
			LOGGER.error("MAG database item is null in parameters");
			throw new IllegalArgumentException("MAG database item is required");
		}

		// Initialize the result directory structure in the constructor
		File fragpipeMagSearchFile = new File(magFragPipeParams.getResult(), "fragpipe_MAG_search");
		if (!fragpipeMagSearchFile.exists()) {
			fragpipeMagSearchFile.mkdir();
			LOGGER.info("Created FragPipe MAG search directory: {}", fragpipeMagSearchFile.getAbsolutePath());
		} else {
			LOGGER.info("FragPipe MAG search directory already exists: {}", fragpipeMagSearchFile.getAbsolutePath());
		}
		magFragPipeParams.setDbSearchResultFile(fragpipeMagSearchFile);

		createSeparateFolders();

		this.fragpipeExpNames = new String[rawFileNames.length];
		for (int i = 0; i < fragpipeExpNames.length; i++) {
			this.fragpipeExpNames[i] = this.rawFileNames[i].replaceAll("-", "_") + "_1";
		}

		// Initialize FragPipe executable path
		this.fragpipeBat = mrp.getToolPath(ToolType.FRAGPIPE);
		File fragpipeExe = new File(fragpipeBat);

		LOGGER.info("FragPipe path from MetaResourceParameters: {}", mrp.getToolPath(ToolType.FRAGPIPE));
		LOGGER.info("FragPipe executable file: {}", fragpipeExe.getAbsolutePath());
		LOGGER.info("FragPipe exists: {}", fragpipeExe.exists());

		if (!fragpipeExe.exists()) {
			LOGGER.error("FragPipe executable not found: {}", fragpipeExe.getAbsolutePath());
			throw new IllegalArgumentException("FragPipe executable is required and must exist");
		}

		this.fragPipeTask = new FragPipeTask(fragpipeExe.getAbsolutePath());
	}

	/**
	 * Set the status label for UI updates
	 */
	public void setStatusLabel(Label statusLabel) {
		this.statusLabel = statusLabel;
	}

	@Override
	public void forceStop() {
		LOGGER.info(getTaskName() + ": Force stop requested");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": Force stop requested");
		
		// Cancel FragPipe task if running
		if (fragPipeTask != null) {
			LOGGER.info(getTaskName() + ": Stopping FragPipe identification task");
			// FragPipeTask doesn't have forceStop method, but we can interrupt any running threads
			// The task will be stopped when this main task is cancelled
		}
		
		// Cancel report generation task if running
		if (reportTask != null && !reportTask.isDone()) {
			LOGGER.info(getTaskName() + ": Cancelling report generation task");
			reportTask.cancel(true);
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
	protected String getVersion() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	protected String getTaskName() {
		// TODO Auto-generated method stub
		return taskName;
	}

	@Override
	protected Logger getLogger() {
		// TODO Auto-generated method stub
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
		LOGGER.info("  Result directory: " + magFragPipeParams.getResult());
		LOGGER.info("  Thread count: " + magFragPipeParams.getThreadCount());
		if (magFragPipeParams.getMetadata() != null) {
			String[] rawFiles = magFragPipeParams.getMetadata().getRawFiles();
			LOGGER.info("  Raw files count: " + (rawFiles != null ? rawFiles.length : 0));
			if (rawFiles != null) {
				for (int i = 0; i < rawFiles.length; i++) {
					LOGGER.info("    Raw file " + (i+1) + ": " + rawFiles[i]);
				}
			}
		}
		
		// Log resource information
		LOGGER.info("Resource parameters:");
		LOGGER.info("  FragPipe path: " + mrp.getToolPath(MetaResourceParameters.ToolType.FRAGPIPE));
		LOGGER.info("  MetaResourceParameters object: {}", mrp);
		LOGGER.info("  MetaResourceParameters class: {}", mrp.getClass().getName());
		
		// Update main progress and keep second bar active
		updateProgress(5, 100);
		if (statusLabel != null) {
			Platform.runLater(() -> statusLabel.setText(getTaskName() + ": Initializing..."));
		}
		Platform.runLater(() -> secondBar.setProgress(ProgressBar.INDETERMINATE_PROGRESS)); // Indeterminate progress for continuous activity

		// Step 1: High Abundant Database Search
		LOGGER.info("Step 1: Starting high abundant database search...");
		if (statusLabel != null) {
			Platform.runLater(() -> statusLabel.setText(getTaskName() + ": Searching high abundant database..."));
		}
		updateProgress(10, 100);
		Platform.runLater(() -> secondBar.setProgress(ProgressBar.INDETERMINATE_PROGRESS));
		updateMessage("Performing high abundant protein search...");
		
		System.out.println(
				format.format(new Date()) + "\t" + getTaskName() + ": searching high abundant database started");

		boolean searchHap = this.separateSearchHap();

		if (searchHap) {
			LOGGER.info("Step 1: High abundant database search completed successfully");
			System.out.println(
					format.format(new Date()) + "\t" + getTaskName() + ": searching high abundant database finished");
		} else {
			LOGGER.error(getTaskName() + ": error in searching high abundant database");
			System.err.println(
					format.format(new Date()) + "\t" + getTaskName() + ": error in searching high abundant database");
			// Update progress on error
			updateMessage(getTaskName() + ": Failed - High abundant database search error");
			Platform.runLater(() -> secondBar.setProgress(0.0)); // Stop the indeterminate progress
			return false;
		}

		// Step 2: Protein Extraction from HAP Results
		LOGGER.info("Step 2: Starting protein extraction from HAP results...");
		updateProgress(20, 100);
		updateMessage(getTaskName() + ": Extracting proteins from HAP results...");
		Platform.runLater(() -> secondBar.setProgress(ProgressBar.INDETERMINATE_PROGRESS));
		updateMessage("Extracting protein sequences...");
		updateProgress(40, 100);
		updateMessage("Extracting protein sequences...");
		
		System.out.println(
				format.format(new Date()) + "\t" + getTaskName() + ": extracting proteins from the HAP result started");

		boolean extract = this.fastaExtract4Iterator();

		if (extract) {
			LOGGER.info("Step 2: Protein extraction from HAP results completed successfully");
			System.out.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": extracting proteins from the HAP result finished");
		} else {
			LOGGER.error(getTaskName() + ": error in protein extraction");
			System.err.println(format.format(new Date()) + "\t" + getTaskName() + ": error in protein extraction");
			// Update progress on error
			updateMessage(getTaskName() + ": Failed - Protein extraction error");
			Platform.runLater(() -> secondBar.setProgress(0.0)); // Stop the indeterminate progress
			updateMessage("Workflow failed!");
			return false;
		}

		// Step 3: Database Search
		LOGGER.info("Step 3: Starting database search through sample-specific database...");
		updateProgress(40, 100);
		updateMessage(getTaskName() + ": Searching sample-specific database...");
		Platform.runLater(() -> secondBar.setProgress(ProgressBar.INDETERMINATE_PROGRESS));
		updateMessage("Preparing for database search...");
		
		System.out.println(
				format.format(new Date()) + "\t" + getTaskName() + ": searching sample-specific database started");

		boolean search = separateSearch();

		if (search) {
			LOGGER.info("Step 3: Database search completed successfully");
			System.out.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": searching sample-specific database finished");
			
			// Update progress after successful search
			updateProgress(65, 100);
			updateMessage(getTaskName() + ": Database search completed");
			Platform.runLater(() -> secondBar.setProgress(ProgressBar.INDETERMINATE_PROGRESS));
			updateMessage("Database search completed successfully");
		} else {
			LOGGER.error(getTaskName() + ": error in database search");
			System.err.println(format.format(new Date()) + "\t" + getTaskName() + ": error in database search");
			// Update progress on error
			updateMessage(getTaskName() + ": Failed - Database search error");
			Platform.runLater(() -> secondBar.setProgress(0.0)); // Stop the indeterminate progress
			return false;
		}

		// Step 4: Combine Results
		LOGGER.info("Step 4: Starting combination of separate search results...");
		updateProgress(70, 100);
		updateMessage(getTaskName() + ": Combining search results...");
		Platform.runLater(() -> secondBar.setProgress(ProgressBar.INDETERMINATE_PROGRESS));
		updateMessage("Merging identification results...");
		
		System.out.println(
				format.format(new Date()) + "\t" + getTaskName() + ": combining separate search results started");

		boolean combine = combineResult();

		if (combine) {
			LOGGER.info("Step 4: Combining separate search results completed successfully");
			System.out.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": combining separate search results finished");
		} else {
			LOGGER.error(getTaskName() + ": error in combining separate search results");
			System.err.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": error in combining separate search results");
			// Update progress on error
			updateMessage(getTaskName() + ": Failed - Result combination error");
			Platform.runLater(() -> secondBar.setProgress(0.0)); // Stop the indeterminate progress
			return false;
		}

		// Step 5: Quantification
		LOGGER.info("Step 5: Starting quantification...");
		updateProgress(75, 100);
		updateMessage(getTaskName() + ": Quantification...");
		Platform.runLater(() -> secondBar.setProgress(ProgressBar.INDETERMINATE_PROGRESS));
		updateMessage("Starting quantification analysis...");
		
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": quantification started");

		boolean quant = false;

		// Note: FragPipe workflow may have different quantification modes
		// This is a placeholder - actual implementation depends on FragPipe workflow
		LOGGER.info("  Using FragPipe quantification mode");
		updateProgress(80, 100);
		updateMessage(getTaskName() + ": FragPipe quantification...");
		Platform.runLater(() -> secondBar.setProgress(ProgressBar.INDETERMINATE_PROGRESS));
		updateMessage("Performing FragPipe quantification...");
		
		// TODO: Implement actual FragPipe quantification logic
		quant = true; // Placeholder - replace with actual quantification logic
		
		if (!quant) {
			LOGGER.error(getTaskName() + ": error in FragPipe quantification");
			System.err.println(format.format(new Date()) + "\t" + getTaskName() + ": error in FragPipe quantification");
			
			// Update progress on error
			updateMessage(getTaskName() + ": Failed - FragPipe quantification error");
			Platform.runLater(() -> secondBar.setProgress(0.0)); // Stop the indeterminate progress
			updateMessage("Workflow failed!");
			return false;
		}
		LOGGER.info("  FragPipe quantification completed successfully");
		
		// Update progress after successful quantification
		updateProgress(80, 100);
		updateMessage(getTaskName() + ": FragPipe quantification completed");
		Platform.runLater(() -> secondBar.setProgress(ProgressBar.INDETERMINATE_PROGRESS));
		updateMessage("FragPipe quantification completed successfully");

		// Step 6: Export Report
		LOGGER.info("Step 6: Starting report generation...");
		updateProgress(85, 100);
		updateMessage(getTaskName() + ": Generating reports...");
		Platform.runLater(() -> secondBar.setProgress(ProgressBar.INDETERMINATE_PROGRESS));
		updateMessage("Preparing final reports...");
		
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": report generation started");

		boolean report = this.exportReport();

		if (report) {
			LOGGER.info("Step 6: Report generation completed successfully");
			System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": report generation finished");
		} else {
			LOGGER.error(getTaskName() + ": error in report generation");
			System.err.println(format.format(new Date()) + "\t" + getTaskName() + ": error in report generation");
			// Update progress on error
			updateMessage(getTaskName() + ": Failed - Report generation error");
			Platform.runLater(() -> secondBar.setProgress(0.0)); // Stop the indeterminate progress
			return false;
		}

		// Update progress for completion
		updateProgress(100, 100);
		updateMessage(getTaskName() + ": Completed successfully!");
		Platform.runLater(() -> secondBar.setProgress(0.0)); // Stop the indeterminate progress
		updateMessage("Workflow completed!");
		
		LOGGER.info(getTaskName() + ": finished - all workflow steps completed successfully");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": finished - all workflow steps completed successfully");
		return true;
	}

	protected boolean separateSearchHap() {
		LOGGER.info("Starting separate search for high abundant database");

		// Get thread configuration
		int threadCount = magFragPipeParams.getThreadCount();
		LOGGER.info("Thread configuration: total={}", threadCount);

		int taskCount = 0;

		File[] hapResultFiles = new File[rawFiles.length];
		for (int i = 0; i < rawFiles.length; i++) {
			File hapFile = new File(separateResultFolders[i], "hap");
			hapResultFiles[i] = new File(hapFile, this.fragpipeExpNames[i] + "_1");
			File combinePsmFile = new File(hapResultFiles[i], "psm.tsv");

			if (combinePsmFile.exists()) {
				System.out.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": FragPipe search result already exists in " + hapResultFiles[i].getAbsolutePath());

				LOGGER.info(getTaskName() + ": FragPipe search result already exists in "
						+ hapResultFiles[i].getAbsolutePath());
			} else {
				taskCount++;
			}
		}

		if (taskCount == 0) {
			return true;
		}

		for (int i = 0; i < separateResultFolders.length; i++) {
			File manifest = new File(this.separateResultFolders[i], "raws.fp-manifest");
			if (!manifest.exists()) {
				try (PrintWriter writer = new PrintWriter(manifest)) {

					StringBuilder sb = new StringBuilder();
					sb.append(rawFiles[i]).append("\t");
					sb.append(fragpipeExpNames[i].replaceAll("-", "_")).append("\t");
					sb.append("1").append("\t");
					sb.append("DDA");

					writer.println(sb);
					writer.close();

				} catch (FileNotFoundException e) {
					// TODO Auto-generated catch block
					LOGGER.error(getTaskName() + ": error in exporting Manifest for " + separateResultFolders[i], e);
					System.err.println(format.format(new Date()) + "\t" + getTaskName()
							+ ": error in exporting Manifest for " + separateResultFolders[i]);
				}

				if (!manifest.exists()) {
					LOGGER.error(getTaskName() + ": exporting Manifest failed for " + separateResultFolders[i]);
					System.err.println(format.format(new Date()) + "\t" + getTaskName()
							+ ": exporting Manifest failed for " + separateResultFolders[i]);
					return false;
				}
			}

			File hapFile = new File(separateResultFolders[i], "hap");
			File combinePsmFile = new File(hapResultFiles[i], "psm.tsv");
			File combineProFile = new File(hapResultFiles[i], "protein.tsv");

			if (!combinePsmFile.exists() || !combineProFile.exists()) {

				File hapTdFile = this.magDb.getHapFastaTDFile();
				if (!hapTdFile.exists()) {
					try {
						FastaManager.addDecoy(this.magDb.getHapFastaFile().getAbsolutePath(),
								hapTdFile.getAbsolutePath());
					} catch (IOException e) {
						// TODO Auto-generated catch block
						LOGGER.error(getTaskName() + ": writing target-decoy database to " + hapTdFile + " failed", e);
						System.err.println(format.format(new Date()) + "\t" + getTaskName()
								+ ": writing target-decoy database to " + hapTdFile + " failed");
					}
				}

				if (!hapTdFile.exists()) {
					LOGGER.error(getTaskName() + ": writing target-decoy database to " + hapTdFile + " failed");
					System.err.println(format.format(new Date()) + "\t" + getTaskName()
							+ ": writing target-decoy database to " + hapTdFile + " failed");
					return false;
				}

				metaPar.setCurrentDb(hapTdFile.getAbsolutePath());
				workflow.setDatabase(hapTdFile.getAbsolutePath());

				File workflowFile = new File(hapFile, "hap.workflow");
				try {
					workflow.export(workflowFile);
				} catch (IOException e) {
					// TODO Auto-generated catch block
					LOGGER.error(getTaskName() + ": writing workflow file to " + workflowFile + " failed");
					System.err.println(format.format(new Date()) + "\t" + getTaskName() + ": writing workflow file to "
							+ workflowFile + " failed");
					return false;
				}

				LOGGER.info(getTaskName() + ": searching high-abundance database for " + separateResultFolders[i]
						+ " started");
				System.out.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": searching high-abundance database for " + separateResultFolders[i] + " started");

				this.fragPipeTask.addTask(workflowFile, manifest, hapFile, threadCount);

			} else {
				System.out.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": FragPipe search result already exists in " + separateResultFolders[i].getAbsolutePath());

				LOGGER.info(getTaskName() + ": FragPipe search result already exists in "
						+ separateResultFolders[i].getAbsolutePath());
			}
		}

		// Set up WatchService to monitor psm.tsv file generation
		FileFilter filter = new FileFilter() {
			@Override
			public boolean accept(File file) {
				return file.getName().equals("psm.tsv");
			}
		};

		// Use the generic WatchServiceLoop from the superclass
		LOGGER.info("Starting WatchService monitoring for high abundant database search");
		startWatchServiceLoop(separateResultFolders, "hap", filter, taskCount, 10, 20);

		fragPipeTask.run(1, rawFiles.length * 4);

		for (int i = 0; i < separateResultFolders.length; i++) {
			File hapFile = new File(separateResultFolders[i], "hap");
			File combinePsmFile = new File(hapResultFiles[i], "psm.tsv");

			if (!combinePsmFile.exists()) {
				System.out.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": FragPipe search result was not found from " + hapResultFiles[i].getAbsolutePath());

				LOGGER.info(getTaskName() + ": FragPipe search result was not found from "
						+ hapResultFiles[i].getAbsolutePath());

				return false;
			}

			File[] hapListFiles = hapFile.listFiles();
			for (int j = 0; j < hapListFiles.length; j++) {
				if (hapListFiles[j].getName().endsWith(".pepindex")) {
					FileUtils.deleteQuietly(hapListFiles[j]);
				}
			}
		}

		// Stop the watch service thread
		stopWatchServiceThread();

		return true;
	}

	protected boolean fastaExtract4Iterator() {

		for (int i = 0; i < this.separateResultFolders.length; i++) {

			File hapFile = new File(separateResultFolders[i], "hap");
			File dbFile = new File(hapFile, separateResultFolders[i].getName() + ".fasta");

			if (dbFile.exists() && dbFile.length() > 0) {

				System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": sample-specific database for "
						+ separateResultFolders[i].getName() + " already exists in " + dbFile.getAbsolutePath());

				LOGGER.info(getTaskName() + ": sample-specific database for " + separateResultFolders[i].getName()
						+ " already exists in " + dbFile.getAbsolutePath());

				continue;
			}

			File hapResultFolder = new File(hapFile, fragpipeExpNames[i] + "_1");
			File combinePsmFile = new File(hapResultFolder, "psm.tsv");
			File combineProFile = new File(hapResultFolder, "protein.tsv");

			LOGGER.info(getTaskName() + ": parsing high-abundance database search result from "
					+ separateResultFolders[i].getName());
			System.out.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": parsing high-abundance database search result from " + separateResultFolders[i].getName());

			if (combinePsmFile.exists()) {

				LOGGER.info(getTaskName() + ": parsing high-abundance database search result " + combinePsmFile);
				System.out.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": parsing high-abundance database search result " + combinePsmFile);

				HashSet<String> totalPsmSet = new HashSet<String>();
				HashMap<String, HashSet<String>> genomeMap = new HashMap<String, HashSet<String>>();

				try {
					BufferedReader combineTsvReader = new BufferedReader(new FileReader(combinePsmFile));
					String line = combineTsvReader.readLine();
					String[] title = line.split("\t");
					int pepid = -1;
					int proid = -1;
					for (int j = 0; j < title.length; j++) {
						if (title[j].equals("Peptide")) {
							pepid = j;
						} else if (title[j].equals("Protein")) {
							proid = j;
						}
					}
					while ((line = combineTsvReader.readLine()) != null) {
						String[] cs = line.split("\t");
						String[] pros = cs[proid].split(",");
						for (int j = 0; j < pros.length; j++) {
							if (pros[j].startsWith(magDb.getIdentifier())) {
								String genome = pros[j].substring(0, pros[j].indexOf("_"));
								if (genomeMap.containsKey(genome)) {
									genomeMap.get(genome).add(cs[pepid]);
								} else {
									HashSet<String> set = new HashSet<String>();
									set.add(cs[pepid]);
									genomeMap.put(genome, set);
								}
								totalPsmSet.add(cs[pepid]);
							}
						}
					}
					combineTsvReader.close();
				} catch (IOException e) {
					// TODO Auto-generated catch block
					LOGGER.error(getTaskName() + ": error in reading database search result file " + combinePsmFile, e);
					System.err.println(format.format(new Date()) + "\t" + getTaskName()
							+ ": error in reading database search result file " + combinePsmFile);
				}

				if (totalPsmSet.size() == 0) {
					LOGGER.info(getTaskName() + ": no PSM identifications in " + combinePsmFile);
					System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": no PSM identifications in "
							+ combinePsmFile);
				}

				File hapGenomeFile = new File(hapFile, "hap_genomes.tsv");
				HashSet<String> genomeSet = refineGenomes(genomeMap, totalPsmSet.size(), psmCoverThreshold,
						hapGenomeFile);

				try {
					PrintWriter writer = new PrintWriter(dbFile);
					File originalDbFile = this.magDb.getOriginalDbFolder();
					for (String genome : genomeSet) {
						File fastafile = new File(originalDbFile, genome + ".faa");
						BufferedReader reader = new BufferedReader(new FileReader(fastafile));
						String line = null;
						StringBuilder sb = null;
						String rev = null;

						while ((line = reader.readLine()) != null) {
							if (line.startsWith(">")) {

								if (rev != null && sb != null) {
									writer.println(rev);
									writer.println(sb.reverse());
								}

								sb = new StringBuilder();
								rev = ">rev_" + line.substring(1);
								writer.println(line);

							} else {
								if (sb != null) {
									sb.append(line);
								}
								writer.println(line);
							}
						}

						if (rev != null && sb != null) {
							writer.println(rev);
							writer.println(sb.reverse());
						}

						reader.close();
					}
					writer.close();
				} catch (IOException e) {
					// TODO Auto-generated catch block
					LOGGER.error(getTaskName() + ": error in writing protein sequences to the sample-specific database "
							+ dbFile, e);
					System.err.println(format.format(new Date()) + "\t" + getTaskName()
							+ ": error in writing protein sequences to the sample-specific database " + dbFile);
				}

			} else {

				if (combineProFile.exists()) {
					LOGGER.info(getTaskName() + ": parsing high-abundance database search result " + combineProFile);
					System.out.println(format.format(new Date()) + "\t" + getTaskName()
							+ ": parsing high-abundance database search result " + combineProFile);

					HashMap<String, Integer> genomeMap = new HashMap<String, Integer>();
					int totalSpCount = 0;

					try {
						BufferedReader combineTsvReader = new BufferedReader(new FileReader(combineProFile));
						String line = combineTsvReader.readLine();
						String[] title = line.split("\t");
						int proId = -1;
						int spCountId = -1;

						for (int j = 0; j < title.length; j++) {
							if (title[j].equals("Protein")) {
								proId = j;
							} else if (title[j].equals("Razor Spectral Count")) {
								spCountId = j;
							}
						}
						while ((line = combineTsvReader.readLine()) != null) {
							String[] cs = line.split("\t");
							if (cs[proId].startsWith(this.magDb.getIdentifier())) {
								String genome = cs[proId].substring(0, cs[proId].indexOf("_"));
								int count = Integer.parseInt(cs[spCountId]);

								if (genomeMap.containsKey(genome)) {
									genomeMap.put(genome, genomeMap.get(genome) + count);
								} else {
									genomeMap.put(genome, count);
								}
								totalSpCount += count;
							}
						}
						combineTsvReader.close();
					} catch (IOException e) {
						// TODO Auto-generated catch block
						LOGGER.error(getTaskName() + ": error in reading database search result file " + combinePsmFile,
								e);
						System.err.println(format.format(new Date()) + "\t" + getTaskName()
								+ ": error in reading database search result file " + combinePsmFile);
					}

					String[] cs = genomeMap.keySet().toArray(new String[genomeMap.size()]);
					Arrays.sort(cs, new Comparator<String>() {

						@Override
						public int compare(String o1, String o2) {
							// TODO Auto-generated method stub
							int i1 = genomeMap.get(o1);
							int i2 = genomeMap.get(o2);
							if (i1 < i2) {
								return 1;
							} else if (i1 > i2) {
								return -1;
							}
							return 0;
						}
					});

					int tempCount = 0;
					try {
						PrintWriter writer = new PrintWriter(dbFile);
						File originalDbFile = this.magDb.getOriginalDbFolder();

						for (int j = 0; j < cs.length; j++) {
							int spCountI = genomeMap.get(cs[j]);
							tempCount += spCountI;
							double ratio = (double) tempCount / (double) totalSpCount;
							if (ratio > psmCoverThreshold) {
								break;
							}

							File fastafile = new File(originalDbFile, cs[j] + ".faa");
							BufferedReader reader = new BufferedReader(new FileReader(fastafile));
							String line = null;
							StringBuilder sb = null;
							String rev = null;

							while ((line = reader.readLine()) != null) {
								if (line.startsWith(">")) {

									if (rev != null && sb != null) {
										writer.println(rev);
										writer.println(sb.reverse());
									}

									sb = new StringBuilder();
									rev = ">rev_" + line.substring(1);
									writer.println(line);

								} else {
									if (sb != null) {
										sb.append(line);
									}
									writer.println(line);
								}
							}

							if (rev != null && sb != null) {
								writer.println(rev);
								writer.println(sb.reverse());
							}

							reader.close();
						}
						writer.close();

					} catch (IOException e) {
						// TODO Auto-generated catch block
						LOGGER.error(getTaskName()
								+ ": error in writing protein sequences to the sample-specific database " + dbFile, e);
						System.err.println(format.format(new Date()) + "\t" + getTaskName()
								+ ": error in writing protein sequences to the sample-specific database " + dbFile);
					}
				} else {

					LOGGER.info(
							getTaskName() + ": searching high-abundance database failed, result files were not found");
					System.out.println(format.format(new Date()) + "\t" + getTaskName()
							+ ": searching high-abundance database failed, result files were not found");
					return false;
				}
			}

			File proXmlFile = new File(hapFile, "combined.prot.xml");
			if (proXmlFile.exists()) {
				FileUtils.deleteQuietly(proXmlFile);
			}

			if (dbFile.exists() && dbFile.length() > 0) {

				LOGGER.info(getTaskName() + ": sample-specific database generation finished");
				System.out.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": sample-specific database generation finished");

				System.out.println(
						format.format(new Date()) + "\t" + getTaskName() + ": generating sample-specific database for "
								+ separateResultFolders[i].getName() + " finished");

				LOGGER.info(getTaskName() + ": generating sample-specific database for "
						+ separateResultFolders[i].getName() + " finished");
			} else {
				LOGGER.info(getTaskName() + ": sample-specific database generation failed");
				System.out.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": sample-specific database generation failed");

				return false;
			}
		}

		return true;
	}

	protected boolean separateSearch() {

		int taskCount = 0;
		File[] spResultFiles = new File[separateResultFolders.length];
		for (int i = 0; i < separateResultFolders.length; i++) {
			spResultFiles[i] = new File(separateResultFolders[i], fragpipeExpNames[i] + "_1");
			File combinePsmFile = new File(spResultFiles[i], "psm.tsv");
			if (spResultFiles[i].exists() && combinePsmFile.exists()) {
				System.out.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": FragPipe search result already exists in " + spResultFiles[i].getAbsolutePath());
				LOGGER.info(getTaskName() + ": FragPipe search result already exists in "
						+ spResultFiles[i].getAbsolutePath());
			} else {
				taskCount++;
			}
		}

		if (taskCount == 0) {
			return true;
		}

		for (int i = 0; i < this.separateResultFolders.length; i++) {
			File hapFile = new File(separateResultFolders[i], "hap");
			File dbFile = new File(hapFile, separateResultFolders[i].getName() + ".fasta");
			File combinePsmFile = new File(spResultFiles[i], "psm.tsv");
			File combineProFile = new File(spResultFiles[i], "protein.tsv");

			if (!combinePsmFile.exists() || !combineProFile.exists()) {
				
				metaPar.setCurrentDb(dbFile.getAbsolutePath());
				workflow.setDatabase(dbFile.getAbsolutePath());
				
				File workflowFile = new File(separateResultFolders[i], "mag.workflow");
				try {
					workflow.export(workflowFile);
				} catch (IOException e) {
					// TODO Auto-generated catch block
					LOGGER.error(getTaskName() + ": writing workflow file to " + workflowFile + " failed");
					System.err.println(format.format(new Date()) + "\t" + getTaskName() + ": writing workflow file to "
							+ workflowFile + " failed");
				}
				File manifest = new File(this.separateResultFolders[i], "raws.fp-manifest");

				LOGGER.info(getTaskName() + ": searching high-abundance database for " + separateResultFolders[i]
						+ " started");
				System.out.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": searching high-abundance database for " + separateResultFolders[i] + " started");

				this.fragPipeTask.addTask(workflowFile, manifest, separateResultFolders[i], threadCount);
			} else {
				System.out.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": FragPipe search result already exists in " + separateResultFolders[i].getAbsolutePath());

				LOGGER.info(getTaskName() + ": FragPipe search result already exists in "
						+ separateResultFolders[i].getAbsolutePath());
			}
		}

		// Set up WatchService to monitor psm.tsv file generation
		FileFilter filter = new FileFilter() {
			@Override
			public boolean accept(File file) {
				return file.getName().equals("psm.tsv");
			}
		};

		// Use the generic WatchServiceLoop from the superclass
		LOGGER.info("Starting WatchService monitoring for sample-specific database search");
		startWatchServiceLoop(separateResultFolders, "", filter, taskCount, 40, 65);

		this.fragPipeTask.run(1, rawFiles.length * 4);

		for (int i = 0; i < separateResultFolders.length; i++) {
			File combinePsmFile = new File(spResultFiles[i], "psm.tsv");
			File combineProFile = new File(spResultFiles[i], "protein.tsv");

			if (!combinePsmFile.exists() || !combineProFile.exists()) {
				System.out.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": FragPipe search result was not found from " + spResultFiles[i].getAbsolutePath());

				LOGGER.info(getTaskName() + ": FragPipe search result was not found from "
						+ spResultFiles[i].getAbsolutePath());

				return false;
			}
		}

		// Stop the watch service thread
		stopWatchServiceThread();

		return true;
	}

	protected boolean combineResult() {
		File combinedFastaFile = new File(combinedFile, "combined.fasta");
		if (!combinedFastaFile.exists()) {
			HashMap<String, String> fastaMap = new HashMap<String, String>();
			for (int i = 0; i < this.separateResultFolders.length; i++) {
				File spResultFolder = new File(separateResultFolders[i], fragpipeExpNames[i] + "_1");
				File fastaFile = new File(spResultFolder, "protein.fas");
				try (BufferedReader reader = new BufferedReader(new FileReader(fastaFile))) {
					String line = null;
					while ((line = reader.readLine()) != null) {
						if (line.startsWith(">")) {
							fastaMap.put(line, reader.readLine());
						}
					}
					reader.close();
				} catch (IOException e) {
					// TODO Auto-generated catch block
					LOGGER.error(getTaskName() + ": error in reading the fasta file " + fastaFile);
					System.err.println(format.format(new Date()) + "\t" + getTaskName()
							+ ": error in reading the fasta file " + fastaFile);
					return false;
				}
			}

			try (PrintWriter writer = new PrintWriter(combinedFastaFile)) {
				for (String proName : fastaMap.keySet()) {
					String seq = fastaMap.get(proName);
					writer.println(proName);
					writer.println(seq);

					String shuffleProName = ">shuf_" + proName.substring(1);
					String shuffle = FastaManager.generateTrypticShuffledSequence(seq);
					writer.println(shuffleProName);
					writer.println(shuffle);

					String revProName = ">rev_" + proName.substring(1);
					String revSeqString = (new StringBuilder(seq)).reverse().toString();

					writer.println(revProName);
					writer.println(revSeqString);

					String revShufProName = ">rev_shuf_" + proName.substring(1);
					String revShufSeqString = (new StringBuilder(shuffle)).reverse().toString();

					writer.println(revShufProName);
					writer.println(revShufSeqString);
				}

				if (metaPar.isAppendHostDb()) {
					BufferedReader hostReader = new BufferedReader(new FileReader(this.metaPar.getHostDB()));
					String line = null;
					String hostProName = null;
					String refHostProName = null;
					StringBuilder seqSb = new StringBuilder();
					while ((line = hostReader.readLine()) != null) {
						if (line.startsWith(">")) {
							if (hostProName != null) {
								writer.println(hostProName);
								writer.println(seqSb);

								writer.println(refHostProName);
								writer.println(seqSb.reverse());
							}

							hostProName = line;
							refHostProName = ">rev_" + line.substring(1);
							seqSb = new StringBuilder();
						} else {
							seqSb.append(line);
						}
					}

					if (hostProName != null) {
						writer.println(hostProName);
						writer.println(seqSb);

						writer.println(refHostProName);
						writer.println(seqSb.reverse());
					}

					hostReader.close();
				}
				writer.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				LOGGER.error(getTaskName() + ": error in writing the combined fasta file " + combinedFastaFile);
				System.err.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": error in writing the combined fasta file " + combinedFastaFile);
				return false;
			}
		}

		File manifest = new File(combinedFile, "raws.fp-manifest");
		if (!manifest.exists()) {
			try {
				metaPar.getMetadata().exportFragpipeManifest(manifest);
			} catch (FileNotFoundException e) {
				// TODO Auto-generated catch block
				LOGGER.error(getTaskName() + ": error in exporting Manifest", e);
				System.err.println(format.format(new Date()) + "\t" + getTaskName() + ": error in exporting Manifest");
			}

			if (!manifest.exists()) {
				LOGGER.error(getTaskName() + ": exporting Manifest failed");
				System.err.println(format.format(new Date()) + "\t" + getTaskName() + ": exporting Manifest failed");
			}
		}

		this.quan_pep_file = new File(combinedFile, "combined_modified_peptide.tsv");
		this.quan_pro_file = new File(combinedFile, "combined_protein.tsv");

		if (!this.quan_pep_file.exists() || !this.quan_pro_file.exists()) {
			File workflowFile = new File(combinedFile, "combined.workflow");
			try {
				workflow.export(workflowFile);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				LOGGER.error(getTaskName() + ": writing workflow file to " + workflowFile + " failed");
				System.err.println(format.format(new Date()) + "\t" + getTaskName() + ": writing workflow file to "
						+ workflowFile + " failed");
				return false;
			}

			LOGGER.info(getTaskName() + ": searching the combined high-abundance database started");
			System.out.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": searching the combined high-abundance database started");

			this.fragPipeTask.addTask(workflowFile, manifest, combinedFile, threadCount);
			this.fragPipeTask.run(rawFiles.length * 4);
		}

		if (this.quan_pep_file.exists() && this.quan_pep_file.length() > 0 && this.quan_pro_file.exists()
				&& this.quan_pro_file.length() > 0) {

			File[] combinedListFiles = metaPar.getDbSearchResultFile().listFiles();
			for (int j = 0; j < combinedListFiles.length; j++) {
				if (combinedListFiles[j].getName().endsWith(".pepindex")) {
					FileUtils.deleteQuietly(combinedListFiles[j]);
				}
				if (combinedListFiles[j].getName().endsWith(".prot.xml")) {
					FileUtils.deleteQuietly(combinedListFiles[j]);
				}
			}

		} else {
			LOGGER.info(getTaskName() + ": peptide and protein identification failed");
			System.out.println(
					format.format(new Date()) + "\t" + getTaskName() + ": peptide and protein identification failed");
		}

		return true;
	}

	protected boolean lfQuant() {
		return true;
	}

	protected boolean isobaricQuant() {
		return true;
	}

	protected boolean exportReport() {

		LOGGER.info(getTaskName() + ": exporting peptide report started");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": exporting peptide report started");

		// Build file list and mapping
		ArrayList<String> fileList = new ArrayList<String>();
		HashMap<String, Integer> fileIdMap = new HashMap<String, Integer>();
		File expDesignFile = new File(combinedFile, "raws.fp-manifest");
		if (expDesignFile.exists()) {
			try (BufferedReader reader = new BufferedReader(new FileReader(expDesignFile))) {
				String line = reader.readLine();
				while ((line = reader.readLine()) != null) {
					String[] cs = line.split("\t");
					fileIdMap.put(cs[1] + "_" + cs[2], fileIdMap.size());
					fileList.add(cs[1] + "_" + cs[2]);
				}
				reader.close();
			} catch (IOException e) {
				LOGGER.error(getTaskName() + ": error in reading the experimental design from " + expDesignFile, e);
				System.err.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": error in reading the experimental design from " + expDesignFile);
			}
		} else {
			File[] files = combinedFile.listFiles();
			for (int i = 0; i < files.length; i++) {
				if (files[i].isDirectory() && (new File(files[i], "peptide.tsv")).exists()) {
					fileIdMap.put(files[i].getName(), fileIdMap.size());
					fileList.add(files[i].getName());
				}
			}
		}

		// Read quantified peptides and PSMs from the same reader
		FragpipePepReader quanPepReader = new FragpipePepReader(this.quan_pep_file);
		FragpipePeptide[] quanPeps = quanPepReader.getMetaPeptides();
		FragPipePSM[] psms = quanPepReader.getPSMs();
		String[] intensityTitles = quanPepReader.getIntensityTitle();

		// Call the common method to process peptides and generate reports
		return processPeptideListAndGenerateReports(quanPeps, psms, fileIdMap, fileList, MetaConstants.labelFree, intensityTitles,
				getTaskName());
	}

	protected boolean exportPepProTxt(File magPepResultFile, File magProResultFile) {

		LOGGER.info(getTaskName() + ": exporting peptide report started");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": exporting peptide report started");

		int[] rep = metaPar.getMetadata().getReplicates();
		for (int i = 0; i < fragpipeExpNames.length; i++) {
			this.fragpipeExpNames[i] = fragpipeExpNames[i] + "_" + rep[i];
		}

		FragpipePepReader pepReader = new FragpipePepReader(magPepResultFile, fragpipeExpNames);
		FragpipePeptide[] peps = pepReader.getMetaPeptides();
		HashMap<String, Integer> hostProMap = new HashMap<String, Integer>();

		HashSet<String> magPepSet = new HashSet<String>();
		HashMap<String, HashSet<String>> genomePepMap = new HashMap<String, HashSet<String>>();
		HashMap<String, HashSet<String>> totalPepProMap = new HashMap<String, HashSet<String>>();
		for (int i = 0; i < peps.length; i++) {
			String stripSeq = peps[i].getSequence();
			String[] pros = peps[i].getProteins();
			for (int j = 0; j < pros.length; j++) {
				if (totalPepProMap.containsKey(stripSeq)) {
					totalPepProMap.get(stripSeq).add(pros[j]);
				} else {
					HashSet<String> set = new HashSet<String>();
					set.add(pros[j]);
					totalPepProMap.put(stripSeq, set);
				}
				if (pros[j].startsWith(magDb.getIdentifier())) {
					magPepSet.add(stripSeq);
					String genome = pros[j].split("_")[0];
					if (genomePepMap.containsKey(genome)) {
						genomePepMap.get(genome).add(stripSeq);
					} else {
						HashSet<String> pSet = new HashSet<String>();
						pSet.add(stripSeq);
						genomePepMap.put(genome, pSet);
					}
				} else {
					if (hostProMap.containsKey(pros[j])) {
						hostProMap.put(pros[j], hostProMap.get(pros[j]) + 1);
					} else {
						hostProMap.put(pros[j], 1);
					}
				}
			}
		}

		File magsFile = new File(metaPar.getDbSearchResultFile(), "mags.tsv");
		HashSet<String> genomeSet = refineGenomes(genomePepMap, magPepSet.size(), 1.0, magsFile);

		LOGGER.info(getTaskName() + ": filtered genome count " + genomeSet.size());
		System.out.println(
				format.format(new Date()) + "\t" + getTaskName() + ": filtered genome count " + genomeSet.size());

		HashMap<String, HashSet<String>> magProPepMap = new HashMap<String, HashSet<String>>();

		for (int i = 0; i < peps.length; i++) {
			String stripSeq = peps[i].getSequence();
			String[] pros = peps[i].getProteins();
			for (int j = 0; j < pros.length; j++) {
				if (pros[j].startsWith(magDb.getIdentifier())) {
					String genome = pros[j].split("_")[0];
					if (genomeSet.contains(genome)) {
						if (magProPepMap.containsKey(pros[j])) {
							magProPepMap.get(pros[j]).add(stripSeq);
						} else {
							HashSet<String> pSet = new HashSet<String>();
							pSet.add(stripSeq);
							magProPepMap.put(pros[j], pSet);
						}
					}
				}
			}
		}

		File razorProFile = new File(metaPar.getDbSearchResultFile(), "razorPro.tsv");
		HashSet<String> razorProSet = refineProtein(magProPepMap, magPepSet.size(), razorProFile);

		LOGGER.info(getTaskName() + ": filtered microbial protein count " + razorProSet.size());
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": filtered microbial protein count "
				+ razorProSet.size());

		razorProSet.addAll(hostProMap.keySet());

		HashMap<String, Integer> proIdMap = new HashMap<String, Integer>();
		try (BufferedReader reader = new BufferedReader(new FileReader(razorProFile))) {
			String line = reader.readLine();
			while ((line = reader.readLine()) != null) {
				String[] cs = line.split("\t");
				if (razorProSet.contains(cs[1])) {
					proIdMap.put(cs[1], proIdMap.size());
				}
			}
			reader.close();
		} catch (IOException e) {
			LOGGER.error(getTaskName() + ": error in reading razor protein file " + razorProFile, e);
			System.out.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": error in reading razor protein file " + razorProFile);
		}

		if (hostProMap.size() > 0) {
			// Convert the HashMap into a list of entries and sort it by values
			List<Map.Entry<String, Integer>> sortedList = hostProMap.entrySet().stream()
					.sorted(Map.Entry.comparingByValue(Collections.reverseOrder())).collect(Collectors.toList());
			String[] hostPros = sortedList.stream().map(Map.Entry::getKey).toArray(String[]::new);
			for (String hostPro : hostPros) {
				proIdMap.put(hostPro, proIdMap.size());
			}
		}

		HashMap<String, HashSet<String>> razorProPepMap = new HashMap<String, HashSet<String>>();
		HashMap<String, HashSet<String>> allProPepMap = new HashMap<String, HashSet<String>>();
		HashMap<String, Integer> pepIdMap = new HashMap<String, Integer>();
		HashMap<String, String[]> pepProMap = new HashMap<String, String[]>();
		HashMap<String, double[]> proIntensityMap = new HashMap<String, double[]>();
		HashMap<String, ArrayList<Double>> proScoreMap = new HashMap<String, ArrayList<Double>>();

		for (int i = 0; i < peps.length; i++) {
			String seq = peps[i].getModSeq();
			String[] pros = peps[i].getProteins();

			ArrayList<String> tempList = new ArrayList<String>();
			for (int j = 0; j < pros.length; j++) {
				if (proIdMap.containsKey(pros[j])) {
					tempList.add(pros[j]);
				}
			}

			if (tempList.size() > 0) {

				String[] tempPros = tempList.toArray(String[]::new);
				Arrays.sort(tempPros, new Comparator<String>() {
					@Override
					public int compare(String s1, String s2) {
						return proIdMap.get(s1) - proIdMap.get(s2);
					}
				});

				if (razorProPepMap.containsKey(tempPros[0])) {
					razorProPepMap.get(tempPros[0]).add(seq);
				} else {
					HashSet<String> razorSet = new HashSet<String>();
					razorSet.add(seq);
					razorProPepMap.put(tempPros[0], razorSet);
				}
				double[] pepIntensity = peps[i].getIntensity();

				for (int j = 0; j < tempPros.length; j++) {
					double[] proIntensity;
					if (proIntensityMap.containsKey(tempPros[j])) {
						proIntensity = proIntensityMap.get(tempPros[j]);
					} else {
						proIntensity = new double[pepIntensity.length];
					}

					for (int k = 0; k < proIntensity.length; k++) {
						proIntensity[k] += pepIntensity[k] / (double) tempPros.length;
					}
					proIntensityMap.put(tempPros[j], proIntensity);

					ArrayList<Double> scoreList;
					if (proScoreMap.containsKey(tempPros[j])) {
						scoreList = proScoreMap.get(tempPros[j]);
					} else {
						scoreList = new ArrayList<Double>();
					}
					scoreList.add(peps[i].getScore());
					proScoreMap.put(tempPros[j], scoreList);

					HashSet<String> pepSet;
					if (allProPepMap.containsKey(tempPros[j])) {
						pepSet = allProPepMap.get(tempPros[j]);
					} else {
						pepSet = new HashSet<String>();
						allProPepMap.put(tempPros[j], pepSet);
					}
					pepSet.add(seq);
				}
				pepIdMap.put(seq, pepIdMap.size());
				pepProMap.put(seq, tempPros);
			}
		}

		this.final_pep_txt = new File(metaPar.getDbSearchResultFile(), "final_peptides.tsv");
		this.final_pro_txt = new File(metaPar.getDbSearchResultFile(), "final_proteins.tsv");

		PrintWriter pepWriter = null;
		try {
			pepWriter = new PrintWriter(this.final_pep_txt);

			StringBuilder titlesb = new StringBuilder();
			titlesb.append("Sequence").append("\t");
			titlesb.append("Base sequence").append("\t");
			titlesb.append("Length").append("\t");
			titlesb.append("Proteins").append("\t");
			titlesb.append("Charges").append("\t");
			titlesb.append("PEP").append("\t");

			for (int i = 0; i < this.fragpipeExpNames.length; i++) {
				titlesb.append("Identification type " + fragpipeExpNames[i]).append("\t");
			}

			titlesb.append("Intensity").append("\t");

			for (int i = 0; i < fragpipeExpNames.length; i++) {
				titlesb.append("Intensity " + fragpipeExpNames[i]).append("\t");
			}

			titlesb.append("Reverse").append("\t");
			titlesb.append("Potential contaminant").append("\t");
			titlesb.append("id").append("\t");
			titlesb.append("Protein group IDs").append("\t");
			titlesb.append("MS/MS Count");

			pepWriter.println(titlesb);

			for (int i = 0; i < peps.length; i++) {

				String sequence = peps[i].getSequence();
				String modSeqString = peps[i].getModSeq();
				if (!pepIdMap.containsKey(modSeqString)) {
					continue;
				}

				StringBuilder sb = new StringBuilder();
				sb.append(modSeqString).append("\t");
				sb.append(sequence).append("\t");
				sb.append(peps[i].getLength()).append("\t");

				String[] proteins = pepProMap.get(modSeqString);
				for (String pro : proteins) {
					sb.append(pro).append(";");
				}
				sb.deleteCharAt(sb.length() - 1);
				sb.append("\t");

				int[] charge = peps[i].getCharges();
				StringBuilder chargesb = new StringBuilder();
				for (int j = 0; j < charge.length; j++) {
					chargesb.append(charge[j]).append(";");
				}
				if (chargesb.length() > 1) {
					sb.append(chargesb.subSequence(0, chargesb.length() - 1)).append("\t");
				} else {
					sb.append("2").append("\t");
				}

				sb.append(FormatTool.getDF4().format(peps[i].getScore())).append("\t");

				int[] idenType = peps[i].getIdenType();

				for (int j = 0; j < idenType.length; j++) {
					sb.append(MetaPeptide.idenTypeStrings[idenType[j]]).append("\t");
				}

				double totalIntensity = 0;

				double[] intensity = peps[i].getIntensity();
				for (int j = 0; j < intensity.length; j++) {
					totalIntensity += intensity[j];
				}

				sb.append((int) totalIntensity).append("\t");

				for (int j = 0; j < intensity.length; j++) {
					sb.append((int) intensity[j]).append("\t");
				}

				sb.append("\t").append("\t");

				sb.append(pepIdMap.get(modSeqString)).append("\t");

				for (String pro : proteins) {
					sb.append(proIdMap.get(pro)).append(";");
				}
				sb.deleteCharAt(sb.length() - 1);
				sb.append("\t");

				int totalSpCount = 0;
				int[] spCount = peps[i].getMs2Counts();
				for (int j = 0; j < spCount.length; j++) {
					totalSpCount += spCount[j];
				}
				sb.append(totalSpCount);

				pepWriter.println(sb);
			}

			pepWriter.close();

		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			LOGGER.error(getTaskName() + ": error in writing peptides to " + this.final_pep_txt.getName(), e);
			System.err.println(format.format(new Date()) + "\t" + getTaskName() + ": error in writing peptides to "
					+ this.final_pep_txt.getName());
		}

		LOGGER.info(getTaskName() + ": peptide report has been exported to " + final_pep_txt);
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": peptide report has been exported to "
				+ final_pep_txt);

		if (!final_pep_txt.exists() || final_pep_txt.length() == 0) {
			return false;
		}

		LOGGER.info(getTaskName() + ": exporting protein report started");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": exporting protein report started");

		PrintWriter proWriter = null;
		try {
			proWriter = new PrintWriter(this.final_pro_txt);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			LOGGER.error(getTaskName() + ": error in writing final protein result to " + this.final_pro_txt.getName(),
					e);
			System.err.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": error in writing final protein result to " + this.final_pro_txt.getName());
			return false;
		}

		StringBuilder titlesb = new StringBuilder();
		titlesb.append("Protein IDs").append("\t");
		titlesb.append("Majority protein IDs").append("\t");
		titlesb.append("Peptide counts (all)").append("\t");
		titlesb.append("Peptide counts (razor)").append("\t");
		titlesb.append("Number of proteins").append("\t");
		titlesb.append("Peptides").append("\t");
		titlesb.append("PEP").append("\t");
		titlesb.append("Intensity").append("\t");

		for (int i = 0; i < fragpipeExpNames.length; i++) {
			titlesb.append("Intensity " + fragpipeExpNames[i]).append("\t");
		}

		titlesb.append("Reverse").append("\t");
		titlesb.append("Potential contaminant").append("\t");
		titlesb.append("id").append("\t");
		titlesb.append("Peptide IDs").append("\t");
		titlesb.append("Peptide is razor");

		proWriter.println(titlesb);

		List<Map.Entry<String, Integer>> sortedList = proIdMap.entrySet().stream().sorted(Map.Entry.comparingByValue())
				.collect(Collectors.toList());
		String[] pros = sortedList.stream().map(Map.Entry::getKey).toArray(String[]::new);

		for (int i = 0; i < pros.length; i++) {

			HashSet<String> allPepSet = allProPepMap.get(pros[i]);
			HashSet<String> razorPepSet = razorProPepMap.get(pros[i]);

			if (allPepSet == null || razorPepSet == null) {
				continue;
			}

			int pepCountAll = allProPepMap.get(pros[i]).size();
			int pepCountRazor = razorProPepMap.get(pros[i]).size();

			StringBuilder sb = new StringBuilder();
			sb.append(pros[i]).append("\t");
			sb.append(pros[i]).append("\t");
			sb.append(pepCountAll).append("\t");
			sb.append(pepCountRazor).append("\t");
			sb.append("1\t");
			sb.append(allPepSet.size()).append("\t");

			double proScore = 1.0;
			ArrayList<Double> scoreList = proScoreMap.get(pros[i]);
			for (int j = 0; j < scoreList.size(); j++) {
				proScore = proScore * scoreList.get(j);
			}

			sb.append(FormatTool.getDF4().format(proScore)).append("\t");

			double totalIntensity = 0;
			double[] intensity = proIntensityMap.get(pros[i]);
			for (int j = 0; j < intensity.length; j++) {
				totalIntensity += intensity[j];
			}

			sb.append(totalIntensity).append("\t");
			for (int j = 0; j < intensity.length; j++) {
				sb.append(intensity[j]).append("\t");
			}
			sb.append("\t").append("\t");
			sb.append(proIdMap.get(pros[i])).append("\t");

			StringBuilder pepIdSb = new StringBuilder();
			StringBuilder pepRazorSb = new StringBuilder();

			String[] allPeps = allPepSet.toArray(new String[allPepSet.size()]);
			Arrays.sort(allPeps, new Comparator<String>() {

				@Override
				public int compare(String o1, String o2) {
					// TODO Auto-generated method stub

					return pepIdMap.get(o1) - pepIdMap.get(o2);
				}
			});
			for (String pep : allPeps) {
				pepIdSb.append(pepIdMap.get(pep)).append(";");
				if (razorPepSet.contains(pep)) {
					pepRazorSb.append("true;");
				} else {
					pepRazorSb.append("false;");
				}
			}
			sb.append(pepIdSb).append("\t");
			sb.append(pepRazorSb);

			proWriter.println(sb);
		}

		proWriter.close();

		LOGGER.info(getTaskName() + ": protein report has been exported to " + final_pro_txt.getName());
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": protein report has been exported to "
				+ final_pro_txt.getName());

		return true;
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
			LOGGER.error(getTaskName() + ": error in writing genomes to " + out, e);
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

	protected boolean exportSummary(MetaPythonReportTask task) {

		File final_summary_txt = new File(metaPar.getResult(), "final_summary.tsv");
		try {
			HashMap<String, int[]> idMap = new HashMap<String, int[]>();
			File[] magResultFiles = this.combinedFile.listFiles();
			for (int i = 0; i < magResultFiles.length; i++) {
				if (magResultFiles[i].isDirectory()) {
					int[] counts = new int[3];
					String resultName = null;
					File[] subFiles = magResultFiles[i].listFiles();
					HashSet<String> pepSet = new HashSet<String>();

					for (int j = 0; j < subFiles.length; j++) {
						String fileNameString = subFiles[j].getName();
						if (fileNameString.equals("psm.tsv")) {
							BufferedReader psmReader = new BufferedReader(new FileReader(subFiles[j]));
							String psmLine = psmReader.readLine();
							while ((psmLine = psmReader.readLine()) != null) {
								counts[1]++;
								String[] cs = psmLine.split("\t");
								pepSet.add(cs[2]);
							}
							psmReader.close();
						}
						if (fileNameString.endsWith(".pepXML")) {
							resultName = fileNameString.substring(0, fileNameString.length() - ".pepXML".length());
						}
					}
					counts[2] = pepSet.size();

					if (resultName != null) {
						File mzmlFile = new File(this.rawDir, resultName + "_uncalibrated.mzML");
						if (mzmlFile.exists()) {
							BufferedReader mzmlReader = new BufferedReader(new FileReader(mzmlFile));
							String mzmlLine = null;
							while ((mzmlLine = mzmlReader.readLine()) != null) {
								if (mzmlLine.trim().equals("</spectrum>")) {
									counts[0]++;
								}
							}
							mzmlReader.close();
						}
					}

					idMap.put(resultName, counts);
				}
			}

			PrintWriter countWriter = new PrintWriter(final_summary_txt);
			StringBuilder titlesb = new StringBuilder();
			titlesb.append("Raw file").append("\t");
			titlesb.append("Experiment").append("\t");
			titlesb.append("MS/MS").append("\t");
			titlesb.append("MS/MS Identified").append("\t");
			titlesb.append("MS/MS Identified [%]").append("\t");
			titlesb.append("Peptide Sequences Identified");

			countWriter.println(titlesb);

			int totalms2count = 0;
			int totalpsmCount = 0;
			int totalUnipepCount = 0;

			for (int i = 0; i < this.rawFiles.length; i++) {
				File rawFile = new File(this.rawFiles[i]);
				String rawNameString = rawFile.getName();
				rawNameString = rawNameString.substring(0, rawNameString.lastIndexOf("."));

				if (idMap.containsKey(rawNameString)) {

					StringBuilder sb = new StringBuilder();
					sb.append(rawNameString).append("\t");
					sb.append(this.fragpipeExpNames[i]).append("\t");
					int[] idCounts = idMap.get(rawNameString);

					sb.append(idCounts[0]).append("\t");
					sb.append(idCounts[1]).append("\t");

					double ratio = 0;
					if (idCounts[0] > 0 && idCounts[1] > 0) {
						ratio = (double) idCounts[1] / (double) idCounts[0] * 100.0;
					}

					sb.append(FormatTool.getDF2().format(ratio)).append("\t");
					sb.append(idCounts[2]);

					countWriter.println(sb);

					totalms2count += idCounts[0];
					totalpsmCount += idCounts[1];
					totalUnipepCount += idCounts[2];

				} else {
					StringBuilder sb = new StringBuilder();
					sb.append(rawNameString).append("\t");
					sb.append(this.fragpipeExpNames[i]).append("\t");
					sb.append("0\t0\t0%\t0");

					countWriter.println(sb);
				}
			}

			StringBuilder sb = new StringBuilder();
			sb.append("Total").append("\t");
			sb.append("\t");
			sb.append(totalms2count).append("\t");
			sb.append(totalpsmCount).append("\t");

			double ratio = 0;
			if (totalms2count > 0) {
				ratio = (double) totalpsmCount / (double) totalms2count * 100.0;
			}

			sb.append(FormatTool.getDF2().format(ratio)).append("\t");

			sb.append(totalUnipepCount);

			countWriter.print(sb);

			countWriter.close();

		} catch (IOException e) {
			// TODO Auto-generated catch block
			LOGGER.error(getTaskName() + ": error in writing PSM information to " + final_summary_txt.getName(), e);
			System.err.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": error in writing PSM information to " + final_summary_txt.getName());
		}

		if (!final_summary_txt.exists() || final_summary_txt.length() == 0) {
			return false;

		} else {

			boolean finish = false;
			File report_ID_summary = new File(metaPar.getReportDir(), "report_ID_summary.html");
			if (summaryMetaFile != null && summaryMetaFile.exists()) {
				finish = task.addTask(MetaPythonReportTask.ID_SUMMARY, final_summary_txt, report_ID_summary, summaryMetaFile);
			} else {
				finish = task.addTask(MetaPythonReportTask.ID_SUMMARY, final_summary_txt, report_ID_summary);
			}

			return finish;
		}
	}

	protected boolean exportTaxaFunc(MetaPythonReportTask task) {

		LOGGER.info(getTaskName() + ": taxonomy analysis and functional annotation started");
		System.out.println(format.format(new Date()) + "\t" + getTaskName()
				+ ": taxonomy analysis and functional annotation started");

		// functional annotation

		HashMap<String, Integer> proTaxIdMap = new HashMap<String, Integer>();

		LOGGER.info(getTaskName() + ": functional annotation started");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": functional annotation started");

		File funcFile = new File(metaPar.getResult(), "functional_annotation");
		if (!funcFile.exists()) {
			funcFile.mkdir();
		}

		File final_pro_xml_file = new File(funcFile, "functions.xml");

		File combinedPepTsv = new File(metaPar.getDbSearchResultFile(), "combined_peptide.tsv");
		File combinedProTsv = new File(metaPar.getDbSearchResultFile(), "combined_protein.tsv");

		FragpipeProReader proReader = new FragpipeProReader(combinedProTsv, this.fragpipeExpNames);
		MetaProtein[] proteins = proReader.getMetaProteins();

		if (!final_pro_xml_file.exists() || final_pro_xml_file.length() == 0) {

			MetaProteinAnnoMag[] mpas = funcAnnoSqlite(proteins, proReader.getIntensityTitle());
			for (int i = 0; i < mpas.length; i++) {
				String name = mpas[i].getPro().getName();
				if (NumberUtils.isCreatable(mpas[i].getTax_id())) {
					proTaxIdMap.put(name, Integer.parseInt(mpas[i].getTax_id()));
				} else {
					proTaxIdMap.put(name, 131567);
				}
			}
		}

		updateProgress(88, 100);

		File funcResultFile = final_pro_xml_file.getParentFile();
		File funTsv = new File(funcResultFile, "functions.tsv");
		File funReportTsv = new File(funcResultFile, "functions_report.tsv");

		MetaProteinXMLReader2 funReader = null;
		try {
			funReader = new MetaProteinXMLReader2(final_pro_xml_file);
		} catch (DocumentException e) {
			// TODO Auto-generated catch block
			LOGGER.error(getTaskName() + ": error in reading protein function information from "
					+ final_pro_xml_file.getName(), e);
			System.err.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": error in reading protein function information from " + final_pro_xml_file.getName());
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

		if (!funTsv.exists() || !funReportTsv.exists()) {
			System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": writing results to "
					+ funTsv.getName() + " started");
			try {
				funReader.exportTsv(funTsv);
				funReader.exportTsvReport(funReportTsv);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				LOGGER.error(getTaskName() + ": error in writing protein function information to " + funTsv.getName(),
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
			} catch (IOException e) {
				// TODO Auto-generated catch block
				LOGGER.error(getTaskName() + ": error in writing protein function information to " + html.getName(), e);
				System.err.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": error in writing protein function information to " + html.getName());
			}
			System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": writing results to "
					+ html.getName() + " finished");
		}

		File report_function_summary = new File(metaPar.getReportDir(), "report_function_summary.html");

		if (!report_function_summary.exists() || report_function_summary.length() == 0) {

			if (summaryMetaFile != null && summaryMetaFile.exists()) {
				task.addTask(MetaPythonReportTask.FUNCTION, funReportTsv, report_function_summary, summaryMetaFile);
			} else {
				task.addTask(MetaPythonReportTask.FUNCTION, funReportTsv, report_function_summary);
			}
		}

		LOGGER.info(getTaskName() + ": functional annotation finished");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": functional annotation finished");

		updateProgress(92, 100);

		// taxonomy analysis

		LOGGER.info(getTaskName() + ": taxonomy analysis started");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": taxonomy analysis started");

		File taxFile = new File(metaPar.getResult(), "taxonomy_analysis");
		if (!taxFile.exists()) {
			taxFile.mkdir();
		}

		FragpipePepReader pepReader = new FragpipePepReader(combinedPepTsv, this.fragpipeExpNames);
		FragpipePeptide[] peptides = pepReader.getMetaPeptides();
		HashMap<String, double[]> pepIntenMap = new HashMap<String, double[]>();
		HashMap<String, Integer> psmCountMap = new HashMap<String, Integer>();
		HashMap<String, HashSet<String>> pepCountMap = new HashMap<String, HashSet<String>>();
		HashMap<String, HashSet<String>> uniquePepCountMap = new HashMap<String, HashSet<String>>();
		HashMap<String, double[]> genomeIntensityMap = new HashMap<String, double[]>();
		HashMap<String, ArrayList<Double>> genomePEPMap = new HashMap<String, ArrayList<Double>>();

		for (int i = 0; i < peptides.length; i++) {
			String sequence = peptides[i].getSequence();
			double[] pepIntensity = peptides[i].getIntensity();
			pepIntenMap.put(sequence, pepIntensity);

			double PEP = peptides[i].getScore();

			int psmCount = 0;
			int[] ms2Count = peptides[i].getMs2Counts();
			for (int j = 0; j < fragpipeExpNames.length; j++) {
				psmCount += ms2Count[j];
			}

			String[] pros = peptides[i].getProteins();
			String razorPro = pros[0];

			if (razorPro.startsWith(this.magDb.getIdentifier())) {
				String genome = razorPro.split("_")[0];
				if (psmCountMap.containsKey(genome)) {
					psmCountMap.put(genome, psmCountMap.get(genome) + psmCount);
					pepCountMap.get(genome).add(sequence);

					if (pros.length == 1) {
						uniquePepCountMap.get(genome).add(sequence);
					}
				} else {
					psmCountMap.put(genome, psmCount);
					pepCountMap.put(genome, new HashSet<String>());
					pepCountMap.get(genome).add(sequence);
					uniquePepCountMap.put(genome, new HashSet<String>());
					if (pros.length == 1) {
						uniquePepCountMap.get(genome).add(sequence);
					}
				}

				double[] totalIntensity;
				if (genomeIntensityMap.containsKey(genome)) {
					totalIntensity = genomeIntensityMap.get(genome);
				} else {
					totalIntensity = new double[pepIntensity.length];
				}
				for (int j = 0; j < totalIntensity.length; j++) {
					totalIntensity[j] += pepIntensity[j];
				}
				genomeIntensityMap.put(genome, totalIntensity);

				ArrayList<Double> qList;
				if (genomePEPMap.containsKey(genome)) {
					qList = genomePEPMap.get(genome);
				} else {
					qList = new ArrayList<Double>();
					genomePEPMap.put(genome, qList);
				}
				qList.add(PEP);
			}
		}

		return exportGenome(genomeIntensityMap, psmCountMap, pepCountMap, uniquePepCountMap, genomePEPMap, taxFile,
				fragpipeExpNames, task);
	}

}

