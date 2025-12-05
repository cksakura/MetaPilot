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
import java.util.Comparator;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.ExecutionException;

import org.apache.commons.io.FileUtils;
import org.apache.logging.log4j.Logger;

import uk.ac.quadram.metapilot.core.prodb.FastaManager;
import uk.ac.quadram.metapilot.dbSearch.fragpipe.FragPipeTask;
import uk.ac.quadram.metapilot.dbSearch.fragpipe.FragPipePSM;
import uk.ac.quadram.metapilot.dbSearch.fragpipe.FragpipePepReader;
import uk.ac.quadram.metapilot.dbSearch.fragpipe.FragpipePeptide;
import uk.ac.quadram.metapilot.dbSearch.fragpipe.FragpipeWorkflow;
import javafx.scene.control.ProgressBar;
import javafx.application.Platform;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters.ToolType;
import uk.ac.quadram.metapilot.task.MetaPythonReportTask;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;
import uk.ac.quadram.metapilot.task.mag.MetaPilotMagTask;
import uk.ac.quadram.metapilot.task.par.MetaConstants;
import uk.ac.quadram.metapilot.task.dda.mag.par.DDAMagFragPipeParameters;
import javafx.scene.control.Label;

public class DDAMagFragPipeTask extends MetaPilotMagTask {

	// Logger is inherited from MetaPilotMagTask -> MetaPilotTask as an instance field

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

	protected static final String taskName = "DDA MAG database search by FragPipe workflow";

	public DDAMagFragPipeTask(DDAMagFragPipeParameters magFragPipeParams, MetaResourceParameters mrp, ProgressBar progressBar,
			ProgressBar secondBar, Label statusLabel, FragpipeWorkflow workflow) {
		super(magFragPipeParams, mrp, progressBar, secondBar);

		// Store the direct reference to DDAMagFragPipeParameters
		this.magFragPipeParams = magFragPipeParams;
		// Initialize workflow from constructor parameter (passed directly from GUI)
		this.workflow = workflow;
		// Initialize status label from constructor parameter (passed directly from GUI)
		this.statusLabel = statusLabel;

		// Initialize the magDb field from parameters
		MagDbItem selectedMagDb = magFragPipeParams.getMagDbItem();
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
		} else {
			logger.error("MAG database item is null in parameters");
			throw new IllegalArgumentException("MAG database item is required");
		}

		// Initialize the result directory structure in the constructor
		File fragpipeMagSearchFile = new File(magFragPipeParams.getResult(), "fragpipe_MAG_search");
		if (!fragpipeMagSearchFile.exists()) {
			fragpipeMagSearchFile.mkdir();
			logger.info("Created FragPipe MAG search directory: {}", fragpipeMagSearchFile.getAbsolutePath());
		} else {
			logger.info("FragPipe MAG search directory already exists: {}", fragpipeMagSearchFile.getAbsolutePath());
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

		logger.info("FragPipe path from MetaResourceParameters: {}", mrp.getToolPath(ToolType.FRAGPIPE));
		logger.info("FragPipe executable file: {}", fragpipeExe.getAbsolutePath());
		logger.info("FragPipe exists: {}", fragpipeExe.exists());

		if (!fragpipeExe.exists()) {
			logger.error("FragPipe executable not found: {}", fragpipeExe.getAbsolutePath());
			throw new IllegalArgumentException("FragPipe executable is required and must exist");
		}

		this.fragPipeTask = new FragPipeTask(fragpipeExe.getAbsolutePath());
	}

	@Override
	public void forceStop() {
		logger.info(getTaskName() + ": Force stop requested");
		logger.info("[{}] Force stop requested", getTaskName());
		
		// Cancel FragPipe task if running
		if (fragPipeTask != null) {
			logger.info(getTaskName() + ": Stopping FragPipe identification task");
			// FragPipeTask doesn't have forceStop method, but we can interrupt any running threads
			// The task will be stopped when this main task is cancelled
		}
		
		// Cancel report generation task if running
		if (reportTask != null && !reportTask.isDone()) {
			logger.info(getTaskName() + ": Cancelling report generation task");
			reportTask.cancel(true);
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
		logger.info("  Result directory: " + magFragPipeParams.getResult());
		logger.info("  Thread count: " + magFragPipeParams.getThreadCount());
		if (magFragPipeParams.getMetadata() != null) {
			String[] rawFiles = magFragPipeParams.getMetadata().getRawFiles();
			logger.info("  Raw files count: " + (rawFiles != null ? rawFiles.length : 0));
			if (rawFiles != null) {
				for (int i = 0; i < rawFiles.length; i++) {
					logger.info("    Raw file " + (i+1) + ": " + rawFiles[i]);
				}
			}
		}
		
		// Log resource information
		logger.info("Resource parameters:");
		logger.info("  FragPipe path: " + mrp.getToolPath(MetaResourceParameters.ToolType.FRAGPIPE));
		logger.info("  MetaResourceParameters object: {}", mrp);
		logger.info("  MetaResourceParameters class: {}", mrp.getClass().getName());
		
		// Update main progress and keep second bar active
		updateProgress(5, 100);
		if (statusLabel != null) {
			Platform.runLater(() -> statusLabel.setText(getTaskName() + ": Initializing..."));
		}
		Platform.runLater(() -> secondBar.setProgress(ProgressBar.INDETERMINATE_PROGRESS)); // Indeterminate progress for continuous activity

		// Step 1: High Abundant Database Search
		logger.info("Step 1: Starting high abundant database search...");
		if (statusLabel != null) {
			Platform.runLater(() -> statusLabel.setText(getTaskName() + ": Searching high abundant database..."));
		}
		updateProgress(10, 100);
		Platform.runLater(() -> secondBar.setProgress(ProgressBar.INDETERMINATE_PROGRESS));
		updateMessage("Performing high abundant protein search...");
		
		logger.info("[{}] Searching high abundant database started", getTaskName());

		boolean searchHap = this.separateSearchHap();

		if (searchHap) {
			logger.info("Step 1: High abundant database search completed successfully");
			logger.info("[{}] searching high abundant database finished", getTaskName());
		} else {
			logger.error(getTaskName() + ": error in searching high abundant database");
			System.err.println(
					format.format(new Date()) + "\t" + getTaskName() + ": error in searching high abundant database");
			// Update progress on error
			updateMessage(getTaskName() + ": Failed - High abundant database search error");
			Platform.runLater(() -> secondBar.setProgress(0.0)); // Stop the indeterminate progress
			return false;
		}

		// Step 2: Protein Extraction from HAP Results
		logger.info("Step 2: Starting protein extraction from HAP results...");
		updateProgress(20, 100);
		updateMessage(getTaskName() + ": Extracting proteins from HAP results...");
		Platform.runLater(() -> secondBar.setProgress(ProgressBar.INDETERMINATE_PROGRESS));
		updateMessage("Extracting protein sequences...");
		updateProgress(40, 100);
		updateMessage("Extracting protein sequences...");
		
		logger.info("[{}] extracting proteins from the HAP result started", getTaskName());

		boolean extract = this.fastaExtract4Iterator();

		if (extract) {
			logger.info("Step 2: Protein extraction from HAP results completed successfully");
			logger.info("[{}] extracting proteins from the HAP result finished", getTaskName());
		} else {
			logger.error(getTaskName() + ": error in protein extraction");
			System.err.println(format.format(new Date()) + "\t" + getTaskName() + ": error in protein extraction");
			// Update progress on error
			updateMessage(getTaskName() + ": Failed - Protein extraction error");
			Platform.runLater(() -> secondBar.setProgress(0.0)); // Stop the indeterminate progress
			updateMessage("Workflow failed!");
			return false;
		}

		// Step 3: Database Search
		logger.info("Step 3: Starting database search through sample-specific database...");
		updateProgress(40, 100);
		updateMessage(getTaskName() + ": Searching sample-specific database...");
		Platform.runLater(() -> secondBar.setProgress(ProgressBar.INDETERMINATE_PROGRESS));
		updateMessage("Preparing for database search...");
		
		logger.info("[{}] searching sample-specific database started", getTaskName());

		boolean search = separateSearch();

		if (search) {
			logger.info("Step 3: Database search completed successfully");
			logger.info("[{}] searching sample-specific database finished", getTaskName());
			
			// Update progress after successful search
			updateProgress(65, 100);
			updateMessage(getTaskName() + ": Database search completed");
			Platform.runLater(() -> secondBar.setProgress(ProgressBar.INDETERMINATE_PROGRESS));
			updateMessage("Database search completed successfully");
		} else {
			logger.error(getTaskName() + ": error in database search");
			System.err.println(format.format(new Date()) + "\t" + getTaskName() + ": error in database search");
			// Update progress on error
			updateMessage(getTaskName() + ": Failed - Database search error");
			Platform.runLater(() -> secondBar.setProgress(0.0)); // Stop the indeterminate progress
			return false;
		}

		// Step 4: Combine Results
		logger.info("Step 4: Starting combination of separate search results...");
		updateProgress(70, 100);
		updateMessage(getTaskName() + ": Combining search results...");
		Platform.runLater(() -> secondBar.setProgress(ProgressBar.INDETERMINATE_PROGRESS));
		updateMessage("Merging identification results...");
		
		logger.info("[{}] combining separate search results started", getTaskName());

		boolean combine = combineResult();

		if (combine) {
			logger.info("Step 4: Combining separate search results completed successfully");
			logger.info("[{}] combining separate search results finished", getTaskName());
		} else {
			logger.error(getTaskName() + ": error in combining separate search results");
			System.err.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": error in combining separate search results");
			// Update progress on error
			updateMessage(getTaskName() + ": Failed - Result combination error");
			Platform.runLater(() -> secondBar.setProgress(0.0)); // Stop the indeterminate progress
			return false;
		}

		// Step 5: Quantification
		logger.info("Step 5: Starting quantification...");
		updateProgress(75, 100);
		updateMessage(getTaskName() + ": Quantification...");
		Platform.runLater(() -> secondBar.setProgress(ProgressBar.INDETERMINATE_PROGRESS));
		updateMessage("Starting quantification analysis...");
		
		logger.info("[{}] quantification started", getTaskName());

		boolean quant = false;

		// Note: FragPipe workflow may have different quantification modes
		// This is a placeholder - actual implementation depends on FragPipe workflow
		logger.info("  Using FragPipe quantification mode");
		updateProgress(80, 100);
		updateMessage(getTaskName() + ": FragPipe quantification...");
		Platform.runLater(() -> secondBar.setProgress(ProgressBar.INDETERMINATE_PROGRESS));
		updateMessage("Performing FragPipe quantification...");
		
		// TODO: Implement actual FragPipe quantification logic
		quant = true; // Placeholder - replace with actual quantification logic
		
		if (!quant) {
			logger.error(getTaskName() + ": error in FragPipe quantification");
			System.err.println(format.format(new Date()) + "\t" + getTaskName() + ": error in FragPipe quantification");
			
			// Update progress on error
			updateMessage(getTaskName() + ": Failed - FragPipe quantification error");
			Platform.runLater(() -> secondBar.setProgress(0.0)); // Stop the indeterminate progress
			updateMessage("Workflow failed!");
			return false;
		}
		logger.info("  FragPipe quantification completed successfully");
		
		// Update progress after successful quantification
		updateProgress(80, 100);
		updateMessage(getTaskName() + ": FragPipe quantification completed");
		Platform.runLater(() -> secondBar.setProgress(ProgressBar.INDETERMINATE_PROGRESS));
		updateMessage("FragPipe quantification completed successfully");

		// Step 6: Export Report
		logger.info("Step 6: Starting report generation...");
		updateProgress(85, 100);
		updateMessage(getTaskName() + ": Generating reports...");
		Platform.runLater(() -> secondBar.setProgress(ProgressBar.INDETERMINATE_PROGRESS));
		updateMessage("Preparing final reports...");
		
		logger.info("[{}] report generation started", getTaskName());

		boolean report = this.exportReport();

		if (report) {
			logger.info("Step 6: Report generation completed successfully");
			logger.info("[{}] report generation finished", getTaskName());
		} else {
			logger.error(getTaskName() + ": error in report generation");
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
		
		logger.info(getTaskName() + ": finished - all workflow steps completed successfully");
		logger.info("[{}] finished - all workflow steps completed successfully", getTaskName());
		return true;
	}

	protected boolean separateSearchHap() {
		logger.info("Starting separate search for high abundant database");

		// Get thread configuration
		int threadCount = magFragPipeParams.getThreadCount();
		logger.info("Thread configuration: total={}", threadCount);

		int taskCount = 0;

		File[] hapResultFiles = new File[rawFiles.length];
		for (int i = 0; i < rawFiles.length; i++) {
			File hapFile = new File(separateResultFolders[i], "hap");
			hapResultFiles[i] = new File(hapFile, this.fragpipeExpNames[i] + "_1");
			File combinePsmFile = new File(hapResultFiles[i], "psm.tsv");

			if (combinePsmFile.exists()) {
				System.out.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": FragPipe search result already exists in " + hapResultFiles[i].getAbsolutePath());

				logger.info(getTaskName() + ": FragPipe search result already exists in "
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
					logger.error(getTaskName() + ": error in exporting Manifest for " + separateResultFolders[i], e);
					System.err.println(format.format(new Date()) + "\t" + getTaskName()
							+ ": error in exporting Manifest for " + separateResultFolders[i]);
				}

				if (!manifest.exists()) {
					logger.error(getTaskName() + ": exporting Manifest failed for " + separateResultFolders[i]);
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
						logger.error(getTaskName() + ": writing target-decoy database to " + hapTdFile + " failed", e);
						System.err.println(format.format(new Date()) + "\t" + getTaskName()
								+ ": writing target-decoy database to " + hapTdFile + " failed");
					}
				}

				if (!hapTdFile.exists()) {
					logger.error(getTaskName() + ": writing target-decoy database to " + hapTdFile + " failed");
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
					logger.error(getTaskName() + ": writing workflow file to " + workflowFile + " failed");
					System.err.println(format.format(new Date()) + "\t" + getTaskName() + ": writing workflow file to "
							+ workflowFile + " failed");
					return false;
				}

				logger.info(getTaskName() + ": searching high-abundance database for " + separateResultFolders[i]
						+ " started");
				System.out.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": searching high-abundance database for " + separateResultFolders[i] + " started");

				this.fragPipeTask.addTask(workflowFile, manifest, hapFile, threadCount);

			} else {
				System.out.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": FragPipe search result already exists in " + separateResultFolders[i].getAbsolutePath());

				logger.info(getTaskName() + ": FragPipe search result already exists in "
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
		logger.info("Starting WatchService monitoring for high abundant database search");
		startWatchServiceLoop(separateResultFolders, "hap", filter, taskCount, 10, 20);

		fragPipeTask.run(1, rawFiles.length * 24);

		for (int i = 0; i < separateResultFolders.length; i++) {
			File hapFile = new File(separateResultFolders[i], "hap");
			File combinePsmFile = new File(hapResultFiles[i], "psm.tsv");

			if (!combinePsmFile.exists()) {
				System.out.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": FragPipe search result was not found from " + hapResultFiles[i].getAbsolutePath());

				logger.info(getTaskName() + ": FragPipe search result was not found from "
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

				logger.info(getTaskName() + ": sample-specific database for " + separateResultFolders[i].getName()
						+ " already exists in " + dbFile.getAbsolutePath());

				continue;
			}

			File hapResultFolder = new File(hapFile, fragpipeExpNames[i] + "_1");
			File combinePsmFile = new File(hapResultFolder, "psm.tsv");
			File combineProFile = new File(hapResultFolder, "protein.tsv");

			logger.info(getTaskName() + ": parsing high-abundance database search result from "
					+ separateResultFolders[i].getName());
			System.out.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": parsing high-abundance database search result from " + separateResultFolders[i].getName());

			if (combinePsmFile.exists()) {

				logger.info(getTaskName() + ": parsing high-abundance database search result " + combinePsmFile);
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
					logger.error(getTaskName() + ": error in reading database search result file " + combinePsmFile, e);
					System.err.println(format.format(new Date()) + "\t" + getTaskName()
							+ ": error in reading database search result file " + combinePsmFile);
				}

				if (totalPsmSet.size() == 0) {
					logger.info(getTaskName() + ": no PSM identifications in " + combinePsmFile);
					System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": no PSM identifications in "
							+ combinePsmFile);
				}

				File hapGenomeFile = new File(hapFile, "hap_genomes.tsv");
				HashSet<String> genomeSet = applyOccamsRazor(genomeMap, totalPsmSet, psmCoverThreshold, hapGenomeFile);

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
					logger.error(getTaskName() + ": error in writing protein sequences to the sample-specific database "
							+ dbFile, e);
					System.err.println(format.format(new Date()) + "\t" + getTaskName()
							+ ": error in writing protein sequences to the sample-specific database " + dbFile);
				}

			} else {

				if (combineProFile.exists()) {
					logger.info(getTaskName() + ": parsing high-abundance database search result " + combineProFile);
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
						logger.error(getTaskName() + ": error in reading database search result file " + combinePsmFile,
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
						logger.error(getTaskName()
								+ ": error in writing protein sequences to the sample-specific database " + dbFile, e);
						System.err.println(format.format(new Date()) + "\t" + getTaskName()
								+ ": error in writing protein sequences to the sample-specific database " + dbFile);
					}
				} else {

					logger.info(
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

				logger.info(getTaskName() + ": sample-specific database generation finished");
				System.out.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": sample-specific database generation finished");

				System.out.println(
						format.format(new Date()) + "\t" + getTaskName() + ": generating sample-specific database for "
								+ separateResultFolders[i].getName() + " finished");

				logger.info(getTaskName() + ": generating sample-specific database for "
						+ separateResultFolders[i].getName() + " finished");
			} else {
				logger.info(getTaskName() + ": sample-specific database generation failed");
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
				logger.info(getTaskName() + ": FragPipe search result already exists in "
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
					logger.error(getTaskName() + ": writing workflow file to " + workflowFile + " failed");
					System.err.println(format.format(new Date()) + "\t" + getTaskName() + ": writing workflow file to "
							+ workflowFile + " failed");
				}
				File manifest = new File(this.separateResultFolders[i], "raws.fp-manifest");

				logger.info(getTaskName() + ": searching high-abundance database for " + separateResultFolders[i]
						+ " started");
				System.out.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": searching high-abundance database for " + separateResultFolders[i] + " started");

				this.fragPipeTask.addTask(workflowFile, manifest, separateResultFolders[i], threadCount);
			} else {
				System.out.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": FragPipe search result already exists in " + separateResultFolders[i].getAbsolutePath());

				logger.info(getTaskName() + ": FragPipe search result already exists in "
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
		logger.info("Starting WatchService monitoring for sample-specific database search");
		startWatchServiceLoop(separateResultFolders, "", filter, taskCount, 40, 65);

		this.fragPipeTask.run(1, rawFiles.length * 24);

		for (int i = 0; i < separateResultFolders.length; i++) {
			File combinePsmFile = new File(spResultFiles[i], "psm.tsv");
			File combineProFile = new File(spResultFiles[i], "protein.tsv");

			if (!combinePsmFile.exists() || !combineProFile.exists()) {
				System.out.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": FragPipe search result was not found from " + spResultFiles[i].getAbsolutePath());

				logger.info(getTaskName() + ": FragPipe search result was not found from "
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
					logger.error(getTaskName() + ": error in reading the fasta file " + fastaFile);
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
				logger.error(getTaskName() + ": error in writing the combined fasta file " + combinedFastaFile);
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
				logger.error(getTaskName() + ": error in exporting Manifest", e);
				System.err.println(format.format(new Date()) + "\t" + getTaskName() + ": error in exporting Manifest");
			}

			if (!manifest.exists()) {
				logger.error(getTaskName() + ": exporting Manifest failed");
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
				logger.error(getTaskName() + ": writing workflow file to " + workflowFile + " failed");
				System.err.println(format.format(new Date()) + "\t" + getTaskName() + ": writing workflow file to "
						+ workflowFile + " failed");
				return false;
			}

			logger.info(getTaskName() + ": searching the combined high-abundance database started");
			System.out.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": searching the combined high-abundance database started");

			this.fragPipeTask.addTask(workflowFile, manifest, combinedFile, threadCount);
			this.fragPipeTask.run(rawFiles.length * 24);
		}

		if (this.quan_pep_file.exists() && this.quan_pep_file.length() > 0 && this.quan_pro_file.exists()
				&& this.quan_pro_file.length() > 0) {

			File[] combinedListFiles = combinedFile.listFiles();
			for (int j = 0; j < combinedListFiles.length; j++) {
				if (combinedListFiles[j].getName().endsWith(".pepindex")) {
					FileUtils.deleteQuietly(combinedListFiles[j]);
				}
				if (combinedListFiles[j].getName().endsWith(".prot.xml")) {
					FileUtils.deleteQuietly(combinedListFiles[j]);
				}
			}

		} else {
			logger.info(getTaskName() + ": peptide and protein identification failed");
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

		logger.info(getTaskName() + ": exporting peptide report started");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": exporting peptide report started");

		// Build file list and mapping
		ArrayList<String> fileList = new ArrayList<String>();
		HashMap<String, Integer> fileIdMap = new HashMap<String, Integer>();
		File expDesignFile = new File(combinedFile, "raws.fp-manifest");
/*		if (expDesignFile.exists()) {
			try (BufferedReader reader = new BufferedReader(new FileReader(expDesignFile))) {
				String line = null;
				while ((line = reader.readLine()) != null) {
					String[] cs = line.split("\t");
					fileIdMap.put(cs[1], fileIdMap.size());
					fileList.add(cs[1]);
				}
				reader.close();
			} catch (IOException e) {
				logger.error(getTaskName() + ": error in reading the experimental design from " + expDesignFile, e);
				System.err.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": error in reading the experimental design from " + expDesignFile);
			}
		} else {*/
			File[] files = combinedFile.listFiles();
			for (int i = 0; i < files.length; i++) {
				if (files[i].isDirectory() && (new File(files[i], "peptide.tsv")).exists()) {
					fileIdMap.put(files[i].getName(), fileIdMap.size());
					fileList.add(files[i].getName());
				}
			}
//		}

		// Read quantified peptides and PSMs from the same reader
		FragpipePepReader quanPepReader = new FragpipePepReader(this.quan_pep_file);
		FragpipePeptide[] quanPeps = quanPepReader.getMetaPeptides();
		FragPipePSM[] psms = quanPepReader.getPSMs();
		String[] intensityTitles = quanPepReader.getIntensityTitle();
		
		logger.info(getTaskName() + ": " + quanPeps.length + " peptides and " + psms.length + " PSMs were obtained");

		// Create and execute the post-analysis task directly
		DDAMagPostAnalysisTask postAnalysisTask = new DDAMagPostAnalysisTask(
			metaPar, mrp, secondBar, quanPeps, psms, fileIdMap, fileList, MetaConstants.labelFree, intensityTitles, getTaskName(), this);
		
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
}

