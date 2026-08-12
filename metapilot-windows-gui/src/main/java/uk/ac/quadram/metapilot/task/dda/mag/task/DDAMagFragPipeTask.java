package uk.ac.quadram.metapilot.task.dda.mag.task;

import java.io.File;
import java.io.IOException;

import org.apache.logging.log4j.Logger;

import uk.ac.quadram.metapilot.core.runner.MagFragPipeHeadlessRunner;
import uk.ac.quadram.metapilot.core.mag.Marker4ResourcePreparer;
import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;
import uk.ac.quadram.metapilot.dbSearch.fragpipe.FragpipeWorkflow;
import javafx.scene.control.ProgressBar;
import javafx.application.Platform;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters.ToolType;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;
import uk.ac.quadram.metapilot.task.mag.MetaPilotMagTask;
import uk.ac.quadram.metapilot.task.dda.mag.par.DDAMagFragPipeParameters;
import uk.ac.quadram.metapilot.task.mag.MagDbItemResolver;
import javafx.scene.control.Label;

public class DDAMagFragPipeTask extends MetaPilotMagTask {

	// Logger is inherited from MetaPilotMagTask -> MetaPilotTask as an instance
	// field

	private String fragpipeBat;

	// Progress tracking components
	private Label statusLabel;

	// Add direct reference to DDAMagFragPipeParameters
	private DDAMagFragPipeParameters magFragPipeParams;
	private FragpipeWorkflow workflow;
	// File references for reports

	protected static final String taskName = "DDA MAG database search by FragPipe workflow";

	/**
	 * Matches pFind's opt-in to the nested {@code search_results/<sample>/}
	 * layout (see {@code MagResultLayoutUtil}). Note this only affects the
	 * {@code separateResultFolders} field this GUI task's base class
	 * pre-creates in its constructor - the actual FragPipe search/resume work
	 * runs entirely inside {@code MagFragPipeHeadlessRunner} (metapilot-core),
	 * which independently resolves and migrates the same
	 * {@code search_results/} location. Flipping this flag just keeps this
	 * pre-creation step consistent with that and avoids a pointless
	 * create-flat-then-immediately-migrate step on every fresh run.
	 */
	@Override
	protected boolean useNestedSearchResultsLayout() {
		return true;
	}

	public DDAMagFragPipeTask(DDAMagFragPipeParameters magFragPipeParams, MetaResourceParameters mrp,
			ProgressBar progressBar, ProgressBar secondBar, Label statusLabel, FragpipeWorkflow workflow) {
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

			// Validate that the MAG database has the required files (use genome selector
			// mode)
			String genomeSelectorMode = magFragPipeParams.getGenomeSelectorMode();
			if (MagDbItem.MARKER4.equalsIgnoreCase(genomeSelectorMode)) {
				try {
					Marker4ResourcePreparer.ensureSelectionResources(selectedMagDb, new ProgressCallback.NoOp());
				} catch (IOException e) {
					logger.error("Failed to prepare Marker4 resources for selected MAG database", e);
					throw new IllegalArgumentException("MAG database Marker4 resources could not be prepared", e);
				}
			}
			File hapFastaFile = selectedMagDb.getHapFastaFile(genomeSelectorMode);
			if (hapFastaFile == null || !hapFastaFile.exists()) {
				logger.error("MAG database HAP FASTA file is null or does not exist: {}",
						hapFastaFile != null ? hapFastaFile.getAbsolutePath() : "null");
				logger.error("Repository file: {}", selectedMagDb.getRepositoryFile());
				logger.error("Current file: {}", selectedMagDb.getCurrentFile());
				throw new IllegalArgumentException("MAG database HAP FASTA file is required and must exist");
			}
			logger.info("MAG database HAP FASTA file validated: {} (genome selector mode: {})",
					hapFastaFile.getAbsolutePath(), genomeSelectorMode);
			if (MagDbItem.MARKER4.equalsIgnoreCase(genomeSelectorMode)) {
				File genomesTaxaTsv = new File(selectedMagDb.getCurrentFile(), "genomes-taxa.tsv");
				if (genomesTaxaTsv == null || !genomesTaxaTsv.exists()) {
					logger.error("Marker4 mode requires genomes-taxa.tsv: {}",
							genomesTaxaTsv != null ? genomesTaxaTsv.getAbsolutePath() : "null");
					throw new IllegalArgumentException("MAG database genomes-taxa.tsv is required for Marker4 mode");
				}
				logger.info("Marker4 genomes-taxa.tsv validated: {}", genomesTaxaTsv.getAbsolutePath());
			}

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

	}

	@Override
	public void forceStop() {
		logger.info(getTaskName() + ": Force stop requested");
		logger.info("[{}] Force stop requested", getTaskName());

		// Cancel this task itself
		this.cancel(true);

		logger.info(getTaskName() + ": Force stop completed");
		logger.info("[{}] Force stop completed", getTaskName());
	}

	@Override
	protected String getVersion() {
		return null;
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
		return callViaCoreRunner();
	}

	private Boolean callViaCoreRunner() throws Exception {
		updateMessage(getTaskName());
		updateProgress(0, 100);

		logger.info(getTaskName() + ": started");
		updateProgress(5, 100);
		if (statusLabel != null) {
			Platform.runLater(() -> statusLabel.setText(getTaskName() + ": Initializing..."));
		}
		Platform.runLater(() -> secondBar.setProgress(ProgressBar.INDETERMINATE_PROGRESS));

		File fragpipeExe = new File(fragpipeBat);
		if (!fragpipeExe.exists()) {
			logger.error("FragPipe executable not found: {}", fragpipeExe.getAbsolutePath());
			Platform.runLater(() -> secondBar.setProgress(0.0));
			return false;
		}

		// Search mode / custom template overrides live entirely on magFragPipeParams
		// (search_mode, open_search_marker_stage, custom_basic_workflow_path,
		// custom_combined_workflow_path) - the same fields the CLI/HPC path reads from
		// the parameter JSON. If the panel already supplied an explicit custom-path
		// override, honor it as-is and skip the GUI-built temp file entirely; otherwise
		// export the GUI-configured (enzyme/mods/ion-mobility) workflow and route it
		// through those same fields rather than passing an explicit workflowTemplate
		// argument, so both stages resolve identically to a CLI/HPC run with the same
		// parameter file.
		File workflowFile = null;
		boolean hasCustomBasicPath = magFragPipeParams.getCustomBasicWorkflowPath() != null
				&& !magFragPipeParams.getCustomBasicWorkflowPath().isBlank();
		if (!hasCustomBasicPath && workflow != null) {
			workflowFile = exportWorkflowTemplateForCore();
			magFragPipeParams.setCustomBasicWorkflowPath(workflowFile.getAbsolutePath());
			boolean hasCustomCombinedPath = magFragPipeParams.getCustomCombinedWorkflowPath() != null
					&& !magFragPipeParams.getCustomCombinedWorkflowPath().isBlank();
			if (!hasCustomCombinedPath && !magFragPipeParams.isOpenSearch()) {
				// Closed mode: reuse the same GUI-built file for the combined stage too,
				// matching the historical single-shared-template behavior exactly. Open
				// mode leaves this unset so the combined stage resolves the Open+LFQ-MBR
				// hybrid instead (the GUI-built file has no LFQ-MBR quant settings baked in).
				magFragPipeParams.setCustomCombinedWorkflowPath(workflowFile.getAbsolutePath());
			}
		}
		File funcDbFile = resolveGuiToolFile(ToolType.FUNCDB);
		try {
			MagFragPipeHeadlessRunner.Result result = MagFragPipeHeadlessRunner.run(
					magFragPipeParams,
					fragpipeExe,
					null,
					MagDbItemResolver.resolveConfigFile(null),
					funcDbFile,
					createGuiProgressCallback());

			this.combinedFile = result.combinedDir;
			this.quan_pep_file = result.combinedModifiedPeptideTsv;
			this.quan_pro_file = result.combinedProteinTsv;
			this.final_pep_txt = new File(result.resultRoot, "final_peptides.tsv");
			this.final_pro_txt = new File(result.resultRoot, "razor_proteins.tsv");

			updateProgress(80, 100);
			updateMessage(getTaskName() + ": Generating visual reports...");
			if (statusLabel != null) {
				Platform.runLater(() -> statusLabel.setText(getTaskName() + ": Generating visual reports..."));
			}
			Platform.runLater(() -> secondBar.setProgress(ProgressBar.INDETERMINATE_PROGRESS));

			boolean report = this.generateFinalReports(getTaskName());
			if (!report) {
				updateMessage(getTaskName() + ": Failed - Report generation error");
				Platform.runLater(() -> secondBar.setProgress(0.0));
				return false;
			}

			updateProgress(100, 100);
			updateMessage(getTaskName() + ": Analysis completed successfully!");
			if (statusLabel != null) {
				Platform.runLater(() -> statusLabel.setText(getTaskName() + ": Completed"));
			}
			Platform.runLater(() -> secondBar.setProgress(1.0));
			return true;
		} finally {
			if (workflowFile != null && workflowFile.exists()) {
				workflowFile.delete();
			}
		}
	}

	private ProgressCallback createGuiProgressCallback() {
		return new ProgressCallback() {
			@Override
			public void onProgress(double workDone, double max) {
				if (max > 0) {
					updateProgress(workDone, max);
				}
			}

			@Override
			public void onMessage(String message) {
				if (message != null && !message.isBlank()) {
					updateMessage(message);
					if (statusLabel != null) {
						Platform.runLater(() -> statusLabel.setText(message));
					}
				}
			}

			@Override
			public void onComplete(Object result) {
			}

			@Override
			public void onError(Exception error) {
				if (error != null) {
					logger.error(getTaskName() + ": core runner error", error);
				}
			}
		};
	}

	private File exportWorkflowTemplateForCore() throws IOException {
		File tmp = File.createTempFile("metapilot-fragpipe-", ".workflow");
		tmp.deleteOnExit();
		if (workflow != null) {
			workflow.export(tmp);
			return tmp;
		}
		throw new IOException("FragPipe workflow template is not available");
	}

	private File resolveGuiToolFile(ToolType toolType) {
		if (mrp == null || toolType == null) {
			return null;
		}
		String path = mrp.getToolPath(toolType);
		if (path == null || path.isBlank()) {
			return null;
		}
		return new File(path);
	}

}
