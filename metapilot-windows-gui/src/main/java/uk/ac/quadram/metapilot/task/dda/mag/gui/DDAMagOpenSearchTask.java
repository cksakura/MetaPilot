package uk.ac.quadram.metapilot.task.dda.mag.gui;

import java.io.File;
import java.util.List;
import java.util.Optional;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.atomic.AtomicBoolean;

import org.apache.logging.log4j.Logger;

import javafx.application.Platform;
import javafx.scene.control.Alert;
import javafx.scene.control.Alert.AlertType;
import javafx.scene.control.ButtonBar;
import javafx.scene.control.ButtonType;
import javafx.scene.control.Label;
import javafx.scene.control.ProgressBar;
import javafx.scene.control.TextArea;
import javafx.scene.layout.GridPane;
import javafx.scene.layout.Priority;

import uk.ac.quadram.metapilot.core.opensearch.OpenSearchEntrapment;
import uk.ac.quadram.metapilot.core.opensearch.OpenSearchLayout;
import uk.ac.quadram.metapilot.core.runner.MagOpenSearchHeadlessRunner;
import uk.ac.quadram.metapilot.core.runner.MagOpenSearchToolSetup;
import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.task.dda.mag.par.DDAMagOpenSearchParameters;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;
import uk.ac.quadram.metapilot.task.mag.MetaPilotMagTask;

/**
 * GUI task for the DDA MAG open-search workflow.
 * <p>
 * Extends {@link MetaPilotMagTask} like the other MAG tasks, rather than a bare
 * {@code javafx.concurrent.Task}, so it inherits the behaviour every MetaPilot
 * workflow is expected to have:
 * <ul>
 * <li><b>Per-run log file.</b> {@code MetaPilotTask}'s constructor calls
 * {@code configLog(taskLogDir())}, which creates {@code logging.log} - here in
 * {@code opensearch_MAG_search/} - with its own non-additive appender
 * and tears it down on completion. An earlier version of this class extended
 * {@code Task} directly and silently wrote no log file at all.</li>
 * <li><b>Progress bar wiring</b> — the {@code progressProperty()} listener that
 * drives the overall bar.</li>
 * <li><b>Shared MAG post-analysis</b> — {@code refineGenomes(...)},
 * {@code refineProteins(...)}, {@code funcAnnoSqlite(...)} and
 * {@code generateFinalReports(...)}, so this workflow's downstream steps are
 * the same code the other engines run.</li>
 * </ul>
 * The workflow logic itself stays in {@link MagOpenSearchHeadlessRunner} so the
 * GUI, CLI and HPC paths execute identically.
 */
public class DDAMagOpenSearchTask extends MetaPilotMagTask {

	protected static final String taskName = "DDA MAG open search (MSFragger standalone)";

	private final DDAMagOpenSearchParameters openSearchParams;
	private final MetaResourceParameters mrp;
	private final Label statusLabel;

	public DDAMagOpenSearchTask(DDAMagOpenSearchParameters params, MetaResourceParameters mrp,
			ProgressBar progressBar, ProgressBar secondBar, Label statusLabel) {
		super(params, mrp, progressBar, secondBar);
		this.openSearchParams = params;
		this.mrp = mrp;
		this.statusLabel = statusLabel;

		// Everything below runs after super(...), so the task logger and its
		// logging.log appender already exist. Anything that fails here must be
		// written through `logger` before it propagates, otherwise it lands only
		// in the application log and the run's own log file shows nothing -
		// which is exactly what happened with the dbSearchResultFile NPE below.
		try {
			MagDbItem selected = params.getMagDbItem();
			if (selected == null) {
				throw new IllegalArgumentException("MAG database item is required");
			}
			this.magDb = selected;
			this.metaPar.setMicroDb(
					selected.getCurrentFile() == null ? "" : selected.getCurrentFile().getAbsolutePath());

			// Must precede createSeparateFolders(): it dereferences
			// getDbSearchResultFile() unguarded. The other MAG tasks set their
			// engine search folder here for the same reason; this workflow sets
			// opensearch_MAG_search/, where its final tables and report go, as
			// MagOpenSearchHeadlessRunner.run() does again.
			File modeRoot = modeRoot(params);
			if (!modeRoot.exists() && !modeRoot.mkdirs()) {
				throw new IllegalStateException("Could not create result folder: " + modeRoot.getAbsolutePath());
			}
			params.setDbSearchResultFile(modeRoot);
			logger.info("{}: result folder {} (per-sample searches in {})", taskName, modeRoot.getAbsolutePath(),
					sampleFoldersRoot().getAbsolutePath());

			createSeparateFolders();

		} catch (RuntimeException e) {
			logger.error("{}: initialization failed - {}", taskName, e.getMessage(), e);
			throw e;
		}
	}

	/** Matches the other MAG engines' {@code search_results/<sample>/} layout. */
	@Override
	protected boolean useNestedSearchResultsLayout() {
		return true;
	}

	/**
	 * Per-sample searches live under {@code opensearch_MAG_search/search_results/}
	 * (or {@code opensearch_entrapment_MAG_search/search_results/} for evaluation),
	 * the same nested layout FragPipe MAG and pFind MAG use.
	 */
	@Override
	protected File sampleFoldersRoot() {
		return OpenSearchLayout.searchResultsRoot(searchRoot((DDAMagOpenSearchParameters) metaPar));
	}

	/** The run's log belongs in the engine folder, with the tables. */
	@Override
	protected File taskLogDir() {
		return modeRoot((DDAMagOpenSearchParameters) metaPar);
	}

	private static File searchRoot(DDAMagOpenSearchParameters params) {
		File result = new File(params.getResult());
		if (params.isEntrapmentEnabled()) {
			return new File(result, OpenSearchEntrapment.RESULT_FOLDER);
		}
		return OpenSearchLayout.engineSearchRoot(result);
	}

	private static File modeRoot(DDAMagOpenSearchParameters params) {
		return OpenSearchLayout.modeRoot(searchRoot(params), params.getSearchEngines());
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
	protected String getVersion() {
		return null;
	}

	@Override
	public void forceStop() {
		logger.info("[{}] Force stop requested", getTaskName());
		this.cancel(true);
	}

	/**
	 * GUI decision after stage-1 when some samples have no refined database
	 * (empty unmodified/rescored set). Mirrors DIA Deep per-file's partial-marker
	 * dialog: continue identifies only successful samples, but FlashLFQ still
	 * quantifies every raw file together.
	 */
	private boolean askContinueWithPartialIdentification(List<String> skippedSampleNames, int successfulSampleCount)
			throws InterruptedException {
		CountDownLatch latch = new CountDownLatch(1);
		AtomicBoolean continueRun = new AtomicBoolean(false);
		Platform.runLater(() -> {
			try {
				ButtonType continueButton = new ButtonType("Continue with successful samples",
						ButtonBar.ButtonData.OK_DONE);
				ButtonType stopButton = new ButtonType("Stop", ButtonBar.ButtonData.CANCEL_CLOSE);
				Alert alert = new Alert(AlertType.CONFIRMATION);
				alert.setTitle("Incomplete stage-1 results");
				alert.setHeaderText("Some samples have no stage-1 genome selection result");
				alert.setContentText(
						"Marker search finished, but " + skippedSampleNames.size()
								+ " sample(s) have too few unmodified targets for genome selection.\n\n"
								+ successfulSampleCount
								+ " sample(s) have complete stage-1 results.\n\n"
								+ "Stop now if you want to replace the failing files, or continue: stage-2 "
								+ "identification uses only the successful samples, and FlashLFQ still "
								+ "quantifies every raw file together (match-between-runs).");

				TextArea details = new TextArea(String.join("\n", skippedSampleNames));
				details.setEditable(false);
				details.setWrapText(true);
				details.setMaxWidth(Double.MAX_VALUE);
				details.setMaxHeight(Double.MAX_VALUE);
				GridPane.setVgrow(details, Priority.ALWAYS);
				GridPane.setHgrow(details, Priority.ALWAYS);
				GridPane expandable = new GridPane();
				expandable.setMaxWidth(Double.MAX_VALUE);
				expandable.add(details, 0, 0);
				alert.getDialogPane().setExpandableContent(expandable);
				alert.getDialogPane().setExpanded(true);
				alert.getButtonTypes().setAll(continueButton, stopButton);

				Optional<ButtonType> answer = alert.showAndWait();
				continueRun.set(answer.isPresent() && answer.get() == continueButton);
			} finally {
				latch.countDown();
			}
		});
		latch.await();
		boolean chosen = continueRun.get();
		logger.info("{}: partial stage-1 decision -> {} (skipped={}, successful={})", getTaskName(),
				chosen ? "continue" : "stop", skippedSampleNames.size(), successfulSampleCount);
		return chosen;
	}

	@Override
	protected Boolean call() throws Exception {
		updateMessage(getTaskName());
		updateProgress(0, 100);
		setSecondBarBusy(true);

		try {
			// Same tool resolution the CLI / HPC entry points use
			// ({@link MagOpenSearchToolSetup}), so a missing Crystal-C / MSBooster /
			// FlashLFQ fails here for the same reason on every surface.
			MagOpenSearchToolSetup.Tools tools = openSearchParams.isEntrapmentEnabled()
					? MagOpenSearchToolSetup.prepareIdentification(openSearchParams, magDb, mrp)
					: MagOpenSearchToolSetup.prepare(openSearchParams, magDb, mrp);
			logger.info("{}: using MSFragger jar {}", getTaskName(), tools.msfraggerJar.getAbsolutePath());
			logger.info("{}: marker FASTA {}", getTaskName(), tools.markerFasta.getAbsolutePath());
			if (tools.flashLfqExe != null) {
				logger.info("{}: using FlashLFQ {}", getTaskName(), tools.flashLfqExe.getAbsolutePath());
			}
			if (openSearchParams.getSageExecutable() != null && !openSearchParams.getSageExecutable().isBlank()) {
				logger.info("{}: using Sage {}", getTaskName(), openSearchParams.getSageExecutable());
			}
			if (openSearchParams.getMsboosterJarPath() != null) {
				logger.info("{}: stage-1 rescoring with MSBooster {} and DIA-NN {}", getTaskName(),
						openSearchParams.getMsboosterJarPath(), openSearchParams.getDiannExecutablePath());
			}
			if (openSearchParams.getCrystalcJarPath() != null) {
				logger.info("{}: Crystal-C precursor refinement with {}", getTaskName(),
						openSearchParams.getCrystalcJarPath());
			}
			if (openSearchParams.getMsconvertExecutable() != null
					&& !openSearchParams.getMsconvertExecutable().isBlank()) {
				logger.info("{}: msconvert for .raw input {}", getTaskName(),
						openSearchParams.getMsconvertExecutable());
			}
			if (tools.funcDbFile != null) {
				logger.info("{}: using functional annotation database {}", getTaskName(),
						tools.funcDbFile.getAbsolutePath());
			} else {
				logger.info("{}: functional annotation database not configured - skipping Step I", getTaskName());
			}

			MagOpenSearchHeadlessRunner.Result result = MagOpenSearchHeadlessRunner.run(openSearchParams,
					tools.msfraggerJar, tools.markerFasta, magDb, tools.flashLfqExe,
					openSearchParams.getThreadCount(), tools.funcDbFile, createProgressCallback(),
					this::askContinueWithPartialIdentification);

			logger.info("{}: {} validated PSMs across {} accepted delta-mass groups", getTaskName(),
					result.validatedPsmCount, result.acceptedModificationCount);

			// Same fields the other MAG tasks populate before reporting, so the
			// shared report generator and the View Report button find their inputs.
			this.combinedFile = OpenSearchLayout.combinedDir(result.modeRoot);
			this.final_pep_txt = OpenSearchLayout.finalPeptides(result.modeRoot);
			this.final_pro_txt = OpenSearchLayout.razorProteins(result.modeRoot);

			if (openSearchParams.isEntrapmentEnabled()) {
				logger.info("{}: entrapment evaluation stops after identification; skipping visual reports",
						getTaskName());
				updateProgress(100, 100);
				setSecondBarBusy(false);
				updateMessage(getTaskName() + ": Entrapment evaluation completed");
				setStatus(getTaskName() + ": Entrapment evaluation completed");
				return true;
			}

			updateProgress(80, 100);
			updateMessage(getTaskName() + ": Generating visual reports...");
			setStatus(getTaskName() + ": Generating visual reports...");
			setSecondBarBusy(true);

			if (!generateFinalReportsIfStale(getTaskName())) {
				logger.error("{}: report generation failed", getTaskName());
				updateMessage(getTaskName() + ": Failed - report generation error");
				setStatus(getTaskName() + ": Failed - report generation error");
				Platform.runLater(() -> {
					if (secondBar != null) {
						secondBar.setProgress(0.0);
					}
				});
				return false;
			}

			updateProgress(100, 100);
			setSecondBarBusy(false);
			updateMessage(getTaskName() + ": Completed");
			setStatus(getTaskName() + ": Completed");
			return true;

		} catch (Exception e) {
			logger.error("{}: failed", getTaskName(), e);
			updateMessage(getTaskName() + ": Failed - " + e.getMessage());
			setStatus(getTaskName() + ": Failed");
			Platform.runLater(() -> {
				if (secondBar != null) {
					secondBar.setProgress(0.0);
				}
			});
			return false;
		}
	}

	private ProgressCallback createProgressCallback() {
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
					// Goes to logging.log via the task logger, and to the console
					// panel via the panel's messageProperty listener.
					logger.info(message);
					updateMessage(message);
					setStatus(message);
				}
			}

			@Override
			public void onComplete(Object result) {
				// Terminal state is reported by call() itself.
			}

			@Override
			public void onError(Exception error) {
				if (error != null) {
					logger.error("{}: runner reported an error", getTaskName(), error);
				}
			}
		};
	}

	/**
	 * The core runner reports progress at step boundaries only, so the second
	 * bar shows activity rather than a fraction - matching
	 * {@code DDAMagFragPipeTask}.
	 */
	private void setSecondBarBusy(boolean busy) {
		if (secondBar == null) {
			return;
		}
		Platform.runLater(() -> secondBar.setProgress(busy ? ProgressBar.INDETERMINATE_PROGRESS : 1.0));
	}

	private void setStatus(String text) {
		if (statusLabel != null) {
			Platform.runLater(() -> statusLabel.setText(text));
		}
	}
}
