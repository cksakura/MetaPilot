package uk.ac.quadram.metapilot.task.dda.mag.gui;

import java.io.File;

import org.apache.logging.log4j.Logger;

import javafx.application.Platform;
import javafx.scene.control.Label;
import javafx.scene.control.ProgressBar;

import uk.ac.quadram.metapilot.core.opensearch.OpenSearchLayout;
import uk.ac.quadram.metapilot.core.runner.MagOpenSearchHeadlessRunner;
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
 * the engine configuration's folder - with its own non-additive appender
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
			// the engine configuration's own folder, where its final tables and
			// report go, as MagOpenSearchHeadlessRunner.run() does again.
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
	 * The per-sample searches are shared by the engine configurations, so they
	 * sit under the engine search root rather than under the mode folder that
	 * {@code dbSearchResultFile} names.
	 */
	@Override
	protected File sampleFoldersRoot() {
		return OpenSearchLayout.searchResultsRoot(OpenSearchLayout.engineSearchRoot(new File(metaPar.getResult())));
	}

	/** The run's log belongs to the engine configuration it describes. */
	@Override
	protected File taskLogDir() {
		return modeRoot((DDAMagOpenSearchParameters) metaPar);
	}

	private static File modeRoot(DDAMagOpenSearchParameters params) {
		return OpenSearchLayout.modeRoot(OpenSearchLayout.engineSearchRoot(new File(params.getResult())),
				params.getSearchEngines());
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

	@Override
	protected Boolean call() throws Exception {
		updateMessage(getTaskName());
		updateProgress(0, 100);
		setSecondBarBusy(true);

		try {
			File msfraggerJar = mrp == null ? null : mrp.resolveMsfraggerJar();
			if (msfraggerJar == null || !msfraggerJar.exists()) {
				throw new IllegalStateException("MSFragger.jar not found inside the configured FragPipe installation. "
						+ "Set Resources > FragPipe to its bin/fragpipe launcher so the bundled jar can be detected.");
			}
			logger.info("{}: using MSFragger jar {}", getTaskName(), msfraggerJar.getAbsolutePath());

			// Resolved here rather than inside the runner so GUI, CLI and HPC all
			// take the executable from the same configured resource.
			if (openSearchParams.usesSage()) {
				String sagePath = mrp == null ? null
						: mrp.getToolPath(uk.ac.quadram.metapilot.resources.MetaResourceParameters.ToolType.SAGE);
				if (sagePath == null || sagePath.isBlank() || !new java.io.File(sagePath).exists()) {
					throw new IllegalStateException("Sage executable not configured or missing. Set Resources > Sage, "
							+ "or select \"" + DDAMagOpenSearchParameters.ENGINES_MSFRAGGER + "\".");
				}
				openSearchParams.setSageExecutable(sagePath);
				logger.info("{}: using Sage {}", getTaskName(), sagePath);
			}

			// Stage-1 rescoring tools, resolved from the same Resources entries the
			// rest of MetaPilot uses: MSBooster ships inside the configured
			// FragPipe installation, and DIA-NN is its prediction backend and is
			// already a configured resource for the DIA workflows. Resolved here,
			// at task setup, so "enabled but the tools are missing" is reported
			// before the search starts rather than after it - and it is reported,
			// not worked around: see MagOpenSearchHeadlessRunner.runStage1MsBooster.
			if (openSearchParams.isStage1RescoringEnabled() && openSearchParams.isStage1MsboosterEnabled()) {
				java.io.File msbooster = mrp == null ? null : mrp.resolveMsboosterJar();
				java.io.File unimod = mrp == null ? null : mrp.resolveUnimodObo();
				// FragPipe's bundled DIA-NN first. MSBooster needs the 1.8.x
				// prediction interface; the standalone DIA-NN configured for the DIA
				// workflows is typically 2.x and will fail.
				// FragPipe's bundled DIA-NN and nothing else. The standalone DIA-NN
				// configured for the DIA workflows is typically 2.x, which MSBooster
				// cannot drive: it was pointed at a 2.5.1 on 2026-09-10 and every
				// one of nine runs produced no predictions while reporting success.
				// There is no fallback here for that reason.
				java.io.File diann = mrp == null ? null : mrp.resolveMsboosterDiann();
				if (msbooster == null || diann == null || !diann.isFile()) {
					throw new IllegalStateException("Stage-1 MSBooster rescoring is enabled but "
							+ (msbooster == null ? "MSBooster.jar was not found in the configured FragPipe "
									+ "installation's tools folder"
									: "FragPipe bundles no DIA-NN 1.8.x under tools/diann/ - MSBooster drives "
											+ "DIA-NN through the 1.8.x prediction interface and a 2.x build "
											+ "silently predicts nothing")
							+ ". Set Resources > FragPipe to the installation that bundles both, or switch "
							+ "stage-1 MSBooster off deliberately.");
				}
				openSearchParams.setMsboosterJarPath(msbooster.getAbsolutePath());
				openSearchParams.setDiannExecutablePath(diann.getAbsolutePath());
				if (unimod != null) {
					openSearchParams.setUnimodOboPath(unimod.getAbsolutePath());
				}
				logger.info("{}: stage-1 rescoring with MSBooster {} and DIA-NN {}", getTaskName(),
						msbooster.getAbsolutePath(), diann.getAbsolutePath());
			}

			// Crystal-C, resolved from the same FragPipe installation. It refines
			// the round-two precursor masses before any delta mass is read;
			// without it the workflow reports isotope-selection and co-isolation
			// artifacts as modifications (measured on 12mix: 316,120 of 785,180
			// validated target PSMs sat at 0.05-5.5 Da, and 335 of 500 fitted
			// peaks lay inside +/-5.5 Da). Resolved here, at task setup, so
			// "enabled but the tool is missing" is reported before the search
			// starts rather than after it.
			if (openSearchParams.isCrystalcEnabled()) {
				java.io.File crystalc = mrp == null ? null : mrp.resolveCrystalcJar();
				if (crystalc == null) {
					throw new IllegalStateException("Crystal-C precursor refinement is enabled but no "
							+ "original-crystalc-*.jar was found in the configured FragPipe installation's "
							+ "tools folder. Set Resources > FragPipe to an installation that bundles it, or "
							+ "switch the step off deliberately.");
				}
				openSearchParams.setCrystalcJarPath(crystalc.getAbsolutePath());
				logger.info("{}: Crystal-C precursor refinement with {}", getTaskName(),
						crystalc.getAbsolutePath());
			}

			// Used only when Thermo .raw files are loaded; the runner converts them
			// to mzML first and reports a missing msconvert itself.
			String msconvertPath = mrp == null ? null
					: mrp.getToolPath(uk.ac.quadram.metapilot.resources.MetaResourceParameters.ToolType.MSCONVERT);
			if (msconvertPath != null && !msconvertPath.isBlank()) {
				openSearchParams.setMsconvertExecutable(msconvertPath);
				logger.info("{}: msconvert for .raw input {}", getTaskName(), msconvertPath);
			}

			File markerFasta = magDb.getHapFastaFile(openSearchParams.getGenomeSelectorMode());
			if (markerFasta == null || !markerFasta.exists()) {
				throw new IllegalStateException("MAG marker FASTA not found for genome selector mode "
						+ openSearchParams.getGenomeSelectorMode());
			}
			logger.info("{}: marker FASTA {}", getTaskName(), markerFasta.getAbsolutePath());

			// Resolved here, at task setup, for the same reason as every other tool
			// above: "enabled but missing" must be reported before the search
			// starts, not after Steps A-D have already run. Same resource entry and
			// getToolPath() call DDAMagPfindTask.lfQuant() uses
			// (mrp.getToolPath(ToolType.FLASHLFQ)), so both workflows' FlashLFQ step
			// is configured from the one place a person would look to fix it.
			String flashLfqPath = mrp == null ? null
					: mrp.getToolPath(uk.ac.quadram.metapilot.resources.MetaResourceParameters.ToolType.FLASHLFQ);
			File flashLfqExe = flashLfqPath == null || flashLfqPath.isBlank() ? null : new File(flashLfqPath);
			if (flashLfqExe == null || !flashLfqExe.exists()) {
				throw new IllegalStateException(
						"FlashLFQ executable not configured or missing. Set Resources > FlashLFQ.");
			}
			logger.info("{}: using FlashLFQ {}", getTaskName(), flashLfqExe.getAbsolutePath());

			// Same resource entry and getToolPath() call DDAMagFragPipeTask's own
			// resolveGuiToolFile(ToolType.FUNCDB) uses. Optional - a missing/blank
			// path resolves to null, and MagOpenSearchHeadlessRunner.run() already
			// skips functional annotation gracefully (with a log message) when
			// funcDbFile is null, exactly as MagFragPipeHeadlessRunner does.
			String funcDbPath = mrp == null ? null
					: mrp.getToolPath(uk.ac.quadram.metapilot.resources.MetaResourceParameters.ToolType.FUNCDB);
			File funcDbFile = funcDbPath == null || funcDbPath.isBlank() ? null : new File(funcDbPath);
			if (funcDbFile != null && funcDbFile.isFile()) {
				logger.info("{}: using functional annotation database {}", getTaskName(),
						funcDbFile.getAbsolutePath());
			} else {
				logger.info("{}: functional annotation database not configured - skipping Step I", getTaskName());
			}

			MagOpenSearchHeadlessRunner.Result result = MagOpenSearchHeadlessRunner.run(openSearchParams,
					msfraggerJar, markerFasta, magDb, flashLfqExe, openSearchParams.getThreadCount(),
					funcDbFile, createProgressCallback());

			logger.info("{}: {} validated PSMs across {} accepted delta-mass groups", getTaskName(),
					result.validatedPsmCount, result.acceptedModificationCount);

			// Same fields the other MAG tasks populate before reporting, so the
			// shared report generator and the View Report button find their inputs.
			this.combinedFile = OpenSearchLayout.combinedDir(result.modeRoot);
			this.final_pep_txt = OpenSearchLayout.finalPeptides(result.modeRoot);
			this.final_pro_txt = OpenSearchLayout.razorProteins(result.modeRoot);

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
