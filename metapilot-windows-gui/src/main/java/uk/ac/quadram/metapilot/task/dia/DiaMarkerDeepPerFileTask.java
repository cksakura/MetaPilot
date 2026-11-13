package uk.ac.quadram.metapilot.task.dia;

import uk.ac.quadram.metapilot.core.dia.PeptideListLibrarySupport;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Optional;
import java.util.Set;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.atomic.AtomicBoolean;

import org.apache.logging.log4j.Logger;

import javafx.application.Platform;
import javafx.scene.control.Alert;
import javafx.scene.control.Alert.AlertType;
import javafx.scene.control.ButtonBar;
import javafx.scene.control.ButtonType;
import javafx.scene.control.ProgressBar;
import javafx.scene.control.TextArea;
import javafx.scene.layout.GridPane;
import javafx.scene.layout.Priority;

import uk.ac.quadram.metapilot.core.dia.DiaDeepProteinMappingService;
import uk.ac.quadram.metapilot.core.mag.MagGenomeSelectionService;
import uk.ac.quadram.metapilot.core.mag.Marker4ResourcePreparer;
import uk.ac.quadram.metapilot.core.runner.DiaMarkerDeepPerFileHeadlessRunner;
import uk.ac.quadram.metapilot.core.runner.DiaNNNoQuantLibrarySearchService;
import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNParquetReader;
import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNPrecursor;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters.ToolType;
import uk.ac.quadram.metapilot.core.prodb.FastaManager;
import uk.ac.quadram.metapilot.task.dia.par.MetaParameterDiaMarker;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;
import uk.ac.quadram.metapilot.task.mag.MagDbItemResolver;
import uk.ac.quadram.metapilot.task.par.MetaData;

/**
 * DIA Deep Marker workflow.
 * <p>
 * High-level flow:
 * <ol>
 * <li>Search marker library against DIA files to get marker evidence parquet.</li>
 * <li>Run marker-based genome ranking to produce {@code genome_ranked.tsv}.</li>
 * <li>Run iterative rounds (top-N genomes per round) where each round builds one
 *     temporary FASTA/speclib and performs non-quant identification.</li>
 * <li>Merge peptides from all rounds and run one final quantifying search to produce
 *     {@code combined.parquet}, {@code combined.pr_matrix.tsv}, {@code combined.pg_matrix.tsv}.</li>
 * </ol>
 *
 * DIA-NN usage map in this class:
 * <ul>
 * <li>{@code addNoQuantLibTask(... marker.parquet ...)}: marker-only identification for genome ranking input.</li>
 * <li>{@code addTask(... combined.fasta -> combined.speclib ...)}: build predicted library from FASTA.</li>
 * <li>{@code addNoQuantLibTask(... roundN.parquet ...)}: round-level non-quant identification.</li>
 * <li>{@code addTask(... combined.parquet ...)}: final quantifying search over all DIA files.</li>
 * </ul>
 */
public class DiaMarkerDeepPerFileTask extends DiaMarkerCoreTask {

	private static final String TASK_NAME = "DIA Deep Marker analysis";
	/** Legacy default genomes per round. Prefer {@link MetaParameterDiaMarker#getGenomesPerRound()} when available. */
	private static final int DEFAULT_GENOMES_PER_ROUND = 5;
	/** Stop iterative rounds when identification count does not increase (like Deep Comprehensive). */
	private static final int MAX_ITERATIVE_ROUNDS = 100;

	/** Last round index from iterative search; used in exportReport to rebuild peptide→protein mapping from round parquets. */
	private int lastRoundForExport;

	/** True when the initial marker search fell back to per-file search. */
	private boolean initialSearchSeparate;
	/** Subset of diaFiles that produced results. Set when initialSearchSeparate is true. */
	private String[] goodDiaFiles;

	public DiaMarkerDeepPerFileTask(MetaParameterDiaMarker markerParams, MetaResourceParameters mrp,
			ProgressBar bar1, ProgressBar bar2) {
		this(withMode(markerParams, DiaRefinementMode.DEEP_PER_FILE), mrp, bar1, bar2, true);
	}

	protected DiaMarkerDeepPerFileTask(MetaParameterDiaMarker markerParams, MetaResourceParameters mrp,
			ProgressBar bar1, ProgressBar bar2, DiaRefinementMode refinementMode) {
		this(withMode(markerParams, refinementMode), mrp, bar1, bar2, true);
	}

	private DiaMarkerDeepPerFileTask(MetaParameterDiaMarker markerParams, MetaResourceParameters mrp,
			ProgressBar bar1, ProgressBar bar2, boolean modeAlreadyApplied) {
		super(markerParams, mrp, bar1, bar2, null);
	}

	private static MetaParameterDiaMarker withMode(MetaParameterDiaMarker markerParams, DiaRefinementMode refinementMode) {
		if (markerParams != null && refinementMode != null) {
			markerParams.setMarkerRefinementMode(refinementMode);
		}
		return markerParams;
	}

	private boolean isDeepPerFileMode() {
		if (!(metaPar instanceof MetaParameterDiaMarker)) {
			return false;
		}
		return ((MetaParameterDiaMarker) metaPar).getMarkerRefinementMode() == DiaRefinementMode.DEEP_PER_FILE;
	}

	private int getGenomesPerRound() {
		if (metaPar instanceof MetaParameterDiaMarker) {
			int v = ((MetaParameterDiaMarker) metaPar).getGenomesPerRound();
			return v > 0 ? v : DEFAULT_GENOMES_PER_ROUND;
		}
		return DEFAULT_GENOMES_PER_ROUND;
	}

	private static String rawBaseName(String path) {
		return DeepComprehensivePaths.runIdFromRawFilePath(path);
	}

	@Override
	protected String getTaskName() {
		return TASK_NAME;
	}

	@Override
	protected Logger getLogger() {
		return logger;
	}

	/**
	 * In DEEP_PER_FILE mode, {@code metaPar.getDbSearchResultFile()} is temporarily set
	 * to {@code <marker_deep>/marker/<rawname>} for the current file, so the iterative folder
	 * should be that directory itself (no extra "/marker" nesting).
	 */
	@Override
	protected File getIterativeResultFolder() {
		if (isDeepPerFileMode()) {
			return metaPar.getDbSearchResultFile();
		}
		return super.getIterativeResultFolder();
	}

	/**
	 * New mode: run Deep Marker independently per DIA raw file.
	 * Outputs per raw file:
	 *   marker_deep/marker/<rawname>/{marker.parquet, genome_ranked.tsv, roundN/...}
	 *
	 * After all per-file rounds finish, a single combined folder is generated (same as original deep mode):
	 *   marker_deep/combined/{combined.fasta, combined.speclib, combined.parquet, matrices...}
	 */
	@Override
	protected Boolean call() throws Exception {
		if (!isDeepPerFileMode()) {
			return super.call();
		}
		return callDeepPerFileViaCoreRunner();
	}

	private Boolean callDeepPerFileViaCoreRunner() throws Exception {
		logger.info(getTaskName() + ": starting DEEP_PER_FILE workflow");
		updateMessage(getTaskName() + " (per raw file)");
		updateProgress(0, 100);
		updateSecondaryProgress(ProgressBar.INDETERMINATE_PROGRESS);

		logger.info(getTaskName() + ": preparing report templates");
		updateMessage(getTaskName() + ": preparing report templates...");
		updateProgress(1, 100);
		if (!copyReportTemplatesInBackground()) {
			logger.error(getTaskName() + ": failed to copy report templates (per-file mode)");
			updateSecondaryProgress(0.0);
			return false;
		}

		logger.info(getTaskName() + ": preparing DIA input files");
		updateMessage(getTaskName() + ": preparing DIA input files...");
		updateProgress(2, 100);
		if (!spectraConvert()) {
			logger.error(getTaskName() + ": spectra conversion failed (per-file mode)");
			updateSecondaryProgress(0.0);
			return false;
		}

		File diannExe = resolveGuiToolFile(ToolType.DIANN);
		if (diannExe == null || !diannExe.exists()) {
			logger.error(getTaskName() + ": DIA-NN executable not found");
			updateSecondaryProgress(0.0);
			return false;
		}

		MetaParameterDiaMarker runnerParams = buildCoreRunnerParams();
		File magDbConfigFile = MagDbItemResolver.resolveConfigFile(null);
		File condaPrefix = inferCondaPrefixFromPython(resolveGuiToolFile(ToolType.PYTHON));
		File funcDbFile = resolveGuiToolFile(ToolType.FUNCDB);

		try {
			logger.info(getTaskName() + ": starting per-file Marker/Deep core runner");
			DiaMarkerDeepPerFileHeadlessRunner.Result result = DiaMarkerDeepPerFileHeadlessRunner.runWithoutPostAnalysis(
					runnerParams,
					diannExe,
					magDbConfigFile,
					createGuiProgressCallback(),
					condaPrefix,
					funcDbFile,
					this::askContinueWithPartialMarkerResults);
			logger.info(getTaskName() + ": per-file Marker/Deep core runner completed");

			this.resultFolderFile = result.markerDeepDir;
			this.metaPar.setDbSearchResultFile(result.markerDeepDir);
			this.combinedFile = result.combinedDir;
			this.quan_pep_file = result.prMatrixTsv;
			this.quan_pro_file = result.pgMatrixTsv;
			this.initialSearchSeparate = false;
			this.goodDiaFiles = this.diaFiles;

			if (runnerParams.isRepresentativeClusteringDryRun()) {
				updateProgress(100, 100);
				updateMessage(getTaskName() + ": grouping dry run completed; clustering reports are ready.");
				logger.info(getTaskName() + ": grouping dry run completed; skipping combined report export");
				updateSecondaryProgress(1.0);
				return true;
			}

			logger.info(getTaskName() + ": generating combined report");
			updateMessage(getTaskName() + ": generating combined report...");
			updateProgress(97, 100);
			if (!exportReport()) {
				logger.error(getTaskName() + ": combined report failed");
				updateSecondaryProgress(0.0);
				return false;
			}
			logger.info(getTaskName() + ": combined report completed");

			updateProgress(100, 100);
			updateMessage(getTaskName() + ": per-file rounds + global combined completed successfully!");
			logger.info(getTaskName() + ": per-file rounds + global combined completed successfully");
			updateSecondaryProgress(1.0);
			return true;
		} catch (Exception e) {
			logger.error(getTaskName() + ": DEEP_PER_FILE core runner failed", e);
			updateMessage(getTaskName() + ": Failed - " + e.getMessage());
			updateSecondaryProgress(0.0);
			throw e;
		}
	}

	/**
	 * GUI decision after the Deep per-file marker phase when some candidate files
	 * (all files, or representatives when clustering is on) lack genome selection outputs.
	 */
	private boolean askContinueWithPartialMarkerResults(List<String> missingRawNames, int successfulFileCount)
			throws InterruptedException {
		CountDownLatch latch = new CountDownLatch(1);
		AtomicBoolean continueRun = new AtomicBoolean(false);
		Platform.runLater(() -> {
			try {
				ButtonType continueButton = new ButtonType("Continue with successful files",
						ButtonBar.ButtonData.OK_DONE);
				ButtonType stopButton = new ButtonType("Stop", ButtonBar.ButtonData.CANCEL_CLOSE);
				Alert alert = new Alert(AlertType.CONFIRMATION);
				alert.setTitle("Incomplete marker results");
				alert.setHeaderText("Some files have no genome selection result");
				alert.setContentText(
						"Marker search finished, but " + missingRawNames.size()
								+ " file(s) are missing genome_evidence.tsv / genome_ranked.tsv / genome_selected.tsv.\n\n"
								+ successfulFileCount
								+ " file(s) have complete results.\n\n"
								+ "Stop now if you want to replace the failing files, or continue and build the "
								+ "combined database from the successful marker results only.");

				TextArea details = new TextArea(String.join("\n", missingRawNames));
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
		logger.info(getTaskName() + ": partial marker result decision -> {} (missing={}, successful={})",
				chosen ? "continue" : "stop", missingRawNames.size(), successfulFileCount);
		return chosen;
	}

	/**
	 * Run per-file iterative phase sequentially (one raw file at a time).
	 * Used when total RAM is below the per-file budget (32 GB) or only one file.
	 *
	 * <p>Each per-file workflow runs in its own child {@link DiaMarkerDeepPerFileTask}
	 * instance using the supplied per-task thread count.
	 */
	private boolean runPerFileSequential(String[] originalDiaFiles,
			java.util.Map<String, File> perFileDirs,
			java.util.Map<String, String> perFileNames,
			java.util.concurrent.ConcurrentHashMap<String, Integer> perFileLastRound,
			int perTaskThreads) {
		int total = originalDiaFiles.length;
		for (int i = 0; i < total; i++) {
			String raw = originalDiaFiles[i];
			File perFileDir = perFileDirs.get(raw);
			String rawName = perFileNames.get(raw);

			int start = 10 + (int) Math.floor((double) i * 80.0 / (double) total);
			int end = 10 + (int) Math.floor((double) (i + 1) * 80.0 / (double) total);
			updateMessage(getTaskName() + ": per-file workflow (" + (i + 1) + "/" + total + "): " + rawName);
			updateProgress(start, 100);

			int lastRound = runOnePerFileWorker(raw, perFileDir, perTaskThreads);
			if (lastRound <= 0) {
				logger.error(getTaskName() + ": per-file workflow failed for {}", raw);
				return false;
			}
			perFileLastRound.put(raw, lastRound);
			updateProgress(end, 100);
		}
		return true;
	}

	/**
	 * Run per-file iterative phase in parallel using a fixed-size thread pool.
	 *
	 * <p>Pool size is computed from total system RAM (32 GB per concurrent task);
	 * each concurrent worker is given {@code totalThreads / poolSize} DIA-NN threads
	 * so the combined CPU/RAM load stays bounded.
	 *
	 * <p>If any worker fails, all submitted workers are awaited so we have a clean
	 * shutdown before returning false.
	 */
	private boolean runPerFileParallel(String[] originalDiaFiles,
			java.util.Map<String, File> perFileDirs,
			java.util.Map<String, String> perFileNames,
			java.util.concurrent.ConcurrentHashMap<String, Integer> perFileLastRound,
			int poolSize, int perTaskThreads) {
		final int total = originalDiaFiles.length;
		logger.info(getTaskName() + ": running {} per-file workflows in parallel (poolSize={}, perTaskThreads={})",
				total, poolSize, perTaskThreads);

		java.util.concurrent.ExecutorService executor = java.util.concurrent.Executors.newFixedThreadPool(poolSize,
				r -> {
					Thread t = new Thread(r, "DeepPerFile-Worker");
					t.setDaemon(true);
					return t;
				});
		java.util.concurrent.atomic.AtomicInteger completed = new java.util.concurrent.atomic.AtomicInteger(0);
		java.util.concurrent.atomic.AtomicBoolean anyFailed = new java.util.concurrent.atomic.AtomicBoolean(false);
		java.util.List<java.util.concurrent.Future<?>> futures = new java.util.ArrayList<>();

		try {
			for (String raw : originalDiaFiles) {
				final File perFileDir = perFileDirs.get(raw);
				final String rawName = perFileNames.get(raw);
				final String rawFile = raw;
				futures.add(executor.submit(() -> {
					try {
						if (anyFailed.get()) {
							return;
						}
						logger.info(getTaskName() + ": [parallel] starting per-file workflow for {} (perTaskThreads={})",
								rawName, perTaskThreads);
						int lastRound = runOnePerFileWorker(rawFile, perFileDir, perTaskThreads);
						if (lastRound <= 0) {
							logger.error(getTaskName() + ": [parallel] per-file workflow failed for {}", rawFile);
							anyFailed.set(true);
							return;
						}
						perFileLastRound.put(rawFile, lastRound);
						int done = completed.incrementAndGet();
						int progress = 10 + (int) Math.floor((double) done * 80.0 / (double) total);
						updateMessage(getTaskName() + ": [parallel] " + done + "/" + total + " complete (latest: "
								+ rawName + ")");
						updateProgress(progress, 100);
						logger.info(getTaskName() + ": [parallel] completed per-file workflow for {} (lastRound={})",
								rawName, lastRound);
					} catch (Throwable t) {
						logger.error(getTaskName() + ": [parallel] exception in per-file worker for " + rawFile, t);
						anyFailed.set(true);
					}
				}));
			}

			for (java.util.concurrent.Future<?> f : futures) {
				try {
					f.get();
				} catch (java.util.concurrent.ExecutionException ee) {
					logger.error(getTaskName() + ": [parallel] worker future failed: {}", ee.toString(), ee);
					anyFailed.set(true);
				} catch (InterruptedException ie) {
					Thread.currentThread().interrupt();
					logger.error(getTaskName() + ": [parallel] interrupted while awaiting worker", ie);
					anyFailed.set(true);
				}
			}
		} finally {
			executor.shutdown();
			try {
				if (!executor.awaitTermination(5, java.util.concurrent.TimeUnit.MINUTES)) {
					executor.shutdownNow();
				}
			} catch (InterruptedException ignored) {
				executor.shutdownNow();
				Thread.currentThread().interrupt();
			}
		}

		if (anyFailed.get()) {
			logger.error(getTaskName() + ": one or more per-file workers failed");
			return false;
		}
		return true;
	}

	/**
	 * Run the iterative phase (marker search + iterative rounds) for a single raw
	 * file in an isolated child {@link DiaMarkerDeepPerFileTask} instance.
	 *
	 * <p>Each child has its own {@code DiaNNTask}, {@code DiannParameter},
	 * {@code DiaLibSearchPar}, and per-file metadata, which makes it safe to run
	 * multiple workers concurrently. The child is configured with
	 * {@code threadCount = perTaskThreads} so DIA-NN sub-processes use the
	 * appropriate share of CPU. The child writes its output under
	 * {@code perFileDir} (the same folder layout used in sequential mode).
	 *
	 * <p>The child is <b>not</b> executed via {@link javafx.concurrent.Task#run()};
	 * we call {@link #iterativeSearch()} and {@link #runRoundsOnly()} directly so we
	 * skip {@code spectraConvert} (already done by the parent), {@code createSSLibrary}
	 * (already done by the parent), and the global combined search (parent runs that
	 * once after all per-file workers finish).
	 *
	 * @return last round index (>0) on success; <=0 on failure.
	 */
	private int runOnePerFileWorker(String raw, File perFileDir, int perTaskThreads) {
		try {
			MetaParameterDiaMarker parentP = (MetaParameterDiaMarker) metaPar;

			String[] childRawFiles = new String[] { raw };
			String[] childExpNames = new String[] { rawBaseName(raw) };
			int[] fractions = new int[] { 1 };
			int[] replicates = new int[] { 1 };
			MetaData childMd = new MetaData(childRawFiles, childExpNames, fractions, replicates);

			MetaParameterDiaMarker childP = new MetaParameterDiaMarker(parentP, parentP.getMagDbItem(),
					parentP.isLibrarySearch(), parentP.isSelfModelSearch(),
					parentP.getLibrary(), parentP.isUseModel(), parentP.getModels());
			// Use DEEP_PER_FILE so the child's getIterativeResultFolder() returns its
			// dbSearchResultFile (perFileDir) directly, matching the per-file folder layout.
			childP.setMarkerRefinementMode(DiaRefinementMode.DEEP_PER_FILE);
			childP.setGenomesPerRound(parentP.getGenomesPerRound());
			childP.setMExcision(parentP.isMExcision());
			childP.setCarbamiC(parentP.isCarbamiC());
			childP.setOxidationM(parentP.isOxidationM());
			childP.setAcetyProN(parentP.isAcetyProN());
			childP.setPhosphoSTY(parentP.isPhosphoSTY());
			childP.setKGG(parentP.isKGG());
			childP.setUsePeptideBasedRound(parentP.isUsePeptideBasedRound());
			childP.setPeptideDetectabilityThreshold(parentP.getPeptideDetectabilityThreshold());
			childP.setMetadata(childMd);
			childP.setDbSearchResultFile(perFileDir);
			childP.setReportDir(null);
			childP.setThreadCount(perTaskThreads);
			childP.setResult(parentP.getResult());

			DiaMarkerDeepPerFileTask child = new DiaMarkerDeepPerFileTask(childP, this.mrp,
					this.progressBar, this.secondBar);
			// Skip spectraConvert in child (parent already converted): set diaFiles directly.
			child.diaFiles = new String[] { raw };
			child.initialSearchSeparate = false;
			child.goodDiaFiles = child.diaFiles;

			if (!child.iterativeSearch()) {
				getLogger().error(getTaskName() + ": child iterativeSearch failed for {}", raw);
				return -1;
			}
			int lastRound = child.runRoundsOnly();
			if (lastRound <= 0) {
				getLogger().error(getTaskName() + ": child runRoundsOnly returned {} for {}", lastRound, raw);
				return -1;
			}
			return lastRound;
		} catch (Throwable t) {
			getLogger().error(getTaskName() + ": runOnePerFileWorker exception for " + raw, t);
			return -1;
		}
	}

	/**
	 * Run only the iterative round loop (no final combined quant search).
	 *
	 * @return lastRound (>0) on success; <=0 on failure.
	 */
	private int runRoundsOnly() {
		File markerFile = getIterativeResultFolder();
		File rankedFile = new File(markerFile, "genome_ranked.tsv");
		if (!rankedFile.exists()) {
			logger.error(getTaskName() + ": genome_ranked.tsv not found: {}", rankedFile.getAbsolutePath());
			return -1;
		}

		String cut = diannPar != null && diannPar.getCut() != null ? diannPar.getCut() : "KR";
		DeepIterativeCommonSupport.IterativeRoundLoopResult loop = DeepIterativeCommonSupport.runIterativeRoundLoop(
				rankedFile, getGenomesPerRound(), MAX_ITERATIVE_ROUNDS, "Round", cut,
				initialSearchSeparate, initialSearchSeparate ? goodDiaFiles.length : 0,
				this::updateMessage,
				(round, genomes) -> runRoundAppendPreviousPeptides(round, markerFile, genomes),
				(inc) -> {
					StopRuleResult local = computeStopRoundFromIncrements(inc);
					return new DeepIterativeStopRule.StopRuleResult(local.stopRound, local.nRounds, local.kChord,
							local.kPiece, local.kKnee, local.k0, local.tMin);
				},
				logger, getTaskName());

		if (loop.lastRound <= 0) {
			logger.error(getTaskName() + ": no round completed");
			return -1;
		}

		StopRuleResult lastStopRuleResult = loop.lastStopRuleResult == null ? null
				: new StopRuleResult(loop.lastStopRuleResult.stopRound, loop.lastStopRuleResult.nRounds,
						loop.lastStopRuleResult.kChord, loop.lastStopRuleResult.kPiece, loop.lastStopRuleResult.kKnee,
						loop.lastStopRuleResult.k0, loop.lastStopRuleResult.tMin);
		writeIterativeRoundsInfo(markerFile, loop.roundsRecord, loop.peptidesInRoundRecord, loop.incrementsPerRound,
				loop.cumulativeRecord, loop.stopReason, lastStopRuleResult);

		return loop.lastRound;
	}

	private boolean runGlobalCombinedFromPerFileRounds(
			String[] originalDiaFiles,
			java.util.LinkedHashMap<String, File> perFileDirs,
			java.util.LinkedHashMap<String, Integer> perFileLastRound,
			File deepRootDir
	) {
		if (originalDiaFiles == null || originalDiaFiles.length == 0) {
			logger.error(getTaskName() + ": no DIA files for global combined search");
			return false;
		}
		if (perFileDirs == null || perFileDirs.isEmpty() || perFileLastRound == null || perFileLastRound.isEmpty()) {
			logger.error(getTaskName() + ": per-file rounds are missing; cannot build global combined");
			return false;
		}

		// 1) Aggregate peptide->proteins across ALL per-file rounds.
		HashMap<String, HashSet<String>> pepProMap = new HashMap<>();
		for (String raw : perFileDirs.keySet()) {
			File perFileDir = perFileDirs.get(raw);
			Integer lastRound = perFileLastRound.get(raw);
			if (perFileDir == null || lastRound == null || lastRound.intValue() <= 0) {
				continue;
			}

			// Temporarily align state for reading round parquets in that folder (important if per-file search outputs exist).
			String[] oldDiaFiles = this.diaFiles;
			boolean oldSeparate = this.initialSearchSeparate;
			String[] oldGood = this.goodDiaFiles;
			try {
				this.diaFiles = new String[] { raw };
				this.goodDiaFiles = this.diaFiles;
				this.initialSearchSeparate = false;

				HashMap<String, HashSet<String>> one = buildPepToProteinsFromRounds(perFileDir, lastRound.intValue());
				for (java.util.Map.Entry<String, HashSet<String>> e : one.entrySet()) {
					if (e.getKey() == null || e.getKey().isEmpty() || e.getValue() == null || e.getValue().isEmpty()) {
						continue;
					}
					HashSet<String> set = pepProMap.computeIfAbsent(e.getKey(), k -> new HashSet<>());
					set.addAll(e.getValue());
				}
			} finally {
				this.diaFiles = oldDiaFiles;
				this.initialSearchSeparate = oldSeparate;
				this.goodDiaFiles = oldGood;
			}
		}

		if (pepProMap.isEmpty()) {
			logger.error(getTaskName() + ": no peptides collected from per-file round parquets for global combined FASTA");
			return false;
		}

		// 2) Write combined FASTA + speclib under the deep root combined/ folder (same as original mode).
		File combinedDir = new File(deepRootDir, "combined");
		if (!combinedDir.exists()) combinedDir.mkdirs();
		File combinedParquet = new File(combinedDir, "combined.parquet");
		File combinedPepMatrix = new File(combinedDir, "combined.pr_matrix.tsv");
		File combinedFasta = new File(combinedDir, "combined.fasta");
		File combinedSpeclib = new File(combinedDir, "combined.speclib");

		if (!DeepIterativeCommonSupport.writeCombinedFastaFromPeptideProteinMap(combinedFasta, pepProMap, logger, getTaskName())) {
			return false;
		}

		if (!DeepIterativeCommonSupport.ensureCombinedSpeclib(combinedDir, combinedSpeclib, () -> {
			diannPar.setThreads(threadCount);
			diannPar.setExportDecoy(true);
			// combined.fasta holds one entry per round peptide: no protein N-terminal rules
			PeptideListLibrarySupport.runWithPeptideEntryRules(diannPar, () -> {
				diaNNTask.addTask(diannPar, combinedFasta.getAbsolutePath(), combinedSpeclib.getAbsolutePath(), false);
				diaNNTask.run(12);
			});
		}, logger, getTaskName())) {
			return false;
		}

		// 3) Final quantifying search across ALL DIA files using the combined speclib (+ host library, if any).
		final String[] finalLibs = finalSearchLibraries(combinedSpeclib, combinedDir);
		if (finalLibs == null) {
			return false;
		}
		boolean combinedOk = DeepIterativeCommonSupport.ensureCombinedParquet(combinedParquet, () -> {
			diaLibSearchPar.setMbr(true);
			diaLibSearchPar.setDiaFiles(originalDiaFiles);
			diaLibSearchPar.setThreadCount(threadCount);
			diaLibSearchPar.setqValue(0.01);
			diannPar.setThreads(threadCount);
			diannPar.setMatrices(true);
			diannPar.setRelaxed_prot_inf(true);
			diannPar.setReanalyse(true);
			diannPar.setExportDecoy(true);
			logger.info(getTaskName() + ": global combined search using combined.speclib {} → {} ({} files)",
					combinedSpeclib.getAbsolutePath(), combinedParquet.getAbsolutePath(), originalDiaFiles.length);
			diaNNTask.addTask(diannPar, diaLibSearchPar, combinedParquet.getAbsolutePath(), finalLibs);
			diaNNTask.run(originalDiaFiles.length * 12);
		});
		if (!combinedOk || !combinedParquet.exists() || !combinedPepMatrix.exists()) {
			logger.error(getTaskName() + ": global combined search failed - combined outputs missing in {}", combinedDir.getAbsolutePath());
			return false;
		}

		this.combinedFile = combinedDir;
		this.quan_pep_file = combinedPepMatrix;
		this.quan_pro_file = new File(combinedDir, "combined.pg_matrix.tsv");
		return true;
	}

	@Override
	protected void initResultFolderIfNeeded() {
		if (metaPar == null) {
			return;
		}
		if (metaPar.getDbSearchResultFile() != null) {
			this.resultFolderFile = metaPar.getDbSearchResultFile();
			return;
		}
		if (!(metaPar instanceof MetaParameterDiaMarker)) {
			super.initResultFolderIfNeeded();
			return;
		}
		MetaParameterDiaMarker markerParams = (MetaParameterDiaMarker) metaPar;
		File base = new File(markerParams.getResult() != null ? markerParams.getResult() : ".");
		File dir = new File(base, "marker_deep");
		if (!dir.exists()) {
			dir.mkdirs();
			logger.info("Created result directory for DIA Deep Marker workflow: {}", dir.getAbsolutePath());
		}
		markerParams.setDbSearchResultFile(dir);
		this.resultFolderFile = dir;
	}

	/**
	 * Same as marker workflow: search marker library → marker.parquet, then run
	 * MarkerGenomeSelector to produce genome_evidence.tsv, genome_ranked.tsv, genome_selected.tsv.
	 * Stops after that; does not build combined.fasta or combined.speclib.
	 */
	@Override
	protected boolean iterativeSearch() {
		logger.info("Starting iterative search for DIA Deep Marker (search + genome selection only)...");
		System.out.println(new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss").format(new Date()) + "\t" + getTaskName()
				+ ": iterative search started");

		File markerLibFile = magDb.resolveMarker4LibraryFile();
		if (markerLibFile == null || !markerLibFile.exists()) {
			logger.error(getTaskName() + ": Marker library file not found for iterative search");
			System.err.println(new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss").format(new Date()) + "\t" + getTaskName()
					+ ": Marker library file not found");
			return false;
		}

		File markerFile = getIterativeResultFolder();
		if (!markerFile.exists()) {
			markerFile.mkdir();
		}

		File tempParquetFile = new File(markerFile, "marker.parquet");
		if (tempParquetFile.exists()) {
			logger.info("Marker library search result already exists in {}", tempParquetFile.getAbsolutePath());
		} else if (hasExistingPerFileParquets(markerFile)) {
			initialSearchSeparate = true;
			goodDiaFiles = collectGoodDiaFilesFromPerFileParquets(markerFile);
			int maxConcurrent = DeepIterativeCommonSupport.computeMaxConcurrentTasks();
			logger.info(getTaskName()
					+ ": per-file parquets already exist from previous run; {} good files; rounds will use per-file search (max concurrent: {})",
					goodDiaFiles.length, maxConcurrent);
		} else {
			DiaNNNoQuantLibrarySearchService.configureNoQuantSearch(
					diannPar, diaLibSearchPar, diaFiles, threadCount, true, false);
			logger.info("Searching marker library against DIA files...");
			logger.info("  Marker library: {}", markerLibFile.getAbsolutePath());
			logger.info("  DIA files count: {}", diaFiles != null ? diaFiles.length : 0);
			diaNNTask.addNoQuantLibTask(diannPar, diaLibSearchPar, tempParquetFile.getAbsolutePath(),
					markerLibFile.getAbsolutePath());
			diaNNTask.run(diaFiles.length * 12);

			if (!tempParquetFile.exists()) {
				logCombinedSearchFailed(tempParquetFile, "marker library search");
				return false;
			}
		}
		if (!initialSearchSeparate) {
			goodDiaFiles = diaFiles;
		}

		// Genome selection only (no combined.fasta / combined.speclib)
		logger.info("Reading PSMs from iterative search results...");
		DiaNNPrecursor[] precursors;
		if (tempParquetFile.exists()) {
			DiaNNParquetReader parquetReader = new DiaNNParquetReader(tempParquetFile.getAbsolutePath(),
					diannPar.getCut());
			precursors = parquetReader.getDiaNNPrecursors(0.01);
		} else {
			logger.warn(getTaskName()
					+ ": Marker library combined parquet not created; using per-file parquet results for genome selection");
			precursors = readPrecursorsFromPerFileParquets(markerFile, 0.01);
		}
		logger.info("Found {} precursors with Global.Q.Value < 0.01", precursors != null ? precursors.length : 0);

		if (precursors == null || precursors.length == 0) {
			logger.warn(getTaskName() + ": No precursors found in iterative search, cannot select genomes");
			return false;
		}

		MetaParameterDiaMarker markerParams = (metaPar instanceof MetaParameterDiaMarker) ? (MetaParameterDiaMarker) metaPar : null;
		// Map markerMode to MarkerGenomeSelector knee mode (FOCUS vs COMPREHENSIVE)
		// DEEP_PER_FILE behaves like DEEP for genome selection (per-file ranking with comprehensive knee).
		DiaRefinementMode markerMode = markerParams != null ? markerParams.getMarkerRefinementMode() : DiaRefinementMode.CORE;
		String mode;
		if (markerMode.usesComprehensiveSelection()) {
			mode = MagGenomeSelectionService.SELECTION_STRATEGY_COMPREHENSIVE;
		} else {
			mode = MagGenomeSelectionService.SELECTION_STRATEGY_FOCUS;
		}
		logger.info("Selecting genomes using MarkerGenomeSelector ({} mode)...", mode);

		File currentFile = magDb.getCurrentFile();
		try {
			Marker4ResourcePreparer.ensureSelectionResources(magDb, createGuiProgressCallback());
		} catch (IOException e) {
			logger.error(getTaskName() + ": failed to prepare Marker4 resources", e);
			return false;
		}
		File markerFasta = new File(currentFile, "Marker4.fasta");
		File genomesTaxaTsv = new File(currentFile, "genomes-taxa.tsv");
		if (!markerFasta.exists() || !genomesTaxaTsv.exists()) {
			logger.error(getTaskName() + ": Marker FASTA or genomes-taxa.tsv not found for genome selection");
			System.err.println(new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss").format(new Date()) + "\t" + getTaskName()
					+ ": Marker FASTA or genomes-taxa.tsv not found");
			return false;
		}

		File evidenceFile = new File(markerFile, "genome_evidence.tsv");
		File rankedFile = new File(markerFile, "genome_ranked.tsv");
		File selectedFile = new File(markerFile, "genome_selected.tsv");
		MagGenomeSelectionService.Request selectionRequest = new MagGenomeSelectionService.Request();
		selectionRequest.magDb = magDb;
		selectionRequest.genomeSelectorMode = MagDbItem.MARKER4;
		selectionRequest.genomeSelectionStrategy = mode;
		selectionRequest.diaPrecursors = precursors;
		selectionRequest.markerFasta = markerFasta;
		selectionRequest.genomesTaxaTsv = genomesTaxaTsv;
		selectionRequest.evidenceFile = evidenceFile;
		selectionRequest.rankedFile = rankedFile;
		selectionRequest.selectedFile = selectedFile;
		selectionRequest.maximumSelected = 300;
		try {
			MagGenomeSelectionService.selectGenomes(selectionRequest, createGuiProgressCallback());
		} catch (Exception e) {
			logger.error(getTaskName() + ": genome selection failed", e);
			return false;
		}

		logger.info(getTaskName() + ": genome selection completed (genome_evidence.tsv, genome_ranked.tsv, genome_selected.tsv); stopping before combined.fasta/speclib");
		System.out.println(new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss").format(new Date()) + "\t" + getTaskName()
				+ ": iterative search finished (genome selection only)");
		return true;
	}

	/**
	 * Reuse the same per-file parquet aggregation helper as the marker workflow.
	 */
	private DiaNNPrecursor[] readPrecursorsFromPerFileParquets(File dir, Double qValue) {
		String cut = diannPar != null && diannPar.getCut() != null ? diannPar.getCut() : "KR";
		return DeepIterativeCommonSupport.readPrecursorsFromPerFileParquets(dir, qValue, diaFiles, cut, logger,
				getTaskName());
	}

	/**
	 * Round controller after {@code genome_ranked.tsv} is available.
	 * Uses shared round-loop logic to:
	 * <ol>
	 * <li>take next genome window from rank list,</li>
	 * <li>run one non-quant identification round,</li>
	 * <li>evaluate stop rule based on incremental peptide gain,</li>
	 * <li>perform final combined quantifying search.</li>
	 * </ol>
	 */
	@Override
	protected boolean iden() {
		File markerFile = getIterativeResultFolder();
		File rankedFile = new File(markerFile, "genome_ranked.tsv");
		if (!rankedFile.exists()) {
			logger.error(getTaskName() + ": genome_ranked.tsv not found: {}", rankedFile.getAbsolutePath());
			return false;
		}

		String cut = diannPar != null && diannPar.getCut() != null ? diannPar.getCut() : "KR";
		DeepIterativeCommonSupport.IterativeRoundLoopResult loop = DeepIterativeCommonSupport.runIterativeRoundLoop(
				rankedFile, getGenomesPerRound(), MAX_ITERATIVE_ROUNDS, "Round", cut,
				initialSearchSeparate, initialSearchSeparate ? goodDiaFiles.length : 0,
				this::updateMessage,
				(round, genomes) -> runRoundAppendPreviousPeptides(round, markerFile, genomes),
				(inc) -> {
					StopRuleResult local = computeStopRoundFromIncrements(inc);
					return new DeepIterativeStopRule.StopRuleResult(local.stopRound, local.nRounds, local.kChord,
							local.kPiece, local.kKnee, local.k0, local.tMin);
				},
				logger, getTaskName());

		if (loop.lastRound <= 0) {
			logger.error(getTaskName() + ": no round completed");
			return false;
		}

		StopRuleResult lastStopRuleResult = loop.lastStopRuleResult == null ? null
				: new StopRuleResult(loop.lastStopRuleResult.stopRound, loop.lastStopRuleResult.nRounds,
						loop.lastStopRuleResult.kChord, loop.lastStopRuleResult.kPiece, loop.lastStopRuleResult.kKnee,
						loop.lastStopRuleResult.k0, loop.lastStopRuleResult.tMin);
		int lastRound = loop.lastRound;
		int prevPeptideCount = loop.previousPeptideCount;
		writeIterativeRoundsInfo(markerFile, loop.roundsRecord, loop.peptidesInRoundRecord, loop.incrementsPerRound,
				loop.cumulativeRecord, loop.stopReason, lastStopRuleResult);

		logger.info(getTaskName() + ": iterative rounds finished at round {} ({} peptides); running combined search using combined FASTA built from all round parquets...", lastRound, prevPeptideCount);

		// Combined folder:
		// - normal mode: under <marker_deep>/combined
		// - per-file mode: under <marker_deep>/marker/<rawname>/combined
		File combinedDir = isDeepPerFileMode()
				? new File(markerFile, "combined")
				: new File(metaPar.getDbSearchResultFile(), "combined");
		if (!combinedDir.exists()) combinedDir.mkdirs();
		File combinedParquet = new File(combinedDir, "combined.parquet");
		File combinedPepMatrix = new File(combinedDir, "combined.pr_matrix.tsv");
		File combinedFasta = new File(combinedDir, "combined.fasta");
		File combinedSpeclib = new File(combinedDir, "combined.speclib");

		if (!combinedParquet.exists() || !combinedPepMatrix.exists()) {
			// 1) Aggregate peptide→all proteins from all round{N}.parquet
			HashMap<String, HashSet<String>> pepProMap = buildPepToProteinsFromRounds(markerFile, lastRound);
			if (pepProMap.isEmpty()) {
				logger.error(getTaskName() + ": no peptides collected from round parquets for combined FASTA");
				return false;
			}
			this.lastRoundForExport = lastRound;

			// 2) Write combined FASTA: one entry per unique peptide, using one of its proteins as the FASTA protein name.
			boolean hadCombinedFasta = combinedFasta.exists() && combinedFasta.length() > 0;
			if (!DeepIterativeCommonSupport.writeCombinedFastaFromPeptideProteinMap(combinedFasta, pepProMap, logger,
					getTaskName())) {
				return false;
			}
			if (!hadCombinedFasta) {
				logger.info(getTaskName() + ": combined FASTA written: {} ({} bytes)", combinedFasta.getAbsolutePath(), combinedFasta.length());
			} else {
				logger.info(getTaskName() + ": using existing combined FASTA at {}", combinedFasta.getAbsolutePath());
			}

			// DIA-NN call purpose: build library for the final quantifying pass
			// from round-aggregated FASTA (combined.fasta -> combined.speclib).
			if (!DeepIterativeCommonSupport.ensureCombinedSpeclib(combinedDir, combinedSpeclib, () -> {
					diannPar.setThreads(threadCount);
					diannPar.setExportDecoy(true);
					// combined.fasta holds one entry per round peptide: no protein N-terminal rules
					PeptideListLibrarySupport.runWithPeptideEntryRules(diannPar, () -> {
						diaNNTask.addTask(diannPar, combinedFasta.getAbsolutePath(), combinedSpeclib.getAbsolutePath(), false);
						diaNNTask.run(12);
					});
				}, logger, getTaskName())) {
				return false;
			}
			logger.info(getTaskName() + ": combined.speclib ready: {} ({} bytes)", combinedSpeclib.getAbsolutePath(), combinedSpeclib.length());

			String[] finalSearchFiles = initialSearchSeparate ? goodDiaFiles : diaFiles;
			final String[] finalLibs = finalSearchLibraries(combinedSpeclib, combinedDir);
			if (finalLibs == null) {
				return false;
			}
			boolean combinedOk = DeepIterativeCommonSupport.ensureCombinedParquet(combinedParquet, () -> {
				diaLibSearchPar.setMbr(true);
				diaLibSearchPar.setDiaFiles(finalSearchFiles);
				diaLibSearchPar.setThreadCount(threadCount);
				diaLibSearchPar.setqValue(0.01);
				diannPar.setThreads(threadCount);
				diannPar.setMatrices(true);
				diannPar.setRelaxed_prot_inf(true);
				diannPar.setReanalyse(true);
				diannPar.setExportDecoy(true);
				logger.info(getTaskName() + ": combined search using combined.speclib {} → {} ({} files)",
						combinedSpeclib.getAbsolutePath(), combinedParquet.getAbsolutePath(), finalSearchFiles.length);
				diaNNTask.addTask(diannPar, diaLibSearchPar, combinedParquet.getAbsolutePath(), finalLibs);
				diaNNTask.run(finalSearchFiles.length * 12);
			});
			if (!combinedOk) {
				logger.error(getTaskName() + ": combined search failed - combined.parquet not created");
				return false;
			}
		} else {
			logger.info(getTaskName() + ": combined result already exists, skipping combined search");
			if (lastRoundForExport <= 0) {
				int gfc = goodDiaFiles != null ? goodDiaFiles.length : (diaFiles != null ? diaFiles.length : 0);
				for (int r = 1; r <= MAX_ITERATIVE_ROUNDS; r++) {
					File rd = new File(markerFile, "round" + r);
					if (rd.isDirectory() && (new File(rd, "round" + r + ".parquet").isFile()
							|| DeepIterativeCommonSupport.hasRoundPerFileParquets(rd, gfc))) {
						lastRoundForExport = r;
					}
				}
			}
		}

		combinedFile = combinedDir;
		this.quan_pep_file = new File(combinedDir, "combined.pr_matrix.tsv");
		this.quan_pro_file = new File(combinedDir, "combined.pg_matrix.tsv");
		logger.info(getTaskName() + ": combined search completed: {} ({} bytes); unified post-analysis will use combined folder", combinedParquet.getAbsolutePath(), combinedParquet.length());
		System.out.println(new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss").format(new Date()) + "\t" + getTaskName() + ": identification finished (combined)");
		return true;
	}

	/**
	 * Writes iterative_rounds_info.tsv so users can see why the round search stopped.
	 * Columns: Round, Peptides_in_round, New_peptides, Cumulative_peptides, Stop_reason (on last row only).
	 * Stop is determined by relative gain: New_peptides/Cumulative_peptides < 1.5% (plateau on total discovery).
	 * If stopRule is non-null (e.g. last evaluation when rounds>=8), writes kChord, kPiece, kKnee, k0, tMin in comment lines.
	 */
	private void writeIterativeRoundsInfo(File resultDir, List<Integer> rounds, List<Integer> peptidesInRound,
			List<Integer> newPeptides, List<Integer> cumulative, String stopReason, StopRuleResult stopRule) {
		DeepIterativeStopRule.StopRuleResult shared = stopRule == null ? null
				: new DeepIterativeStopRule.StopRuleResult(stopRule.stopRound, stopRule.nRounds, stopRule.kChord,
						stopRule.kPiece, stopRule.kKnee, stopRule.k0, stopRule.tMin);
		DeepIterativeCommonSupport.writeIterativeRoundsInfo(resultDir, rounds, peptidesInRound, newPeptides, cumulative,
				stopReason, shared, CONSECUTIVE_ROUNDS_FOR_PLATEAU, logger, getTaskName());
	}

	/**
	 * Build peptide→proteins map by scanning all round{N}/round{N}.parquet files (q≤0.01, non-decoy).
	 */
	private HashMap<String, HashSet<String>> buildPepToProteinsFromRounds(File markerFile, int lastRound) {
		if (markerFile == null || lastRound <= 0) {
			return new HashMap<>();
		}
		return DiaDeepProteinMappingService.buildPeptideProteinMap(
				DiaDeepProteinMappingService.Request.deepTogether(markerFile, lastRound,
						initialSearchSeparate ? goodDiaFiles : diaFiles,
						diannPar != null ? diannPar.getCut() : "KR"),
				logger, getTaskName());
	}

	/**
	 * Rebuild peptide→protein mapping from round parquets in the report phase, so that
	 * each precursor's ProGroups reflect all proteins from iterative rounds (not just
	 * the single protein from the combined FASTA). Called from exportReport path only.
	 */
	@Override
	protected void enrichPrecursorsWithPeptideProteinMapping(DiaNNPrecursor[] precursors, File parquetFolder) {
		File markerFile = getIterativeResultFolder();
		DiaDeepProteinMappingService.Request request = isDeepPerFileMode()
				? DiaDeepProteinMappingService.Request.deepPerFile(markerFile, goodDiaFiles != null ? goodDiaFiles : diaFiles,
						diannPar != null ? diannPar.getCut() : "KR")
				: DiaDeepProteinMappingService.Request.deepTogether(markerFile, lastRoundForExport,
						goodDiaFiles != null ? goodDiaFiles : diaFiles, diannPar != null ? diannPar.getCut() : "KR");
		DiaDeepProteinMappingService.remapPrecursors(precursors, request, logger, getTaskName());
	}

	/**
	 * Round strategy used by deep marker: one library per round.
	 * FASTA = [proteins from 5 genomes] + [identified peptides from previous round as single-protein entries].
	 * Returns identified peptide count or -1 on failure.
	 */
	private int runRoundAppendPreviousPeptides(int round, File markerFile, List<String> genomes) {
		String[] searchFiles = initialSearchSeparate ? goodDiaFiles : diaFiles;
		return DeepIterativeCommonSupport.runAppendPreviousPeptidesRound(round, markerFile, genomes, searchFiles, true,
				"append strategy", initialSearchSeparate,
				this::writeCombinedFastaFromGenomes,
				this::writeCombinedFastaFromGenomesAndPreviousPeptides,
				(combinedFasta, combinedSpeclib) -> {
					diannPar.setThreads(threadCount);
					diaNNTask.addTask(diannPar, combinedFasta.getAbsolutePath(), combinedSpeclib.getAbsolutePath(), false);
					diaNNTask.run(12);
				},
				(roundParquet, combinedSpeclibPath, files, reanalyse) -> {
					if (initialSearchSeparate) {
						File roundDir = roundParquet.getParentFile();
						int maxConcurrent = DeepIterativeCommonSupport.computeMaxConcurrentTasks();
						for (int i = 0; i < files.length; i++) {
							File pf = new File(roundDir, "per_file_" + i + ".parquet");
						if (pf.exists() && pf.length() > 0) {
							continue;
						}
						DiaNNNoQuantLibrarySearchService.configureNoQuantSearch(
								diannPar, diaLibSearchPar, new String[] { files[i] }, threadCount, reanalyse, false);
						diaNNTask.addNoQuantLibTask(diannPar, diaLibSearchPar, pf.getAbsolutePath(), combinedSpeclibPath);
					}
					diaNNTask.run(maxConcurrent, files.length * 12);
				} else {
					DiaNNNoQuantLibrarySearchService.configureNoQuantSearch(
							diannPar, diaLibSearchPar, files, threadCount, reanalyse, false);
					diaNNTask.addNoQuantLibTask(diannPar, diaLibSearchPar, roundParquet.getAbsolutePath(), combinedSpeclibPath);
					diaNNTask.run(files.length * 12);
					}
				},
				this::countIdentifiedPeptides, logger, getTaskName());
	}

	/**
	 * Write combined FASTA: first genome proteins from the given genomes, then append one entry per unique identified peptide from previousRoundParquet (q≤0.01, non-decoy).
	 * Each peptide is written as ">proteinId\nsequence\n" using the first non-decoy protein from that peptide's proGroups so DIA-NN can link peptide to protein.
	 */
	private boolean writeCombinedFastaFromGenomesAndPreviousPeptides(List<String> genomes, File previousRoundParquet, File outputFasta) {
		File dbFolder = magDb.getOriginalDbFolder();
		if (dbFolder == null || !dbFolder.isDirectory()) {
			logger.error(getTaskName() + ": MAG original DB folder not available");
			return false;
		}
		String cut = diannPar != null && diannPar.getCut() != null ? diannPar.getCut() : "KR";
		DiaNNPrecursor[] precursors;
		if (previousRoundParquet.exists()) {
			DiaNNParquetReader reader = new DiaNNParquetReader(previousRoundParquet.getAbsolutePath(), cut);
			precursors = reader.getDiaNNPrecursors(0.01);
		} else {
			File prevDir = previousRoundParquet.getParentFile();
			String[] filesForRead = initialSearchSeparate ? goodDiaFiles : diaFiles;
			precursors = DeepIterativeCommonSupport.readPrecursorsFromPerFileParquets(
					prevDir, 0.01, filesForRead, cut, logger, getTaskName());
			logger.info(getTaskName() + ": read {} precursors from per-file parquets in {}",
					precursors != null ? precursors.length : 0, prevDir.getAbsolutePath());
		}
		Set<String> addedPeptideSequences = new HashSet<>();
		try (PrintWriter w = new PrintWriter(outputFasta)) {
			for (String genome : genomes) {
				File genomeFaa = new File(dbFolder, genome + ".faa");
				if (!genomeFaa.exists()) genomeFaa = new File(dbFolder, genome + ".fasta");
				if (!genomeFaa.exists()) {
					logger.warn(getTaskName() + ": genome file not found: {}", genome);
					continue;
				}
				try (BufferedReader reader = new BufferedReader(new FileReader(genomeFaa))) {
					String line;
					while ((line = reader.readLine()) != null) w.println(line);
				}
			}
			if (precursors != null) {
				for (DiaNNPrecursor p : precursors) {
					if (p == null || p.getIsDecoy()) continue;
					String seq = p.getSeqString();
					if (seq == null || seq.isEmpty()) continue;
					if (addedPeptideSequences.contains(seq)) continue;
					String[] proGroups = p.getProGroups();
					String proteinId = null;
					if (proGroups != null) {
						for (String pro : proGroups) {
							if (pro == null || pro.isEmpty()) continue;
							if (pro.startsWith(FastaManager.shuffle) || pro.startsWith("REV_") || pro.startsWith("rev_")) continue;
							proteinId = pro;
							break;
						}
					}
					if (proteinId == null) proteinId = "pep_" + seq.length();
					w.println(">" + proteinId);
					w.println(seq);
					addedPeptideSequences.add(seq);
				}
			}
			logger.info(getTaskName() + ": appended {} unique peptides from previous round to combined.fasta", addedPeptideSequences.size());
		} catch (IOException e) {
			logger.error(getTaskName() + ": error writing combined FASTA (genomes + peptides)", e);
			return false;
		}
		return outputFasta.exists() && outputFasta.length() > 0;
	}

	/** Relative gain threshold: stop when New_peptides/Cumulative_peptides < this (e.g. 0.015 = 1.5%). */
	private static final double RELATIVE_GAIN_THRESHOLD = 0.015;
	/** Require this many consecutive rounds below threshold before stopping (trend-based, not single dip). */
	private static final int CONSECUTIVE_ROUNDS_FOR_PLATEAU = 2;

	/** Result of stop-rule evaluation: stop round (0 = don't stop) plus knee indices for iterative_rounds_info. */
	public static class StopRuleResult {
		public final int stopRound;
		public final int nRounds;
		public final int kChord;
		public final int kPiece;
		public final int kKnee;
		public final int k0;
		public final int tMin;

		public StopRuleResult(int stopRound, int nRounds, int kChord, int kPiece, int kKnee, int k0, int tMin) {
			this.stopRound = stopRound;
			this.nRounds = nRounds;
			this.kChord = kChord;
			this.kPiece = kPiece;
			this.kKnee = kKnee;
			this.k0 = k0;
			this.tMin = tMin;
		}
	}

	/**
	 * Iterative DIA stop rule: plateau on cumulative discovery by trend, not a single dip.
	 * Stop only when relative gain (New_peptides/Cumulative_peptides) is below 1.5% for CONSECUTIVE_ROUNDS_FOR_PLATEAU consecutive rounds,
	 * so we do not stop while the trend is still rising (e.g. one round at 4.7% but previous round at 5.3%).
	 * inc is 0-based list of length n (rounds 1..n). Returns result with stopRound in 1..n, or 0 if no plateau detected yet; knee indices set when n>=8.
	 */
	private StopRuleResult computeStopRoundFromIncrements(List<Integer> inc) {
		DeepIterativeStopRule.StopRuleResult shared = DeepIterativeStopRule.computeStopRoundFromIncrements(inc,
				RELATIVE_GAIN_THRESHOLD, CONSECUTIVE_ROUNDS_FOR_PLATEAU,
				DeepIterativeStopRule.AlgorithmProfile.DEEP_MARKER);
		return new StopRuleResult(shared.stopRound, shared.nRounds, shared.kChord, shared.kPiece, shared.kKnee,
				shared.k0, shared.tMin);
	}

	/** Count unique peptide sequences from parquet (q≤0.01, exclude decoys). */
	private int countIdentifiedPeptides(File parquetFile) {
		String cut = diannPar != null && diannPar.getCut() != null ? diannPar.getCut() : "KR";
		return DeepIterativeCommonSupport.countIdentifiedPeptides(parquetFile, cut, logger, getTaskName());
	}

	/**
	 * Check whether the marker folder already has per-file parquets from a previous separate search.
	 */
	private boolean hasExistingPerFileParquets(File markerDir) {
		if (diaFiles == null || diaFiles.length == 0 || markerDir == null) {
			return false;
		}
		for (int i = 0; i < diaFiles.length; i++) {
			File pf = new File(markerDir, "per_file_" + i + ".parquet");
			if (pf.exists() && pf.length() > 0) {
				return true;
			}
		}
		return false;
	}

	/**
	 * Collect the subset of diaFiles whose per-file parquets exist (non-damaged files).
	 */
	private String[] collectGoodDiaFilesFromPerFileParquets(File markerDir) {
		java.util.ArrayList<String> goodList = new java.util.ArrayList<>();
		if (diaFiles != null && markerDir != null) {
			for (int i = 0; i < diaFiles.length; i++) {
				File pf = new File(markerDir, "per_file_" + i + ".parquet");
				if (pf.exists() && pf.length() > 0) {
					goodList.add(diaFiles[i]);
				}
			}
		}
		return goodList.toArray(new String[0]);
	}

	/** Write combined FASTA from genome IDs using MAG DB original folder (genome.faa or genome.fasta). */
	private boolean writeCombinedFastaFromGenomes(List<String> genomes, File outputFasta) {
		File dbFolder = magDb.getOriginalDbFolder();
		if (dbFolder == null || !dbFolder.isDirectory()) {
			logger.error(getTaskName() + ": MAG original DB folder not available");
			return false;
		}
		try (PrintWriter w = new PrintWriter(outputFasta)) {
			for (String genome : genomes) {
				File genomeFaa = new File(dbFolder, genome + ".faa");
				if (!genomeFaa.exists()) {
					genomeFaa = new File(dbFolder, genome + ".fasta");
				}
				if (!genomeFaa.exists()) {
					logger.warn(getTaskName() + ": genome file not found: {}", genome);
					continue;
				}
				try (BufferedReader reader = new BufferedReader(new FileReader(genomeFaa))) {
					String line;
					while ((line = reader.readLine()) != null) {
						w.println(line);
					}
				}
			}
		} catch (IOException e) {
			logger.error(getTaskName() + ": error writing combined FASTA", e);
			return false;
		}
		return outputFasta.exists() && outputFasta.length() > 0;
	}
}
