package uk.ac.quadram.metapilot.core.runner;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicInteger;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import uk.ac.quadram.metapilot.core.dia.DiaDeepProteinMappingService;
import uk.ac.quadram.metapilot.core.diann.DiaNNHeadlessRunner;
import uk.ac.quadram.metapilot.core.mag.MagGenomeSelectionService;
import uk.ac.quadram.metapilot.core.report.MetaPilotHtmlReportService;
import uk.ac.quadram.metapilot.core.resources.MetaPilotResourceConfig;
import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;
import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNParquetReader;
import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNPrecursor;
import uk.ac.quadram.metapilot.dbSearch.diann.DiannParameter;
import uk.ac.quadram.metapilot.task.dia.DeepComprehensivePaths;
import uk.ac.quadram.metapilot.task.dia.DeepIterativeCommonSupport;
import uk.ac.quadram.metapilot.task.dia.DeepIterativeStopRule;
import uk.ac.quadram.metapilot.task.dia.DiaLibSearchPar;
import uk.ac.quadram.metapilot.task.dia.par.MetaParameterDiaMarker;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;

/**
 * Headless DEEP_PER_FILE DIA Marker workflow runner (HPC version).
 *
 * <p>Mirrors the behavior of {@code DiaMarkerDeepPerFileTask} when {@code markerMode == DEEP_PER_FILE},
 * but without any JavaFX dependency. For each raw DIA file:
 * <ol>
 *   <li>run a marker no-quant search ({@code marker.parquet}) using only that raw file,</li>
 *   <li>rank/select genomes per file ({@code genome_evidence/ranked/selected.tsv}),</li>
 *   <li>run iterative rounds (one library per round = 5 genomes + previous-round peptides),</li>
 *   <li>stop on the same trend-based plateau / knee+plateau rule as the GUI workflow.</li>
 * </ol>
 *
 * <p>After all per-file workers finish, the runner builds a single global combined
 * library from the peptides identified across <i>all</i> per-file rounds (for low
 * missingness) and runs one final quantifying DIA-NN search across all DIA files.
 *
 * <p>Per-file parallelism follows the DIA-NN RAM-budget rule:
 * {@code poolSize = totalRAM_GB / 30} (sequential when total RAM &lt; 30 GB);
 * each concurrent worker is given {@code totalThreads / poolSize} DIA-NN threads.
 *
 * <p>Outputs (under {@code <resultDir>/marker_deep/}):
 * <pre>
 *   marker/&lt;rawname&gt;/marker.parquet
 *   marker/&lt;rawname&gt;/genome_evidence.tsv, genome_ranked.tsv, genome_selected.tsv
 *   marker/&lt;rawname&gt;/round1/{combined.fasta, combined.speclib, round1.parquet}
 *   marker/&lt;rawname&gt;/roundN/...
 *   combined/{combined.fasta, combined.speclib,
 *             combined.parquet, combined.pr_matrix.tsv, combined.pg_matrix.tsv}
 * </pre>
 *
 * <p>This runner uses the shared final HTML report service after post-analysis.
 */
public final class DiaMarkerDeepPerFileHeadlessRunner {

    private static final Logger LOGGER = LogManager.getLogger(DiaMarkerDeepPerFileHeadlessRunner.class);

    private static final String TASK_NAME = "DIA Deep Marker (per-file, HPC)";
    private static final int DEFAULT_GENOMES_PER_ROUND = 5;
    private static final int MAX_ITERATIVE_ROUNDS = 100;
    private static final double RELATIVE_GAIN_THRESHOLD = 0.015;
    private static final int CONSECUTIVE_ROUNDS_FOR_PLATEAU = 2;

    // Overall progress is reported on the same 0-100 scale the GUI task uses around this
    // runner (updateProgress(2,100) right before calling in, updateProgress(97,100) right
    // after it returns), so these bands sit inside [2, 97] and give real, monotonically
    // increasing feedback instead of the bar sitting frozen for the whole run.
    private static final double PROGRESS_SETUP_END = 5.0;
    private static final double PROGRESS_INITIAL_EVIDENCE_START = 5.0;
    private static final double PROGRESS_INITIAL_EVIDENCE_END = 10.0;
    private static final double PROGRESS_PERFILE_START = 10.0;
    private static final double PROGRESS_PERFILE_END = 85.0;
    private static final double PROGRESS_COMBINED_START = 85.0;
    private static final double PROGRESS_COMBINED_FASTA = 88.0;
    private static final double PROGRESS_COMBINED_SPECLIB = 91.0;
    private static final double PROGRESS_COMBINED_END = 95.0;
    // Kept below the GUI's own updateProgress(97,100) checkpoint (set right after this
    // runner returns, for "generating combined report") so progress stays monotonic.
    private static final double PROGRESS_FINAL = 96.0;

    private DiaMarkerDeepPerFileHeadlessRunner() {
    }

    /** Report progress as a 0-100 percentage (matches the surrounding GUI updateProgress(x,100) calls). */
    private static void reportProgress(ProgressCallback cb, double pct) {
        if (cb == null) {
            return;
        }
        cb.onProgress(Math.max(0.0, Math.min(100.0, pct)), 100.0);
    }

    /** Report progress for {@code completed}/{@code total} discrete units of work within [bandStart, bandEnd]. */
    private static void reportBandProgress(ProgressCallback cb, int completed, int total, double bandStart,
            double bandEnd) {
        if (cb == null || total <= 0) {
            return;
        }
        double fraction = Math.max(0.0, Math.min(1.0, (double) completed / (double) total));
        reportProgress(cb, bandStart + (bandEnd - bandStart) * fraction);
    }

    public static final class Result {
        public final File resultDir;
        public final File markerDeepDir;
        public final File combinedDir;
        public final File combinedFasta;
        public final File combinedSpeclib;
        public final File combinedParquet;
        public final File prMatrixTsv;
        public final File pgMatrixTsv;

        Result(File resultDir, File markerDeepDir, File combinedDir, File combinedFasta, File combinedSpeclib,
               File combinedParquet, File prMatrixTsv, File pgMatrixTsv) {
            this.resultDir = resultDir;
            this.markerDeepDir = markerDeepDir;
            this.combinedDir = combinedDir;
            this.combinedFasta = combinedFasta;
            this.combinedSpeclib = combinedSpeclib;
            this.combinedParquet = combinedParquet;
            this.prMatrixTsv = prMatrixTsv;
            this.pgMatrixTsv = pgMatrixTsv;
        }
    }

    public static Result run(MetaParameterDiaMarker params, File diannExe, File magDbConfigJson,
                             ProgressCallback cb) throws Exception {
        return run(params, diannExe, magDbConfigJson, cb, null, null);
    }

    public static Result run(MetaParameterDiaMarker params, File diannExe, File magDbConfigJson,
                             ProgressCallback cb, File condaPrefix) throws Exception {
        return run(params, diannExe, magDbConfigJson, cb, condaPrefix, null);
    }

    public static Result run(MetaParameterDiaMarker params, File diannExe, File magDbConfigJson,
                             ProgressCallback cb, File condaPrefix, File funcDbFile) throws Exception {
        return run(params, diannExe, magDbConfigJson, cb, condaPrefix, funcDbFile, true, null);
    }

    public static Result run(MetaParameterDiaMarker params, File diannExe, File magDbConfigJson,
                             ProgressCallback cb, File condaPrefix, File funcDbFile,
                             boolean generateHtmlReport, String reportPythonHint) throws Exception {
        return run(params, diannExe, magDbConfigJson, cb, condaPrefix, funcDbFile, generateHtmlReport,
                reportPythonHint, null);
    }

    public static Result run(MetaParameterDiaMarker params, File diannExe, File magDbConfigJson,
                             ProgressCallback cb, File condaPrefix, File funcDbFile,
                             boolean generateHtmlReport, String reportPythonHint,
                             MetaPilotResourceConfig reportResourceConfig) throws Exception {
        return run(params, diannExe, magDbConfigJson, cb, condaPrefix, funcDbFile, generateHtmlReport,
                reportPythonHint, reportResourceConfig, false);
    }

    /**
     * Same as the 9-argument {@link #run}, with an additional opt-in
     * {@code exportMetaAtlasParquet} flag (see {@code --export-metaatlas} in the CLI). This
     * is a brand new overload, not a change to any existing signature, so it cannot affect
     * any existing caller.
     */
    public static Result run(MetaParameterDiaMarker params, File diannExe, File magDbConfigJson,
                             ProgressCallback cb, File condaPrefix, File funcDbFile,
                             boolean generateHtmlReport, String reportPythonHint,
                             MetaPilotResourceConfig reportResourceConfig, boolean exportMetaAtlasParquet)
            throws Exception {
        return runInternal(params, diannExe, magDbConfigJson, cb, condaPrefix, funcDbFile, true,
                generateHtmlReport, reportPythonHint, reportResourceConfig, exportMetaAtlasParquet);
    }

    public static Result runWithoutPostAnalysis(MetaParameterDiaMarker params, File diannExe, File magDbConfigJson,
                                                ProgressCallback cb, File condaPrefix, File funcDbFile) throws Exception {
        return runInternal(params, diannExe, magDbConfigJson, cb, condaPrefix, funcDbFile, false,
                false, null, null, false);
    }

    private static Result runInternal(MetaParameterDiaMarker params, File diannExe, File magDbConfigJson,
                                      ProgressCallback cb, File condaPrefix, File funcDbFile,
                                      boolean exportPostAnalysis, boolean generateHtmlReport,
                                      String reportPythonHint, MetaPilotResourceConfig reportResourceConfig,
                                      boolean exportMetaAtlasParquet)
            throws Exception {
        if (cb == null) cb = new ProgressCallback.NoOp();
        if (params == null) throw new IllegalArgumentException("params is null");
        if (diannExe == null || !diannExe.exists()) {
            throw new IOException("DIA-NN not found: " + (diannExe == null ? "null" : diannExe.getAbsolutePath()));
        }

        File resultDir = new File(params.getResult());
        if (!resultDir.exists()) resultDir.mkdirs();
        if (!resultDir.isDirectory()) {
            throw new IOException("Result dir is not a directory: " + resultDir.getAbsolutePath());
        }

        String[] diaFiles = params.getMetadata() != null ? params.getMetadata().getRawFiles() : null;
        if (diaFiles == null || diaFiles.length == 0) {
            throw new IOException("No DIA files found in params metadata (input_output.experiments)");
        }

        MagDbItem magDb = DiaMarkerHeadlessRunner.resolveMagDb(params, magDbConfigJson);
        if (magDb == null) throw new IOException("MAG DB could not be resolved from params + magdb config");

        // Ensure Marker4 FASTA, genomes-taxa.tsv and Marker4.speclib exist (shared across all per-file workers).
        DiaMarkerHeadlessRunner.ensureMarker4Files(magDb, diannExe, params, cb);
        cb.onMessage(TASK_NAME + ": preparing Marker4 library");
        DiaMarkerHeadlessRunner.ensureMarker4Library(magDb, diannExe, params, cb, condaPrefix);
        reportProgress(cb, PROGRESS_SETUP_END);

        // Per-file folder layout:  <resultDir>/marker_deep/marker/<rawname>/...
        File markerDeepDir = new File(resultDir, "marker_deep");
        if (!markerDeepDir.exists()) markerDeepDir.mkdirs();
        File markerRoot = new File(markerDeepDir, "marker");
        if (!markerRoot.exists()) markerRoot.mkdirs();

        LinkedHashMap<String, File> perFileDirs = new LinkedHashMap<>();
        LinkedHashMap<String, String> perFileNames = new LinkedHashMap<>();
        for (int i = 0; i < diaFiles.length; i++) {
            String raw = diaFiles[i];
            String rawName = DeepComprehensivePaths.runIdFromRawFilePath(raw);
            File perFileDir = new File(markerRoot, rawName);
            perFileDir.mkdirs();
            perFileDirs.put(raw, perFileDir);
            perFileNames.put(raw, rawName);
        }

        int totalThreads = Math.max(1, params.getThreadCount());
        DeepIterativeCommonSupport.PoolDecision decision =
                DeepIterativeCommonSupport.computeDeepPerFilePoolSizeHpc(diaFiles.length, totalThreads);
        int poolSize = decision.poolSize;
        int perTaskThreads = decision.perTaskThreads;
        String decisionMsg = String.format(
                "per-file processing budget -> maxConcurrentTasks=%d, perTaskThreads=%d, totalThreads=%d, files=%d, "
                        + "rule=%s, memorySource=%s (%s), effectiveMemoryGB=%d, ramPerTaskGB=%d, explicitOverride=%s",
                poolSize, perTaskThreads, totalThreads, diaFiles.length,
                "min(files, memoryGB/ramPerTaskGB, totalThreads/" + DeepIterativeCommonSupport.MIN_THREADS_PER_TASK + ")",
                decision.memorySource, decision.memorySourceDetail,
                decision.effectiveMemoryGB, decision.ramPerTaskGB, decision.explicitOverride);
        LOGGER.info("{}: {}", TASK_NAME, decisionMsg);
        cb.onMessage(TASK_NAME + ": " + decisionMsg);
        cb.onMessage(TASK_NAME + ": (override pool with METAPILOT_DEEP_PER_FILE_POOL=N or"
                + " ram-per-task with METAPILOT_DEEP_PER_FILE_RAM_PER_TASK_GB=N)");

        ConcurrentHashMap<String, Integer> perFileLastRound = new ConcurrentHashMap<>();

        boolean representativeClustering = params.isRepresentativeClusteringEnabled();
        String[] filesForIterativeRounds = diaFiles;
        if (representativeClustering) {
            cb.onMessage(TASK_NAME + ": representative clustering enabled (preset="
                    + params.getRepresentativeClusteringPreset() + ")");
            boolean initialOk;
            if (poolSize <= 1) {
                initialOk = runInitialEvidenceSequential(diaFiles, perFileDirs, perFileNames,
                        diannExe, params, magDb, condaPrefix, perTaskThreads, cb);
            } else {
                initialOk = runInitialEvidenceParallel(diaFiles, perFileDirs, perFileNames,
                        diannExe, params, magDb, condaPrefix, poolSize, perTaskThreads, cb);
            }
            if (!initialOk) {
                if (params.isRepresentativeClusteringDryRun()) {
                    LOGGER.warn("{}: representative clustering initial evidence failed during grouping dry run",
                            TASK_NAME);
                    throw new IOException(TASK_NAME + ": grouping dry run failed before clustering reports could be written");
                }
                LOGGER.warn("{}: representative clustering initial evidence failed; falling back to full per-file mode",
                        TASK_NAME);
                cb.onMessage("WARNING: representative clustering failed before grouping; falling back to full per-file mode");
                representativeClustering = false;
            } else {
                RepresentativePeptideClusterer.Result grouping = buildRepresentativeGrouping(diaFiles, perFileDirs,
                        perFileNames, params);
                writeRepresentativeGroupsTsv(new File(markerDeepDir, "representative_groups.tsv"), grouping);
                writeRepresentativePairwiseTsv(new File(markerDeepDir, "representative_pairwise_similarity.tsv"),
                        grouping, params);
                writeRepresentativeClusteringSummaryTsv(new File(markerDeepDir,
                        "representative_clustering_summary.tsv"), grouping, params);
                if (params.isRepresentativeClusteringDryRun()) {
                    cb.onMessage(TASK_NAME + ": grouping dry run completed; reports written under "
                            + markerDeepDir.getAbsolutePath());
                    return new Result(resultDir, markerDeepDir, new File(markerDeepDir, "combined"),
                            null, null, null, null, null);
                }
                if (grouping.representativeFiles.isEmpty()) {
                    LOGGER.warn("{}: representative clustering produced no representatives; falling back to full per-file mode",
                            TASK_NAME);
                    cb.onMessage("WARNING: representative clustering produced no representatives; falling back to full per-file mode");
                    representativeClustering = false;
                } else {
                    filesForIterativeRounds = grouping.representativeFiles.toArray(new String[0]);
                    cb.onMessage(TASK_NAME + ": representative clustering selected "
                            + filesForIterativeRounds.length + " representatives from " + diaFiles.length
                            + " raw files");
                }
            }
        }

        boolean ok;
        if (poolSize <= 1 || filesForIterativeRounds.length <= 1) {
            ok = runPerFileSequential(filesForIterativeRounds, perFileDirs, perFileNames, perFileLastRound,
                    diannExe, params, magDb, condaPrefix, perTaskThreads, cb);
        } else {
            int iterativePoolSize = Math.min(poolSize, filesForIterativeRounds.length);
            ok = runPerFileParallel(filesForIterativeRounds, perFileDirs, perFileNames, perFileLastRound,
                    diannExe, params, magDb, condaPrefix, iterativePoolSize, perTaskThreads, cb);
        }
        if (!ok) {
            throw new IOException(TASK_NAME + ": one or more per-file workers failed");
        }

        // Order perFileLastRound by original DIA file order for deterministic downstream processing.
        LinkedHashMap<String, Integer> orderedLastRound = new LinkedHashMap<>();
        for (String raw : diaFiles) {
            Integer v = perFileLastRound.get(raw);
            orderedLastRound.put(raw, v == null ? 0 : v);
        }

        // Build global combined.fasta / combined.speclib from all per-file round parquets and run final search.
        cb.onMessage(TASK_NAME + ": building global combined search from all per-file rounds");
        File combinedDir = new File(markerDeepDir, "combined");
        if (!combinedDir.exists()) combinedDir.mkdirs();
        File combinedFasta = new File(combinedDir, "combined.fasta");
        File combinedSpeclib = new File(combinedDir, "combined.speclib");
        File combinedParquet = new File(combinedDir, "combined.parquet");
        File prMatrix = new File(combinedDir, "combined.pr_matrix.tsv");
        File pgMatrix = new File(combinedDir, "combined.pg_matrix.tsv");

        if (!runGlobalCombinedFromPerFileRounds(diaFiles, perFileDirs, orderedLastRound,
                combinedDir, combinedFasta, combinedSpeclib, combinedParquet,
                diannExe, params, condaPrefix, representativeClustering, cb)) {
            throw new IOException(TASK_NAME + ": global combined search failed");
        }

        if (!combinedParquet.exists()) throw new IOException("Combined quant did not produce: " + combinedParquet.getAbsolutePath());
        if (!prMatrix.exists()) throw new IOException("Combined quant did not produce: " + prMatrix.getAbsolutePath());
        if (!pgMatrix.exists()) throw new IOException("Combined quant did not produce: " + pgMatrix.getAbsolutePath());

        DiaMarkerHeadlessRunner.runHostSearchIfRequested(markerDeepDir, diannExe, params,
                params.getMetadata() != null ? params.getMetadata().getRawFiles() : null, cb, condaPrefix);

        // TSV/report exports go under marker_deep/, matching the Windows task layout.
        if (exportPostAnalysis) {
            try {
                File gtdbNcbi = null;
                File maybe = new File("Resources", "gtdb_ncbi");
                if (maybe.exists()) gtdbNcbi = maybe;
                DiaMagPostAnalysisTsvExporter.export(
                        markerDeepDir,
                        magDb,
                        combinedParquet,
                        prMatrix,
                        gtdbNcbi,
                        params.getMarkerMode(),
                        funcDbFile,
                        params.getMicroDb(),
                        params.getQuantitativeMode(),
                        cb,
                        DiaDeepProteinMappingService.Request.deepPerFile(markerDeepDir, diaFiles, "KR"),
                        exportMetaAtlasParquet
                );
                runHtmlReportIfRequested(markerDeepDir, params, generateHtmlReport, reportPythonHint,
                        reportResourceConfig, cb);
            } catch (Exception e) {
                cb.onMessage("WARNING: DIA marker post-analysis export failed: " + e.getMessage());
                throw e;
            }
        }

        reportProgress(cb, PROGRESS_FINAL);
        cb.onComplete("DIA Deep Marker (per-file) complete: " + resultDir.getAbsolutePath());
        return new Result(resultDir, markerDeepDir, combinedDir, combinedFasta, combinedSpeclib,
                combinedParquet, prMatrix, pgMatrix);
    }

    private static void runHtmlReportIfRequested(File resultDir, MetaParameterDiaMarker params,
                                                 boolean generateHtmlReport, String reportPythonHint,
                                                 MetaPilotResourceConfig reportResourceConfig,
                                                 ProgressCallback cb) {
        if (!generateHtmlReport) {
            cb.onMessage("HTML report generation skipped (disabled)");
            return;
        }
        cb.onMessage(TASK_NAME + ": generating shared HTML report");
        MetaPilotHtmlReportService.Result report = MetaPilotHtmlReportService.runStandardReportWarningOnly(
                resultDir,
                "DIA_MAG",
                params != null ? params.getThreadCount() : 1,
                reportPythonHint,
                reportResourceConfig,
                cb
        );
        if (!report.warnings.isEmpty()) {
            cb.onMessage("WARNING: HTML report generation completed with warnings: "
                    + String.join("; ", report.warnings));
        }
    }

    // ----------------------- per-file orchestration -----------------------

    private static boolean runPerFileSequential(String[] diaFiles,
            LinkedHashMap<String, File> perFileDirs, LinkedHashMap<String, String> perFileNames,
            ConcurrentHashMap<String, Integer> perFileLastRound, File diannExe, MetaParameterDiaMarker params,
            MagDbItem magDb, File condaPrefix, int perTaskThreads, ProgressCallback cb) {
        int total = diaFiles.length;
        for (int i = 0; i < total; i++) {
            String raw = diaFiles[i];
            File perFileDir = perFileDirs.get(raw);
            String rawName = perFileNames.get(raw);
            cb.onMessage(TASK_NAME + ": per-file workflow (" + (i + 1) + "/" + total + "): " + rawName);
            int lastRound = runOnePerFileWorker(raw, perFileDir, diannExe, params, magDb, condaPrefix, perTaskThreads, cb);
            if (lastRound <= 0) {
                LOGGER.error("{}: per-file workflow failed for {}", TASK_NAME, raw);
                return false;
            }
            perFileLastRound.put(raw, lastRound);
            reportBandProgress(cb, i + 1, total, PROGRESS_PERFILE_START, PROGRESS_PERFILE_END);
        }
        return true;
    }

    private static boolean runPerFileParallel(String[] diaFiles,
            LinkedHashMap<String, File> perFileDirs, LinkedHashMap<String, String> perFileNames,
            ConcurrentHashMap<String, Integer> perFileLastRound, File diannExe, MetaParameterDiaMarker params,
            MagDbItem magDb, File condaPrefix, int poolSize, int perTaskThreads, ProgressCallback cb) {
        int total = diaFiles.length;
        LOGGER.info("{}: running {} per-file workflows in parallel (poolSize={}, perTaskThreads={})",
                TASK_NAME, total, poolSize, perTaskThreads);

        ExecutorService executor = Executors.newFixedThreadPool(poolSize, r -> {
            Thread t = new Thread(r, "DeepPerFile-Worker");
            t.setDaemon(true);
            return t;
        });
        AtomicBoolean anyFailed = new AtomicBoolean(false);
        AtomicInteger completed = new AtomicInteger(0);
        List<Future<?>> futures = new ArrayList<>();

        try {
            for (String raw : diaFiles) {
                final File perFileDir = perFileDirs.get(raw);
                final String rawName = perFileNames.get(raw);
                final String rawFile = raw;
                futures.add(executor.submit(() -> {
                    try {
                        if (anyFailed.get()) return;
                        LOGGER.info("{}: [parallel] starting per-file workflow for {} (perTaskThreads={})",
                                TASK_NAME, rawName, perTaskThreads);
                        int lastRound = runOnePerFileWorker(rawFile, perFileDir, diannExe, params, magDb,
                                condaPrefix, perTaskThreads, cb);
                        if (lastRound <= 0) {
                            LOGGER.error("{}: [parallel] per-file workflow failed for {}", TASK_NAME, rawFile);
                            anyFailed.set(true);
                            return;
                        }
                        perFileLastRound.put(rawFile, lastRound);
                        int done = completed.incrementAndGet();
                        cb.onMessage(TASK_NAME + ": [parallel] " + done + "/" + total + " complete (latest: "
                                + rawName + ")");
                        reportBandProgress(cb, done, total, PROGRESS_PERFILE_START, PROGRESS_PERFILE_END);
                        LOGGER.info("{}: [parallel] completed per-file workflow for {} (lastRound={})",
                                TASK_NAME, rawName, lastRound);
                    } catch (Throwable t) {
                        LOGGER.error("{}: [parallel] exception in per-file worker for {}", TASK_NAME, rawFile, t);
                        anyFailed.set(true);
                    }
                }));
            }
            for (Future<?> f : futures) {
                try {
                    f.get();
                } catch (Exception e) {
                    LOGGER.error("{}: [parallel] worker future failed", TASK_NAME, e);
                    anyFailed.set(true);
                }
            }
        } finally {
            executor.shutdown();
            try {
                if (!executor.awaitTermination(5, TimeUnit.MINUTES)) {
                    executor.shutdownNow();
                }
            } catch (InterruptedException ignored) {
                executor.shutdownNow();
                Thread.currentThread().interrupt();
            }
        }

        return !anyFailed.get();
    }

    private static boolean runInitialEvidenceSequential(String[] diaFiles,
            LinkedHashMap<String, File> perFileDirs, LinkedHashMap<String, String> perFileNames,
            File diannExe, MetaParameterDiaMarker params, MagDbItem magDb, File condaPrefix,
            int perTaskThreads, ProgressCallback cb) {
        int total = diaFiles.length;
        for (int i = 0; i < total; i++) {
            String raw = diaFiles[i];
            File perFileDir = perFileDirs.get(raw);
            String rawName = perFileNames.get(raw);
            cb.onMessage(TASK_NAME + ": initial marker evidence (" + (i + 1) + "/" + total + "): " + rawName);
            if (!prepareInitialEvidenceForOneFile(raw, perFileDir, diannExe, params, magDb,
                    condaPrefix, perTaskThreads, cb)) {
                LOGGER.error("{}: initial evidence failed for {}", TASK_NAME, raw);
                return false;
            }
            reportBandProgress(cb, i + 1, total, PROGRESS_INITIAL_EVIDENCE_START, PROGRESS_INITIAL_EVIDENCE_END);
        }
        return true;
    }

    private static boolean runInitialEvidenceParallel(String[] diaFiles,
            LinkedHashMap<String, File> perFileDirs, LinkedHashMap<String, String> perFileNames,
            File diannExe, MetaParameterDiaMarker params, MagDbItem magDb, File condaPrefix,
            int poolSize, int perTaskThreads, ProgressCallback cb) {
        LOGGER.info("{}: running initial marker evidence for {} files in parallel (poolSize={}, perTaskThreads={})",
                TASK_NAME, diaFiles.length, poolSize, perTaskThreads);
        ExecutorService executor = Executors.newFixedThreadPool(poolSize, r -> {
            Thread t = new Thread(r, "DeepPerFile-Initial");
            t.setDaemon(true);
            return t;
        });
        AtomicBoolean anyFailed = new AtomicBoolean(false);
        AtomicInteger completed = new AtomicInteger(0);
        int total = diaFiles.length;
        List<Future<?>> futures = new ArrayList<>();
        try {
            for (String raw : diaFiles) {
                final File perFileDir = perFileDirs.get(raw);
                final String rawName = perFileNames.get(raw);
                final String rawFile = raw;
                futures.add(executor.submit(() -> {
                    try {
                        if (anyFailed.get()) return;
                        cb.onMessage(TASK_NAME + ": initial marker evidence for " + rawName);
                        if (!prepareInitialEvidenceForOneFile(rawFile, perFileDir, diannExe, params, magDb,
                                condaPrefix, perTaskThreads, cb)) {
                            anyFailed.set(true);
                            return;
                        }
                        reportBandProgress(cb, completed.incrementAndGet(), total,
                                PROGRESS_INITIAL_EVIDENCE_START, PROGRESS_INITIAL_EVIDENCE_END);
                    } catch (Throwable t) {
                        LOGGER.error("{}: initial evidence exception for {}", TASK_NAME, rawFile, t);
                        anyFailed.set(true);
                    }
                }));
            }
            for (Future<?> f : futures) {
                try {
                    f.get();
                } catch (Exception e) {
                    LOGGER.error("{}: initial evidence worker future failed", TASK_NAME, e);
                    anyFailed.set(true);
                }
            }
        } finally {
            executor.shutdown();
            try {
                if (!executor.awaitTermination(5, TimeUnit.MINUTES)) {
                    executor.shutdownNow();
                }
            } catch (InterruptedException ignored) {
                executor.shutdownNow();
                Thread.currentThread().interrupt();
            }
        }
        return !anyFailed.get();
    }

    /**
     * Run the iterative phase (marker search + iterative rounds) for a single raw file.
     * Returns lastRound (>0) on success, <=0 on failure.
     */
    private static int runOnePerFileWorker(String raw, File perFileDir, File diannExe,
                                           MetaParameterDiaMarker params, MagDbItem magDb, File condaPrefix,
                                           int perTaskThreads, ProgressCallback cb) {
        try {
            if (!prepareInitialEvidenceForOneFile(raw, perFileDir, diannExe, params, magDb,
                    condaPrefix, perTaskThreads, cb)) {
                return -1;
            }
            return runRoundsForPreparedFile(raw, perFileDir, diannExe, params, magDb, condaPrefix,
                    perTaskThreads, cb);
        } catch (Throwable t) {
            LOGGER.error("{}: runOnePerFileWorker exception for {}", TASK_NAME, raw, t);
            return -1;
        }
    }

    private static boolean prepareInitialEvidenceForOneFile(String raw, File perFileDir, File diannExe,
                                                           MetaParameterDiaMarker params, MagDbItem magDb,
                                                           File condaPrefix, int perTaskThreads,
                                                           ProgressCallback cb) {
        try {
            File markerParquet = new File(perFileDir, "marker.parquet");
            File markerLib = magDb.resolveMarker4LibraryFile();
            if (markerLib == null || !markerLib.exists()) {
                LOGGER.error("{}: Marker4 library not found for per-file marker search: {}",
                        TASK_NAME, markerLib == null ? "null" : markerLib.getAbsolutePath());
                return false;
            }
            if (!markerParquet.exists()) {
                cb.onMessage(TASK_NAME + ": marker no-quant search for " + perFileDir.getName());
                DiaNNNoQuantLibrarySearchService.run(diannExe, params, new String[] { raw }, markerLib,
                        markerParquet, perTaskThreads, true, false, cb, condaPrefix, LOGGER, TASK_NAME);
            }
            if (!markerParquet.exists()) {
                LOGGER.error("{}: marker.parquet not produced for {}", TASK_NAME, perFileDir.getName());
                return false;
            }

            File evidenceTsv = new File(perFileDir, "genome_evidence.tsv");
            File rankedTsv = new File(perFileDir, "genome_ranked.tsv");
            File selectedTsv = new File(perFileDir, "genome_selected.tsv");
            if (!rankedTsv.exists() || !selectedTsv.exists()) {
                cb.onMessage(TASK_NAME + ": genome selection for " + perFileDir.getName());
                DiaNNPrecursor[] precursors = DiaMarkerHeadlessRunner.readPrecursors(markerParquet, params);
                if (precursors == null || precursors.length == 0) {
                    LOGGER.error("{}: no precursors in marker.parquet for {}", TASK_NAME, perFileDir.getName());
                    return false;
                }
                File currentFile = magDb.getCurrentFile();
                File markerFasta = new File(currentFile, "Marker4.fasta");
                File genomesTaxaTsv = new File(currentFile, "genomes-taxa.tsv");
                MagGenomeSelectionService.Request selectionRequest = new MagGenomeSelectionService.Request();
                selectionRequest.magDb = magDb;
                selectionRequest.genomeSelectorMode = MagDbItem.MARKER4;
                selectionRequest.genomeSelectionStrategy = MagGenomeSelectionService.SELECTION_STRATEGY_COMPREHENSIVE;
                selectionRequest.diaPrecursors = precursors;
                selectionRequest.markerFasta = markerFasta;
                selectionRequest.genomesTaxaTsv = genomesTaxaTsv;
                selectionRequest.evidenceFile = evidenceTsv;
                selectionRequest.rankedFile = rankedTsv;
                selectionRequest.selectedFile = selectedTsv;
                selectionRequest.maximumSelected = 300;
                MagGenomeSelectionService.selectGenomes(selectionRequest, cb);
            }
            if (!rankedTsv.exists()) {
                LOGGER.error("{}: genome_ranked.tsv missing for {}", TASK_NAME, perFileDir.getName());
                return false;
            }
            return true;
        } catch (Throwable t) {
            LOGGER.error("{}: prepareInitialEvidenceForOneFile exception for {}", TASK_NAME, raw, t);
            return false;
        }
    }

    private static int runRoundsForPreparedFile(String raw, File perFileDir, File diannExe,
                                                MetaParameterDiaMarker params, MagDbItem magDb, File condaPrefix,
                                                int perTaskThreads, ProgressCallback cb) {
        File rankedTsv = new File(perFileDir, "genome_ranked.tsv");
        if (!rankedTsv.exists()) {
            LOGGER.error("{}: genome_ranked.tsv missing for {}", TASK_NAME, perFileDir.getName());
            return -1;
        }
        int genomesPerRound = params.getGenomesPerRound() > 0 ? params.getGenomesPerRound() : DEFAULT_GENOMES_PER_ROUND;
        String cut = "KR";
        String[] singleRawArr = new String[] { raw };
        DeepIterativeCommonSupport.IterativeRoundLoopResult loop = DeepIterativeCommonSupport.runIterativeRoundLoop(
                rankedTsv, genomesPerRound, MAX_ITERATIVE_ROUNDS, "Round", cut,
                cb::onMessage,
                (round, genomes) -> runRoundAppendPreviousPeptides(round, perFileDir, genomes, singleRawArr,
                        magDb, diannExe, params, condaPrefix, perTaskThreads, cb),
                DiaMarkerDeepPerFileHeadlessRunner::computeStopRule,
                LOGGER, TASK_NAME);

        if (loop.lastRound <= 0) {
            LOGGER.error("{}: no round completed for {}", TASK_NAME, perFileDir.getName());
            return -1;
        }
        DeepIterativeCommonSupport.writeIterativeRoundsInfo(perFileDir, loop.roundsRecord, loop.peptidesInRoundRecord,
                loop.incrementsPerRound, loop.cumulativeRecord, loop.stopReason, loop.lastStopRuleResult,
                CONSECUTIVE_ROUNDS_FOR_PLATEAU, LOGGER, TASK_NAME);
        return loop.lastRound;
    }

    // ----------------------- representative clustering -----------------------

    private static RepresentativePeptideClusterer.Result buildRepresentativeGrouping(String[] diaFiles,
            LinkedHashMap<String, File> perFileDirs, LinkedHashMap<String, String> perFileNames,
            MetaParameterDiaMarker params) {
        List<RepresentativePeptideClusterer.Input> signatures = new ArrayList<>();
        for (int i = 0; i < diaFiles.length; i++) {
            String raw = diaFiles[i];
            File dir = perFileDirs.get(raw);
            String rawName = perFileNames.get(raw);
            RepresentativePeptideClusterer.Input signature = readMarkerSignature(signatures.size(), raw, rawName,
                    dir, params);
            if (signature != null) {
                signatures.add(signature);
            }
        }
        int minPeptides = Math.max(0, params.getRepresentativeMinMarkerPeptides());
        double minJaccard = params.getRepresentativePeptideJaccardThreshold();
        double minContainment = params.getRepresentativePeptideContainmentThreshold();
        RepresentativePeptideClusterer.Result result = RepresentativePeptideClusterer.cluster(signatures,
                minPeptides, minJaccard, minContainment);
        LOGGER.info("{}: representative peptide clustering algorithm={}, groups={}, representatives={}, files={}, comparedPairs={}",
                TASK_NAME, result.algorithmMode, result.groups.size(), result.representativeFiles.size(),
                result.assignments.size(), result.pairSimilarities.size());
        return result;
    }

    private static RepresentativePeptideClusterer.Input readMarkerSignature(int index, String raw, String rawName,
                                                                            File perFileDir,
                                                                            MetaParameterDiaMarker params) {
        if (perFileDir == null) {
            return null;
        }
        LinkedHashSet<String> peptides = new LinkedHashSet<>();
        File markerParquet = new File(perFileDir, "marker.parquet");
        if (markerParquet.exists() && markerParquet.length() > 0) {
            try {
                DiaNNParquetReader reader = new DiaNNParquetReader(markerParquet.getAbsolutePath(), "KR");
                DiaNNPrecursor[] precursors = reader.getDiaNNPrecursors(0.01);
                if (precursors != null) {
                    for (DiaNNPrecursor p : precursors) {
                        if (p == null || p.getIsDecoy()) continue;
                        String seq = p.getSeqString();
                        if (seq != null && !seq.isEmpty()) {
                            peptides.add(seq);
                        }
                    }
                }
            } catch (Exception e) {
                LOGGER.warn("{}: failed to read marker signature from {}", TASK_NAME, markerParquet.getAbsolutePath(), e);
            }
        }
        List<String> topGenomes = readTopGenomes(new File(perFileDir, "genome_ranked.tsv"),
                params.getRepresentativeGenomeTopN());
        return new RepresentativePeptideClusterer.Input(index, raw, rawName, peptides, topGenomes);
    }

    private static List<String> readTopGenomes(File rankedTsv, int topN) {
        if (rankedTsv == null || !rankedTsv.isFile() || topN <= 0) {
            return Collections.emptyList();
        }
        ArrayList<String> genomes = new ArrayList<>();
        try (BufferedReader reader = new BufferedReader(new FileReader(rankedTsv))) {
            String line;
            int genomeColumn = 0;
            boolean headerChecked = false;
            while ((line = reader.readLine()) != null && genomes.size() < topN) {
                if (line.isBlank()) continue;
                String[] cols = line.split("\t", -1);
                if (!headerChecked) {
                    headerChecked = true;
                    for (int i = 0; i < cols.length; i++) {
                        if ("Genome".equalsIgnoreCase(cols[i].trim())
                                || "GenomeId".equalsIgnoreCase(cols[i].trim())
                                || "Genome_ID".equalsIgnoreCase(cols[i].trim())) {
                            genomeColumn = i;
                            break;
                        }
                    }
                    if (cols.length > 0 && (cols[0].equalsIgnoreCase("rank")
                            || cols[0].equalsIgnoreCase("genome"))) {
                        continue;
                    }
                }
                if (genomeColumn >= cols.length) {
                    continue;
                }
                String genome = cols[genomeColumn].trim();
                if (genome.isEmpty() || genome.equalsIgnoreCase("genome") || genome.equalsIgnoreCase("genome_id")) {
                    continue;
                }
                genomes.add(genome);
            }
        } catch (IOException e) {
            LOGGER.warn("{}: failed to read top genomes from {}", TASK_NAME, rankedTsv.getAbsolutePath(), e);
        }
        return genomes;
    }

    private static void writeRepresentativeGroupsTsv(File out, RepresentativePeptideClusterer.Result grouping) {
        if (out == null || grouping == null) return;
        try (PrintWriter w = new PrintWriter(out)) {
            w.println("RawFile\tRawName\tGroupId\tRepresentativeRawFile\tRepresentativeRawName\tRole\t"
                    + "MarkerPeptides\tRepresentativeMarkerPeptides\tClusterSize\tPeptideJaccard\t"
                    + "PeptideContainment\tPeptideDistance\tGenomeTopNOverlapDiagnostic\t"
                    + "BestMatchRawFile\tBestMatchRawName\tBestPeptideJaccard\tBestPeptideContainment\t"
                    + "BestGenomeTopNOverlapDiagnostic\tReason");
            for (RepresentativePeptideClusterer.Assignment assignment : grouping.assignments) {
                RepresentativePeptideClusterer.Input file = assignment.file;
                RepresentativePeptideClusterer.Input rep = assignment.group.representative;
                RepresentativePeptideClusterer.Input best = assignment.bestMatch;
                String role = file.raw.equals(rep.raw) ? "representative" : "member";
                w.println(file.raw + "\t" + file.rawName + "\t" + assignment.group.groupId + "\t"
                        + rep.raw + "\t" + rep.rawName + "\t" + role + "\t"
                        + file.peptides.size() + "\t" + rep.peptides.size() + "\t"
                        + assignment.group.size() + "\t"
                        + formatDouble(assignment.peptideJaccardToRepresentative) + "\t"
                        + formatDouble(assignment.peptideContainmentToRepresentative) + "\t"
                        + formatDouble(1.0 - assignment.peptideJaccardToRepresentative) + "\t"
                        + formatDouble(assignment.genomeOverlapToRepresentative) + "\t"
                        + (best == null ? "" : best.raw) + "\t"
                        + (best == null ? "" : best.rawName) + "\t"
                        + formatDouble(assignment.bestPeptideJaccard) + "\t"
                        + formatDouble(assignment.bestPeptideContainment) + "\t"
                        + formatDouble(assignment.bestGenomeOverlap) + "\t" + assignment.reason);
            }
        } catch (IOException e) {
            LOGGER.warn("{}: failed to write representative groups TSV {}", TASK_NAME, out.getAbsolutePath(), e);
        }
    }

    private static void writeRepresentativePairwiseTsv(File out, RepresentativePeptideClusterer.Result grouping,
                                                       MetaParameterDiaMarker params) {
        if (out == null || grouping == null) return;
        try (PrintWriter w = new PrintWriter(out)) {
            w.println("RawFileA\tRawNameA\tRawFileB\tRawNameB\tMarkerPeptidesA\tMarkerPeptidesB\t"
                    + "IntersectionPeptides\tUnionPeptides\tSmallerPeptideSet\tLargerPeptideSet\t"
                    + "PeptideJaccard\tPeptideContainment\tPeptideDistance\tGenomeTopNOverlapDiagnostic\t"
                    + "PassesPeptideThresholds\tComparedFor");
            for (RepresentativePeptideClusterer.PairSimilarity sim : grouping.pairSimilarities) {
                boolean passes = sim.passes(params.getRepresentativePeptideJaccardThreshold(),
                        params.getRepresentativePeptideContainmentThreshold());
                w.println(sim.a.raw + "\t" + sim.a.rawName + "\t"
                        + sim.b.raw + "\t" + sim.b.rawName + "\t"
                        + sim.a.peptideCount() + "\t" + sim.b.peptideCount() + "\t"
                        + sim.intersectionPeptides + "\t" + sim.unionPeptides + "\t"
                        + sim.smallerPeptideSet + "\t" + sim.largerPeptideSet + "\t"
                        + formatDouble(sim.jaccard) + "\t"
                        + formatDouble(sim.containment) + "\t"
                        + formatDouble(1.0 - sim.jaccard) + "\t"
                        + formatDouble(sim.genomeTopNOverlap) + "\t"
                        + passes + "\t" + sim.source);
            }
        } catch (IOException e) {
            LOGGER.warn("{}: failed to write representative pairwise TSV {}", TASK_NAME, out.getAbsolutePath(), e);
        }
    }

    private static void writeRepresentativeClusteringSummaryTsv(File out,
            RepresentativePeptideClusterer.Result grouping, MetaParameterDiaMarker params) {
        if (out == null || grouping == null) return;
        try (PrintWriter w = new PrintWriter(out)) {
            w.println("Metric\tValue");
            w.println("Preset\t" + params.getRepresentativeClusteringPreset());
            w.println("PeptideJaccardThreshold\t" + formatDouble(params.getRepresentativePeptideJaccardThreshold()));
            w.println("PeptideContainmentThreshold\t"
                    + formatDouble(params.getRepresentativePeptideContainmentThreshold()));
            w.println("MinMarkerPeptides\t" + params.getRepresentativeMinMarkerPeptides());
            w.println("GenomeTopNDiagnostic\t" + params.getRepresentativeGenomeTopN());
            w.println("GenomeOverlapUsedForClustering\tfalse");
            w.println("GroupingDryRun\t" + params.isRepresentativeClusteringDryRun());
            w.println("AlgorithmMode\t" + grouping.algorithmMode);
            w.println("ExactPairwiseLimit\t" + RepresentativePeptideClusterer.EXACT_PAIRWISE_LIMIT);
            w.println("MinHashHashes\t" + RepresentativePeptideClusterer.MINHASH_HASHES);
            w.println("MinHashBands\t" + RepresentativePeptideClusterer.MINHASH_BANDS);
            w.println("MinHashRowsPerBand\t" + RepresentativePeptideClusterer.MINHASH_ROWS_PER_BAND);
            w.println("TotalRawFiles\t" + grouping.inputs.size());
            w.println("CandidatePairs\t" + grouping.candidatePairCount);
            w.println("ComparedPairs\t" + grouping.pairSimilarities.size());
            w.println("Clusters\t" + grouping.groups.size());
            w.println("Representatives\t" + grouping.representativeFiles.size());
            w.println("Singletons\t" + grouping.singletonCount);
        } catch (IOException e) {
            LOGGER.warn("{}: failed to write representative clustering summary TSV {}", TASK_NAME,
                    out.getAbsolutePath(), e);
        }
    }

    private static String formatDouble(double value) {
        return String.format(java.util.Locale.ROOT, "%.4f", value);
    }

    // ----------------------- DIA-NN per-round and per-file searches -----------------------

    /**
     * One iterative round for a single raw file:
     * <pre>
     *   roundN/combined.fasta = [proteins from N genomes] [+ identified peptides from previous round]
     *   roundN/combined.speclib = DIA-NN library predicted from the FASTA
     *   roundN/roundN.parquet  = DIA-NN no-quant identification against the new library
     * </pre>
     */
    private static int runRoundAppendPreviousPeptides(int round, File perFileDir, List<String> genomes,
                                                      String[] searchFiles, MagDbItem magDb, File diannExe,
                                                      MetaParameterDiaMarker params, File condaPrefix, int threads,
                                                      ProgressCallback cb) {
        return DeepIterativeCommonSupport.runAppendPreviousPeptidesRound(round, perFileDir, genomes, searchFiles, true,
                "append strategy",
                (gs, fastaOut) -> writeFastaFromGenomes(magDb, gs, fastaOut),
                (gs, prevParquet, fastaOut) -> writeFastaFromGenomesAndPreviousPeptides(magDb, gs, prevParquet, fastaOut),
                (combinedFasta, combinedSpeclib) -> {
                    try {
                        DiannParameter diannPar = DiaMarkerHeadlessRunner.buildDiannPar(params);
                        diannPar.setThreads(threads);
                        diannPar.setExportDecoy(true);
                        String cmd = diannPar.generateFasta(combinedFasta.getAbsolutePath(),
                                combinedSpeclib.getAbsolutePath(), false);
                        int exit = DiaNNHeadlessRunner.run(diannExe, cmd, combinedFasta.getParentFile(), cb, condaPrefix);
                        if (exit != 0) LOGGER.warn("{}: DIA-NN round speclib build exit code {}", TASK_NAME, exit);
                        File alt = new File(combinedFasta.getParentFile(), "combined.predicted.speclib");
                        if (!combinedSpeclib.exists() && alt.exists()) {
                            try { Files.move(alt.toPath(), combinedSpeclib.toPath()); } catch (Exception ignored) {}
                        }
                    } catch (Exception e) {
                        LOGGER.error("{}: round {} build speclib failed", TASK_NAME, round, e);
                    }
                },
                (roundParquet, combinedSpeclibPath, files, reanalyse) -> {
                    try {
                        DiaNNNoQuantLibrarySearchService.run(diannExe, params, files, new File(combinedSpeclibPath),
                                roundParquet, threads, reanalyse, false, cb, condaPrefix, LOGGER, TASK_NAME);
                    } catch (Exception e) {
                        LOGGER.error("{}: round {} no-quant search failed", TASK_NAME, round, e);
                    }
                },
                (parquetFile) -> DeepIterativeCommonSupport.countIdentifiedPeptides(parquetFile, "KR", LOGGER, TASK_NAME),
                LOGGER, TASK_NAME);
    }

    // ----------------------- global combined search -----------------------

    private static boolean runGlobalCombinedFromPerFileRounds(String[] diaFiles,
            LinkedHashMap<String, File> perFileDirs, LinkedHashMap<String, Integer> perFileLastRound,
            File combinedDir, File combinedFasta, File combinedSpeclib, File combinedParquet,
            File diannExe, MetaParameterDiaMarker params, File condaPrefix,
            boolean includeInitialMarkerEvidence, ProgressCallback cb) {

        // 1) Aggregate peptide -> proteins across ALL per-file rounds.
        HashMap<String, HashSet<String>> pepProMap = new HashMap<>();
        for (java.util.Map.Entry<String, File> entry : perFileDirs.entrySet()) {
            String raw = entry.getKey();
            File perFileDir = entry.getValue();
            Integer lastRound = perFileLastRound.get(raw);
            if (perFileDir == null || lastRound == null || lastRound <= 0) continue;
            HashMap<String, HashSet<String>> one = buildPepToProteinsFromRounds(perFileDir, lastRound);
            for (java.util.Map.Entry<String, HashSet<String>> e : one.entrySet()) {
                if (e.getKey() == null || e.getKey().isEmpty() || e.getValue() == null || e.getValue().isEmpty()) continue;
                HashSet<String> set = pepProMap.computeIfAbsent(e.getKey(), k -> new HashSet<>());
                set.addAll(e.getValue());
            }
        }

        if (includeInitialMarkerEvidence) {
            int before = pepProMap.size();
            for (Map.Entry<String, File> entry : perFileDirs.entrySet()) {
                File perFileDir = entry.getValue();
                if (perFileDir == null) continue;
                mergePeptideProteinMapFromParquet(new File(perFileDir, "marker.parquet"), pepProMap);
            }
            LOGGER.info("{}: added initial marker evidence to global combined FASTA input ({} -> {} peptides)",
                    TASK_NAME, before, pepProMap.size());
        }

        if (pepProMap.isEmpty()) {
            LOGGER.error("{}: no peptides collected from per-file round parquets for global combined FASTA", TASK_NAME);
            return false;
        }
        LOGGER.info("{}: aggregated {} unique peptides across all per-file rounds", TASK_NAME, pepProMap.size());
        reportProgress(cb, PROGRESS_COMBINED_START);

        // 2) Write combined.fasta
        if (!DeepIterativeCommonSupport.writeCombinedFastaFromPeptideProteinMap(combinedFasta, pepProMap, LOGGER, TASK_NAME)) {
            return false;
        }
        reportProgress(cb, PROGRESS_COMBINED_FASTA);

        // 3) Build combined.speclib
        if (!DeepIterativeCommonSupport.ensureCombinedSpeclib(combinedDir, combinedSpeclib, () -> {
            try {
                DiannParameter diannPar = DiaMarkerHeadlessRunner.buildDiannPar(params);
                diannPar.setThreads(params.getThreadCount());
                diannPar.setExportDecoy(true);
                String cmd = diannPar.generateFasta(combinedFasta.getAbsolutePath(), combinedSpeclib.getAbsolutePath(), false);
                int exit = DiaNNHeadlessRunner.run(diannExe, cmd, combinedDir, cb, condaPrefix);
                if (exit != 0) LOGGER.warn("{}: DIA-NN combined library build exit code {}", TASK_NAME, exit);
                File alt = new File(combinedDir, "combined.predicted.speclib");
                if (!combinedSpeclib.exists() && alt.exists()) {
                    try { Files.move(alt.toPath(), combinedSpeclib.toPath()); } catch (Exception ignored) {}
                }
            } catch (Exception e) {
                LOGGER.error("{}: combined speclib build failed", TASK_NAME, e);
            }
        }, LOGGER, TASK_NAME)) {
            return false;
        }
        reportProgress(cb, PROGRESS_COMBINED_SPECLIB);

        // 4) Final quantifying search across ALL DIA files
        boolean ok = DeepIterativeCommonSupport.ensureCombinedParquet(combinedParquet, () -> {
            try {
                DiannParameter diannPar = DiaMarkerHeadlessRunner.buildDiannPar(params);
                diannPar.setThreads(params.getThreadCount());
                diannPar.setMatrices(true);
                diannPar.setRelaxed_prot_inf(true);
                diannPar.setReanalyse(true);
                diannPar.setExportDecoy(true);

                DiaLibSearchPar libSearchPar = new DiaLibSearchPar();
                libSearchPar.setMbr(true);
                libSearchPar.setDiaFiles(diaFiles);
                libSearchPar.setThreadCount(params.getThreadCount());
                libSearchPar.setqValue(0.01);

                String cmd = diannPar.generateLibSearch(libSearchPar, combinedParquet.getAbsolutePath(),
                        combinedSpeclib.getAbsolutePath());
                cmd = cmd + DiaMarkerHeadlessRunner.buildFArgs(diaFiles);
                int exit = DiaNNHeadlessRunner.run(diannExe, cmd, combinedDir, cb, condaPrefix);
                if (exit != 0) LOGGER.warn("{}: DIA-NN combined quant exit code {}", TASK_NAME, exit);
            } catch (Exception e) {
                LOGGER.error("{}: combined quant failed", TASK_NAME, e);
            }
        });
        if (ok) {
            reportProgress(cb, PROGRESS_COMBINED_END);
        }
        return ok;
    }

    /** Build peptide -> proteins map by scanning roundN.parquet (q<=0.01, non-decoy) for one per-file folder. */
    private static HashMap<String, HashSet<String>> buildPepToProteinsFromRounds(File perFileDir, int lastRound) {
        if (perFileDir == null || lastRound <= 0) return new HashMap<>();
        return DiaDeepProteinMappingService.buildPeptideProteinMap(
                DiaDeepProteinMappingService.Request.deepTogether(perFileDir, lastRound, null, "KR"),
                LOGGER, TASK_NAME);
    }

    private static void mergePeptideProteinMapFromParquet(File parquet, HashMap<String, HashSet<String>> target) {
        if (parquet == null || !parquet.isFile() || parquet.length() <= 0 || target == null) {
            return;
        }
        try {
            DiaNNParquetReader reader = new DiaNNParquetReader(parquet.getAbsolutePath(), "KR");
            DiaNNPrecursor[] precursors = reader.getDiaNNPrecursors(0.01);
            if (precursors == null) return;
            for (DiaNNPrecursor p : precursors) {
                if (p == null || p.getIsDecoy()) continue;
                String seq = p.getSeqString();
                if (seq == null || seq.isEmpty()) continue;
                String[] proteins = p.getProGroups();
                if (proteins == null || proteins.length == 0) continue;
                HashSet<String> set = target.computeIfAbsent(seq, k -> new HashSet<>());
                for (String protein : proteins) {
                    if (protein == null || protein.isEmpty()) continue;
                    if (protein.startsWith("REV_") || protein.startsWith("rev_")
                            || protein.startsWith("DECOY_") || protein.startsWith("decoy_")) {
                        continue;
                    }
                    set.add(protein);
                }
            }
        } catch (Exception e) {
            LOGGER.warn("{}: failed to merge initial marker parquet {}", TASK_NAME, parquet.getAbsolutePath(), e);
        }
    }

    // ----------------------- per-round FASTA helpers -----------------------

    private static boolean writeFastaFromGenomes(MagDbItem magDb, List<String> genomes, File outputFasta) {
        File dbFolder = magDb.getOriginalDbFolder();
        if (dbFolder == null || !dbFolder.isDirectory()) {
            LOGGER.error("{}: MAG original DB folder not available", TASK_NAME);
            return false;
        }
        try (PrintWriter w = new PrintWriter(outputFasta)) {
            for (String genome : genomes) {
                File genomeFaa = resolveGenomeFile(dbFolder, genome);
                if (genomeFaa == null) {
                    LOGGER.warn("{}: genome file not found: {}", TASK_NAME, genome);
                    continue;
                }
                try (BufferedReader reader = new BufferedReader(new FileReader(genomeFaa))) {
                    String line;
                    while ((line = reader.readLine()) != null) w.println(line);
                }
            }
        } catch (IOException e) {
            LOGGER.error("{}: error writing FASTA from genomes", TASK_NAME, e);
            return false;
        }
        return outputFasta.exists() && outputFasta.length() > 0;
    }

    /**
     * FASTA = [genome proteins] + one entry per unique identified peptide from the previous round
     * (q <= 0.01, non-decoy). Each peptide entry uses a non-decoy protein id so DIA-NN can link
     * peptide -> protein.
     */
    private static boolean writeFastaFromGenomesAndPreviousPeptides(MagDbItem magDb, List<String> genomes,
            File previousRoundParquet, File outputFasta) {
        File dbFolder = magDb.getOriginalDbFolder();
        if (dbFolder == null || !dbFolder.isDirectory()) {
            LOGGER.error("{}: MAG original DB folder not available", TASK_NAME);
            return false;
        }
        DiaNNPrecursor[] precursors = null;
        if (previousRoundParquet.exists()) {
            try {
                DiaNNParquetReader reader = new DiaNNParquetReader(previousRoundParquet.getAbsolutePath(), "KR");
                precursors = reader.getDiaNNPrecursors(0.01);
            } catch (Exception e) {
                LOGGER.warn("{}: failed to read previous round parquet {}",
                        TASK_NAME, previousRoundParquet.getAbsolutePath(), e);
            }
        }
        Set<String> addedPeptideSequences = new HashSet<>();
        try (PrintWriter w = new PrintWriter(outputFasta)) {
            for (String genome : genomes) {
                File genomeFaa = resolveGenomeFile(dbFolder, genome);
                if (genomeFaa == null) {
                    LOGGER.warn("{}: genome file not found: {}", TASK_NAME, genome);
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
                            // Skip common decoy prefixes
                            if (pro.startsWith("REV_") || pro.startsWith("rev_") || pro.startsWith("DECOY_")
                                    || pro.startsWith("decoy_")) continue;
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
        } catch (IOException e) {
            LOGGER.error("{}: error writing FASTA (genomes + previous peptides)", TASK_NAME, e);
            return false;
        }
        return outputFasta.exists() && outputFasta.length() > 0;
    }

    private static File resolveGenomeFile(File dbFolder, String genome) {
        File f = new File(dbFolder, genome + ".faa");
        if (f.exists()) return f;
        f = new File(dbFolder, genome + ".fasta");
        if (f.exists()) return f;
        return null;
    }

    // ----------------------- stop rule (matches DEEP_MARKER profile) -----------------------

    private static DeepIterativeStopRule.StopRuleResult computeStopRule(List<Integer> increments) {
        return DeepIterativeStopRule.computeStopRoundFromIncrements(increments,
                RELATIVE_GAIN_THRESHOLD, CONSECUTIVE_ROUNDS_FOR_PLATEAU,
                DeepIterativeStopRule.AlgorithmProfile.DEEP_MARKER);
    }

    // ----------------------- misc -----------------------

}
