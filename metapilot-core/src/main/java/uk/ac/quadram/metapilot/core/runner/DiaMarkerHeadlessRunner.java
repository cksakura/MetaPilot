package uk.ac.quadram.metapilot.core.runner;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import uk.ac.quadram.metapilot.core.dia.Marker4LibrarySupport;
import uk.ac.quadram.metapilot.core.diann.DiaNNHeadlessRunner;
import uk.ac.quadram.metapilot.core.mag.MagGenomeSelectionService;
import uk.ac.quadram.metapilot.core.mag.Marker4ResourcePreparer;
import uk.ac.quadram.metapilot.core.report.MetaPilotHtmlReportService;
import uk.ac.quadram.metapilot.core.resources.MetaPilotResourceConfig;
import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;
import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNParquetReader;
import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNPrecursor;
import uk.ac.quadram.metapilot.dbSearch.diann.DiannParameter;
import uk.ac.quadram.metapilot.task.dia.DiaLibSearchPar;
import uk.ac.quadram.metapilot.task.dia.DiaRefinementMode;
import uk.ac.quadram.metapilot.task.dia.par.MetaParameterDiaMarker;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;
import uk.ac.quadram.metapilot.core.MGYG.HapExtractor;

import java.io.*;
import java.nio.file.Files;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Headless DIA Marker workflow runner (Marker4 -> genome selection -> combined library -> quant matrices).
 *
 * Outputs (under resultDir):
 * - marker/marker.parquet
 * - marker/genome_evidence.tsv, marker/genome_ranked.tsv, marker/genome_selected.tsv
 * - marker/combined.fasta
 * - marker/combined.speclib
 * - combined/combined.parquet + combined.pr_matrix.tsv + combined.pg_matrix.tsv
 *
 * This runner uses the shared final HTML report service after post-analysis.
 */
public final class DiaMarkerHeadlessRunner {

    private static final Logger LOGGER = LogManager.getLogger(DiaMarkerHeadlessRunner.class);

    private DiaMarkerHeadlessRunner() {
    }

    public static final class Result {
        public final File resultDir;
        public final File markerDir;
        public final File combinedDir;
        public final File markerParquet;
        public final File rankedTsv;
        public final File selectedTsv;
        public final File combinedFasta;
        public final File combinedSpeclib;
        public final File combinedParquet;
        public final File prMatrixTsv;
        public final File pgMatrixTsv;

        Result(File resultDir, File markerDir, File combinedDir, File markerParquet, File rankedTsv, File selectedTsv,
               File combinedFasta, File combinedSpeclib, File combinedParquet, File prMatrixTsv, File pgMatrixTsv) {
            this.resultDir = resultDir;
            this.markerDir = markerDir;
            this.combinedDir = combinedDir;
            this.markerParquet = markerParquet;
            this.rankedTsv = rankedTsv;
            this.selectedTsv = selectedTsv;
            this.combinedFasta = combinedFasta;
            this.combinedSpeclib = combinedSpeclib;
            this.combinedParquet = combinedParquet;
            this.prMatrixTsv = prMatrixTsv;
            this.pgMatrixTsv = pgMatrixTsv;
        }
    }

    public static Result run(
            MetaParameterDiaMarker params,
            File diannExe,
            File magDbConfigJson,
            ProgressCallback cb
    ) throws Exception {
        return run(params, diannExe, magDbConfigJson, cb, null);
    }

    /**
     * @param condaPrefix optional conda env prefix used to inject {@code <prefix>/lib} into LD_LIBRARY_PATH for DIA-NN (Linux)
     */
    public static Result run(
            MetaParameterDiaMarker params,
            File diannExe,
            File magDbConfigJson,
            ProgressCallback cb,
            File condaPrefix
    ) throws Exception {
        return run(params, diannExe, magDbConfigJson, cb, condaPrefix, null);
    }

    public static Result run(
            MetaParameterDiaMarker params,
            File diannExe,
            File magDbConfigJson,
            ProgressCallback cb,
            File condaPrefix,
            File funcDbFile
    ) throws Exception {
        return run(params, diannExe, magDbConfigJson, cb, condaPrefix, funcDbFile, true, null);
    }

    public static Result run(
            MetaParameterDiaMarker params,
            File diannExe,
            File magDbConfigJson,
            ProgressCallback cb,
            File condaPrefix,
            File funcDbFile,
            boolean generateHtmlReport,
            String reportPythonHint
    ) throws Exception {
        return run(params, diannExe, magDbConfigJson, cb, condaPrefix, funcDbFile, generateHtmlReport,
                reportPythonHint, null);
    }

    public static Result run(
            MetaParameterDiaMarker params,
            File diannExe,
            File magDbConfigJson,
            ProgressCallback cb,
            File condaPrefix,
            File funcDbFile,
            boolean generateHtmlReport,
            String reportPythonHint,
            MetaPilotResourceConfig reportResourceConfig
    ) throws Exception {
        return runInternal(params, diannExe, magDbConfigJson, cb, condaPrefix, funcDbFile, true,
                generateHtmlReport, reportPythonHint, reportResourceConfig);
    }

    public static Result runWithoutPostAnalysis(
            MetaParameterDiaMarker params,
            File diannExe,
            File magDbConfigJson,
            ProgressCallback cb,
            File condaPrefix,
            File funcDbFile
    ) throws Exception {
        return runInternal(params, diannExe, magDbConfigJson, cb, condaPrefix, funcDbFile, false,
                false, null, null);
    }

    private static Result runInternal(
            MetaParameterDiaMarker params,
            File diannExe,
            File magDbConfigJson,
            ProgressCallback cb,
            File condaPrefix,
            File funcDbFile,
            boolean exportPostAnalysis,
            boolean generateHtmlReport,
            String reportPythonHint,
            MetaPilotResourceConfig reportResourceConfig
    ) throws Exception {
        if (cb == null) cb = new ProgressCallback.NoOp();
        if (params == null) throw new IllegalArgumentException("params is null");
        if (diannExe == null || !diannExe.exists()) throw new IOException("DIA-NN not found: " + (diannExe == null ? "null" : diannExe.getAbsolutePath()));

        DiaRefinementMode markerMode = params.getMarkerRefinementMode();
        if (markerMode == DiaRefinementMode.DIRECT) {
            throw new IOException("DIA Marker does not support DIRECT mode. DIRECT is only valid for DIA Seed workflows.");
        }

        // Dispatch deep marker modes to their dedicated headless runners. Returns are adapted
        // to this runner's Result type so existing callers keep working.
        if (markerMode == DiaRefinementMode.DEEP_TOGETHER) {
            cb.onMessage("DIA Marker: dispatching to DEEP_TOGETHER headless runner");
            DiaMarkerDeepTogetherHeadlessRunner.Result deepTogetherResult = exportPostAnalysis
                    ? DiaMarkerDeepTogetherHeadlessRunner.run(params, diannExe, magDbConfigJson, cb, condaPrefix,
                            funcDbFile, generateHtmlReport, reportPythonHint, reportResourceConfig)
                    : DiaMarkerDeepTogetherHeadlessRunner.runWithoutPostAnalysis(params, diannExe, magDbConfigJson, cb,
                            condaPrefix, funcDbFile);
            return new Result(
                    deepTogetherResult.markerDeepDir,
                    deepTogetherResult.markerDir,
                    deepTogetherResult.combinedDir,
                    deepTogetherResult.markerParquet,
                    deepTogetherResult.rankedTsv,
                    deepTogetherResult.selectedTsv,
                    deepTogetherResult.combinedFasta,
                    deepTogetherResult.combinedSpeclib,
                    deepTogetherResult.combinedParquet,
                    deepTogetherResult.prMatrixTsv,
                    deepTogetherResult.pgMatrixTsv
            );
        }
        if (markerMode == DiaRefinementMode.DEEP_PER_FILE) {
            cb.onMessage("DIA Marker: dispatching to DEEP_PER_FILE headless runner");
            DiaMarkerDeepPerFileHeadlessRunner.Result perFileResult = exportPostAnalysis
                    ? DiaMarkerDeepPerFileHeadlessRunner.run(params, diannExe, magDbConfigJson, cb, condaPrefix,
                            funcDbFile, generateHtmlReport, reportPythonHint, reportResourceConfig)
                    : DiaMarkerDeepPerFileHeadlessRunner.runWithoutPostAnalysis(params, diannExe, magDbConfigJson, cb,
                            condaPrefix, funcDbFile);
            return new Result(
                    perFileResult.markerDeepDir,
                    new File(perFileResult.markerDeepDir, "marker"),
                    perFileResult.combinedDir,
                    null,
                    null,
                    null,
                    perFileResult.combinedFasta,
                    perFileResult.combinedSpeclib,
                    perFileResult.combinedParquet,
                    perFileResult.prMatrixTsv,
                    perFileResult.pgMatrixTsv
            );
        }

        File studyRoot = new File(params.getResult());
        if (!studyRoot.exists()) studyRoot.mkdirs();
        if (!studyRoot.isDirectory()) throw new IOException("Result dir is not a directory: " + studyRoot.getAbsolutePath());

        // Mirror the Windows task: export under marker_focus / marker_comprehensive,
        // matching the folder name chosen by MetaDiaTask for non-DEEP_PER_FILE modes.
        String workFolderName;
        if (markerMode.usesComprehensiveSelection()) {
            workFolderName = "marker_comprehensive";
        } else {
            workFolderName = "marker_focus";
        }
        File resultDir = new File(studyRoot, workFolderName);
        if (!resultDir.exists()) resultDir.mkdirs();
        cb.onMessage("DIA Marker: result folder -> " + resultDir.getAbsolutePath());

        String[] diaFiles = params.getMetadata() != null ? params.getMetadata().getRawFiles() : null;
        if (diaFiles == null || diaFiles.length == 0) {
            throw new IOException("No DIA files found in params metadata (input_output.experiments)");
        }

        MagDbItem magDb = resolveMagDb(params, magDbConfigJson);
        if (magDb == null) throw new IOException("MAG DB could not be resolved from params + magdb config");

        // Ensure Marker4 FASTA and genomes-taxa.tsv exist.
        ensureMarker4Files(magDb, diannExe, params, cb);

        cb.onMessage("DIA Marker: preparing Marker4 library");
        ensureMarker4Library(magDb, diannExe, params, cb, condaPrefix);

        // Marker no-quant search -> genome selection -> combined.fasta + combined.speclib
        File markerDir = new File(resultDir, "marker");
        if (!markerDir.exists()) markerDir.mkdirs();

        File markerParquet = new File(markerDir, "marker.parquet");
        if (!markerParquet.exists()) {
            cb.onMessage("DIA Marker: searching Marker4 library (no quant) to select genomes");
            runMarkerNoQuantSearch(diannExe, params, diaFiles, magDb.resolveMarker4LibraryFile(), markerParquet, cb, condaPrefix);
        }
        if (!markerParquet.exists()) throw new IOException("Marker search did not produce: " + markerParquet.getAbsolutePath());

        File currentFile = magDb.getCurrentFile();
        File markerFasta = new File(currentFile, "Marker4.fasta");
        File genomesTaxaTsv = new File(currentFile, "genomes-taxa.tsv");
        if (!markerFasta.exists() || !genomesTaxaTsv.exists()) {
            throw new IOException("Marker4.fasta or genomes-taxa.tsv missing under: " + currentFile.getAbsolutePath());
        }

        File evidenceTsv = new File(markerDir, "genome_evidence.tsv");
        File rankedTsv = new File(markerDir, "genome_ranked.tsv");
        File selectedTsv = new File(markerDir, "genome_selected.tsv");
        HashSet<String> selectedGenomes = new HashSet<>();

        if (!selectedTsv.exists()) {
            cb.onMessage("DIA Marker: reading precursors for genome ranking");
            DiaNNPrecursor[] precursors = readPrecursors(markerParquet, params);
            if (precursors == null || precursors.length == 0) {
                throw new IOException("No precursors found in marker parquet: " + markerParquet.getAbsolutePath());
            }

            String kneeMode = markerModeToKneeMode(params.getMarkerMode());
            cb.onMessage("DIA Marker: selecting genomes (mode=" + kneeMode + ")");
            MagGenomeSelectionService.Request selectionRequest = new MagGenomeSelectionService.Request();
            selectionRequest.magDb = magDb;
            selectionRequest.genomeSelectorMode = MagDbItem.MARKER4;
            selectionRequest.genomeSelectionStrategy = kneeMode;
            selectionRequest.diaPrecursors = precursors;
            selectionRequest.markerFasta = markerFasta;
            selectionRequest.genomesTaxaTsv = genomesTaxaTsv;
            selectionRequest.evidenceFile = evidenceTsv;
            selectionRequest.rankedFile = rankedTsv;
            selectionRequest.selectedFile = selectedTsv;
            selectionRequest.maximumSelected = 300;
            MagGenomeSelectionService.selectGenomes(selectionRequest, cb);
        }

        readSelectedGenomesIntoSet(selectedTsv, selectedGenomes);
        if (selectedGenomes.isEmpty()) {
            throw new IOException("No genomes selected (empty genome_selected.tsv): " + selectedTsv.getAbsolutePath());
        }

        File combinedFasta = new File(markerDir, "combined.fasta");
        if (!combinedFasta.exists() || combinedFasta.length() == 0) {
            cb.onMessage("DIA Marker: building combined.fasta from " + selectedGenomes.size() + " genomes");
            buildCombinedFasta(magDb, selectedGenomes, combinedFasta);
        }
        if (!combinedFasta.exists() || combinedFasta.length() == 0) {
            throw new IOException("combined.fasta was not created: " + combinedFasta.getAbsolutePath());
        }

        File combinedSpeclib = new File(markerDir, "combined.speclib");
        if (!combinedSpeclib.exists()) {
            cb.onMessage("DIA Marker: generating combined.speclib from combined.fasta");
            DiannParameter diannPar = buildDiannPar(params);
            String cmd = diannPar.generateFasta(combinedFasta.getAbsolutePath(), combinedSpeclib.getAbsolutePath(), false);
            int exit = DiaNNHeadlessRunner.run(diannExe, cmd, markerDir, cb, condaPrefix);
            if (exit != 0) LOGGER.warn("DIA-NN library build exit code {}", exit);

            // DIA-NN may create *.predicted.speclib in some setups; try to adopt it.
            File alt = new File(markerDir, "combined.predicted.speclib");
            if (!combinedSpeclib.exists() && alt.exists()) {
                Files.move(alt.toPath(), combinedSpeclib.toPath());
            }
        }
        if (!combinedSpeclib.exists()) throw new IOException("combined.speclib not created: " + combinedSpeclib.getAbsolutePath());

        // Combined quantifying search
        File combinedDir = new File(resultDir, "combined");
        if (!combinedDir.exists()) combinedDir.mkdirs();

        File combinedParquet = new File(combinedDir, "combined.parquet");
        File prMatrix = new File(combinedDir, "combined.pr_matrix.tsv");
        File pgMatrix = new File(combinedDir, "combined.pg_matrix.tsv");

        if (!combinedParquet.exists() || !prMatrix.exists() || !pgMatrix.exists()) {
            cb.onMessage("DIA Marker: running DIA-NN quant with combined library");
            runCombinedQuantSearch(diannExe, params, diaFiles, combinedSpeclib, combinedParquet, cb, condaPrefix);
        }

        if (!combinedParquet.exists()) throw new IOException("Combined quant did not produce: " + combinedParquet.getAbsolutePath());
        if (!prMatrix.exists()) throw new IOException("Combined quant did not produce: " + prMatrix.getAbsolutePath());
        if (!pgMatrix.exists()) throw new IOException("Combined quant did not produce: " + pgMatrix.getAbsolutePath());

        runHostSearchIfRequested(resultDir, diannExe, params, diaFiles, cb, condaPrefix);

        // Post-analysis + TSV/report exports under resultDir (marker_focus / marker_comprehensive),
        // matching the Windows task layout.
        if (exportPostAnalysis) {
            try {
                File gtdbNcbi = null;
                File maybe = new File("Resources", "gtdb_ncbi");
                if (maybe.exists()) gtdbNcbi = maybe;
                DiaMagPostAnalysisTsvExporter.export(
                        resultDir,
                        magDb,
                        combinedParquet,
                        prMatrix,
                        gtdbNcbi,
                        params.getMarkerMode(),
                        funcDbFile,
                        params.getMicroDb(),
                        params.getQuantitativeMode(),
                        cb
                );
                runHtmlReportIfRequested(resultDir, params, generateHtmlReport, reportPythonHint,
                        reportResourceConfig, cb);
            } catch (Exception e) {
                cb.onMessage("WARNING: DIA marker post-analysis export failed: " + e.getMessage());
                throw e;
            }
        }

        cb.onComplete("DIA marker complete: " + resultDir.getAbsolutePath());
        return new Result(resultDir, markerDir, combinedDir, markerParquet, rankedTsv, selectedTsv,
                combinedFasta, combinedSpeclib, combinedParquet, prMatrix, pgMatrix);
    }

    static MagDbItem resolveMagDb(MetaParameterDiaMarker params, File magDbConfigJson) throws Exception {
        MagDbItem fromParams = params.getMagDbItem();
        if (fromParams != null && fromParams.getCurrentFile() != null && fromParams.getCurrentFile().exists()) {
            return fromParams;
        }
        MagDbItem resolved = uk.ac.quadram.metapilot.task.mag.MagDbItemResolver.resolve(
                params.getMicroDb(),
                params.getMagDbIdentifier(),
                params.getMagDbCatalogueVersion(),
                magDbConfigJson);
        if (resolved != null) {
            params.setMagDbItem(resolved);
            return resolved;
        }
        return fromParams;
    }

    static void ensureMarker4Library(MagDbItem magDb, File diannExe, MetaParameterDiaMarker params, ProgressCallback cb, File condaPrefix) throws Exception {
        // Migrate legacy layout if present.
        magDb.moveLegacyMarker4SpeclibToCanonical();
        if (!magDb.ensureMarker4LibraryFolderExists()) {
            throw new IOException("Could not create Marker4 library folder under MAG DB: " + magDb.getCurrentFile());
        }

        File fastaFile = magDb.getMarker4FastaForLibrary();
        File libFile = magDb.getMarker4LibFile();
        if (fastaFile == null || !fastaFile.exists()) {
            throw new IOException("Marker4 FASTA not found: " + (fastaFile == null ? "null" : fastaFile.getAbsolutePath()));
        }
        if (libFile == null) throw new IOException("Marker4 lib path is null");

        if (!libFile.exists()) {
            cb.onMessage("Marker4: generating libraries/Marker4/Marker4.speclib from " + fastaFile.getName());
            DiannParameter diannPar = buildDiannPar(params);
            String cmd = diannPar.generateFasta(fastaFile.getAbsolutePath(), libFile.getAbsolutePath(), false);
            int exit = DiaNNHeadlessRunner.run(diannExe, cmd, libFile.getParentFile(), cb, condaPrefix);
            if (exit != 0) LOGGER.warn("DIA-NN Marker4 library build exit code {}", exit);

            // Some DIA-NN builds write <name>.predicted.speclib.
            File predicted = new File(libFile.getParentFile(), "Marker4.predicted.speclib");
            if (!libFile.exists() && predicted.exists()) {
                Files.move(predicted.toPath(), libFile.toPath());
            }
        }

        if (!libFile.exists()) {
            throw new IOException("Failed to create Marker4 library speclib: " + libFile.getAbsolutePath());
        }

        Marker4LibrarySupport.ensureMarker4LibraryMetadata(
                diannExe,
                condaPrefix,
                magDb,
                params,
                cb,
                LOGGER
        );
    }

    static void ensureMarker4Files(MagDbItem magDb, File diannExe, MetaParameterDiaMarker params, ProgressCallback cb) throws Exception {
        Marker4ResourcePreparer.ensureSelectionResources(magDb, cb);
    }

    private static void runMarkerNoQuantSearch(File diannExe, MetaParameterDiaMarker params, String[] diaFiles,
                                              File markerLibFile, File outParquet, ProgressCallback cb, File condaPrefix) throws Exception {
        DiaNNNoQuantLibrarySearchService.run(diannExe, params, diaFiles, markerLibFile, outParquet,
                params.getThreadCount(), true, true, cb, condaPrefix, LOGGER, "DIA Marker");
    }

    private static void runCombinedQuantSearch(File diannExe, MetaParameterDiaMarker params, String[] diaFiles,
                                              File combinedLibFile, File outParquet, ProgressCallback cb, File condaPrefix) throws Exception {
        if (combinedLibFile == null || !combinedLibFile.exists()) {
            throw new IOException("Combined library not found: " + (combinedLibFile == null ? "null" : combinedLibFile.getAbsolutePath()));
        }
        DiannParameter diannPar = buildDiannPar(params);
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

        String cmd = diannPar.generateLibSearch(libSearchPar, outParquet.getAbsolutePath(), combinedLibFile.getAbsolutePath());
        cmd = cmd + buildFArgs(diaFiles);
        int exit = DiaNNHeadlessRunner.run(diannExe, cmd, outParquet.getParentFile(), cb, condaPrefix);
        if (exit != 0) LOGGER.warn("DIA-NN combined quant exit code {}", exit);
    }

    private static void runHtmlReportIfRequested(File resultDir, MetaParameterDiaMarker params,
                                                 boolean generateHtmlReport, String reportPythonHint,
                                                 MetaPilotResourceConfig reportResourceConfig,
                                                 ProgressCallback cb) {
        if (!generateHtmlReport) {
            cb.onMessage("HTML report generation skipped (disabled)");
            return;
        }
        cb.onMessage("DIA Marker: generating shared HTML report");
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

    public static void runHostSearchIfRequested(File resultDir, File diannExe, MetaParameterDiaMarker params,
                                                String[] diaFiles, ProgressCallback cb, File condaPrefix) throws Exception {
        if (cb == null) cb = new ProgressCallback.NoOp();
        if (params == null || !params.isAppendHostDb()) {
            return;
        }
        String hostDb = params.getHostDB();
        if (hostDb == null || hostDb.isBlank()) {
            throw new IOException("Append host DB is enabled but host database path is empty");
        }
        File hostInput = resolveUserPath(hostDb);
        if (!hostInput.exists()) {
            throw new IOException("Host database file not found: " + hostInput.getAbsolutePath());
        }
        if (resultDir == null) {
            throw new IOException("Result directory is null for host search");
        }
        String[] files = diaFiles != null ? diaFiles
                : (params.getMetadata() != null ? params.getMetadata().getRawFiles() : null);
        if (files == null || files.length == 0) {
            throw new IOException("No DIA files available for host search");
        }

        File hostDir = new File(resultDir, "host");
        if (!hostDir.exists() && !hostDir.mkdirs()) {
            throw new IOException("Failed to create host search directory: " + hostDir.getAbsolutePath());
        }
        cb.onMessage("DIA Marker: running optional host DB search -> " + hostDir.getAbsolutePath());

        File hostLibrary = resolveHostLibraryForSearch(hostInput, hostDir, diannExe, params, cb, condaPrefix);

        File hostParquet = new File(hostDir, "host.parquet");
        File hostPrMatrix = new File(hostDir, "host.pr_matrix.tsv");
        File hostPgMatrix = new File(hostDir, "host.pg_matrix.tsv");
        if (!hostParquet.exists() || !hostPrMatrix.exists() || !hostPgMatrix.exists()) {
            cb.onMessage("DIA Marker: running DIA-NN quant with host library");
            runCombinedQuantSearch(diannExe, params, files, hostLibrary, hostParquet, cb, condaPrefix);
        }
        if (!hostParquet.exists()) {
            throw new IOException("Host search did not produce: " + hostParquet.getAbsolutePath());
        }
        if (!hostPrMatrix.exists()) {
            throw new IOException("Host search did not produce: " + hostPrMatrix.getAbsolutePath());
        }
        if (!hostPgMatrix.exists()) {
            throw new IOException("Host search did not produce: " + hostPgMatrix.getAbsolutePath());
        }
    }

    private static File resolveHostLibraryForSearch(File hostInput, File hostDir, File diannExe,
                                                    MetaParameterDiaMarker params, ProgressCallback cb,
                                                    File condaPrefix) throws Exception {
        String name = hostInput.getName().toLowerCase();
        if (name.endsWith(".speclib") || name.endsWith(".parquet")) {
            cb.onMessage("DIA Marker: using host library directly: " + hostInput.getAbsolutePath());
            return hostInput;
        }
        if (!name.endsWith(".fasta") && !name.endsWith(".fa")) {
            throw new IOException("Unsupported DIA host database file type: " + hostInput.getAbsolutePath()
                    + " (expected .fasta, .fa, .speclib, or .parquet)");
        }

        File hostSpeclib = new File(hostDir, "host.speclib");
        if (!hostSpeclib.exists()) {
            cb.onMessage("DIA Marker: generating host library from FASTA: " + hostInput.getAbsolutePath());
            DiannParameter diannPar = buildDiannPar(params);
            String cmd = diannPar.generateFasta(hostInput.getAbsolutePath(), hostSpeclib.getAbsolutePath(), false);
            int exit = DiaNNHeadlessRunner.run(diannExe, cmd, hostDir, cb, condaPrefix);
            if (exit != 0) LOGGER.warn("DIA-NN host library build exit code {}", exit);

            File predicted = new File(hostDir, "host.predicted.speclib");
            if (!hostSpeclib.exists() && predicted.exists()) {
                Files.move(predicted.toPath(), hostSpeclib.toPath());
            }
        }
        if (!hostSpeclib.exists()) {
            throw new IOException("Host speclib was not created: " + hostSpeclib.getAbsolutePath());
        }
        return hostSpeclib;
    }

    private static File resolveUserPath(String path) {
        if (path == null) {
            return null;
        }
        String trimmed = path.trim();
        if (trimmed.equals("~")) {
            String home = System.getProperty("user.home", "");
            return home == null || home.isBlank() ? new File(trimmed) : new File(home);
        }
        if (trimmed.startsWith("~/") || trimmed.startsWith("~\\")) {
            String home = System.getProperty("user.home", "");
            if (home == null || home.isBlank()) {
                return new File(trimmed);
            }
            return new File(home, trimmed.substring(2));
        }
        String home = System.getenv("HOME");
        if (home != null && !home.isBlank()) {
            if (trimmed.equals("$HOME")) {
                return new File(home);
            }
            if (trimmed.startsWith("$HOME/") || trimmed.startsWith("$HOME\\")) {
                return new File(home, trimmed.substring(6));
            }
        }
        return new File(trimmed);
    }

    static String buildFArgs(String[] diaFiles) {
        StringBuilder sb = new StringBuilder();
        for (String f : diaFiles) {
            if (f == null || f.isBlank()) continue;
            sb.append(" --f \"").append(f.trim()).append("\"");
        }
        return sb.toString();
    }

    static DiaNNPrecursor[] readPrecursors(File parquet, MetaParameterDiaMarker params) throws Exception {
        DiannParameter diannPar = buildDiannPar(params);
        DiaNNParquetReader parquetReader = new DiaNNParquetReader(parquet.getAbsolutePath(), diannPar.getCut());
        return parquetReader.getDiaNNPrecursors(0.01);
    }

    static DiannParameter buildDiannPar(MetaParameterDiaMarker params) {
        DiannParameter diannPar = new DiannParameter();
        diannPar.setThreads(params.getThreadCount());
        diannPar.setEnzyme("Trypsin/P");
        diannPar.setMissed_cleavages(1);

        diannPar.setM_excision(params.isMExcision());
        diannPar.setCarbami_C(params.isCarbamiC());
        diannPar.setOxidation_M(params.isOxidationM());
        diannPar.setAcety_proN(params.isAcetyProN());
        diannPar.setPhospho_STY(params.isPhosphoSTY());
        diannPar.setK_GG(params.isKGG());

        return diannPar;
    }

    static String markerModeToKneeMode(String markerMode) {
        DiaRefinementMode mode = DiaRefinementMode.fromMarkerMode(markerMode);
        if (mode.usesComprehensiveSelection()) {
            return MagGenomeSelectionService.SELECTION_STRATEGY_COMPREHENSIVE;
        }
        return MagGenomeSelectionService.SELECTION_STRATEGY_FOCUS;
    }

    static void readSelectedGenomesIntoSet(File selectedFile, HashSet<String> genomeSet) {
        if (selectedFile == null || !selectedFile.exists() || genomeSet == null) return;
        try (BufferedReader reader = new BufferedReader(new FileReader(selectedFile))) {
            reader.readLine(); // header
            String line;
            while ((line = reader.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty() || line.startsWith("#")) break;
                String[] parts = line.split("\t");
                if (parts.length > 0 && !parts[0].trim().isEmpty()) genomeSet.add(parts[0].trim());
            }
        } catch (IOException e) {
            LOGGER.warn("Error reading selected genomes from {}", selectedFile.getAbsolutePath(), e);
        }
    }

    static void buildCombinedFasta(MagDbItem magDb, HashSet<String> genomes, File out) throws IOException {
        File dbFolder = magDb.getOriginalDbFolder();
        if (dbFolder == null || !dbFolder.exists()) {
            throw new IOException("Original DB folder not found for combined FASTA build");
        }
        try (PrintWriter writer = new PrintWriter(out)) {
            for (String genome : genomes) {
                File genomeDbFile = new File(dbFolder, genome + ".faa");
                if (!genomeDbFile.exists()) genomeDbFile = new File(dbFolder, genome + ".fasta");
                if (!genomeDbFile.exists()) {
                    LOGGER.warn("Genome file not found: {}", genome);
                    continue;
                }
                try (BufferedReader reader = new BufferedReader(new FileReader(genomeDbFile))) {
                    String line;
                    while ((line = reader.readLine()) != null) {
                        writer.println(line);
                    }
                }
            }
        }
        if (!out.exists() || out.length() == 0) {
            throw new IOException("Combined FASTA empty after build: " + out.getAbsolutePath());
        }
    }
}

