package uk.ac.quadram.metapilot.core.runner;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import uk.ac.quadram.metapilot.core.dia.DiaDeepProteinMappingService;
import uk.ac.quadram.metapilot.core.dia.DiaHostLibrarySupport;
import uk.ac.quadram.metapilot.core.dia.MagPredictedPeptideProteinMapper;
import uk.ac.quadram.metapilot.core.diann.DiaNNHeadlessRunner;
import uk.ac.quadram.metapilot.core.mag.MagGenomeSelectionService;
import uk.ac.quadram.metapilot.core.mag.Marker4ResourcePreparer;
import uk.ac.quadram.metapilot.core.report.MetaPilotHtmlReportService;
import uk.ac.quadram.metapilot.core.resources.MetaPilotResourceConfig;
import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;
import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNParquetReader;
import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNPrecursor;
import uk.ac.quadram.metapilot.dbSearch.diann.DiannParameter;
import uk.ac.quadram.metapilot.task.dia.DeepIterativeCommonSupport;
import uk.ac.quadram.metapilot.task.dia.DeepIterativeStopRule;
import uk.ac.quadram.metapilot.task.dia.DiaLibSearchPar;
import uk.ac.quadram.metapilot.task.dia.DiaRefinementMode;
import uk.ac.quadram.metapilot.task.dia.par.MetaParameterDia;
import uk.ac.quadram.metapilot.task.dia.par.MetaParameterDiaMarker;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;
import uk.ac.quadram.metapilot.task.mag.MagDbItem.MagPepLib;
import uk.ac.quadram.metapilot.task.mag.MarkerGenomeSelector;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;

/**
 * Headless DIA Seed workflow, following the Windows tasks step by step.
 * <ul>
 *   <li>Core / Focus and Comprehensive ({@code DiaSeedCoreTask}): seed library search -> peptide-to-MAG-protein
 *   mapping via {@code predicted/} -> Marker4 genome selection -> combined library -> quant search ->
 *   post-analysis.</li>
 *   <li>Direct ({@code DiaSeedDirectTask}): one quantifying search against the seed library, then
 *   post-analysis; the exporter maps the seed peptides to MAG proteins via {@code predicted/}. Output under
 *   {@code seed_direct}.</li>
 *   <li>Deep together ({@code DiaSeedDeepTogetherTask}): the same seed library search and mapping, comprehensive
 *   genome ranking, iterative rounds of {@code genomes_per_round} genomes across all DIA files (stop rule
 *   DEEP_COMPREHENSIVE), one final combined search from the round peptides, post-analysis with the peptide
 *   proteins remapped across rounds. Output under {@code seed_deep} ({@code seed/} holds the seed search,
 *   ranking, {@code roundN/} and {@code iterative_rounds_info.tsv}).</li>
 * </ul>
 * Deep per file is not supported (the GUI does not support it for Seed either).
 *
 * <p>Libraries: exactly one seed library, registered in the MAG database as Microbiome
 * (created in the GUI), plus any Standalone (host) libraries, which must be {@code .parquet} libraries. A library given as an existing file path is Standalone, as in the GUI.
 * The host can also be given as {@code host_database} with {@code append_host_db}, as in DIA Marker.
 *
 * <p>Host: in every mode the host library is searched together with the microbial library only in the final
 * combined search, and the export removes host peptides (host takes shared peptides); see
 * {@link DiaHostLibrarySupport}. The seed search, genome selection and Deep rounds use the seed library only.
 *
 * <p>Outputs (under {@code <output_directory>/seed_focus} or {@code seed_comprehensive}):
 * <ul>
 *   <li>seed/seed.parquet</li>
 *   <li>seed/genome_evidence.tsv, seed/genome_ranked.tsv, seed/genome_selected.tsv</li>
 *   <li>seed/combined.fasta, seed/combined.speclib (+ seed/combined.lib.parquet with a host library)</li>
 *   <li>combined/host_library.txt with a host library; host_peptides.tsv and host_summary.tsv</li>
 *   <li>combined/combined.parquet + combined.pr_matrix.tsv + combined.pg_matrix.tsv</li>
 *   <li>post-analysis TSVs and the HTML report</li>
 * </ul>
 */
public final class DiaSeedHeadlessRunner {

    private static final Logger LOGGER = LogManager.getLogger(DiaSeedHeadlessRunner.class);
    private static final String TASK_NAME = "DIA Seed";
    /** Deep together: as in DiaSeedDeepPerFileTask (used by DiaSeedDeepTogetherTask). */
    private static final int DEFAULT_GENOMES_PER_ROUND = 5;
    private static final int MAX_ITERATIVE_ROUNDS = 100;
    private static final double RELATIVE_GAIN_THRESHOLD = 0.015;
    private static final int CONSECUTIVE_ROUNDS_FOR_PLATEAU = 2;
    private static final String CUT = "KR";

    private DiaSeedHeadlessRunner() {
    }

    public static final class Result {
        public final File resultDir;
        public final File seedDir;
        public final File combinedDir;
        public final File seedParquet;
        public final File selectedTsv;
        public final File combinedFasta;
        public final File combinedSpeclib;
        public final File combinedParquet;
        public final File prMatrixTsv;
        public final File pgMatrixTsv;

        Result(File resultDir, File seedDir, File combinedDir, File seedParquet, File selectedTsv,
               File combinedFasta, File combinedSpeclib, File combinedParquet, File prMatrixTsv, File pgMatrixTsv) {
            this.resultDir = resultDir;
            this.seedDir = seedDir;
            this.combinedDir = combinedDir;
            this.seedParquet = seedParquet;
            this.selectedTsv = selectedTsv;
            this.combinedFasta = combinedFasta;
            this.combinedSpeclib = combinedSpeclib;
            this.combinedParquet = combinedParquet;
            this.prMatrixTsv = prMatrixTsv;
            this.pgMatrixTsv = pgMatrixTsv;
        }
    }

    public static final class SeedLibraries {
        public final File seedLibrary;
        /** Standalone (host) libraries; empty when only a seed library is used. */
        public final List<File> hostLibraries;

        SeedLibraries(File seedLibrary, List<File> hostLibraries) {
            this.seedLibrary = seedLibrary;
            this.hostLibraries = hostLibraries;
        }
    }

    /**
     * Run parameters for a DIA Seed parameter file, built as the GUI Seed panel builds them
     * for Core / Focus and Comprehensive runs.
     */
    public static MetaParameterDiaMarker toRunParameters(MetaParameterDia seedParams) throws IOException {
        if (seedParams == null) throw new IllegalArgumentException("seedParams is null");
        DiaRefinementMode mode = seedParams.getSeedRefinementMode();
        requireSupportedMode(mode);
        MetaParameterDiaMarker params = new MetaParameterDiaMarker(seedParams, seedParams.getMagDbItem(),
                false, false, seedParams.getLibrary(), false, new String[0]);
        params.setMarkerRefinementMode(mode);
        return params;
    }

    /**
     * Throws unless {@code mode} is one of the Seed modes this runner supports.
     */
    public static void requireSupportedMode(DiaRefinementMode mode) throws IOException {
        if (mode != DiaRefinementMode.DIRECT && mode != DiaRefinementMode.CORE
                && mode != DiaRefinementMode.COMPREHENSIVE && mode != DiaRefinementMode.DEEP_TOGETHER) {
            throw new IOException("DIA Seed in the CLI supports Direct, Core / Focus, Comprehensive and Deep together"
                    + " modes (got " + (mode == null ? "null" : mode.getDisplayLabel()) + ")");
        }
    }

    /**
     * Classifies the selected libraries as in the GUI ({@link MagPepLib#isMatchToMag()}) and
     * resolves their files.
     */
    public static SeedLibraries resolveLibraries(String[] libraries, MagDbItem magDb) throws IOException {
        if (magDb == null) {
            throw new IOException("MAG database is required to resolve DIA Seed libraries");
        }
        if (libraries == null || libraries.length == 0) {
            throw new IOException("No library selected (dia_settings.libraries). DIA Seed needs one seed library;"
                    + " create it in the GUI first.");
        }
        File seed = null;
        List<File> hosts = new ArrayList<>();
        for (String entry : libraries) {
            if (entry == null || entry.isBlank()) {
                throw new IOException("Empty library entry in dia_settings.libraries");
            }
            String name = entry.trim();
            File direct = new File(name);
            boolean microbiome;
            File libFile;
            if (direct.exists()) {
                microbiome = false;
                libFile = direct.getAbsoluteFile();
            } else {
                MagPepLib pepLib = findPepLib(magDb, name);
                if (pepLib == null) {
                    throw new IOException("Library '" + name + "' is neither an existing file nor a library of MAG database "
                            + magDb.getCurrentFile() + ". Create the seed library in the GUI first, and copy its folder "
                            + MagDbItem.LIBRARIES_FOLDER_NAME + "/" + name + " together with the MAG database.");
                }
                microbiome = pepLib.isMatchToMag();
                libFile = microbiome ? magDb.resolveSpectralLibraryFile(name) : registeredParquet(magDb, name);
                if (libFile == null || !libFile.exists()) {
                    throw new IOException("No library file found for library '" + name + "' under "
                            + new File(magDb.getLibrariesRootFolder(), name)
                            + (microbiome ? " (.speclib/.tsv/.parquet)" : " (" + name + ".parquet)"));
                }
            }
            if (microbiome) {
                if (seed != null) {
                    throw new IOException("More than one seed (Microbiome) library selected; the CLI accepts one");
                }
                seed = libFile;
            } else {
                hosts.add(libFile);
            }
        }
        if (seed == null) {
            throw new IOException("No seed library: one library must be registered as Microbiome in the MAG database."
                    + " Create the seed library in the GUI first.");
        }
        for (File host : hosts) {
            if (!host.getName().toLowerCase(Locale.ROOT).endsWith(".parquet")) {
                throw new IOException("A host (Standalone) library must be a .parquet library: " + host.getAbsolutePath());
            }
        }
        return new SeedLibraries(seed, hosts);
    }

    private static MagPepLib findPepLib(MagDbItem magDb, String name) {
        MagPepLib lib = registeredPepLib(magDb, name);
        if (lib != null) {
            return lib;
        }
        File root = magDb.getLibrariesRootFolder();
        if (root == null) {
            return null;
        }
        File folder = new File(root, name);
        if (!new File(folder, name + "_info.txt").isFile()) {
            return null;
        }
        magDb.addPepLibFromFile(folder);
        return registeredPepLib(magDb, name);
    }

    private static MagPepLib registeredPepLib(MagDbItem magDb, String name) {
        MagPepLib[] libs = magDb.getMagPepLibs();
        if (libs == null) {
            return null;
        }
        for (MagPepLib lib : libs) {
            if (name.equals(lib.getLibName())) {
                return lib;
            }
        }
        return null;
    }

    private static File registeredParquet(MagDbItem magDb, String name) {
        File root = magDb.getLibrariesRootFolder();
        return root == null ? null : new File(new File(root, name), name + ".parquet");
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
            MetaPilotResourceConfig reportResourceConfig,
            boolean exportMetaAtlasParquet
    ) throws Exception {
        if (cb == null) cb = new ProgressCallback.NoOp();
        if (params == null) throw new IllegalArgumentException("params is null");
        if (diannExe == null || !diannExe.exists()) throw new IOException("DIA-NN not found: " + (diannExe == null ? "null" : diannExe.getAbsolutePath()));

        DiaRefinementMode mode = params.getMarkerRefinementMode();
        requireSupportedMode(mode);
        if (mode == DiaRefinementMode.DIRECT) {
            return runDirect(params, diannExe, magDbConfigJson, cb, condaPrefix, funcDbFile, generateHtmlReport,
                    reportPythonHint, reportResourceConfig, exportMetaAtlasParquet);
        }
        if (mode == DiaRefinementMode.DEEP_TOGETHER) {
            return runDeepTogether(params, diannExe, magDbConfigJson, cb, condaPrefix, funcDbFile, generateHtmlReport,
                    reportPythonHint, reportResourceConfig, exportMetaAtlasParquet);
        }

        File studyRoot = new File(params.getResult());
        if (!studyRoot.exists()) studyRoot.mkdirs();
        if (!studyRoot.isDirectory()) throw new IOException("Result dir is not a directory: " + studyRoot.getAbsolutePath());

        File resultDir = new File(studyRoot, mode.usesComprehensiveSelection() ? "seed_comprehensive" : "seed_focus");
        if (!resultDir.exists()) resultDir.mkdirs();
        cb.onMessage(TASK_NAME + ": result folder -> " + resultDir.getAbsolutePath());

        String[] diaFiles = params.getMetadata() != null ? params.getMetadata().getRawFiles() : null;
        if (diaFiles == null || diaFiles.length == 0) {
            throw new IOException("No DIA files found in params metadata (input_output.experiments)");
        }

        MagDbItem magDb = DiaMarkerHeadlessRunner.resolveMagDb(params, magDbConfigJson);
        if (magDb == null) throw new IOException("MAG DB could not be resolved from params + magdb config");

        SeedLibraries libs = resolveLibraries(params.getLibrary(), magDb);
        cb.onMessage(TASK_NAME + ": seed library -> " + libs.seedLibrary.getAbsolutePath());
        for (File host : libs.hostLibraries) {
            cb.onMessage(TASK_NAME + ": host library (final combined search only) -> " + host.getAbsolutePath());
        }

        int threads = params.getThreadCount();
        File seedDir = new File(resultDir, "seed");
        if (!seedDir.exists()) seedDir.mkdirs();

        // Step 1: seed library search (no quant); the host library is not part of it.
        File seedParquet = runSeedLibrarySearch(params, libs, diaFiles, seedDir, threads, diannExe, cb, condaPrefix);

        // Steps 2-3: map precursors to MAG proteins, then select genomes.
        File evidenceTsv = new File(seedDir, "genome_evidence.tsv");
        File rankedTsv = new File(seedDir, "genome_ranked.tsv");
        File selectedTsv = new File(seedDir, "genome_selected.tsv");
        if (!selectedTsv.exists()) {
            selectGenomes(params, magDb, mode, seedParquet, evidenceTsv, rankedTsv, selectedTsv, cb);
        }
        HashSet<String> selectedGenomes = new HashSet<>();
        DiaMarkerHeadlessRunner.readSelectedGenomesIntoSet(selectedTsv, selectedGenomes);
        if (selectedGenomes.isEmpty()) {
            throw new IOException("No genomes selected (empty genome_selected.tsv): " + selectedTsv.getAbsolutePath());
        }

        // Step 4: combined FASTA (peptide-based or whole genomes).
        File combinedFasta = new File(seedDir, "combined.fasta");
        if (!combinedFasta.exists() || combinedFasta.length() == 0) {
            if (params.isUsePeptideBasedRound()) {
                double threshold = params.getPeptideDetectabilityThreshold();
                cb.onMessage(TASK_NAME + ": building peptide-based combined.fasta (threshold=" + threshold + ")");
                if (!DeepIterativeCommonSupport.writeCombinedPeptideFastaFromSelectedGenomes(
                        new File(magDb.getCurrentFile(), "predicted"), selectedGenomes, threshold, combinedFasta,
                        LOGGER, TASK_NAME)) {
                    throw new IOException("Failed to write peptide-based combined FASTA: " + combinedFasta.getAbsolutePath());
                }
            } else {
                cb.onMessage(TASK_NAME + ": building combined.fasta from " + selectedGenomes.size() + " genomes");
                DiaMarkerHeadlessRunner.buildCombinedFasta(magDb, selectedGenomes, combinedFasta);
            }
        }
        if (!combinedFasta.exists() || combinedFasta.length() == 0) {
            throw new IOException("combined.fasta was not created: " + combinedFasta.getAbsolutePath());
        }

        // Step 5: combined library. DIA-NN settings carry over from the seed search, as in the GUI.
        File combinedSpeclib = new File(seedDir, "combined.speclib");
        if (!combinedSpeclib.exists()) {
            File alt = new File(seedDir, "combined.predicted.speclib");
            if (!alt.exists()) {
                cb.onMessage(TASK_NAME + ": generating combined.speclib from combined.fasta");
                DiannParameter diannPar = seedSearchDiannPar(params, threads);
                if (params.isUsePeptideBasedRound()) {
                    // peptide-based combined.fasta: one entry per predicted peptide, no protein N-terminal rules
                    uk.ac.quadram.metapilot.core.dia.PeptideListLibrarySupport.applyPeptideEntryRules(diannPar);
                }
                String cmd = diannPar.generateFasta(combinedFasta.getAbsolutePath(), combinedSpeclib.getAbsolutePath(), false);
                int exit = DiaNNHeadlessRunner.run(diannExe, cmd, seedDir, cb, condaPrefix);
                if (exit != 0) LOGGER.warn("DIA-NN combined library build exit code {}", exit);
            }
            if (!combinedSpeclib.exists() && alt.exists()) {
                Files.move(alt.toPath(), combinedSpeclib.toPath());
            }
        }
        if (!combinedSpeclib.exists()) throw new IOException("combined.speclib not created: " + combinedSpeclib.getAbsolutePath());

        // Step 6: quant search into combined/.
        File combinedDir = new File(resultDir, "combined");
        if (!combinedDir.exists()) combinedDir.mkdirs();
        File combinedParquet = new File(combinedDir, "combined.parquet");
        File prMatrix = new File(combinedDir, "combined.pr_matrix.tsv");
        File pgMatrix = new File(combinedDir, "combined.pg_matrix.tsv");

        List<File> hostParquets = resolveHostParquets(libs, params, resultDir, diannExe, cb, condaPrefix);
        String[] finalLibs = DiaHostLibrarySupport.finalSearchLibraries(combinedSpeclib, hostParquets, combinedDir,
                DiaMarkerHeadlessRunner.libraryBuildPar(params), diannExe, cb, condaPrefix);

        if (!prMatrix.exists() || !combinedParquet.exists()) {
            cb.onMessage(TASK_NAME + ": running DIA-NN quant with the combined library"
                    + (hostParquets.isEmpty() ? "" : " and the host library"));
            DiaLibSearchPar libSearchPar = new DiaLibSearchPar();
            libSearchPar.setMbr(true);
            libSearchPar.setDiaFiles(diaFiles);
            libSearchPar.setThreadCount(threads);
            libSearchPar.setqValue(0.01);
            String cmd = quantDiannPar(params, threads).generateLibSearch(libSearchPar,
                    combinedParquet.getAbsolutePath(), finalLibs);
            cmd = cmd + DiaMarkerHeadlessRunner.buildFArgs(diaFiles);
            int exit = DiaNNHeadlessRunner.run(diannExe, cmd, combinedDir, cb, condaPrefix);
            if (exit != 0) LOGGER.warn("DIA-NN combined quant exit code {}", exit);
        }
        if (!combinedParquet.exists()) throw new IOException("Combined quant did not produce: " + combinedParquet.getAbsolutePath());
        if (!prMatrix.exists()) throw new IOException("Combined quant did not produce: " + prMatrix.getAbsolutePath());
        if (!pgMatrix.exists()) throw new IOException("Combined quant did not produce: " + pgMatrix.getAbsolutePath());

        // Step 7: post-analysis + report under seed_focus / seed_comprehensive.
        File gtdbNcbi = null;
        File maybe = new File("Resources", "gtdb_ncbi");
        if (maybe.exists()) gtdbNcbi = maybe;
        DiaMagPostAnalysisTsvExporter.exportSeed(resultDir, magDb, combinedParquet, prMatrix,
                gtdbNcbi, funcDbFile, params.getMicroDb(), params.getQuantitativeMode(), cb,
                selectedGenomes, exportMetaAtlasParquet);
        runHtmlReportIfRequested(resultDir, params, generateHtmlReport, reportPythonHint, reportResourceConfig, cb);

        cb.onComplete("DIA seed complete: " + resultDir.getAbsolutePath());
        return new Result(resultDir, seedDir, combinedDir, seedParquet, selectedTsv, combinedFasta, combinedSpeclib,
                combinedParquet, prMatrix, pgMatrix);
    }

    /**
     * Seed library search without quantification into {@code seedDir/seed.parquet} (skipped when it exists).
     * Shared by Core / Focus, Comprehensive and Deep together, which all start with this search in the GUI
     * ({@code DiaSeedCoreTask.runFirstSearchUserLibrary}).
     */
    private static File runSeedLibrarySearch(MetaParameterDiaMarker params, SeedLibraries libs, String[] diaFiles,
                                             File seedDir, int threads, File diannExe, ProgressCallback cb,
                                             File condaPrefix) throws Exception {
        File seedParquet = new File(seedDir, "seed.parquet");
        if (!seedParquet.exists()) {
            DiannParameter diannPar = DiaMarkerHeadlessRunner.buildDiannPar(params);
            DiaLibSearchPar libSearchPar = new DiaLibSearchPar();
            DiaNNNoQuantLibrarySearchService.configureNoQuantSearch(diannPar, libSearchPar, diaFiles, threads, true, true);
            cb.onMessage(TASK_NAME + ": searching seed library against DIA files");
            String cmd = diannPar.generateLibSearchFastLib(libSearchPar, seedParquet.getAbsolutePath(),
                    libs.seedLibrary.getAbsolutePath());
            cmd = cmd + DiaMarkerHeadlessRunner.buildFArgs(diaFiles);
            int exit = DiaNNHeadlessRunner.run(diannExe, cmd, seedDir, cb, condaPrefix);
            if (exit != 0) LOGGER.warn("DIA-NN seed library search exit code {}", exit);
        }
        if (!seedParquet.exists()) throw new IOException("Seed library search did not produce: " + seedParquet.getAbsolutePath());
        return seedParquet;
    }

    private static File gtdbNcbiMapping() {
        File maybe = new File("Resources", "gtdb_ncbi");
        return maybe.exists() ? maybe : null;
    }

    /**
     * Seed Direct ({@code DiaSeedDirectTask}): one quantifying search of the DIA files against the seed library
     * into {@code seed_direct/combined/}, then post-analysis. The seed library proteins are not MAG proteins; the
     * exporter maps the seed peptides to MAG proteins via {@code predicted/} before genome and protein scoring.
     */
    private static Result runDirect(MetaParameterDiaMarker params, File diannExe, File magDbConfigJson,
                                    ProgressCallback cb, File condaPrefix, File funcDbFile, boolean generateHtmlReport,
                                    String reportPythonHint, MetaPilotResourceConfig reportResourceConfig,
                                    boolean exportMetaAtlasParquet) throws Exception {
        final String task = TASK_NAME + " Direct";
        File studyRoot = new File(params.getResult());
        if (!studyRoot.exists()) studyRoot.mkdirs();
        File resultDir = new File(studyRoot, "seed_direct");
        if (!resultDir.exists()) resultDir.mkdirs();
        cb.onMessage(task + ": result folder -> " + resultDir.getAbsolutePath());

        String[] diaFiles = params.getMetadata() != null ? params.getMetadata().getRawFiles() : null;
        if (diaFiles == null || diaFiles.length == 0) {
            throw new IOException("No DIA files found in params metadata (input_output.experiments)");
        }
        MagDbItem magDb = DiaMarkerHeadlessRunner.resolveMagDb(params, magDbConfigJson);
        if (magDb == null) throw new IOException("MAG DB could not be resolved from params + magdb config");
        SeedLibraries libs = resolveLibraries(params.getLibrary(), magDb);
        cb.onMessage(task + ": seed library -> " + libs.seedLibrary.getAbsolutePath());

        int threads = params.getThreadCount();
        File combinedDir = new File(resultDir, "combined");
        if (!combinedDir.exists()) combinedDir.mkdirs();
        File combinedParquet = new File(combinedDir, "combined.parquet");
        File prMatrix = new File(combinedDir, "combined.pr_matrix.tsv");
        File pgMatrix = new File(combinedDir, "combined.pg_matrix.tsv");
        // Direct has a single search, which is the final combined search: the host library joins it here.
        List<File> hostParquets = resolveHostParquets(libs, params, resultDir, diannExe, cb, condaPrefix);
        File seedForSearch = libs.seedLibrary;
        if (!hostParquets.isEmpty()) {
            // multi-library search needs .parquet: the library folder holds <name>.parquet next to <name>.speclib
            seedForSearch = DiaHostLibrarySupport.parquetInLibraryFolder(seedForSearch);
            if (!seedForSearch.getName().toLowerCase(Locale.ROOT).endsWith(".parquet")) {
                seedForSearch = DiaHostLibrarySupport.ensureParquetLibrary(seedForSearch, combinedDir, "seed",
                        DiaMarkerHeadlessRunner.libraryBuildPar(params), diannExe, cb, condaPrefix);
            }
        }
        String[] directLibs = DiaHostLibrarySupport.finalSearchLibraries(seedForSearch, hostParquets, combinedDir,
                DiaMarkerHeadlessRunner.libraryBuildPar(params), diannExe, cb, condaPrefix);
        if (!prMatrix.exists() || !combinedParquet.exists()) {
            cb.onMessage(task + ": searching the seed library against DIA files (quantification)"
                    + (hostParquets.isEmpty() ? "" : " together with the host library"));
            DiaLibSearchPar libSearchPar = new DiaLibSearchPar();
            libSearchPar.setMbr(true);
            libSearchPar.setDiaFiles(diaFiles);
            libSearchPar.setThreadCount(threads);
            libSearchPar.setqValue(0.01);
            String cmd = quantDiannPar(params, threads).generateLibSearch(libSearchPar,
                    combinedParquet.getAbsolutePath(), directLibs);
            cmd = cmd + DiaMarkerHeadlessRunner.buildFArgs(diaFiles);
            int exit = DiaNNHeadlessRunner.run(diannExe, cmd, combinedDir, cb, condaPrefix);
            if (exit != 0) LOGGER.warn("DIA-NN seed direct search exit code {}", exit);
        }
        if (!combinedParquet.exists()) throw new IOException("Seed direct search did not produce: " + combinedParquet.getAbsolutePath());
        if (!prMatrix.exists()) throw new IOException("Seed direct search did not produce: " + prMatrix.getAbsolutePath());
        if (!pgMatrix.exists()) throw new IOException("Seed direct search did not produce: " + pgMatrix.getAbsolutePath());

        DiaMagPostAnalysisTsvExporter.exportSeed(resultDir, magDb, combinedParquet, prMatrix,
                gtdbNcbiMapping(), funcDbFile, params.getMicroDb(), params.getQuantitativeMode(), cb,
                null, exportMetaAtlasParquet);
        runHtmlReportIfRequested(resultDir, params, generateHtmlReport, reportPythonHint, reportResourceConfig, cb);

        cb.onComplete("DIA seed direct complete: " + resultDir.getAbsolutePath());
        return new Result(resultDir, null, combinedDir, null, null, null, null, combinedParquet, prMatrix, pgMatrix);
    }

    /**
     * Seed Deep together ({@code DiaSeedDeepTogetherTask} / {@code DiaSeedDeepPerFileTask} in DEEP_TOGETHER mode):
     * seed library search -> mapping to MAG proteins -> comprehensive genome ranking -> iterative rounds across all
     * DIA files -> final combined search from the round peptides -> post-analysis.
     */
    private static Result runDeepTogether(MetaParameterDiaMarker params, File diannExe, File magDbConfigJson,
                                          ProgressCallback cb, File condaPrefix, File funcDbFile,
                                          boolean generateHtmlReport, String reportPythonHint,
                                          MetaPilotResourceConfig reportResourceConfig,
                                          boolean exportMetaAtlasParquet) throws Exception {
        final String task = TASK_NAME + " Deep";
        final ProgressCallback callback = cb;
        File studyRoot = new File(params.getResult());
        if (!studyRoot.exists()) studyRoot.mkdirs();
        File resultDir = new File(studyRoot, "seed_deep");
        if (!resultDir.exists()) resultDir.mkdirs();
        cb.onMessage(task + ": result folder -> " + resultDir.getAbsolutePath());

        String[] diaFiles = params.getMetadata() != null ? params.getMetadata().getRawFiles() : null;
        if (diaFiles == null || diaFiles.length == 0) {
            throw new IOException("No DIA files found in params metadata (input_output.experiments)");
        }
        MagDbItem magDb = DiaMarkerHeadlessRunner.resolveMagDb(params, magDbConfigJson);
        if (magDb == null) throw new IOException("MAG DB could not be resolved from params + magdb config");
        SeedLibraries libs = resolveLibraries(params.getLibrary(), magDb);
        cb.onMessage(task + ": seed library -> " + libs.seedLibrary.getAbsolutePath());

        int threads = params.getThreadCount();
        File seedDir = new File(resultDir, "seed");
        if (!seedDir.exists()) seedDir.mkdirs();
        File seedParquet = runSeedLibrarySearch(params, libs, diaFiles, seedDir, threads, diannExe, cb, condaPrefix);

        File evidenceTsv = new File(seedDir, "genome_evidence.tsv");
        File rankedTsv = new File(seedDir, "genome_ranked.tsv");
        File selectedTsv = new File(seedDir, "genome_selected.tsv");
        if (!rankedTsv.exists() || !selectedTsv.exists()) {
            selectGenomesDeep(magDb, seedParquet, evidenceTsv, rankedTsv, selectedTsv, cb);
        }
        if (!rankedTsv.exists()) throw new IOException("Genome ranked TSV missing: " + rankedTsv.getAbsolutePath());

        int genomesPerRound = params.getGenomesPerRound() > 0 ? params.getGenomesPerRound() : DEFAULT_GENOMES_PER_ROUND;
        cb.onMessage(task + ": iterative rounds of " + genomesPerRound + " genomes");
        DeepIterativeCommonSupport.IterativeRoundLoopResult loop = DeepIterativeCommonSupport.runIterativeRoundLoop(
                rankedTsv, genomesPerRound, MAX_ITERATIVE_ROUNDS, "Deep seed round", CUT,
                callback::onMessage,
                (round, genomes) -> runSeedDeepRound(round, seedDir, genomes, diaFiles, magDb, diannExe, params,
                        condaPrefix, threads, callback),
                DiaSeedHeadlessRunner::computeDeepStopRule,
                LOGGER, task);
        if (loop.lastRound <= 0) throw new IOException(task + ": no iterative round completed");
        DeepIterativeCommonSupport.writeIterativeRoundsInfo(seedDir, loop.roundsRecord, loop.peptidesInRoundRecord,
                loop.incrementsPerRound, loop.cumulativeRecord, loop.stopReason, loop.lastStopRuleResult,
                CONSECUTIVE_ROUNDS_FOR_PLATEAU, LOGGER, task);

        File combinedDir = new File(resultDir, "combined");
        if (!combinedDir.exists()) combinedDir.mkdirs();
        File combinedFasta = new File(combinedDir, "combined.fasta");
        File combinedSpeclib = new File(combinedDir, "combined.speclib");
        File combinedParquet = new File(combinedDir, "combined.parquet");
        File prMatrix = new File(combinedDir, "combined.pr_matrix.tsv");
        File pgMatrix = new File(combinedDir, "combined.pg_matrix.tsv");
        cb.onMessage(task + ": final combined search from the peptides of " + loop.lastRound + " rounds");
        List<File> hostParquets = resolveHostParquets(libs, params, resultDir, diannExe, cb, condaPrefix);
        if (!DiaMarkerDeepTogetherHeadlessRunner.runGlobalCombinedFromRounds(diaFiles, seedDir, loop.lastRound,
                combinedDir, combinedFasta, combinedSpeclib, combinedParquet, hostParquets, diannExe, params,
                condaPrefix, cb)) {
            throw new IOException(task + ": final combined search failed");
        }
        if (!combinedParquet.exists()) throw new IOException("Combined quant did not produce: " + combinedParquet.getAbsolutePath());
        if (!prMatrix.exists()) throw new IOException("Combined quant did not produce: " + prMatrix.getAbsolutePath());
        if (!pgMatrix.exists()) throw new IOException("Combined quant did not produce: " + pgMatrix.getAbsolutePath());

        DiaMagPostAnalysisTsvExporter.export(resultDir, magDb, combinedParquet, prMatrix, gtdbNcbiMapping(),
                params.getMarkerMode(), funcDbFile, params.getMicroDb(), params.getQuantitativeMode(), cb,
                DiaDeepProteinMappingService.Request.deepTogether(seedDir, loop.lastRound, diaFiles, CUT),
                exportMetaAtlasParquet);
        runHtmlReportIfRequested(resultDir, params, generateHtmlReport, reportPythonHint, reportResourceConfig, cb);

        cb.onComplete("DIA seed deep complete: " + resultDir.getAbsolutePath());
        return new Result(resultDir, seedDir, combinedDir, seedParquet, selectedTsv, combinedFasta, combinedSpeclib,
                combinedParquet, prMatrix, pgMatrix);
    }

    /** Comprehensive genome ranking for Seed Deep, as in DiaSeedDeepPerFileTask.iterativeSearch(). */
    private static void selectGenomesDeep(MagDbItem magDb, File seedParquet, File evidenceTsv, File rankedTsv,
                                          File selectedTsv, ProgressCallback cb) throws Exception {
        cb.onMessage(TASK_NAME + " Deep: reading precursors (q<=0.01) and mapping them to MAG proteins via predicted TSVs");
        DiaNNPrecursor[] precursors = new DiaNNParquetReader(seedParquet.getAbsolutePath(), CUT).getDiaNNPrecursors(0.01);
        if (precursors == null || precursors.length == 0) {
            throw new IOException("No precursors found in seed library search: " + seedParquet.getAbsolutePath());
        }
        File predictedFolder = new File(magDb.getCurrentFile(), "predicted");
        if (!predictedFolder.isDirectory()) {
            throw new IOException("MAG predicted folder not found: " + predictedFolder.getAbsolutePath());
        }
        List<DiaNNPrecursor> withPro = MagPredictedPeptideProteinMapper.mapPrecursorsToProteins(precursors, predictedFolder);
        if (withPro == null || withPro.isEmpty()) {
            throw new IOException("No precursors could be mapped to MAG proteins; check predicted TSVs under "
                    + predictedFolder.getAbsolutePath());
        }
        cb.onMessage(TASK_NAME + " Deep: mapped " + withPro.size() + " / " + precursors.length + " precursors to MAG proteins");
        Marker4ResourcePreparer.ensureSelectionResources(magDb, cb);
        File markerFasta = new File(magDb.getCurrentFile(), "Marker4.fasta");
        File genomesTaxaTsv = new File(magDb.getCurrentFile(), "genomes-taxa.tsv");
        if (!markerFasta.exists() || !genomesTaxaTsv.exists()) {
            throw new IOException("Marker4.fasta or genomes-taxa.tsv missing under: " + magDb.getCurrentFile().getAbsolutePath());
        }
        File effectiveMarkerFasta = markerFastaForSelection(withPro, markerFasta, TASK_NAME + " Deep", cb);
        MagGenomeSelectionService.Request request = new MagGenomeSelectionService.Request();
        request.magDb = magDb;
        request.genomeSelectorMode = MagDbItem.MARKER4;
        request.genomeSelectionStrategy = MagGenomeSelectionService.SELECTION_STRATEGY_COMPREHENSIVE;
        request.diaPrecursors = withPro.toArray(new DiaNNPrecursor[0]);
        request.markerFasta = effectiveMarkerFasta;
        request.disableMarkerFasta = effectiveMarkerFasta == null;
        request.genomesTaxaTsv = genomesTaxaTsv;
        request.evidenceFile = evidenceTsv;
        request.rankedFile = rankedTsv;
        request.selectedFile = selectedTsv;
        request.maximumSelected = 300;
        MagGenomeSelectionService.selectGenomes(request, cb);
    }

    /**
     * One Seed Deep round, as DiaSeedDeepPerFileTask.runSeedRoundAppendPreviousPeptides (rounds without
     * --reanalyse; round library built with the settings in force after the seed search).
     */
    private static int runSeedDeepRound(int round, File seedDir, List<String> genomes, String[] diaFiles,
                                        MagDbItem magDb, File diannExe, MetaParameterDiaMarker params,
                                        File condaPrefix, int threads, ProgressCallback cb) {
        final String task = TASK_NAME + " Deep";
        return DeepIterativeCommonSupport.runAppendPreviousPeptidesRound(round, seedDir, genomes, diaFiles, false,
                "deep seed",
                (gs, fastaOut) -> DiaMarkerDeepTogetherHeadlessRunner.writeFastaFromGenomes(magDb, gs, fastaOut),
                (gs, prevParquet, fastaOut) -> DiaMarkerDeepTogetherHeadlessRunner
                        .writeFastaFromGenomesAndPreviousPeptides(magDb, gs, prevParquet, fastaOut),
                (combinedFasta, combinedSpeclib) -> {
                    try {
                        String cmd = seedSearchDiannPar(params, threads).generateFasta(combinedFasta.getAbsolutePath(),
                                combinedSpeclib.getAbsolutePath(), false);
                        int exit = DiaNNHeadlessRunner.run(diannExe, cmd, combinedFasta.getParentFile(), cb, condaPrefix);
                        if (exit != 0) LOGGER.warn("{}: DIA-NN round speclib build exit code {}", task, exit);
                        File alt = new File(combinedFasta.getParentFile(), "combined.predicted.speclib");
                        if (!combinedSpeclib.exists() && alt.exists()) {
                            Files.move(alt.toPath(), combinedSpeclib.toPath());
                        }
                    } catch (Exception e) {
                        LOGGER.error("{}: round {} build speclib failed", task, round, e);
                    }
                },
                (roundParquet, combinedSpeclibPath, files, reanalyse) -> {
                    try {
                        DiaNNNoQuantLibrarySearchService.run(diannExe, params, files, new File(combinedSpeclibPath),
                                roundParquet, threads, reanalyse, false, cb, condaPrefix, LOGGER, task);
                    } catch (Exception e) {
                        LOGGER.error("{}: round {} no-quant search failed", task, round, e);
                    }
                },
                (parquetFile) -> DeepIterativeCommonSupport.countIdentifiedPeptides(parquetFile, CUT, LOGGER, task),
                LOGGER, task);
    }

    /** Seed Deep stop rule: relative gain < 1.5% in 2 consecutive rounds, profile DEEP_COMPREHENSIVE (as the GUI). */
    private static DeepIterativeStopRule.StopRuleResult computeDeepStopRule(List<Integer> increments) {
        return DeepIterativeStopRule.computeStopRoundFromIncrements(increments, RELATIVE_GAIN_THRESHOLD,
                CONSECUTIVE_ROUNDS_FOR_PLATEAU, DeepIterativeStopRule.AlgorithmProfile.DEEP_COMPREHENSIVE);
    }

    private static void selectGenomes(MetaParameterDiaMarker params, MagDbItem magDb, DiaRefinementMode mode,
                                      File seedParquet, File evidenceTsv, File rankedTsv, File selectedTsv,
                                      ProgressCallback cb) throws Exception {
        cb.onMessage(TASK_NAME + ": reading precursors and mapping them to MAG proteins via predicted TSVs");
        DiaNNParquetReader parquetReader = new DiaNNParquetReader(seedParquet.getAbsolutePath(),
                DiaMarkerHeadlessRunner.buildDiannPar(params).getCut());
        DiaNNPrecursor[] precursors = parquetReader.getDiaNNPrecursors();
        if (precursors == null || precursors.length == 0) {
            throw new IOException("No precursors found in seed library search: " + seedParquet.getAbsolutePath());
        }
        File predictedFolder = new File(magDb.getCurrentFile(), "predicted");
        if (!predictedFolder.isDirectory()) {
            throw new IOException("MAG predicted folder not found: " + predictedFolder.getAbsolutePath());
        }
        List<DiaNNPrecursor> precursorsWithPro = MagPredictedPeptideProteinMapper.mapPrecursorsToProteins(precursors,
                predictedFolder);
        if (precursorsWithPro == null || precursorsWithPro.isEmpty()) {
            throw new IOException("No precursors could be mapped to MAG proteins; check predicted TSVs under "
                    + predictedFolder.getAbsolutePath());
        }
        cb.onMessage(TASK_NAME + ": mapped " + precursorsWithPro.size() + " / " + precursors.length
                + " precursors to MAG proteins");

        Marker4ResourcePreparer.ensureSelectionResources(magDb, cb);
        File markerFasta = new File(magDb.getCurrentFile(), "Marker4.fasta");
        File genomesTaxaTsv = new File(magDb.getCurrentFile(), "genomes-taxa.tsv");
        if (!markerFasta.exists() || !genomesTaxaTsv.exists()) {
            throw new IOException("Marker4.fasta or genomes-taxa.tsv missing under: " + magDb.getCurrentFile().getAbsolutePath());
        }

        File effectiveMarkerFasta = markerFastaForSelection(precursorsWithPro, markerFasta, TASK_NAME, cb);

        String strategy = mode.usesComprehensiveSelection()
                ? MarkerGenomeSelector.KNEE_MODE_COMPREHENSIVE
                : MarkerGenomeSelector.KNEE_MODE_FOCUS;
        cb.onMessage(TASK_NAME + ": selecting genomes (mode=" + strategy + ")");
        MagGenomeSelectionService.Request request = new MagGenomeSelectionService.Request();
        request.magDb = magDb;
        request.genomeSelectorMode = MagDbItem.MARKER4;
        request.genomeSelectionStrategy = strategy;
        request.diaPrecursors = precursorsWithPro.toArray(new DiaNNPrecursor[0]);
        request.markerFasta = effectiveMarkerFasta;
        request.disableMarkerFasta = effectiveMarkerFasta == null;
        request.genomesTaxaTsv = genomesTaxaTsv;
        request.evidenceFile = evidenceTsv;
        request.rankedFile = rankedTsv;
        request.selectedFile = selectedTsv;
        request.maximumSelected = 300;
        MagGenomeSelectionService.selectGenomes(request, cb);
    }

    /**
     * Host libraries of the final combined search: the Standalone libraries of {@code dia_settings.libraries} plus
     * {@code host_database} with {@code append_host_db} (as in DIA Marker), without duplicates; empty without a host.
     */
    private static List<File> resolveHostParquets(SeedLibraries libs, MetaParameterDiaMarker params, File resultDir,
                                                  File diannExe, ProgressCallback cb, File condaPrefix) throws Exception {
        List<File> hosts = new ArrayList<>(libs.hostLibraries);
        File fromParams = DiaMarkerHeadlessRunner.resolveHostParquet(params, new File(resultDir, "host_library"),
                diannExe, cb, condaPrefix);
        if (fromParams != null) {
            boolean listed = false;
            for (File h : hosts) listed |= h.getCanonicalFile().equals(fromParams.getCanonicalFile());
            if (!listed) hosts.add(fromParams);
        }
        return hosts;
    }

    /**
     * A seed library need not contain marker proteins: below the marker peptide ratio threshold, genomes are selected
     * from all peptides (greedy only, returns null); otherwise from marker peptides (returns {@code markerFasta}).
     */
    private static File markerFastaForSelection(List<DiaNNPrecursor> precursorsWithPro, File markerFasta, String label,
                                                ProgressCallback cb) {
        MarkerGenomeSelector.MarkerPeptideRatioResult ratioResult =
                MarkerGenomeSelector.computeMarkerPeptideRatio(precursorsWithPro, markerFasta);
        String ratioPercent = String.format(Locale.ROOT, "%.1f", ratioResult.ratio * 100);
        if (ratioResult.ratio < MarkerGenomeSelector.MARKER_PEPTIDE_RATIO_THRESHOLD) {
            cb.onMessage(label + ": marker peptides (" + ratioResult.markerUniquePeptides + ") are " + ratioPercent
                    + "% of identified peptides (" + ratioResult.totalUniquePeptides
                    + "); using all peptides for genome selection (greedy only)");
            return null;
        }
        cb.onMessage(label + ": marker peptide ratio " + ratioResult.markerUniquePeptides + " / "
                + ratioResult.totalUniquePeptides + " = " + ratioPercent + "% (marker-based genome selection)");
        return markerFasta;
    }

    /** DIA-NN settings in force after the seed search (no-quant search configuration). */
    private static DiannParameter seedSearchDiannPar(MetaParameterDiaMarker params, int threads) {
        DiannParameter diannPar = DiaMarkerHeadlessRunner.buildDiannPar(params);
        diannPar.setThreads(threads);
        diannPar.setMatrices(false);
        diannPar.setRelaxed_prot_inf(true);
        diannPar.setReanalyse(true);
        diannPar.setExportDecoy(true);
        return diannPar;
    }

    /** DIA-NN settings of the final quant search. */
    private static DiannParameter quantDiannPar(MetaParameterDiaMarker params, int threads) {
        DiannParameter diannPar = DiaMarkerHeadlessRunner.buildDiannPar(params);
        diannPar.setThreads(threads);
        diannPar.setMatrices(true);
        diannPar.setRelaxed_prot_inf(true);
        diannPar.setReanalyse(true);
        diannPar.setExportDecoy(true);
        return diannPar;
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
                cb,
                params != null ? params.getMetadata() : null
        );
        if (!report.warnings.isEmpty()) {
            cb.onMessage("WARNING: HTML report generation completed with warnings: "
                    + String.join("; ", report.warnings));
        }
    }
}
