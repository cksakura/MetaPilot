package uk.ac.quadram.metapilot.core.runner;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import uk.ac.quadram.metapilot.core.dia.Marker4LibrarySupportHeadless;
import uk.ac.quadram.metapilot.core.diann.DiaNNHeadlessRunner;
import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;
import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNParquetReader;
import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNPrecursor;
import uk.ac.quadram.metapilot.dbSearch.diann.DiannParameter;
import uk.ac.quadram.metapilot.task.dia.DiaLibSearchPar;
import uk.ac.quadram.metapilot.task.dia.par.MetaParameterDiaMarker;
import uk.ac.quadram.metapilot.task.mag.MagDbConfigIO;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;
import uk.ac.quadram.metapilot.task.mag.MarkerGenomeSelector;
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
 * This runner intentionally avoids HTML report generation (HPC-friendly).
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
        if (cb == null) cb = new ProgressCallback.NoOp();
        if (params == null) throw new IllegalArgumentException("params is null");
        if (diannExe == null || !diannExe.exists()) throw new IOException("DIA-NN not found: " + (diannExe == null ? "null" : diannExe.getAbsolutePath()));

        // Dispatch to the per-file runner when markerMode == DEEP_PER_FILE so the CLI/HPC entrypoint
        // matches the GUI's DEEP_PER_FILE workflow (per-raw-file iterative rounds + global combined search).
        // Returns are adapted to this runner's Result type so existing callers keep working.
        if ("DEEP_PER_FILE".equalsIgnoreCase(params.getMarkerMode())) {
            cb.onMessage("DIA Marker: dispatching to DEEP_PER_FILE headless runner");
            DiaDeepPerFileMarkerHeadlessRunner.Result perFileResult =
                    DiaDeepPerFileMarkerHeadlessRunner.run(params, diannExe, magDbConfigJson, cb, condaPrefix, funcDbFile);
            return new Result(
                    perFileResult.resultDir,
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
        String markerMode = params.getMarkerMode();
        String workFolderName;
        if ("COMPREHENSIVE".equalsIgnoreCase(markerMode) || "DEEP".equalsIgnoreCase(markerMode)) {
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
            MarkerGenomeSelector selector = new MarkerGenomeSelector();
            selector.setKneeSelectionMode(kneeMode);
            selector.selectMarkerGenomesWithKnee(
                    precursors,
                    markerFasta,
                    genomesTaxaTsv,
                    evidenceTsv.getAbsolutePath(),
                    rankedTsv.getAbsolutePath(),
                    selectedTsv.getAbsolutePath(),
                    300,
                    true
            );
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

        // Post-analysis + TSV exports (no HTML) under resultDir (marker_focus / marker_comprehensive),
        // matching the Windows task layout.
        try {
            File gtdbNcbi = null;
            File maybe = new File("Resources", "gtdb_ncbi");
            if (maybe.exists()) gtdbNcbi = maybe;
            File rankedForLabels = rankedTsv != null && rankedTsv.exists() ? rankedTsv : null;
            DiaMarkerPostAnalysisTsvExporter.export(
                    resultDir,
                    magDb,
                    combinedParquet,
                    prMatrix,
                    rankedForLabels,
                    gtdbNcbi,
                    params.getMarkerMode(),
                    funcDbFile,
                    params.getMicroDb(),
                    params.getQuantitativeMode(),
                    cb
            );
        } catch (Exception e) {
            cb.onMessage("WARNING: DIA marker post-analysis export failed: " + e.getMessage());
            throw e;
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
        if (magDbConfigJson == null || !magDbConfigJson.exists()) {
            return fromParams;
        }
        MagDbItem[][] parsed = MagDbConfigIO.parse(magDbConfigJson);
        if (parsed == null) return fromParams;
        java.util.ArrayList<MagDbItem> items = new java.util.ArrayList<>();
        for (MagDbItem[] group : parsed) {
            if (group == null) continue;
            for (MagDbItem it : group) {
                if (it != null) items.add(it);
            }
        }
        if (items.isEmpty()) return fromParams;

        // Prefer matching by catalogue+identifier from params, otherwise fall back to first item.
        String cat = fromParams != null ? fromParams.getCatalogueID() : null;
        String id = fromParams != null ? fromParams.getIdentifier() : null;
        if (cat != null && id != null) {
            for (MagDbItem it : items) {
                if (cat.equals(it.getCatalogueID()) && id.equals(it.getIdentifier())) {
                    if (fromParams != null && fromParams.getUsedVersion() != null) {
                        it.setUsedVersion(fromParams.getUsedVersion());
                    }
                    return it;
                }
            }
        }
        return items.get(0);
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

        Marker4LibrarySupportHeadless.ensureMarker4LibraryMetadata(
                diannExe,
                condaPrefix,
                magDb,
                params,
                cb,
                LOGGER
        );
    }

    static void ensureMarker4Files(MagDbItem magDb, File diannExe, MetaParameterDiaMarker params, ProgressCallback cb) throws Exception {
        if (magDb == null) throw new IOException("MAG DB is null");
        File currentFile = magDb.getCurrentFile();
        if (currentFile == null || !currentFile.exists()) throw new IOException("MAG DB current folder missing");

        File marker4Fasta = new File(currentFile, "Marker4.fasta");
        if (!marker4Fasta.exists()) {
            cb.onMessage("Marker4.fasta missing; extracting Marker4 proteins from original DB");
            File originalDbFolder = magDb.getOriginalDbFolder();
            if (originalDbFolder == null || !originalDbFolder.exists()) {
                throw new IOException("Original DB folder not available to extract Marker4.fasta");
            }
            HapExtractor extractor = new HapExtractor();
            extractor.extractFourCategories(originalDbFolder, marker4Fasta);
            if (!marker4Fasta.exists()) throw new IOException("Marker4.fasta was not created: " + marker4Fasta.getAbsolutePath());
        }

        // Ensure genomes-taxa.tsv exists (needed by MarkerGenomeSelector).
        File genomeTaxa = new File(currentFile, "genomes-taxa.tsv");
        if (!genomeTaxa.exists()) {
            cb.onMessage("genomes-taxa.tsv missing; generating from MAG taxa file");
            generateGenomesTaxaTsv(magDb, genomeTaxa);
        }
        if (!genomeTaxa.exists()) throw new IOException("genomes-taxa.tsv was not created: " + genomeTaxa.getAbsolutePath());
    }

    private static void generateGenomesTaxaTsv(MagDbItem magDb, File genomeTaxaFile) throws IOException {
        File taxa = magDb.getOriginalTaxaFile();
        if (taxa == null || !taxa.exists()) {
            throw new IOException("Original taxa file not found for MAG DB; cannot generate genomes-taxa.tsv");
        }
        try (BufferedReader reader = new BufferedReader(new FileReader(taxa))) {
            String header = reader.readLine();
            if (header == null) throw new IOException("Empty taxa file: " + taxa.getAbsolutePath());
            String[] heads = header.split("\t");
            int repId = -1, lineageId = -1;
            for (int i = 0; i < heads.length; i++) {
                String h = heads[i];
                if ("Rep_ID".equalsIgnoreCase(h) || "Rep".equalsIgnoreCase(h) || "Genome".equalsIgnoreCase(h)) repId = i;
                if ("Lineage".equalsIgnoreCase(h) || "GTDB_lineage".equalsIgnoreCase(h)) lineageId = i;
            }
            if (repId < 0) repId = 0;
            if (lineageId < 0) lineageId = Math.min(1, heads.length - 1);

            HashMap<String, String> map = new HashMap<>();
            String line;
            while ((line = reader.readLine()) != null) {
                String[] cs = line.split("\t");
                if (cs.length > Math.max(repId, lineageId) && !map.containsKey(cs[repId])) {
                    map.put(cs[repId], cs[lineageId]);
                }
            }
            try (PrintWriter writer = new PrintWriter(genomeTaxaFile)) {
                writer.println("Genome\tLineage");
                for (String key : map.keySet()) {
                    writer.println(key + "\t" + map.get(key));
                }
            }
        }
    }

    private static void runMarkerNoQuantSearch(File diannExe, MetaParameterDiaMarker params, String[] diaFiles,
                                              File markerLibFile, File outParquet, ProgressCallback cb, File condaPrefix) throws Exception {
        if (markerLibFile == null || !markerLibFile.exists()) {
            throw new IOException("Marker library not found: " + (markerLibFile == null ? "null" : markerLibFile.getAbsolutePath()));
        }
        DiannParameter diannPar = buildDiannPar(params);
        diannPar.setMatrices(false);
        diannPar.setRelaxed_prot_inf(true);
        diannPar.setReanalyse(true);
        diannPar.setExportDecoy(true);

        DiaLibSearchPar libSearchPar = new DiaLibSearchPar();
        libSearchPar.setMbr(false);
        libSearchPar.setDiaFiles(diaFiles);
        libSearchPar.setThreadCount(params.getThreadCount());
        libSearchPar.setqValue(0.01);

        String cmd = diannPar.generateLibSearchNoQuant(libSearchPar, outParquet.getAbsolutePath(), markerLibFile.getAbsolutePath());
        cmd = cmd + buildFArgs(diaFiles);
        int exit = DiaNNHeadlessRunner.run(diannExe, cmd, outParquet.getParentFile(), cb, condaPrefix);
        if (exit != 0) LOGGER.warn("DIA-NN marker no-quant exit code {}", exit);
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
        if (markerMode == null) return MarkerGenomeSelector.KNEE_MODE_FOCUS;
        if ("COMPREHENSIVE".equalsIgnoreCase(markerMode) || "DEEP".equalsIgnoreCase(markerMode)
                || "DEEP_PER_FILE".equalsIgnoreCase(markerMode)) {
            return MarkerGenomeSelector.KNEE_MODE_COMPREHENSIVE;
        }
        return MarkerGenomeSelector.KNEE_MODE_FOCUS;
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

