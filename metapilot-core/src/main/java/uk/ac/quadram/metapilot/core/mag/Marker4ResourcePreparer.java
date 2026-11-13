package uk.ac.quadram.metapilot.core.mag;

import uk.ac.quadram.metapilot.core.MGYG.HapExtractor;
import uk.ac.quadram.metapilot.core.deepdetect.DeepDetectPythonRunner;
import uk.ac.quadram.metapilot.core.deepdetect.DeepDetectWindowsRunner;
import uk.ac.quadram.metapilot.core.resources.MetaPilotResourceConfig;
import uk.ac.quadram.metapilot.core.resources.ToolKey;
import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;
import uk.ac.quadram.metapilot.dbSearch.deepDetect.DeepDetectParameter;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;
import uk.ac.quadram.metapilot.task.mag.MarkerGenomeSelector;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.LinkedHashMap;
import java.util.Locale;
import java.util.Map;
import java.util.Optional;
import java.util.concurrent.TimeUnit;

/**
 * Prepares Marker4 resource files needed by genome selection and DIA Marker library build.
 *
 * <p>Selection resources ({@link #ensureSelectionResources}) ensure {@code Marker4.fasta} and
 * {@code genomes-taxa.tsv}. Library resources ({@link #ensureLibraryResources}) additionally
 * enforce the GUI highpep rule: catalogues with more than {@link #HIGHPEP_PROTEIN_THRESHOLD}
 * proteins require DeepDetect and {@code Marker4_highpep.fasta} at detectability
 * {@link #HIGHPEP_DETECTABILITY_THRESHOLD}.
 */
public final class Marker4ResourcePreparer {

    public static final String MARKER4_FASTA_NAME = "Marker4.fasta";
    public static final String MARKER4_HIGHPEP_FASTA_NAME = "Marker4_highpep.fasta";
    public static final String MARKER4_PREDICTED_TSV_NAME = "Marker4.predicted.tsv";
    public static final String GENOMES_TAXA_TSV_NAME = "genomes-taxa.tsv";

    /** Same threshold as the GUI Marker4 library path. */
    public static final int HIGHPEP_PROTEIN_THRESHOLD = 100_000;

    /** Same export threshold as the GUI Marker4 library path. */
    public static final double HIGHPEP_DETECTABILITY_THRESHOLD = 0.8;

    private static final int DEEPDETECT_TIMEOUT_HOURS = 4;

    /**
     * Runs DeepDetect prediction for a single FASTA → predicted TSV.
     */
    @FunctionalInterface
    public interface DeepDetectPredictor {
        void predict(File inputFasta, File outputPredictedTsv, ProgressCallback cb) throws IOException;
    }

    private Marker4ResourcePreparer() {
    }

    /**
     * Ensures {@code Marker4.fasta} and {@code genomes-taxa.tsv} exist (selection path).
     * Does not require DeepDetect / highpep.
     */
    public static SelectionResources ensureSelectionResources(MagDbItem magDb, ProgressCallback cb) throws IOException {
        if (cb == null) {
            cb = new ProgressCallback.NoOp();
        }
        if (magDb == null) {
            throw new IOException("MAG DB is null");
        }
        File currentFile = magDb.getCurrentFile();
        if (currentFile == null || !currentFile.exists()) {
            throw new IOException("MAG DB current folder missing");
        }

        File marker4Fasta = new File(currentFile, MARKER4_FASTA_NAME);
        ensureMarker4Fasta(magDb, marker4Fasta, cb);

        File genomesTaxaTsv = new File(currentFile, GENOMES_TAXA_TSV_NAME);
        ensureGenomesTaxaTsv(magDb, genomesTaxaTsv, cb);

        return new SelectionResources(marker4Fasta, genomesTaxaTsv);
    }

    /**
     * Same as {@link #ensureSelectionResources}, then applies the DIA Marker highpep rule
     * used by the GUI: if {@code Marker4.fasta} has more than
     * {@link #HIGHPEP_PROTEIN_THRESHOLD} proteins, ensure {@code Marker4_highpep.fasta}
     * exists (DeepDetect + export at {@link #HIGHPEP_DETECTABILITY_THRESHOLD}).
     *
     * <p>Fails with a clear {@link IOException} when the catalogue is above the threshold,
     * highpep is missing, and DeepDetect is not configured or fails — no silent fallback
     * to the full FASTA.
     *
     * @param predictor may be null only when highpep is already present or protein count ≤ threshold
     */
    public static SelectionResources ensureLibraryResources(MagDbItem magDb, ProgressCallback cb,
            DeepDetectPredictor predictor) throws IOException {
        SelectionResources base = ensureSelectionResources(magDb, cb);
        ensureHighpepIfNeeded(base.marker4Fasta, magDb.getCurrentFile(), cb, predictor);
        return base;
    }

    /**
     * Build a DeepDetect predictor from CLI/headless resource config.
     * Returns {@code null} when DEEPDETECT is not configured or missing on disk.
     */
    public static DeepDetectPredictor predictorFromResources(MetaPilotResourceConfig resources,
            ProgressCallback cb) {
        if (resources == null) {
            return null;
        }
        ProgressCallback progress = cb == null ? new ProgressCallback.NoOp() : cb;
        Optional<Path> deepDetectPath = resources.resolveToolPath(ToolKey.DEEPDETECT, progress);
        if (deepDetectPath.isEmpty() || !Files.exists(deepDetectPath.get())) {
            return null;
        }
        File deepDetect = deepDetectPath.get().toFile();
        File python = resources.resolveToolPath(ToolKey.PYTHON, progress).map(Path::toFile).orElse(null);
        return predictorFromTool(deepDetect, python);
    }

    /**
     * Build a DeepDetect predictor from an explicit tool path (GUI {@code DeepDetect.exe}
     * or HPC {@code main.py}).
     */
    public static DeepDetectPredictor predictorFromTool(File deepDetectTool, File pythonExe) {
        if (deepDetectTool == null || !deepDetectTool.exists()) {
            return null;
        }
        final File tool = deepDetectTool;
        final File python = pythonExe;
        return (inputFasta, outputPredictedTsv, cb) -> runDeepDetect(tool, python, inputFasta, outputPredictedTsv, cb);
    }

    public static void ensureHighpepIfNeeded(File marker4Fasta, File currentFile, ProgressCallback cb,
            DeepDetectPredictor predictor) throws IOException {
        if (cb == null) {
            cb = new ProgressCallback.NoOp();
        }
        if (marker4Fasta == null || !marker4Fasta.exists()) {
            throw new IOException("Marker4.fasta missing; cannot evaluate highpep requirement");
        }
        if (currentFile == null) {
            throw new IOException("MAG DB current folder is null");
        }

        int proteinCount = countFastaProteins(marker4Fasta);
        if (proteinCount < 0) {
            throw new IOException("Failed to read protein count from " + marker4Fasta.getAbsolutePath());
        }

        File highpep = new File(currentFile, MARKER4_HIGHPEP_FASTA_NAME);
        if (proteinCount <= HIGHPEP_PROTEIN_THRESHOLD) {
            cb.onMessage(String.format(Locale.ROOT,
                    "Marker4.fasta has %,d proteins (≤%,d); using full FASTA for library",
                    proteinCount, HIGHPEP_PROTEIN_THRESHOLD));
            return;
        }

        cb.onMessage(String.format(Locale.ROOT,
                "Marker4.fasta has %,d proteins (>%s); Marker4_highpep.fasta required",
                proteinCount, String.format(Locale.ROOT, "%,d", HIGHPEP_PROTEIN_THRESHOLD)));

        if (highpep.exists() && highpep.length() > 0) {
            cb.onMessage("Marker4_highpep.fasta already present: " + highpep.getAbsolutePath());
            return;
        }

        if (predictor == null) {
            throw new IOException(String.format(Locale.ROOT,
                    "Marker4.fasta has %,d proteins (>%s) but DeepDetect is not configured. "
                            + "Configure DEEPDETECT (and PYTHON when using main.py) and retry, "
                            + "or pre-generate %s. Refusing to build the Marker4 library from the full FASTA.",
                    proteinCount,
                    String.format(Locale.ROOT, "%,d", HIGHPEP_PROTEIN_THRESHOLD),
                    MARKER4_HIGHPEP_FASTA_NAME));
        }

        File predictedTsv = new File(currentFile, MARKER4_PREDICTED_TSV_NAME);
        if (!predictedTsv.exists() || predictedTsv.length() == 0) {
            cb.onMessage("Running DeepDetect for Marker4 peptide detectability → " + predictedTsv.getName());
            predictor.predict(marker4Fasta, predictedTsv, cb);
            if (!predictedTsv.exists() || predictedTsv.length() == 0) {
                throw new IOException("DeepDetect did not produce " + predictedTsv.getAbsolutePath());
            }
        } else {
            cb.onMessage("Reusing existing " + predictedTsv.getName());
        }

        cb.onMessage(String.format(Locale.ROOT,
                "Exporting peptides with detectability > %.1f → %s",
                HIGHPEP_DETECTABILITY_THRESHOLD, MARKER4_HIGHPEP_FASTA_NAME));
        MarkerGenomeSelector.exportPeptideFasta(predictedTsv.getAbsolutePath(), highpep.getAbsolutePath(),
                HIGHPEP_DETECTABILITY_THRESHOLD);
        if (!highpep.exists() || highpep.length() == 0) {
            throw new IOException("Failed to write " + highpep.getAbsolutePath()
                    + " (no peptides above detectability " + HIGHPEP_DETECTABILITY_THRESHOLD + "?)");
        }
        cb.onMessage("Generated " + highpep.getAbsolutePath());
    }

    private static void runDeepDetect(File deepDetectTool, File pythonExe, File inputFasta,
            File outputPredictedTsv, ProgressCallback cb) throws IOException {
        ProgressCallback progress = cb == null ? new ProgressCallback.NoOp() : cb;
        DeepDetectParameter par = new DeepDetectParameter();
        par.setInput(inputFasta);
        par.setOutput(outputPredictedTsv);

        File outDir = outputPredictedTsv.getParentFile();
        if (outDir != null && !outDir.exists() && !outDir.mkdirs()) {
            throw new IOException("Could not create DeepDetect output folder: " + outDir.getAbsolutePath());
        }

        try {
            if (isPythonEntrypoint(deepDetectTool)) {
                if (pythonExe == null || !pythonExe.exists()) {
                    throw new IOException("PYTHON executable is required when DEEPDETECT points to "
                            + deepDetectTool.getName() + " (" + deepDetectTool.getAbsolutePath() + ")");
                }
                int exit = DeepDetectPythonRunner.run(pythonExe, deepDetectTool, par, null, progress);
                if (exit != 0) {
                    throw new IOException("DeepDetect (python) exited with code " + exit);
                }
                return;
            }

            File workDir = outDir != null ? outDir : inputFasta.getParentFile();
            File prepared = DeepDetectPythonRunner.ensureFastaHeadersHaveProteinName(inputFasta, workDir);
            if (prepared != null) {
                par.setInput(prepared);
            }
            File batFile = new File(workDir != null ? workDir : new File(System.getProperty("java.io.tmpdir")),
                    "Marker4_deepDetect.bat");
            DeepDetectWindowsRunner.createBatchFile(deepDetectTool, par, batFile);
            progress.onMessage("Running DeepDetect (Windows): " + batFile.getAbsolutePath());

            ProcessBuilder pb = new ProcessBuilder("cmd.exe", "/c", batFile.getAbsolutePath());
            pb.redirectErrorStream(true);
            Process process = pb.start();
            try (BufferedReader reader = new BufferedReader(
                    new java.io.InputStreamReader(process.getInputStream()))) {
                String line;
                while ((line = reader.readLine()) != null) {
                    if (Thread.currentThread().isInterrupted()) {
                        process.destroyForcibly();
                        throw new IOException("Interrupted while running DeepDetect");
                    }
                    String t = line.trim();
                    if (!t.isEmpty()) {
                        progress.onMessage(t);
                    }
                }
            }
            boolean finished = process.waitFor(DEEPDETECT_TIMEOUT_HOURS, TimeUnit.HOURS);
            if (!finished) {
                process.destroyForcibly();
                throw new IOException("DeepDetect timed out after " + DEEPDETECT_TIMEOUT_HOURS + " hour(s)");
            }
            int exit = process.exitValue();
            try {
                Files.deleteIfExists(batFile.toPath());
            } catch (Exception ignored) {
                // best-effort cleanup
            }
            if (exit != 0) {
                throw new IOException("DeepDetect (Windows) exited with code " + exit);
            }
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            throw new IOException("Interrupted while running DeepDetect", e);
        }
    }

    private static boolean isPythonEntrypoint(File tool) {
        if (tool == null) {
            return false;
        }
        String name = tool.getName().toLowerCase(Locale.ROOT);
        return name.endsWith(".py") || name.equals("main.py") || name.contains("deepdetect_pred");
    }

    public static int countFastaProteins(File fastaFile) {
        if (fastaFile == null || !fastaFile.exists()) {
            return -1;
        }
        int count = 0;
        try (BufferedReader r = new BufferedReader(new FileReader(fastaFile))) {
            String line;
            while ((line = r.readLine()) != null) {
                if (line.startsWith(">")) {
                    count++;
                }
            }
        } catch (IOException e) {
            return -1;
        }
        return count;
    }

    public static File ensureGenomesTaxaTsv(MagDbItem magDb, File genomeTaxaFile, ProgressCallback cb)
            throws IOException {
        if (cb == null) {
            cb = new ProgressCallback.NoOp();
        }
        if (genomeTaxaFile != null && genomeTaxaFile.exists() && genomeTaxaFile.length() > 0
                && !startsWithHeader(genomeTaxaFile)) {
            return genomeTaxaFile;
        }
        if (genomeTaxaFile == null) {
            throw new IOException("genomes-taxa.tsv path is null");
        }
        File taxa = magDb != null ? magDb.getOriginalTaxaFile() : null;
        if (taxa == null || !taxa.exists()) {
            throw new IOException("Original taxa file not found for MAG DB; cannot generate genomes-taxa.tsv");
        }

        cb.onMessage("genomes-taxa.tsv missing or outdated; generating from MAG taxa file");
        writeGenomesTaxaTsv(taxa, genomeTaxaFile);
        if (!genomeTaxaFile.exists() || genomeTaxaFile.length() == 0) {
            throw new IOException("genomes-taxa.tsv was not created: " + genomeTaxaFile.getAbsolutePath());
        }
        return genomeTaxaFile;
    }

    private static void ensureMarker4Fasta(MagDbItem magDb, File marker4Fasta, ProgressCallback cb)
            throws IOException {
        if (marker4Fasta.exists() && marker4Fasta.length() > 0) {
            return;
        }
        cb.onMessage("Marker4.fasta missing; extracting Marker4 proteins from original DB");
        File originalDbFolder = magDb.getOriginalDbFolder();
        if (originalDbFolder == null || !originalDbFolder.exists()) {
            throw new IOException("Original DB folder not available to extract Marker4.fasta");
        }
        HapExtractor extractor = new HapExtractor();
        extractor.extractFourCategories(originalDbFolder, marker4Fasta);
        if (!marker4Fasta.exists() || marker4Fasta.length() == 0) {
            throw new IOException("Marker4.fasta was not created: " + marker4Fasta.getAbsolutePath());
        }
    }

    private static void writeGenomesTaxaTsv(File originalTaxaFile, File outputFile) throws IOException {
        File parent = outputFile.getParentFile();
        if (parent != null && !parent.exists() && !parent.mkdirs()) {
            throw new IOException("Could not create folder for genomes-taxa.tsv: " + parent.getAbsolutePath());
        }

        try (BufferedReader reader = new BufferedReader(new FileReader(originalTaxaFile))) {
            String header = reader.readLine();
            if (header == null) {
                throw new IOException("Empty taxa file: " + originalTaxaFile.getAbsolutePath());
            }
            String[] heads = header.split("\t", -1);
            int genomeId = findColumn(heads, "Genome", "Rep_ID", "Rep", "Species_rep");
            int lineageId = findColumn(heads, "Lineage", "GTDB_lineage");
            if (genomeId < 0) {
                throw new IOException("Cannot generate genomes-taxa.tsv: Genome column not found in "
                        + originalTaxaFile.getAbsolutePath());
            }
            if (lineageId < 0) {
                throw new IOException("Cannot generate genomes-taxa.tsv: Lineage column not found in "
                        + originalTaxaFile.getAbsolutePath());
            }

            Map<String, String> genomeToLineage = new LinkedHashMap<>();
            String line;
            while ((line = reader.readLine()) != null) {
                if (line.trim().isEmpty()) {
                    continue;
                }
                String[] cols = line.split("\t", -1);
                if (cols.length <= Math.max(genomeId, lineageId)) {
                    continue;
                }
                String genome = cols[genomeId].trim();
                String lineage = cols[lineageId].trim();
                if (!genome.isEmpty() && !lineage.isEmpty() && !genomeToLineage.containsKey(genome)) {
                    genomeToLineage.put(genome, lineage);
                }
            }

            try (PrintWriter writer = new PrintWriter(outputFile)) {
                for (Map.Entry<String, String> entry : genomeToLineage.entrySet()) {
                    writer.println(entry.getKey() + "\t" + entry.getValue());
                }
            }
        }
    }

    private static int findColumn(String[] headers, String... names) {
        for (String name : names) {
            for (int i = 0; i < headers.length; i++) {
                if (name.equalsIgnoreCase(headers[i].trim())) {
                    return i;
                }
            }
        }
        return -1;
    }

    private static boolean startsWithHeader(File file) throws IOException {
        try (BufferedReader reader = new BufferedReader(new FileReader(file))) {
            String line;
            while ((line = reader.readLine()) != null) {
                line = line.trim();
                if (!line.isEmpty()) {
                    String[] cols = line.split("\t", -1);
                    return cols.length >= 2
                            && "Genome".equalsIgnoreCase(cols[0].trim())
                            && "Lineage".equalsIgnoreCase(cols[1].trim());
                }
            }
        }
        return false;
    }

    public static final class SelectionResources {
        public final File marker4Fasta;
        public final File genomesTaxaTsv;

        private SelectionResources(File marker4Fasta, File genomesTaxaTsv) {
            this.marker4Fasta = marker4Fasta;
            this.genomesTaxaTsv = genomesTaxaTsv;
        }
    }
}
