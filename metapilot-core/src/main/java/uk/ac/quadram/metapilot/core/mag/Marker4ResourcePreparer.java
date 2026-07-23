package uk.ac.quadram.metapilot.core.mag;

import uk.ac.quadram.metapilot.core.MGYG.HapExtractor;
import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 * Prepares the small Marker4 resource files needed by genome selection.
 */
public final class Marker4ResourcePreparer {

    public static final String MARKER4_FASTA_NAME = "Marker4.fasta";
    public static final String GENOMES_TAXA_TSV_NAME = "genomes-taxa.tsv";

    private Marker4ResourcePreparer() {
    }

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
