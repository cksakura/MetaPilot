package uk.ac.quadram.metapilot.core.runner;

import uk.ac.quadram.metapilot.core.fragpipe.FragPipeHeadlessRunner;
import uk.ac.quadram.metapilot.core.fragpipe.FragpipeWorkflowFile;
import uk.ac.quadram.metapilot.core.prodb.DecoyShuffler;
import uk.ac.quadram.metapilot.core.prodb.FastaManager;
import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;
import uk.ac.quadram.metapilot.dbSearch.fragpipe.FragPipePSM;
import uk.ac.quadram.metapilot.dbSearch.fragpipe.FragpipePSMTsvReader;
import uk.ac.quadram.metapilot.dbSearch.fragpipe.FragpipePepReader;
import uk.ac.quadram.metapilot.dbSearch.fragpipe.FragpipePeptide;
import uk.ac.quadram.metapilot.task.dda.mag.par.DDAMagFragPipeParameters;
import uk.ac.quadram.metapilot.task.mag.MagDbConfigIO;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;
import uk.ac.quadram.metapilot.task.mag.MarkerGenomeSelector;
import uk.ac.quadram.metapilot.task.par.MetaData;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Headless implementation of the MAG+FragPipe workflow (first CLI-ready slice).
 *
 * Implements:
 * - HAP search (Marker2/Marker4) per sample
 * - genome selection (MarkerGenomeSelector) and sample-specific FASTA creation
 * - sample-specific search per sample
 * - combined search for quant outputs (combined_modified_peptide.tsv / combined_protein.tsv)
 *
 * Skips GUI-only report generation and downstream visualization steps.
 */
public final class MagFragPipeHeadlessRunner {

    private MagFragPipeHeadlessRunner() {
    }

    public static void run(
            DDAMagFragPipeParameters params,
            File resourcesFragpipe,
            File workflowTemplate,
            File magdbConfigJson,
            ProgressCallback cb
    ) throws Exception {
        if (cb == null) cb = new ProgressCallback.NoOp();
        if (params == null) throw new IllegalArgumentException("params is null");

        MetaData md = params.getMetadata();
        if (md == null || md.getRawFiles() == null || md.getRawFiles().length == 0) {
            throw new IllegalArgumentException("No raw files in parameters metadata");
        }

        // Resolve MAG DB item from MAGDB_1_0.json using catalogueId + version
        MagDbItem mag = resolveMagDb(params.getMagDbCatalogueId(), params.getMagDbVersion(), magdbConfigJson);
        if (mag == null) {
            throw new IllegalArgumentException("MAG database not found in " + magdbConfigJson + " for catalogue="
                    + params.getMagDbCatalogueId() + " version=" + params.getMagDbVersion());
        }

        cb.onMessage("Using MAG DB: " + mag.getCatalogueID() + " / " + mag.getUsedVersion().getCatalogueVersion());
        cb.onMessage("MAG DB folder: " + (mag.getCurrentFile() != null ? mag.getCurrentFile().getAbsolutePath() : "null"));

        // Prepare output folders
        File resultRoot = new File(params.getResult(), "fragpipe_MAG_search");
        resultRoot.mkdirs();
        cb.onMessage("Result root: " + resultRoot.getAbsolutePath());

        String[] rawFiles = md.getRawFiles();
        String[] expNames = md.getExpNames();
        int threads = Math.max(1, params.getThreadCount());

        // Parse template workflow once; we will write per-step workflow files with updated DB.
        FragpipeWorkflowFile wfTemplate = FragpipeWorkflowFile.parse(workflowTemplate.toPath());

        // Step A: per-sample HAP searches
        cb.onMessage("Step A: running HAP searches (per sample) ...");
        for (int i = 0; i < rawFiles.length; i++) {
            File sampleDir = new File(resultRoot, baseNameNoExt(rawFiles[i]));
            File hapDir = new File(sampleDir, "hap");
            hapDir.mkdirs();

            File manifest = new File(sampleDir, "raws.fp-manifest");
            if (!manifest.exists()) {
                writeManifest(manifest, rawFiles[i], safeExpName(expNames, i));
            }

            String genomeSelectorMode = params.getGenomeSelectorMode();
            File hapFasta = mag.getHapFastaFile(genomeSelectorMode);
            File hapTdFasta = mag.getHapFastaTDFile(genomeSelectorMode);
            if (hapFasta == null || !hapFasta.exists()) {
                throw new IOException("HAP FASTA not found for selectorMode=" + genomeSelectorMode + ": " + hapFasta);
            }
            if (hapTdFasta == null) {
                throw new IOException("HAP TD FASTA path is null (MAG DB not configured correctly)");
            }
            if (!hapTdFasta.exists()) {
                cb.onMessage("Creating target-decoy FASTA: " + hapTdFasta.getAbsolutePath());
                FastaManager.addDecoy(hapFasta.getAbsolutePath(), hapTdFasta.getAbsolutePath());
            }
            if (!hapTdFasta.exists()) {
                throw new IOException("Failed to create target-decoy FASTA: " + hapTdFasta.getAbsolutePath());
            }

            File wfOut = new File(hapDir, "hap.workflow");
            wfTemplate.setDatabase(hapTdFasta.getAbsolutePath()).write(wfOut.toPath());

            cb.onMessage("HAP search [" + (i + 1) + "/" + rawFiles.length + "]: " + rawFiles[i]);
            int exit = FragPipeHeadlessRunner.run(resourcesFragpipe, wfOut, manifest, hapDir, threads, cb);
            if (exit != 0) {
                throw new IOException("FragPipe HAP search failed (exit=" + exit + ") for " + rawFiles[i]);
            }
        }

        // Step B: select genomes and build sample-specific FASTAs
        cb.onMessage("Step B: selecting genomes and building sample-specific FASTAs ...");
        for (int i = 0; i < rawFiles.length; i++) {
            File sampleDir = new File(resultRoot, baseNameNoExt(rawFiles[i]));
            File hapDir = new File(sampleDir, "hap");

            // FragPipe output folder name: expName_1 (same as GUI: replace '-' with '_' + "_1")
            String expFolder = safeExpName(expNames, i).replaceAll("-", "_") + "_1";
            File fpOut = new File(hapDir, expFolder);
            File psmTsv = new File(fpOut, "psm.tsv");
            if (!psmTsv.exists()) {
                throw new IOException("Missing FragPipe psm.tsv after HAP search: " + psmTsv.getAbsolutePath());
            }

            boolean useMarker4 = MagDbItem.MARKER4.equalsIgnoreCase(params.getGenomeSelectorMode());
            String prefix = useMarker4 ? "marker_genome" : "genome";
            File evidence = new File(hapDir, prefix + "_evidence.tsv");
            File ranked = new File(hapDir, prefix + "_ranked.tsv");
            File selected = new File(hapDir, prefix + "_selected.tsv");

            HashSet<String> genomeSet = selectGenomes(mag, params, psmTsv, evidence, ranked, selected, cb);

            File sampleFasta = new File(hapDir, sampleDir.getName() + ".fasta");
            if (!sampleFasta.exists() || sampleFasta.length() == 0) {
                cb.onMessage("Writing sample-specific FASTA: " + sampleFasta.getAbsolutePath() + " (" + genomeSet.size() + " genomes)");
                try (PrintWriter writer = new PrintWriter(sampleFasta)) {
                    File originalDbFolder = mag.getOriginalDbFolder();
                    if (originalDbFolder == null || !originalDbFolder.exists()) {
                        throw new IOException("Original DB folder not found in MAG DB: " + originalDbFolder);
                    }
                    for (String genome : genomeSet) {
                        File faa = new File(originalDbFolder, genome + ".faa");
                        File fa = new File(originalDbFolder, genome + ".fasta");
                        File src = faa.exists() ? faa : (fa.exists() ? fa : null);
                        if (src == null) continue;
                        appendFastaWithDecoys(src, writer);
                    }
                }
            }
            if (!sampleFasta.exists() || sampleFasta.length() == 0) {
                throw new IOException("Sample-specific FASTA is empty: " + sampleFasta.getAbsolutePath());
            }
        }

        // Step C: sample-specific searches (per sample)
        cb.onMessage("Step C: running sample-specific searches (per sample) ...");
        for (int i = 0; i < rawFiles.length; i++) {
            File sampleDir = new File(resultRoot, baseNameNoExt(rawFiles[i]));
            File hapDir = new File(sampleDir, "hap");
            File sampleFasta = new File(hapDir, sampleDir.getName() + ".fasta");

            File sampleSearchDir = new File(sampleDir, "sample_db");
            sampleSearchDir.mkdirs();

            File manifest = new File(sampleDir, "raws.fp-manifest");
            if (!manifest.exists()) {
                writeManifest(manifest, rawFiles[i], safeExpName(expNames, i));
            }

            File wfOut = new File(sampleSearchDir, "sample_db.workflow");
            wfTemplate.setDatabase(sampleFasta.getAbsolutePath()).write(wfOut.toPath());

            cb.onMessage("Sample-specific search [" + (i + 1) + "/" + rawFiles.length + "]: " + rawFiles[i]);
            int exit = FragPipeHeadlessRunner.run(resourcesFragpipe, wfOut, manifest, sampleSearchDir, threads, cb);
            if (exit != 0) {
                throw new IOException("FragPipe sample-specific search failed (exit=" + exit + ") for " + rawFiles[i]);
            }
        }

        // Step D/E/F (CLI slice): build combined FASTA and run combined search for quant outputs.
        cb.onMessage("Step D/E/F: building combined FASTA and running combined search ...");
        File combinedDir = new File(resultRoot, "combined");
        combinedDir.mkdirs();

        File combinedFasta = new File(combinedDir, "combined.fasta");
        if (!combinedFasta.exists() || combinedFasta.length() == 0) {
            cb.onMessage("Building combined FASTA: " + combinedFasta.getAbsolutePath());
            HashMap<String, String> fastaMap = new HashMap<>();

            for (int i = 0; i < rawFiles.length; i++) {
                File sampleDir = new File(resultRoot, baseNameNoExt(rawFiles[i]));
                File sampleSearchDir = new File(sampleDir, "sample_db");
                String expFolder = safeExpName(expNames, i).replaceAll("-", "_") + "_1";
                File fpOut = new File(sampleSearchDir, expFolder);
                File proteinFas = new File(fpOut, "protein.fas");
                if (!proteinFas.exists()) {
                    throw new IOException("Missing FragPipe protein.fas for combined FASTA: " + proteinFas.getAbsolutePath());
                }
                try (BufferedReader reader = new BufferedReader(new FileReader(proteinFas))) {
                    String line;
                    while ((line = reader.readLine()) != null) {
                        if (!line.startsWith(">")) continue;
                        int space = line.indexOf(' ');
                        String header = (space > 0) ? line.substring(0, space) : line.trim();
                        String seq = reader.readLine();
                        if (seq == null) break;
                        fastaMap.put(header, seq.trim());
                    }
                }
            }

            String cleavage = getCleavageResiduesOrDefault(wfTemplate);
            try (PrintWriter writer = new PrintWriter(combinedFasta)) {
                for (String proName : fastaMap.keySet()) {
                    String seq = fastaMap.get(proName);
                    writer.println(proName);
                    writer.println(seq);

                    String shuffleProName = ">shuf_" + proName.substring(1);
                    String shuffle = DecoyShuffler.generateShuffledDecoy(seq, cleavage, 42);
                    writer.println(shuffleProName);
                    writer.println(shuffle);

                    String revProName = ">rev_" + proName.substring(1);
                    writer.println(revProName);
                    writer.println(new StringBuilder(seq).reverse().toString());

                    String revShufProName = ">rev_shuf_" + proName.substring(1);
                    writer.println(revShufProName);
                    writer.println(new StringBuilder(shuffle).reverse().toString());
                }

                // Optional host DB append
                if (params.isAppendHostDb() && params.getHostDB() != null && !params.getHostDB().isBlank()) {
                    File hostFasta = new File(params.getHostDB());
                    if (hostFasta.exists()) {
                        cb.onMessage("Appending host DB to combined FASTA: " + hostFasta.getAbsolutePath());
                        appendHostFastaWithDecoys(hostFasta, writer, cleavage);
                    } else {
                        cb.onMessage("WARNING: host DB enabled but not found: " + hostFasta.getAbsolutePath());
                    }
                }
            }
        }

        File combinedManifest = new File(combinedDir, "raws.fp-manifest");
        if (!combinedManifest.exists()) {
            // Reuse the same metadata that the GUI uses to build the manifest.
            params.getMetadata().exportFragpipeManifest(combinedManifest);
        }

        File combinedWorkflow = new File(combinedDir, "combined.workflow");
        wfTemplate.setDatabase(combinedFasta.getAbsolutePath()).write(combinedWorkflow.toPath());

        File quanPep = new File(combinedDir, "combined_modified_peptide.tsv");
        File quanPro = new File(combinedDir, "combined_protein.tsv");
        if (!quanPep.exists() || !quanPro.exists() || quanPep.length() == 0 || quanPro.length() == 0) {
            cb.onMessage("Running combined FragPipe search (quant outputs) ...");
            int exit = FragPipeHeadlessRunner.run(resourcesFragpipe, combinedWorkflow, combinedManifest, combinedDir, threads, cb);
            if (exit != 0) {
                throw new IOException("FragPipe combined search failed (exit=" + exit + ")");
            }
        } else {
            cb.onMessage("Combined quant outputs already exist: " + quanPep.getName() + ", " + quanPro.getName());
        }

        // Step G (headless): post-analysis and TSV exports (no HTML)
        cb.onMessage("Step G: post-analysis + TSV export (no HTML) ...");
        FragpipePostAnalysisTsvExporter.Result post = FragpipePostAnalysisTsvExporter.exportFinalTsvs(
                mag,
                quanPep,
                resultRoot,
                "DDA MAG FragPipe",
                cb
        );

        // Step H (headless): taxonomy TSV exports (no HTML)
        cb.onMessage("Step H: taxonomy TSV export (no HTML) ...");
        exportTaxonomyTsvs(params, quanPep, mag, resultRoot, post, cb);

        cb.onMessage("MAG+FragPipe headless run completed.");
        cb.onComplete(resultRoot.getAbsolutePath());
    }

    private static MagDbItem resolveMagDb(String catalogueId, String version, File magdbConfigJson) {
        if (catalogueId == null || catalogueId.isBlank()) return null;
        if (magdbConfigJson == null || !magdbConfigJson.exists()) return null;

        MagDbItem[][] all = MagDbConfigIO.parse(magdbConfigJson);
        if (all == null) return null;
        for (MagDbItem[] group : all) {
            if (group == null) continue;
            for (MagDbItem item : group) {
                if (item == null) continue;
                if (!catalogueId.equals(item.getCatalogueID())) continue;
                // Find matching version in versionList
                if (version != null && !version.isBlank()) {
                    for (MagDbItem.MagVersion mv : item.getVersionList()) {
                        if (mv != null && version.equals(mv.getCatalogueVersion())) {
                            item.setUsedVersion(mv);
                            return item;
                        }
                    }
                }
                // fallback: use the default localAvailable version
                MagDbItem.MagVersion used = item.getUsedVersion();
                if (used != null) {
                    item.setUsedVersion(used);
                    return item;
                }
            }
        }
        return null;
    }

    private static void writeManifest(File manifest, String rawFile, String expName) throws IOException {
        manifest.getParentFile().mkdirs();
        try (PrintWriter pw = new PrintWriter(manifest)) {
            // file_path \t experiment_name \t fraction \t DDA
            pw.println(rawFile + "\t" + expName.replaceAll("-", "_") + "\t1\tDDA");
        }
    }

    private static String safeExpName(String[] expNames, int idx) {
        if (expNames == null || idx < 0 || idx >= expNames.length) return "sample_" + (idx + 1);
        String s = expNames[idx];
        if (s == null || s.isBlank()) return "sample_" + (idx + 1);
        return s.trim();
    }

    private static String baseNameNoExt(String path) {
        if (path == null) return "sample";
        String name = new File(path).getName();
        int dot = name.lastIndexOf('.');
        return dot > 0 ? name.substring(0, dot) : name;
    }

    private static HashSet<String> selectGenomes(
            MagDbItem mag,
            DDAMagFragPipeParameters params,
            File psmTsv,
            File evidenceFile,
            File rankedFile,
            File selectedFile,
            ProgressCallback cb
    ) throws Exception {
        boolean useMarker4 = MagDbItem.MARKER4.equalsIgnoreCase(params.getGenomeSelectorMode());
        var psmList = new java.util.ArrayList<>(new FragpipePSMTsvReader(psmTsv).read());
        HashSet<String> genomeSet = new HashSet<>();

        if (!psmList.isEmpty()) {
            var selector = new MarkerGenomeSelector();
            selector.setKneeSelectionMode(params.getGenomeSelectionStrategy());
            int maximumSelected = 300;
            boolean enableAutoRep = useMarker4;

            if (useMarker4) {
                File current = mag.getCurrentFile();
                File markerFasta = new File(current, "Marker4.fasta");
                if (!markerFasta.exists()) {
                    markerFasta = mag.getMarker4FastaForLibrary();
                }
                File genomesTaxaTsv = new File(current, "genomes-taxa.tsv");
                if (markerFasta == null || !markerFasta.exists() || !genomesTaxaTsv.exists()) {
                    throw new IOException("Marker4 mode requires Marker4.fasta and genomes-taxa.tsv");
                }
                selector.selectMarkerGenomesWithKnee(
                        psmList.toArray(new FragPipePSM[0]),
                        markerFasta,
                        genomesTaxaTsv,
                        evidenceFile.getAbsolutePath(),
                        rankedFile.getAbsolutePath(),
                        selectedFile.getAbsolutePath(),
                        maximumSelected,
                        enableAutoRep
                );
            } else {
                selector.selectMarkerGenomesWithKnee(
                        psmList.toArray(new FragPipePSM[0]),
                        evidenceFile.getAbsolutePath(),
                        rankedFile.getAbsolutePath(),
                        selectedFile.getAbsolutePath(),
                        maximumSelected,
                        enableAutoRep
                );
            }
        }

        if (selectedFile.exists()) {
            try (BufferedReader reader = new BufferedReader(new FileReader(selectedFile))) {
                reader.readLine(); // header
                String line;
                while ((line = reader.readLine()) != null) {
                    line = line.trim();
                    if (line.isEmpty() || line.startsWith("#")) break;
                    String[] parts = line.split("\t");
                    if (parts.length > 0 && !parts[0].trim().isEmpty()) {
                        genomeSet.add(parts[0].trim());
                    }
                }
            }
        }

        // fallback: infer genomes from PSM protein names
        if (genomeSet.isEmpty() && !psmList.isEmpty()) {
            for (FragPipePSM psm : psmList) {
                String[] proteins = psm.getProteins();
                if (proteins == null) continue;
                for (String pro : proteins) {
                    if (pro == null) continue;
                    pro = pro.trim();
                    if (!pro.startsWith(mag.getIdentifier())) continue;
                    int idx = pro.indexOf("_");
                    if (idx > 0) {
                        genomeSet.add(pro.substring(0, idx));
                    }
                }
            }
            if (!genomeSet.isEmpty()) {
                cb.onMessage("WARNING: marker selector produced no genomes; falling back to " + genomeSet.size() + " genomes inferred from proteins");
            }
        }

        if (genomeSet.isEmpty()) {
            cb.onMessage("WARNING: no genomes selected from " + psmTsv.getAbsolutePath());
        }
        return genomeSet;
    }

    /**
     * Append all sequences from {@code fastaFile} and create simple reverse decoys as >rev_ entries.
     * (Matches the existing GUI behavior used for sample-specific DB construction.)
     */
    private static void appendFastaWithDecoys(File fastaFile, PrintWriter writer) throws IOException {
        try (BufferedReader reader = new BufferedReader(new FileReader(fastaFile))) {
            String line;
            StringBuilder seq = null;
            String revHeader = null;

            while ((line = reader.readLine()) != null) {
                if (line.startsWith(">")) {
                    if (revHeader != null && seq != null) {
                        writer.println(revHeader);
                        writer.println(seq.reverse());
                    }
                    seq = new StringBuilder();
                    revHeader = ">rev_" + line.substring(1);
                    writer.println(line);
                } else {
                    if (seq != null) seq.append(line.trim());
                    writer.println(line);
                }
            }
            if (revHeader != null && seq != null) {
                writer.println(revHeader);
                writer.println(seq.reverse());
            }
        }
    }

    private static String getCleavageResiduesOrDefault(FragpipeWorkflowFile wf) {
        try {
            String cut = wf.get("msfragger.search_enzyme_cut_1");
            if (cut != null && !cut.trim().isEmpty()) {
                return cut.trim();
            }
        } catch (Exception ignored) {
        }
        return "KR";
    }

    private static void appendHostFastaWithDecoys(File hostFasta, PrintWriter writer, String cleavage) throws IOException {
        try (BufferedReader reader = new BufferedReader(new FileReader(hostFasta))) {
            String line;
            String header = null;
            StringBuilder seqSb = new StringBuilder();
            while ((line = reader.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty()) continue;
                if (line.startsWith(">")) {
                    if (header != null) {
                        String seq = seqSb.toString();
                        writer.println(header);
                        writer.println(seq);

                        String shuffle = DecoyShuffler.generateShuffledDecoy(seq, cleavage, 42);
                        writer.println(">shuf_" + header.substring(1));
                        writer.println(shuffle);

                        writer.println(">rev_" + header.substring(1));
                        writer.println(new StringBuilder(seq).reverse().toString());

                        writer.println(">rev_shuf_" + header.substring(1));
                        writer.println(new StringBuilder(shuffle).reverse().toString());
                    }
                    header = line;
                    seqSb.setLength(0);
                } else {
                    seqSb.append(line);
                }
            }
            if (header != null) {
                String seq = seqSb.toString();
                writer.println(header);
                writer.println(seq);

                String shuffle = DecoyShuffler.generateShuffledDecoy(seq, cleavage, 42);
                writer.println(">shuf_" + header.substring(1));
                writer.println(shuffle);

                writer.println(">rev_" + header.substring(1));
                writer.println(new StringBuilder(seq).reverse().toString());

                writer.println(">rev_shuf_" + header.substring(1));
                writer.println(new StringBuilder(shuffle).reverse().toString());
            }
        }
    }

    private static void exportTaxonomyTsvs(DDAMagFragPipeParameters params,
                                          File combinedModifiedPeptideTsv,
                                          MagDbItem mag,
                                          File resultRoot,
                                          FragpipePostAnalysisTsvExporter.Result post,
                                          ProgressCallback cb) throws Exception {
        if (cb == null) cb = new ProgressCallback.NoOp();
        if (mag == null) throw new IllegalArgumentException("mag is null");
        if (post == null) throw new IllegalArgumentException("post-analysis result is null");

        File taxDir = new File(resultRoot, "taxonomy_analysis");
        taxDir.mkdirs();

        // Ensure combined ranked.tsv exists so we can label core/mid/tail.
        File combinedDir = new File(resultRoot, "combined");
        File rankedFile = new File(combinedDir, "ranked.tsv");
        if (!rankedFile.exists()) {
            generateCombinedRankedFiles(params, combinedModifiedPeptideTsv, mag, combinedDir, cb);
        }
        RankedTsvLabelParser.RankedGenomeSelectionResult labels =
                RankedTsvLabelParser.buildExportSelectionFromRankedTsv(rankedFile, post.genomeIntensityMap.keySet());

        File taxaTsv = mag.getOriginalTaxaFile();
        File gtdbNcbi = new File("Resources", "gtdb_ncbi");
        if (!gtdbNcbi.exists()) {
            gtdbNcbi = null;
        }
        MagTaxonomyTsvExporter.exportGenomeAndTaxaReports(
                taxDir,
                taxaTsv,
                gtdbNcbi,
                post.genomeIntensityMap,
                post.genomeScoreMap,
                post.genomeDiaStats,
                post.psmCountMap,
                labels.genomeSelectionLabel,
                post.totalPepCountMap,
                post.razorPepCountMap,
                post.uniquePepCountMap,
                post.totalPepWithILMap,
                post.razorPepWithILMap,
                post.uniquePepWithILMap,
                post.intensityTitles,
                cb
        );
    }

    private static void generateCombinedRankedFiles(DDAMagFragPipeParameters params,
                                                   File combinedModifiedPeptideTsv,
                                                   MagDbItem mag,
                                                   File combinedDir,
                                                   ProgressCallback cb) throws Exception {
        if (cb == null) cb = new ProgressCallback.NoOp();
        combinedDir.mkdirs();
        cb.onMessage("Generating combined evidence/ranked/selected TSVs for taxonomy labels ...");

        FragpipePepReader reader = new FragpipePepReader(combinedModifiedPeptideTsv);
        FragpipePeptide[] peps = reader.getMetaPeptides();
        FragPipePSM[] psms = reader.getPSMs();
        if (psms == null || psms.length == 0) {
            // build pseudo PSMs from peptides (same as GUI fallback)
            ArrayList<FragPipePSM> list = new ArrayList<>();
            if (peps != null) {
                for (FragpipePeptide p : peps) {
                    String[] prots = p.getProteins();
                    if (prots == null || prots.length == 0) continue;
                    list.add(new FragPipePSM("", 0, "", p.getSequence(),
                            p.getModSeq() != null ? p.getModSeq() : p.getSequence(),
                            "", 0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, 0.0, prots));
                }
            }
            psms = list.toArray(new FragPipePSM[0]);
        }

        File evidenceTsv = new File(combinedDir, "evidence.tsv");
        File rankedTsv = new File(combinedDir, "ranked.tsv");
        File selectedTsv = new File(combinedDir, "selected.tsv");

        boolean useMarker4 = MagDbItem.MARKER4.equalsIgnoreCase(params.getGenomeSelectorMode());
        MarkerGenomeSelector selector = new MarkerGenomeSelector();
        selector.setKneeSelectionMode(MarkerGenomeSelector.KNEE_MODE_COMPREHENSIVE);
        int maximumSelected = MarkerGenomeSelector.UNLIMITED_MAXIMUM_SELECTED;
        boolean enableAutoRep = useMarker4;

        if (useMarker4) {
            File current = mag.getCurrentFile();
            File markerFasta = new File(current, "Marker4.fasta");
            if (!markerFasta.exists()) markerFasta = mag.getMarker4FastaForLibrary();
            File genomesTaxaTsv = new File(current, "genomes-taxa.tsv");
            if (markerFasta != null && markerFasta.exists() && genomesTaxaTsv.exists()) {
                selector.selectMarkerGenomesWithKnee(psms, markerFasta, genomesTaxaTsv,
                        evidenceTsv.getAbsolutePath(), rankedTsv.getAbsolutePath(), selectedTsv.getAbsolutePath(),
                        maximumSelected, enableAutoRep);
                return;
            }
        }
        selector.selectMarkerGenomesWithKnee(psms,
                evidenceTsv.getAbsolutePath(), rankedTsv.getAbsolutePath(), selectedTsv.getAbsolutePath(),
                maximumSelected, enableAutoRep);
    }
}

