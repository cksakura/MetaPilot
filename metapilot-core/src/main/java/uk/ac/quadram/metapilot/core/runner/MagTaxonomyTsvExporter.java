package uk.ac.quadram.metapilot.core.runner;

import uk.ac.quadram.metapilot.core.mag.DdaMagGenomeScoringService;
import uk.ac.quadram.metapilot.core.mag.DdaMagProteinInferenceService;
import uk.ac.quadram.metapilot.core.mag.GenomeEvidenceTest;
import uk.ac.quadram.metapilot.core.taxonomy.TaxonomyRanks;
import uk.ac.quadram.metapilot.core.tools.FormatTool;
import uk.ac.quadram.metapilot.core.tools.ResourceLoader;
import uk.ac.quadram.metapilot.core.tools.TypedParquetTableWriter;
import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

/**
 * Headless taxonomy report exporter.
 *
 * Writes:
 * - taxonomy_analysis/Genome.html
 * - taxonomy_analysis/Genome.tsv
 * - taxonomy_analysis/Genome_unique_peptides.tsv (same columns; intensity from unique peptides)
 * - taxonomy_analysis/Taxa.tsv (DIA-style scoring)
 * - taxonomy_analysis/Taxa_report.tsv
 * - taxonomy_analysis/genome_rank_scores.tsv (DIA-style scoring), or genome_scores.tsv and
 *   genome_scores_runinfo.txt when the genomes were tested by GenomeEvidenceTest (open search)
 *
 * This intentionally avoids the GUI-only dependency on a packaged "gtdb_ncbi" mapping resource.
 * Instead, it uses the MAG DB's original taxa TSV directly (from MAGDB_1_0.json version metadata).
 */
public final class MagTaxonomyTsvExporter {

    static final String GENOME_SCORES_TSV = "genome_scores.tsv";
    static final String GENOME_SCORES_RUNINFO = "genome_scores_runinfo.txt";
    static final String GENOME_UNIQUE_PEPTIDES_TSV = "Genome_unique_peptides.tsv";

    private MagTaxonomyTsvExporter() {
    }

    public static void exportGenomeAndTaxaReports(File taxonomyAnalysisDir,
                                                  File originalTaxaTsv,
                                                  File gtdbNcbiMappingTsv,
                                                  Map<String, double[]> genomeIntensityMap,
                                                  Map<String, Double> genomeScoreMap,
                                                  Map<String, double[]> genomeDiaStats,
                                                  Map<String, Integer> genomePsmCountMap,
                                                  Map<String, String> genomeSelectionLabel,
                                                  Map<String, ? extends Set<String>> totalPepCountMap,
                                                  Map<String, ? extends Set<String>> razorPepCountMap,
                                                  Map<String, ? extends Set<String>> uniquePepCountMap,
                                                  Map<String, ? extends Set<String>> totalPepWithILMap,
                                                  Map<String, ? extends Set<String>> razorPepWithILMap,
                                                  Map<String, ? extends Set<String>> uniquePepWithILMap,
                                                  String[] intensityTitles,
                                                  ProgressCallback cb) throws IOException {
        exportGenomeAndTaxaReports(taxonomyAnalysisDir, originalTaxaTsv, gtdbNcbiMappingTsv, genomeIntensityMap,
                genomeScoreMap, genomeDiaStats, null, genomePsmCountMap, genomeSelectionLabel, totalPepCountMap,
                razorPepCountMap, uniquePepCountMap, totalPepWithILMap, razorPepWithILMap, uniquePepWithILMap,
                intensityTitles, null, cb);
    }

    public static void exportGenomeAndTaxaReports(File taxonomyAnalysisDir,
                                                  File originalTaxaTsv,
                                                  File gtdbNcbiMappingTsv,
                                                  Map<String, double[]> genomeIntensityMap,
                                                  Map<String, Double> genomeScoreMap,
                                                  Map<String, double[]> genomeDiaStats,
                                                  List<String> orderedGenomes,
                                                  Map<String, Integer> genomePsmCountMap,
                                                  Map<String, String> genomeSelectionLabel,
                                                  Map<String, ? extends Set<String>> totalPepCountMap,
                                                  Map<String, ? extends Set<String>> razorPepCountMap,
                                                  Map<String, ? extends Set<String>> uniquePepCountMap,
                                                  Map<String, ? extends Set<String>> totalPepWithILMap,
                                                  Map<String, ? extends Set<String>> razorPepWithILMap,
                                                  Map<String, ? extends Set<String>> uniquePepWithILMap,
                                                  String[] intensityTitles,
                                                  DdaMagGenomeScoringService.DdaGenomeScoreRuninfo runinfo,
                                                  ProgressCallback cb) throws IOException {
        exportGenomeAndTaxaReports(taxonomyAnalysisDir, originalTaxaTsv, gtdbNcbiMappingTsv, genomeIntensityMap,
                genomeScoreMap, genomeDiaStats, orderedGenomes, genomePsmCountMap, genomeSelectionLabel,
                totalPepCountMap, razorPepCountMap, uniquePepCountMap, totalPepWithILMap, razorPepWithILMap,
                uniquePepWithILMap, intensityTitles, runinfo, null, null, cb);
    }

    /**
     * @param genomeEvidence when not null, {@code genomeDiaStats} holds its
     *                       {kept peptides, p, q}: Genome.tsv names the columns
     *                       accordingly and genome_scores.tsv with its run-info
     *                       replaces genome_rank_scores.tsv ({@code runinfo} is
     *                       then ignored)
     */
    public static void exportGenomeAndTaxaReports(File taxonomyAnalysisDir,
                                                  File originalTaxaTsv,
                                                  File gtdbNcbiMappingTsv,
                                                  Map<String, double[]> genomeIntensityMap,
                                                  Map<String, Double> genomeScoreMap,
                                                  Map<String, double[]> genomeDiaStats,
                                                  List<String> orderedGenomes,
                                                  Map<String, Integer> genomePsmCountMap,
                                                  Map<String, String> genomeSelectionLabel,
                                                  Map<String, ? extends Set<String>> totalPepCountMap,
                                                  Map<String, ? extends Set<String>> razorPepCountMap,
                                                  Map<String, ? extends Set<String>> uniquePepCountMap,
                                                  Map<String, ? extends Set<String>> totalPepWithILMap,
                                                  Map<String, ? extends Set<String>> razorPepWithILMap,
                                                  Map<String, ? extends Set<String>> uniquePepWithILMap,
                                                  String[] intensityTitles,
                                                  DdaMagGenomeScoringService.DdaGenomeScoreRuninfo runinfo,
                                                  GenomeEvidenceTest.Result genomeEvidence,
                                                  Map<String, double[]> genomeUniqueIntensityMap,
                                                  ProgressCallback cb) throws IOException {
        if (cb == null) cb = new ProgressCallback.NoOp();
        if (taxonomyAnalysisDir == null) throw new IOException("taxonomyAnalysisDir is null");
        taxonomyAnalysisDir.mkdirs();

        if (originalTaxaTsv == null || !originalTaxaTsv.exists()) {
            throw new IOException("MAG taxa TSV not found: " + (originalTaxaTsv == null ? "null" : originalTaxaTsv.getAbsolutePath()));
        }

        if (genomeIntensityMap == null) {
            throw new IOException("genomeIntensityMap is null");
        }

        if (intensityTitles == null) intensityTitles = new String[0];
        if (genomeScoreMap == null) genomeScoreMap = new HashMap<>();
        if (genomePsmCountMap == null) genomePsmCountMap = new HashMap<>();

        cb.onMessage("Loading MAG taxa TSV: " + originalTaxaTsv.getAbsolutePath());
        HashMap<String, Lineage> lineageMap = loadLineageByGenome(originalTaxaTsv);

        HashMap<String, String> gtdbToNcbi = null;
        if (gtdbNcbiMappingTsv != null && gtdbNcbiMappingTsv.exists()) {
            cb.onMessage("Loading GTDB->NCBI mapping: " + gtdbNcbiMappingTsv.getAbsolutePath());
            gtdbToNcbi = loadGtdbNcbiMap(gtdbNcbiMappingTsv);
        } else {
            gtdbToNcbi = loadBundledGtdbNcbiMap();
        }

        File genomeHtml = new File(taxonomyAnalysisDir, "Genome.html");
        cb.onMessage("Writing taxonomy: " + genomeHtml.getName());
        writeGenomeHtml(genomeHtml, genomeIntensityMap.keySet(), lineageMap);

        File genomeTsv = new File(taxonomyAnalysisDir, "Genome.tsv");
        cb.onMessage("Writing taxonomy: " + genomeTsv.getName());
        writeGenomeTsv(genomeTsv, genomeIntensityMap, genomeScoreMap, genomeDiaStats, orderedGenomes,
                genomePsmCountMap, lineageMap, gtdbToNcbi,
                genomeSelectionLabel,
                totalPepCountMap, razorPepCountMap, uniquePepCountMap,
                totalPepWithILMap, razorPepWithILMap, uniquePepWithILMap,
                intensityTitles, genomeEvidence != null);

        if (genomeUniqueIntensityMap != null) {
            File uniqueTsv = new File(taxonomyAnalysisDir, GENOME_UNIQUE_PEPTIDES_TSV);
            cb.onMessage("Writing taxonomy: " + uniqueTsv.getName());
            writeGenomeTsv(uniqueTsv, alignIntensityToGenomes(genomeUniqueIntensityMap, genomeIntensityMap),
                    genomeScoreMap, genomeDiaStats, orderedGenomes,
                    genomePsmCountMap, lineageMap, gtdbToNcbi,
                    genomeSelectionLabel,
                    totalPepCountMap, razorPepCountMap, uniquePepCountMap,
                    totalPepWithILMap, razorPepWithILMap, uniquePepWithILMap,
                    intensityTitles, genomeEvidence != null);
        }

        if (genomeEvidence != null) {
            File scoresFile = new File(taxonomyAnalysisDir, GENOME_SCORES_TSV);
            cb.onMessage("Writing taxonomy: " + scoresFile.getName());
            writeGenomeScores(scoresFile, genomeEvidence, genusByGenome(lineageMap));
            File runinfoFile = new File(taxonomyAnalysisDir, GENOME_SCORES_RUNINFO);
            cb.onMessage("Writing taxonomy: " + runinfoFile.getName());
            writeGenomeScoresRuninfo(runinfoFile, scoresFile, genomeEvidence);
        }

        boolean useDiaStyle = genomeDiaStats != null && !genomeDiaStats.isEmpty();
        if (useDiaStyle) {
            File taxaTsv = new File(taxonomyAnalysisDir, "Taxa.tsv");
            File taxaReportTsv = new File(taxonomyAnalysisDir, "Taxa_report.tsv");
            cb.onMessage("Writing taxonomy: " + taxaTsv.getName());
            writeDiaStyleTaxaTables(taxaTsv, taxaReportTsv, genomeIntensityMap, genomeScoreMap, genomeDiaStats,
                    lineageMap, intensityTitles);
            if (genomeEvidence == null) {
                File rankScoresFile = new File(taxonomyAnalysisDir, "genome_rank_scores.tsv");
                cb.onMessage("Writing taxonomy: " + rankScoresFile.getName());
                writeGenomeRankScores(rankScoresFile, genomeDiaStats, totalPepCountMap, razorPepCountMap);
                if (runinfo != null) {
                    File runinfoFile = new File(taxonomyAnalysisDir, "genome_rank_scores_runinfo.txt");
                    cb.onMessage("Writing taxonomy: " + runinfoFile.getName());
                    writeGenomeRankScoresRuninfo(runinfoFile, rankScoresFile, runinfo);
                }
            }
        } else {
            File taxaReportTsv = new File(taxonomyAnalysisDir, "Taxa_report.tsv");
            cb.onMessage("Writing taxonomy: " + taxaReportTsv.getName());
            writeTaxaReport(taxaReportTsv, genomeIntensityMap, lineageMap, intensityTitles);
        }
    }

    /**
     * Opt-in MetaAtlas export: writes taxon_abundance.parquet (see
     * {@link #writeTaxonAbundanceParquet}) alongside the standard TSV/HTML reports written
     * by {@link #exportGenomeAndTaxaReports}, but is never called by that method - callers
     * (CLI/headless runners only; not GUI tasks) must invoke this explicitly and separately,
     * gated by their own opt-in flag (e.g. {@code --export-metaatlas}).
     */
    public static void exportMetaAtlasTaxonAbundance(File taxonomyAnalysisDir,
                                                     File originalTaxaTsv,
                                                     Map<String, double[]> genomeIntensityMap,
                                                     Map<String, double[]> genomeDiaStats,
                                                     Map<String, Double> genomeScoreMap,
                                                     String[] intensityTitles,
                                                     ProgressCallback cb) throws IOException {
        if (cb == null) cb = new ProgressCallback.NoOp();
        if (taxonomyAnalysisDir == null) throw new IOException("taxonomyAnalysisDir is null");
        taxonomyAnalysisDir.mkdirs();

        if (originalTaxaTsv == null || !originalTaxaTsv.exists()) {
            throw new IOException(
                    "MAG taxa TSV not found: " + (originalTaxaTsv == null ? "null" : originalTaxaTsv.getAbsolutePath()));
        }
        if (genomeIntensityMap == null) {
            throw new IOException("genomeIntensityMap is null");
        }
        if (intensityTitles == null) intensityTitles = new String[0];

        HashMap<String, Lineage> lineageMap = loadLineageByGenome(originalTaxaTsv);
        File out = new File(taxonomyAnalysisDir, "taxon_abundance.parquet");
        cb.onMessage("Writing MetaAtlas export: " + out.getName());
        writeTaxonAbundanceParquet(out, genomeIntensityMap, genomeDiaStats, genomeScoreMap, lineageMap,
                intensityTitles);
    }

    private static HashMap<String, Lineage> loadLineageByGenome(File taxaTsv) throws IOException {
        HashMap<String, Lineage> out = new HashMap<>();
        try (BufferedReader br = new BufferedReader(new FileReader(taxaTsv))) {
            String header = br.readLine();
            if (header == null) return out;
            String[] title = header.split("\t");

            int genomeId = 0; // usually first column is genome id
            int lineageId = -1;
            for (int i = 0; i < title.length; i++) {
                String t = title[i].trim();
                if (t.equalsIgnoreCase("Lineage") || t.equalsIgnoreCase("Species_designation")) {
                    lineageId = i;
                }
            }
            if (lineageId < 0) {
                // fallback: try second column if present
                lineageId = Math.min(1, title.length - 1);
            }

            String line;
            while ((line = br.readLine()) != null) {
                if (line.isBlank()) continue;
                String[] cs = line.split("\t");
                if (cs.length <= Math.max(genomeId, lineageId)) continue;
                String genome = cs[genomeId].trim();
                String lineage = cs[lineageId].trim();
                if (genome.isEmpty()) continue;
                out.put(genome, Lineage.parse(lineage));
            }
        }
        return out;
    }

    private static HashMap<String, String> loadGtdbNcbiMap(File mappingTsv) throws IOException {
        HashMap<String, String> out = new HashMap<>();
        try (BufferedReader br = new BufferedReader(new FileReader(mappingTsv))) {
            String header = br.readLine(); // gtdb_taxonomy\tncbi_taxonomy
            if (header == null) return out;
            String line;
            while ((line = br.readLine()) != null) {
                if (line.isBlank()) continue;
                String[] cs = line.split("\t");
                if (cs.length < 2) continue;
                out.put(cs[0].trim(), cs[1].trim());
            }
        }
        return out;
    }

    private static HashMap<String, String> loadBundledGtdbNcbiMap() throws IOException {
        HashMap<String, String> out = new HashMap<>();
        try (InputStream input = ResourceLoader.load("gtdb_ncbi")) {
            if (input == null) return out;
            try (BufferedReader br = new BufferedReader(new InputStreamReader(input, "UTF-8"))) {
                br.readLine();
                String line;
                while ((line = br.readLine()) != null) {
                    if (line.isBlank()) continue;
                    String[] cs = line.split("\t");
                    if (cs.length < 2) continue;
                    out.put(cs[0].trim(), cs[1].trim());
                }
            }
        } catch (IOException e) {
            return out;
        }
        return out;
    }

    private static void writeGenomeTsv(File out,
                                      Map<String, double[]> genomeIntensityMap,
                                      Map<String, Double> genomeScoreMap,
                                      Map<String, double[]> genomeDiaStats,
                                      List<String> orderedGenomes,
                                      Map<String, Integer> genomePsmCountMap,
                                      Map<String, Lineage> lineageMap,
                                      Map<String, String> gtdbToNcbi,
                                      Map<String, String> genomeSelectionLabel,
                                      Map<String, ? extends Set<String>> totalPepCountMap,
                                      Map<String, ? extends Set<String>> razorPepCountMap,
                                      Map<String, ? extends Set<String>> uniquePepCountMap,
                                      Map<String, ? extends Set<String>> totalPepWithILMap,
                                      Map<String, ? extends Set<String>> razorPepWithILMap,
                                      Map<String, ? extends Set<String>> uniquePepWithILMap,
                                      String[] intensityTitles,
                                      boolean evidenceStats) throws IOException {
        TaxonomyRanks[] mainRanks = TaxonomyRanks.getMainRanks();
        String[] genomes = orderGenomes(genomeIntensityMap, orderedGenomes);

        try (PrintWriter pw = new PrintWriter(out)) {
            StringBuilder h = new StringBuilder();
            h.append("Genome\t");
            boolean useDiaStyleHeader = (genomeDiaStats != null && !genomeDiaStats.isEmpty());
            if (useDiaStyleHeader && evidenceStats) {
                h.append("kept_unique_peptides\t");
                h.append("genome_p_value\t");
                h.append("genome_q_value\t");
            } else if (useDiaStyleHeader) {
                h.append("genome_score_S\t");
                h.append("genome_p_emp\t");
                h.append("genome_q_value\t");
            } else {
                h.append("Minus_log(p-value)\t");
            }
            h.append("Selection_label\t");
            h.append("PSM count\t");
            h.append("Peptide count\t");
            h.append("Razor peptide count\t");
            h.append("Unique peptide count\t");
            h.append("Peptide I/L ratio\t");
            h.append("Razor peptide I/L ratio\t");
            h.append("Unique peptide I/L ratio\t");
            h.append("GTDB_Name\t");
            h.append("GTDB_Rank\t");
            for (TaxonomyRanks r : mainRanks) {
                h.append("GTDB_").append(r.getName()).append("\t");
            }
            h.append("NCBI_Name\t");
            h.append("NCBI_Rank");
            for (TaxonomyRanks r : mainRanks) {
                h.append("\t").append("NCBI_").append(r.getName());
            }
            for (String exp : intensityTitles) {
                if (exp == null || exp.isBlank()) continue;
                h.append("\t").append(exp.startsWith("Intensity ") ? exp : ("Intensity " + exp));
            }
            pw.println(h);

            for (String g : genomes) {
                double[] inten = genomeIntensityMap.get(g);
                double score = genomeScoreMap.getOrDefault(g, 1.0);
                int psm = genomePsmCountMap.getOrDefault(g, 0);

                Lineage lin = lineageMap.get(g);
                if (lin == null) lin = Lineage.empty();

                Lineage ncbi = Lineage.empty();
                if (gtdbToNcbi != null && lin.rawLineage != null && !lin.rawLineage.isBlank()) {
                    String ncbiLineage = gtdbToNcbi.get(lin.rawLineage);
                    if (ncbiLineage != null && !ncbiLineage.isBlank()) {
                        ncbi = Lineage.parse(ncbiLineage);
                    }
                }

                StringBuilder sb = new StringBuilder();
                sb.append(g).append("\t");
                Map<String, double[]> dia = genomeDiaStats;
                double[] st = (dia != null) ? dia.get(g) : null;
                if (st != null && st.length >= 3 && evidenceStats) {
                    sb.append((long) st[0]).append("\t");
                    sb.append(FormatTool.getDFE4().format(st[1])).append("\t");
                    sb.append(FormatTool.getDFE4().format(st[2])).append("\t");
                } else if (st != null && st.length >= 3) {
                    sb.append(FormatTool.getDFE4().format(st[0])).append("\t");
                    sb.append(FormatTool.getDF4().format(st[1])).append("\t");
                    sb.append(FormatTool.getDF4().format(st[2])).append("\t");
                } else {
                    double minusLog = score > 0 ? -Math.log10(score) : 0.0;
                    sb.append(FormatTool.getDF2().format(minusLog)).append("\t");
                }
                sb.append(genomeSelectionLabel != null ? genomeSelectionLabel.getOrDefault(g, "tail") : "tail").append("\t");
                sb.append(psm).append("\t");
                int totalPep = totalPepCountMap != null && totalPepCountMap.containsKey(g) ? totalPepCountMap.get(g).size() : 0;
                int razorPep = razorPepCountMap != null && razorPepCountMap.containsKey(g) ? razorPepCountMap.get(g).size() : 0;
                int uniqPep = uniquePepCountMap != null && uniquePepCountMap.containsKey(g) ? uniquePepCountMap.get(g).size() : 0;
                sb.append(totalPep).append("\t");
                sb.append(razorPep).append("\t");
                sb.append(uniqPep).append("\t");
                double totalIL = (totalPep > 0 && totalPepWithILMap != null && totalPepWithILMap.containsKey(g))
                        ? ((double) totalPepWithILMap.get(g).size() / (double) totalPep) : 0.0;
                double razorIL = (razorPep > 0 && razorPepWithILMap != null && razorPepWithILMap.containsKey(g))
                        ? ((double) razorPepWithILMap.get(g).size() / (double) razorPep) : 0.0;
                double uniqIL = (uniqPep > 0 && uniquePepWithILMap != null && uniquePepWithILMap.containsKey(g))
                        ? ((double) uniquePepWithILMap.get(g).size() / (double) uniqPep) : 0.0;
                sb.append(FormatTool.getDF4().format(totalIL)).append("\t");
                sb.append(FormatTool.getDF4().format(razorIL)).append("\t");
                sb.append(FormatTool.getDF4().format(uniqIL)).append("\t");

                sb.append(lin.name).append("\t");
                sb.append(lin.rank).append("\t");
                for (TaxonomyRanks r : mainRanks) {
                    sb.append(lin.byRank.getOrDefault(r.getName(), "")).append("\t");
                }

                sb.append(ncbi.name).append("\t");
                sb.append(ncbi.rank);
                for (TaxonomyRanks r : mainRanks) {
                    sb.append("\t").append(ncbi.byRank.getOrDefault(r.getName(), ""));
                }

                if (inten != null) {
                    for (double v : inten) {
                        sb.append("\t").append(FormatTool.getDFE4().format(v));
                    }
                }
                pw.println(sb);
            }
        }
    }

    /**
     * Writes taxon_abundance.parquet: one row per taxon (at each main rank) per sample,
     * with real types (abundance as DOUBLE) instead of a TSV-mirrored string table -
     * shaped to match the MetaAtlas Parquet Export Upgrade Guide's taxon_abundance view
     * (section 9.8) so it can be queried/joined by sample in DuckDB.
     *
     * Known deviations from the guide, since this pivots existing per-genome maps rather
     * than adding a new taxonomy-ID resolution layer:
     * - There is no numeric/accession taxon_id anywhere in this codebase's taxonomy path
     *   (only free-text names derived from GTDB lineage strings), so taxon_name is used as
     *   part of the natural key instead of a stable taxon_id. The guide explicitly warns
     *   against using display names as primary keys (section 4.3) - flagging this rather
     *   than pretending otherwise.
     * - taxonomy_source is always "GTDB" here because the aggregation upstream (see
     *   writeDiaStyleTaxaTables/writeTaxaReport) only aggregates against the GTDB-derived
     *   Lineage, not the mapped NCBI lineage.
     * - sample_id is the raw experiment/run name from the intensity titles, not a resolved
     *   SDRF-backed stable sample_accession (guide section 12) - no sample metadata layer
     *   exists yet to do that resolution.
     * - supporting_peptide_count / unique_peptide_count / minimum_evidence_passed (guide
     *   fields) are omitted: peptide counts are only aggregated per-genome upstream, not
     *   per-taxon, so there's nothing correct to put there yet.
     */
    private static void writeTaxonAbundanceParquet(File out,
                                                   Map<String, double[]> genomeIntensityMap,
                                                   Map<String, double[]> genomeDiaStats,
                                                   Map<String, Double> genomeScoreMap,
                                                   Map<String, Lineage> lineageMap,
                                                   String[] intensityTitles) throws IOException {
        TaxonomyRanks[] mainRanks = TaxonomyRanks.getMainRanks();
        String[] sampleIds = TypedParquetTableWriter.sampleIdsFromIntensityTitles(intensityTitles);
        boolean useDia = genomeDiaStats != null && !genomeDiaStats.isEmpty();

        // key = rank + "\t" + taxon name
        LinkedHashMap<String, double[]> taxonIntensity = new LinkedHashMap<>();
        LinkedHashMap<String, Double> taxonScore = new LinkedHashMap<>();
        LinkedHashMap<String, String> taxonRankOf = new LinkedHashMap<>();
        LinkedHashMap<String, String> taxonNameOf = new LinkedHashMap<>();

        for (String genome : genomeIntensityMap.keySet()) {
            double[] gi = genomeIntensityMap.get(genome);
            if (gi == null) continue;
            Lineage lin = lineageMap.get(genome);
            if (lin == null) continue;

            for (TaxonomyRanks r : mainRanks) {
                String name = lin.byRank.get(r.getName());
                if (name == null || name.isBlank()) continue;
                String key = r.getName() + "\t" + name;

                double[] acc = taxonIntensity.computeIfAbsent(key, k -> new double[gi.length]);
                for (int i = 0; i < Math.min(acc.length, gi.length); i++) {
                    acc[i] += gi[i];
                }
                taxonRankOf.put(key, r.getName());
                taxonNameOf.put(key, name);

                double contribution;
                if (useDia) {
                    double[] st = genomeDiaStats.get(genome);
                    contribution = (st != null && st.length > 0) ? st[0] : 0.0;
                } else {
                    double score = (genomeScoreMap != null) ? genomeScoreMap.getOrDefault(genome, 1.0) : 1.0;
                    contribution = score > 0 ? -Math.log10(score) : 0.0;
                }
                taxonScore.merge(key, contribution, Double::sum);
            }
        }

        List<TypedParquetTableWriter.Column> columns = Arrays.asList(
                new TypedParquetTableWriter.Column("taxon_name", TypedParquetTableWriter.ColumnType.STRING),
                new TypedParquetTableWriter.Column("taxon_rank", TypedParquetTableWriter.ColumnType.STRING),
                new TypedParquetTableWriter.Column("taxonomy_source", TypedParquetTableWriter.ColumnType.STRING),
                new TypedParquetTableWriter.Column("sample_id", TypedParquetTableWriter.ColumnType.STRING),
                new TypedParquetTableWriter.Column("abundance", TypedParquetTableWriter.ColumnType.DOUBLE),
                new TypedParquetTableWriter.Column("abundance_type", TypedParquetTableWriter.ColumnType.STRING),
                new TypedParquetTableWriter.Column("evidence_score", TypedParquetTableWriter.ColumnType.DOUBLE));

        List<Object[]> rows = new ArrayList<>();
        for (String key : taxonIntensity.keySet()) {
            double[] inten = taxonIntensity.get(key);
            String rank = taxonRankOf.get(key);
            String name = taxonNameOf.get(key);
            Double score = taxonScore.get(key);
            for (int i = 0; i < sampleIds.length && i < inten.length; i++) {
                rows.add(new Object[] { name, rank, "GTDB", sampleIds[i], inten[i], "raw", score });
            }
        }

        TypedParquetTableWriter.write(out, columns, rows);
    }

    private static void writeDiaStyleTaxaTables(File taxaTsv,
                                              File taxaReportTsv,
                                              Map<String, double[]> genomeIntensityMap,
                                              Map<String, Double> genomeScoreMap,
                                              Map<String, double[]> genomeDiaStats,
                                              Map<String, Lineage> lineageMap,
                                              String[] intensityTitles) throws IOException {
        TaxonomyRanks[] mainRanks = TaxonomyRanks.getMainRanks();
        @SuppressWarnings("unchecked")
        HashMap<String, String[]>[] taxRankMaps = new HashMap[mainRanks.length];
        @SuppressWarnings("unchecked")
        HashMap<String, double[]>[] taxIntensityMaps = new HashMap[mainRanks.length];
        @SuppressWarnings("unchecked")
        HashMap<String, Double>[] taxScoreMaps = new HashMap[mainRanks.length];
        for (int i = 0; i < mainRanks.length; i++) {
            taxRankMaps[i] = new HashMap<>();
            taxIntensityMaps[i] = new HashMap<>();
            taxScoreMaps[i] = new HashMap<>();
        }

        for (String genome : genomeIntensityMap.keySet()) {
            double[] genomeIntensity = genomeIntensityMap.get(genome);
            if (genomeIntensity == null) continue;
            Lineage lin = lineageMap.get(genome);
            if (lin == null) lin = Lineage.empty();

            for (int j = 0; j < mainRanks.length; j++) {
                String taxon = lin.byRank.getOrDefault(mainRanks[j].getName(), "");
                if (taxon.isEmpty()) continue;

                if (!taxRankMaps[j].containsKey(taxon)) {
                    String[] lineage = new String[j + 1];
                    for (int k = 0; k <= j; k++) {
                        lineage[k] = lin.byRank.getOrDefault(mainRanks[k].getName(), "");
                    }
                    taxRankMaps[j].put(taxon, lineage);
                }

                double[] taxIntensity = taxIntensityMaps[j].get(taxon);
                if (taxIntensity == null) {
                    taxIntensity = new double[genomeIntensity.length];
                }
                for (int k = 0; k < taxIntensity.length; k++) {
                    taxIntensity[k] += genomeIntensity[k];
                }
                taxIntensityMaps[j].put(taxon, taxIntensity);

                double[] gStats = genomeDiaStats.get(genome);
                if (gStats != null && gStats.length > 0) {
                    double s = gStats[0];
                    taxScoreMaps[j].put(taxon, taxScoreMaps[j].getOrDefault(taxon, 0.0) + s);
                } else {
                    double score = genomeScoreMap.getOrDefault(genome, 1.0);
                    taxScoreMaps[j].put(taxon, taxScoreMaps[j].getOrDefault(taxon, 0.0) - Math.log10(score));
                }
            }
        }

        try (PrintWriter taxaWriter = new PrintWriter(taxaTsv);
             PrintWriter taxaReportWriter = new PrintWriter(taxaReportTsv)) {
            StringBuilder taxaTitlesb1 = new StringBuilder();
            taxaTitlesb1.append("Name\tScore\tRank");
            for (TaxonomyRanks r : mainRanks) taxaTitlesb1.append("\t").append(r.getName());
            for (String exp : intensityTitles) {
                if (exp == null || exp.isBlank()) continue;
                taxaTitlesb1.append("\t").append(exp.startsWith("Intensity ") ? exp : ("Intensity " + exp));
            }

            StringBuilder taxaTitlesb2 = new StringBuilder();
            taxaTitlesb2.append("Name\tRank");
            for (TaxonomyRanks r : mainRanks) taxaTitlesb2.append("\t").append(r.getName());
            for (String exp : intensityTitles) {
                if (exp == null || exp.isBlank()) continue;
                taxaTitlesb2.append("\t").append(exp.startsWith("Intensity ") ? exp : ("Intensity " + exp));
            }

            taxaWriter.println(taxaTitlesb1);
            taxaReportWriter.println(taxaTitlesb2);

            for (int i = 0; i < taxRankMaps.length; i++) {
                for (String taxon : taxRankMaps[i].keySet()) {
                    String[] taxaContent = taxRankMaps[i].get(taxon);
                    double[] intensity = taxIntensityMaps[i].get(taxon);
                    double taxonScore = taxScoreMaps[i].getOrDefault(taxon, 0.0);

                    StringBuilder taxaSb = new StringBuilder();
                    taxaSb.append(taxon).append("\t");
                    taxaSb.append(FormatTool.getDFE4().format(taxonScore)).append("\t");
                    taxaSb.append(mainRanks[i].getName());
                    for (int j = 0; j < mainRanks.length; j++) {
                        taxaSb.append("\t").append(j < taxaContent.length ? taxaContent[j] : "");
                    }
                    if (intensity != null) {
                        for (double v : intensity) taxaSb.append("\t").append(v);
                    }
                    taxaWriter.println(taxaSb);
                }
            }

            for (int i = 0; i < taxRankMaps.length; i++) {
                for (String taxon : taxRankMaps[i].keySet()) {
                    String[] taxaContent = taxRankMaps[i].get(taxon);
                    double[] intensity = taxIntensityMaps[i].get(taxon);

                    StringBuilder taxaSb = new StringBuilder();
                    taxaSb.append(taxon).append("\t").append(mainRanks[i].getName());
                    for (int j = 0; j < mainRanks.length; j++) {
                        taxaSb.append("\t").append(j < taxaContent.length ? taxaContent[j] : "");
                    }
                    if (intensity != null) {
                        for (double v : intensity) taxaSb.append("\t").append(v);
                    }
                    taxaReportWriter.println(taxaSb);
                }
            }
        }
    }

    private static void writeGenomeRankScores(File out,
                                              Map<String, double[]> genomeDiaStats,
                                              Map<String, ? extends Set<String>> totalPepCountMap,
                                              Map<String, ? extends Set<String>> razorPepCountMap) throws IOException {
        String[] statsGenomes = genomeDiaStats.keySet().toArray(new String[0]);
        Arrays.sort(statsGenomes, Comparator.comparingDouble(g -> -genomeDiaStats.get(g)[0]));
        try (PrintWriter pw = new PrintWriter(out)) {
            pw.println("Rank\tGenome\tgenome_score_S\tn_peptides\tgenome_p_emp\tgenome_q_value");
            for (int r = 0; r < statsGenomes.length; r++) {
                String g = statsGenomes[r];
                double[] st = genomeDiaStats.get(g);
                int nPep = 0;
                if (razorPepCountMap != null && razorPepCountMap.containsKey(g)) {
                    nPep = razorPepCountMap.get(g).size();
                } else if (totalPepCountMap != null && totalPepCountMap.containsKey(g)) {
                    nPep = totalPepCountMap.get(g).size();
                }
                pw.println((r + 1) + "\t" + g + "\t" + FormatTool.getDFE4().format(st[0]) + "\t" + nPep
                        + "\t" + FormatTool.getDF4().format(st[1]) + "\t" + FormatTool.getDF4().format(st[2]));
            }
        }
    }

    private static HashMap<String, String> genusByGenome(Map<String, Lineage> lineageMap) {
        HashMap<String, String> out = new HashMap<>();
        String genus = TaxonomyRanks.Genus.getName();
        for (Map.Entry<String, Lineage> e : lineageMap.entrySet()) {
            String name = e.getValue().byRank.get(genus);
            if (name != null && !name.isBlank()) {
                out.put(e.getKey(), name);
            }
        }
        return out;
    }

    /**
     * Accepted target genome -> the accepted genome with the most kept peptides
     * in its GTDB genus (ties by name), when that is another genome. A label, not
     * a rejection: the database genome closest to a present strain wins its
     * peptides, and this names the stronger same-genus neighbour.
     */
    static HashMap<String, String> sameGenusAs(List<GenomeEvidenceTest.GenomeRow> rows,
                                               Map<String, String> genusByGenome) {
        HashMap<String, GenomeEvidenceTest.GenomeRow> topByGenus = new HashMap<>();
        for (GenomeEvidenceTest.GenomeRow row : rows) {
            String genus = genusByGenome.get(row.genome);
            if (row.decoy || !row.accepted || genus == null) continue;
            GenomeEvidenceTest.GenomeRow top = topByGenus.get(genus);
            if (top == null || row.keptPeptides > top.keptPeptides
                    || (row.keptPeptides == top.keptPeptides && row.genome.compareTo(top.genome) < 0)) {
                topByGenus.put(genus, row);
            }
        }
        HashMap<String, String> out = new HashMap<>();
        for (GenomeEvidenceTest.GenomeRow row : rows) {
            String genus = genusByGenome.get(row.genome);
            if (row.decoy || !row.accepted || genus == null) continue;
            GenomeEvidenceTest.GenomeRow top = topByGenus.get(genus);
            if (!top.genome.equals(row.genome)) {
                out.put(row.genome, top.genome);
            }
        }
        return out;
    }

    static void writeGenomeScores(File out, GenomeEvidenceTest.Result evidence, Map<String, String> genusByGenome)
            throws IOException {
        HashMap<String, String> sameGenus = sameGenusAs(evidence.rows, genusByGenome);
        try (PrintWriter pw = new PrintWriter(out)) {
            pw.println("genome\tis_decoy\truns_searched\tunique_peptides\texcluded_near_neighbour\tkept_peptides"
                    + "\tdecoy_unique_peptides\texpected_false\tlog_p\tp_value\tq_value\taccepted\tacceptance_order"
                    + "\tresolved_rank\tsame_genus_as");
            for (GenomeEvidenceTest.GenomeRow row : evidence.rows) {
                pw.println(row.genome + "\t" + row.decoy + "\t" + row.runsSearched + "\t" + row.uniquePeptides
                        + "\t" + row.excludedNearNeighbour + "\t" + row.keptPeptides + "\t" + row.decoyUniquePeptides
                        + "\t" + FormatTool.getDFE4().format(row.expectedFalse)
                        + "\t" + FormatTool.getDFE4().format(row.logP)
                        + "\t" + FormatTool.getDFE4().format(Math.exp(row.logP))
                        + "\t" + FormatTool.getDFE4().format(row.qValue)
                        + "\t" + row.accepted + "\t" + row.acceptanceOrder
                        + "\tnot_assessed\t" + sameGenus.getOrDefault(row.genome, ""));
            }
        }
    }

    static void writeGenomeScoresRuninfo(File out, File scoresFile, GenomeEvidenceTest.Result evidence)
            throws IOException {
        GenomeEvidenceTest.Calibration c = evidence.calibration;
        try (PrintWriter pw = new PrintWriter(out)) {
            pw.println("=================================================================================");
            pw.println("genome_scores run info (DDA MAG open search, expected false counts)");
            pw.println("=================================================================================");
            pw.println();
            pw.println("output_scores\t" + scoresFile.getAbsolutePath());
            pw.println("runs\t" + c.runs);
            pw.println("candidate_genomes\t" + c.candidateGenomes);
            pw.println();
            pw.println("--- Step 1: peptide FDR per class (spectrum winners; unit = base peptide, I = L) ---");
            pw.println("peptide_fdr\t" + c.peptideFdr);
            pw.println("class\ttarget_psms\tdecoy_psms\ttarget_peptides_accepted\tdecoy_peptides_accepted");
            for (GenomeEvidenceTest.PeptideClass cls : GenomeEvidenceTest.PeptideClass.values()) {
                int k = cls.ordinal();
                pw.println(cls.name().toLowerCase() + "\t" + c.psms[k][0] + "\t" + c.psms[k][1] + "\t"
                        + c.acceptedUnits[k][0] + "\t" + c.acceptedUnits[k][1]);
            }
            pw.println();
            pw.println("--- Step 2: uniqueness against each run's database ---");
            pw.println("status\ttarget_evidence_peptides\tdecoy_peptides");
            Set<String> statuses = new TreeSet<>(c.targetUniqueness.keySet());
            statuses.addAll(c.decoyUniqueness.keySet());
            for (String status : statuses) {
                pw.println(status + "\t" + c.targetUniqueness.getOrDefault(status, 0) + "\t"
                        + c.decoyUniqueness.getOrDefault(status, 0));
            }
            pw.println("F_decoy_evidence_peptides_unique\t" + c.decoyEvidencePeptides);
            pw.println("window_K\t" + c.windowK);
            pw.println("window_size\t" + c.windowSize + "\t(of " + c.windowAvailable
                    + " unique open-search decoy peptides of any class)");
            pw.println("window_edge_score\t" + formatOrNa(c.windowEdgeScore));
            pw.println();
            pw.println("--- Step 3: sequential near-neighbour exclusion ---");
            pw.println("max_substitutions\t" + c.maxSubstitutions);
            pw.println("genome_fdr\t" + c.genomeFdr);
            pw.println("target_genomes_accepted\t" + c.targetGenomesAccepted);
            pw.println("running_vs_final_q_mismatches\t" + c.runningVersusFinalMismatches);
            pw.println();
            pw.println("--- Calibration: reversed genomes through the same test ---");
            pw.println("decoy_genomes_accepted\t" + c.decoyGenomesAccepted + "\t(expected near 0)");
            pw.println("chi_square_accepted_decoys_vs_F_pi\t" + formatOrNa(c.chiSquare)
                    + "\tdf\t" + c.chiSquareDf + "\tp\t" + formatOrNa(c.chiSquareP));
            if (c.psms[GenomeEvidenceTest.PeptideClass.CLOSED.ordinal()][0]
                    + c.psms[GenomeEvidenceTest.PeptideClass.CLOSED.ordinal()][1] > 0) {
                // F includes closed-stage decoys, pi_g comes from open-search decoys: only a
                // failure of the open-stage reading points at pi_g.
                pw.println("chi_square_open_stage_decoys_vs_F_pi\t" + formatOrNa(c.chiSquareOpen)
                        + "\tdf\t" + c.chiSquareOpenDf + "\tp\t" + formatOrNa(c.chiSquareOpenP)
                        + "\tF_open\t" + c.openDecoyEvidencePeptides);
            }
            pw.println();
            pw.println("--- Step 5: picked protein FDR within accepted genomes (reported, not a filter) ---");
            pw.println("target_proteins_scored\t" + c.targetProteinsScored);
            pw.println("decoy_proteins_scored\t" + c.decoyProteinsScored);
            pw.println("target_proteins_q<=" + c.proteinFdr + "\t" + c.targetProteinsAtFdr);
        }
    }

    private static String formatOrNa(double v) {
        return Double.isNaN(v) ? "NA" : FormatTool.getDFE4().format(v);
    }

    private static void writeGenomeHtml(File out, Set<String> genomes, Map<String, Lineage> lineageMap)
            throws IOException {
        HashMap<String, String[]> genomeMap = new HashMap<>();
        TaxonomyRanks[] mainRanks = TaxonomyRanks.getMainRanks();
        for (String genome : genomes) {
            Lineage lineage = lineageMap.get(genome);
            if (lineage == null) lineage = Lineage.empty();
            ArrayList<String> content = new ArrayList<>();
            content.add(genome);
            content.add(lineage.name);
            content.add(lineage.rank);
            for (TaxonomyRanks rank : mainRanks) {
                content.add(lineage.byRank.getOrDefault(rank.getName(), ""));
            }
            genomeMap.put(genome, content.toArray(new String[0]));
        }
        GenomeHtmlWriter writer = new GenomeHtmlWriter(out);
        writer.write(genomeMap);
        writer.close();
    }

    private static void writeGenomeRankScoresRuninfo(File out, File rankScoresFile,
                                                     DdaMagGenomeScoringService.DdaGenomeScoreRuninfo runinfo)
            throws IOException {
        try (PrintWriter pw = new PrintWriter(out)) {
            pw.println("=================================================================================");
            pw.println("genome_rank_scores run info - diagnostic report (DDA MAG, DIA-style scoring)");
            pw.println("=================================================================================");
            pw.println();
            pw.println("--- Input / Output ---");
            pw.println("workflow\tDDA MAG post-analysis");
            pw.println("output_scores\t" + rankScoresFile.getAbsolutePath());
            pw.println();
            pw.println("--- Peptide counts ---");
            pw.println("n_peptides_target\t" + runinfo.nTargetPeptides + "\t(after PSM collapse, unique peptides only m=1)");
            pw.println("n_peptides_decoy\t" + runinfo.nDecoyPeptides);
            pw.println();
            pw.println("--- Peptide sharing (target) ---");
            pw.println("n_unique_peptides\t" + runinfo.nUniquePeptides + "\t(m=1 genome)");
            pw.println("n_shared_peptides\t" + runinfo.nSharedPeptides + "\t(m>1 genomes, skipped in unique-only mode)");
            pw.print("m_distribution\t");
            List<Integer> mKeys = new ArrayList<>(runinfo.mDistribution.keySet());
            Collections.sort(mKeys);
            for (int i = 0; i < mKeys.size(); i++) {
                int m = mKeys.get(i);
                if (i > 0) pw.print(";");
                pw.print(m + ":" + runinfo.mDistribution.get(m));
            }
            pw.println();
            pw.println();
            pw.println("--- Scoring parameters ---");
            pw.println("top_k\t" + runinfo.k);
            pw.println("decoy_resampling_B\t" + runinfo.B);
            pw.println("score_cap\t" + runinfo.cap);
            pw.println("intensity_log_min\t" + runinfo.minLog);
            pw.println("intensity_log_max\t" + runinfo.maxLog);
            pw.println("intensity_log_range\t" + runinfo.logRange);
            pw.println();
            pw.println("--- Genome counts ---");
            pw.println("n_target_genomes\t" + (runinfo.SMap != null ? runinfo.SMap.size() : 0));
            pw.println("n_decoy_pool_entries\t" + runinfo.decoyPoolSize);
            pw.println();
            double threshold = DdaMagProteinInferenceService.DEFAULT_GENOME_Q_VALUE_THRESHOLD;
            pw.println("--- Decoy genome calibration (reversed proteomes through the same test) ---");
            pw.println("Estimated genome FDR at a p threshold = decoy genomes passing / target genomes passing.");
            int targetScorePass = 0;
            if (runinfo.pempArr != null) {
                for (double p : runinfo.pempArr) {
                    if (p <= threshold) targetScorePass++;
                }
            }
            int decoyScorePass = countAtOrBelow(runinfo.decoyGenomeScorePMap, threshold);
            pw.println("score_test\ttarget_genomes_p<=" + threshold + "\t" + targetScorePass
                    + "\tdecoy_genomes_p<=" + threshold + "\t" + decoyScorePass + " of "
                    + runinfo.decoyGenomeScorePMap.size());
            pw.println("reversed_genome\tn_peptides\tp");
            List<String> decoyGenomes = new ArrayList<>(runinfo.decoyGenomeScorePMap.keySet());
            decoyGenomes.sort(Comparator.comparingDouble(runinfo.decoyGenomeScorePMap::get));
            for (String g : decoyGenomes) {
                pw.println(g + "\t" + runinfo.decoyGenomeScoreNMap.getOrDefault(g, 0) + "\t"
                        + FormatTool.getDF4().format(runinfo.decoyGenomeScorePMap.get(g)));
            }
        }
    }

    private static int countAtOrBelow(Map<String, Double> pValues, double threshold) {
        int n = 0;
        for (double p : pValues.values()) {
            if (p <= threshold) n++;
        }
        return n;
    }

    private static Map<String, double[]> alignIntensityToGenomes(Map<String, double[]> values,
            Map<String, double[]> genomeIntensityMap) {
        HashMap<String, double[]> aligned = new HashMap<>();
        for (String genome : genomeIntensityMap.keySet()) {
            double[] template = genomeIntensityMap.get(genome);
            int n = template == null ? 0 : template.length;
            double[] src = values == null ? null : values.get(genome);
            if (src == null) {
                aligned.put(genome, new double[n]);
                continue;
            }
            if (src.length == n) {
                aligned.put(genome, src);
                continue;
            }
            double[] copy = new double[n];
            System.arraycopy(src, 0, copy, 0, Math.min(n, src.length));
            aligned.put(genome, copy);
        }
        return aligned;
    }

    private static String[] orderGenomes(Map<String, double[]> genomeIntensityMap, List<String> orderedGenomes) {
        ArrayList<String> genomes = new ArrayList<>();
        if (orderedGenomes != null) {
            for (String genome : orderedGenomes) {
                if (genomeIntensityMap.containsKey(genome) && !genomes.contains(genome)) {
                    genomes.add(genome);
                }
            }
        }
        for (String genome : genomeIntensityMap.keySet()) {
            if (!genomes.contains(genome)) {
                genomes.add(genome);
            }
        }
        if (orderedGenomes == null || orderedGenomes.isEmpty()) {
            Collections.sort(genomes);
        }
        return genomes.toArray(new String[0]);
    }

    private static void writeTaxaReport(File out,
                                       Map<String, double[]> genomeIntensityMap,
                                       Map<String, Lineage> lineageMap,
                                       String[] intensityTitles) throws IOException {
        // Aggregate intensities at main ranks using GTDB-derived lineage from taxa TSV.
        TaxonomyRanks[] mainRanks = TaxonomyRanks.getMainRanks();
        HashMap<String, double[]> taxaIntensity = new HashMap<>();
        HashMap<String, String> taxaRank = new HashMap<>();

        for (String g : genomeIntensityMap.keySet()) {
            double[] gi = genomeIntensityMap.get(g);
            Lineage lin = lineageMap.get(g);
            if (lin == null) continue;

            for (TaxonomyRanks r : mainRanks) {
                String name = lin.byRank.get(r.getName());
                if (name == null || name.isBlank()) continue;
                String key = r.getName() + "\t" + name;
                double[] acc = taxaIntensity.computeIfAbsent(key, k -> new double[gi.length]);
                for (int i = 0; i < gi.length; i++) acc[i] += gi[i];
                taxaRank.put(key, r.getName());
            }
        }

        try (PrintWriter pw = new PrintWriter(out)) {
            StringBuilder h = new StringBuilder();
            h.append("Taxon\tRank");
            for (String exp : intensityTitles) {
                if (exp == null || exp.isBlank()) continue;
                h.append("\t").append(exp.startsWith("Intensity ") ? exp : ("Intensity " + exp));
            }
            pw.println(h);

            String[] keys = taxaIntensity.keySet().toArray(new String[0]);
            Arrays.sort(keys, String::compareToIgnoreCase);
            for (String k : keys) {
                String[] parts = k.split("\t", 2);
                String rank = parts[0];
                String taxon = parts.length > 1 ? parts[1] : "";
                double[] inten = taxaIntensity.get(k);

                StringBuilder sb = new StringBuilder();
                sb.append(taxon).append("\t").append(rank);
                if (inten != null) {
                    for (double v : inten) {
                        sb.append("\t").append(FormatTool.getDFE4().format(v));
                    }
                }
                pw.println(sb);
            }
        }
    }

    private static final class Lineage {
        final String name; // most specific label
        final String rank; // most specific rank name
        final HashMap<String, String> byRank; // rankName -> label
        final String rawLineage;

        private Lineage(String name, String rank, HashMap<String, String> byRank, String rawLineage) {
            this.name = name;
            this.rank = rank;
            this.byRank = byRank;
            this.rawLineage = rawLineage;
        }

        static Lineage empty() {
            return new Lineage("", "", new HashMap<>(), "");
        }

        static Lineage parse(String lineage) {
            if (lineage == null) return empty();
            String[] parts = lineage.split(";");
            HashMap<String, String> map = new HashMap<>();
            String bestName = "";
            String bestRank = "";
            for (String p : parts) {
                p = p.trim();
                if (p.length() < 3) continue;
                String label = p;
                // handle GTDB-style prefixes like d__/p__/...
                if (p.contains("__")) {
                    int idx = p.indexOf("__");
                    if (idx >= 0 && idx + 2 < p.length()) {
                        label = p.substring(idx + 2).trim();
                    }
                }
                if (label.isEmpty()) continue;

                char r = Character.toLowerCase(p.charAt(0));
                String rankName = rankFromPrefix(r);
                if (!rankName.isEmpty()) {
                    map.put(rankName, label);
                    bestName = label;
                    bestRank = rankName;
                }
            }
            if (!map.containsKey(TaxonomyRanks.Superkingdom.getName()) && hasLowerRank(map)) {
                map.put(TaxonomyRanks.Superkingdom.getName(), "Bacteria");
            }
            return new Lineage(bestName, bestRank, map, lineage);
        }

        private static boolean hasLowerRank(HashMap<String, String> map) {
            return map.containsKey(TaxonomyRanks.Phylum.getName())
                    || map.containsKey(TaxonomyRanks.Class.getName())
                    || map.containsKey(TaxonomyRanks.Order.getName())
                    || map.containsKey(TaxonomyRanks.Family.getName())
                    || map.containsKey(TaxonomyRanks.Genus.getName())
                    || map.containsKey(TaxonomyRanks.Species.getName());
        }

        private static String rankFromPrefix(char r) {
            switch (r) {
                case 'd': return TaxonomyRanks.Superkingdom.getName();
                case 'p': return TaxonomyRanks.Phylum.getName();
                case 'c': return TaxonomyRanks.Class.getName();
                case 'o': return TaxonomyRanks.Order.getName();
                case 'f': return TaxonomyRanks.Family.getName();
                case 'g': return TaxonomyRanks.Genus.getName();
                case 's': return TaxonomyRanks.Species.getName();
                default: return "";
            }
        }
    }
}

