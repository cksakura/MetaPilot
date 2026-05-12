package uk.ac.quadram.metapilot.core.runner;

import uk.ac.quadram.metapilot.core.taxonomy.TaxonomyRanks;
import uk.ac.quadram.metapilot.core.tools.FormatTool;
import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

/**
 * Headless taxonomy TSV exporter (no HTML).
 *
 * Writes:
 * - taxonomy_analysis/Genome.tsv
 * - taxonomy_analysis/Taxa.tsv (DIA-style scoring)
 * - taxonomy_analysis/Taxa_report.tsv
 * - taxonomy_analysis/genome_rank_scores.tsv (DIA-style scoring)
 *
 * This intentionally avoids the GUI-only dependency on a packaged "gtdb_ncbi" mapping resource.
 * Instead, it uses the MAG DB's original taxa TSV directly (from MAGDB_1_0.json version metadata).
 */
public final class MagTaxonomyTsvExporter {

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
        if (cb == null) cb = new ProgressCallback.NoOp();
        if (taxonomyAnalysisDir == null) throw new IOException("taxonomyAnalysisDir is null");
        taxonomyAnalysisDir.mkdirs();

        if (originalTaxaTsv == null || !originalTaxaTsv.exists()) {
            throw new IOException("MAG taxa TSV not found: " + (originalTaxaTsv == null ? "null" : originalTaxaTsv.getAbsolutePath()));
        }

        if (genomeIntensityMap == null || genomeIntensityMap.isEmpty()) {
            throw new IOException("genomeIntensityMap is empty (no genomes to export)");
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
        }

        File genomeTsv = new File(taxonomyAnalysisDir, "Genome.tsv");
        cb.onMessage("Writing taxonomy: " + genomeTsv.getName());
        writeGenomeTsv(genomeTsv, genomeIntensityMap, genomeScoreMap, genomeDiaStats, genomePsmCountMap, lineageMap, gtdbToNcbi,
                genomeSelectionLabel,
                totalPepCountMap, razorPepCountMap, uniquePepCountMap,
                totalPepWithILMap, razorPepWithILMap, uniquePepWithILMap,
                intensityTitles);

        boolean useDiaStyle = genomeDiaStats != null && !genomeDiaStats.isEmpty();
        if (useDiaStyle) {
            File taxaTsv = new File(taxonomyAnalysisDir, "Taxa.tsv");
            File taxaReportTsv = new File(taxonomyAnalysisDir, "Taxa_report.tsv");
            cb.onMessage("Writing taxonomy: " + taxaTsv.getName());
            writeDiaStyleTaxaTables(taxaTsv, taxaReportTsv, genomeIntensityMap, genomeScoreMap, genomeDiaStats,
                    lineageMap, intensityTitles);
            File rankScoresFile = new File(taxonomyAnalysisDir, "genome_rank_scores.tsv");
            cb.onMessage("Writing taxonomy: " + rankScoresFile.getName());
            writeGenomeRankScores(rankScoresFile, genomeDiaStats, totalPepCountMap, razorPepCountMap);
        } else {
            File taxaReportTsv = new File(taxonomyAnalysisDir, "Taxa_report.tsv");
            cb.onMessage("Writing taxonomy: " + taxaReportTsv.getName());
            writeTaxaReport(taxaReportTsv, genomeIntensityMap, lineageMap, intensityTitles);
        }
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

    private static void writeGenomeTsv(File out,
                                      Map<String, double[]> genomeIntensityMap,
                                      Map<String, Double> genomeScoreMap,
                                      Map<String, double[]> genomeDiaStats,
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
                                      String[] intensityTitles) throws IOException {
        TaxonomyRanks[] mainRanks = TaxonomyRanks.getMainRanks();
        String[] genomes = genomeIntensityMap.keySet().toArray(new String[0]);
        Arrays.sort(genomes);

        try (PrintWriter pw = new PrintWriter(out)) {
            StringBuilder h = new StringBuilder();
            h.append("Genome\t");
            boolean useDiaStyleHeader = (genomeDiaStats != null && !genomeDiaStats.isEmpty());
            if (useDiaStyleHeader) {
                h.append("S\t");
                h.append("p_emp\t");
                h.append("q_value\t");
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
                if (st != null && st.length >= 3) {
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
            pw.println("Rank\tGenome\tS\tn_peptides\tp_emp\tq_value");
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
            return new Lineage(bestName, bestRank, map, lineage);
        }

        private static String rankFromPrefix(char r) {
            switch (r) {
                case 'd': return TaxonomyRanks.Subkingdom.getName();
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

