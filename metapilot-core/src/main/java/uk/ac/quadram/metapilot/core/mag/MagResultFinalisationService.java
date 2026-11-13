package uk.ac.quadram.metapilot.core.mag;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.*;
import uk.ac.quadram.metapilot.core.runner.MagFinalGenomeRankingService;
import uk.ac.quadram.metapilot.task.io.pep.MetaPeptide;

/** One frozen genome decision shared by DIA/DDA and GUI/CLI. Scores are never
 * recalibrated after selection: pruning changes attribution, not the hypothesis test. */
public final class MagResultFinalisationService {
    private MagResultFinalisationService() {}

    public static void writeCatalogueContext(File resultDir, uk.ac.quadram.metapilot.task.mag.MagDbItem catalogue) throws IOException {
        Files.createDirectories(resultDir.toPath());
        File file = new File(resultDir, "list_consistency_context.json");
        com.google.gson.JsonObject context = file.isFile()
                ? com.google.gson.JsonParser.parseString(Files.readString(file.toPath())).getAsJsonObject()
                : new com.google.gson.JsonObject();
        context.addProperty("catalogue_identifier", catalogue.getIdentifier());
        context.addProperty("catalogue_id", catalogue.getCatalogueID());
        Files.writeString(file.toPath(), new com.google.gson.GsonBuilder().setPrettyPrinting().create().toJson(context), StandardCharsets.UTF_8);
    }

    public static final class Selection {
        public final Set<String> genomes;
        public final Map<String, String> reasons;
        public final Map<String, double[]> stats;
        private final MagFinalGenomeRankingService.Request ranking;
        public final List<String[]> excludedPeptides = new ArrayList<>();
        private Selection(Set<String> genomes, Map<String, String> reasons,
                          Map<String, double[]> stats, MagFinalGenomeRankingService.Request ranking) {
            this.genomes = genomes;
            this.reasons = Collections.unmodifiableMap(reasons);
            Map<String, double[]> snapshot = new HashMap<>();
            if (stats != null) stats.forEach((g, s) -> snapshot.put(g, s == null ? null : s.clone()));
            this.stats = Collections.unmodifiableMap(snapshot);
            this.ranking = ranking;
        }
        public void exclude(String sequence, String modified, String[] proteins) {
            excludedPeptides.add(new String[]{sequence, modified,
                    proteins == null ? "" : String.join(";", proteins), "NO_RETAINED_PROTEIN"});
        }
        public void write(File resultDir) throws IOException {
            File tax = new File(resultDir, "taxonomy_analysis");
            Files.createDirectories(tax.toPath());
            ranking.taxonomyAnalysisDir = tax;
            MagFinalGenomeRankingService.write(ranking);
            try (PrintWriter out = new PrintWriter(new File(tax, "genome_final_decisions.tsv"), StandardCharsets.UTF_8)) {
                out.println("Genome\tSelected\tReason\tp_emp\tq_value");
                for (String g : new TreeSet<>(reasons.keySet())) {
                    double[] s = stats == null ? null : stats.get(g);
                    out.println(g + "\t" + genomes.contains(g) + "\t" + reasons.get(g) + "\t"
                            + (s != null && s.length > 1 ? s[1] : "") + "\t"
                            + (s != null && s.length > 2 ? s[2] : ""));
                }
            }
            File evidence = new File(resultDir, MagResultLayoutUtil.EXCLUDED_EVIDENCE_DIR_NAME);
            Files.createDirectories(evidence.toPath());
            try (PrintWriter out = new PrintWriter(new File(evidence, "peptides.tsv"), StandardCharsets.UTF_8)) {
                out.println("Base sequence\tModified sequence\tOriginal proteins\tReason");
                for (String[] row : excludedPeptides) out.println(String.join("\t", Arrays.stream(row)
                        .map(v -> v == null ? "" : v.replace('\t', ' ').replace('\n', ' ')).toArray(String[]::new)));
            }
        }
    }

    /** Carries the frozen decision with the rebuilt abundance maps through existing export wrappers. */
    public static final class FinalGenomeIntensityMap extends HashMap<String, double[]> {
        public final Selection selection;
        public FinalGenomeIntensityMap(Map<String, double[]> values, Selection selection) {
            super(values);
            keySet().retainAll(selection.genomes);
            for (String genome : selection.genomes) {
                double[] provisional = selection.ranking.genomeIntensityMap.get(genome);
                putIfAbsent(genome, new double[provisional == null ? 0 : provisional.length]);
            }
            this.selection = selection;
        }
    }

    public static Selection select(Map<String, double[]> intensities, Map<String, double[]> stats,
            Set<String> acceptedWithoutStats, Map<String, Integer> psms,
            Map<String, ? extends Set<String>> total, Map<String, ? extends Set<String>> razor,
            Map<String, ? extends Set<String>> unique) {
        Map<String, String> reasons = new TreeMap<>();
        Set<String> candidates = new TreeSet<>(intensities.keySet());
        if (stats != null) candidates.addAll(stats.keySet());
        if (acceptedWithoutStats != null) candidates.addAll(acceptedWithoutStats);
        Map<String, double[]> eligible = new HashMap<>();
        for (String g : candidates) {
            double[] s = stats == null ? null : stats.get(g);
            boolean confident = stats != null && !stats.isEmpty()
                    ? s != null && s.length > 2 && Double.isFinite(s[2]) && s[2] <= 0.01
                    : acceptedWithoutStats != null && acceptedWithoutStats.contains(g);
            reasons.put(g, !confident ? "GENOME_CONFIDENCE_GATE" : "NO_QUANTIFIED_SUPPORT");
            if (confident && intensities.containsKey(g)) eligible.put(g, intensities.get(g).clone());
        }
        MagFinalGenomeRankingService.Request r = new MagFinalGenomeRankingService.Request();
        r.genomeIntensityMap = eligible;
        r.genomePsmCountMap = psms == null ? null : new HashMap<>(psms);
        r.totalPepCountMap = total;
        r.razorPepCountMap = razor;
        r.uniquePepCountMap = unique;
        Set<String> selected = MagFinalGenomeRankingService.selectGenomes(r);
        for (String g : eligible.keySet()) reasons.put(g, selected.contains(g) ? "SELECTED" : "NO_INCREMENTAL_COVERAGE");
        return new Selection(selected, reasons, stats, r);
    }

    public static String genome(String protein) {
        int i = protein.indexOf('_');
        return i < 0 ? protein : protein.substring(0, i);
    }
    public static String[] retainedProteins(String[] proteins, String prefix, String decoy, Set<String> genomes) {
        if (proteins == null) return new String[0];
        return Arrays.stream(proteins).filter(p -> p != null && !p.isBlank())
                .filter(p -> decoy == null || decoy.isEmpty() || !p.startsWith(decoy))
                .filter(p -> prefix == null || prefix.isEmpty() || !p.startsWith(prefix) || genomes.contains(genome(p)))
                .distinct().sorted().toArray(String[]::new);
    }
    /** Preserves array positions for existing DDA razor maps. */
    public static MetaPeptide[] prune(MetaPeptide[] peptides, String prefix, String decoy, Selection selection) {
        MetaPeptide[] retained = peptides.clone();
        for (int i = 0; i < retained.length; i++) {
            MetaPeptide p = retained[i];
            if (p == null) continue;
            String[] proteins = retainedProteins(p.getProteins(), prefix, decoy, selection.genomes);
            if (proteins.length == 0) {
                selection.exclude(p.getSequence(), p.getModSeq(), p.getProteins());
                retained[i] = null;
            } else p.setProteins(proteins);
        }
        return retained;
    }
}
