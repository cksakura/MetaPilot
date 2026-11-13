package uk.ac.quadram.metapilot.core.dia;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Set;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNParquetReader;
import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNPrecursor;
import uk.ac.quadram.metapilot.task.dia.DeepIterativeCommonSupport;

/**
 * Shared DIA deep post-processing support.
 *
 * <p>Deep workflows build a non-redundant final combined FASTA from round peptides.
 * Before final scoring/export, the peptides from {@code combined.parquet} must be
 * mapped back to every protein observed for that peptide across the iterative rounds.
 *
 * <p>That mapping can either be derived directly from the round parquet files
 * ({@link Request#deepTogether} / {@link Request#deepPerFile}, reading potentially many
 * files under a {@code marker} folder), or from a single precomputed
 * {@code peptide -> proteins} map ({@link Request#prebuiltMap}) - typically produced once,
 * locally, via {@link #writeMapTsv} from the very same round parquets, then transferred as
 * one small file instead of the whole round-parquet folder tree. See
 * {@link #readMapTsv(File)} for the read side.
 */
public final class DiaDeepProteinMappingService {

    private static final Logger DEFAULT_LOGGER = LogManager.getLogger(DiaDeepProteinMappingService.class);
    private static final int DEFAULT_MAX_ROUNDS = 100;

    private DiaDeepProteinMappingService() {
    }

    public static final class Request {
        public final File roundRoot;
        public final boolean perFileLayout;
        public final int lastRoundHint;
        public final int maxRounds;
        public final String[] diaFiles;
        public final String cut;
        /** When non-null, used as-is instead of reading anything under {@link #roundRoot}. */
        public final HashMap<String, HashSet<String>> prebuiltMap;

        private Request(File roundRoot, boolean perFileLayout, int lastRoundHint, int maxRounds,
                        String[] diaFiles, String cut, HashMap<String, HashSet<String>> prebuiltMap) {
            this.roundRoot = roundRoot;
            this.perFileLayout = perFileLayout;
            this.lastRoundHint = lastRoundHint;
            this.maxRounds = maxRounds > 0 ? maxRounds : DEFAULT_MAX_ROUNDS;
            this.diaFiles = diaFiles;
            this.cut = cut != null && !cut.isBlank() ? cut : "KR";
            this.prebuiltMap = prebuiltMap;
        }

        public static Request deepTogether(File markerDir, int lastRoundHint, String[] diaFiles, String cut) {
            return new Request(markerDir, false, lastRoundHint, DEFAULT_MAX_ROUNDS, diaFiles, cut, null);
        }

        public static Request deepPerFile(File markerDeepDir, String[] diaFiles, String cut) {
            return new Request(markerDeepDir, true, 0, DEFAULT_MAX_ROUNDS, diaFiles, cut, null);
        }

        /**
         * Wraps an already-computed peptide-to-proteins map (e.g. read back via
         * {@link DiaDeepProteinMappingService#readMapTsv(File)}) so {@link #remapPrecursors}
         * uses it directly - no round parquet files are read or required.
         */
        public static Request prebuiltMap(HashMap<String, HashSet<String>> map) {
            return new Request(null, false, 0, DEFAULT_MAX_ROUNDS, null, "KR", map);
        }
    }

    public static final class RemapResult {
        public final DiaNNPrecursor[] precursors;
        public final int mappedPeptides;
        public final int proteinsBefore;
        public final int proteinsAfter;

        private RemapResult(DiaNNPrecursor[] precursors, int mappedPeptides, int proteinsBefore, int proteinsAfter) {
            this.precursors = precursors;
            this.mappedPeptides = mappedPeptides;
            this.proteinsBefore = proteinsBefore;
            this.proteinsAfter = proteinsAfter;
        }
    }

    public static RemapResult remapPrecursors(DiaNNPrecursor[] precursors, Request request,
                                              Logger logger, String taskName) {
        Logger log = logger != null ? logger : DEFAULT_LOGGER;
        String name = taskName != null && !taskName.isBlank() ? taskName : "DIA deep protein remapping";
        if (precursors == null) {
            return new RemapResult(new DiaNNPrecursor[0], 0, 0, 0);
        }
        boolean hasPrebuiltMap = request != null && request.prebuiltMap != null;
        if (request == null || (!hasPrebuiltMap && (request.roundRoot == null || !request.roundRoot.isDirectory()))) {
            throw new IllegalStateException(name + ": recorded round protein mappings are unavailable");
        }

        HashMap<String, HashSet<String>> pepProMap = buildPeptideProteinMap(request, log, name);
        if (pepProMap.isEmpty()) throw new IllegalStateException(name + ": no recorded round protein mappings");
        pepProMap = mapFinalPeptides(precursors, pepProMap);

        int proteinsBefore = countDistinctProteins(precursors);
        int mappedPrecursors = 0;
        for (DiaNNPrecursor precursor : precursors) {
            if (precursor == null) {
                continue;
            }
            if (precursor.getIsDecoy()) continue; // preserve the original scoring background
            String sequence = precursor.getSeqString();
            if (sequence == null || sequence.isEmpty()) {
                continue;
            }
            HashSet<String> proteins = pepProMap.get(sequence);
            if (proteins == null || proteins.isEmpty()) {
                if (!precursor.getIsDecoy() && precursor.getGlobalQvalue() < 0.01)
                    throw new IllegalStateException(name + ": missing round protein mapping for " + sequence);
                continue;
            }
            precursor.setProGroups(proteins.stream().sorted().toArray(String[]::new));
            mappedPrecursors++;
        }
        int proteinsAfter = countDistinctProteins(precursors);
        log.info("{}: deep peptide-protein remap applied; mapping peptides={}, mapped precursors={}, proteins before={}, proteins after={}",
                name, pepProMap.size(), mappedPrecursors, proteinsBefore, proteinsAfter);
        return new RemapResult(precursors, pepProMap.size(), proteinsBefore, proteinsAfter);
    }

    /** The final peptide FASTA can yield shorter cleavage products. Restore all
     * associations from containing recorded peptides, including when an exact
     * match also exists. This never falls back to arbitrary combined protein IDs. */
    static HashMap<String, HashSet<String>> mapFinalPeptides(DiaNNPrecursor[] precursors,
            Map<String, HashSet<String>> recorded) {
        HashSet<String> queries = new HashSet<>();
        java.util.TreeSet<Integer> lengths = new java.util.TreeSet<>();
        for (DiaNNPrecursor p : precursors) {
            if (p != null && !p.getIsDecoy() && p.getSeqString() != null && !p.getSeqString().isEmpty()) {
                queries.add(p.getSeqString());
                lengths.add(p.getSeqString().length());
            }
        }
        HashMap<String, HashSet<String>> result = new HashMap<>();
        for (Map.Entry<String, HashSet<String>> e : recorded.entrySet()) {
            String source = e.getKey();
            for (int length : lengths) {
                if (length > source.length()) break;
                for (int start = 0; start <= source.length() - length; start++) {
                    String peptide = source.substring(start, start + length);
                    if (queries.contains(peptide))
                        result.computeIfAbsent(peptide, k -> new HashSet<>()).addAll(e.getValue());
                }
            }
        }
        return result;
    }

    public static HashMap<String, HashSet<String>> buildPeptideProteinMap(Request request,
                                                                          Logger logger,
                                                                          String taskName) {
        Logger log = logger != null ? logger : DEFAULT_LOGGER;
        String name = taskName != null && !taskName.isBlank() ? taskName : "DIA deep protein remapping";
        if (request != null && request.prebuiltMap != null) {
            log.info("{}: using prebuilt peptide-protein map ({} peptides)", name, request.prebuiltMap.size());
            return request.prebuiltMap;
        }
        HashMap<String, HashSet<String>> pepProMap = new HashMap<>();
        if (request == null || request.roundRoot == null || !request.roundRoot.isDirectory()) {
            return pepProMap;
        }

        if (request.perFileLayout) {
            File markerRoot = new File(request.roundRoot, "marker");
            if (!markerRoot.isDirectory()) {
                log.warn("{}: DEEP_PER_FILE remap expected marker folder under {}",
                        name, request.roundRoot.getAbsolutePath());
                return pepProMap;
            }
            File[] perFileDirs = markerRoot.listFiles(File::isDirectory);
            if (perFileDirs == null || perFileDirs.length == 0) {
                return pepProMap;
            }
            Arrays.sort(perFileDirs, (a, b) -> a.getName().compareToIgnoreCase(b.getName()));
            LinkedHashMap<File, HashSet<String>> genomesByDir = new LinkedHashMap<>();
            HashSet<String> searchedInAllRounds = new HashSet<>();
            for (File perFileDir : perFileDirs) {
                int lastRound = detectLastRound(perFileDir, request);
                HashMap<String, HashSet<String>> rounds = lastRound > 0
                        ? buildPeptideProteinMapFromRounds(perFileDir, lastRound, request, log, name)
                        : new HashMap<>();
                HashSet<String> genomes = genomesOf(rounds);
                genomesByDir.put(perFileDir, genomes);
                searchedInAllRounds.addAll(genomes);
                mergeInto(pepProMap, rounds);
            }
            // BALANCED clustering can leave a sample with initial Marker
            // evidence but no iterative rounds of its own.
            for (Map.Entry<File, HashSet<String>> e : genomesByDir.entrySet()) {
                Set<String> searched = e.getValue().isEmpty() ? searchedInAllRounds : e.getValue();
                mergeInto(pepProMap, markerEvidence(e.getKey(), searched, request, log, name));
            }
            log.info("{}: DEEP_PER_FILE remap collected {} peptides from all per-file round parquets",
                    name, pepProMap.size());
            return pepProMap;
        }

        int lastRound = request.lastRoundHint > 0 ? request.lastRoundHint : detectLastRound(request.roundRoot, request);
        HashMap<String, HashSet<String>> rounds = lastRound > 0
                ? buildPeptideProteinMapFromRounds(request.roundRoot, lastRound, request, log, name)
                : new HashMap<>();
        mergeInto(pepProMap, markerEvidence(request.roundRoot, rounds, request, log, name));
        mergeInto(pepProMap, rounds);
        log.info("{}: deep remap collected {} peptides from round parquets", name, pepProMap.size());
        return pepProMap;
    }

    /**
     * Marker-search evidence ({@code marker.parquet}) of one round folder, restricted to the genomes searched in its
     * rounds. The marker search runs against the Marker4 library of the whole catalogue, so without this restriction its
     * peptides would bring proteins of genomes that were never selected into the final library and the export.
     */
    public static HashMap<String, HashSet<String>> markerEvidence(File dir, Map<String, HashSet<String>> rounds,
                                                           Request request, Logger log, String name) {
        return markerEvidence(dir, genomesOf(rounds), request, log, name);
    }

    /**
     * Marker-search evidence of one folder, restricted to {@code searchedGenomes}. A sample without rounds of its own
     * (BALANCED clustering) is searched through the rounds of its representative, so callers pass the genomes searched
     * in all rounds. An empty set (no rounds at all) falls back to the folder's selected genomes.
     */
    public static HashMap<String, HashSet<String>> markerEvidence(File dir, Set<String> searchedGenomes,
                                                           Request request, Logger log, String name) {
        HashMap<String, HashSet<String>> marker = new HashMap<>();
        File parquet = new File(dir, "marker.parquet");
        if (!parquet.isFile()) {
            return marker;
        }
        addRecordedParquet(marker, parquet, request, name);
        HashSet<String> genomes = new HashSet<>(searchedGenomes);
        if (genomes.isEmpty()) {
            genomes.addAll(readSelectedGenomes(new File(dir, "genome_selected.tsv"), name));
        }
        int before = marker.size();
        restrictToGenomes(marker, genomes);
        log.info("{}: marker evidence of {} kept for {} searched genomes: {} of {} peptides",
                name, dir.getName(), genomes.size(), marker.size(), before);
        return marker;
    }

    public static HashSet<String> genomesOf(Map<String, HashSet<String>> pepProMap) {
        HashSet<String> genomes = new HashSet<>();
        for (HashSet<String> proteins : pepProMap.values()) {
            for (String protein : proteins) {
                genomes.add(genomeOf(protein));
            }
        }
        return genomes;
    }

    /** Drops proteins of other genomes, and peptides left without a protein. */
    static void restrictToGenomes(Map<String, HashSet<String>> pepProMap, java.util.Set<String> genomes) {
        pepProMap.values().forEach(proteins -> proteins.removeIf(p -> !genomes.contains(genomeOf(p))));
        pepProMap.values().removeIf(HashSet::isEmpty);
    }

    /** Same genome rule as DdaMagPostAnalysisService / DdaMagProteinInferenceService. */
    static String genomeOf(String protein) {
        int idx = protein.indexOf('_');
        return idx > 0 ? protein.substring(0, idx) : protein;
    }

    private static HashSet<String> readSelectedGenomes(File tsv, String name) {
        HashSet<String> genomes = new HashSet<>();
        if (!tsv.isFile()) {
            return genomes;
        }
        try (BufferedReader reader = Files.newBufferedReader(tsv.toPath(), StandardCharsets.UTF_8)) {
            String line = reader.readLine(); // header: Genome ...
            while ((line = reader.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty() || line.startsWith("#")) {
                    break;
                }
                String genome = line.split("\t", 2)[0].trim();
                if (!genome.isEmpty()) {
                    genomes.add(genome);
                }
            }
        } catch (IOException e) {
            throw new IllegalStateException(name + ": failed to read " + tsv, e);
        }
        return genomes;
    }

    /**
     * Writes a peptide-to-proteins map as one compact TSV: {@code Peptide\tProtein1;Protein2;...}
     * per line. Intended to be produced once, locally, from the round parquet files (via
     * {@link #buildPeptideProteinMap}) and transferred as a single small file instead of the
     * whole round-parquet folder tree - see {@link #readMapTsv(File)} for the read side.
     */
    public static void writeMapTsv(File out, Map<String, HashSet<String>> pepProMap) throws IOException {
        if (out == null) {
            throw new IOException("output file is null");
        }
        File parent = out.getParentFile();
        if (parent != null && !parent.exists()) {
            Files.createDirectories(parent.toPath());
        }
        try (PrintWriter w = new PrintWriter(out, StandardCharsets.UTF_8)) {
            w.println("Peptide\tProteins");
            if (pepProMap != null) {
                for (Map.Entry<String, HashSet<String>> entry : pepProMap.entrySet()) {
                    String peptide = entry.getKey();
                    HashSet<String> proteins = entry.getValue();
                    if (peptide == null || peptide.isEmpty() || proteins == null || proteins.isEmpty()) {
                        continue;
                    }
                    w.print(peptide);
                    w.print('\t');
                    w.println(String.join(";", proteins));
                }
            }
        }
    }

    /**
     * Reads back a map written by {@link #writeMapTsv}. Pass the result to
     * {@link Request#prebuiltMap(HashMap)} to use it in {@link #remapPrecursors} without
     * touching any round parquet files.
     */
    public static HashMap<String, HashSet<String>> readMapTsv(File in) throws IOException {
        if (in == null || !in.isFile()) {
            throw new IOException("peptide-protein map file not found: " + (in == null ? "null" : in.getAbsolutePath()));
        }
        HashMap<String, HashSet<String>> map = new HashMap<>();
        try (BufferedReader reader = new BufferedReader(new java.io.InputStreamReader(
                new java.io.FileInputStream(in), StandardCharsets.UTF_8))) {
            String line = reader.readLine(); // header, discarded
            while ((line = reader.readLine()) != null) {
                if (line.isEmpty()) {
                    continue;
                }
                int tab = line.indexOf('\t');
                if (tab < 0) {
                    continue;
                }
                String peptide = line.substring(0, tab);
                String rest = line.substring(tab + 1);
                if (peptide.isEmpty() || rest.isEmpty()) {
                    continue;
                }
                HashSet<String> proteins = new HashSet<>(Arrays.asList(rest.split(";")));
                proteins.removeIf(p -> p == null || p.isEmpty());
                if (!proteins.isEmpty()) {
                    map.put(peptide, proteins);
                }
            }
        }
        return map;
    }

    private static HashMap<String, HashSet<String>> buildPeptideProteinMapFromRounds(File markerDir, int lastRound,
                                                                                     Request request, Logger logger,
                                                                                     String taskName) {
        HashMap<String, HashSet<String>> pepProMap = new HashMap<>();
        for (int round = 1; round <= lastRound; round++) {
            File roundDir = new File(markerDir, "round" + round);
            DiaNNPrecursor[] precursors = null;
            File roundParquet = new File(roundDir, "round" + round + ".parquet");
            if (roundParquet.isFile()) {
                addRecordedParquet(pepProMap, roundParquet, request, taskName);
            } else if (request.diaFiles != null && request.diaFiles.length > 0) {
                precursors = DeepIterativeCommonSupport.readPrecursorsFromPerFileParquets(
                        roundDir, null, request.diaFiles, request.cut, logger, taskName);
            }
            addPrecursors(pepProMap, precursors);
        }
        return pepProMap;
    }

    private static void addRecordedParquet(HashMap<String, HashSet<String>> map, File parquet,
                                            Request request, String taskName) {
        if (!parquet.isFile()) return;
        try {
            // This restores sequence-to-protein evidence, not peptide acceptance.
            // Round and final searches have independent q-values. Final scoring
            // and the genome gate still control the exported identifications.
            DiaNNParquetReader reader = new DiaNNParquetReader(parquet.getAbsolutePath(), request.cut);
            addPrecursors(map, reader.getDiaNNPrecursors());
        } catch (Exception e) {
            throw new IllegalStateException(taskName + ": failed to read recorded mapping " + parquet, e);
        }
    }

    private static void addPrecursors(HashMap<String, HashSet<String>> pepProMap, DiaNNPrecursor[] precursors) {
        if (precursors == null) {
            return;
        }
        for (DiaNNPrecursor precursor : precursors) {
            if (precursor == null || precursor.getIsDecoy()) {
                continue;
            }
            String peptide = precursor.getSeqString();
            if (peptide == null || peptide.isEmpty()) {
                continue;
            }
            String[] proteins = precursor.getProGroups();
            if (proteins == null || proteins.length == 0) {
                continue;
            }
            HashSet<String> set = pepProMap.computeIfAbsent(peptide, k -> new HashSet<>());
            for (String protein : proteins) {
                if (protein != null && !protein.isEmpty()) {
                    set.add(protein);
                }
            }
        }
    }

    private static void mergeInto(HashMap<String, HashSet<String>> target, HashMap<String, HashSet<String>> source) {
        if (source == null || source.isEmpty()) {
            return;
        }
        for (java.util.Map.Entry<String, HashSet<String>> entry : source.entrySet()) {
            if (entry.getKey() == null || entry.getKey().isEmpty()
                    || entry.getValue() == null || entry.getValue().isEmpty()) {
                continue;
            }
            target.computeIfAbsent(entry.getKey(), k -> new HashSet<>()).addAll(entry.getValue());
        }
    }

    private static int detectLastRound(File markerDir, Request request) {
        if (markerDir == null || !markerDir.isDirectory()) {
            return -1;
        }
        int lastRound = -1;
        int fileCount = request.diaFiles != null ? request.diaFiles.length : 0;
        for (int round = 1; round <= request.maxRounds; round++) {
            File roundDir = new File(markerDir, "round" + round);
            if (!roundDir.isDirectory()) {
                continue;
            }
            File roundParquet = new File(roundDir, "round" + round + ".parquet");
            if (roundParquet.isFile()
                    || (fileCount > 0 && DeepIterativeCommonSupport.hasRoundPerFileParquets(roundDir, fileCount))) {
                lastRound = round;
            }
        }
        return lastRound;
    }

    private static int countDistinctProteins(DiaNNPrecursor[] precursors) {
        HashSet<String> proteins = new HashSet<>();
        if (precursors == null) {
            return 0;
        }
        for (DiaNNPrecursor precursor : precursors) {
            if (precursor == null) {
                continue;
            }
            String[] proGroups = precursor.getProGroups();
            if (proGroups == null) {
                continue;
            }
            for (String protein : proGroups) {
                if (protein != null && !protein.isEmpty()) {
                    proteins.add(protein);
                }
            }
        }
        return proteins.size();
    }
}
