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
import java.util.Map;

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
            log.debug("{}: deep peptide-protein remap skipped; round root is unavailable", name);
            return new RemapResult(precursors, 0, countDistinctProteins(precursors), countDistinctProteins(precursors));
        }

        HashMap<String, HashSet<String>> pepProMap = buildPeptideProteinMap(request, log, name);
        if (pepProMap.isEmpty()) {
            log.warn("{}: deep peptide-protein remap found no peptide mappings under {}", name,
                    hasPrebuiltMap ? "(prebuilt map)" : request.roundRoot.getAbsolutePath());
            int proteins = countDistinctProteins(precursors);
            return new RemapResult(precursors, 0, proteins, proteins);
        }

        int proteinsBefore = countDistinctProteins(precursors);
        int mappedPrecursors = 0;
        for (DiaNNPrecursor precursor : precursors) {
            if (precursor == null) {
                continue;
            }
            String sequence = precursor.getSeqString();
            if (sequence == null || sequence.isEmpty()) {
                continue;
            }
            HashSet<String> proteins = pepProMap.get(sequence);
            if (proteins == null || proteins.isEmpty()) {
                continue;
            }
            precursor.setProGroups(proteins.toArray(new String[0]));
            mappedPrecursors++;
        }
        int proteinsAfter = countDistinctProteins(precursors);
        log.info("{}: deep peptide-protein remap applied; mapping peptides={}, mapped precursors={}, proteins before={}, proteins after={}",
                name, pepProMap.size(), mappedPrecursors, proteinsBefore, proteinsAfter);
        return new RemapResult(precursors, pepProMap.size(), proteinsBefore, proteinsAfter);
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
            for (File perFileDir : perFileDirs) {
                int lastRound = detectLastRound(perFileDir, request);
                if (lastRound <= 0) {
                    continue;
                }
                mergeInto(pepProMap, buildPeptideProteinMapFromRounds(perFileDir, lastRound, request, log, name));
            }
            log.info("{}: DEEP_PER_FILE remap collected {} peptides from all per-file round parquets",
                    name, pepProMap.size());
            return pepProMap;
        }

        int lastRound = request.lastRoundHint > 0 ? request.lastRoundHint : detectLastRound(request.roundRoot, request);
        if (lastRound <= 0) {
            return pepProMap;
        }
        pepProMap = buildPeptideProteinMapFromRounds(request.roundRoot, lastRound, request, log, name);
        log.info("{}: deep remap collected {} peptides from round parquets", name, pepProMap.size());
        return pepProMap;
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
                try {
                    DiaNNParquetReader reader = new DiaNNParquetReader(roundParquet.getAbsolutePath(), request.cut);
                    precursors = reader.getDiaNNPrecursors(0.01);
                } catch (Exception e) {
                    logger.warn("{}: failed to read round parquet {}", taskName, roundParquet.getAbsolutePath(), e);
                }
            } else if (request.diaFiles != null && request.diaFiles.length > 0) {
                precursors = DeepIterativeCommonSupport.readPrecursorsFromPerFileParquets(
                        roundDir, 0.01, request.diaFiles, request.cut, logger, taskName);
            }
            addPrecursors(pepProMap, precursors);
        }
        return pepProMap;
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
