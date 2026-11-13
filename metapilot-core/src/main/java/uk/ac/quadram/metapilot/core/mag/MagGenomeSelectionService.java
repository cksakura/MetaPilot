package uk.ac.quadram.metapilot.core.mag;

import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;
import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNPrecursor;
import uk.ac.quadram.metapilot.task.io.pep.MetaPSM;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;
import uk.ac.quadram.metapilot.task.mag.MarkerGenomeSelector;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Locale;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Shared MAG genome selection from PSM evidence.
 *
 * Search-engine-specific code should only parse native result files into
 * {@link MetaPSM} implementations, then delegate marker selection here.
 */
public final class MagGenomeSelectionService {

    public static final int DEFAULT_MAXIMUM_SELECTED = 300;
    public static final int UNLIMITED_MAXIMUM_SELECTED = MarkerGenomeSelector.UNLIMITED_MAXIMUM_SELECTED;
    public static final String SELECTION_STRATEGY_FOCUS = MarkerGenomeSelector.KNEE_MODE_FOCUS;
    public static final String SELECTION_STRATEGY_COMPREHENSIVE = MarkerGenomeSelector.KNEE_MODE_COMPREHENSIVE;
    public static final String SELECTION_STRATEGY_PSM_COVERAGE = MarkerGenomeSelector.KNEE_MODE_PSM_COVERAGE;

    private MagGenomeSelectionService() {
    }

    public static final class Request {
        public MagDbItem magDb;
        public String genomeSelectorMode;
        public String genomeSelectionStrategy;
        public double genomeSelectionPsmCoverage = uk.ac.quadram.metapilot.task.mag.GenomeSelectionCut.DEFAULT_PSM_COVERAGE;
        public double genomeSelectionStepGuard = uk.ac.quadram.metapilot.task.mag.GenomeSelectionCut.DEFAULT_STEP_GUARD;
        public int genomeSelectionMaxGenomes = uk.ac.quadram.metapilot.task.mag.GenomeSelectionCut.DEFAULT_MAX_GENOMES;
        public Collection<? extends MetaPSM> psms;
        public DiaNNPrecursor[] diaPrecursors;
        public File inputResultFile;
        public File markerFasta;
        public File genomesTaxaTsv;
        public boolean disableMarkerFasta;
        public File evidenceFile;
        public File rankedFile;
        public File selectedFile;
        public int maximumSelected = DEFAULT_MAXIMUM_SELECTED;
        public String representativeLevel;
        public boolean inferFromProteinNamesWhenEmpty;
        public boolean allowMarker4LibraryFallback;
        public boolean logDiagnostics;
        public String sourceDescription;
        /**
         * Files the selection is derived from, for the staleness check below.
         * <p>
         * Optional: leave null and the resume check is existence-only, as it was.
         * Supply the PSM file a caller just rewrote and a cached selection older
         * than it is recomputed instead of reused - which is what a caller means
         * when it rewrites its evidence and calls again.
         */
        public File[] inputs;
    }

    public static HashSet<String> selectGenomes(Request request, ProgressCallback cb) throws Exception {
        if (cb == null) cb = new ProgressCallback.NoOp();
        validate(request);

        if (hasContent(request.evidenceFile) && hasContent(request.rankedFile) && hasContent(request.selectedFile)) {
            // Existence is not enough. This resume check sits INSIDE a step whose
            // own caller may already have decided to re-run: the open-search
            // workflow rewrites unmodified_psms.tsv with a rescored PSM set, calls
            // here, and - before this check looked at timestamps - got the previous
            // run's genomes back and rebuilt a byte-identical FASTA from them. The
            // rescoring was computed, written to disk, and then discarded, with the
            // selection audit files left describing evidence that no longer
            // existed. A caller that supplies request.inputs gets the selection
            // recomputed whenever its evidence is newer than the cached answer.
            File stale = newerThan(request.selectedFile, request.inputs);
            if (stale == null) {
                cb.onMessage("Genome selection already exists, skipping: " + request.selectedFile.getAbsolutePath());
                return readSelectedGenomes(request.selectedFile);
            }
            cb.onMessage("Genome selection is older than " + stale.getName() + "; recomputing it from the current "
                    + "evidence rather than reusing " + request.selectedFile.getName());
        }

        List<MetaPSM> psmList = new ArrayList<>();
        if (request.psms != null) {
            for (MetaPSM psm : request.psms) {
                if (psm != null) {
                    psmList.add(psm);
                }
            }
        }

        HashSet<String> genomeSet = new HashSet<>();
        boolean useMarker4 = MagDbItem.MARKER4.equalsIgnoreCase(request.genomeSelectorMode);
        SelectionResources resources = resolveResources(request, useMarker4, cb);
        boolean enableAutoRep = useMarker4;
        MarkerGenomeSelector selectorUsed = null;
        if (!psmList.isEmpty()) {
            MarkerGenomeSelector selector = newSelector(request);
            selectorUsed = selector;
            MetaPSM[] psmArray = psmList.toArray(new MetaPSM[0]);
            if (resources.markerFasta != null && request.logDiagnostics) {
                logMarkerSelectionDiagnostics(psmList, resources.markerFasta, cb);
            }
            selector.selectMarkerGenomesWithKnee(psmArray, resources.markerFasta, resources.genomesTaxaTsv,
                    request.evidenceFile.getAbsolutePath(),
                    request.rankedFile.getAbsolutePath(),
                    request.selectedFile.getAbsolutePath(),
                    request.maximumSelected,
                    enableAutoRep);
        } else if (request.diaPrecursors != null && request.diaPrecursors.length > 0) {
            MarkerGenomeSelector selector = newSelector(request);
            selectorUsed = selector;
            selector.selectMarkerGenomesWithKnee(request.diaPrecursors, resources.markerFasta, resources.genomesTaxaTsv,
                    request.evidenceFile.getAbsolutePath(),
                    request.rankedFile.getAbsolutePath(),
                    request.selectedFile.getAbsolutePath(),
                    request.maximumSelected,
                    enableAutoRep);
        } else if (hasContent(request.inputResultFile)) {
            MarkerGenomeSelector selector = newSelector(request);
            selectorUsed = selector;
            selector.selectMarkerGenomesWithKnee(request.inputResultFile.getAbsolutePath(),
                    resources.markerFasta, resources.genomesTaxaTsv,
                    request.evidenceFile.getAbsolutePath(),
                    request.rankedFile.getAbsolutePath(),
                    request.selectedFile.getAbsolutePath(),
                    request.maximumSelected,
                    enableAutoRep);
        }

        if (selectorUsed != null) {
            publishCut(selectorUsed, request, cb);
        }

        genomeSet.addAll(readSelectedGenomes(request.selectedFile));
        if (genomeSet.isEmpty() && request.inferFromProteinNamesWhenEmpty && !psmList.isEmpty()) {
            genomeSet.addAll(inferGenomesFromProteinNames(psmList, request.magDb));
            if (!genomeSet.isEmpty()) {
                cb.onMessage("WARNING: marker selector produced no genomes; falling back to " + genomeSet.size()
                        + " genomes inferred from proteins");
            }
        }

        if (genomeSet.isEmpty()) {
            String source = request.sourceDescription == null ? String.valueOf(request.selectedFile) : request.sourceDescription;
            cb.onMessage("WARNING: no genomes selected from " + source);
        }
        return genomeSet;
    }

    public static HashSet<String> readSelectedGenomes(File selectedFile) throws IOException {
        HashSet<String> genomeSet = new HashSet<>();
        if (!hasContent(selectedFile)) {
            return genomeSet;
        }
        try (BufferedReader reader = new BufferedReader(new FileReader(selectedFile))) {
            reader.readLine(); // header
            String line;
            while ((line = reader.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty() || line.startsWith("#")) {
                    break;
                }
                String[] parts = line.split("\t");
                if (parts.length > 0 && !parts[0].trim().isEmpty()) {
                    genomeSet.add(parts[0].trim());
                }
            }
        }
        return genomeSet;
    }

    private static void validate(Request request) {
        if (request == null) {
            throw new IllegalArgumentException("request is null");
        }
        if (request.magDb == null) {
            throw new IllegalArgumentException("magDb is null");
        }
        if (request.evidenceFile == null || request.rankedFile == null || request.selectedFile == null) {
            throw new IllegalArgumentException("selection output files are required");
        }
        if (request.maximumSelected <= 0) {
            request.maximumSelected = DEFAULT_MAXIMUM_SELECTED;
        }
    }

    private static void publishCut(MarkerGenomeSelector selector, Request request, ProgressCallback cb) {
        MarkerGenomeSelector.SelectionResult result = selector.getLastSelectionResult();
        if (result == null || request.selectedFile == null) {
            return;
        }
        String raw = rawFileName(request);
        int genomes = result.selected == null ? 0 : result.selected.size();
        String coverage = Double.isNaN(result.psmCoverageReached) ? "n/a"
                : String.format(Locale.ROOT, "%.1f%%", result.psmCoverageReached * 100.0);
        cb.onMessage(raw + ": " + result.stopReason + ", PSM coverage " + coverage + ", " + genomes + " genomes");
        if (uk.ac.quadram.metapilot.task.mag.GenomeSelectionCut.isPsmCoverage(result.selectionStrategy)) {
            double target = uk.ac.quadram.metapilot.task.mag.GenomeSelectionCut.DEFAULT_PSM_COVERAGE;
            try {
                if (result.psmCoverageThresholdText != null && !result.psmCoverageThresholdText.isEmpty()) {
                    target = Double.parseDouble(result.psmCoverageThresholdText);
                }
            } catch (NumberFormatException ignored) {
            }
            int maxGenomes = request.genomeSelectionMaxGenomes;
            String missed = uk.ac.quadram.metapilot.task.mag.GenomeSelectionCut.missedTargetMessage(target,
                    result.psmCoverageReached, result.stopReason, maxGenomes);
            if (!missed.isEmpty()) {
                cb.onMessage(raw + ": " + missed);
            }
        }
        try {
            uk.ac.quadram.metapilot.task.mag.GenomeSelectionSummary.writeSample(
                    new File(request.selectedFile.getParentFile(),
                            uk.ac.quadram.metapilot.task.mag.GenomeSelectionSummary.FILE_NAME),
                    raw, result);
        } catch (IOException e) {
            cb.onMessage("WARNING: could not write genome selection summary for " + raw + ": " + e.getMessage());
        }
    }

    private static String rawFileName(Request request) {
        File selected = request.selectedFile;
        if (selected != null && selected.getParentFile() != null && selected.getParentFile().getParentFile() != null) {
            return selected.getParentFile().getParentFile().getName();
        }
        if (request.sourceDescription != null) {
            return new File(request.sourceDescription).getName();
        }
        return "";
    }

    private static String normalizeSelectionStrategy(String strategy) {
        if (strategy == null || strategy.trim().isEmpty()) {
            return SELECTION_STRATEGY_FOCUS;
        }
        return strategy.trim();
    }

    private static MarkerGenomeSelector newSelector(Request request) {
        MarkerGenomeSelector selector = new MarkerGenomeSelector();
        selector.setKneeSelectionMode(normalizeSelectionStrategy(request.genomeSelectionStrategy));
        selector.setPsmCoverageCut(request.genomeSelectionPsmCoverage, request.genomeSelectionStepGuard,
                request.genomeSelectionMaxGenomes);
        if (request.representativeLevel != null && !request.representativeLevel.trim().isEmpty()) {
            selector.setRepresentativeLevel(request.representativeLevel);
        }
        return selector;
    }

    private static SelectionResources resolveResources(Request request, boolean useMarker4, ProgressCallback cb) throws IOException {
        if (!useMarker4) {
            return new SelectionResources(null, request.genomesTaxaTsv);
        }
        File current = request.magDb.getCurrentFile();
        if (current == null) {
            throw new IOException("Marker4 mode requires a resolved MAG DB current folder");
        }
        File defaultMarkerFasta = new File(current, Marker4ResourcePreparer.MARKER4_FASTA_NAME);
        File defaultGenomesTaxaTsv = new File(current, Marker4ResourcePreparer.GENOMES_TAXA_TSV_NAME);
        File markerFasta = request.disableMarkerFasta ? null
                : request.markerFasta != null ? request.markerFasta : defaultMarkerFasta;
        File genomesTaxaTsv = request.genomesTaxaTsv != null ? request.genomesTaxaTsv : defaultGenomesTaxaTsv;

        boolean usesDefaultMarkerFasta = !request.disableMarkerFasta && samePath(markerFasta, defaultMarkerFasta);
        boolean usesDefaultGenomesTaxa = samePath(genomesTaxaTsv, defaultGenomesTaxaTsv);
        if (usesDefaultMarkerFasta && !hasContent(markerFasta)) {
            Marker4ResourcePreparer.SelectionResources prepared =
                    Marker4ResourcePreparer.ensureSelectionResources(request.magDb, cb);
            markerFasta = prepared.marker4Fasta;
            if (usesDefaultGenomesTaxa) {
                genomesTaxaTsv = prepared.genomesTaxaTsv;
            }
        } else if (usesDefaultGenomesTaxa) {
            genomesTaxaTsv = Marker4ResourcePreparer.ensureGenomesTaxaTsv(request.magDb, defaultGenomesTaxaTsv, cb);
        }

        if (!request.disableMarkerFasta && !hasContent(markerFasta) && request.allowMarker4LibraryFallback) {
            File libraryMarkerFasta = request.magDb.getMarker4FastaForLibrary();
            if (hasContent(libraryMarkerFasta)) {
                markerFasta = libraryMarkerFasta;
            }
        }
        if (!request.disableMarkerFasta && !hasContent(markerFasta)) {
            throw new IOException("Marker4 mode requires Marker4.fasta for " + current.getAbsolutePath());
        }
        if (!hasContent(genomesTaxaTsv)) {
            throw new IOException("Marker4 mode requires genomes-taxa.tsv for " + current.getAbsolutePath());
        }
        return new SelectionResources(markerFasta, genomesTaxaTsv);
    }

    private static boolean samePath(File a, File b) {
        if (a == null || b == null) {
            return false;
        }
        return a.getAbsoluteFile().equals(b.getAbsoluteFile());
    }

    private static final class SelectionResources {
        final File markerFasta;
        final File genomesTaxaTsv;

        SelectionResources(File markerFasta, File genomesTaxaTsv) {
            this.markerFasta = markerFasta;
            this.genomesTaxaTsv = genomesTaxaTsv;
        }
    }

    private static HashSet<String> inferGenomesFromProteinNames(List<? extends MetaPSM> psms, MagDbItem magDb) {
        HashSet<String> genomeSet = new HashSet<>();
        for (MetaPSM psm : psms) {
            if (psm == null) {
                continue;
            }
            String[] proteins = psm.getProteins();
            if (proteins == null) {
                continue;
            }
            for (String protein : proteins) {
                if (protein == null) {
                    continue;
                }
                String pro = protein.trim();
                if (!pro.startsWith(magDb.getIdentifier())) {
                    continue;
                }
                int idx = pro.indexOf("_");
                if (idx > 0) {
                    genomeSet.add(pro.substring(0, idx));
                }
            }
        }
        return genomeSet;
    }

    private static boolean hasContent(File file) {
        return file != null && file.exists() && file.length() > 0;
    }

    /** The first input newer than {@code artefact}, or null when none is. */
    private static File newerThan(File artefact, File[] inputs) {
        if (inputs == null || artefact == null) {
            return null;
        }
        for (File input : inputs) {
            if (input != null && input.exists() && input.lastModified() > artefact.lastModified()) {
                return input;
            }
        }
        return null;
    }

    private static void logMarkerSelectionDiagnostics(List<? extends MetaPSM> psmList, File markerFasta,
                                                      ProgressCallback cb) {
        try {
            Map<String, String> markerCategories = loadMarkerCategories(markerFasta);
            Map<String, Set<String>> peptideToGenomes = new HashMap<>();
            Map<String, Set<String>> genomeToUniquePeptides = new HashMap<>();
            HashSet<String> allPeptides = new HashSet<>();
            HashSet<String> markerPeptides = new HashSet<>();
            HashSet<String> genomesWithMarkerEvidence = new HashSet<>();

            int emptyProteinRows = 0;
            int targetProteinRefs = 0;
            int decoyProteinRefs = 0;
            int markerProteinRefs = 0;

            for (MetaPSM psm : psmList) {
                if (psm == null) {
                    continue;
                }
                String peptide = psm.getSequence();
                if (peptide == null || peptide.isEmpty()) {
                    continue;
                }
                allPeptides.add(peptide);
                String[] proteins = psm.getProteins();
                if (proteins == null || proteins.length == 0) {
                    emptyProteinRows++;
                    continue;
                }
                for (String protein : proteins) {
                    String proteinId = normalizeProteinId(protein);
                    if (proteinId.isEmpty()) {
                        continue;
                    }
                    if (isDecoyProteinId(proteinId)) {
                        decoyProteinRefs++;
                        continue;
                    }
                    targetProteinRefs++;
                    String markerProteinId = resolveMarkerProteinId(proteinId, markerCategories);
                    if (markerProteinId == null) {
                        continue;
                    }
                    markerProteinRefs++;
                    String genomeId = extractGenomeId(markerProteinId);
                    if (genomeId == null || genomeId.isEmpty()) {
                        continue;
                    }
                    markerPeptides.add(peptide);
                    genomesWithMarkerEvidence.add(genomeId);
                    peptideToGenomes.computeIfAbsent(peptide, k -> new HashSet<>()).add(genomeId);
                }
            }

            int uniqueGenomePeptides = 0;
            int sharedGenomePeptides = 0;
            int m2 = 0;
            int m3to5 = 0;
            int m6to20 = 0;
            int mOver20 = 0;
            ArrayList<String> sharedExamples = new ArrayList<>();
            for (Map.Entry<String, Set<String>> entry : peptideToGenomes.entrySet()) {
                int genomeCount = entry.getValue().size();
                if (genomeCount == 1) {
                    uniqueGenomePeptides++;
                    String genome = entry.getValue().iterator().next();
                    genomeToUniquePeptides.computeIfAbsent(genome, k -> new HashSet<>()).add(entry.getKey());
                } else if (genomeCount > 1) {
                    sharedGenomePeptides++;
                    if (genomeCount == 2) {
                        m2++;
                    } else if (genomeCount <= 5) {
                        m3to5++;
                    } else if (genomeCount <= 20) {
                        m6to20++;
                    } else {
                        mOver20++;
                    }
                    if (sharedExamples.size() < 3) {
                        sharedExamples.add(entry.getKey() + "->" + genomeCount + " genomes");
                    }
                }
            }

            int maxUnique = 0;
            for (Set<String> peptides : genomeToUniquePeptides.values()) {
                maxUnique = Math.max(maxUnique, peptides.size());
            }
            cb.onMessage("Marker diagnostics: markerProteins=" + markerCategories.size()
                    + ", uniquePeptides(all)=" + allPeptides.size()
                    + ", markerPeptides=" + markerPeptides.size()
                    + ", genomeUniqueMarkerPeptides=" + uniqueGenomePeptides
                    + ", sharedMarkerPeptides=" + sharedGenomePeptides);
            cb.onMessage("Marker diagnostics: genomesWithMarkerEvidence=" + genomesWithMarkerEvidence.size()
                    + ", genomesWithUniquePeptides=" + genomeToUniquePeptides.size()
                    + ", maxUniquePeptidesPerGenome=" + maxUnique
                    + ", emptyProteinRows=" + emptyProteinRows);
            cb.onMessage("Marker diagnostics: proteinRefs target=" + targetProteinRefs
                    + ", markerMatched=" + markerProteinRefs
                    + ", decoyIgnored=" + decoyProteinRefs
                    + ", sharedMultiplicity m2=" + m2
                    + ", m3-5=" + m3to5
                    + ", m6-20=" + m6to20
                    + ", m>20=" + mOver20);
            if (!sharedExamples.isEmpty()) {
                cb.onMessage("Marker diagnostics: shared peptide examples " + String.join("; ", sharedExamples));
            }
            if (uniqueGenomePeptides == 0 && markerPeptides.size() > 0) {
                cb.onMessage("Marker diagnostics: no genome-unique marker peptides; all matched marker peptides map to multiple genomes.");
            }
        } catch (Exception e) {
            cb.onMessage("Marker diagnostics unavailable: " + e.getMessage());
        }
    }

    private static Map<String, String> loadMarkerCategories(File markerFasta) throws IOException {
        Map<String, String> categories = new HashMap<>();
        try (BufferedReader reader = new BufferedReader(new FileReader(markerFasta))) {
            String line;
            while ((line = reader.readLine()) != null) {
                if (!line.startsWith(">")) {
                    continue;
                }
                String header = line.substring(1).trim();
                if (header.isEmpty()) {
                    continue;
                }
                String proteinId = header.split("\\s+", 2)[0];
                String category = MarkerGenomeSelector.categorizeMarkerHeader(proteinId, line);
                if (category == null) {
                    continue;
                }
                categories.put(proteinId, category);
            }
        }
        return categories;
    }

    private static String resolveMarkerProteinId(String proteinId, Map<String, String> markerCategories) {
        if (markerCategories.containsKey(proteinId)) {
            return proteinId;
        }
        int last = proteinId.lastIndexOf("_");
        if (last > 0) {
            String baseId = proteinId.substring(0, last);
            if (markerCategories.containsKey(baseId)) {
                return baseId;
            }
        }
        return null;
    }

    private static String normalizeProteinId(String protein) {
        if (protein == null) {
            return "";
        }
        String id = protein.trim();
        if (id.startsWith(">")) {
            id = id.substring(1).trim();
        }
        int whitespace = -1;
        for (int i = 0; i < id.length(); i++) {
            if (Character.isWhitespace(id.charAt(i))) {
                whitespace = i;
                break;
            }
        }
        if (whitespace > 0) {
            id = id.substring(0, whitespace);
        }
        return id;
    }

    private static boolean isDecoyProteinId(String proteinId) {
        if (proteinId == null) {
            return false;
        }
        String lower = proteinId.toLowerCase();
        return lower.startsWith("rev_") || lower.startsWith("decoy_");
    }

    private static String extractGenomeId(String proteinId) {
        if (proteinId == null) {
            return "";
        }
        int underscore = proteinId.indexOf('_');
        if (underscore > 0) {
            return proteinId.substring(0, underscore);
        }
        return proteinId;
    }
}
