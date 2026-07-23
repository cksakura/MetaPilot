package uk.ac.quadram.metapilot.core.runner;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;

final class RepresentativePeptideClusterer {

    static final int EXACT_PAIRWISE_LIMIT = 300;
    static final int MINHASH_HASHES = 128;
    static final int MINHASH_BANDS = 32;
    static final int MINHASH_ROWS_PER_BAND = 4;

    private static final String SOURCE_EXACT_PAIRWISE = "exact_pairwise";
    private static final String SOURCE_LSH_CANDIDATE = "lsh_candidate";
    private static final String SOURCE_COMPLETE_LINK_ON_DEMAND = "complete_link_on_demand";

    private RepresentativePeptideClusterer() {
    }

    static final class Input {
        final int index;
        final String raw;
        final String rawName;
        final LinkedHashSet<String> peptides;
        final List<String> topGenomes;

        Input(int index, String raw, String rawName, Collection<String> peptides, Collection<String> topGenomes) {
            this.index = index;
            this.raw = raw;
            this.rawName = rawName;
            this.peptides = new LinkedHashSet<>();
            if (peptides != null) {
                this.peptides.addAll(peptides);
            }
            this.topGenomes = topGenomes == null ? Collections.emptyList() : new ArrayList<>(topGenomes);
        }

        int peptideCount() {
            return peptides.size();
        }
    }

    static final class PairSimilarity {
        final Input a;
        final Input b;
        final int intersectionPeptides;
        final int unionPeptides;
        final int smallerPeptideSet;
        final int largerPeptideSet;
        final double jaccard;
        final double containment;
        final double genomeTopNOverlap;
        final String source;

        PairSimilarity(Input a, Input b, int intersectionPeptides, String source) {
            this.a = a;
            this.b = b;
            this.intersectionPeptides = intersectionPeptides;
            this.unionPeptides = a.peptideCount() + b.peptideCount() - intersectionPeptides;
            this.smallerPeptideSet = Math.min(a.peptideCount(), b.peptideCount());
            this.largerPeptideSet = Math.max(a.peptideCount(), b.peptideCount());
            this.jaccard = unionPeptides <= 0 ? 0.0 : (double) intersectionPeptides / (double) unionPeptides;
            this.containment = smallerPeptideSet <= 0 ? 0.0 : (double) intersectionPeptides / (double) smallerPeptideSet;
            this.genomeTopNOverlap = overlap(a.topGenomes, b.topGenomes);
            this.source = source;
        }

        boolean passes(double minJaccard, double minContainment) {
            return jaccard >= minJaccard && containment >= minContainment;
        }

        Input other(Input input) {
            if (input == a) return b;
            if (input == b) return a;
            return null;
        }
    }

    static final class Group {
        String groupId;
        final List<Input> members;
        Input representative;
        double representativeMeanDistance;

        Group(List<Input> members) {
            this.members = new ArrayList<>(members);
        }

        int size() {
            return members.size();
        }
    }

    static final class Assignment {
        final Input file;
        final Group group;
        final double peptideJaccardToRepresentative;
        final double peptideContainmentToRepresentative;
        final double genomeOverlapToRepresentative;
        final Input bestMatch;
        final double bestPeptideJaccard;
        final double bestPeptideContainment;
        final double bestGenomeOverlap;
        final String reason;

        Assignment(Input file, Group group, double peptideJaccardToRepresentative,
                   double peptideContainmentToRepresentative, double genomeOverlapToRepresentative,
                   Input bestMatch, double bestPeptideJaccard, double bestPeptideContainment,
                   double bestGenomeOverlap, String reason) {
            this.file = file;
            this.group = group;
            this.peptideJaccardToRepresentative = peptideJaccardToRepresentative;
            this.peptideContainmentToRepresentative = peptideContainmentToRepresentative;
            this.genomeOverlapToRepresentative = genomeOverlapToRepresentative;
            this.bestMatch = bestMatch;
            this.bestPeptideJaccard = bestPeptideJaccard;
            this.bestPeptideContainment = bestPeptideContainment;
            this.bestGenomeOverlap = bestGenomeOverlap;
            this.reason = reason;
        }
    }

    static final class Result {
        final List<Input> inputs;
        final List<Group> groups;
        final List<Assignment> assignments;
        final List<PairSimilarity> pairSimilarities;
        final List<String> representativeFiles;
        final String algorithmMode;
        final int candidatePairCount;
        final int singletonCount;

        Result(List<Input> inputs, List<Group> groups, List<Assignment> assignments,
               List<PairSimilarity> pairSimilarities, List<String> representativeFiles,
               String algorithmMode, int candidatePairCount, int singletonCount) {
            this.inputs = inputs;
            this.groups = groups;
            this.assignments = assignments;
            this.pairSimilarities = pairSimilarities;
            this.representativeFiles = representativeFiles;
            this.algorithmMode = algorithmMode;
            this.candidatePairCount = candidatePairCount;
            this.singletonCount = singletonCount;
        }
    }

    static Result cluster(List<Input> inputs, int minMarkerPeptides, double minJaccard, double minContainment) {
        List<Input> orderedInputs = new ArrayList<>();
        if (inputs != null) {
            for (Input input : inputs) {
                orderedInputs.add(new Input(orderedInputs.size(), input.raw, input.rawName,
                        input.peptides, input.topGenomes));
            }
        }
        int minPeptides = Math.max(0, minMarkerPeptides);
        String algorithmMode = orderedInputs.size() <= EXACT_PAIRWISE_LIMIT ? "exact_pairwise" : "minhash_lsh";

        Map<PairKey, PairSimilarity> pairMap = new HashMap<>();
        Set<PairKey> candidatePairs = "exact_pairwise".equals(algorithmMode)
                ? allPairs(orderedInputs)
                : minHashCandidatePairs(orderedInputs, minPeptides);
        String initialSource = "exact_pairwise".equals(algorithmMode) ? SOURCE_EXACT_PAIRWISE : SOURCE_LSH_CANDIDATE;
        for (PairKey key : candidatePairs) {
            getOrComputePair(orderedInputs.get(key.left), orderedInputs.get(key.right), pairMap, initialSource);
        }

        Set<Integer> eligible = new HashSet<>();
        for (Input input : orderedInputs) {
            if (input.peptideCount() >= minPeptides) {
                eligible.add(input.index);
            }
        }

        Map<Integer, Set<Integer>> passGraph = new HashMap<>();
        for (PairSimilarity sim : pairMap.values()) {
            if (eligible.contains(sim.a.index) && eligible.contains(sim.b.index)
                    && sim.passes(minJaccard, minContainment)) {
                passGraph.computeIfAbsent(sim.a.index, k -> new HashSet<>()).add(sim.b.index);
                passGraph.computeIfAbsent(sim.b.index, k -> new HashSet<>()).add(sim.a.index);
            }
        }

        List<Group> groups = new ArrayList<>();
        Set<Integer> seenEligible = new HashSet<>();
        for (Input input : orderedInputs) {
            if (!eligible.contains(input.index) || seenEligible.contains(input.index)) {
                continue;
            }
            List<Input> component = connectedComponent(input, orderedInputs, eligible, passGraph, seenEligible);
            groups.addAll(completeLinkClusters(component, pairMap, minJaccard, minContainment));
        }

        for (Input input : orderedInputs) {
            if (!eligible.contains(input.index)) {
                groups.add(new Group(Collections.singletonList(input)));
            }
        }

        for (Group group : groups) {
            chooseRepresentative(group, pairMap);
        }

        groups.sort((a, b) -> {
            int cmp = Integer.compare(b.representative.peptideCount(), a.representative.peptideCount());
            if (cmp != 0) return cmp;
            cmp = a.representative.rawName.compareToIgnoreCase(b.representative.rawName);
            if (cmp != 0) return cmp;
            return a.representative.raw.compareToIgnoreCase(b.representative.raw);
        });
        for (int i = 0; i < groups.size(); i++) {
            groups.get(i).groupId = String.format(Locale.ROOT, "G%03d", i + 1);
        }

        List<Assignment> assignments = new ArrayList<>();
        List<String> representativeFiles = new ArrayList<>();
        int singletonCount = 0;
        for (Group group : groups) {
            if (group.size() == 1) {
                singletonCount++;
            }
            representativeFiles.add(group.representative.raw);
            List<Input> members = new ArrayList<>(group.members);
            members.sort((a, b) -> {
                if (a == group.representative) return -1;
                if (b == group.representative) return 1;
                int cmp = Integer.compare(b.peptideCount(), a.peptideCount());
                if (cmp != 0) return cmp;
                return a.rawName.compareToIgnoreCase(b.rawName);
            });
            for (Input member : members) {
                assignments.add(createAssignment(member, group, pairMap, minPeptides));
            }
        }

        List<PairSimilarity> pairSimilarities = new ArrayList<>(pairMap.values());
        pairSimilarities.sort(Comparator
                .comparingInt((PairSimilarity s) -> s.a.index)
                .thenComparingInt(s -> s.b.index));

        return new Result(orderedInputs, groups, assignments, pairSimilarities, representativeFiles,
                algorithmMode, candidatePairs.size(), singletonCount);
    }

    private static Assignment createAssignment(Input member, Group group, Map<PairKey, PairSimilarity> pairMap,
                                               int minPeptides) {
        double jaccardToRep = 1.0;
        double containmentToRep = 1.0;
        double genomeOverlapToRep = 1.0;
        if (member != group.representative) {
            PairSimilarity toRep = getOrComputePair(member, group.representative, pairMap,
                    SOURCE_COMPLETE_LINK_ON_DEMAND);
            jaccardToRep = toRep.jaccard;
            containmentToRep = toRep.containment;
            genomeOverlapToRep = toRep.genomeTopNOverlap;
        }

        PairSimilarity best = bestComparedMatch(member, pairMap);
        String reason;
        if (member == group.representative) {
            if (member.peptideCount() < minPeptides) {
                reason = "singleton_low_evidence";
            } else if (group.size() == 1) {
                reason = "singleton_no_match";
            } else {
                reason = "representative_medoid";
            }
        } else {
            reason = "member_complete_link";
        }

        return new Assignment(member, group, jaccardToRep, containmentToRep, genomeOverlapToRep,
                best == null ? null : best.other(member),
                best == null ? 0.0 : best.jaccard,
                best == null ? 0.0 : best.containment,
                best == null ? 0.0 : best.genomeTopNOverlap,
                reason);
    }

    private static PairSimilarity bestComparedMatch(Input input, Map<PairKey, PairSimilarity> pairMap) {
        PairSimilarity best = null;
        double bestScore = -1.0;
        for (PairSimilarity sim : pairMap.values()) {
            if (sim.a != input && sim.b != input) {
                continue;
            }
            double score = sim.jaccard + sim.containment;
            if (best == null || score > bestScore
                    || (score == bestScore && sim.other(input).rawName.compareToIgnoreCase(best.other(input).rawName) < 0)) {
                best = sim;
                bestScore = score;
            }
        }
        return best;
    }

    private static List<Group> completeLinkClusters(List<Input> component, Map<PairKey, PairSimilarity> pairMap,
                                                    double minJaccard, double minContainment) {
        if (component.size() <= 1) {
            return Collections.singletonList(new Group(component));
        }

        List<Group> clusters = new ArrayList<>();
        for (Input input : component) {
            clusters.add(new Group(Collections.singletonList(input)));
        }

        while (true) {
            MergeCandidate best = null;
            for (int i = 0; i < clusters.size(); i++) {
                for (int j = i + 1; j < clusters.size(); j++) {
                    MergeCandidate candidate = evaluateCompleteLinkMerge(clusters.get(i), clusters.get(j),
                            pairMap, minJaccard, minContainment);
                    if (candidate == null) {
                        continue;
                    }
                    if (best == null || candidate.score > best.score
                            || (candidate.score == best.score && candidate.tieBreakKey().compareTo(best.tieBreakKey()) < 0)) {
                        candidate.leftIndex = i;
                        candidate.rightIndex = j;
                        best = candidate;
                    }
                }
            }
            if (best == null) {
                break;
            }
            Group merged = new Group(best.mergedMembers());
            clusters.remove(best.rightIndex);
            clusters.remove(best.leftIndex);
            clusters.add(merged);
        }
        return clusters;
    }

    private static MergeCandidate evaluateCompleteLinkMerge(Group left, Group right,
                                                            Map<PairKey, PairSimilarity> pairMap,
                                                            double minJaccard, double minContainment) {
        double sumJaccard = 0.0;
        double sumContainment = 0.0;
        int count = 0;
        for (Input a : left.members) {
            for (Input b : right.members) {
                PairSimilarity sim = getOrComputePair(a, b, pairMap, SOURCE_COMPLETE_LINK_ON_DEMAND);
                if (!sim.passes(minJaccard, minContainment)) {
                    return null;
                }
                sumJaccard += sim.jaccard;
                sumContainment += sim.containment;
                count++;
            }
        }
        double score = count == 0 ? 0.0 : (sumJaccard + sumContainment) / count;
        return new MergeCandidate(left, right, score);
    }

    private static void chooseRepresentative(Group group, Map<PairKey, PairSimilarity> pairMap) {
        if (group.members.isEmpty()) {
            return;
        }
        Input best = null;
        double bestMeanDistance = Double.POSITIVE_INFINITY;
        for (Input candidate : group.members) {
            double meanDistance = meanDistance(candidate, group.members, pairMap);
            if (best == null || meanDistance < bestMeanDistance
                    || (meanDistance == bestMeanDistance && candidate.peptideCount() > best.peptideCount())
                    || (meanDistance == bestMeanDistance && candidate.peptideCount() == best.peptideCount()
                    && candidate.rawName.compareToIgnoreCase(best.rawName) < 0)) {
                best = candidate;
                bestMeanDistance = meanDistance;
            }
        }
        group.representative = best;
        group.representativeMeanDistance = Double.isFinite(bestMeanDistance) ? bestMeanDistance : 0.0;
    }

    private static double meanDistance(Input candidate, List<Input> members, Map<PairKey, PairSimilarity> pairMap) {
        if (members.size() <= 1) {
            return 0.0;
        }
        double sum = 0.0;
        int count = 0;
        for (Input other : members) {
            if (candidate == other) {
                continue;
            }
            PairSimilarity sim = getOrComputePair(candidate, other, pairMap, SOURCE_COMPLETE_LINK_ON_DEMAND);
            sum += 1.0 - sim.jaccard;
            count++;
        }
        return count == 0 ? 0.0 : sum / count;
    }

    private static List<Input> connectedComponent(Input start, List<Input> inputs, Set<Integer> eligible,
                                                  Map<Integer, Set<Integer>> passGraph, Set<Integer> seen) {
        List<Input> component = new ArrayList<>();
        ArrayDeque<Integer> queue = new ArrayDeque<>();
        queue.add(start.index);
        seen.add(start.index);
        while (!queue.isEmpty()) {
            int idx = queue.removeFirst();
            component.add(inputs.get(idx));
            for (Integer next : passGraph.getOrDefault(idx, Collections.emptySet())) {
                if (eligible.contains(next) && seen.add(next)) {
                    queue.addLast(next);
                }
            }
        }
        return component;
    }

    private static Set<PairKey> allPairs(List<Input> inputs) {
        Set<PairKey> pairs = new HashSet<>();
        for (int i = 0; i < inputs.size(); i++) {
            for (int j = i + 1; j < inputs.size(); j++) {
                pairs.add(new PairKey(i, j));
            }
        }
        return pairs;
    }

    private static Set<PairKey> minHashCandidatePairs(List<Input> inputs, int minPeptides) {
        Set<PairKey> pairs = new HashSet<>();
        if (inputs.isEmpty()) {
            return pairs;
        }
        long[][] signatures = new long[inputs.size()][MINHASH_HASHES];
        for (long[] row : signatures) {
            Arrays.fill(row, -1L);
        }
        for (Input input : inputs) {
            if (input.peptideCount() < minPeptides) {
                continue;
            }
            for (String peptide : input.peptides) {
                for (int i = 0; i < MINHASH_HASHES; i++) {
                    long hash = hash64(peptide, i);
                    if (Long.compareUnsigned(hash, signatures[input.index][i]) < 0) {
                        signatures[input.index][i] = hash;
                    }
                }
            }
        }

        Map<String, List<Integer>> buckets = new HashMap<>();
        for (Input input : inputs) {
            if (input.peptideCount() < minPeptides) {
                continue;
            }
            for (int band = 0; band < MINHASH_BANDS; band++) {
                StringBuilder key = new StringBuilder();
                key.append(band);
                int offset = band * MINHASH_ROWS_PER_BAND;
                for (int row = 0; row < MINHASH_ROWS_PER_BAND; row++) {
                    key.append(':').append(Long.toUnsignedString(signatures[input.index][offset + row]));
                }
                buckets.computeIfAbsent(key.toString(), k -> new ArrayList<>()).add(input.index);
            }
        }

        for (List<Integer> bucket : buckets.values()) {
            if (bucket.size() < 2) {
                continue;
            }
            for (int i = 0; i < bucket.size(); i++) {
                for (int j = i + 1; j < bucket.size(); j++) {
                    pairs.add(new PairKey(bucket.get(i), bucket.get(j)));
                }
            }
        }
        return pairs;
    }

    private static PairSimilarity getOrComputePair(Input a, Input b, Map<PairKey, PairSimilarity> pairMap,
                                                   String source) {
        PairKey key = new PairKey(a.index, b.index);
        PairSimilarity existing = pairMap.get(key);
        if (existing != null) {
            return existing;
        }
        PairSimilarity sim = new PairSimilarity(key.left == a.index ? a : b,
                key.left == a.index ? b : a, intersectionSize(a.peptides, b.peptides), source);
        pairMap.put(key, sim);
        return sim;
    }

    private static int intersectionSize(Set<String> a, Set<String> b) {
        Set<String> small = a.size() <= b.size() ? a : b;
        Set<String> large = a.size() <= b.size() ? b : a;
        int count = 0;
        for (String value : small) {
            if (large.contains(value)) {
                count++;
            }
        }
        return count;
    }

    private static double overlap(List<String> a, List<String> b) {
        if (a == null || b == null || a.isEmpty() || b.isEmpty()) {
            return 0.0;
        }
        HashSet<String> bs = new HashSet<>(b);
        int count = 0;
        for (String value : a) {
            if (bs.contains(value)) {
                count++;
            }
        }
        int denominator = Math.min(a.size(), b.size());
        return denominator <= 0 ? 0.0 : (double) count / (double) denominator;
    }

    private static long hash64(String value, int seed) {
        long hash = 0xcbf29ce484222325L ^ (0x9e3779b97f4a7c15L * (seed + 1L));
        for (int i = 0; i < value.length(); i++) {
            hash ^= value.charAt(i);
            hash *= 0x100000001b3L;
        }
        hash ^= (hash >>> 33);
        hash *= 0xff51afd7ed558ccdL;
        hash ^= (hash >>> 33);
        hash *= 0xc4ceb9fe1a85ec53L;
        hash ^= (hash >>> 33);
        return hash;
    }

    private static final class PairKey {
        final int left;
        final int right;

        PairKey(int a, int b) {
            this.left = Math.min(a, b);
            this.right = Math.max(a, b);
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (!(o instanceof PairKey)) return false;
            PairKey pairKey = (PairKey) o;
            return left == pairKey.left && right == pairKey.right;
        }

        @Override
        public int hashCode() {
            return 31 * left + right;
        }
    }

    private static final class MergeCandidate {
        int leftIndex;
        int rightIndex;
        final Group left;
        final Group right;
        final double score;

        MergeCandidate(Group left, Group right, double score) {
            this.left = left;
            this.right = right;
            this.score = score;
        }

        List<Input> mergedMembers() {
            List<Input> merged = new ArrayList<>(left.members);
            merged.addAll(right.members);
            return merged;
        }

        String tieBreakKey() {
            List<String> names = new ArrayList<>();
            for (Input input : mergedMembers()) {
                names.add(input.rawName);
            }
            Collections.sort(names, String.CASE_INSENSITIVE_ORDER);
            return String.join("\t", names);
        }
    }
}
