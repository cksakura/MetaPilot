package uk.ac.quadram.metapilot.core.runner;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import org.junit.jupiter.api.Test;

class RepresentativePeptideClustererTest {

    @Test
    void clustersThreeTechnicalReplicateGroupsWithThresholds() {
        List<RepresentativePeptideClusterer.Input> inputs = new ArrayList<>();
        inputs.add(input(0, "A1", range(0, 100), "g1"));
        inputs.add(input(1, "A2", union(range(0, 80), range(100, 120)), "g2"));
        inputs.add(input(2, "A3", union(range(0, 79), range(120, 141)), "g3"));
        inputs.add(input(3, "B1", range(1000, 1100), "g4"));
        inputs.add(input(4, "B2", union(range(1000, 1080), range(1100, 1120)), "g5"));
        inputs.add(input(5, "B3", union(range(1000, 1079), range(1120, 1141)), "g6"));
        inputs.add(input(6, "C1", range(2000, 2100), "g7"));
        inputs.add(input(7, "C2", union(range(2000, 2080), range(2100, 2120)), "g8"));
        inputs.add(input(8, "C3", union(range(2000, 2079), range(2120, 2141)), "g9"));

        RepresentativePeptideClusterer.Result result = RepresentativePeptideClusterer.cluster(inputs, 10, 0.60, 0.78);

        assertEquals(3, result.groups.size());
        assertTrue(result.groups.stream().allMatch(group -> group.size() == 3));
    }

    @Test
    void completeLinkPreventsChainMerge() {
        List<RepresentativePeptideClusterer.Input> inputs = List.of(
                input(0, "A", range(0, 100), "g1"),
                input(1, "B", range(20, 120), "g1"),
                input(2, "C", range(40, 140), "g1"));

        RepresentativePeptideClusterer.Result result = RepresentativePeptideClusterer.cluster(inputs, 10, 0.60, 0.80);

        assertEquals(2, result.groups.size());
        assertTrue(result.groups.stream().anyMatch(group -> group.size() == 2));
        assertTrue(result.groups.stream().anyMatch(group -> group.size() == 1));
    }

    @Test
    void lowEvidenceFileStaysSingleton() {
        List<RepresentativePeptideClusterer.Input> inputs = List.of(
                input(0, "A1", range(0, 100), "g1"),
                input(1, "A2", union(range(0, 80), range(100, 120)), "g1"),
                input(2, "Low", range(0, 5), "g1"));

        RepresentativePeptideClusterer.Result result = RepresentativePeptideClusterer.cluster(inputs, 10, 0.60, 0.78);

        assertEquals(2, result.groups.size());
        assertTrue(result.assignments.stream()
                .anyMatch(a -> a.file.rawName.equals("Low") && a.reason.equals("singleton_low_evidence")));
    }

    @Test
    void representativeTieBreakUsesHighestMarkerPeptides() {
        List<RepresentativePeptideClusterer.Input> inputs = List.of(
                input(0, "Large", range(0, 100), "g1"),
                input(1, "Small", range(0, 80), "g1"));

        RepresentativePeptideClusterer.Result result = RepresentativePeptideClusterer.cluster(inputs, 10, 0.60, 0.78);

        assertEquals(1, result.groups.size());
        assertEquals("Large", result.groups.get(0).representative.rawName);
    }

    @Test
    void genomeOverlapDoesNotControlMembership() {
        List<RepresentativePeptideClusterer.Input> inputs = List.of(
                input(0, "A1", range(0, 100), "genome_a"),
                input(1, "A2", union(range(0, 80), range(100, 120)), "genome_b"));

        RepresentativePeptideClusterer.Result result = RepresentativePeptideClusterer.cluster(inputs, 10, 0.60, 0.78);

        assertEquals(1, result.groups.size());
        assertEquals(2, result.groups.get(0).size());
    }

    private static RepresentativePeptideClusterer.Input input(int index, String name, Set<String> peptides,
                                                              String genome) {
        return new RepresentativePeptideClusterer.Input(index, name + ".mzML", name, peptides,
                Collections.singletonList(genome));
    }

    private static Set<String> range(int startInclusive, int endExclusive) {
        return IntStream.range(startInclusive, endExclusive)
                .mapToObj(i -> "PEP" + i)
                .collect(Collectors.toSet());
    }

    private static Set<String> union(Set<String> a, Set<String> b) {
        Set<String> result = a.stream().collect(Collectors.toSet());
        result.addAll(b);
        return result;
    }
}
