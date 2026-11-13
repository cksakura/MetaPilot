package uk.ac.quadram.metapilot.core.mag;

import static org.junit.jupiter.api.Assertions.*;
import java.nio.file.*;
import java.util.*;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import org.dom4j.Element;
import uk.ac.quadram.metapilot.task.io.pep.MetaPeptide;

class MagResultFinalisationServiceTest {
    @TempDir Path temp;
    private static MetaPeptide peptide(String seq, double amount, String... proteins) {
        return new MetaPeptide(seq, seq, proteins, new double[]{amount}, new int[]{1}) {
            public Object[] getTableObjects() { return null; }
            public Element getXmlPepElement() { return null; }
        };
    }
    @Test void gateThenRankDropsRedundantGenomeAndReassignsSharedEvidence() throws Exception {
        HashMap<String,double[]> intensities = new HashMap<>();
        intensities.put("CAT1",new double[]{200});
        intensities.put("CAT2",new double[]{100});
        intensities.put("CAT3",new double[]{1000});
        Map<String,double[]> stats = Map.of("CAT1",new double[]{10,.001,.001},
                "CAT2",new double[]{9,.001,.001}, "CAT3",new double[]{1,.2,.3});
        Map<String,Set<String>> support = Map.of("CAT1",Set.of("SHARED"),"CAT2",Set.of("SHARED"),"CAT3",Set.of("BAD"));
        var selection = MagResultFinalisationService.select(intensities,stats,null,Map.of(),support,Map.of(),Map.of());
        assertEquals(Set.of("CAT1"),selection.genomes);
        MetaPeptide shared = peptide("SHARED",100,"CAT2_p2","CAT1_p1","CAT3_p3","shuf_CAT1_p1");
        MetaPeptide rejected = peptide("BAD",1000,"CAT3_p3");
        MetaPeptide host = peptide("HOST",50,"HOST1");
        MetaPeptide[] kept = MagResultFinalisationService.prune(new MetaPeptide[]{shared,rejected,host},"CAT","shuf_",selection);
        assertArrayEquals(new String[]{"CAT1_p1"},kept[0].getProteins());
        assertNull(kept[1]);
        assertArrayEquals(new String[]{"HOST1"},kept[2].getProteins());
        assertEquals(1,selection.excludedPeptides.size());
        var rebuilt = new MagResultFinalisationService.FinalGenomeIntensityMap(Map.of(), selection);
        assertEquals(selection.genomes, rebuilt.keySet());
        assertArrayEquals(new double[]{0}, rebuilt.get("CAT1"));
        selection.write(temp.toFile());
        String decisions = Files.readString(temp.resolve("taxonomy_analysis/genome_final_decisions.tsv"));
        assertTrue(decisions.contains("CAT3\tfalse\tGENOME_CONFIDENCE_GATE"));
        assertTrue(decisions.contains("CAT2\tfalse\tNO_INCREMENTAL_COVERAGE"));
        assertTrue(Files.readString(temp.resolve("excluded_evidence/peptides.tsv")).contains("BAD"));
    }
    @Test void ddaRebuildsAllProteinQuantitiesFromRetainedLinksOnly() {
        MetaPeptide shared = peptide("SHARED",120,"CAT1_a","CAT2_b");
        MetaPeptide unique = peptide("UNIQUE",30,"CAT1_a");
        MetaPeptide rejected = peptide("REJECTED",900,"CAT3_c");
        var result = DdaMagPostAnalysisService.analyzeWithAcceptedGenomes(
                new MetaPeptide[]{shared,unique,rejected},1,"CAT",Set.of("CAT1","CAT2"));
        assertEquals(Set.of("CAT1"),result.finalSelection.genomes);
        assertEquals(2,result.peptideList.size());
        assertArrayEquals(new String[]{"CAT1_a"},shared.getProteins());
        assertArrayEquals(new String[]{"CAT1_a"},result.allProteins);
        assertEquals(150,result.genomeIntensityMap.get("CAT1")[0],1e-8);
        assertEquals(150,result.genomeUniqueIntensityMap.get("CAT1")[0],1e-8);
        assertEquals(150,result.allProteinIntensityMap.get("CAT1_a")[0],1e-8);
        assertEquals(Set.of("SHARED","UNIQUE"),result.fileSequenceSets[0]);
        assertEquals(1,result.finalSelection.excludedPeptides.size());
    }

    @Test void uniquePeptideIntensityExcludesSharedPeptides() {
        MetaPeptide cat1 = peptide("UNIQ1", 50, "CAT1_a");
        MetaPeptide cat2 = peptide("UNIQ2", 40, "CAT2_b");
        MetaPeptide shared = peptide("SHARED", 100, "CAT1_a", "CAT2_b");
        var result = DdaMagPostAnalysisService.analyzeWithAcceptedGenomes(
                new MetaPeptide[]{cat1, cat2, shared}, 1, "CAT", Set.of("CAT1", "CAT2"));
        assertEquals(Set.of("CAT1", "CAT2"), result.finalSelection.genomes);
        assertEquals(50, result.genomeUniqueIntensityMap.get("CAT1")[0], 1e-8);
        assertEquals(40, result.genomeUniqueIntensityMap.get("CAT2")[0], 1e-8);
        assertEquals(100, result.genomeIntensityMap.get("CAT1")[0], 1e-8);
        assertEquals(90, result.genomeIntensityMap.get("CAT2")[0], 1e-8);
    }
    @Test void nanConfidenceAndAbsentStatisticsCannotPassTheGate() {
        var selection = MagResultFinalisationService.select(Map.of("CAT1",new double[]{100}),
                Map.of("CAT1",new double[]{1,.1,Double.NaN}),Set.of("CAT1"),Map.of(),Map.of(),Map.of(),Map.of());
        assertTrue(selection.genomes.isEmpty());
    }
}
