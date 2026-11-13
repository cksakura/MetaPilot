package uk.ac.quadram.metapilot.core.runner;

import static org.junit.jupiter.api.Assertions.*;
import java.io.File;
import java.nio.file.*;
import java.util.*;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import uk.ac.quadram.metapilot.core.dia.*;
import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;
import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNPrecursor;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;

class DiaFinalisationRegressionTest {
    @TempDir Path temp;
    private static DiaNNPrecursor precursor(String seq, String[] proteins, double score, boolean decoy) {
        DiaNNPrecursor p = new DiaNNPrecursor("S1",proteins,seq,seq,2,0,800,.001,.01,.001,.001,.001,0,score);
        p.setIsDecoy(decoy);
        return p;
    }
    @Test void seedSelectedScopeAndDirectScopePreserveDecoyBackground() throws Exception {
        Path predicted = temp.resolve("predicted");
        Files.createDirectories(predicted);
        Files.writeString(predicted.resolve("CAT1.predicted.tsv"),"Protein\tSequence\nCAT1_p\tPEPTIDE\n");
        Files.writeString(predicted.resolve("CAT2.predicted.tsv"),"Protein\tSequence\nCAT2_p\tPEPTIDE\n");
        DiaNNPrecursor p = precursor("PEPTIDE",new String[]{"library_entry"},10,false);
        DiaNNPrecursor decoy = precursor("DECOYPEP",new String[]{"shuf_library_entry"},1,true);
        var mapped = MagPredictedPeptideProteinMapper.mapForFinalisation(new DiaNNPrecursor[]{p,decoy}, predicted.toFile(),Set.of("CAT1"));
        assertEquals(Set.of("CAT1_p"),Set.of(mapped[0].getProGroups()));
        assertSame(decoy,mapped[1]);
        assertArrayEquals(new String[]{"shuf_library_entry"},decoy.getProGroups());
        mapped = MagPredictedPeptideProteinMapper.mapForFinalisation(new DiaNNPrecursor[]{p,decoy},predicted.toFile(),null);
        assertEquals(Set.of("CAT1_p","CAT2_p"),Set.of(mapped[0].getProGroups()));
    }
    @Test void deepUsesRecordedRoundMappingsAndRejectsMissingEvidence() {
        DiaNNPrecursor p = precursor("PEPTIDE",new String[]{"arbitrary_combined_entry"},10,false);
        HashMap<String,HashSet<String>> map = new HashMap<>();
        map.put("PEPTIDE",new HashSet<>(Set.of("CAT1_p","CAT2_p")));
        var result = DiaDeepProteinMappingService.remapPrecursors(new DiaNNPrecursor[]{p},
                DiaDeepProteinMappingService.Request.prebuiltMap(map),null,"test");
        assertEquals(Set.of("CAT1_p","CAT2_p"),Set.of(result.precursors[0].getProGroups()));
        assertThrows(IllegalStateException.class,() -> DiaDeepProteinMappingService.remapPrecursors(
                new DiaNNPrecursor[]{p},DiaDeepProteinMappingService.Request.prebuiltMap(new HashMap<>()),null,"test"));
    }
    @Test void sharedDiaExporterPrunesWeakGenomeAndNeverReadsUnrelatedPredictedFolder() throws Exception {
        Files.createDirectories(temp.resolve("predicted"));
        Files.writeString(temp.resolve("predicted/CAT9.predicted.tsv"),"Protein\tSequence\nCAT9_p\tPEPTIDEA\n");
        Path taxa = temp.resolve("taxa.tsv");
        Files.writeString(taxa,"Genome\tLineage\nCAT1\td__Bacteria;p__P;c__C;o__O;f__F;g__G;s__S\nCAT2\td__Bacteria\n");
        MagDbItem db = new MagDbItem("test","CAT") {
            @Override public File getOriginalTaxaFile() { return taxa.toFile(); }
            @Override public File getCurrentFile() { return temp.toFile(); }
        };
        List<DiaNNPrecursor> pp = new ArrayList<>();
        StringBuilder matrix = new StringBuilder("Stripped.Sequence\tModified.Sequence\tPrecursor.Charge\tProtein.Group\tS1.raw\n");
        for (int i=0;i<6;i++) {
            String seq="PEPTIDE"+(char)('A'+i);
            pp.add(precursor(seq,new String[]{"CAT1_p"},10,false));
            matrix.append(seq+"\t"+seq+"\t2\tCAT1_p\t100\n");
        }
        pp.add(precursor("SHAREDPEP",new String[]{"CAT2_p","CAT1_p"},10,false));
        matrix.append("SHAREDPEP\tSHAREDPEP\t2\tCAT2_p;CAT1_p\t100\n");
        pp.add(precursor("WEAKPEP",new String[]{"CAT2_p"},2,false));
        matrix.append("WEAKPEP\tWEAKPEP\t2\tCAT2_p\t900\n");
        for(int i=0;i<30;i++) pp.add(precursor("DECOY"+i,new String[]{"shuf_CAT1_p"},i % 2 == 0 ? 1 : 3,true));
        Path pr=temp.resolve("combined.pr_matrix.tsv");
        Files.writeString(pr,matrix.toString());
        Path out=temp.resolve("out"); Files.createDirectories(out);
        var result=DiaMagPostAnalysisTsvExporter.exportPrepared(out.toFile(),db,pr.toFile(),null,null,"test","Label-free",
                new ProgressCallback.NoOp(),pp.toArray(new DiaNNPrecursor[0]),false);
        assertEquals(Set.of("CAT1"),result.genomeIntensityMap.keySet());
        assertTrue(result.genomeDiaStats.get("CAT2")[1] < 1.0);
        assertTrue(result.genomeDiaStats.get("CAT2")[2] > 0.01);
        String peptides=Files.readString(out.resolve("final_peptides.tsv"));
        assertFalse(peptides.contains("CAT2_p"));
        assertFalse(peptides.contains("CAT9_p"));
        assertFalse(peptides.contains("WEAKPEP"));
        assertTrue(peptides.contains("SHAREDPEP"));
        assertEquals(8,Files.readAllLines(out.resolve("final_peptides.tsv")).size());
        assertEquals(700,result.genomeIntensityMap.get("CAT1")[0],1e-8);
        assertTrue(Files.readString(out.resolve("excluded_evidence/peptides.tsv")).contains("WEAKPEP"));
        String genome=Files.readString(out.resolve("taxonomy_analysis/Genome.tsv"));
        assertFalse(genome.contains("CAT2\t"));
        Path uniqueTsv=out.resolve("taxonomy_analysis/Genome_unique_peptides.tsv");
        assertTrue(Files.exists(uniqueTsv));
        List<String> genomeLines=Files.readAllLines(out.resolve("taxonomy_analysis/Genome.tsv"));
        List<String> uniqueLines=Files.readAllLines(uniqueTsv);
        assertEquals(genomeLines.get(0), uniqueLines.get(0));
        assertEquals(genomeLines.size(), uniqueLines.size());
        assertTrue(uniqueLines.get(1).startsWith("CAT1\t"));
        assertFalse(Files.readString(uniqueTsv).contains("CAT2\t"));
    }
    @Test void emptyFinalSelectionOverwritesEarlierGenomeAndPeptideTables() throws Exception {
        Path taxa=temp.resolve("taxa.tsv");
        Files.writeString(taxa,"Genome\tLineage\nCAT1\td__Bacteria\n");
        MagDbItem db=new MagDbItem("test","CAT") {
            @Override public File getOriginalTaxaFile() { return taxa.toFile(); }
        };
        Path matrix=temp.resolve("combined.pr_matrix.tsv");
        Files.writeString(matrix,"Stripped.Sequence\tModified.Sequence\tPrecursor.Charge\tProtein.Group\tS1.raw\nWEAKPEP\tWEAKPEP\t2\tCAT1_p\t100\n");
        Path out=temp.resolve("out"); Files.createDirectories(out.resolve("taxonomy_analysis"));
        Files.writeString(out.resolve("taxonomy_analysis/Genome.tsv"),"Genome\nSTALE\n");
        Files.writeString(out.resolve("taxonomy_analysis/genome_final_selected.tsv"),"Genome\nSTALE\n");
        DiaMagPostAnalysisTsvExporter.exportPrepared(out.toFile(),db,matrix.toFile(),null,null,"test","Label-free",
                new ProgressCallback.NoOp(), new DiaNNPrecursor[]{
                    precursor("WEAKPEP",new String[]{"CAT1_p"},.1,false),
                    precursor("DECOYPEP",new String[]{"shuf_CAT1_p"},1,true)},false);
        assertEquals(1, Files.readAllLines(out.resolve("final_peptides.tsv")).size());
        assertEquals(1, Files.readAllLines(out.resolve("taxonomy_analysis/Genome.tsv")).size());
        assertFalse(Files.readString(out.resolve("taxonomy_analysis/genome_final_selected.tsv")).contains("STALE"));
    }

}
