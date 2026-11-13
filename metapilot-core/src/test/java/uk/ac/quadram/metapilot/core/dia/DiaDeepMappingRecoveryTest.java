package uk.ac.quadram.metapilot.core.dia;

import static org.junit.jupiter.api.Assertions.*;
import java.nio.file.*;
import java.io.*;
import java.util.*;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import org.apache.parquet.example.data.simple.SimpleGroupFactory;
import org.apache.parquet.hadoop.example.ExampleParquetWriter;
import org.apache.parquet.schema.MessageTypeParser;
import org.apache.parquet.io.OutputFile;
import org.apache.parquet.io.PositionOutputStream;

class DiaDeepMappingRecoveryTest {
    @TempDir Path temp;

    @Test void finalCleavageProductsInheritAllContainingPeptideMappings() {
        var p = new uk.ac.quadram.metapilot.dbSearch.diann.DiaNNPrecursor("sample",
                new String[]{"arbitrary_combined_entry"}, "PEPTIDE", "PEPTIDE", 2, 0, 800,
                .001, .001, .001, .001, .001, 0, 10);
        var decoy = new uk.ac.quadram.metapilot.dbSearch.diann.DiaNNPrecursor("sample",
                new String[]{"shuf_CAT4_p"}, "DECOYPEP", "DECOYPEP", 2, 0, 800,
                .001, .001, .001, .001, .001, 0, 1);
        decoy.setIsDecoy(true);
        HashMap<String,HashSet<String>> map = new HashMap<>();
        map.put("KPEPTIDE", new HashSet<>(Set.of("CAT1_p")));
        map.put("PEPTIDE", new HashSet<>(Set.of("CAT2_p")));
        map.put("PEPTIDER", new HashSet<>(Set.of("CAT3_p")));
        var request = DiaDeepProteinMappingService.Request.prebuiltMap(map);
        var result = DiaDeepProteinMappingService.remapPrecursors(
                new uk.ac.quadram.metapilot.dbSearch.diann.DiaNNPrecursor[]{p,decoy}, request,null,"test");
        assertEquals(Set.of("CAT1_p","CAT2_p","CAT3_p"), Set.of(result.precursors[0].getProGroups()));
        assertEquals("shuf_CAT4_p", result.precursors[1].getProGroups()[0]);
        assertFalse(map.containsKey("DECOYPEP"));
    }

    private void parquet(Path file, String seq, String proteins, float q, boolean decoy) throws Exception {
        Files.createDirectories(file.getParent());
        var schema = MessageTypeParser.parseMessageType("message diann { "
                + "required binary Run (UTF8); required binary Protein.Ids (UTF8); "
                + "required binary Modified.Sequence (UTF8); required binary Stripped.Sequence (UTF8); "
                + "required int64 Precursor.Charge; required int64 Decoy; "
                + "required float Precursor.Mz; required float Precursor.Normalised; required float Precursor.Quantity; "
                + "required float Ms1.Normalised; required float Q.Value; required float PEP; "
                + "required float Global.Q.Value; required float PG.Q.Value; required float Global.PG.Q.Value; }");
        OutputFile output = new OutputFile() {
            public PositionOutputStream create(long hint) throws IOException { return stream(); }
            public PositionOutputStream createOrOverwrite(long hint) throws IOException { return stream(); }
            public boolean supportsBlockSize() { return false; }
            public long defaultBlockSize() { return 0; }
            private PositionOutputStream stream() throws IOException {
                var out = Files.newOutputStream(file);
                return new PositionOutputStream() {
                    long position;
                    public long getPos() { return position; }
                    public void write(int b) throws IOException { out.write(b); position++; }
                    public void write(byte[] b, int offset, int length) throws IOException {
                        out.write(b, offset, length); position += length;
                    }
                    public void close() throws IOException { out.close(); }
                };
            }
        };
        try (var writer = ExampleParquetWriter.builder(output).withType(schema)
                .withCompressionCodec(org.apache.parquet.hadoop.metadata.CompressionCodecName.UNCOMPRESSED).build()) {
            var row = new SimpleGroupFactory(schema).newGroup().append("Run", "sample")
                    .append("Protein.Ids", proteins).append("Modified.Sequence", seq)
                    .append("Stripped.Sequence", seq).append("Precursor.Charge", 2L)
                    .append("Decoy", decoy ? 1L : 0L);
            for (String f : List.of("Precursor.Mz", "Precursor.Normalised", "Precursor.Quantity", "Ms1.Normalised")) row.append(f, 100f);
            for (String f : List.of("Q.Value", "PEP", "Global.Q.Value", "PG.Q.Value", "Global.PG.Q.Value")) row.append(f, q);
            writer.write(row);
        }
    }

    @Test void roundMappingDoesNotApplyAnIndependentIdentificationThreshold() throws Exception {
        parquet(temp.resolve("round1/round1.parquet"), "QMAQDVTSMVLALR", "CAT1_p", .02148f, false);
        parquet(temp.resolve("round2/round2.parquet"), "QMAQDVTSMVLALR", "CAT2_p", .01355f, false);
        parquet(temp.resolve("round3/round3.parquet"), "QMAQDVTSMVLALR", "shuf_CAT3_p", .001f, true);
        var map = DiaDeepProteinMappingService.buildPeptideProteinMap(
                DiaDeepProteinMappingService.Request.deepTogether(temp.toFile(), 0, null, "K*,R*"), null, "test");
        assertEquals(Set.of("CAT1_p", "CAT2_p"), map.get("QMAQDVTSMVLALR"));
    }

    @Test void perFileMappingIncludesInitialEvidenceFromSamplesWithoutRounds() throws Exception {
        parquet(temp.resolve("marker/nonrepresentative/marker.parquet"), "CDMVDDPEMLELVEMEMR", "CAT1_p", .00987f, false);
        parquet(temp.resolve("marker/representative/round1/round1.parquet"), "OTHERPEPTIDE", "CAT1_q", .001f, false);
        // Genome CAT9 was never searched in any round: its marker evidence must not come back.
        parquet(temp.resolve("marker/unsearched/marker.parquet"), "NEVERSEARCHEDK", "CAT9_p", .001f, false);
        var map = DiaDeepProteinMappingService.buildPeptideProteinMap(
                DiaDeepProteinMappingService.Request.deepPerFile(temp.toFile(), null, "K*,R*"), null, "test");
        assertEquals(Set.of("CAT1_p"), map.get("CDMVDDPEMLELVEMEMR"));
        assertEquals(Set.of("CAT1_q"), map.get("OTHERPEPTIDE"));
        assertFalse(map.containsKey("NEVERSEARCHEDK"));
    }
}
