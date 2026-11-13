package uk.ac.quadram.metapilot.core.opensearch;

import static org.junit.jupiter.api.Assertions.*;

import java.io.File;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.attribute.FileTime;
import java.util.List;
import java.util.zip.GZIPOutputStream;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import uk.ac.quadram.metapilot.core.dia.EntrapmentShuffle;
import uk.ac.quadram.metapilot.core.prodb.TargetDecoyFasta;

class OpenSearchEntrapmentTest {
    @TempDir Path dir;
    private static final String SEQ = "MACDEFGHILMNPQSTVWYRKACDEFGHILMNPQSTVWYK";

    private File source() throws Exception {
        Path source = dir.resolve("real.fasta");
        Files.writeString(source, ">MGYG000000001_00001 description\n" + SEQ + "\n");
        return source.toFile();
    }

    @Test void bothPopulationsReceiveDecoysAndRealInputIsUntouched() throws Exception {
        File source = source(), searched = dir.resolve("searched.fasta").toFile();
        String original = Files.readString(source.toPath());
        OpenSearchEntrapment.prepare(source, searched, new OpenSearchEntrapment.Config(1, 8));
        assertEquals(original, Files.readString(source.toPath()));
        List<String> lines = Files.readAllLines(searched.toPath());
        assertEquals(4, lines.stream().filter(l -> l.startsWith(">")).count());
        for (String id : List.of("MGYG000000001_00001", "MGYG800000001_00001",
                "rev_MGYG000000001_00001", "rev_MGYG800000001_00001")) {
            assertTrue(lines.stream().anyMatch(l -> l.equals(">" + id) || l.startsWith(">" + id + " ")), id);
        }
        assertTrue(lines.contains(SEQ));
        assertTrue(lines.contains(EntrapmentShuffle.shuffle(SEQ, "MGYG800000001_00001", 1)));
        TargetDecoyFasta.Counts n = TargetDecoyFasta.count(searched);
        assertEquals(2, n.targets);
        assertEquals(2, n.decoys);
        assertTrue(Files.readString(Path.of(searched + ".entrapment.tsv")).contains("shuffled_proteins\t1"));
    }

    @Test void resumeKeepsTimestampsButSeedAndTamperingInvalidate() throws Exception {
        File source = source(), searched = dir.resolve("searched.fasta").toFile();
        OpenSearchEntrapment.prepare(source, searched, new OpenSearchEntrapment.Config(1, 8));
        FileTime first = Files.getLastModifiedTime(searched.toPath());
        String firstBytes = Files.readString(searched.toPath());
        OpenSearchEntrapment.prepare(source, searched, new OpenSearchEntrapment.Config(1, 8));
        assertEquals(first, Files.getLastModifiedTime(searched.toPath()));
        assertEquals(firstBytes, Files.readString(searched.toPath()));
        OpenSearchEntrapment.prepare(source, searched, new OpenSearchEntrapment.Config(2, 8));
        assertNotEquals(firstBytes, Files.readString(searched.toPath()));
        String secondBytes = Files.readString(searched.toPath());
        Files.writeString(searched.toPath(), secondBytes.replace(SEQ, "A".repeat(SEQ.length())));
        OpenSearchEntrapment.prepare(source, searched, new OpenSearchEntrapment.Config(2, 8));
        assertEquals(secondBytes, Files.readString(searched.toPath()));
    }

    @Test void changedSourceAndNamespaceRebuildWithoutAccumulatingControls() throws Exception {
        File source = source(), searched = dir.resolve("searched.fasta").toFile();
        OpenSearchEntrapment.prepare(source, searched, new OpenSearchEntrapment.Config(1, 8));
        Files.writeString(source.toPath(), ">MGYG000000001_00002\n" + SEQ + "\n");
        OpenSearchEntrapment.prepare(source, searched, new OpenSearchEntrapment.Config(1, 7));
        String content = Files.readString(searched.toPath());
        assertTrue(content.contains(">MGYG700000001_00002"));
        assertFalse(content.contains("MGYG800"));
        assertFalse(content.contains("00001 description"));
        assertEquals(2, TargetDecoyFasta.count(searched).targets);
    }

    @Test void invalidInputsFailBeforeReplacingSearchedDatabase() throws Exception {
        File source = source(), searched = dir.resolve("searched.fasta").toFile();
        Files.writeString(searched.toPath(), "prior database");
        for (String header : List.of("MGYG800000001_00001", "rev_MGYG000000001_00001", "sp|P12345|HUMAN")) {
            Files.writeString(source.toPath(), ">" + header + "\n" + SEQ + "\n");
            assertThrows(java.io.IOException.class,
                    () -> OpenSearchEntrapment.prepare(source, searched, new OpenSearchEntrapment.Config(1, 8)));
            assertEquals("prior database", Files.readString(searched.toPath()));
        }
        assertThrows(IllegalArgumentException.class, () -> new OpenSearchEntrapment.Config(1, 0));
        assertThrows(IllegalArgumentException.class, () -> new OpenSearchEntrapment.Config(1, 10));
        Files.writeString(source.toPath(), ">MGYG000000001_00001\n\n");
        assertThrows(java.io.IOException.class,
                () -> OpenSearchEntrapment.prepare(source, searched, new OpenSearchEntrapment.Config(1, 8)));
        assertEquals("prior database", Files.readString(searched.toPath()));
    }

    @Test void summarySeparatesControlsSharedMappingsAndDecoysWithoutEstimatingFdp() throws Exception {
        gzip("identification_arm_peptides.tsv.gz", "peptide\tproteins\n"
                + "PEPTIDEK\tMGYG000000001_00001\nENTRAPK\tMGYG800000001_00001\n"
                + "SHAREDK\tMGYG000000001_00001;MGYG800000001_00001\n");
        gzip("identification_arms.tsv.gz", "arm\tfile\tscan\tpeptide\tis_target\tclass\tq_value\taccepted\n"
                + "joint\trun\t1\tPEPTIDEK\ttrue\texact\t0.001\ttrue\n"
                + "joint\trun\t2\tENTRAPK\ttrue\texact\t0.001\ttrue\n"
                + "joint\trun\t3\tSHAREDK\ttrue\texact\t0.001\ttrue\n"
                + "joint\trun\t4\tDECOYK\tfalse\texact\t0.001\ttrue\n"
                + "joint\trun\t5\tENTRAPK\ttrue\texact\t0.04\tfalse\n");
        OpenSearchEntrapment.writeSummary(dir.toFile(), new OpenSearchEntrapment.Config(1, 8));
        String summary = Files.readString(dir.resolve("entrapment_summary.tsv"));
        for (String label : List.of("real_accessions", "entrapment_accessions", "shared_accessions", "decoy")) {
            assertTrue(summary.contains("joint\t" + label + "\t1\t1"));
        }
        assertTrue(summary.contains("NOT an FDP estimate"));
    }

    private void gzip(String name, String text) throws Exception {
        try (PrintWriter writer = new PrintWriter(new GZIPOutputStream(Files.newOutputStream(dir.resolve(name))),
                false, StandardCharsets.UTF_8)) { writer.print(text); }
    }

    @Test void decoyMethodSwitchRebuildsTheSearchedDatabase() throws Exception {
        File source = source(), searched = dir.resolve("searched.fasta").toFile();
        OpenSearchEntrapment.Config config = new OpenSearchEntrapment.Config(1, 8);
        OpenSearchEntrapment.prepare(source, searched, config);
        assertEquals(TargetDecoyFasta.DecoyMethod.REVERSE, TargetDecoyFasta.methodOf(searched));
        OpenSearchEntrapment.prepare(source, searched, config, TargetDecoyFasta.DecoyMethod.PSEUDO_REVERSE);
        assertEquals(TargetDecoyFasta.DecoyMethod.PSEUDO_REVERSE, TargetDecoyFasta.methodOf(searched));
        assertTrue(Files.readString(Path.of(searched.getPath() + ".entrapment.tsv")).contains("decoy_method\tpseudo_reverse\n"));
        long stamp = searched.lastModified();
        OpenSearchEntrapment.prepare(source, searched, config, TargetDecoyFasta.DecoyMethod.PSEUDO_REVERSE);
        assertEquals(stamp, searched.lastModified(), "same method reuses the database");
        OpenSearchEntrapment.prepare(source, searched, config);
        assertEquals(TargetDecoyFasta.DecoyMethod.REVERSE, TargetDecoyFasta.methodOf(searched));
    }
}
