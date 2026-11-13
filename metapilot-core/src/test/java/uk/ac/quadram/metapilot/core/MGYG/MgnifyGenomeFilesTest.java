package uk.ac.quadram.metapilot.core.MGYG;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.net.InetSocketAddress;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.zip.GZIPOutputStream;

import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import com.sun.net.httpserver.HttpServer;

class MgnifyGenomeFilesTest {

    private static final String PARENT = "MGYG0000000";
    private static final String ID = "MGYG000000001";
    private static final byte[] FAA = ">MGYG000000001_1\nMPEPTIDEK\n".getBytes(StandardCharsets.US_ASCII);
    private static final byte[] EGG = "MGYG000000001_1\tCOG1234\n".getBytes(StandardCharsets.US_ASCII);

    @TempDir
    File tmp;

    private HttpServer server;
    private String base;

    @BeforeEach
    void start() throws IOException {
        server = HttpServer.create(new InetSocketAddress("127.0.0.1", 0), 0);
        server.start();
        base = "http://127.0.0.1:" + server.getAddress().getPort();
    }

    @AfterEach
    void stop() {
        if (server != null) {
            server.stop(0);
        }
    }

    @Test
    @DisplayName("v3-style: gzipped remote files decompress to the uncompressed local names")
    void gzippedCatalogue() throws IOException {
        serve(MgnifyGenomeFiles.fileUrl("", PARENT, ID, MgnifyGenomeFiles.FAA) + ".gz", gzip(FAA));
        serve(MgnifyGenomeFiles.fileUrl("", PARENT, ID, MgnifyGenomeFiles.EGGNOG) + ".gz", gzip(EGG));
        serveNotFound(MgnifyGenomeFiles.fileUrl("", PARENT, ID, MgnifyGenomeFiles.FAA));

        assertTrue(MgnifyGenomeFiles.remoteFilesGzipped(base, PARENT, ID, null));

        File faa = new File(tmp, ID + ".faa");
        File egg = new File(tmp, ID + "_eggNOG.tsv");
        assertTrue(MgnifyGenomeFiles.download(
                MgnifyGenomeFiles.fileUrl(base, PARENT, ID, MgnifyGenomeFiles.FAA), faa, true, true, null));
        assertTrue(MgnifyGenomeFiles.download(
                MgnifyGenomeFiles.fileUrl(base, PARENT, ID, MgnifyGenomeFiles.EGGNOG), egg, true, true, null));
        assertEquals(new String(FAA, StandardCharsets.US_ASCII), Files.readString(faa.toPath()));
        assertEquals(new String(EGG, StandardCharsets.US_ASCII), Files.readString(egg.toPath()));
        assertFalse(new File(faa.getAbsolutePath() + ".gz").exists());
    }

    @Test
    @DisplayName("v2-style: uncompressed remote files are stored as-is")
    void uncompressedCatalogue() throws IOException {
        serve(MgnifyGenomeFiles.fileUrl("", PARENT, ID, MgnifyGenomeFiles.FAA), FAA);
        serveNotFound(MgnifyGenomeFiles.fileUrl("", PARENT, ID, MgnifyGenomeFiles.FAA) + ".gz");

        assertFalse(MgnifyGenomeFiles.remoteFilesGzipped(base, PARENT, ID, null));

        File faa = new File(tmp, ID + ".faa");
        assertTrue(MgnifyGenomeFiles.download(
                MgnifyGenomeFiles.fileUrl(base, PARENT, ID, MgnifyGenomeFiles.FAA), faa, false, true, null));
        assertEquals(new String(FAA, StandardCharsets.US_ASCII), Files.readString(faa.toPath()));
    }

    @Test
    @DisplayName("Resume skips a complete local file and does not hit the server")
    void resumeSkipsCompleteFile() throws IOException {
        File faa = new File(tmp, ID + ".faa");
        Files.write(faa.toPath(), FAA);
        serveNotFound(MgnifyGenomeFiles.fileUrl("", PARENT, ID, MgnifyGenomeFiles.FAA) + ".gz");
        assertTrue(MgnifyGenomeFiles.download(
                MgnifyGenomeFiles.fileUrl(base, PARENT, ID, MgnifyGenomeFiles.FAA), faa, true, true, null));
        assertEquals(new String(FAA, StandardCharsets.US_ASCII), Files.readString(faa.toPath()));
    }

    @Test
    @DisplayName("Probe fails clearly when neither .faa nor .faa.gz exists")
    void missingBothLayouts() {
        serveNotFound(MgnifyGenomeFiles.fileUrl("", PARENT, ID, MgnifyGenomeFiles.FAA));
        serveNotFound(MgnifyGenomeFiles.fileUrl("", PARENT, ID, MgnifyGenomeFiles.FAA) + ".gz");
        IOException e = assertThrows(IOException.class,
                () -> MgnifyGenomeFiles.remoteFilesGzipped(base, PARENT, ID, null));
        assertTrue(e.getMessage().contains(".faa.gz"));
        assertTrue(e.getMessage().contains("HTTP 404"));
    }

    @Test
    @DisplayName("If uncompressed 404s, the gzipped file is fetched (v3 layout, no probe)")
    void fallbackFromPlainToGzip() throws IOException {
        serve(MgnifyGenomeFiles.fileUrl("", PARENT, ID, MgnifyGenomeFiles.FAA) + ".gz", gzip(FAA));
        serveNotFound(MgnifyGenomeFiles.fileUrl("", PARENT, ID, MgnifyGenomeFiles.FAA));
        File faa = new File(tmp, ID + ".faa");
        assertTrue(MgnifyGenomeFiles.download(
                MgnifyGenomeFiles.fileUrl(base, PARENT, ID, MgnifyGenomeFiles.FAA), faa, false, true, null));
        assertEquals(new String(FAA, StandardCharsets.US_ASCII), Files.readString(faa.toPath()));
    }

    @Test
    @DisplayName("HTTP 404 is reported as such and is not retried")
    void notFoundMessage() {
        String path = MgnifyGenomeFiles.fileUrl("", PARENT, ID, MgnifyGenomeFiles.FAA);
        serveNotFound(path);
        serveNotFound(path + ".gz");
        IOException e = assertThrows(IOException.class,
                () -> MgnifyGenomeFiles.download(base + path, new File(tmp, ID + ".faa"), false, true, null));
        assertTrue(e.getMessage().contains("HTTP 404"), e.getMessage());
        assertTrue(e.getMessage().contains(".faa.gz"), e.getMessage());
    }

    private void serve(String path, byte[] body) {
        server.createContext(path, ex -> {
            if ("HEAD".equalsIgnoreCase(ex.getRequestMethod())) {
                ex.sendResponseHeaders(200, body.length);
                ex.close();
                return;
            }
            ex.sendResponseHeaders(200, body.length);
            try (OutputStream os = ex.getResponseBody()) {
                os.write(body);
            }
        });
    }

    private void serveNotFound(String path) {
        server.createContext(path, ex -> {
            ex.sendResponseHeaders(404, -1);
            ex.close();
        });
    }

    private static byte[] gzip(byte[] raw) throws IOException {
        ByteArrayOutputStream bos = new ByteArrayOutputStream();
        try (GZIPOutputStream gz = new GZIPOutputStream(bos)) {
            gz.write(raw);
        }
        return bos.toByteArray();
    }
}
