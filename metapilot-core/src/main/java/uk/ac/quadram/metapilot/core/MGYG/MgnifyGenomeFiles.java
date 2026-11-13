package uk.ac.quadram.metapilot.core.MGYG;

import uk.ac.quadram.metapilot.core.tools.ResumableDownload;
import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.zip.GZIPInputStream;

/**
 * Per-genome MAG files as MGnify publishes them.
 * <p>
 * Human-gut v2.x serves uncompressed {@code .faa} and {@code _eggNOG.tsv}.
 * v3.0 (and later catalogues that follow it) gzip those same files. The local
 * MAG database always stores the uncompressed names, which is what the rest of
 * MetaPilot reads.
 */
final class MgnifyGenomeFiles {

    static final String FAA = ".faa";
    static final String EGGNOG = "_eggNOG.tsv";

    private MgnifyGenomeFiles() {
    }

    static String fileUrl(String httpBase, String parent, String genomeId, String suffix) {
        return httpBase + "/species_catalogue/" + parent + "/" + genomeId + "/genome/" + genomeId + suffix;
    }

    /**
     * True when the catalogue version gzips genome files. Probes the first
     * genome's protein FASTA once so every later download uses the matching URL.
     */
    static boolean remoteFilesGzipped(String httpBase, String parent, String genomeId, ProgressCallback cb)
            throws IOException {
        String plain = fileUrl(httpBase, parent, genomeId, FAA);
        String gz = plain + ".gz";
        int gzCode = ResumableDownload.httpStatus(gz);
        if (gzCode >= 200 && gzCode < 300) {
            if (cb != null) {
                cb.onMessage("Remote genome files are gzipped (" + genomeId + FAA + ".gz)");
            }
            return true;
        }
        int plainCode = ResumableDownload.httpStatus(plain);
        if (plainCode >= 200 && plainCode < 300) {
            if (cb != null) {
                cb.onMessage("Remote genome files are uncompressed (" + genomeId + FAA + ")");
            }
            return false;
        }
        throw new IOException("Neither " + gz + " (HTTP " + gzCode + ") nor " + plain
                + " (HTTP " + plainCode + ") is available");
    }

    /**
     * Downloads {@code uncompressedUrl} or {@code uncompressedUrl + ".gz"} into
     * {@code dest} as the uncompressed file. If the preferred layout 404s, the
     * other layout is tried immediately.
     */
    static boolean download(String uncompressedUrl, File dest, boolean gzipped, boolean resume, ProgressCallback cb)
            throws IOException {
        if (resume && dest.isFile() && dest.length() > 0) {
            File leftoverGz = new File(dest.getAbsolutePath() + ".gz");
            if (leftoverGz.isFile()) {
                leftoverGz.delete();
            }
            if (cb != null) {
                cb.onMessage("Already downloaded: " + dest.getName());
            }
            return true;
        }
        boolean useGz = gzipped || ResumableDownload.remoteExists(uncompressedUrl + ".gz");
        if (useGz) {
            try {
                return downloadGzipped(uncompressedUrl, dest, resume, cb);
            } catch (IOException e) {
                if (!ResumableDownload.isNotFound(e)) {
                    throw e;
                }
                if (cb != null) {
                    cb.onMessage("Gzipped file missing, trying uncompressed: " + uncompressedUrl);
                }
                return ResumableDownload.download(uncompressedUrl, dest, resume, cb);
            }
        }
        try {
            return ResumableDownload.download(uncompressedUrl, dest, resume, cb);
        } catch (IOException e) {
            if (!ResumableDownload.isNotFound(e)) {
                throw e;
            }
            if (cb != null) {
                cb.onMessage("Uncompressed file missing (HTTP 404), trying " + uncompressedUrl + ".gz");
            }
            return downloadGzipped(uncompressedUrl, dest, resume, cb);
        }
    }

    private static boolean downloadGzipped(String uncompressedUrl, File dest, boolean resume, ProgressCallback cb)
            throws IOException {
        File gz = new File(dest.getAbsolutePath() + ".gz");
        ResumableDownload.download(uncompressedUrl + ".gz", gz, resume, cb);
        gunzip(gz, dest);
        if (!gz.delete() && gz.isFile() && cb != null) {
            cb.onMessage("WARNING: could not delete " + gz.getAbsolutePath());
        }
        return true;
    }

    static void gunzip(File gz, File dest) throws IOException {
        File part = new File(dest.getAbsolutePath() + ".part");
        try (InputStream in = new GZIPInputStream(new BufferedInputStream(new FileInputStream(gz)));
             OutputStream out = new BufferedOutputStream(new FileOutputStream(part))) {
            in.transferTo(out);
        } catch (IOException e) {
            part.delete();
            gz.delete();
            throw new IOException("Failed to decompress " + gz.getName() + ": " + e.getMessage(), e);
        }
        Files.move(part.toPath(), dest.toPath(), StandardCopyOption.REPLACE_EXISTING);
    }
}
