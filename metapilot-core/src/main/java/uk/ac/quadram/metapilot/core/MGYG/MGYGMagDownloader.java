package uk.ac.quadram.metapilot.core.MGYG;

import org.apache.commons.net.ftp.FTPClient;
import org.apache.commons.net.ftp.FTPFile;
import uk.ac.quadram.metapilot.core.tools.ResumableDownload;
import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Objects;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

/**
 * Headless MGnify (MGYG) catalogue downloader for CLI/HPC.
 *
 * Download strategy:
 * - Use FTP only for directory listing (species catalogue and genome IDs)
 * - Use HTTP for actual file downloads (supports resume)
 * - Generate Marker2.fasta and Marker4.fasta after download
 */
public final class MGYGMagDownloader {

    private static final String SERVER = "ftp.ebi.ac.uk";
    private static final int PORT = 21;
    private static final String USER = "anonymous";
    private static final String PASS = "";

    private static final String MGNIFY_ROOT = "/pub/databases/metagenomics/mgnify_genomes";

    private static final String META_FILE = "genomes-all_metadata.tsv";

    private MGYGMagDownloader() {
    }

    public static void downloadCatalogue(String catalogueId,
                                         String version,
                                         File outputDir,
                                         boolean resume,
                                         int threads,
                                         int waitHours,
                                         ProgressCallback cb) throws IOException {
        Objects.requireNonNull(catalogueId, "catalogueId");
        Objects.requireNonNull(version, "version");
        Objects.requireNonNull(outputDir, "outputDir");
        if (cb == null) cb = new ProgressCallback.NoOp();

        if (threads <= 0) {
            threads = Math.max(1, Runtime.getRuntime().availableProcessors() / 2);
        }
        if (waitHours <= 0) {
            waitHours = 6;
        }

        File versionDir = new File(new File(outputDir, catalogueId), version);
        File originalDbDir = new File(versionDir, "original_db");
        File eggnogDir = new File(versionDir, "eggNOG");

        versionDir.mkdirs();
        originalDbDir.mkdirs();
        eggnogDir.mkdirs();

        String httpBase = "http://" + SERVER + MGNIFY_ROOT + "/" + catalogueId + "/" + version;

        cb.onMessage("MGYG download: " + catalogueId + " / " + version);
        cb.onMessage("Output: " + versionDir.getAbsolutePath());

        List<String> genomeEntries = listGenomeEntriesOverFtp(catalogueId, version, cb);
        if (genomeEntries.isEmpty()) {
            throw new IOException("No genomes found for " + catalogueId + "/" + version + " (FTP listing returned 0)");
        }

        cb.onMessage("Genomes to download: " + genomeEntries.size());

        ExecutorService exec = Executors.newFixedThreadPool(threads);
        List<Future<Boolean>> futures = new ArrayList<>();

        // 1) metadata file
        File localMeta = new File(versionDir, META_FILE);
        futures.add(exec.submit(downloadJob(httpBase + "/" + META_FILE, localMeta, resume, cb)));

        // 2) per-genome files
        for (String entry : genomeEntries) {
            if (Thread.currentThread().isInterrupted()) {
                throw new IOException("Interrupted before scheduling downloads");
            }
            String parent = parentFolder(entry);
            String genomeId = genomeId(entry);

            String faaUrl = httpBase + "/species_catalogue/" + parent + "/" + genomeId + "/genome/" + genomeId + ".faa";
            String eggUrl = httpBase + "/species_catalogue/" + parent + "/" + genomeId + "/genome/" + genomeId + "_eggNOG.tsv";

            futures.add(exec.submit(downloadJob(faaUrl, new File(originalDbDir, genomeId + ".faa"), resume, cb)));
            futures.add(exec.submit(downloadJob(eggUrl, new File(eggnogDir, genomeId + "_eggNOG.tsv"), resume, cb)));
        }

        exec.shutdown();

        cb.onMessage("Waiting for downloads to finish (" + waitHours + "h timeout)...");
        boolean finished;
        try {
            finished = exec.awaitTermination(waitHours, TimeUnit.HOURS);
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            throw new IOException("Interrupted while waiting for downloads", e);
        }
        if (!finished) {
            exec.shutdownNow();
            throw new IOException("Downloads did not finish within " + waitHours + " hours");
        }

        // Ensure all futures succeeded (will rethrow exceptions)
        for (Future<Boolean> f : futures) {
            try {
                Boolean ok = f.get();
                if (ok == null || !ok) {
                    throw new IOException("One or more downloads returned false");
                }
            } catch (Exception e) {
                throw new IOException("Download failed: " + e.getMessage(), e);
            }
        }

        cb.onMessage("Generating Marker2.fasta and Marker4.fasta...");
        try {
            HapExtractor extractor = new HapExtractor();
            extractor.extractTwoMarkers(originalDbDir, new File(versionDir, "Marker2.fasta"));
            extractor.extractFourCategories(originalDbDir, new File(versionDir, "Marker4.fasta"));
        } catch (IOException e) {
            throw new IOException("Downloaded, but marker FASTA generation failed: " + e.getMessage(), e);
        }

        cb.onMessage("MGYG download complete.");
        cb.onComplete(versionDir);
    }

    private static Callable<Boolean> downloadJob(String url, File dest, boolean resume, ProgressCallback cb) {
        return () -> ResumableDownload.download(url, dest, resume, cb);
    }

    /**
     * List genome IDs under species_catalogue for a catalogueId/version.
     *
     * The directory structure is:
     * /.../{catalogueId}/{version}/species_catalogue/{parent}/{genomeId}/
     *
     * We avoid storing the full 2-level tree and just emit genome IDs.
     */
    private static List<String> listGenomeEntriesOverFtp(String catalogueId, String version, ProgressCallback cb) throws IOException {
        FTPClient ftp = new FTPClient();
        List<String> entries = new ArrayList<>();
        String speciesRoot = MGNIFY_ROOT + "/" + catalogueId + "/" + version + "/species_catalogue";

        try {
            cb.onMessage("Connecting FTP for listing: " + SERVER);
            ftp.connect(SERVER, PORT);
            ftp.login(USER, PASS);
            ftp.enterLocalPassiveMode();

            FTPFile[] parents = ftp.listFiles(speciesRoot);
            if (parents == null) return entries;

            for (FTPFile parent : parents) {
                if (Thread.currentThread().isInterrupted()) {
                    throw new IOException("Interrupted during FTP listing");
                }
                if (parent == null || !parent.isDirectory()) continue;
                String parentName = parent.getName();
                FTPFile[] children = ftp.listFiles(speciesRoot + "/" + parentName);
                if (children == null) continue;
                for (FTPFile child : children) {
                    if (Thread.currentThread().isInterrupted()) {
                        throw new IOException("Interrupted during FTP listing");
                    }
                    if (child != null && child.isDirectory()) {
                        entries.add(parentName + "/" + child.getName());
                    }
                }
            }

            cb.onMessage("FTP listing done.");
            return entries;
        } catch (IOException e) {
            throw new IOException("FTP listing failed: " + e.getMessage(), e);
        } finally {
            try {
                if (ftp.isConnected()) {
                    ftp.logout();
                    ftp.disconnect();
                }
            } catch (Exception ignored) {
            }
        }
    }

    /**
     * MGYG uses a parent folder under species_catalogue. In the GUI task this is tracked explicitly.
     * For CLI download we reconstruct it from FTP listing and keep it as the immediate parent in the path.
     *
     * Current implementation assumes the parent folder is the first 2 characters of the genome id group
     * is NOT reliable, so we instead compute it lazily by mapping using the FTP listing path itself.
     *
     * Here we keep a simple fallback: parent folder == genomeId.substring(0, 2) is NOT valid for MGYG ids,
     * so callers should ensure listGenomeIdsOverFtp returns ids and we still need the parent folder.
     *
     * For now, parent folder is determined by scanning the FTP structure; we encode it into the id string
     * as "parent/genomeId" via {@link #parentFolder(String)}.
     */
    private static String parentFolder(String genomeEntry) {
        int slash = genomeEntry.indexOf('/');
        if (slash > 0) {
            return genomeEntry.substring(0, slash);
        }
        return genomeEntry;
    }

    private static String genomeId(String genomeEntry) {
        int slash = genomeEntry.indexOf('/');
        if (slash > 0 && slash + 1 < genomeEntry.length()) {
            return genomeEntry.substring(slash + 1);
        }
        return genomeEntry;
    }
}

