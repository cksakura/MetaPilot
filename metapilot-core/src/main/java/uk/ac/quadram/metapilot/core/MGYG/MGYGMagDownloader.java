package uk.ac.quadram.metapilot.core.MGYG;

import org.apache.commons.net.ftp.FTPClient;
import org.apache.commons.net.ftp.FTPFile;
import uk.ac.quadram.metapilot.core.tools.ResumableDownload;
import uk.ac.quadram.metapilot.core.variant.ProteinCatalogueInstaller;
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
 * - Accept both uncompressed genome files (human-gut v2.x) and gzipped ones
 *   (v3.0+); store uncompressed {@code .faa} / {@code _eggNOG.tsv} locally
 * - Generate Marker2.fasta and Marker4.fasta after download
 * - Optionally download and prepare the UHGP-100 protein catalogue (human-gut only)
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
        downloadCatalogue(catalogueId, version, outputDir, resume, threads, waitHours, false, true, cb);
    }

    /**
     * @param proteinCatalogue     also download and prepare the UHGP-100 protein
     *                             catalogue (human-gut v2.0.2 or v3.0) for the
     *                             pangenome variant check. A failure there is
     *                             reported and leaves the MAG database usable.
     * @param deleteProteinArchive delete {@code uhgp-100.tar.gz} once prepared
     */
    public static void downloadCatalogue(String catalogueId,
                                         String version,
                                         File outputDir,
                                         boolean resume,
                                         int threads,
                                         int waitHours,
                                         boolean proteinCatalogue,
                                         boolean deleteProteinArchive,
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

        String httpBase = httpBase(catalogueId, version);

        cb.onMessage("MGYG download: " + catalogueId + " / " + version);
        cb.onMessage("Output: " + versionDir.getAbsolutePath());

        List<String> genomeEntries = listGenomeEntriesOverFtp(catalogueId, version, cb);
        if (genomeEntries.isEmpty()) {
            throw new IOException("No genomes found for " + catalogueId + "/" + version + " (FTP listing returned 0)");
        }

        cb.onMessage("Genomes to download: " + genomeEntries.size());

        String probe = genomeEntries.get(0);
        boolean gzipped = MgnifyGenomeFiles.remoteFilesGzipped(
                httpBase, parentFolder(probe), genomeId(probe), cb);

        ExecutorService exec = Executors.newFixedThreadPool(threads);
        List<Future<Boolean>> futures = new ArrayList<>();

        // 1) metadata file
        File localMeta = new File(versionDir, META_FILE);
        futures.add(exec.submit(downloadJob(httpBase + "/" + META_FILE, localMeta, resume, cb)));

        // 2) per-genome files (v3.0+ are gzipped remotely; stored uncompressed locally)
        for (String entry : genomeEntries) {
            if (Thread.currentThread().isInterrupted()) {
                throw new IOException("Interrupted before scheduling downloads");
            }
            String parent = parentFolder(entry);
            String genomeId = genomeId(entry);
            String faaUrl = MgnifyGenomeFiles.fileUrl(httpBase, parent, genomeId, MgnifyGenomeFiles.FAA);
            String eggUrl = MgnifyGenomeFiles.fileUrl(httpBase, parent, genomeId, MgnifyGenomeFiles.EGGNOG);

            futures.add(exec.submit(genomeFileJob(
                    faaUrl, new File(originalDbDir, genomeId + ".faa"), gzipped, resume, cb)));
            futures.add(exec.submit(genomeFileJob(
                    eggUrl, new File(eggnogDir, genomeId + "_eggNOG.tsv"), gzipped, resume, cb)));
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

        if (proteinCatalogue) {
            if (!ProteinCatalogueInstaller.supports(catalogueId)) {
                cb.onMessage("UHGP-100 protein catalogue skipped: only available for the "
                        + ProteinCatalogueInstaller.SUPPORTED_CATALOGUE + " catalogue");
            } else {
                try {
                    ProteinCatalogueInstaller.install(httpBase, versionDir, deleteProteinArchive, cb);
                } catch (IOException | RuntimeException e) {
                    cb.onMessage("WARNING: the MAG database is complete and usable, but the UHGP-100 protein "
                            + "catalogue failed (" + e.getMessage() + "). Retry with 'Prepare protein catalogue' "
                            + "(GUI) or mag-download --protein-catalogue-only (CLI); the archive download resumes.");
                }
            }
        }

        cb.onMessage("MGYG download complete.");
        cb.onComplete(versionDir);
    }

    /**
     * Downloads and prepares only the UHGP-100 protein catalogue, for a MAG
     * database that is already downloaded.
     *
     * @param outputDir the same base directory as {@link #downloadCatalogue}; the
     *                  database is {@code <outputDir>/<catalogueId>/<version>/}
     */
    public static void prepareProteinCatalogue(String catalogueId,
                                               String version,
                                               File outputDir,
                                               boolean deleteArchive,
                                               ProgressCallback cb) throws IOException {
        Objects.requireNonNull(catalogueId, "catalogueId");
        Objects.requireNonNull(version, "version");
        Objects.requireNonNull(outputDir, "outputDir");
        if (cb == null) cb = new ProgressCallback.NoOp();
        if (!ProteinCatalogueInstaller.supports(catalogueId)) {
            throw new IOException("The UHGP-100 protein catalogue is only available for the "
                    + ProteinCatalogueInstaller.SUPPORTED_CATALOGUE + " catalogue, not " + catalogueId);
        }
        File versionDir = new File(new File(outputDir, catalogueId), version);
        if (!versionDir.isDirectory()) {
            throw new IOException("MAG database not found: " + versionDir.getAbsolutePath());
        }
        File dir = ProteinCatalogueInstaller.install(httpBase(catalogueId, version), versionDir, deleteArchive, cb);
        cb.onComplete(dir);
    }

    private static String httpBase(String catalogueId, String version) {
        return "http://" + SERVER + MGNIFY_ROOT + "/" + catalogueId + "/" + version;
    }

    private static Callable<Boolean> downloadJob(String url, File dest, boolean resume, ProgressCallback cb) {
        return () -> ResumableDownload.download(url, dest, resume, cb);
    }

    private static Callable<Boolean> genomeFileJob(String uncompressedUrl,
                                                   File dest,
                                                   boolean gzipped,
                                                   boolean resume,
                                                   ProgressCallback cb) {
        return () -> MgnifyGenomeFiles.download(uncompressedUrl, dest, gzipped, resume, cb);
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

