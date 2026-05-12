package uk.ac.quadram.metapilot.core.MGYG;

import org.apache.commons.net.ftp.FTPClient;
import org.apache.commons.net.ftp.FTPFile;
import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Objects;

/**
 * Headless MGYG repository discovery/compare service (FTP listing only).
 *
 * This is the shared core counterpart to the JavaFX MGYG FTP task's:
 * - repository discovery (catalogues + latest version + species list)
 * - compare remote latest vs local configured items
 */
public final class MGYGRepositoryService {

    private static final String SERVER = "ftp.ebi.ac.uk";
    private static final int PORT = 21;
    private static final String USER = "anonymous";
    private static final String PASS = "";

    private static final String MGNIFY_ROOT = "/pub/databases/metagenomics/mgnify_genomes";

    private MGYGRepositoryService() {
    }

    public static MagDbItem[] discoverRepositories(ProgressCallback cb) throws IOException {
        if (cb == null) cb = new ProgressCallback.NoOp();

        FTPClient ftp = new FTPClient();
        List<MagDbItem> out = new ArrayList<>();

        try {
            cb.onMessage("Connecting to " + SERVER + " ...");
            ftp.connect(SERVER, PORT);
            ftp.login(USER, PASS);
            ftp.enterLocalPassiveMode();
            cb.onMessage("Connected.");

            FTPFile[] catalogues = ftp.listFiles(MGNIFY_ROOT);
            if (catalogues == null) {
                return new MagDbItem[0];
            }

            int total = catalogues.length;
            int processed = 0;

            for (FTPFile cat : catalogues) {
                if (Thread.currentThread().isInterrupted()) {
                    throw new IOException("Interrupted during repository discovery");
                }
                if (cat == null || !cat.isDirectory()) {
                    continue;
                }

                String catalogueId = cat.getName();
                String cataloguePath = MGNIFY_ROOT + "/" + catalogueId;

                processed++;
                cb.onMessage("Processing catalog: " + catalogueId);

                LatestVersion lv = findLatestVersion(ftp, cataloguePath);
                if (lv == null) {
                    cb.onMessage("Skipping (no versions): " + catalogueId);
                    continue;
                }

                SpeciesInfo info = listSpeciesCatalogue(ftp, catalogueId, lv.version, cataloguePath);

                MagDbItem item = new MagDbItem(catalogueId, "MGYG");
                item.resetMgygMap();
                if (info.parentToChildren != null) {
                    for (String parent : info.parentToChildren.keySet()) {
                        for (String child : info.parentToChildren.get(parent)) {
                            item.addMgygName(parent, child);
                        }
                    }
                }
                item.addVersion(lv.version, info.speciesCount, lv.date, false);
                out.add(item);

                cb.onProgress(processed, total);
            }

            cb.onMessage("Discovery completed (" + out.size() + " catalogues).");
            return out.toArray(new MagDbItem[0]);
        } catch (IOException e) {
            cb.onError(e);
            throw e;
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
     * Compare remote latest versions against local configured items.
     * Mirrors legacy behavior:
     * - If local catalogue+version matches latest: keep local item unchanged
     * - If catalogue exists but version differs: reuse local item object but reset map and add new version with localAvailable=true
     * - If new catalogue: create new MagDbItem and add new version with localAvailable=true
     */
    public static MagDbItem[] compareRepositories(MagDbItem[] localItems, ProgressCallback cb) throws IOException {
        Objects.requireNonNull(localItems, "localItems");
        if (cb == null) cb = new ProgressCallback.NoOp();

        HashMap<String, MagDbItem> localByCatalogue = new HashMap<>();
        for (MagDbItem li : localItems) {
            if (li != null && li.getCatalogueID() != null) {
                localByCatalogue.put(li.getCatalogueID(), li);
            }
        }

        FTPClient ftp = new FTPClient();
        List<MagDbItem> out = new ArrayList<>();

        try {
            cb.onMessage("Connecting to " + SERVER + " ...");
            ftp.connect(SERVER, PORT);
            ftp.login(USER, PASS);
            ftp.enterLocalPassiveMode();
            cb.onMessage("Connected.");

            FTPFile[] catalogues = ftp.listFiles(MGNIFY_ROOT);
            if (catalogues == null) {
                return new MagDbItem[0];
            }

            int total = catalogues.length;
            int processed = 0;

            for (FTPFile cat : catalogues) {
                if (Thread.currentThread().isInterrupted()) {
                    throw new IOException("Interrupted during repository comparison");
                }
                if (cat == null || !cat.isDirectory()) {
                    continue;
                }

                String catalogueId = cat.getName();
                String cataloguePath = MGNIFY_ROOT + "/" + catalogueId;
                processed++;

                cb.onMessage("Checking catalog: " + catalogueId);

                LatestVersion lv = findLatestVersion(ftp, cataloguePath);
                if (lv == null) {
                    continue;
                }

                MagDbItem local = localByCatalogue.get(catalogueId);
                if (local != null && safeEq(local.getCatalogueVersion(), lv.version)) {
                    out.add(local);
                    cb.onProgress(processed, total);
                    continue;
                }

                MagDbItem item = (local != null) ? local : new MagDbItem(catalogueId, "MGYG");
                item.resetMgygMap();

                SpeciesInfo info = listSpeciesCatalogue(ftp, catalogueId, lv.version, cataloguePath);
                if (info.parentToChildren != null) {
                    for (String parent : info.parentToChildren.keySet()) {
                        for (String child : info.parentToChildren.get(parent)) {
                            item.addMgygName(parent, child);
                        }
                    }
                }

                // localAvailable=true here means "update available" in the current GUI semantics
                item.addVersion(lv.version, info.speciesCount, lv.date, true);
                out.add(item);

                cb.onProgress(processed, total);
            }

            cb.onMessage("Comparison completed.");
            return out.toArray(new MagDbItem[0]);
        } catch (IOException e) {
            cb.onError(e);
            throw e;
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

    private static boolean safeEq(String a, String b) {
        if (a == null && b == null) return true;
        if (a == null || b == null) return false;
        return a.equals(b);
    }

    private static LatestVersion findLatestVersion(FTPClient ftp, String cataloguePath) throws IOException {
        FTPFile[] versions = ftp.listFiles(cataloguePath);
        if (versions == null || versions.length == 0) return null;

        String bestName = null;
        Calendar bestTs = null;
        Date bestDate = null;

        for (FTPFile f : versions) {
            if (f == null) continue;
            String name = f.getName();
            if (name == null) continue;

            // Skip obvious non-version entries if any (defensive)
            String lower = name.toLowerCase(Locale.ROOT);
            if (lower.equals(".") || lower.equals("..")) continue;

            Calendar ts = f.getTimestamp();
            if (bestTs == null || (ts != null && ts.after(bestTs))) {
                bestTs = ts;
                bestName = name;
                bestDate = ts != null ? ts.getTime() : null;
            }
        }

        if (bestName == null) return null;
        return new LatestVersion(bestName, bestDate);
    }

    private static SpeciesInfo listSpeciesCatalogue(FTPClient ftp, String catalogueId, String version, String cataloguePath) throws IOException {
        String versionPath = cataloguePath + "/" + version;
        FTPFile[] versionEntries = ftp.listFiles(versionPath);
        if (versionEntries == null) {
            return new SpeciesInfo(0, new HashMap<>());
        }

        HashMap<String, List<String>> parentToChildren = new HashMap<>();
        int speciesCount = 0;

        for (FTPFile ve : versionEntries) {
            if (Thread.currentThread().isInterrupted()) {
                throw new IOException("Interrupted during species catalogue listing");
            }
            if (ve == null) continue;
            if (!"species_catalogue".equals(ve.getName())) continue;

            String speciesRoot = versionPath + "/species_catalogue";
            FTPFile[] parents = ftp.listFiles(speciesRoot);
            if (parents == null) {
                break;
            }

            for (FTPFile p : parents) {
                if (Thread.currentThread().isInterrupted()) {
                    throw new IOException("Interrupted during species catalogue listing");
                }
                if (p == null || !p.isDirectory()) continue;
                String parent = p.getName();
                FTPFile[] children = ftp.listFiles(speciesRoot + "/" + parent);
                if (children == null) continue;
                for (FTPFile c : children) {
                    if (c == null || !c.isDirectory()) continue;
                    String child = c.getName();
                    parentToChildren.computeIfAbsent(parent, k -> new ArrayList<>()).add(child);
                    speciesCount++;
                }
            }
        }

        return new SpeciesInfo(speciesCount, parentToChildren);
    }

    private static final class LatestVersion {
        final String version;
        final Date date;

        LatestVersion(String version, Date date) {
            this.version = version;
            this.date = date;
        }
    }

    private static final class SpeciesInfo {
        final int speciesCount;
        final HashMap<String, List<String>> parentToChildren;

        SpeciesInfo(int speciesCount, HashMap<String, List<String>> parentToChildren) {
            this.speciesCount = speciesCount;
            this.parentToChildren = parentToChildren;
        }
    }
}

