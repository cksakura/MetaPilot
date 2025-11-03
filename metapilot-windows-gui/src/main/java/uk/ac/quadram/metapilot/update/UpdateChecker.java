package uk.ac.quadram.metapilot.update;

import java.awt.Desktop;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.net.HttpURLConnection;
import java.net.URI;
import java.net.URL;
import java.nio.charset.StandardCharsets;
import java.time.Duration;
import java.util.Optional;
import java.util.logging.Level;
import java.util.logging.Logger;

import javafx.application.HostServices;
import javafx.application.Platform;
import javafx.scene.control.Alert;
import javafx.scene.control.Alert.AlertType;
import javafx.scene.control.ButtonType;
import javafx.scene.control.Hyperlink;
import javafx.scene.control.Label;
import javafx.scene.layout.VBox;

/**
 * Update checker for MetaPilot (GitHub-hosted .7z build).
 * Works with: https://github.com/cksakura/MetaPilot-Releases/releases
 */
public class UpdateChecker {
    private static final Logger LOGGER = Logger.getLogger(UpdateChecker.class.getName());

    // === Configuration ===
    private static final String GITHUB_OWNER = "cksakura";
    private static final String GITHUB_REPO  = "MetaPilot-Releases";
    private static final String GITHUB_API   = "https://api.github.com/repos/" + GITHUB_OWNER + "/" + GITHUB_REPO + "/releases?per_page=5";
    private static final String GITHUB_PAGE  = "https://github.com/" + GITHUB_OWNER + "/" + GITHUB_REPO + "/releases";
    private static final String WEBSITE_PAGE = "https://your-site.com/download"; // optional

    private static final int CONNECT_TIMEOUT_MS = (int) Duration.ofSeconds(6).toMillis();
    private static final int READ_TIMEOUT_MS    = (int) Duration.ofSeconds(8).toMillis();

    private final HostServices hostServices;
    private volatile String latestUrl;

    public UpdateChecker(HostServices hostServices) {
        this.hostServices = hostServices;
    }

    /** Check for newer version asynchronously */
    public void checkForUpdates(boolean showNoUpdatesMessage) {
        new Thread(() -> {
            try {
                String current = getCurrentVersion();
                ReleaseInfo latest = getLatestRelease();
                if (latest == null) {
                    if (showNoUpdatesMessage) showUpdateError();
                    return;
                }

                if (isNewerVersion(current, latest.version)) {
                    showUpdateAvailable(current, latest.version, latest.downloadUrl);
                } else if (showNoUpdatesMessage) {
                    showNoUpdatesAvailable(current);
                }
            } catch (Exception e) {
                LOGGER.log(Level.WARNING, "Update check failed", e);
                if (showNoUpdatesMessage) showUpdateError();
            }
        }, "MetaPilot-UpdateChecker").start();
    }

    // === Core ===

    private ReleaseInfo getLatestRelease() {
        try {
            URL url = new URL(GITHUB_API);
            HttpURLConnection conn = (HttpURLConnection) url.openConnection();
            conn.setRequestMethod("GET");
            conn.setConnectTimeout(CONNECT_TIMEOUT_MS);
            conn.setReadTimeout(READ_TIMEOUT_MS);
            conn.setRequestProperty("Accept", "application/vnd.github+json");
            conn.setRequestProperty("User-Agent", "MetaPilot-Updater");

            if (conn.getResponseCode() != HttpURLConnection.HTTP_OK) {
                LOGGER.warning("GitHub API status: " + conn.getResponseCode());
                return null;
            }

            String body;
            try (BufferedReader br = new BufferedReader(
                    new InputStreamReader(conn.getInputStream(), StandardCharsets.UTF_8))) {
                StringBuilder sb = new StringBuilder();
                String line;
                while ((line = br.readLine()) != null) sb.append(line);
                body = sb.toString();
            }

            // Find first non-draft release (pre-releases allowed)
            int idx = body.indexOf("{\"url\"");
            if (idx < 0) return null;
            String obj = body.substring(idx);

            String tag = extractJson(body, "\"tag_name\"");
            tag = stripV(tag);

            // Find first .7z asset
            String assetName = findAsset(body, ".7z");
            String url7z = (assetName != null) ? extractJsonNear(body, assetName, "\"browser_download_url\"") : null;
            if (url7z == null) url7z = GITHUB_PAGE;

            latestUrl = url7z;
            return new ReleaseInfo(tag, url7z);
        } catch (Exception e) {
            LOGGER.log(Level.WARNING, "GitHub query failed", e);
            return null;
        }
    }

    // === JSON helpers (minimal parser) ===

    private static String extractJson(String json, String key) {
        int i = json.indexOf(key);
        if (i < 0) return null;
        int colon = json.indexOf(':', i);
        int q1 = json.indexOf('"', colon + 1);
        int q2 = json.indexOf('"', q1 + 1);
        if (q1 < 0 || q2 < 0) return null;
        return json.substring(q1 + 1, q2);
    }

    private static String extractJsonNear(String json, String anchor, String key) {
        int pos = json.indexOf(anchor);
        if (pos < 0) return null;
        int sliceStart = Math.max(0, pos - 400);
        int sliceEnd = Math.min(json.length(), pos + 400);
        String slice = json.substring(sliceStart, sliceEnd);
        return extractJson(slice, key);
    }

    private static String findAsset(String json, String suffix) {
        int idx = 0;
        while ((idx = json.indexOf("\"name\"", idx)) >= 0) {
            String val = extractJson(json.substring(idx), "\"name\"");
            if (val != null && val.endsWith(suffix)) return val;
            idx += 6;
        }
        return null;
    }

    private static String stripV(String v) {
        return (v != null && v.matches("^[vV].*")) ? v.substring(1) : v;
    }

    // === Version comparison (supports pre-release) ===

    private static boolean isNewerVersion(String current, String latest) {
        SemVer c = SemVer.parse(current);
        SemVer l = SemVer.parse(latest);
        return l.compareTo(c) > 0;
    }

    private static class SemVer implements Comparable<SemVer> {
        int major, minor, patch;
        String[] pre;
        static SemVer parse(String v) {
            if (v == null) return new SemVer(0, 0, 0, new String[0]);
            v = v.trim();
            if (v.startsWith("v") || v.startsWith("V")) v = v.substring(1);
            String[] parts = v.split("-", 2);
            String[] nums = parts[0].split("\\.");
            int M = nums.length > 0 ? parseInt(nums[0]) : 0;
            int m = nums.length > 1 ? parseInt(nums[1]) : 0;
            int p = nums.length > 2 ? parseInt(nums[2]) : 0;
            String[] pre = (parts.length > 1) ? parts[1].split("\\.") : new String[0];
            return new SemVer(M, m, p, pre);
        }
        SemVer(int M, int m, int p, String[] pre) { this.major = M; this.minor = m; this.patch = p; this.pre = pre; }
        @Override public int compareTo(SemVer o) {
            if (major != o.major) return major - o.major;
            if (minor != o.minor) return minor - o.minor;
            if (patch != o.patch) return patch - o.patch;
            if (pre.length == 0 && o.pre.length > 0) return 1;
            if (pre.length > 0 && o.pre.length == 0) return -1;
            for (int i = 0; i < Math.max(pre.length, o.pre.length); i++) {
                String a = i < pre.length ? pre[i] : "";
                String b = i < o.pre.length ? o.pre[i] : "";
                if (a.equals(b)) continue;
                boolean ad = a.matches("\\d+"), bd = b.matches("\\d+");
                if (ad && bd) return Integer.compare(Integer.parseInt(a), Integer.parseInt(b));
                if (ad != bd) return ad ? -1 : 1;
                return a.compareTo(b);
            }
            return 0;
        }
        private static int parseInt(String s) { try { return Integer.parseInt(s); } catch (Exception e) { return 0; } }
    }

    // === Dialogs ===

    private void showUpdateAvailable(String current, String latest, String url) {
        Platform.runLater(() -> {
            Alert alert = new Alert(AlertType.INFORMATION);
            alert.setTitle("Update Available");
            alert.setHeaderText("A new version of MetaPilot is available!");

            VBox box = new VBox(10);
            Label info = new Label(
                "Current version: " + current + "\n" +
                "Latest version: " + latest + "\n\n" +
                "Would you like to download it?");
            Hyperlink link = new Hyperlink("Download Update (.7z)");
            link.setOnAction(e -> openExternal(url));
            Hyperlink site = new Hyperlink("Visit MetaPilot website");
            site.setOnAction(e -> openExternal(WEBSITE_PAGE));

            box.getChildren().addAll(info, link, site);
            alert.getDialogPane().setContent(box);
            alert.getButtonTypes().setAll(ButtonType.CLOSE);
            alert.showAndWait();
        });
    }

    private void showNoUpdatesAvailable(String current) {
        Platform.runLater(() -> {
            Alert a = new Alert(AlertType.INFORMATION);
            a.setTitle("Up to Date");
            a.setHeaderText("Your MetaPilot is up to date");
            a.setContentText("Current version: " + current);
            a.showAndWait();
        });
    }

    private void showUpdateError() {
        Platform.runLater(() -> {
            Alert a = new Alert(AlertType.WARNING);
            a.setTitle("Update Check Failed");
            a.setHeaderText("Could not contact GitHub");
            a.setContentText("Please check your connection or try again later.\nCurrent version: " + getCurrentVersion());
            a.showAndWait();
        });
    }

    private void openExternal(String url) {
        try {
            if (hostServices != null) {
                hostServices.showDocument(url);
            } else if (Desktop.isDesktopSupported()) {
                Desktop.getDesktop().browse(URI.create(url));
            }
        } catch (Exception e) {
            LOGGER.log(Level.WARNING, "Failed to open " + url, e);
        }
    }

    public static String getCurrentVersion() {
        String v = Optional.ofNullable(UpdateChecker.class.getPackage())
                .map(Package::getImplementationVersion).orElse(null);
        return (v != null && !v.isBlank()) ? v : "1.0.0-beta.2";
    }

    private static class ReleaseInfo {
        final String version;
        final String downloadUrl;
        ReleaseInfo(String v, String u) { version = v; downloadUrl = u; }
    }
}
