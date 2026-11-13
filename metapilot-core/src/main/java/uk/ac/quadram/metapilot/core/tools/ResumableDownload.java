package uk.ac.quadram.metapilot.core.tools;

import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.net.HttpURLConnection;
import java.net.URI;
import java.net.URL;
import java.net.URLConnection;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.Objects;

import org.apache.commons.io.FileUtils;

/**
 * Headless resumable downloader (no JavaFX/Swing).
 *
 * Mirrors the behavior of GUI download helpers but reports progress via {@link ProgressCallback}.
 */
public final class ResumableDownload {

    private static final int DEFAULT_CONNECT_TIMEOUT_MS = 30000;
    private static final int DEFAULT_READ_TIMEOUT_MS = 60000;
    private static final int MAX_RETRY_ATTEMPTS = 5;
    private static final long PROGRESS_UPDATE_INTERVAL_MS = 500;
    private static final long SMALL_FILE_THRESHOLD = 1000;
    private static final int DEFAULT_BUFFER_SIZE = 8 * 1024 * 1024;
    private static final int SMALL_FILE_BUFFER_SIZE = 1024 * 1024;

    private ResumableDownload() {
    }

    public static boolean download(String url, File dest, boolean resume, ProgressCallback cb) throws IOException {
        return download(url, dest, resume, DEFAULT_CONNECT_TIMEOUT_MS, DEFAULT_READ_TIMEOUT_MS, MAX_RETRY_ATTEMPTS, cb);
    }

    /**
     * HTTP status of {@code url} (GET with a 1-byte Range, following redirects
     * including http→https). Negative when the request cannot be made.
     */
    public static int httpStatus(String url) {
        return httpStatus(url, 5);
    }

    public static boolean remoteExists(String url) {
        int code = httpStatus(url);
        return code >= 200 && code < 300;
    }

    /**
     * Remote file size in bytes, following redirects. Negative when the server
     * does not report it.
     */
    public static long contentLength(String url) {
        return getContentLength(url, DEFAULT_CONNECT_TIMEOUT_MS, DEFAULT_READ_TIMEOUT_MS);
    }

    public static boolean isNotFound(IOException e) {
        for (Throwable t = e; t != null; t = t.getCause()) {
            if (t instanceof FileNotFoundException) {
                return true;
            }
            String m = t.getMessage();
            if (m != null && m.contains("HTTP 404")) {
                return true;
            }
        }
        return false;
    }

    private static int httpStatus(String url, int redirectsLeft) {
        Objects.requireNonNull(url, "url");
        HttpURLConnection http = null;
        try {
            http = open(url, "GET");
            http.setInstanceFollowRedirects(false);
            http.setRequestProperty("Range", "bytes=0-0");
            int code = http.getResponseCode();
            if (redirectsLeft > 0 && isRedirect(code)) {
                String loc = http.getHeaderField("Location");
                if (loc != null && !loc.isBlank()) {
                    String next = URI.create(url).resolve(loc).toString();
                    http.disconnect();
                    http = null;
                    return httpStatus(next, redirectsLeft - 1);
                }
            }
            return code;
        } catch (FileNotFoundException e) {
            return HttpURLConnection.HTTP_NOT_FOUND;
        } catch (Exception e) {
            return -1;
        } finally {
            if (http != null) {
                http.disconnect();
            }
        }
    }

    private static boolean isRedirect(int code) {
        return code == HttpURLConnection.HTTP_MOVED_PERM
                || code == HttpURLConnection.HTTP_MOVED_TEMP
                || code == HttpURLConnection.HTTP_SEE_OTHER
                || code == 307
                || code == 308;
    }

    static IOException notFound(String url, IOException cause) {
        String hint = url.endsWith(".gz")
                ? ""
                : " (MGnify v3.0+ genome files are gzipped; expected " + url + ".gz)";
        return new IOException("HTTP 404 Not Found: " + url + hint, cause);
    }

    private static HttpURLConnection open(String url, String method) throws IOException {
        URLConnection connection;
        try {
            connection = new URI(url).toURL().openConnection();
        } catch (Exception e) {
            throw new IOException("Invalid URL: " + url, e);
        }
        if (!(connection instanceof HttpURLConnection)) {
            throw new IOException("Not an HTTP URL: " + url);
        }
        HttpURLConnection http = (HttpURLConnection) connection;
        http.setRequestMethod(method);
        http.setConnectTimeout(DEFAULT_CONNECT_TIMEOUT_MS);
        http.setReadTimeout(DEFAULT_READ_TIMEOUT_MS);
        http.setInstanceFollowRedirects(true);
        http.setRequestProperty("User-Agent", "MetaPilot/cli");
        return http;
    }

    public static boolean download(String url,
                                   File dest,
                                   boolean resume,
                                   int connectTimeoutMs,
                                   int readTimeoutMs,
                                   int maxRetries,
                                   ProgressCallback cb) throws IOException {
        Objects.requireNonNull(url, "url");
        Objects.requireNonNull(dest, "dest");
        if (cb == null) cb = new ProgressCallback.NoOp();

        File parent = dest.getParentFile();
        if (parent != null && !parent.exists()) {
            parent.mkdirs();
        }

        int attempt = 0;
        while (true) {
            try {
                cb.onMessage("Downloading: " + dest.getName());
                return attemptDownload(url, dest, resume, connectTimeoutMs, readTimeoutMs, cb);
            } catch (IOException e) {
                IOException wrapped = isNotFound(e) ? notFound(url, e) : e;
                attempt++;
                if (isNotFound(e) || attempt > maxRetries) {
                    cb.onError(wrapped);
                    throw wrapped;
                }
                cb.onMessage("Download failed (" + attempt + "/" + maxRetries + "): " + wrapped.getMessage());
                try {
                    Thread.sleep(1000L * attempt);
                } catch (InterruptedException ie) {
                    Thread.currentThread().interrupt();
                    throw new IOException("Interrupted during download retry", ie);
                }
            }
        }
    }

    private static boolean attemptDownload(String url,
                                           File dest,
                                           boolean resume,
                                           int connectTimeoutMs,
                                           int readTimeoutMs,
                                           ProgressCallback cb) throws IOException {
        URLConnection connection;
        try {
            connection = new URI(url).toURL().openConnection();
        } catch (Exception e) {
            throw new IOException("Invalid URL: " + url, e);
        }

        connection.setConnectTimeout(connectTimeoutMs);
        connection.setReadTimeout(readTimeoutMs);
        connection.setRequestProperty("User-Agent", "MetaPilot/cli");

        if (!(connection instanceof HttpURLConnection)) {
            FileUtils.copyURLToFile(new URL(url), dest, connectTimeoutMs, readTimeoutMs);
            cb.onProgress(1, 1);
            return true;
        }

        long fileLength = getContentLength(url, connectTimeoutMs, readTimeoutMs);
        if (fileLength > 0 && fileLength <= SMALL_FILE_THRESHOLD) {
            FileUtils.copyURLToFile(new URL(url), dest, connectTimeoutMs, readTimeoutMs);
            cb.onProgress(1, 1);
            return true;
        }

        File tmp = new File(dest.getAbsolutePath() + ".tmp");
        File downloadFile = resume ? dest : tmp;

        long startPosition = 0;
        if (resume && dest.exists() && fileLength > 0) {
            startPosition = dest.length();
            if (startPosition >= fileLength) {
                cb.onMessage("Already downloaded: " + dest.getName());
                cb.onProgress(1, 1);
                return true;
            }
        }

        HttpURLConnection http = (HttpURLConnection) connection;
        if (startPosition > 0) {
            http.setRequestProperty("Range", "bytes=" + startPosition + "-");
            cb.onMessage("Resuming at byte " + startPosition + " (" + dest.getName() + ")");
        }

        long bytesDownloaded = transferData(http, downloadFile, resume && startPosition > 0, startPosition, fileLength, cb);
        long totalBytes = startPosition + bytesDownloaded;

        if (fileLength > 0 && totalBytes < fileLength) {
            throw new IOException("Download incomplete: " + totalBytes + "/" + fileLength + " bytes for " + dest.getName());
        }

        if (!resume) {
            Files.move(tmp.toPath(), dest.toPath(), StandardCopyOption.REPLACE_EXISTING);
        }

        cb.onProgress(1, 1);
        cb.onMessage("Finished: " + dest.getName());
        return true;
    }

    private static long getContentLength(String url, int connectTimeoutMs, int readTimeoutMs) {
        HttpURLConnection head = null;
        try {
            head = (HttpURLConnection) new URI(url).toURL().openConnection();
            head.setRequestMethod("HEAD");
            head.setConnectTimeout(connectTimeoutMs);
            head.setReadTimeout(readTimeoutMs);
            head.setRequestProperty("User-Agent", "MetaPilot/cli");
            return head.getContentLengthLong();
        } catch (Exception e) {
            return -1;
        } finally {
            if (head != null) {
                head.disconnect();
            }
        }
    }

    private static long transferData(URLConnection conn,
                                     File outputFile,
                                     boolean append,
                                     long startPosition,
                                     long fileLength,
                                     ProgressCallback cb) throws IOException {
        long bytesDownloaded = 0;

        int bufferSize = (fileLength > 0 && fileLength > 100L * 1024L * 1024L) ? DEFAULT_BUFFER_SIZE : SMALL_FILE_BUFFER_SIZE;
        byte[] buffer = new byte[bufferSize];

        long lastUpdate = System.currentTimeMillis();

        try (BufferedInputStream bis = new BufferedInputStream(conn.getInputStream());
             BufferedOutputStream bos = new BufferedOutputStream(new FileOutputStream(outputFile, append))) {

            int bytesRead;
            while ((bytesRead = bis.read(buffer)) > 0) {
                if (Thread.currentThread().isInterrupted()) {
                    throw new IOException("Download interrupted");
                }
                bos.write(buffer, 0, bytesRead);
                bytesDownloaded += bytesRead;

                long now = System.currentTimeMillis();
                if (now - lastUpdate >= PROGRESS_UPDATE_INTERVAL_MS) {
                    if (fileLength > 0) {
                        long total = startPosition + bytesDownloaded;
                        cb.onProgress(total, fileLength);
                    }
                    lastUpdate = now;
                }
            }
            bos.flush();
        }

        if (fileLength > 0) {
            cb.onProgress(startPosition + bytesDownloaded, fileLength);
        }

        return bytesDownloaded;
    }
}

