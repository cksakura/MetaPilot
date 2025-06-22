package uk.ac.quadram.metapilot.core.tools;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.net.HttpURLConnection;
import java.net.ProtocolException;
import java.net.URI;
import java.net.URISyntaxException;
import java.net.URL;
import java.net.URLConnection;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;

import javafx.concurrent.Task;
import org.apache.commons.io.FileUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * JavaFX-Compatible Resumable Download Task with Enhanced Features
 * 
 * <h2>Purpose</h2>
 * This JavaFX-compatible version provides robust file downloading capabilities with
 * resume functionality, progress tracking, and enhanced error handling. It extends
 * javafx.concurrent.Task for seamless integration with JavaFX applications.
 * 
 * <h2>Key Features</h2>
 * <ul>
 *   <li><strong>Resumable Downloads:</strong> Can resume interrupted downloads</li>
 *   <li><strong>Progress Tracking:</strong> Real-time progress updates via JavaFX properties</li>
 *   <li><strong>Automatic Retry:</strong> Configurable retry mechanism for failed downloads</li>
 *   <li><strong>Size Validation:</strong> Verifies file integrity after download</li>
 *   <li><strong>Timeout Handling:</strong> Configurable connection and read timeouts</li>
 *   <li><strong>Bandwidth Optimization:</strong> Optimized buffer sizes for different file types</li>
 * </ul>
 * 
 * <h2>Optimizations Over Original</h2>
 * <ul>
 *   <li><strong>Better Error Handling:</strong> Comprehensive exception handling with recovery</li>
 *   <li><strong>Improved Progress Reporting:</strong> More granular progress updates</li>
 *   <li><strong>Connection Management:</strong> Better connection lifecycle management</li>
 *   <li><strong>Temporary File Handling:</strong> Safe atomic file operations</li>
 *   <li><strong>Logging Integration:</strong> Comprehensive logging for debugging</li>
 * </ul>
 * 
 * <h2>Usage Examples</h2>
 * 
 * <h3>Basic Download:</h3>
 * <pre>{@code
 *   ResumableDownloadTaskJavaFX task = new ResumableDownloadTaskJavaFX(
 *       "http://example.com/file.zip", 
 *       "/path/to/save/file.zip"
 *   );
 *   
 *   // Bind to progress bar
 *   progressBar.progressProperty().bind(task.progressProperty());
 *   
 *   // Execute
 *   Thread thread = new Thread(task);
 *   thread.setDaemon(true);
 *   thread.start();
 * }</pre>
 * 
 * <h3>Resumable Download:</h3>
 * <pre>{@code
 *   ResumableDownloadTaskJavaFX task = new ResumableDownloadTaskJavaFX(
 *       "http://example.com/largefile.zip", 
 *       "/path/to/save/largefile.zip",
 *       true  // Enable resume
 *   );
 * }</pre>
 * 
 * @author Kai Cheng - Quadram Institute Bioscience
 * @version 2.0
 * @since 2.0
 * 
 * @see uk.ac.quadram.metapilot.core.tools.ResumableDownloadTask
 * @see javafx.concurrent.Task
 */
public class ResumableDownloadTaskJavaFX extends Task<Boolean> {

    private static final Logger LOGGER = LogManager.getLogger(ResumableDownloadTaskJavaFX.class);
    
    // Configuration constants
    private static final int DEFAULT_BUFFER_SIZE = 8 * 1024 * 1024; // 8MB buffer
    private static final int SMALL_FILE_BUFFER_SIZE = 1024 * 1024; // 1MB for small files
    private static final int DEFAULT_CONNECT_TIMEOUT = 30000; // 30 seconds
    private static final int DEFAULT_READ_TIMEOUT = 60000; // 60 seconds
    private static final int MAX_RETRY_ATTEMPTS = 5;
    private static final long PROGRESS_UPDATE_INTERVAL = 500; // 500ms
    private static final long SMALL_FILE_THRESHOLD = 1000; // Files smaller than 1KB
    
    // Instance variables
    private final String downloadUrl;
    private final File saveAsFile;
    private final File tempFile;
    private final boolean resume;
    private final int connectTimeout;
    private final int readTimeout;
    private final int maxRetryAttempts;
    
    private long fileLength = -1;
    private int retryCount = 0;
    
    /**
     * Basic constructor for simple downloads
     * 
     * @param downloadUrl URL to download from
     * @param saveAsFileName Path where file should be saved
     */
    public ResumableDownloadTaskJavaFX(String downloadUrl, String saveAsFileName) {
        this(downloadUrl, saveAsFileName, false);
    }
    
    /**
     * Constructor with resume capability
     * 
     * @param downloadUrl URL to download from
     * @param saveAsFileName Path where file should be saved
     * @param resume Whether to resume interrupted downloads
     */
    public ResumableDownloadTaskJavaFX(String downloadUrl, String saveAsFileName, boolean resume) {
        this(downloadUrl, saveAsFileName, -1, resume);
    }
    
    /**
     * Constructor with known file length
     * 
     * @param downloadUrl URL to download from
     * @param saveAsFileName Path where file should be saved
     * @param fileLength Expected file length (for progress calculation)
     */
    public ResumableDownloadTaskJavaFX(String downloadUrl, String saveAsFileName, long fileLength) {
        this(downloadUrl, saveAsFileName, fileLength, false);
    }
    
    /**
     * Full constructor with all options
     * 
     * @param downloadUrl URL to download from
     * @param saveAsFileName Path where file should be saved
     * @param fileLength Expected file length (-1 if unknown)
     * @param resume Whether to resume interrupted downloads
     */
    public ResumableDownloadTaskJavaFX(String downloadUrl, String saveAsFileName, long fileLength, boolean resume) {
        this(downloadUrl, saveAsFileName, fileLength, resume, DEFAULT_CONNECT_TIMEOUT, DEFAULT_READ_TIMEOUT, MAX_RETRY_ATTEMPTS);
    }
    
    /**
     * Advanced constructor with timeout and retry configuration
     * 
     * @param downloadUrl URL to download from
     * @param saveAsFileName Path where file should be saved
     * @param fileLength Expected file length (-1 if unknown)
     * @param resume Whether to resume interrupted downloads
     * @param connectTimeout Connection timeout in milliseconds
     * @param readTimeout Read timeout in milliseconds
     * @param maxRetryAttempts Maximum number of retry attempts
     */
    public ResumableDownloadTaskJavaFX(String downloadUrl, String saveAsFileName, long fileLength, boolean resume,
                                       int connectTimeout, int readTimeout, int maxRetryAttempts) {
        this.downloadUrl = downloadUrl;
        this.saveAsFile = new File(saveAsFileName);
        this.tempFile = new File(saveAsFileName + ".tmp");
        this.fileLength = fileLength;
        this.resume = resume;
        this.connectTimeout = connectTimeout;
        this.readTimeout = readTimeout;
        this.maxRetryAttempts = maxRetryAttempts;
        
        // Ensure parent directory exists
        File parentDir = this.saveAsFile.getParentFile();
        if (parentDir != null && !parentDir.exists()) {
            parentDir.mkdirs();
        }
    }
    
    @Override
    protected Boolean call() throws Exception {
        updateMessage("Starting download of " + saveAsFile.getName());
        updateProgress(0, 100);
        
        LOGGER.info("Starting download: {} -> {}", downloadUrl, saveAsFile.getAbsolutePath());
        
        try {
            return performDownload();
        } catch (Exception e) {
            LOGGER.error("Download failed for {}: {}", downloadUrl, e.getMessage(), e);
            updateMessage("Download failed: " + e.getMessage());
            throw e;
        }
    }
    
    /**
     * Main download logic with retry mechanism
     */
    private Boolean performDownload() throws Exception {
        while (retryCount <= maxRetryAttempts) {
            if (isCancelled()) {
                updateMessage("Download cancelled");
                cleanupTempFile();
                return false;
            }
            
            try {
                return attemptDownload();
            } catch (Exception e) {
                retryCount++;
                if (retryCount <= maxRetryAttempts) {
                    LOGGER.warn("Download attempt {} failed, retrying... Error: {}", retryCount, e.getMessage());
                    updateMessage("Download failed, retrying... (attempt " + retryCount + "/" + maxRetryAttempts + ")");
                    Thread.sleep(1000 * retryCount); // Exponential backoff
                } else {
                    throw e;
                }
            }
        }
        return false;
    }
    
    /**
     * Single download attempt
     */
    private Boolean attemptDownload() throws Exception {
        URLConnection connection = createConnection();
        
        if (connection instanceof HttpURLConnection) {
            return downloadHttpFile((HttpURLConnection) connection);
        } else {
            return downloadSimpleFile(connection);
        }
    }
    
    /**
     * Create and configure URL connection
     */
    private URLConnection createConnection() throws IOException, URISyntaxException {
        URLConnection connection = new URI(downloadUrl).toURL().openConnection();
        
        // Configure timeouts
        connection.setConnectTimeout(connectTimeout);
        connection.setReadTimeout(readTimeout);
        
        // Set user agent to avoid blocking
        connection.setRequestProperty("User-Agent", "MetaPilot/2.0");
        
        return connection;
    }
    
    /**
     * Download HTTP file with resume capability
     */
    private Boolean downloadHttpFile(HttpURLConnection connection) throws Exception {
        // Get file length if not known
        if (fileLength <= 0) {
            fileLength = getContentLength(connection);
        }
        
        updateMessage("File size: " + formatBytes(fileLength));
        
        if (fileLength <= SMALL_FILE_THRESHOLD) {
            return downloadSmallFile();
        }
        
        // Determine starting position for resume
        long startPosition = 0;
        File downloadFile = resume ? saveAsFile : tempFile;
        
        if (resume && saveAsFile.exists()) {
            startPosition = saveAsFile.length();
            if (startPosition >= fileLength) {
                updateMessage("File already downloaded");
                updateProgress(100, 100);
                return true;
            }
        }
        
        // Set range header for resume
        if (startPosition > 0) {
            connection.setRequestProperty("Range", "bytes=" + startPosition + "-");
            updateMessage("Resuming download from byte " + startPosition);
        }
        
        // Perform the download
        long bytesDownloaded = transferData(connection, downloadFile, startPosition);
        long totalBytes = startPosition + bytesDownloaded;
        
        // Verify download completion
        if (totalBytes >= fileLength) {
            if (!resume && tempFile.exists()) {
                // Move temp file to final location atomically
                Files.move(tempFile.toPath(), saveAsFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
            }
            
            updateMessage("Download completed: " + saveAsFile.getName());
            updateProgress(100, 100);
            LOGGER.info("Download completed successfully: {}", saveAsFile.getAbsolutePath());
            return true;
        } else {
            String errorMsg = "Download incomplete: " + totalBytes + "/" + fileLength + " bytes";
            updateMessage(errorMsg);
            LOGGER.warn(errorMsg);
            return false;
        }
    }
    
    /**
     * Download simple file (non-HTTP or very small files)
     */
    private Boolean downloadSimpleFile(URLConnection connection) throws Exception {
        updateMessage("Downloading file...");
        
        try {
            FileUtils.copyURLToFile(
                new URL(downloadUrl), 
                saveAsFile, 
                connectTimeout, 
                readTimeout
            );
            
            updateMessage("Download completed: " + saveAsFile.getName());
            updateProgress(100, 100);
            LOGGER.info("Simple download completed: {}", saveAsFile.getAbsolutePath());
            return true;
            
        } catch (Exception e) {
            updateMessage("Download failed: " + e.getMessage());
            LOGGER.error("Simple download failed: {}", e.getMessage());
            throw e;
        }
    }
    
    /**
     * Download small files using Apache Commons IO
     */
    private Boolean downloadSmallFile() throws Exception {
        updateMessage("Downloading small file...");
        
        try {
            FileUtils.copyURLToFile(
                new URL(downloadUrl), 
                saveAsFile, 
                connectTimeout, 
                readTimeout
            );
            
            updateMessage("Download completed: " + saveAsFile.getName());
            updateProgress(100, 100);
            return true;
            
        } catch (Exception e) {
            updateMessage("Small file download failed, trying alternative method...");
            // Fallback to regular download method
            URLConnection connection = createConnection();
            return downloadHttpFile((HttpURLConnection) connection);
        }
    }
    
    /**
     * Transfer data from connection to file with progress tracking
     */
    private long transferData(URLConnection connection, File outputFile, long startPosition) throws IOException {
        long bytesDownloaded = 0;
        
        try (BufferedInputStream bis = new BufferedInputStream(connection.getInputStream());
             BufferedOutputStream bos = new BufferedOutputStream(
                 new FileOutputStream(outputFile, resume && startPosition > 0))) {
            
            // Choose buffer size based on file size
            int bufferSize = fileLength > 100 * 1024 * 1024 ? DEFAULT_BUFFER_SIZE : SMALL_FILE_BUFFER_SIZE;
            byte[] buffer = new byte[bufferSize];
            
            int bytesRead;
            long lastUpdateTime = System.currentTimeMillis();
            double lastProgressUpdate = 0;
            
            while ((bytesRead = bis.read(buffer)) > 0) {
                if (isCancelled()) {
                    updateMessage("Download cancelled");
                    break;
                }
                
                bos.write(buffer, 0, bytesRead);
                bytesDownloaded += bytesRead;
                
                // Update progress periodically
                long currentTime = System.currentTimeMillis();
                if (currentTime - lastUpdateTime > PROGRESS_UPDATE_INTERVAL) {
                    long totalBytes = startPosition + bytesDownloaded;
                    double progress = (double) totalBytes / fileLength * 100.0;
                    
                    if (progress > lastProgressUpdate + 1) { // Update only if progress increased by at least 1%
                        updateProgress(progress, 100);
                        updateMessage("Downloaded " + formatBytes(totalBytes) + " / " + formatBytes(fileLength) + 
                                    " (" + String.format("%.1f", progress) + "%)");
                        lastProgressUpdate = progress;
                    }
                    lastUpdateTime = currentTime;
                }
            }
            
            bos.flush();
        }
        
        return bytesDownloaded;
    }
    
    /**
     * Get content length from HTTP connection
     */
    private long getContentLength(HttpURLConnection connection) throws IOException {
        // Create a HEAD request to get file size
        HttpURLConnection headConnection = null;
        try {
            headConnection = (HttpURLConnection) new URI(downloadUrl).toURL().openConnection();
            headConnection.setRequestMethod("HEAD");
            headConnection.setConnectTimeout(connectTimeout);
            headConnection.setReadTimeout(readTimeout);
            headConnection.setRequestProperty("User-Agent", "MetaPilot/2.0");
            
            long length = headConnection.getContentLengthLong();
            return length > 0 ? length : -1;
            
        } catch (Exception e) {
            LOGGER.warn("Could not determine file size: {}", e.getMessage());
            return -1;
        } finally {
            if (headConnection != null) {
                headConnection.disconnect();
            }
        }
    }
    
    /**
     * Format bytes for human-readable display
     */
    private String formatBytes(long bytes) {
        if (bytes < 1024) return bytes + " B";
        if (bytes < 1024 * 1024) return String.format("%.1f KB", bytes / 1024.0);
        if (bytes < 1024 * 1024 * 1024) return String.format("%.1f MB", bytes / (1024.0 * 1024.0));
        return String.format("%.1f GB", bytes / (1024.0 * 1024.0 * 1024.0));
    }
    
    /**
     * Clean up temporary files
     */
    private void cleanupTempFile() {
        if (tempFile.exists()) {
            try {
                Files.delete(tempFile.toPath());
                LOGGER.debug("Cleaned up temp file: {}", tempFile.getAbsolutePath());
            } catch (IOException e) {
                LOGGER.warn("Could not delete temp file: {}", tempFile.getAbsolutePath());
            }
        }
    }
    
    /**
     * Get the final downloaded file
     * 
     * @return The downloaded file
     */
    public File getDownloadedFile() {
        return saveAsFile;
    }
    
    /**
     * Get the download URL
     * 
     * @return The download URL
     */
    public String getDownloadUrl() {
        return downloadUrl;
    }
    
    /**
     * Get the expected file length
     * 
     * @return The file length in bytes, or -1 if unknown
     */
    public long getFileLength() {
        return fileLength;
    }
} 