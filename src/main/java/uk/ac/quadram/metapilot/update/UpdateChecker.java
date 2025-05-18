package uk.ac.quadram.metapilot.update;

import java.io.BufferedReader;
import java.io.File;
import java.io.InputStreamReader;
import java.net.HttpURLConnection;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.logging.Level;
import java.util.logging.Logger;

import javafx.application.Platform;
import javafx.scene.control.Alert;
import javafx.scene.control.ButtonType;
import javafx.scene.control.Hyperlink;
import javafx.scene.control.Label;
import javafx.scene.control.Alert.AlertType;
import javafx.scene.layout.VBox;
import javafx.application.HostServices;

/**
 * Checks for updates to the MetaPilot application
 */
public class UpdateChecker {
    private static final Logger LOGGER = Logger.getLogger(UpdateChecker.class.getName());
    
    // Current version
    private static final String CURRENT_VERSION = "1.0.0";
    
    // URL to check for updates - can be changed to your server when available
    // For now, we'll use a file hosting service (example URL)
    private static final String UPDATE_CHECK_URL = "https://example.com/metapilot/version.txt";
    
    // Download URL for the latest version
    private static final String DOWNLOAD_URL = "https://example.com/metapilot/download";
    
    private HostServices hostServices;
    
    public UpdateChecker(HostServices hostServices) {
        this.hostServices = hostServices;
    }
    
    /**
     * Check for updates and return latest version if available
     */
    public void checkForUpdates(boolean showNoUpdatesMessage) {
        new Thread(() -> {
            try {
                // First try to check online
                String latestVersion = getLatestVersionFromServer();
                
                // If server check fails, try to check local file
                if (latestVersion == null) {
                    latestVersion = getLatestVersionFromFile();
                }
                
                // If still null, we couldn't check for updates
                if (latestVersion == null) {
                    if (showNoUpdatesMessage) {
                        showUpdateError();
                    }
                    return;
                }
                
                // Compare versions
                if (isNewerVersion(latestVersion)) {
                    showUpdateAvailable(latestVersion);
                } else if (showNoUpdatesMessage) {
                    showNoUpdatesAvailable();
                }
            } catch (Exception e) {
                LOGGER.log(Level.WARNING, "Error checking for updates", e);
                if (showNoUpdatesMessage) {
                    showUpdateError();
                }
            }
        }).start();
    }
    
    /**
     * Try to get the latest version from the server
     */
    private String getLatestVersionFromServer() {
        try {
            URL url = new URL(UPDATE_CHECK_URL);
            HttpURLConnection connection = (HttpURLConnection) url.openConnection();
            connection.setRequestMethod("GET");
            connection.setConnectTimeout(5000);
            connection.setReadTimeout(5000);
            
            if (connection.getResponseCode() == HttpURLConnection.HTTP_OK) {
                try (BufferedReader reader = new BufferedReader(
                        new InputStreamReader(connection.getInputStream()))) {
                    String version = reader.readLine();
                    if (version != null && !version.trim().isEmpty()) {
                        return version.trim();
                    }
                }
            }
        } catch (Exception e) {
            LOGGER.log(Level.WARNING, "Failed to get version from server", e);
        }
        
        return null;
    }
    
    /**
     * Try to get the latest version from a local file
     */
    private String getLatestVersionFromFile() {
        try {
            File versionFile = new File("config/latest_version.txt");
            if (versionFile.exists()) {
                return Files.readString(Paths.get(versionFile.getPath())).trim();
            }
        } catch (Exception e) {
            LOGGER.log(Level.WARNING, "Failed to get version from local file", e);
        }
        
        return null;
    }
    
    /**
     * Check if the latest version is newer than the current version
     */
    private boolean isNewerVersion(String latestVersion) {
        try {
            // Compare versions assuming format is x.y.z
            String[] currentParts = CURRENT_VERSION.split("\\.");
            String[] latestParts = latestVersion.split("\\.");
            
            for (int i = 0; i < Math.min(currentParts.length, latestParts.length); i++) {
                int current = Integer.parseInt(currentParts[i]);
                int latest = Integer.parseInt(latestParts[i]);
                
                if (latest > current) {
                    return true;
                } else if (latest < current) {
                    return false;
                }
            }
            
            // If we get here, the versions are equal up to the common parts
            // The longer version is considered newer
            return latestParts.length > currentParts.length;
        } catch (Exception e) {
            LOGGER.log(Level.WARNING, "Error comparing versions", e);
            return false;
        }
    }
    
    /**
     * Show a dialog that updates are available
     */
    private void showUpdateAvailable(String latestVersion) {
        Platform.runLater(() -> {
            Alert alert = new Alert(AlertType.INFORMATION);
            alert.setTitle("Update Available");
            alert.setHeaderText("A new version of MetaPilot is available!");
            
            VBox content = new VBox(10);
            Label infoLabel = new Label(
                "Your version: " + CURRENT_VERSION + "\n" +
                "Latest version: " + latestVersion + "\n\n" +
                "Would you like to download the new version?"
            );
            
            Hyperlink downloadLink = new Hyperlink("Download Update");
            downloadLink.setOnAction(e -> {
                if (hostServices != null) {
                    hostServices.showDocument(DOWNLOAD_URL);
                }
                alert.close();
            });
            
            content.getChildren().addAll(infoLabel, downloadLink);
            alert.getDialogPane().setContent(content);
            
            alert.getButtonTypes().setAll(ButtonType.CLOSE);
            alert.showAndWait();
        });
    }
    
    /**
     * Show a dialog that no updates are available
     */
    private void showNoUpdatesAvailable() {
        Platform.runLater(() -> {
            Alert alert = new Alert(AlertType.INFORMATION);
            alert.setTitle("No Updates Available");
            alert.setHeaderText("Your version of MetaPilot is up to date");
            alert.setContentText("You are running the latest version (" + CURRENT_VERSION + ").");
            alert.showAndWait();
        });
    }
    
    /**
     * Show an error dialog when update check fails
     */
    private void showUpdateError() {
        Platform.runLater(() -> {
            Alert alert = new Alert(AlertType.WARNING);
            alert.setTitle("Update Check Failed");
            alert.setHeaderText("Could not check for updates");
            alert.setContentText(
                "Unable to connect to the update server.\n" +
                "Please check your internet connection or try again later.\n\n" +
                "Current version: " + CURRENT_VERSION
            );
            alert.showAndWait();
        });
    }
    
    /**
     * Get the current version
     */
    public static String getCurrentVersion() {
        return CURRENT_VERSION;
    }
} 