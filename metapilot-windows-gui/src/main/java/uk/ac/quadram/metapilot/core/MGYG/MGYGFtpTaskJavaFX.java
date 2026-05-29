package uk.ac.quadram.metapilot.core.MGYG;

import java.io.File;
import java.io.IOException;
import javafx.concurrent.Task;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;

/**
 * MGnify Genome Database FTP Download and Management System (JavaFX).
 * 
 * <h2>Purpose</h2>
 * Downloads, updates, and manages MGnify (formerly EBI Metagenomics) genome databases.
 * Extends javafx.concurrent.Task for integration with JavaFX progress and UI.
 * 
 * <h2>JavaFX Integration</h2>
 * <ul>
 *   <li><strong>Progress Property:</strong> Bindable to JavaFX ProgressBar</li>
 *   <li><strong>Message Property:</strong> Real-time status updates</li>
 *   <li><strong>Value Property:</strong> Task result as Observable</li>
 *   <li><strong>State Property:</strong> Task state (READY, RUNNING, SUCCEEDED, etc.)</li>
 * </ul>
 * 
 * <h2>Usage</h2>
 * <pre>{@code
 *   MGYGFtpTaskJavaFX task = new MGYGFtpTaskJavaFX();
 *   progressBar.progressProperty().bind(task.progressProperty());
 *   statusLabel.textProperty().bind(task.messageProperty());
 *   new Thread(task).setDaemon(true).start();
 * }</pre>
 * 
 * <h2>Functionality</h2>
 * <ul>
 *   <li>Repository discovery and catalog enumeration</li>
 *   <li>Version detection and metadata extraction</li>
 *   <li>Full catalog and selective downloads</li>
 *   <li>Resumable and multi-threaded downloads</li>
 *   <li>Version comparison and incremental updates</li>
 * </ul>
 * 
 * @see uk.ac.quadram.metapilot.task.mag.MagDbItem
 * @see javafx.concurrent.Task
 */
public class MGYGFtpTaskJavaFX extends Task<Boolean> {

	private MagDbItem item;
	private MagDbItem[] repositoryItems;

	private int waitHours = 4;
	
	private static String taskName = "MGnify database FTP task";
	private static final Logger LOGGER = LogManager.getLogger(MGYGFtpTaskJavaFX.class);
	
	/**
	 * Default constructor for repository discovery
	 */
	public MGYGFtpTaskJavaFX() {
	}
	
	/**
	 * Constructor for downloading a specific catalog item
	 * 
	 * @param item The catalog item to download
	 */
	public MGYGFtpTaskJavaFX(MagDbItem item) {
		this.item = item;
	}
	
	/**
	 * Constructor for comparing repositories with local items
	 * 
	 * @param repositoryItems Local repository items to compare
	 */
	public MGYGFtpTaskJavaFX(MagDbItem[] repositoryItems) {
		this.repositoryItems = repositoryItems;
	}

	/**
	 * Get the discovered repository items
	 * 
	 * @return Array of discovered MagDbItem objects
	 */
	public MagDbItem[] getRepositoryItems() {
		return repositoryItems;
	}

	@Override
	protected Boolean call() throws Exception {
		boolean finish = false;
		
		updateMessage("Initializing MGnify FTP connection...");
		updateProgress(0, 100);
		
		if (item == null) {
			if (repositoryItems == null || repositoryItems.length == 0) {
				updateMessage("Discovering MGnify repositories...");
				getRepositories();
			} else {
				updateMessage("Comparing with local repositories...");
				compareRepositories(repositoryItems);
			}

			if (this.repositoryItems != null) {
				finish = true;
			}
		} else {
			updateMessage("Downloading catalog: " + item.getCatalogueID());
			finish = download();
		}

		if (finish) {
			updateMessage("Task completed successfully");
			updateProgress(100, 100);
		} else {
			updateMessage("Task failed");
		}

		return finish;
	}
	
	/**
	 * Discover all available MGnify repositories
	 */
	private void getRepositories() {
		try {
			ProgressCallback cb = new ProgressCallback() {
				@Override
				public void onProgress(double workDone, double max) {
					if (isCancelled()) {
						Thread.currentThread().interrupt();
						return;
					}
					updateProgress(workDone, max);
				}

				@Override
				public void onMessage(String message) {
					if (isCancelled()) {
						Thread.currentThread().interrupt();
						return;
					}
					updateMessage(message);
				}

				@Override
				public void onComplete(Object result) {
					// no-op
				}

				@Override
				public void onError(Exception error) {
					if (error != null) {
						updateMessage("Connection error: " + error.getMessage());
					}
				}
			};

			this.repositoryItems = MGYGRepositoryService.discoverRepositories(cb);

		} catch (IOException ex) {
			LOGGER.error("FTP connection error", ex);
			updateMessage("Connection error: " + ex.getMessage());
		}
	}
	
	/**
	 * Compare remote repositories with local items
	 * 
	 * @param localItems Local repository items to compare
	 */
	private void compareRepositories(MagDbItem[] localItems) {
		try {
			ProgressCallback cb = new ProgressCallback() {
				@Override
				public void onProgress(double workDone, double max) {
					if (isCancelled()) {
						Thread.currentThread().interrupt();
						return;
					}
					updateProgress(workDone, max);
				}

				@Override
				public void onMessage(String message) {
					if (isCancelled()) {
						Thread.currentThread().interrupt();
						return;
					}
					updateMessage(message);
				}

				@Override
				public void onComplete(Object result) {
					// no-op
				}

				@Override
				public void onError(Exception error) {
					if (error != null) {
						updateMessage("Connection error: " + error.getMessage());
					}
				}
			};

			this.repositoryItems = MGYGRepositoryService.compareRepositories(localItems, cb);

		} catch (IOException ex) {
			LOGGER.error("FTP connection error during comparison", ex);
			updateMessage("Connection error: " + ex.getMessage());
		}
	}
	
	/**
	 * Download the specified catalog item with complete FTP functionality
	 * 
	 * @return true if download was successful
	 */
	private boolean download() {
		if (item == null) {
			updateMessage("No item specified for download");
			return false;
		}

		String catalogueId = item.getCatalogueID();
		String version = item.getCatalogueVersion();

		// The GUI stores repositoryFile as the catalogue folder (…/<catalogueId>).
		// The core downloader writes to <base>/<catalogueId>/<version>, so we pass the parent directory as base.
		File catalogueDir = item.getRepositoryFile();
		if (catalogueDir == null) {
			updateMessage("ERROR: Local repository is not set for catalogue: " + catalogueId);
			return false;
		}
		File baseDir = catalogueDir.getParentFile();
		if (baseDir == null) {
			updateMessage("ERROR: Invalid local repository: " + catalogueDir.getAbsolutePath());
			return false;
		}

		updateMessage("Starting MGYG download (headless core implementation)...");
		updateProgress(0, 100);

		final long[] lastProgressUpdateMs = new long[] { 0L };
		ProgressCallback cb = new ProgressCallback() {
			@Override
			public void onProgress(double workDone, double max) {
				if (isCancelled()) {
					Thread.currentThread().interrupt();
					return;
				}
				// Throttle JavaFX updates
				long now = System.currentTimeMillis();
				if (now - lastProgressUpdateMs[0] < 250) {
					return;
				}
				lastProgressUpdateMs[0] = now;

				double pct = 0.0;
				if (max > 0) {
					pct = (workDone / max) * 100.0;
				}
				if (pct < 0) pct = 0;
				if (pct > 100) pct = 100;
				updateProgress(pct, 100);
			}

			@Override
			public void onMessage(String message) {
				if (isCancelled()) {
					Thread.currentThread().interrupt();
					return;
				}
				updateMessage(message);
			}

			@Override
			public void onComplete(Object result) {
				updateProgress(100, 100);
			}

			@Override
			public void onError(Exception error) {
				// Surface the error message via JavaFX message property; exception will still be thrown to Task
				if (error != null) {
					updateMessage("ERROR: " + error.getMessage());
				}
			}
		};

		try {
			// Resume is safer for large downloads and matches expected HPC behavior.
			boolean resume = true;
			int threads = Math.max(1, Runtime.getRuntime().availableProcessors() / 2);

			MGYGMagDownloader.downloadCatalogue(catalogueId, version, baseDir, resume, threads, waitHours, cb);

			if (isCancelled()) {
				updateMessage("Download cancelled by user");
				return false;
			}

			item.setAvailable();
			updateMessage("Download completed successfully!");
			updateProgress(100, 100);
			return true;

		} catch (IOException e) {
			if (isCancelled() || Thread.currentThread().isInterrupted()) {
				updateMessage("Download cancelled");
				return false;
			}
			String errorMsg = "Download failed: " + e.getMessage();
			updateMessage(errorMsg);
			LOGGER.error(taskName + ": failed", e);
			return false;
		}
	}
	
	/**
	 * Get the task name for logging
	 * 
	 * @return The task name
	 */
	public String getTaskName() {
		return taskName;
	}
	
	/**
	 * Get the logger instance
	 * 
	 * @return The logger
	 */
	public Logger getLogger() {
		return LOGGER;
	}
} 