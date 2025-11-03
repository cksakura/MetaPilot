package uk.ac.quadram.metapilot.core.MGYG;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Date;
import java.util.HashMap;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

import javafx.concurrent.Task;

import java.text.SimpleDateFormat;

import org.apache.commons.net.ftp.FTPClient;
import org.apache.commons.net.ftp.FTPFile;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import uk.ac.quadram.metapilot.core.tools.ResumableDownloadTaskJavaFX;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;

/**
 * JavaFX-Compatible MGnify Genome Database FTP Download and Management System
 * 
 * <h2>Purpose</h2>
 * This JavaFX-compatible version provides the same comprehensive functionality as MGYGFtpTask
 * for downloading, updating, and managing MGnify (formerly EBI Metagenomics) genome databases.
 * It extends javafx.concurrent.Task for seamless integration with JavaFX applications and
 * progress monitoring systems.
 * 
 * <h2>Key Differences from Original MGYGFtpTask</h2>
 * <ul>
 *   <li><strong>Threading Model:</strong> Uses JavaFX Task instead of SwingWorker</li>
 *   <li><strong>Progress Reporting:</strong> Native JavaFX progress property binding</li>
 *   <li><strong>UI Thread Safety:</strong> Automatic Platform.runLater() for UI updates</li>
 *   <li><strong>Property Binding:</strong> Observable properties for real-time UI updates</li>
 * </ul>
 * 
 * <h2>JavaFX Integration Features</h2>
 * <ul>
 *   <li><strong>Progress Property:</strong> Bindable to JavaFX ProgressBar</li>
 *   <li><strong>Message Property:</strong> Real-time status updates</li>
 *   <li><strong>Value Property:</strong> Task result as Observable</li>
 *   <li><strong>State Property:</strong> Task state monitoring (READY, RUNNING, SUCCEEDED, etc.)</li>
 * </ul>
 * 
 * <h2>Usage with JavaFX Components</h2>
 * <pre>{@code
 *   // Create task
 *   MGYGFtpTaskJavaFX task = new MGYGFtpTaskJavaFX();
 *   
 *   // Bind to JavaFX ProgressBar
 *   progressBar.progressProperty().bind(task.progressProperty());
 *   
 *   // Bind to status label
 *   statusLabel.textProperty().bind(task.messageProperty());
 *   
 *   // Execute on background thread
 *   Thread thread = new Thread(task);
 *   thread.setDaemon(true);
 *   thread.start();
 * }</pre>
 * 
 * <h2>Maintained Functionality</h2>
 * <p>All core functionality from the original MGYGFtpTask is preserved:</p>
 * <ul>
 *   <li>Repository discovery and catalog enumeration</li>
 *   <li>Version detection and metadata extraction</li>
 *   <li>Full catalog and selective downloads</li>
 *   <li>Resumable and multi-threaded downloads</li>
 *   <li>Version comparison and incremental updates</li>
 *   <li>Error handling and network resilience</li>
 * </ul>
 * 
 * @author Kai Cheng - Quadram Institute Bioscience
 * @version 2.0
 * @since 2.0
 * 
 * @see uk.ac.quadram.metapilot.core.MGYG.MGYGFtpTask
 * @see uk.ac.quadram.metapilot.task.mag.MagDbItem
 * @see javafx.concurrent.Task
 */
public class MGYGFtpTaskJavaFX extends Task<Boolean> {

	private static final String server = "ftp.ebi.ac.uk";
	private static final int port = 21;
	private static final String user = "anonymous";
	private static final String pass = "";

	private static final String mgnify = "/pub/databases/metagenomics/mgnify_genomes";

	private MagDbItem item;
	private MagDbItem[] repositoryItems;
	
	private long missLength = 146;
	private int waitHours = 4;
	
	private static String taskName = "MGnify database FTP task";
	private static final Logger LOGGER = LogManager.getLogger(MGYGFtpTaskJavaFX.class);
	private SimpleDateFormat format = new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss");
	
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
		ArrayList<MagDbItem> list = new ArrayList<MagDbItem>();
		FTPClient ftpClient = new FTPClient();

		try {
			// connect and login to the server
			ftpClient.connect(server, port);
			ftpClient.login(user, pass);

			// use passive mode to pass firewall
			ftpClient.enterLocalPassiveMode();

			updateMessage("Connected to " + server);
			System.out.println("---Connected to " + server + "---");

			FTPFile[] mgnifyFiles = ftpClient.listFiles(mgnify);
			int finishCount = 0;
			int totalFiles = mgnifyFiles.length;

			for (FTPFile mgnifyFile : mgnifyFiles) {
				if (isCancelled()) {
					break;
				}
				
				if (mgnifyFile.isDirectory()) {
					try {
						String catalogueID = mgnifyFile.getName();
						String path = mgnify + "/" + catalogueID;

						updateMessage("Processing catalog: " + catalogueID);

						MagDbItem magDbItem = new MagDbItem(catalogueID, "MGYG");

						FTPFile[] files = ftpClient.listFiles(path);

						String version = "";
						Date date = null;
						if (files.length == 1) {
							Calendar calendar = files[0].getTimestamp();
							date = calendar.getTime();
							version = files[0].getName();
						} else {
							Calendar calendar = files[0].getTimestamp();
							version = files[0].getName();
							date = calendar.getTime();

							for (int j = 1; j < files.length; j++) {
								if (files[j].getTimestamp().after(calendar)) {
									version = files[j].getName();
									calendar = files[j].getTimestamp();
									date = calendar.getTime();
								}
							}
						}
						ftpClient.changeWorkingDirectory(path);

						FTPFile[] genomeFile = ftpClient.listFiles(path + "/" + version);
						int speciesCount = 0;

						for (int j = 0; j < genomeFile.length; j++) {
							String nameJ = genomeFile[j].getName();
							if (nameJ.equals("species_catalogue")) {

								String pathString = path + "/" + version + "/" + nameJ;
								FTPFile[] speciesFiles = ftpClient.listFiles(pathString);

								for (FTPFile speciesFile : speciesFiles) {

									String parent = speciesFile.getName();
									FTPFile[] subSpeciesFiles = ftpClient.listFiles(pathString + "/" + parent);
									speciesCount += subSpeciesFiles.length;

									for (FTPFile subSpeciesFile : subSpeciesFiles) {
										String child = subSpeciesFile.getName();
										magDbItem.addMgygName(parent, child);
									}
								}
							}
						}

						magDbItem.addVersion(version, speciesCount, date, false);
						list.add(magDbItem);

						finishCount++;

						// Update progress
						double progress = ((double) finishCount / (double) totalFiles) * 100.0;
						updateProgress(progress, 100);

						String message = catalogueID + " detected, " + speciesCount + " species found";
						updateMessage(message);
						System.out.println(format.format(new Date()) + "\t" + message);
						
					} catch (IOException e) {
						LOGGER.error("Error processing catalog: " + mgnifyFile.getName(), e);
						e.printStackTrace();
					}
				}
			}

			ftpClient.logout();
			ftpClient.disconnect();

			updateMessage("Disconnected from server");
			System.out.println("---Disconnected---");

		} catch (IOException ex) {
			LOGGER.error("FTP connection error", ex);
			ex.printStackTrace();
			updateMessage("Connection error: " + ex.getMessage());
		}

		this.repositoryItems = list.toArray(new MagDbItem[list.size()]);
	}
	
	/**
	 * Compare remote repositories with local items
	 * 
	 * @param localItems Local repository items to compare
	 */
	private void compareRepositories(MagDbItem[] localItems) {
		ArrayList<MagDbItem> list = new ArrayList<MagDbItem>();
		FTPClient ftpClient = new FTPClient();

		try {
			// connect and login to the server
			ftpClient.connect(server, port);
			ftpClient.login(user, pass);

			// use passive mode to pass firewall
			ftpClient.enterLocalPassiveMode();

			updateMessage("Connected to " + server);
			System.out.println("---Connected to " + server + "---");

			FTPFile[] mgnifyFiles = ftpClient.listFiles(mgnify);
			int totalFiles = mgnifyFiles.length;
			int processedFiles = 0;

			L: for (FTPFile mgnifyFile : mgnifyFiles) {
				if (isCancelled()) {
					break;
				}
				
				if (mgnifyFile.isDirectory()) {
					try {
						String catalogueId = mgnifyFile.getName();
						String path = mgnify + "/" + catalogueId;
						FTPFile[] files = ftpClient.listFiles(path);

						updateMessage("Checking catalog: " + catalogueId);

						String version = "";
						Date date = null;
						if (files.length == 1) {
							Calendar calendar = files[0].getTimestamp();
							version = files[0].getName();
							date = calendar.getTime();
						} else {
							Calendar calendar = files[0].getTimestamp();
							version = files[0].getName();
							date = calendar.getTime();

							for (int j = 1; j < files.length; j++) {
								if (files[j].getTimestamp().after(calendar)) {
									version = files[j].getName();
									calendar = files[j].getTimestamp();
									date = calendar.getTime();
								}
							}
						}

						MagDbItem magDbItem = null;
						for (int i = 0; i < localItems.length; i++) {
							if (localItems[i].getCatalogueID().equals(catalogueId)) {
								if (localItems[i].getCatalogueVersion().equals(version)) {
									list.add(localItems[i]);
									String message = localItems[i].getCatalogueID() + " " + 
													localItems[i].getCatalogueVersion() + " already exists";
									updateMessage(message);
									System.out.println(format.format(new Date()) + "\t" + message);
									
									processedFiles++;
									updateProgress((double) processedFiles / totalFiles * 100, 100);
									continue L;
								} else {
									magDbItem = localItems[i];
								}
							}
						}
						if (magDbItem == null) {
							magDbItem = new MagDbItem(catalogueId, "MGYG");
						} else {
							magDbItem.resetMgygMap();
						}
						
						String message = "Parsing new MAG database " + catalogueId;
						updateMessage(message);
						System.out.println(format.format(new Date()) + "\t" + message);

						ftpClient.changeWorkingDirectory(path);

						FTPFile[] genomeFile = ftpClient.listFiles(path + "/" + version);
						int speciesCount = 0;

						for (int j = 0; j < genomeFile.length; j++) {
							String nameJ = genomeFile[j].getName();
							if (nameJ.equals("species_catalogue")) {

								String pathString = path + "/" + version + "/" + nameJ;
								FTPFile[] speciesFiles = ftpClient.listFiles(pathString);

								for (FTPFile speciesFile : speciesFiles) {

									String parent = speciesFile.getName();
									FTPFile[] subSpeciesFiles = ftpClient.listFiles(pathString + "/" + parent);
									speciesCount += subSpeciesFiles.length;

									for (FTPFile subSpeciesFile : subSpeciesFiles) {
										String child = subSpeciesFile.getName();
										magDbItem.addMgygName(parent, child);
									}
								}
							}
						}

						magDbItem.addVersion(version, speciesCount, date, true);
						list.add(magDbItem);

						processedFiles++;
						updateProgress((double) processedFiles / totalFiles * 100, 100);

					} catch (IOException e) {
						LOGGER.error("Error comparing catalog: " + mgnifyFile.getName(), e);
						e.printStackTrace();
					}
				}
			}

			ftpClient.logout();
			ftpClient.disconnect();

			updateMessage("Comparison completed");
			System.out.println("---Disconnected---");

		} catch (IOException ex) {
			LOGGER.error("FTP connection error during comparison", ex);
			ex.printStackTrace();
			updateMessage("Connection error: " + ex.getMessage());
		}

		this.repositoryItems = list.toArray(new MagDbItem[list.size()]);
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
		HashMap<String, ArrayList<String>> map = item.getMgygNameMap();
		File dbFile = new File(item.getRepositoryFile(), item.getCatalogueVersion());
		
		updateMessage("Setting up download directories for " + catalogueId);
		
		if (!dbFile.exists()) {
			dbFile.mkdirs();
		}

		File originalDbFile = new File(dbFile, "original_db");
		if (!originalDbFile.exists()) {
			originalDbFile.mkdir();
		}

		File eggnogDbFile = new File(dbFile, "eggNOG");
		if (!eggnogDbFile.exists()) {
			eggnogDbFile.mkdir();
		}

		String baseLink = "http://" + server + "/" + mgnify + "/" + catalogueId + "/" + version;

		ExecutorService executor = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors() / 2);
		ArrayList<Future<?>> flist = new ArrayList<Future<?>>();

		updateMessage("Downloading metadata file...");
		updateProgress(5, 100);

		// Download metadata file
		File localMetaFile = new File(dbFile, "genomes-all_metadata.tsv");
		if (!localMetaFile.exists() || localMetaFile.length() <= missLength) {
			ResumableDownloadTaskJavaFX task = new ResumableDownloadTaskJavaFX(baseLink + "/genomes-all_metadata.tsv",
					localMetaFile.getAbsolutePath(), false);
			Future<?> f = executor.submit(task);
			flist.add(f);
		}

		updateMessage("Preparing genome file downloads...");
		updateProgress(10, 100);

		// Download genome files
		for (String parent : map.keySet()) {
			if (isCancelled()) {
				executor.shutdownNow();
				updateMessage("Download cancelled by user");
				return false;
			}
			
			ArrayList<String> list = map.get(parent);
			for (String child : list) {
				if (isCancelled()) {
					executor.shutdownNow();
					updateMessage("Download cancelled by user");
					return false;
				}
				
				String faaLink = baseLink + "/species_catalogue/" + parent + "/" + child + "/genome/" + child + ".faa";
				String eggnogLink = baseLink + "/species_catalogue/" + parent + "/" + child + "/genome/" + child
						+ "_eggNOG.tsv";

				File faaFile = new File(originalDbFile, child + ".faa");
				File eggnogFile = new File(eggnogDbFile, child + "_eggNOG.tsv");

				if (!faaFile.exists() || faaFile.length() <= missLength) {
					ResumableDownloadTaskJavaFX task = new ResumableDownloadTaskJavaFX(faaLink, faaFile.getAbsolutePath(), false);
					Future<?> f = executor.submit(task);
					flist.add(f);
				}

				if (!eggnogFile.exists() || eggnogFile.length() <= missLength) {
					ResumableDownloadTaskJavaFX task = new ResumableDownloadTaskJavaFX(eggnogLink, eggnogFile.getAbsolutePath(),
							false);
					Future<?> f = executor.submit(task);
					flist.add(f);
				}
			}
		}

		updateMessage("Downloading " + flist.size() + " files...");
		updateProgress(15, 100);

		// Monitor download progress
		for (int i = 0; i < flist.size(); i++) {
			if (isCancelled()) {
				executor.shutdownNow();
				updateMessage("Download cancelled by user");
				return false;
			}
			
			Future<?> f = flist.get(i);
			try {
				f.get();
			} catch (InterruptedException | ExecutionException e) {
				LOGGER.error("Error downloading file", e);
				e.printStackTrace();
			}

			// Update progress from 15% to 85% during downloads
			double progress = 15 + (i * 70.0 / (double) flist.size());
			updateProgress(progress, 100);
			updateMessage("Downloaded " + (i + 1) + " of " + flist.size() + " files");
		}

		boolean finish = false;

		try {
			updateMessage("Waiting for all downloads to complete...");
			executor.shutdown();

			finish = executor.awaitTermination(waitHours, TimeUnit.HOURS);

			updateProgress(90, 100);

			if (finish) {
				if (isCancelled()) {
					updateMessage("Download cancelled by user");
					return false;
				}
				
				updateMessage("Extracting high abundant proteins...");
				LOGGER.info(taskName + ": extracting high abundant proteins...");
				System.out.println(
						format.format(new Date()) + "\t" + taskName + ": extracting high abundant proteins...");

				File hapFile = new File(dbFile, "rib_elon.fasta");
				HapExtractor extractor = new HapExtractor(true);
				extractor.extract(originalDbFile, hapFile);
				
				updateProgress(95, 100);
				
				item.setAvailable();

				updateMessage("Download completed successfully!");
				LOGGER.info(taskName + ": finished");
				System.out.println(format.format(new Date()) + "\t" + taskName + ": finished");

				updateProgress(100, 100);
				finish = true;

			} else {
				String errorMsg = "Task didn't finish in " + waitHours + " hours, please try again";
				updateMessage(errorMsg);
				LOGGER.info(taskName + ": " + errorMsg);
				System.out.println(format.format(new Date()) + "\t" + taskName + ": " + errorMsg);
			}

		} catch (InterruptedException | IOException e) {
			String errorMsg = "Download failed: " + e.getMessage();
			updateMessage(errorMsg);
			LOGGER.error(taskName + ": failed", e);
			System.out.println(format.format(new Date()) + "\t" + taskName + ": failed");
		}

		return finish;
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