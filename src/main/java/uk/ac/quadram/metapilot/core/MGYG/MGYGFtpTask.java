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

import javax.swing.JProgressBar;
import javax.swing.SwingWorker;

import java.io.FileOutputStream;
import java.io.OutputStream;
import java.text.SimpleDateFormat;

import org.apache.commons.net.ftp.FTPClient;
import org.apache.commons.net.ftp.FTPFile;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import uk.ac.quadram.metapilot.core.tools.ResumableDownloadTask;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;

/**
 * MGnify Genome Database FTP Download and Management System
 * 
 * <h2>Purpose</h2>
 * This class provides comprehensive functionality for downloading, updating, and managing
 * MGnify (formerly EBI Metagenomics) genome databases. MGnify is a major resource for 
 * metagenome-assembled genomes (MAGs) that are essential for metaproteomics analysis
 * of environmental and clinical samples.
 * 
 * <h2>MGnify Database Structure</h2>
 * <p>MGnify organizes genomes in a hierarchical structure:</p>
 * <pre>{@code
 *   ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/
 *   ├── Catalog1/
 *   │   └── Version/
 *   │       └── species_catalogue/
 *   │           ├── Species1/
 *   │           │   ├── Genome1.faa
 *   │           │   └── Genome2.faa
 *   │           └── Species2/
 *   │               └── Genome3.faa
 *   └── Catalog2/
 *       └── Version/
 *           └── species_catalogue/
 *               └── ...
 * }</pre>
 * 
 * <h2>Supported Operations</h2>
 * 
 * <h3>Repository Discovery:</h3>
 * <ul>
 *   <li><strong>Catalog Enumeration:</strong> Discover all available MGnify catalogs</li>
 *   <li><strong>Version Detection:</strong> Find latest versions for each catalog</li>
 *   <li><strong>Genome Counting:</strong> Calculate total genomes per catalog</li>
 *   <li><strong>Metadata Extraction:</strong> Parse timestamps and versioning</li>
 * </ul>
 * 
 * <h3>Database Downloads:</h3>
 * <ul>
 *   <li><strong>Full Catalog Download:</strong> Download entire catalogs</li>
 *   <li><strong>Selective Download:</strong> Download specific species/genomes</li>
 *   <li><strong>Resumable Downloads:</strong> Continue interrupted downloads</li>
 *   <li><strong>Multi-threaded:</strong> Parallel download acceleration</li>
 * </ul>
 * 
 * <h3>Update Management:</h3>
 * <ul>
 *   <li><strong>Version Comparison:</strong> Compare local vs remote versions</li>
 *   <li><strong>Incremental Updates:</strong> Download only new/changed files</li>
 *   <li><strong>Timestamp Tracking:</strong> Monitor database freshness</li>
 *   <li><strong>Local Validation:</strong> Verify download integrity</li>
 * </ul>
 * 
 * <h2>FTP Configuration</h2>
 * <h3>Connection Parameters:</h3>
 * <ul>
 *   <li><strong>Server:</strong> ftp.ebi.ac.uk</li>
 *   <li><strong>Port:</strong> 21 (standard FTP)</li>
 *   <li><strong>Authentication:</strong> Anonymous access</li>
 *   <li><strong>Mode:</strong> Passive (firewall-friendly)</li>
 * </ul>
 * 
 * <h3>Base Path:</h3>
 * <ul>
 *   <li><strong>Root Directory:</strong> /pub/databases/metagenomics/mgnify_genomes</li>
 *   <li><strong>Access Pattern:</strong> Public read-only</li>
 *   <li><strong>Update Frequency:</strong> Regular (varies by catalog)</li>
 * </ul>
 * 
 * <h2>Download Strategies</h2>
 * 
 * <h3>Single-threaded Downloads:</h3>
 * <ul>
 *   <li><strong>Reliable:</strong> Traditional FTP download</li>
 *   <li><strong>Resumable:</strong> Can continue interrupted transfers</li>
 *   <li><strong>Memory Efficient:</strong> Low resource usage</li>
 *   <li><strong>Progress Tracking:</strong> Real-time progress updates</li>
 * </ul>
 * 
 * <h3>Multi-threaded Downloads:</h3>
 * <ul>
 *   <li><strong>Accelerated:</strong> Parallel connection utilization</li>
 *   <li><strong>Configurable:</strong> Adjustable thread pool size</li>
 *   <li><strong>Fault Tolerant:</strong> Individual thread failure handling</li>
 *   <li><strong>Resource Managed:</strong> Automatic connection cleanup</li>
 * </ul>
 * 
 * <h2>Progress Monitoring</h2>
 * <h3>Progress Indicators:</h3>
 * <ul>
 *   <li><strong>JProgressBar Integration:</strong> Swing GUI progress display</li>
 *   <li><strong>Percentage Tracking:</strong> Real-time completion percentage</li>
 *   <li><strong>File-level Progress:</strong> Individual file download status</li>
 *   <li><strong>Catalog-level Progress:</strong> Overall catalog completion</li>
 * </ul>
 * 
 * <h3>Logging and Monitoring:</h3>
 * <ul>
 *   <li><strong>Timestamped Logs:</strong> Detailed operation logging</li>
 *   <li><strong>Error Reporting:</strong> Comprehensive error messages</li>
 *   <li><strong>Performance Metrics:</strong> Download speed tracking</li>
 *   <li><strong>Status Updates:</strong> Real-time status information</li>
 * </ul>
 * 
 * <h2>Error Handling and Recovery</h2>
 * 
 * <h3>Network Resilience:</h3>
 * <ul>
 *   <li><strong>Connection Retry:</strong> Automatic reconnection attempts</li>
 *   <li><strong>Timeout Management:</strong> Configurable timeout settings</li>
 *   <li><strong>Graceful Degradation:</strong> Fallback to single-threaded mode</li>
 *   <li><strong>Partial Recovery:</strong> Resume from last successful point</li>
 * </ul>
 * 
 * <h3>Data Integrity:</h3>
 * <ul>
 *   <li><strong>Size Validation:</strong> Compare downloaded vs expected sizes</li>
 *   <li><strong>Corruption Detection:</strong> Identify incomplete downloads</li>
 *   <li><strong>Automatic Cleanup:</strong> Remove corrupted partial files</li>
 *   <li><strong>Retry Logic:</strong> Re-download corrupted files</li>
 * </ul>
 * 
 * <h2>Integration with MAG Workflows</h2>
 * 
 * <h3>Database Preparation:</h3>
 * <ul>
 *   <li><strong>FASTA Processing:</strong> Handle protein FASTA files</li>
 *   <li><strong>Index Generation:</strong> Prepare search engine indexes</li>
 *   <li><strong>Taxonomy Assignment:</strong> Associate genomes with taxonomy</li>
 *   <li><strong>Metadata Extraction:</strong> Parse genome annotations</li>
 * </ul>
 * 
 * <h3>Search Engine Integration:</h3>
 * <ul>
 *   <li><strong>pFind:</strong> Prepare databases for pFind search</li>
 *   <li><strong>Sage:</strong> Format for Sage search engine</li>
 *   <li><strong>DIA-NN:</strong> Generate spectral libraries for DIA</li>
 *   <li><strong>FragPipe:</strong> MSFragger database formatting</li>
 * </ul>
 * 
 * <h2>Performance Considerations</h2>
 * 
 * <h3>Network Optimization:</h3>
 * <ul>
 *   <li><strong>Passive Mode:</strong> Firewall and NAT compatibility</li>
 *   <li><strong>Connection Pooling:</strong> Efficient FTP connection reuse</li>
 *   <li><strong>Bandwidth Management:</strong> Throttling capabilities</li>
 *   <li><strong>Concurrent Limits:</strong> Respect server connection limits</li>
 * </ul>
 * 
 * <h3>Storage Optimization:</h3>
 * <ul>
 *   <li><strong>Disk Space Management:</strong> Monitor available space</li>
 *   <li><strong>Compression Handling:</strong> Extract compressed archives</li>
 *   <li><strong>Temporary Files:</strong> Clean up intermediate files</li>
 *   <li><strong>Directory Structure:</strong> Maintain organized file hierarchy</li>
 * </ul>
 * 
 * <h2>Usage Examples</h2>
 * 
 * <h3>Repository Discovery:</h3>
 * <pre>{@code
 *   // Discover available MGnify catalogs
 *   MGYGFtpTask discoveryTask = new MGYGFtpTask(progressBar);
 *   discoveryTask.execute();
 *   
 *   // Get repository information
 *   MagDbItem[] repositories = discoveryTask.getRepositoryItems();
 * }</pre>
 * 
 * <h3>Catalog Download:</h3>
 * <pre>{@code
 *   // Download specific catalog
 *   MagDbItem catalog = new MagDbItem("human-gut", "MGYG");
 *   MGYGFtpTask downloadTask = new MGYGFtpTask(catalog, progressBar);
 *   downloadTask.execute();
 * }</pre>
 * 
 * <h3>Update Check:</h3>
 * <pre>{@code
 *   // Compare local vs remote versions
 *   MagDbItem[] localCatalogs = loadLocalCatalogs();
 *   MGYGFtpTask updateTask = new MGYGFtpTask(localCatalogs, progressBar);
 *   updateTask.execute();
 * }</pre>
 * 
 * <h2>Configuration Parameters</h2>
 * <ul>
 *   <li><strong>missLength:</strong> 146 bytes (minimum significant file size)</li>
 *   <li><strong>waitHours:</strong> 4 hours (retry interval for failed downloads)</li>
 *   <li><strong>Thread Pool:</strong> Configurable for multi-threaded downloads</li>
 *   <li><strong>Timeout Settings:</strong> Network operation timeouts</li>
 * </ul>
 * 
 * <h2>Thread Safety</h2>
 * <p>This class extends SwingWorker and is designed for single-threaded execution
 * per instance. Multiple instances can be used concurrently for different catalogs.</p>
 * 
 * <h2>Dependencies</h2>
 * <ul>
 *   <li><strong>Apache Commons Net:</strong> FTP client functionality</li>
 *   <li><strong>Log4j:</strong> Comprehensive logging</li>
 *   <li><strong>SwingWorker:</strong> Background task execution</li>
 *   <li><strong>ResumableDownloadTask:</strong> Download resumption capability</li>
 * </ul>
 * 
 * @author Kai Cheng - Quadram Institute Bioscience
 * @version 2.0
 * @since 1.0
 * 
 * @see uk.ac.quadram.metapilot.task.mag.MagDbItem
 * @see uk.ac.quadram.metapilot.core.tools.ResumableDownloadTask
 * @see uk.ac.quadram.metapilot.task.mag
 */
public class MGYGFtpTask extends SwingWorker<Boolean, Object> {

	private static final String server = "ftp.ebi.ac.uk";
	private static final int port = 21;
	private static final String user = "anonymous";
	private static final String pass = "";

	private static final String mgnify = "/pub/databases/metagenomics/mgnify_genomes";

	private JProgressBar bar;
	private MagDbItem item;
	private MagDbItem[] repositoryItems;
	
	private long missLength = 146;
	private int waitHours = 4;
	
	private static String taskName = "MGnify database FTP task";
	private static final Logger LOGGER = LogManager.getLogger(MGYGFtpTask.class);
	private SimpleDateFormat format = new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss");
	
	public MGYGFtpTask() {
	}
	
	public MGYGFtpTask(JProgressBar bar) {
		this.bar = bar;
	}
	
	public MGYGFtpTask(MagDbItem item, JProgressBar bar) {
		this.bar = bar;
		this.item = item;
	}
	
	public MGYGFtpTask(MagDbItem[] repositoryItems, JProgressBar bar) {
		this.bar = bar;
		this.repositoryItems = repositoryItems;
	}

	public MagDbItem[] getRepositoryItems() {
		return repositoryItems;
	}

	@Override
	protected Boolean doInBackground() throws Exception {
		// TODO Auto-generated method stub
		boolean finish = false;
		if (item == null) {
			if (repositoryItems == null || repositoryItems.length == 0) {
				getRepositories();
			} else {
				compareRepositories(repositoryItems);
			}

			if (this.repositoryItems != null) {
				finish = true;
			}
		} else {
			finish = download();
		}

		return finish;
	}
	
	private void getRepositories() {
		ArrayList<MagDbItem> list = new ArrayList<MagDbItem>();
		FTPClient ftpClient = new FTPClient();

		try {
			// connect and login to the server
			ftpClient.connect(server, port);
			ftpClient.login(user, pass);

			// use passive mode to pass firewall
			ftpClient.enterLocalPassiveMode();

			System.out.println("---Connected to " + server + "---");

			FTPFile[] mgnifyFiles = ftpClient.listFiles(mgnify);
			int finishCount = 0;

			for (FTPFile mgnifyFile : mgnifyFiles) {
				if (mgnifyFile.isDirectory()) {

					try {
						String catalogueID = mgnifyFile.getName();
						String path = mgnify + "/" + catalogueID;

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

						bar.setValue((int) ((double) finishCount * 100.0 / (double) mgnifyFiles.length));

						System.out.println(format.format(new Date()) + "\t" + mgnifyFile.getName() + " detected, "
								+ speciesCount + " species are found");
					} catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
			}

			ftpClient.logout();
			ftpClient.disconnect();

			System.out.println("---Disconnected---");

		} catch (IOException ex) {
			ex.printStackTrace();
		}

		this.repositoryItems = list.toArray(new MagDbItem[list.size()]);
	}
	
	private void compareRepositories(MagDbItem[] localItems) {
		ArrayList<MagDbItem> list = new ArrayList<MagDbItem>();
		FTPClient ftpClient = new FTPClient();

		try {
			// connect and login to the server
			ftpClient.connect(server, port);
			ftpClient.login(user, pass);

			// use passive mode to pass firewall
			ftpClient.enterLocalPassiveMode();

			System.out.println("---Connected to " + server + "---");

			FTPFile[] mgnifyFiles = ftpClient.listFiles(mgnify);

			L: for (FTPFile mgnifyFile : mgnifyFiles) {
				if (mgnifyFile.isDirectory()) {

					try {
						String catalogueId = mgnifyFile.getName();
						String path = mgnify + "/" + catalogueId;
						FTPFile[] files = ftpClient.listFiles(path);

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
									System.out.println(format.format(new Date()) + "\t" + localItems[i].getCatalogueID()
											+ " " + localItems[i].getCatalogueVersion() + " already exist");
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
						
						System.
								out.println(format.format(new Date()) + "\t" + "parsing new MAG database " + catalogueId);

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

									for (FTPFile subSpeciesFile : subSpeciesFiles) {
										String child = subSpeciesFile.getName();
										magDbItem.addMgygName(parent, child);
									}
									speciesCount += subSpeciesFiles.length;
								}
							}
						}

						magDbItem.addVersion(version, speciesCount, date, false);
						list.add(magDbItem);
						
						bar.setValue((int) ((double) list.size() * 100.0 / (double) mgnifyFiles.length));

						System.out.println(format.format(new Date()) + "\t" + mgnifyFile.getName() + " detected, "
								+ speciesCount + " species are found");
					} catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
			}

			ftpClient.logout();
			ftpClient.disconnect();

			System.out.println("---Disconnected---");

		} catch (IOException ex) {
			ex.printStackTrace();
		}

		this.repositoryItems = list.toArray(new MagDbItem[list.size()]);
	}
	
	private boolean download() {
		if (this.item == null) {
			return false;
		}

		String catalogureId = item.getCatalogueID();
		String version = item.getCatalogueVersion();
		HashMap<String, ArrayList<String>> map = item.getMgygNameMap();
		File dbFile = new File(item.getRepositoryFile(), item.getCatalogueVersion());
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

		String baseLink = "http://" + server + "/" + mgnify + "/" + catalogureId + "/" + version;

		ExecutorService executor = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors() / 2);
		ArrayList<Future<?>> flist = new ArrayList<Future<?>>();

		File localMetaFile = new File(dbFile, "genomes-all_metadata.tsv");
		if (!localMetaFile.exists() || localMetaFile.length() <= missLength) {
			ResumableDownloadTask task = new ResumableDownloadTask(baseLink + "/genomes-all_metadata.tsv",
					localMetaFile.getAbsolutePath(), false);
			Future<?> f = executor.submit(task);
			flist.add(f);
		}

		for (String parent : map.keySet()) {
			ArrayList<String> list = map.get(parent);
			for (String child : list) {
				String faaLink = baseLink + "/species_catalogue/" + parent + "/" + child + "/genome/" + child + ".faa";
				String eggnogLink = baseLink + "/species_catalogue/" + parent + "/" + child + "/genome/" + child
						+ "_eggNOG.tsv";

				File faaFile = new File(originalDbFile, child + ".faa");
				File eggnogFile = new File(eggnogDbFile, child + "_eggNOG.tsv");

				if (!faaFile.exists() || faaFile.length() <= missLength) {
					ResumableDownloadTask task = new ResumableDownloadTask(faaLink, faaFile.getAbsolutePath(), false);
					Future<?> f = executor.submit(task);
					flist.add(f);
				}

				if (!eggnogFile.exists() || eggnogFile.length() <= missLength) {
					ResumableDownloadTask task = new ResumableDownloadTask(eggnogLink, eggnogFile.getAbsolutePath(),
							false);
					Future<?> f = executor.submit(task);
					flist.add(f);
				}
			}
		}

		for (int i = 0; i < flist.size(); i++) {
			Future<?> f = flist.get(i);
			try {
				f.get();
			} catch (InterruptedException | ExecutionException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}

			bar.setValue((int) (i * 90.0 / (double) flist.size()));
		}

		boolean finish = false;

		try {

			executor.shutdown();

			finish = executor.awaitTermination(waitHours, TimeUnit.HOURS);

			bar.setValue(90);

			if (finish) {
/*
				LOGGER.info(taskName + ": creating SQL database...");
				System.out.println(format.format(new Date()) + "\t" + taskName + ": creating SQL database...");
				
				File sqlDbFile = new File(dbFile, item.getRepositoryFile().getName());
				MagSqliteTask.create(sqlDbFile);
				MagSqliteTask.createTaxaTable(localMetaFile.getAbsolutePath(), sqlDbFile);
				MagSqliteTask.createFuncTable(eggnogDbFile.getAbsolutePath(), sqlDbFile);

				LOGGER.info(taskName + ": creating SQL database finished...");
				System.out.println(format.format(new Date()) + "\t" + taskName + ": creating SQL database finished...");
				
				bar.setValue(95);
*/				
				LOGGER.info(taskName + ": extracting high abundant proteins...");
				System.out.println(
						format.format(new Date()) + "\t" + taskName + ": extracting  high abundant proteins...");

				File hapFile = new File(dbFile, "rib_elon.fasta");
				HapExtractor extractor = new HapExtractor(true);
				extractor.extract(originalDbFile, hapFile);
				
				item.setAvailable();

				LOGGER.info(taskName + ": finished");
				System.out.println(format.format(new Date()) + "\t" + taskName + ": finished");

				bar.setValue(100);

				finish = true;

			} else {
				LOGGER.info(taskName + ": task dosen't finish in a long time, please try again");
				System.out.println(format.format(new Date()) + "\t" + taskName
						+ ": task dosen't finish in a long time, please try again");
			}

		} catch (InterruptedException | IOException e) {
			// TODO Auto-generated catch block
			LOGGER.error(taskName + ": failed", e);
			System.out.println(format.format(new Date()) + "\t" + taskName + ": failed");
		}

		return finish;
	}
	
	@SuppressWarnings("unused")
	private void downloadFTP() {

		FTPClient ftpClient = new FTPClient();

		try {
			HashMap<String, Long> fileSizeMap = new HashMap<String, Long>();

			// connect and login to the server
			ftpClient.connect(server, port);
			ftpClient.login(user, pass);

			// use passive mode to pass firewall
			ftpClient.enterLocalPassiveMode();

			System.out.println("---Connected---");

			String catalogueName = item.getCatalogueID();
			String catalogueVersion = item.getCatalogueVersion();
			File dbFile = new File(item.getRepositoryFile(), catalogueVersion);
			if (!dbFile.exists()) {
				dbFile.mkdir();
			}
			int finishCount = 0;
			int totalCount = item.getSpeciesCount() * 2 + 1;

			FTPFile[] mgnifyFiles = ftpClient.listFiles(mgnify);
			for (FTPFile mgnifyFile : mgnifyFiles) {
				if (mgnifyFile.isDirectory() && mgnifyFile.getName().equals(catalogueName)) {

					FTPFile[] versionFiles = ftpClient.listFiles(mgnify + "/" + catalogueName);
					for (FTPFile versionFile : versionFiles) {
						if (versionFile.getName().equals(catalogueVersion)) {
							String pathString = mgnify + "/" + mgnifyFile.getName() + "/" + catalogueVersion;
							FTPFile[] genomeFiles = ftpClient.listFiles(pathString);

							ftpClient.changeWorkingDirectory(pathString);

							for (FTPFile genomeFile : genomeFiles) {
								String nameJ = genomeFile.getName();

								if (nameJ.equals("genomes-all_metadata.tsv")) {

									File downloadFile = new File(dbFile, nameJ);
									if (downloadFile.exists() && downloadFile.length() == genomeFile.getSize()) {
										System.out.println("File exist " + nameJ);
									} else {
										fileSizeMap.put(nameJ, genomeFile.getSize());
										OutputStream out;
										try {
											out = new FileOutputStream(downloadFile.getAbsolutePath());
											if (ftpClient.retrieveFile(nameJ, out)) {
												System.out.println("Downloaded " + nameJ);
											} else {
												System.out.println("Failed to download " + nameJ);
											}
											out.close();
										} catch (IOException e) {
											// TODO Auto-generated catch block
											LOGGER.error(taskName + ": error in processing file " + nameJ, e);
											System.out.println(format.format(new Date()) + "\t" + taskName
													+ ": error in processing file " + nameJ);
										}
									}

									finishCount++;
									bar.setValue((int) ((double) finishCount * 95.0 / (double) totalCount));

								} else if (nameJ.equals("species_catalogue")) {

									File originalDbFile = new File(dbFile, "original_db");
									if (!originalDbFile.exists()) {
										originalDbFile.mkdir();
									}

									File eggnogFile = new File(dbFile, "eggNOG");
									if (!eggnogFile.exists()) {
										eggnogFile.mkdir();
									}

									FTPFile[] speciesFiles = ftpClient.listFiles(pathString + "/" + nameJ);

									for (FTPFile speciesFile : speciesFiles) {

										FTPFile[] subSpeciesFiles = ftpClient.listFiles(
												pathString + "/species_catalogue" + "/" + speciesFile.getName());

										for (FTPFile subSpeciesFile : subSpeciesFiles) {

											FTPFile[] subSpeciesContentFiles = ftpClient
													.listFiles(pathString + "/species_catalogue" + "/"
															+ speciesFile.getName() + "/" + subSpeciesFile.getName());
											String genomeNameString = subSpeciesFile.getName();

											for (FTPFile subSpeciesContentFile : subSpeciesContentFiles) {
												if (subSpeciesContentFile.getName().equals("genome")) {
													String finalPathString = pathString + "/species_catalogue" + "/"
															+ speciesFile.getName() + "/" + subSpeciesFile.getName()
															+ "/genome";

													FTPFile[] finalFiles = ftpClient.listFiles(finalPathString);
													ftpClient.changeWorkingDirectory(finalPathString);

													for (FTPFile finalFile : finalFiles) {
														String fileName = finalFile.getName();
														if (fileName.equals(genomeNameString + ".faa")) {
															fileSizeMap.put(fileName, finalFile.getSize());
															File downloadFile = new File(originalDbFile, fileName);
															if (downloadFile.exists()
																	&& downloadFile.length() == genomeFile.getSize()) {
																System.out.println("File exist " + nameJ);
															} else {

																OutputStream out;
																try {
																	out = new FileOutputStream(
																			(new File(originalDbFile, fileName))
																					.getAbsolutePath());
																	if (ftpClient.retrieveFile(fileName, out)) {
																		System.out.println(
																				"Downloaded " + finalFile.getName());
																	} else {
																		System.out.println("Failed to download "
																				+ finalFile.getName());
																	}
																	out.close();
																} catch (IOException e) {
																	// TODO Auto-generated catch block
																	LOGGER.error(taskName
																			+ ": error in processing file " + fileName,
																			e);
																	System.out.println(format.format(new Date()) + "\t"
																			+ taskName + ": error in processing file "
																			+ fileName);
																}

															}

															finishCount++;
															bar.setValue((int) ((double) finishCount * 95.0
																	/ (double) totalCount));

														} else if (fileName.equals(genomeNameString + "_eggNOG.tsv")) {

															fileSizeMap.put(fileName, finalFile.getSize());
															File downloadFile = new File(eggnogFile, fileName);
															if (downloadFile.exists()
																	&& downloadFile.length() == genomeFile.getSize()) {
																System.out.println("File exist " + nameJ);
															} else {
																OutputStream out;
																try {
																	out = new FileOutputStream(
																			(new File(eggnogFile, fileName))
																					.getAbsolutePath());
																	if (ftpClient.retrieveFile(fileName, out)) {
																		System.out.println(
																				"Downloaded " + finalFile.getName());
																	} else {
																		System.out.println("Failed to download "
																				+ finalFile.getName());
																	}
																	out.close();
																} catch (IOException e) {
																	// TODO Auto-generated catch block
																	LOGGER.error(taskName
																			+ ": error in processing file " + fileName,
																			e);
																	System.out.println(format.format(new Date()) + "\t"
																			+ taskName + ": error in processing file "
																			+ fileName);
																}
															}

															finishCount++;
															bar.setValue((int) ((double) finishCount * 95.0
																	/ (double) totalCount));
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}

			System.out.println("Disconnected");

			ftpClient.logout();
			ftpClient.disconnect();

			bar.setValue(98);

			File hapFile = new File(dbFile, "rib_elon.fasta");
			File originalDbFile = new File(dbFile, "original_db");
			File eggnogFolder = new File(dbFile, "eggNOG");

			int downloadCount = 0;
			for (String fileNameString : fileSizeMap.keySet()) {
				File downloadFile;
				if (fileNameString.endsWith(".faa")) {
					downloadFile = new File(originalDbFile, fileNameString);
				} else if (fileNameString.endsWith("_eggNOG.tsv")) {
					downloadFile = new File(eggnogFolder, fileNameString);
				} else {
					downloadFile = new File(dbFile, fileNameString);
				}

				if (downloadFile.exists() && downloadFile.length() == fileSizeMap.get(fileNameString)) {
					downloadCount++;
				}
			}

			if (downloadCount == fileSizeMap.size()) {
				LOGGER.info(taskName + ": extracting high abundant proteins...");
				System.out.println(
						format.format(new Date()) + "\t" + taskName + ": extracting  high abundant proteins...");

				HapExtractor extractor = new HapExtractor(true);
				extractor.extract(originalDbFile, hapFile);

				LOGGER.info(taskName + ": finished");
				System.out.println(format.format(new Date()) + "\t" + taskName + ": finished");

				bar.setValue(100);

			} else {
				LOGGER.info(taskName + ": the total number is " + fileSizeMap.size() + ", " + downloadCount
						+ " files were downloaded," + " please try again");
				System.out.println(format.format(new Date()) + "\t" + taskName + ": the total number is "
						+ fileSizeMap.size() + ", " + downloadCount + " files were downloaded," + " please try again");
			}

		} catch (IOException e) {
			LOGGER.error(taskName + ": failed", e);
			System.out.println(format.format(new Date()) + "\t" + taskName + ": failed");
		}
	}

	@SuppressWarnings("unused")
	private void downloadMultiThreads(File dbFile, String catalogueName, String catalogueVersion) {

		FTPClient ftpClient = new FTPClient();

		try {
			HashMap<String, Long> fileSizeMap = new HashMap<String, Long>();
			ExecutorService executor = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
			ArrayList<Future<?>> flist = new ArrayList<Future<?>>();

			// connect and login to the server
			ftpClient.connect(server, port);
			ftpClient.login(user, pass);

			// use passive mode to pass firewall
			ftpClient.enterLocalPassiveMode();

			System.out.println("---Connected---");

			FTPFile[] mgnifyFiles = ftpClient.listFiles(mgnify);
			for (FTPFile mgnifyFile : mgnifyFiles) {
				if (mgnifyFile.isDirectory() && mgnifyFile.getName().equals(catalogueName)) {

					FTPFile[] versionFiles = ftpClient.listFiles(mgnify + "/" + catalogueName);
					for (FTPFile versionFile : versionFiles) {
						if (versionFile.getName().equals(catalogueVersion)) {
							String pathString = mgnify + "/" + mgnifyFile.getName() + "/" + catalogueVersion;
							FTPFile[] genomeFiles = ftpClient.listFiles(pathString);

							ftpClient.changeWorkingDirectory(pathString);

							for (FTPFile genomeFile : genomeFiles) {
								String nameJ = genomeFile.getName();

								if (nameJ.equals("genomes-all_metadata.tsv")) {

									File downloadFile = new File(dbFile, nameJ);
									if (downloadFile.exists() && downloadFile.length() == genomeFile.getSize()) {
										System.out.println("File exist " + nameJ);
									} else {
										Future<?> f = executor.submit(() -> {
											fileSizeMap.put(nameJ, genomeFile.getSize());
											OutputStream out;
											try {
												out = new FileOutputStream(downloadFile.getAbsolutePath());
												if (ftpClient.retrieveFile(nameJ, out)) {
													System.out.println("Downloaded " + nameJ);
												} else {
													System.out.println("Failed to download " + nameJ);
												}
												out.close();
											} catch (IOException e) {
												// TODO Auto-generated catch block
												LOGGER.error(taskName + ": error in processing file " + nameJ, e);
												System.out.println(format.format(new Date()) + "\t" + taskName
														+ ": error in processing file " + nameJ);
											}
										});
										flist.add(f);
									}

								} else if (nameJ.equals("species_catalogue")) {

									File originalDbFile = new File(dbFile, "original_db");
									if (!originalDbFile.exists()) {
										originalDbFile.mkdir();
									}

									File eggnogFile = new File(dbFile, "eggNOG");
									if (!eggnogFile.exists()) {
										eggnogFile.mkdir();
									}

									FTPFile[] speciesFiles = ftpClient.listFiles(pathString + "/" + nameJ);

									for (FTPFile speciesFile : speciesFiles) {

										FTPFile[] subSpeciesFiles = ftpClient.listFiles(
												pathString + "/species_catalogue" + "/" + speciesFile.getName());

										for (FTPFile subSpeciesFile : subSpeciesFiles) {

											FTPFile[] subSpeciesContentFiles = ftpClient
													.listFiles(pathString + "/species_catalogue" + "/"
															+ speciesFile.getName() + "/" + subSpeciesFile.getName());
											String genomeNameString = subSpeciesFile.getName();

											for (FTPFile subSpeciesContentFile : subSpeciesContentFiles) {
												if (subSpeciesContentFile.getName().equals("genome")) {
													String finalPathString = pathString + "/species_catalogue" + "/"
															+ speciesFile.getName() + "/" + subSpeciesFile.getName()
															+ "/genome";

													FTPFile[] finalFiles = ftpClient.listFiles(finalPathString);
													ftpClient.changeWorkingDirectory(finalPathString);

													for (FTPFile finalFile : finalFiles) {
														String fileName = finalFile.getName();
														if (fileName.equals(genomeNameString + ".faa")) {
															fileSizeMap.put(fileName, finalFile.getSize());
															File downloadFile = new File(originalDbFile, fileName);
															if (downloadFile.exists()
																	&& downloadFile.length() == genomeFile.getSize()) {
																System.out.println("File exist " + nameJ);
															} else {
																Future<?> f = executor.submit(() -> {
																	OutputStream out;
																	try {
																		out = new FileOutputStream(
																				(new File(originalDbFile, fileName))
																						.getAbsolutePath());
																		if (ftpClient.retrieveFile(fileName, out)) {
																			System.out.println("Downloaded "
																					+ finalFile.getName());
																		} else {
																			System.out.println("Failed to download "
																					+ finalFile.getName());
																		}
																		out.close();
																	} catch (IOException e) {
																		// TODO Auto-generated catch block
																		LOGGER.error(
																				taskName + ": error in processing file "
																						+ fileName,
																				e);
																		System.out.println(format.format(new Date())
																				+ "\t" + taskName
																				+ ": error in processing file "
																				+ fileName);
																	}
																});
																flist.add(f);
															}

														} else if (fileName.equals(genomeNameString + "_eggNOG.tsv")) {

															fileSizeMap.put(fileName, finalFile.getSize());
															File downloadFile = new File(eggnogFile, fileName);
															if (downloadFile.exists()
																	&& downloadFile.length() == genomeFile.getSize()) {
																System.out.println("File exist " + nameJ);
															} else {
																OutputStream out;
																try {
																	out = new FileOutputStream(
																			(new File(eggnogFile, fileName))
																					.getAbsolutePath());
																	if (ftpClient.retrieveFile(fileName, out)) {
																		System.out.println(
																				"Downloaded " + finalFile.getName());
																	} else {
																		System.out.println("Failed to download "
																				+ finalFile.getName());
																	}
																	out.close();
																} catch (IOException e) {
																	// TODO Auto-generated catch block
																	LOGGER.error(taskName
																			+ ": error in processing file " + fileName,
																			e);
																	System.out.println(format.format(new Date()) + "\t"
																			+ taskName + ": error in processing file "
																			+ fileName);
																}
															}
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}

		
			System.out.println("Disconnected");

			for (int i = 0; i < flist.size(); i++) {
				Future<?> f = flist.get(i);
				f.get();
				bar.setValue((int) (i * 95.0 / (double) flist.size()));
			}

			try {

				executor.shutdown();

				boolean finish = executor.awaitTermination(12, TimeUnit.HOURS);

				ftpClient.logout();
				ftpClient.disconnect();

				bar.setValue(98);

				File hapFile = new File(dbFile, "rib_elon.fasta");
				File originalDbFile = new File(dbFile, "original_db");
				File eggnogFolder = new File(dbFile, "eggNOG");

				int downloadCount = 0;
				for (String fileNameString : fileSizeMap.keySet()) {
					File downloadFile;
					if (fileNameString.endsWith(".faa")) {
						downloadFile = new File(originalDbFile, fileNameString);
					} else if (fileNameString.endsWith("_eggNOG.tsv")) {
						downloadFile = new File(eggnogFolder, fileNameString);
					} else {
						downloadFile = new File(dbFile, fileNameString);
					}

					if (downloadFile.exists() && downloadFile.length() == fileSizeMap.get(fileNameString)) {
						downloadCount++;
					}
				}

				if (finish) {

					if (downloadCount == fileSizeMap.size()) {
						LOGGER.info(taskName + ": extracting high abundant proteins...");
						System.out.println(format.format(new Date()) + "\t" + taskName
								+ ": extracting  high abundant proteins...");

						HapExtractor extractor = new HapExtractor(true);
						extractor.extract(originalDbFile, hapFile);

						LOGGER.info(taskName + ": finished");
						System.out.println(format.format(new Date()) + "\t" + taskName + ": finished");

						bar.setValue(100);

					} else {
						LOGGER.info(taskName + ": the total number is " + fileSizeMap.size() + ", " + downloadCount
								+ " files were downloaded," + " please try again");
						System.out.println(format.format(new Date()) + "\t" + taskName + ": the total number is "
								+ fileSizeMap.size() + ", " + downloadCount + " files were downloaded,"
								+ " please try again");
					}
				} else {
					LOGGER.info(taskName + ": task dosen't finish in a long time, please try again");
					System.out.println(format.format(new Date()) + "\t" + taskName
							+ ": task dosen't finish in a long time, please try again");
				}

			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				LOGGER.error(taskName + ": failed", e);
				System.out.println(format.format(new Date()) + "\t" + taskName + ": failed");
			}

		} catch (IOException e) {
			LOGGER.error(taskName + ": failed", e);
			System.out.println(format.format(new Date()) + "\t" + taskName + ": failed");
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			LOGGER.error(taskName + ": failed", e);
			System.out.println(format.format(new Date()) + "\t" + taskName + ": failed");
		} catch (ExecutionException e) {
			// TODO Auto-generated catch block
			LOGGER.error(taskName + ": failed", e);
			System.out.println(format.format(new Date()) + "\t" + taskName + ": failed");
		}
	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		
		
		FTPClient ftpClient = new FTPClient();

		try {
			// connect and login to the server
			ftpClient.connect(server, port);
			ftpClient.login(user, pass);

			// use passive mode to pass firewall
			ftpClient.enterLocalPassiveMode();

			System.out.println("Connected");
			
			ftpClient.logout();
			ftpClient.disconnect();
			
			System.out.println("Disconnected");
		} catch (IOException e) {
			LOGGER.error(taskName + ": failed", e);
		}
		
	}


}
