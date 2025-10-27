/**
 * 
 */
package uk.ac.quadram.metapilot.task;

import java.io.File;
import java.io.FileFilter;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.net.URL;
import java.nio.file.FileSystems;
import java.nio.file.Path;
import java.nio.file.StandardWatchEventKinds;
import java.nio.file.WatchEvent;
import java.nio.file.WatchKey;
import java.nio.file.WatchService;
import java.text.SimpleDateFormat;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

import javax.swing.UIManager;
import javax.swing.UnsupportedLookAndFeelException;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.apache.logging.log4j.core.config.builder.api.AppenderComponentBuilder;
import org.apache.logging.log4j.core.config.builder.api.ConfigurationBuilder;
import org.apache.logging.log4j.core.config.builder.api.ConfigurationBuilderFactory;
import org.apache.logging.log4j.core.config.builder.api.LayoutComponentBuilder;
import org.apache.logging.log4j.core.config.builder.api.RootLoggerComponentBuilder;
import org.apache.logging.log4j.core.config.builder.impl.BuiltConfiguration;

import javafx.scene.control.ProgressBar;
import javafx.application.Platform;
import uk.ac.quadram.metapilot.license.LicenseVerifier;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.task.par.MetaData;
import uk.ac.quadram.metapilot.task.par.BaseParameterIO;
import uk.ac.quadram.metapilot.task.par.MetaBaseParameter;

/**
 * MetaPilot provides a super class of MetaPilotTask. Different versions of
 * MetaPilotTask could inherit this class to inherit the method
 * {@code exportMetadata()} and {@code copyReport()}
 * 
 * (Migrated from SwingWorker to JavaFX Task)
 * 
 * @author Kai Cheng
 * @version 1.0.0
 */
public abstract class MetaPilotTask extends MetaAbstractTask {

	protected ProgressBar secondBar;
	
	// Instance logger - initialized in constructor AFTER configLog() to ensure proper configuration
	protected Logger logger;
	
	protected SimpleDateFormat format = new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss");
	
	// WatchService thread for monitoring file generation progress
	protected Thread watchServiceThread;

	protected ConfigurationBuilder<BuiltConfiguration> builder;

	protected File final_pep_txt;
	protected File final_pro_txt;
	protected File reportHtmlDir;

	private static final String taskName = "MetaPilot task";
	public static final String version = "1.0.0";
	public static final String versionFile = "_1_0_0.json";

	public MetaPilotTask(MetaBaseParameter metaPar, MetaResourceParameters mrp, ProgressBar progressBar,
			ProgressBar secondBar) {
		super(metaPar, mrp, progressBar);
		this.secondBar = secondBar;
		
		// Configure logging for this task's result folder and get the task-specific logger name
		String taskLoggerName = this.configLog(new File(metaPar.getResult()));
		
		// Initialize logger with the task-specific logger name (not the class name)
		// This ensures the logger uses the file appender we just configured
		if (taskLoggerName != null) {
			this.logger = LogManager.getLogger(taskLoggerName);
		} else {
			// Fallback to class-based logger if configuration failed
			this.logger = LogManager.getLogger(this.getClass());
			System.err.println("Warning: Using fallback logger configuration");
		}
		
		// Log task initialization
		logger.info("========================================");
		logger.info("Task initialized: " + getTaskName());
		logger.info("Result folder: " + metaPar.getResult());
		logger.info("========================================");
		
		this.progressProperty().addListener((observable, oldValue, newValue) -> {
			if (this.progressBar != null) {
				this.progressBar.setProgress(newValue.doubleValue());
			}
		});
	}

	public abstract void forceStop();

	protected abstract String getVersion();

	// Helper method to update progress on the Swing progress bar
	protected void setProgressValue(double value) {
		updateProgress(value, 100);
	}

	// Helper method to send a message from the background thread
	protected void publishMessage(String message) {
		updateMessage(message);
	}

	/**
	 * Compatibility method to start the task execution. Equivalent to SwingWorker's
	 * execute() method.
	 */
	public void execute() {
		Thread thread = new Thread(this);
		thread.setDaemon(true);
		thread.start();
	}

	protected void exportMetadata(File result) throws FileNotFoundException {
		File metafile = new File(result, "metainfo.txt");
		if (metafile.exists()) {
			return;
		}

		MetaData metadata = metaPar.getMetadata();
		int metaTypeCount = metadata.getMetaTypeCount();
		if (metaTypeCount == 0) {
			return;
		}

		PrintWriter writer = new PrintWriter(metafile);
		StringBuilder title = new StringBuilder();
		title.append("Raw").append("\t");
		title.append("Experiment");
		for (int i = 0; i < metaTypeCount; i++) {
			title.append("\t").append("meta" + (i + 1));
		}
		writer.println(title);

		if (metadata.getLabelTitle().length == 0) {
			String[] raws = metadata.getRawFiles();
			String[] expNames = metadata.getExpNames();
			String[][] metainfo = metadata.getMetaInfo();

			StringBuilder[] sbs = new StringBuilder[expNames.length];

			for (int i = 0; i < sbs.length; i++) {
				String raw = raws[i].substring(raws[i].lastIndexOf("\\"), raws[i].length() - ".raw".length());
				sbs[i] = new StringBuilder();
				sbs[i].append(raw).append("\t");
				sbs[i].append(expNames[i]);
				for (int j = 0; j < metainfo[i].length; j++) {
					sbs[i].append("\t").append(metainfo[i][j]);
				}
				writer.println(sbs[i]);
			}
		} else {
			String[] labelExpNames = metadata.getLabelExpNames();
			String[][] metainfo = metadata.getMetaInfo();

			StringBuilder[] sbs = new StringBuilder[labelExpNames.length];
			for (int i = 0; i < sbs.length; i++) {
				sbs[i] = new StringBuilder(labelExpNames[i]);
				for (int j = 0; j < metainfo[i].length; j++) {
					sbs[i].append("\t").append(metainfo[i][j]);
				}
				writer.println(sbs[i]);
			}
		}
		writer.close();

		metaPar.setMetafile(metafile);
	}

	/**
	 * Copy report templates in background thread with progress updates. This method
	 * should be called from the call() method after setting up the result
	 * directory.
	 * 
	 * @return true if copy was successful, false otherwise
	 */
	protected boolean copyReportTemplatesInBackground() {
		try {
			updateMessage("Setting up report templates...");
			// Get the result directory from metaPar
			File resultDir = metaPar.getDbSearchResultFile();
			if (resultDir == null) {
				logger.error(getTaskName() + ": Result directory is null, cannot copy report templates");
				return false;
			}

			// Create report directory
			File reportDir = new File(resultDir, "report");
			if (!reportDir.exists()) {
				reportDir.mkdirs();
			}
			metaPar.setReportDir(reportDir);

			// Copy report templates using MetaPythonReportCopyTask
			MetaPythonReportCopyTask reportCopyTask = new MetaPythonReportCopyTask(resultDir);
			Thread reportCopyThread = new Thread(reportCopyTask);
			reportCopyThread.start();

			updateMessage("Copying report templates...");

			try {
				reportCopyThread.join();
				if (!reportCopyTask.get()) {
					logger.error(getTaskName() + ": error in copying report template");
					updateMessage("Error: Failed to copy report templates");
					return false;
				}
			} catch (Exception e) {
				logger.error(getTaskName() + ": error in copying report template", e);
				updateMessage("Error: " + e.getMessage());
				return false;
			}

			reportHtmlDir = new File(reportDir + "\\reports");

			updateMessage("Report templates ready");

			logger.info(getTaskName() + ": Report templates copied successfully");
			return true;

		} catch (Exception e) {
			logger.error(getTaskName() + ": error in copying report templates", e);
			updateMessage("Error: " + e.getMessage());
			return false;
		}
	}

	protected String configLog(File logDir) {
		// Ensure log directory exists
		if (!logDir.exists()) {
			boolean created = logDir.mkdirs();
			if (!created) {
				System.err.println("Warning: Could not create log directory: " + logDir.getAbsolutePath());
				return null;
			}
		}

		try {
			// Get the current LoggerContext
			org.apache.logging.log4j.core.LoggerContext context = 
				(org.apache.logging.log4j.core.LoggerContext) LogManager.getContext(false);
			org.apache.logging.log4j.core.config.Configuration config = context.getConfiguration();
			
			// Create pattern layout
			String pattern = "%d{yyyy-MM-dd HH:mm:ss} %-5level %logger{36} - %msg%n";
			org.apache.logging.log4j.core.layout.PatternLayout layout = 
				org.apache.logging.log4j.core.layout.PatternLayout.newBuilder()
					.withPattern(pattern)
					.withConfiguration(config)
					.build();
			
			// Create file appender for task logging
			File logFile = new File(logDir, "logging.log");
			org.apache.logging.log4j.core.appender.FileAppender fileAppender = 
				org.apache.logging.log4j.core.appender.FileAppender.newBuilder()
					.setName("TaskFileAppender_" + System.currentTimeMillis())
					.withFileName(logFile.getAbsolutePath())
					.withAppend(true)
					.withImmediateFlush(true)
					.setLayout(layout)
					.setConfiguration(config)
					.build();
			
			fileAppender.start();
			config.addAppender(fileAppender);
			
			// Create unique logger name for this task
			String taskLoggerName = "TaskLogger_" + System.currentTimeMillis();
			
			// Create new LoggerConfig for this task (not using root logger)
			org.apache.logging.log4j.core.config.LoggerConfig taskLoggerConfig = 
				org.apache.logging.log4j.core.config.LoggerConfig.createLogger(
					false,  // additive = false, don't pass events to parent
					Level.DEBUG,
					taskLoggerName,
					"true",
					new org.apache.logging.log4j.core.config.AppenderRef[0],
					null,
					config,
					null
				);
			
			// Add the file appender to the task logger config
			taskLoggerConfig.addAppender(fileAppender, Level.DEBUG, null);
			
			// Explicitly set additive to false to prevent duplication in metapilot.log
			taskLoggerConfig.setAdditive(false);
			
			// Add the logger config to the configuration
			config.addLogger(taskLoggerName, taskLoggerConfig);
			
			// CRITICAL: Update loggers to activate the new configuration
			context.updateLoggers();
			
			System.out.println("Task logger configured successfully to: " + logFile.getAbsolutePath());
			System.out.println("Task logger name: " + taskLoggerName);
			
			// Redirect System.out and System.err to system_output.log in task folder
			File systemOutFile = new File(logDir, "system_output.log");
			java.io.PrintStream fileOut = new java.io.PrintStream(
				new java.io.FileOutputStream(systemOutFile, true), true);
			
			// Create TeeOutputStream to write to both original console and file
			TeeOutputStream teeOut = new TeeOutputStream(System.out, fileOut);
			TeeOutputStream teeErr = new TeeOutputStream(System.err, fileOut);
			
			System.setOut(new java.io.PrintStream(teeOut, true));
			System.setErr(new java.io.PrintStream(teeErr, true));
			
			System.out.println("========================================");
			System.out.println("System output/error redirected to: " + systemOutFile.getAbsolutePath());
			System.out.println("Task logger file: " + logFile.getAbsolutePath());
			System.out.println("========================================");
			
			// Return the task logger name so caller can get the right logger
			return taskLoggerName;
			
		} catch (Exception e) {
			System.err.println("Warning: Could not configure logging: " + e.getMessage());
			e.printStackTrace();
			return null;
		}
	}

	private static void run(String[] args) {

		CommandLineParser parser = new DefaultParser();
		CommandLine commandLine = null;
		try {
			commandLine = parser.parse(RunOption.getOptions(), args);
		} catch (ParseException e) {
			// TODO Auto-generated catch block
			System.err.println("Error parsing command line: " + e.getMessage());
			e.printStackTrace();
			return;
		}

		if (commandLine.hasOption("help")) {
			printUsage();
			return;
		}

		if (commandLine.hasOption("gui")) {

			try {
				UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
			} catch (ClassNotFoundException | InstantiationException | IllegalAccessException
					| UnsupportedLookAndFeelException e) {
				// TODO Auto-generated catch block
				System.err.println("Error in creating the GUI: " + e.getMessage());
				e.printStackTrace();
			}

			javax.swing.SwingUtilities.invokeLater(new Runnable() {
				public void run() {
					try {

					} catch (Exception e) {
						System.err.println("Error in creating the GUI: " + e.getMessage());
						e.printStackTrace();
					}
				}
			});

			return;
		}

		if (commandLine.hasOption("par")) {

			LicenseVerifier verifier = new LicenseVerifier();
			boolean verified = verifier.verify();
			if (!verified) {

				System.out.println("Thank you for using MetaPilot. A license is needed to run MetaPilot, "
						+ "please go to Tools -> Verification for the license.");
				return;
			}

			String par = commandLine.getOptionValue("par");
			File parafile = new File(par);
			if (parafile.exists() && parafile.length() > 0) {
				Object parameters = BaseParameterIO.loadParametersFromFilePath(par);

			}
		}
	}

	static class RunOption {

		private static final Options options = new Options();

		static {

			Option help = Option.builder("h").longOpt("help").desc("print this message").build();

			options.addOption(help);

			Option gui = Option.builder("g").longOpt("gui").desc("show the GUI").build();

			options.addOption(gui);

			Option paraPath = Option.builder("p").argName("file").longOpt("par").hasArg()
					.desc("path to the parameter file, which is used to run this program in cmd").build();

			options.addOption(paraPath);

		};

		public static Options getOptions() {
			return options;
		}

	}

	private static void printUsage() {
		HelpFormatter formatter = new HelpFormatter();
		formatter.printHelp("Metaproteomics data analysis platform\n",
				"\nBoth GUI and cmd running modes are supported, in cmd mode a parameter file is required.\n",
				RunOption.getOptions(), "\nOnline version see http://imetalab.ca/", true);
	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub

		MetaPilotTask.run(args);

	}

	/**
	 * Gets the application base directory using multiple fallback mechanisms
	 * 
	 * @return A File representing the application base directory
	 */
	public static File getApplicationBaseDirectory() {
		File baseDir = null;

		// First attempt: using ClassLoader.getSystemClassLoader().getResource(".")
		URL resourceUrl = ClassLoader.getSystemClassLoader().getResource(".");
		if (resourceUrl != null) {
			try {
				baseDir = new File(resourceUrl.getPath());
				System.out.println("Using ClassLoader.getSystemClassLoader().getResource(\".\") to find base directory: "
						+ baseDir);
				return baseDir;
			} catch (Exception e) {
				System.err.println(
						"Could not use ClassLoader.getSystemClassLoader().getResource(\".\") to find base directory: "
						+ e.getMessage());
			}
		}

		// Second attempt: using
		// Thread.currentThread().getContextClassLoader().getResource(".")
		ClassLoader cl = Thread.currentThread().getContextClassLoader();
		if (cl != null) {
			URL contextResource = cl.getResource(".");
			if (contextResource != null) {
				try {
					baseDir = new File(contextResource.getPath());
					System.out.println(
							"Using Thread.currentThread().getContextClassLoader().getResource(\".\") to find base directory: "
									+ baseDir);
					return baseDir;
				} catch (Exception e) {
					System.err.println(
							"Could not use Thread.currentThread().getContextClassLoader().getResource(\".\") to find base directory: "
							+ e.getMessage());
				}
			}
		}

		// Third attempt: using current working directory
		baseDir = new File(System.getProperty("user.dir"));
		System.out.println("Using System.getProperty(\"user.dir\") to find base directory: " + baseDir);

		// Make sure Resources directory exists
		File resourcesDir = new File(baseDir, "Resources");
		if (!resourcesDir.exists()) {
			System.err.println("WARNING: Resources directory not found at: " + resourcesDir.getAbsolutePath());

			// If Resources directory doesn't exist in the current directory,
			// check if it exists one directory up (happens sometimes with projects run from
			// IDEs)
			File parentDir = baseDir.getParentFile();
			if (parentDir != null) {
				File parentResourcesDir = new File(parentDir, "Resources");
				if (parentResourcesDir.exists()) {
					baseDir = parentDir;
					System.out.println("Found Resources directory in parent directory: " + baseDir);
				}
			}
		}

		return baseDir;
	}

	/**
	 * Generic WatchServiceLoop that can be used by any MetaPilotTask subclass
	 * to monitor file generation progress. This eliminates code duplication
	 * across different workflow tasks.
	 */
	public class WatchServiceLoop implements Runnable {

		private WatchService watchService;
		private FileFilter filter;
		private int totalCount = 0;
		private int currentCount = 0;
		private int barStart;
		private int barEnd;

		public WatchServiceLoop(WatchService watchService, FileFilter filter, int totalCount, int currentCount,
				int barStart, int barEnd) {
			this.watchService = watchService;
			this.filter = filter;
			this.totalCount = totalCount;
			this.currentCount = currentCount;
			this.barStart = barStart;
			this.barEnd = barEnd;
		}

		@SuppressWarnings("unchecked")
		@Override
		public void run() {
			while (currentCount < totalCount) {
				try {
					// Wait for a watch key to be available
					WatchKey key = watchService.take();
					// Iterate over the watch events
					for (WatchEvent<?> event : key.pollEvents()) {
						// Get the event kind and context
						WatchEvent.Kind<?> kind = event.kind();
						WatchEvent<Path> pathEvent = (WatchEvent<Path>) event;
						Path fileName = pathEvent.context();
						
						// If the event is a creation of a file matching the filter, increment the count and update progress
						if (kind == StandardWatchEventKinds.ENTRY_CREATE && filter.accept(fileName.toFile())) {
							currentCount++;
							int percentage = (int) ((double) currentCount / (double) totalCount
									* (double) (barEnd - barStart));

							// Update overall progress bar using the task's updateProgress method
							MetaPilotTask.this.updateProgress(barStart + percentage, 100);
							
							// Update current task bar (second bar) to show progress within current step
							if (secondBar != null) {
								double currentTaskProgress = (double) currentCount / (double) totalCount;
								Platform.runLater(() -> secondBar.setProgress(currentTaskProgress));
							}

							System.out.println(format.format(new Date()) + "\t" + getTaskName()
									+ ": WatchService found one finished task at " + fileName.toFile());

							logger.info("WatchService found one finished task at {}", fileName.toFile());
						}
					}
					// Reset the key
					key.reset();
				} catch (InterruptedException e) {
					logger.info("WatchService thread interrupted, stopping monitoring");
					Thread.currentThread().interrupt();
					break;
				} catch (Exception e) {
					logger.error("Error in WatchService loop", e);
					System.out.println(format.format(new Date()) + "\t" + getTaskName()
							+ ": error in WatchService loop: " + e.getMessage());
				}
			}
			logger.info("WatchService monitoring completed. Found {}/{} files", currentCount, totalCount);
		}
	}

	/**
	 * Helper method to create and start a WatchServiceLoop for monitoring file generation.
	 * This method handles the setup of WatchService and starts the monitoring thread.
	 * 
	 * @param folders Array of folders to monitor
	 * @param subFolderName Subfolder name within each folder to monitor (e.g., "hap", "result")
	 * @param filter FileFilter to determine which files to count
	 * @param totalCount Total number of files expected
	 * @param barStart Progress bar start percentage
	 * @param barEnd Progress bar end percentage
	 * @return true if WatchService was set up successfully, false otherwise
	 */
	protected boolean startWatchServiceLoop(File[] folders, String subFolderName, FileFilter filter, 
			int totalCount, int barStart, int barEnd) {
		WatchService watchService = null;
		
		try {
			logger.info("Setting up WatchService to monitor {} folders", folders.length);
			watchService = FileSystems.getDefault().newWatchService();
			
			// Register each folder with the watch service
			for (File folder : folders) {
				File subFolder = new File(folder, subFolderName);
				if (subFolder.exists()) {
					Path path = subFolder.toPath();
					path.register(watchService, StandardWatchEventKinds.ENTRY_CREATE);
					logger.debug("Registered folder for monitoring: {}", subFolder.getAbsolutePath());
				} else {
					logger.warn("Subfolder does not exist, cannot monitor: {}", subFolder.getAbsolutePath());
				}
			}
			logger.info("WatchService setup completed successfully");
		} catch (IOException e) {
			logger.error("Error creating WatchService", e);
			System.out.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": error in creating WatchService: " + e.getMessage());
			return false;
		}

		logger.info("Starting WatchService monitoring thread");
		watchServiceThread = new Thread(new WatchServiceLoop(watchService, filter, totalCount, 0, barStart, barEnd));
		watchServiceThread.start();
		
		return true;
	}

	/**
	 * Helper method to stop the WatchService thread if it's running.
	 * This should be called at the end of methods that use WatchService monitoring.
	 */
	protected void stopWatchServiceThread() {
		if (watchServiceThread != null && watchServiceThread.isAlive()) {
			logger.info("Stopping WatchService thread");
			watchServiceThread.interrupt();
		}
	}

	/**
	 * Sophisticated Occam's razor implementation for selecting items (proteins,
	 * genomes, etc.) based on peptide coverage. This method implements the advanced
	 * algorithm with iterative re-sorting and multiple breakpoints used by the MAG
	 * workflow's refineGenomes method.
	 * 
	 * @param itemToPeptides Map from item names to their associated peptides
	 * @param totalPeptides  Total set of unique peptides (for threshold
	 *                       calculation)
	 * @param threshold      Coverage threshold (e.g., 0.95 for 95% coverage)
	 * @param stepThreshold  Increment threshold for breakpoint detection (e.g.,
	 *                       0.01 for 1%)
	 * @param maximum        Maximum number of items to select
	 * @param decoyFilter    Whether to apply decoy filtering logic
	 * @param out            Output file for ranking statistics (optional, can be
	 *                       null)
	 * @return Set of selected item names that together cover the threshold
	 *         percentage of peptides
	 */
	protected HashSet<String> applyOccamsRazor(HashMap<String, HashSet<String>> itemToPeptides,
			Set<String> totalPeptides, double threshold, double stepThreshold, int maximum, boolean decoyFilter,
			File out) {

		int totalPepCount = totalPeptides.size();
		String[] items = itemToPeptides.keySet().toArray(new String[itemToPeptides.size()]);

		// Initial sort by peptide count (descending)
		Arrays.sort(items, (i1, i2) -> {
			int s1 = itemToPeptides.get(i1).size();
			int s2 = itemToPeptides.get(i2).size();
			return s2 - s1;
		});

		// Iterative re-sorting based on actual contributions (5 iterations for
		// optimization)
		for (int iteratorCount = 0; iteratorCount < 5; iteratorCount++) {
			HashSet<String> totalSet = new HashSet<String>();
			HashMap<String, Integer> orderMap = new HashMap<String, Integer>();
			for (int i = 0; i < items.length; i++) {
				int currentSize = totalSet.size();
				totalSet.addAll(itemToPeptides.get(items[i]));
				int increaseSize = totalSet.size() - currentSize;
				orderMap.put(items[i], increaseSize);
			}

			// Re-sort by actual increment contributions
			Arrays.sort(items, (i1, i2) -> {
				return orderMap.get(i2) - orderMap.get(i1);
			});
		}

		double[] percentage = new double[items.length];
		HashSet<String> selectedItems = new HashSet<String>();
		HashSet<String> totalSet = new HashSet<String>();
		int breakPoint1 = -1;
		int breakPoint2 = -1;

		// Write ranking statistics if output file is provided
		if (out != null) {
			try (PrintWriter writer = new PrintWriter(out)) {
				writer.println(
						"Rank\tItem\tPeptide_count\tTotal_covered_peptide_count\tPercentage\tIncrement\tIncrement ratio");
				for (int i = 0; i < items.length; i++) {

					int currentSize = totalSet.size();
					totalSet.addAll(itemToPeptides.get(items[i]));
					int add = totalSet.size() - currentSize;
					double additional = ((double) add / (double) totalPepCount);

					percentage[i] = ((double) totalSet.size() / (double) totalPepCount);

					writer.println(i + "\t" + items[i] + "\t" + itemToPeptides.get(items[i]).size() + "\t"
							+ totalSet.size() + "\t" + percentage[i] + "\t" + add + "\t" + additional);

					if (add > 1) {
						breakPoint1 = i;
					}
					if (additional > stepThreshold) {
						breakPoint2 = i;
					}
				}
				writer.close();
			} catch (IOException e) {
				logger.error(getTaskName() + ": error in writing ranking to " + out, e);
				System.out.println(
						format.format(new Date()) + "\t" + getTaskName() + ": error in writing ranking to " + out);
			}
		} else {
			// Calculate breakpoints without writing to file
			for (int i = 0; i < items.length; i++) {
				int currentSize = totalSet.size();
				totalSet.addAll(itemToPeptides.get(items[i]));
				int add = totalSet.size() - currentSize;
				double additional = ((double) add / (double) totalPepCount);

				percentage[i] = ((double) totalSet.size() / (double) totalPepCount);

				if (add > 1) {
					breakPoint1 = i;
				}
				if (additional > stepThreshold) {
					breakPoint2 = i;
				}
			}
		}

		// Reset totalSet for actual selection
		totalSet.clear();

		int decoyCount = 0;
		for (int i = 0; i < items.length; i++) {
			selectedItems.add(items[i]);

			// Check for decoy items (common decoy patterns)
			if (items[i].startsWith("REV_") || items[i].contains("DECOY") || items[i].startsWith("REV")
					|| items[i].endsWith("_REV")) {
				decoyCount++;
			}

			if (decoyFilter) {
				if (decoyCount >= 3) {
					if (percentage[i] >= threshold) {
						break;
					}
					if (i == breakPoint1 || i == breakPoint2) {
						break;
					}
					if (selectedItems.size() > maximum) {
						break;
					}
				}
			} else {
				if (percentage[i] >= threshold) {
					break;
				}
				if (threshold < 1.0) {
					if (i == breakPoint1 || i == breakPoint2) {
						break;
					}
					if (selectedItems.size() > maximum) {
						break;
					}
				}
			}
		}

		return selectedItems;
	}

	/**
	 * Simplified version of applyOccamsRazor for backward compatibility. Uses
	 * default values for advanced parameters.
	 * 
	 * @param itemToPeptides Map from item names to their associated peptides
	 * @param totalPeptides  Total set of unique peptides (for threshold
	 *                       calculation)
	 * @param threshold      Coverage threshold (e.g., 0.95 for 95% coverage)
	 * @param out            Output file for ranking statistics
	 * @return Set of selected item names that together cover the threshold
	 *         percentage of peptides
	 */
	protected HashSet<String> applyOccamsRazor(HashMap<String, HashSet<String>> itemToPeptides,
			Set<String> totalPeptides, double threshold, File out) {
		return applyOccamsRazor(itemToPeptides, totalPeptides, threshold, 0.01, Integer.MAX_VALUE, false, out);
	}

	/**
	 * Get the logger for this task
	 * 
	 * @return The logger instance
	 */
	@Override
	protected Logger getLogger() {
		return logger;
	}

}
