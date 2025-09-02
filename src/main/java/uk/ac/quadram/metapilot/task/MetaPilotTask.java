/**
 * 
 */
package uk.ac.quadram.metapilot.task;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.net.JarURLConnection;
import java.net.URL;
import java.net.URLConnection;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.Enumeration;
import java.util.jar.JarEntry;
import java.util.jar.JarFile;

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
import org.json.JSONException;
import org.json.JSONObject;

import uk.ac.quadram.metapilot.task.pfind.gui.MetaPilotPFindMainFrame;
import uk.ac.quadram.metapilot.task.pfind.par.MetaParameterPFind;
import uk.ac.quadram.metapilot.gui.JavaFXProgressBar;
import uk.ac.quadram.metapilot.license.LicenseVerifier;
import uk.ac.quadram.metapilot.resources.MetaResourceIO;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.task.par.MetaData;
import uk.ac.quadram.metapilot.task.par.MetaBaseParameter;
import uk.ac.quadram.metapilot.task.par.MetaParameterIOHandler;

/**
 * MetaPilot provides a super class of MetaPilotTask.
 * Different versions of MetaPilotTask could inherit this class to inherit the
 * method {@code exportMetadata()} and {@code copyReport()}
 * 
 * (Migrated from SwingWorker to JavaFX Task)
 * 
 * @author Kai Cheng
 * @version 1.0.0
 */
public abstract class MetaPilotTask extends MetaAbstractTask {

	protected JavaFXProgressBar secondBar;
	protected static Logger LOGGER = LogManager.getLogger(MetaPilotTask.class);
	protected SimpleDateFormat format = new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss");
	
	protected ConfigurationBuilder<BuiltConfiguration> builder;
	
	protected File final_pep_txt;
	protected File final_pro_txt;

	private static final String taskName = "MetaPilot task";
	public static final String version = "1.0.0";
	public static final String versionFile = "_1_0_0.json";

	public MetaPilotTask(MetaBaseParameter metaPar, MetaResourceParameters mrp, JavaFXProgressBar progressBar,
			JavaFXProgressBar secondBar) {
		super(metaPar, mrp, progressBar);
		this.secondBar = secondBar;
		this.configLog(new File(metaPar.getResult()));
		this.progressProperty().addListener((observable, oldValue, newValue) -> {
			if (this.progressBar != null) {
				this.progressBar.setValue(newValue.intValue() * 100);
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
	 * Compatibility method to start the task execution.
	 * Equivalent to SwingWorker's execute() method.
	 */
	public void execute() {
	    Thread thread = new Thread(this);
	    thread.setDaemon(true);
	    thread.start();
	}

	protected void exportMetadata(File result) throws FileNotFoundException {
		if (metaPar.getWorkflowType() == MetaLabWorkflowType.TaxonomyAnalysis
				|| metaPar.getWorkflowType() == MetaLabWorkflowType.FunctionalAnnotation) {
			return;
		}

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
	 * @param result
	 * @throws IOException
	 */
	protected void copyReport(File result) throws IOException {

		File destDir = new File(result, "report");
		if (!destDir.exists()) {
			destDir.mkdirs();
		}

		try {
			loadRecourseFromJarByFolder("/uk/ac/quadram/metapilot/report", destDir);
		} catch (IOException e) {
			LOGGER.error(taskName + ": error in creating the report folder", e);
			System.err
					.println(format.format(new Date()) + "\t" + taskName + ": error in creating the report folder");
		}

		metaPar.setReportDir(destDir);
	}

	protected void loadRecourseFromJarByFolder(String fromPath, File toPath) throws IOException {

		URL url = getClass().getResource(fromPath);

		LOGGER.info(taskName + ": copying report folder from " + url + " started");

		URLConnection urlConnection = url.openConnection();
		if (url.getProtocol().equals("jar")) {
			copyJarResources((JarURLConnection) urlConnection, toPath);
		} else if (url.getProtocol().equals("file")) {
			copyFileResources(url, fromPath, toPath);
		} else {
			LOGGER.info(taskName + ": unknown resource file type " + url.getProtocol() + " for " + url);
		}

		LOGGER.info(taskName + ": copying report folder from " + url + " finished");
	}

	protected void copyFileResources(URL url, String fromPath, File toPath) throws IOException {
		File root = new File(url.getPath());
		if (root.isDirectory()) {
			File[] files = root.listFiles();
			for (File file : files) {
				if (file.isDirectory()) {
					File newToFile = new File(toPath, file.getName());
					if (!newToFile.exists()) {
						newToFile.mkdirs();
					}
					loadRecourseFromJarByFolder(fromPath + "/" + file.getName(), newToFile);
				} else {
					copyFromFile(fromPath + "/" + file.getName(), toPath);
				}
			}
		}
	}

	protected void copyJarResources(JarURLConnection jarURLConnection, File toPath) throws IOException {
		JarFile jarFile = jarURLConnection.getJarFile();
		Enumeration<JarEntry> entrys = jarFile.entries();
		while (entrys.hasMoreElements()) {
			JarEntry entry = entrys.nextElement();
			if (entry.getName().startsWith(jarURLConnection.getEntryName()) && !entry.getName().endsWith("/")) {
				copyFromJar("/" + entry.getName(), toPath);
			}
		}
		jarFile.close();
	}

	protected void copyFromFile(String fromPath, File toPath) throws IOException {

		int index = fromPath.lastIndexOf('/');

		String filename = fromPath.substring(index + 1);
		File file = new File(toPath, filename);

		if (!file.exists() && !file.createNewFile()) {
			return;
		}

		byte[] buffer = new byte[1024];
		int readBytes;

		URL fromUrl = getClass().getResource(fromPath);
		URLConnection urlConnection = fromUrl.openConnection();
		InputStream is = urlConnection.getInputStream();

		if (is == null) {
			throw new FileNotFoundException("File " + fromPath + " was not found inside JAR.");
		}

		OutputStream os = new FileOutputStream(file);
		try {
			while ((readBytes = is.read(buffer)) != -1) {
				os.write(buffer, 0, readBytes);
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			LOGGER.error("Error in writing report file to " + file, e);
			System.err.println(format.format(new Date()) + "\t" + "Error in writing report file to " + file);
		} finally {
			os.close();
			is.close();
		}
	}

	protected void copyFromJar(String fromPath, File toPath) throws IOException {

		String toPathName = toPath.getName();
		int index = fromPath.indexOf("/" + toPathName + "/");

		File file = new File(toPath.getParent() + fromPath.substring(index));

		File parent = file.getParentFile();

		if (!parent.exists()) {
			parent.mkdirs();
		}

		byte[] buffer = new byte[1024];
		int readBytes;

		URL fromUrl = getClass().getResource(fromPath);
		URLConnection urlConnection = fromUrl.openConnection();
		InputStream is = urlConnection.getInputStream();

		if (is == null) {
			throw new FileNotFoundException("File " + fromPath + " was not found inside JAR.");
		}

		OutputStream os = new FileOutputStream(file);
		try {
			while ((readBytes = is.read(buffer)) != -1) {
				os.write(buffer, 0, readBytes);
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			LOGGER.error("Error in writing report file to " + file, e);
			System.err.println(format.format(new Date()) + "\t" + "Error in writing report file to " + file);
		} finally {
			os.close();
			is.close();
		}
	}

	protected void configLog(File logDir) {
		String pattern = "%d %p %c [%t] %m%n";

		this.builder = ConfigurationBuilderFactory.newConfigurationBuilder();

		builder.setStatusLevel(Level.INFO);
		builder.setConfigurationName("DefaultFileLogger");

		// set the pattern layout and pattern
		LayoutComponentBuilder layoutBuilder = builder.newLayout("PatternLayout").addAttribute("pattern", pattern);

		// Ensure log directory exists
		if (!logDir.exists()) {
			boolean created = logDir.mkdirs();
			if (!created) {
				System.err.println("Warning: Could not create log directory: " + logDir.getAbsolutePath());
			}
		}

		// create a file appender
		File logFile = new File(logDir, "logging.log");
		AppenderComponentBuilder appenderBuilder = builder.newAppender("LogToFile", "File")
				.addAttribute("fileName", logFile.getAbsolutePath()).add(layoutBuilder);

		builder.add(appenderBuilder);

		RootLoggerComponentBuilder rootLogger = builder.newRootLogger(Level.INFO);
		rootLogger.add(builder.newAppenderRef("LogToFile"));
		builder.add(rootLogger);

		try {
			Configurator.reconfigure(builder.build());
		} catch (Exception e) {
			System.err.println("Warning: Could not configure logging to file: " + e.getMessage());
			// Continue without file logging - console logging will still work
		}
	}

	private static void run(String[] args) {

		CommandLineParser parser = new DefaultParser();
		CommandLine commandLine = null;
		try {
			commandLine = parser.parse(RunOption.getOptions(), args);
		} catch (ParseException e) {
			// TODO Auto-generated catch block
			LOGGER.error(e);
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
				LOGGER.error("Error in creating the GUI", e);
			}

			javax.swing.SwingUtilities.invokeLater(new Runnable() {
				public void run() {
					try {

						MetaPilotPFindMainFrame frame = new MetaPilotPFindMainFrame();
						frame.setVisible(true);

					} catch (Exception e) {
						LOGGER.error("Error in creating the GUI", e);
					}
				}
			});

			return;
		}

		if (commandLine.hasOption("par")) {

			LicenseVerifier verifier = new LicenseVerifier();
			boolean verified = verifier.verify();
			if (!verified) {

				LOGGER.info("Thank you for using MetaPilot. A license is needed to run MetaPilot, "
						+ "please go to Tools -> Verification for the license.");
				return;
			}

			String par = commandLine.getOptionValue("par");
			File parafile = new File(par);
			if (parafile.exists() && parafile.length() > 0) {
				StringBuilder sb = new StringBuilder();
				String line = null;
				try (BufferedReader reader = new BufferedReader(new FileReader(parafile))) {
					while ((line = reader.readLine()) != null) {
						sb.append(line);
					}
					reader.close();
				} catch (IOException e) {
					// TODO Auto-generated catch block
					LOGGER.error("Error in reading MetaLab parameter file " + parafile, e);
				}

				JSONObject obj = null;
				try {
					obj = new JSONObject(sb.toString());
				} catch (JSONException e) {
					LOGGER.error("Error in reading MetaLab parameter file " + parafile, e);
				}

				MetaBaseParameter loadPar = MetaParameterIOHandler.parse(obj);

				if (loadPar.getWorkflowType() == MetaLabWorkflowType.pFindWorkflow) {

					MetaParameterPFind metaPar = (MetaParameterPFind) loadPar;

					JavaFXProgressBar bar1 = new JavaFXProgressBar();
					JavaFXProgressBar bar2 = new JavaFXProgressBar();

					File file = getApplicationBaseDirectory();
					File setting = new File(file, "resources" + versionFile);

					MetaResourceIO mrio = new MetaResourceIO(setting);
					MetaResourceParameters mrp = mrio.loadParameters();


					File report = new File(metaPar.getReportDir(), "index.html");
					try {
						java.awt.Desktop.getDesktop().browse(report.toURI());
					} catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}

				} else {

				}
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
	 * @return A File representing the application base directory
	 */
	public static File getApplicationBaseDirectory() {
		File baseDir = null;
		
		// First attempt: using ClassLoader.getSystemClassLoader().getResource(".")
		URL resourceUrl = ClassLoader.getSystemClassLoader().getResource(".");
		if (resourceUrl != null) {
			try {
				baseDir = new File(resourceUrl.getPath());
				LOGGER.info("Using ClassLoader.getSystemClassLoader().getResource(\".\") to find base directory: " + baseDir);
				return baseDir;
			} catch (Exception e) {
				LOGGER.warn("Could not use ClassLoader.getSystemClassLoader().getResource(\".\") to find base directory", e);
			}
		}
		
		// Second attempt: using Thread.currentThread().getContextClassLoader().getResource(".")
		ClassLoader cl = Thread.currentThread().getContextClassLoader();
		if (cl != null) {
			URL contextResource = cl.getResource(".");
			if (contextResource != null) {
				try {
					baseDir = new File(contextResource.getPath());
					LOGGER.info("Using Thread.currentThread().getContextClassLoader().getResource(\".\") to find base directory: " + baseDir);
					return baseDir;
				} catch (Exception e) {
					LOGGER.warn("Could not use Thread.currentThread().getContextClassLoader().getResource(\".\") to find base directory", e);
				}
			}
		}
		
		// Third attempt: using current working directory
		baseDir = new File(System.getProperty("user.dir"));
		LOGGER.info("Using System.getProperty(\"user.dir\") to find base directory: " + baseDir);
		
		// Make sure Resources directory exists
		File resourcesDir = new File(baseDir, "Resources");
		if (!resourcesDir.exists()) {
			LOGGER.warn("Resources directory not found at: " + resourcesDir.getAbsolutePath());
			
			// If Resources directory doesn't exist in the current directory, 
			// check if it exists one directory up (happens sometimes with projects run from IDEs)
			File parentDir = baseDir.getParentFile();
			if (parentDir != null) {
				File parentResourcesDir = new File(parentDir, "Resources");
				if (parentResourcesDir.exists()) {
					baseDir = parentDir;
					LOGGER.info("Found Resources directory in parent directory: " + baseDir);
				}
			}
		}
		
		return baseDir;
	}

}
