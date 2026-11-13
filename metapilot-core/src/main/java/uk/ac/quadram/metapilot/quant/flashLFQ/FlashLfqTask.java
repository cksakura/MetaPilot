package uk.ac.quadram.metapilot.quant.flashLFQ;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileFilter;
import java.io.FileNotFoundException;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;
import java.util.concurrent.TimeUnit;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import uk.ac.quadram.metapilot.task.par.MetaData;

/**
 * Before using FlashLFQ, please make sure to modify the file
 * LicenceAgreements.toml with HasAcceptedThermoLicence = true
 * 
 * <p>Launch is OS-aware:
 * <ul>
 * <li>Windows: {@code cmd.exe /c call CMD.exe ...} (existing GUI behaviour)</li>
 * <li>Linux/macOS: {@code dotnet CMD.dll ...} (HPC pattern; uses {@code DOTNET_ROOT}
 * when set, otherwise {@code dotnet} on PATH)</li>
 * </ul>
 * {@code FLASHLFQ} may point at {@code CMD.exe}, {@code CMD.dll}, or the FlashLFQ
 * install directory.
 *
 * @author Kai Cheng
 *
 */
public class FlashLfqTask {

	private String flashLFQ;
	private MetaData metadata;
	private String idt;
	private String rep;
	private boolean mbr;
	private boolean nor;
	private int thr;
	
	// Task state management
	private volatile boolean isRunning = false;
	private volatile boolean isCancelled = false;
	private volatile boolean isDone = false;
	private Process flashLFQProcess = null;

	protected Logger LOGGER = LogManager.getLogger(FlashLfqTask.class);
	protected String taskName = "Label free quantification";
	protected SimpleDateFormat format = new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss");

	private static final String DOTNET_STRING_WINDOWS = "https://dotnet.microsoft.com/en-us/download/dotnet/8.0/runtime?cid=getdotnetcore&os=windows&arch=x64";
	private static final String DOTNET_STRING_LINUX = "https://dotnet.microsoft.com/en-us/download/dotnet/8.0";

	public FlashLfqTask(String flashLFQ, int thr) {
		this.flashLFQ = flashLFQ;
		this.thr = thr;
		this.thermoLicence();
	}

	private void thermoLicence() {
		File flashLFQFile = resolveFlashLfqArtifact(new File(flashLFQ), isWindowsOs());
		File parent = flashLFQFile.getParentFile();
		if (parent == null) {
			LOGGER.error(taskName + ": error in creating the license file (no parent for " + flashLFQFile + ")");
			return;
		}
		File rawLicenceFile = new File(parent, "LicenceAgreements.toml");
		try (PrintWriter writer = new PrintWriter(rawLicenceFile)) {
			writer.println("HasAcceptedThermoLicence = true");
		} catch (IOException e) {
			LOGGER.error(taskName + ": error in creating the license file", e);
			System.out.println(format.format(new Date()) + "\t" + taskName + ": error in creating the license file");
		}
	}

	/** True when the JVM is running on Windows. */
	static boolean isWindowsOs() {
		return isWindowsOs(System.getProperty("os.name"));
	}

	static boolean isWindowsOs(String osName) {
		return osName != null && osName.toLowerCase(Locale.ROOT).startsWith("windows");
	}

	/**
	 * Resolve {@code FLASHLFQ} config to the concrete artifact to launch
	 * ({@code CMD.exe} / {@code CMD.dll} / wrapper script).
	 */
	static File resolveFlashLfqArtifact(File configured, boolean windows) {
		if (configured == null) {
			throw new IllegalArgumentException("FlashLFQ path is null");
		}
		if (configured.isDirectory()) {
			File exe = new File(configured, "CMD.exe");
			File dll = new File(configured, "CMD.dll");
			if (windows) {
				if (exe.isFile()) {
					return exe;
				}
				if (dll.isFile()) {
					return dll;
				}
			} else {
				if (dll.isFile()) {
					return dll;
				}
				if (exe.isFile()) {
					return exe;
				}
			}
			return configured;
		}
		if (!windows && configured.getName().toLowerCase(Locale.ROOT).endsWith(".exe")) {
			File siblingDll = new File(configured.getParentFile(), "CMD.dll");
			if (siblingDll.isFile()) {
				return siblingDll;
			}
		}
		return configured;
	}

	/**
	 * Prefer {@code DOTNET_ROOT/dotnet} when set (HPC installs), else bare {@code dotnet}.
	 */
	static String resolveDotnetCommand() {
		String root = System.getenv("DOTNET_ROOT");
		if (root != null && !root.isBlank()) {
			File unix = new File(root.trim(), "dotnet");
			if (unix.isFile()) {
				return unix.getAbsolutePath();
			}
			File win = new File(root.trim(), "dotnet.exe");
			if (win.isFile()) {
				return win.getAbsolutePath();
			}
		}
		return "dotnet";
	}

	/**
	 * Build the process command for FlashLFQ for the given OS and artifact.
	 * Package-visible for unit tests.
	 */
	static List<String> buildFlashLfqCommand(File artifact, boolean windows, String idt, String rep, String out,
			boolean mbr, boolean nor, int thr, String dotnetCommand) {
		List<String> list = new ArrayList<>();
		String name = artifact.getName().toLowerCase(Locale.ROOT);
		boolean isDll = name.endsWith(".dll");
		boolean isExe = name.endsWith(".exe");

		if (windows && (isExe || !isDll)) {
			// Historical Windows GUI launch: cmd.exe /c call CMD.exe ...
			list.add("cmd.exe");
			list.add("/c");
			list.add("call");
			list.add(artifact.getAbsolutePath());
			list.add("--idt");
			list.add("\"" + idt + "\"");
			list.add("--rep");
			list.add("\"" + rep + "\"");
			list.add("--out");
			list.add("\"" + out + "\"");
			addSwitches(list, mbr, nor);
			list.add("--thr");
			list.add(String.valueOf(thr));
			return list;
		}

		if (isDll) {
			String dotnet = (dotnetCommand == null || dotnetCommand.isBlank()) ? "dotnet" : dotnetCommand;
			list.add(dotnet);
			list.add(artifact.getAbsolutePath());
		} else {
			// Shell wrapper or native binary
			list.add(artifact.getAbsolutePath());
		}
		list.add("--idt");
		list.add(idt);
		list.add("--rep");
		list.add(rep);
		list.add("--out");
		list.add(out);
		addSwitches(list, mbr, nor);
		list.add("--thr");
		list.add(String.valueOf(thr));
		return list;
	}

	/**
	 * {@code --mbr} and {@code --nor} are switches in FlashLFQ's command line (2.0.0.308): the flag
	 * alone turns the option on and a following {@code false} is ignored, so {@code --nor false}
	 * normalised. Each is passed only when wanted. {@code --mbr} defaults to true and no syntax turns
	 * it off (checked 2026-11-09 with {@code --mbr false}, {@code --mbr=false} and no flag: every one
	 * wrote {@code MatchBetweenRuns = true}), so a caller that wants no match-between-runs must drop
	 * the rows FlashLFQ reports with detection type {@code MBR}.
	 */
	private static void addSwitches(List<String> list, boolean mbr, boolean nor) {
		if (mbr) {
			list.add("--mbr");
		}
		if (nor) {
			list.add("--nor");
		}
	}

	private static String dotnetHelpUrl(boolean windows) {
		return windows ? DOTNET_STRING_WINDOWS : DOTNET_STRING_LINUX;
	}

	public boolean run(MetaData metadata, String idt, String rep, boolean mbr, boolean nor) {
		// Reset task state
		isRunning = true;
		isCancelled = false;
		isDone = false;
		
		this.metadata = metadata;
		this.idt = idt;
		this.rep = rep;
		this.mbr = mbr;
		this.nor = nor;

		try {
			if (!this.exportExpDesign()) {
				isDone = true;
				isRunning = false;
				return false;
			}

			if (isCancelled) {
				isDone = true;
				isRunning = false;
				return false;
			}

			if (!this.flashLFQ()) {
				isDone = true;
				isRunning = false;
				return false;
			}

			isDone = true;
			isRunning = false;
			return true;
		} catch (Exception e) {
			LOGGER.error(taskName + ": Unexpected error during execution", e);
			isDone = true;
			isRunning = false;
			return false;
		}
	}

	/**
	 * FlashLFQ CMD version has a bug that if the raw file folder contains other raw
	 * files which are not shown in the ExperimentalDesign.tsv, an Exception will
	 * thrown
	 * 
	 * @param rawDir
	 * @throws FileNotFoundException
	 */
	private boolean exportExpDesign() {

		File output = new File(this.rep, "ExperimentalDesign.tsv");
		PrintWriter writer = null;
		try {
			writer = new PrintWriter(output);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			LOGGER.error(taskName + ": error in writing experimental design file", e);
			System.out.println(
					format.format(new Date()) + "\t" + taskName + ": error in writing experimental design file");
			return false;
		}
		StringBuilder title = new StringBuilder();
		title.append("FileName").append("\t");
		title.append("Condition").append("\t");
		title.append("Biorep").append("\t");
		title.append("Fraction").append("\t");
		title.append("Techrep");

		writer.println(title);

		String[] raws = metadata.getRawFiles();
		String[] expNames = metadata.getExpNames();
		int[] fractions = metadata.getFractions();
		int[] replicates = metadata.getReplicates();

		HashSet<String> rawSet = new HashSet<String>();
		for (int i = 0; i < raws.length; i++) {
			StringBuilder sb = new StringBuilder();
			File rawFile = new File(raws[i]);
			String fileName = rawFile.getName();
			int dot = fileName.lastIndexOf('.');
			String name = dot > 0 ? fileName.substring(0, dot) : fileName;
			sb.append(name).append("\t");
			sb.append(expNames[i]).append("\t");

			// biorep
			sb.append("1").append("\t");

			sb.append(fractions[i]).append("\t");
			sb.append(replicates[i]);

			writer.println(sb);

			rawSet.add(fileName);
		}

		writer.close();

		File rawDir = new File(this.rep);
		File[] listFiles = rawDir.listFiles(new FilenameFilter() {

			@Override
			public boolean accept(File dir, String name) {
				// TODO Auto-generated method stub

				if (name.endsWith(".raw") || name.endsWith(".RAW") || name.endsWith(".mzML") || name.endsWith(".MZML")
						|| name.endsWith(".mzml"))
					return true;

				return false;
			}
		});

		if (listFiles.length != raws.length) {
			LOGGER.info(taskName + ": the number of the raw files in " + rawDir
					+ " didn't exactly matched to the expected raw file number in this task, " + raws.length
					+ " raw files were added in this task but " + listFiles.length + " raw files in this folder "
					+ rawDir);
			System.out.println(format.format(new Date()) + "\t" + taskName + ": the number of the raw files in "
					+ rawDir + " didn't exactly matched to the expected raw file number in this task, " + raws.length
					+ " raw files were added in this task but " + listFiles.length + " raw files in this folder "
					+ rawDir);

			return false;
		}

		for (int i = 0; i < listFiles.length; i++) {
			if (!rawSet.contains(listFiles[i].getName())) {
				LOGGER.info(taskName + ": this raw file " + listFiles[i].getName()
						+ " was not in this task, please move it out of this folder");
				System.out.println(format.format(new Date()) + "\t" + taskName + ": this raw file "
						+ listFiles[i].getName() + " was not in this task, please move it out of this folder");

				return false;
			}
		}

		LOGGER.info(taskName + ": the experimental design has been exported to " + output);
		System.out.println(format.format(new Date()) + "\t" + taskName
				+ ": the experimental design has been exported to " + output);

		return true;
	}

	private boolean flashLFQ() {
		File outFile = (new File(idt)).getParentFile();
		boolean windows = isWindowsOs();
		File artifact;
		try {
			artifact = resolveFlashLfqArtifact(new File(flashLFQ), windows);
		} catch (IllegalArgumentException e) {
			LOGGER.error(taskName + ": " + e.getMessage());
			return false;
		}
		if (!artifact.isFile()) {
			LOGGER.error(taskName + ": FlashLFQ artifact not found: " + artifact.getAbsolutePath()
					+ " (configured FLASHLFQ=" + flashLFQ + ")");
			System.out.println(format.format(new Date()) + "\t" + taskName + ": FlashLFQ artifact not found: "
					+ artifact.getAbsolutePath());
			return false;
		}
		if (!windows && artifact.getName().toLowerCase(Locale.ROOT).endsWith(".exe")) {
			LOGGER.error(taskName + ": FlashLFQ CMD.exe cannot run natively on "
					+ System.getProperty("os.name") + "; set FLASHLFQ to CMD.dll (or the FlashLFQ folder) "
					+ "and ensure DOTNET_ROOT or dotnet is on PATH");
			System.out.println(format.format(new Date()) + "\t" + taskName
					+ ": FlashLFQ CMD.exe cannot run natively on this OS; use CMD.dll with dotnet");
			return false;
		}

		if (windows) {
			writeWindowsBatHint(outFile, artifact);
		}

		String helpUrl = dotnetHelpUrl(windows);
		try {
			List<String> list = buildFlashLfqCommand(artifact, windows, idt, rep, outFile.getAbsolutePath(), mbr, nor,
					thr, resolveDotnetCommand());
			String[] args = list.toArray(new String[0]);

			StringBuilder sb = new StringBuilder();
			for (int i = 0; i < args.length; i++) {
				if (i > 0) {
					sb.append(' ');
				}
				sb.append(args[i]);
			}

			System.out.println(format.format(new Date()) + "\t" + sb);
			LOGGER.info(taskName + ": " + sb);

			ProcessBuilder pb = new ProcessBuilder(args);
			File workDir = artifact.getParentFile();
			if (workDir != null && workDir.isDirectory()) {
				pb.directory(workDir);
			}
			pb.redirectErrorStream(true);
			flashLFQProcess = pb.start();

			BufferedInputStream in = new BufferedInputStream(flashLFQProcess.getInputStream());
			BufferedReader inBr = new BufferedReader(new InputStreamReader(in));
			String lineStr;
			while ((lineStr = inBr.readLine()) != null) {
				if (isCancelled) {
					LOGGER.info(taskName + ": Task cancelled, terminating FlashLFQ process");
					System.out.println(format.format(new Date()) + "\t" + taskName
							+ ": Task cancelled, terminating FlashLFQ process");
					flashLFQProcess.destroyForcibly();
					inBr.close();
					in.close();
					return false;
				}

				LOGGER.info(lineStr);
				System.out.println(format.format(new Date()) + "\t" + taskName + ": " + lineStr);
			}

			int exitCode = flashLFQProcess.waitFor();
			if (exitCode != 0) {
				LOGGER.error(taskName + ": Process exited with code " + exitCode);
				System.out.println(
						format.format(new Date()) + "\t" + taskName + ": Process exited with code " + exitCode);

				LOGGER.info(taskName + ": try to fix this problem by installing the .NET runtime, please download from "
						+ helpUrl);
				System.out.println(format.format(new Date()) + "\t" + taskName
						+ ": try to fix this problem by installing the .NET runtime, please download from " + helpUrl);

				inBr.close();
				in.close();
				return false;
			}
			inBr.close();
			in.close();
		} catch (Exception e) {
			LOGGER.error(taskName + ": error in quantifying PSMs from " + idt, e);
			System.out.println(format.format(new Date()) + "\t" + taskName + ": error in quantifying PSMs from " + idt);

			LOGGER.info(taskName + ": try to fix this problem by installing the .NET runtime, please download from "
					+ helpUrl);
			System.out.println(format.format(new Date()) + "\t" + taskName
					+ ": try to fix this problem by installing the .NET runtime, please download from " + helpUrl);

			return false;
		}

		LOGGER.info(taskName + ": quantifying PSMs from " + idt + " finished");
		System.out
				.println(format.format(new Date()) + "\t" + taskName + ": quantifying PSMs from " + idt + " finished");

		return true;
	}

	/** Best-effort LFQ.bat for Windows debugging; not used to launch the process. */
	private void writeWindowsBatHint(File outFile, File artifact) {
		File batFile = new File(outFile, "LFQ.bat");
		try {
			String abs = artifact.getAbsolutePath();
			int id = abs.indexOf(':');
			if (id < 0) {
				LOGGER.warn(taskName + ": skipping LFQ.bat (no Windows drive letter in " + abs + ")");
				return;
			}
			PrintWriter writer = new PrintWriter(batFile);
			writer.println(abs.substring(0, id) + ":");
			writer.println("cd " + artifact.getParent());
			writer.println("CMD.exe --idt " + "\"" + idt + "\"" + " --rep " + "\"" + rep + "\"" + " --out " + "\""
					+ outFile.getAbsolutePath() + "\"" + " --mbr " + mbr + " --nor " + nor + " --thr " + thr);
			writer.print("exit");
			writer.close();
			LOGGER.info(artifact.getAbsolutePath() + " --idt \"" + idt + "\" --rep \"" + rep + "\" --out \""
					+ outFile.getAbsolutePath() + "\" --mbr " + mbr + " --nor " + nor + " --thr " + thr);
		} catch (FileNotFoundException e) {
			LOGGER.error(taskName + ": error in writing the command file " + batFile, e);
			System.out.println(
					format.format(new Date()) + "\t" + taskName + ": error in writing the command file " + batFile);
		}
	}

	protected static void exportExpDesign(File repFile) throws FileNotFoundException {

		File output = new File(repFile, "ExperimentalDesign.tsv");
		PrintWriter writer = new PrintWriter(output);

		StringBuilder title = new StringBuilder();
		title.append("FileName").append("\t");
		title.append("Condition").append("\t");
		title.append("Biorep").append("\t");
		title.append("Fraction").append("\t");
		title.append("Techrep");

		writer.println(title);

		HashSet<String> rawSet = new HashSet<String>();

		File[] rawFiles = repFile.listFiles(new FileFilter() {

			@Override
			public boolean accept(File pathname) {
				// TODO Auto-generated method stub
				if (pathname.getName().endsWith("raw"))
					return true;
				return false;
			}
		});
		for (int i = 0; i < rawFiles.length; i++) {
			StringBuilder sb = new StringBuilder();
			String name = rawFiles[i].getName().substring(0, rawFiles[i].getName().lastIndexOf("."));

			sb.append(name).append("\t");
			sb.append(name).append("\t");

			// biorep
			sb.append("1").append("\t");

			sb.append("1").append("\t");
			sb.append("1");

			writer.println(sb);

			rawSet.add(name + ".raw");
		}

		writer.close();
	}

	/**
	 * Check if the task is done (completed, failed, or cancelled)
	 * @return true if the task has finished execution
	 */
	public boolean isDone() {
		return isDone;
	}
	
	/**
	 * Check if the task is currently running
	 * @return true if the task is currently executing
	 */
	public boolean isRunning() {
		return isRunning;
	}
	
	/**
	 * Cancel the running task
	 * This will terminate the FlashLFQ process if it's running
	 */
	public void cancel() {
		if (isRunning && !isCancelled) {
			LOGGER.info(taskName + ": Cancellation requested");
			System.out.println(format.format(new Date()) + "\t" + taskName + ": Cancellation requested");
			
			isCancelled = true;
			
			// Terminate the FlashLFQ process if it's running
			if (flashLFQProcess != null && flashLFQProcess.isAlive()) {
				LOGGER.info(taskName + ": Terminating FlashLFQ process");
				System.out.println(format.format(new Date()) + "\t" + taskName + ": Terminating FlashLFQ process");
				
				flashLFQProcess.destroyForcibly();
				
				try {
					// Wait a bit for the process to terminate
					if (!flashLFQProcess.waitFor(5, TimeUnit.SECONDS)) {
						LOGGER.warn(taskName + ": FlashLFQ process did not terminate gracefully within 5 seconds");
						System.out.println(format.format(new Date()) + "\t" + taskName + ": FlashLFQ process did not terminate gracefully within 5 seconds");
					}
				} catch (InterruptedException e) {
					LOGGER.warn(taskName + ": Interrupted while waiting for process termination", e);
					Thread.currentThread().interrupt();
				}
			}
			
			isDone = true;
			isRunning = false;
		}
	}

}
