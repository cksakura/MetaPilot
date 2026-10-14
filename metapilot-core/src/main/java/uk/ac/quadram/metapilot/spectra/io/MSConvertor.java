/**
 * 
 */
package uk.ac.quadram.metapilot.spectra.io;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.function.Consumer;

/**
 * @author Kai Cheng
 *
 */
public class MSConvertor {
	
	private static String tofString = "--filter \"analyzer TOF\"";
	private static String polarityString = "--filter \"polarity  positive\"";
	private static String scanSumString = "--filter \"scanSumming precursorTol=0.05 scanTimeTol=5 ionMobilityTol=0.1\"";
	private static String titleMakerString = "--filter \"titleMaker <RunId>.<ScanNumber>.<ScanNumber>.<ChargeState>\"";

	private static final String convertCmd = "Resources//ProteoWizard//msconvert.exe";
	
	public static void d2mgf(String convertCmd, String raw, String resultDir) {

		File file = new File(raw);
		String rawName = file.getName();
		rawName = rawName.substring(0, rawName.lastIndexOf("."));

		File batFile = new File(resultDir, rawName + ".bat");
		File msconvertFile = new File(convertCmd);
		int id = convertCmd.indexOf(":");
		if (id < 0) {
			return;
		}

		try {

			PrintWriter writer = new PrintWriter(batFile);
			writer.println(convertCmd.substring(0, id) + ":");
			writer.println("cd " + msconvertFile.getParent());
			writer.println("msconvert.exe " + "\"" + file + "\"" + " --mgf " + "-o " + "\"" + resultDir + "\" "
					+ tofString + " " + polarityString + " " + scanSumString + " " + titleMakerString);
			writer.print("exit");
			writer.close();

			ArrayList<String> list = new ArrayList<String>();
			list.add("cmd.exe");
			list.add("/c");
			list.add("start");
			list.add(batFile.getAbsolutePath());
			Process p = Runtime.getRuntime().exec(list.toArray(new String[list.size()]));

			BufferedInputStream in = new BufferedInputStream(p.getInputStream());
			BufferedReader inBr = new BufferedReader(new InputStreamReader(in));
			String lineStr;
			while ((lineStr = inBr.readLine()) != null)
				System.out.println(lineStr);
			if (p.waitFor() != 0) {
				if (p.exitValue() == 1)
					System.err.println("false");
			}
			inBr.close();
			in.close();

		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public static String d2mzML(String convertCmd, String raw, String resultDir) {

		File file = new File(raw);
		String rawName = file.getName();
		rawName = rawName.substring(0, rawName.lastIndexOf("."));

		File batFile = new File(resultDir, rawName + ".bat");
		File msconvertFile = new File(convertCmd);
		int id = convertCmd.indexOf(":");
		if (id < 0) {
			return null;
		}

		try {

			PrintWriter writer = new PrintWriter(batFile);
			writer.println(convertCmd.substring(0, id) + ":");
			writer.println("cd " + msconvertFile.getParent());
			writer.println("msconvert.exe " + "\"" + file + "\"" + " --mzML " + "-o " + "\"" + resultDir + "\" "
					+ tofString + " " + polarityString + " " + scanSumString + " " + titleMakerString);
			writer.print("exit");
			writer.close();

			ArrayList<String> list = new ArrayList<String>();
			list.add("cmd.exe");
			list.add("/c");
			list.add("start");
			list.add(batFile.getAbsolutePath());
			Process p = Runtime.getRuntime().exec(list.toArray(new String[list.size()]));

			BufferedInputStream in = new BufferedInputStream(p.getInputStream());
			BufferedReader inBr = new BufferedReader(new InputStreamReader(in));
			String lineStr;
			while ((lineStr = inBr.readLine()) != null)
				System.out.println(lineStr);
			if (p.waitFor() != 0) {
				if (p.exitValue() == 1)
					System.err.println("false");
			}
			inBr.close();
			in.close();

		} catch (Exception e) {
			e.printStackTrace();
			return null;
		}

		String name = file.getName();
		name = name.substring(0, name.lastIndexOf(".")) + ".mzML";
		return new File(resultDir, name).getAbsolutePath();
	}
	
	
	public static void rawMzml2mgf(String convertCmd, String raw, String resultDir) {

		File file = new File(raw);
		String rawName = file.getName();
		rawName = rawName.substring(0, rawName.lastIndexOf("."));

		File batFile = new File(resultDir, rawName + ".bat");
		File msconvertFile = new File(convertCmd);
		int id = convertCmd.indexOf(":");
		if (id < 0) {
			return;
		}

		try {

			PrintWriter writer = new PrintWriter(batFile);
			writer.println(convertCmd.substring(0, id) + ":");
			writer.println("cd " + msconvertFile.getParent());
			writer.println("msconvert.exe " + "\"" + file + "\"" + " --mgf " + "-o " + "\"" + resultDir + "\" "
					+ polarityString + " " + scanSumString + " " + titleMakerString);
			writer.print("exit");
			writer.close();

			ArrayList<String> list = new ArrayList<String>();
			list.add("cmd.exe");
			list.add("/c");
			list.add("start");
			list.add(batFile.getAbsolutePath());
			Process p = Runtime.getRuntime().exec(list.toArray(new String[list.size()]));

			BufferedInputStream in = new BufferedInputStream(p.getInputStream());
			BufferedReader inBr = new BufferedReader(new InputStreamReader(in));
			String lineStr;
			while ((lineStr = inBr.readLine()) != null)
				System.out.println(lineStr);
			if (p.waitFor() != 0) {
				if (p.exitValue() == 1)
					System.err.println("false");
			}
			inBr.close();
			in.close();

		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public static void tofMzml2mgf(String convertCmd, String raw, String resultDir) {

		File file = new File(raw);
		String rawName = file.getName();
		rawName = rawName.substring(0, rawName.lastIndexOf("."));

		File batFile = new File(resultDir, rawName + ".bat");
		File msconvertFile = new File(convertCmd);
		int id = convertCmd.indexOf(":");
		if (id < 0) {
			return;
		}

		try {

			PrintWriter writer = new PrintWriter(batFile);
			writer.println(convertCmd.substring(0, id) + ":");
			writer.println("cd " + msconvertFile.getParent());
			writer.println("msconvert.exe " + "\"" + file + "\"" + " --mgf " + "-o " + "\"" + resultDir + "\" "
					+ tofString + " " + polarityString + " " + scanSumString + " " + titleMakerString);
			writer.print("exit");
			writer.close();

			ArrayList<String> list = new ArrayList<String>();
			list.add("cmd.exe");
			list.add("/c");
			list.add("start");
			list.add(batFile.getAbsolutePath());
			Process p = Runtime.getRuntime().exec(list.toArray(new String[list.size()]));

			BufferedInputStream in = new BufferedInputStream(p.getInputStream());
			BufferedReader inBr = new BufferedReader(new InputStreamReader(in));
			String lineStr;
			while ((lineStr = inBr.readLine()) != null)
				System.out.println(lineStr);
			if (p.waitFor() != 0) {
				if (p.exitValue() == 1)
					System.err.println("false");
			}
			inBr.close();
			in.close();

		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public static void raw2mgf(String convertCmd, String raw, String resultDir) {

		File file = new File(raw);
		String cmd = convertCmd + " " + file + " --mgf " + "-o " + resultDir
				+ " --zlib --filter \"peakPicking true 1-\"";

		Runtime run = Runtime.getRuntime();
		try {
			Process p = run.exec(cmd);
			BufferedInputStream in = new BufferedInputStream(p.getInputStream());
			BufferedReader inBr = new BufferedReader(new InputStreamReader(in));
			String lineStr;
			while ((lineStr = inBr.readLine()) != null)
				System.out.println(lineStr);
			if (p.waitFor() != 0) {
				if (p.exitValue() == 1)
					System.err.println("false");
			}
			inBr.close();
			in.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public static File raw2mgf(String raw) {
		return raw2mgf(new File(raw));
	}

	public static File raw2mgf(File raw) {
		String parent = raw.getParent();
		String cmd = convertCmd + " " + raw + " --mgf " + "-o " + parent + " --zlib --filter \"peakPicking true 1-\"";

		Runtime run = Runtime.getRuntime();
		try {
			Process p = run.exec(cmd);
			BufferedInputStream in = new BufferedInputStream(p.getInputStream());
			BufferedReader inBr = new BufferedReader(new InputStreamReader(in));
			String lineStr;
			while ((lineStr = inBr.readLine()) != null)
				System.out.println(lineStr);
			if (p.waitFor() != 0) {
				if (p.exitValue() == 1)
					System.err.println("false");
			}
			inBr.close();
			in.close();
		} catch (Exception e) {
			e.printStackTrace();
		}

		String name = raw.getName();
		name = name.substring(0, name.lastIndexOf(".")) + ".mgf";
		return new File(parent, name);
	}

	public static String raw2Mzxml(String raw) {
		File file = new File(raw);
		String parent = file.getParent();
		String cmd = convertCmd + " " + file + " --mzXML " + "-o " + parent
				+ " --zlib --filter \"peakPicking true 1-\"";

		Runtime run = Runtime.getRuntime();
		try {
			Process p = run.exec(cmd);
			BufferedInputStream in = new BufferedInputStream(p.getInputStream());
			BufferedReader inBr = new BufferedReader(new InputStreamReader(in));
			String lineStr;
			while ((lineStr = inBr.readLine()) != null)
				System.out.println(lineStr);
			if (p.waitFor() != 0) {
				if (p.exitValue() == 1)
					System.err.println("false");
			}
			inBr.close();
			in.close();
		} catch (Exception e) {
			e.printStackTrace();
		}

		String name = file.getName();
		name = name.substring(0, name.lastIndexOf(".")) + ".mzXML";
		return parent + "\\" + name;
	}
	
	public static String raw2mzML(String raw) {
		File file = new File(raw);
		String parent = file.getParent();
		String cmd = convertCmd + " " + file + " --mzML " + "-o " + parent
				+ " --zlib --filter \"peakPicking true 1-\"";

		Runtime run = Runtime.getRuntime();
		try {
			Process p = run.exec(cmd);
			BufferedInputStream in = new BufferedInputStream(p.getInputStream());
			BufferedReader inBr = new BufferedReader(new InputStreamReader(in));
			String lineStr;
			while ((lineStr = inBr.readLine()) != null)
				System.out.println(lineStr);
			if (p.waitFor() != 0) {
				if (p.exitValue() == 1)
					System.err.println("false");
			}
			inBr.close();
			in.close();
		} catch (Exception e) {
			e.printStackTrace();
		}

		String name = file.getName();
		name = name.substring(0, name.lastIndexOf(".")) + ".mzML";
		return parent + "\\" + name;
	}
	
	public static String raw2mzML(String convertCmd, String raw) {
		File file = new File(raw);
		String parent = file.getParent();
		String cmd = convertCmd + " " + file + " --mzML " + "-o " + parent
				+ " --zlib --filter \"peakPicking true 1-\"";

		Runtime run = Runtime.getRuntime();
		try {
			Process p = run.exec(cmd);
			BufferedInputStream in = new BufferedInputStream(p.getInputStream());
			BufferedReader inBr = new BufferedReader(new InputStreamReader(in));
			String lineStr;
			while ((lineStr = inBr.readLine()) != null)
				System.out.println(lineStr);
			if (p.waitFor() != 0) {
				if (p.exitValue() == 1)
					System.err.println("false");
			}
			inBr.close();
			in.close();
		} catch (Exception e) {
			e.printStackTrace();
		}

		String name = file.getName();
		name = name.substring(0, name.lastIndexOf(".")) + ".mzML";
		return parent + "\\" + name;
	}

	/**
	 * The conversion {@link #raw2mzML(String, String)} runs, written into
	 * {@code outDir} instead of beside the raw file. Every line msconvert prints,
	 * errors included, goes to {@code output}, and the exit code is returned
	 * rather than swallowed.
	 *
	 * @return msconvert's exit code
	 */
	public static int raw2mzML(String convertCmd, File raw, File outDir, Consumer<String> output)
			throws IOException, InterruptedException {
		ProcessBuilder pb = new ProcessBuilder(raw2mzMLCommand(convertCmd, raw, outDir));
		pb.redirectErrorStream(true);
		Process p = pb.start();
		try (BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()))) {
			String line;
			while ((line = br.readLine()) != null) {
				if (Thread.currentThread().isInterrupted()) {
					throw new InterruptedException("Interrupted while converting " + raw.getName());
				}
				if (output != null && !line.isBlank()) {
					output.accept(line.trim());
				}
			}
			return p.waitFor();
		} finally {
			if (p.isAlive()) {
				p.destroyForcibly();
			}
		}
	}

	/**
	 * Arguments for {@link #raw2mzML(String, File, File, Consumer)}, one list
	 * element each, so paths containing spaces stay single arguments. Vendor
	 * centroiding of every MS level: Crystal-C and FlashLFQ read MS1 too.
	 */
	public static List<String> raw2mzMLCommand(String convertCmd, File raw, File outDir) {
		return Arrays.asList(convertCmd, raw.getAbsolutePath(), "--mzML", "-o", outDir.getAbsolutePath(), "--zlib",
				"--filter", "peakPicking true 1-");
	}

	/**
	 * The {@code ProteoWizard release: ...} line of msconvert's help, or null
	 * when it cannot be read.
	 */
	public static String version(String convertCmd) {
		try {
			ProcessBuilder pb = new ProcessBuilder(convertCmd, "--help");
			pb.redirectErrorStream(true);
			Process p = pb.start();
			String version = null;
			try (BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()))) {
				String line;
				while ((line = br.readLine()) != null) {
					if (version == null && line.trim().startsWith("ProteoWizard release:")) {
						version = line.trim();
					}
				}
			} finally {
				p.waitFor();
			}
			return version;
		} catch (IOException e) {
			return null;
		} catch (InterruptedException e) {
			Thread.currentThread().interrupt();
			return null;
		}
	}

	public static void help() {
		String cmd = convertCmd + " --help";
		Runtime run = Runtime.getRuntime();
		try {
			Process p = run.exec(cmd);
			BufferedInputStream in = new BufferedInputStream(p.getInputStream());
			BufferedReader inBr = new BufferedReader(new InputStreamReader(in));
			String lineStr;
			while ((lineStr = inBr.readLine()) != null)
				System.out.println(lineStr);
			if (p.waitFor() != 0) {
				if (p.exitValue() == 1)
					System.err.println("false");
			}
			inBr.close();
			in.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public static void main(String[] args) {
		File[] files = (new File("Z:\\PRIDE_dataset\\PXD051104_metaExpertPro\\MetaPilot_DDA_MPW")).listFiles();
		for (int i = 0; i < files.length; i++) {
			if (files[i].getName().endsWith(".raw")) {
				MSConvertor.raw2mzML("C:\\Users\\daz25hak\\MetaProAI\\Resources\\ProteoWizard\\msconvert.exe",
						files[i].getAbsolutePath());
			}
		}
	}
}
