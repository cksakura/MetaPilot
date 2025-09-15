package uk.ac.quadram.metapilot.dbSearch.sage;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import javax.swing.JProgressBar;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class SageTask {

	private static String taskName = "Sage task";
	private static final Logger LOGGER = LogManager.getLogger(SageTask.class);
	private SimpleDateFormat format = new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss");

	private File sageFile;
	private ArrayList<String> parList;
	private StringBuilder headsb;
	private StringBuilder tailsb;

	public SageTask(String sage) {
		this(new File(sage));
	}

	public SageTask(File sage) {
		this.sageFile = sage;
		this.parList = new ArrayList<String>();
		this.headsb = new StringBuilder();
		this.headsb.append("@echo off\n");
		this.headsb.append("cd /d \"" + sageFile.getParent() + "\"\n");

		this.tailsb = new StringBuilder();
		this.tailsb.append("if errorlevel 1 (\n");
		this.tailsb.append("\techo Error: Failed to start the application.\n");
		this.tailsb.append(") else (\n");
		this.tailsb.append("\techo Application started successfully.\n");
		this.tailsb.append(")\n");
		this.tailsb.append("exit");
	}

	public void addTask(File configFile) {
		File batFile = new File(configFile.getParentFile(), "sage.bat");
		try (PrintWriter writer = new PrintWriter(batFile)) {
			writer.println(headsb);
			writer.println("start /B " + sageFile.getName() + " " + configFile.getAbsolutePath());
			writer.println();
			writer.print("exit");
			writer.close();
		} catch (IOException e) {
			LOGGER.error(taskName + ": error in writing the .bat file to " + batFile, e);
			System.out.println(
					format.format(new Date()) + "\t" + taskName + ": error in writing the .bat file to " + batFile);
		}
		this.parList.add(batFile.getAbsolutePath());
	}

	/**
	 * Search mzML files against a FASTA database using Sage
	 * 
	 * @param mzMLFiles Array of mzML file paths
	 * @param result Output directory for results
	 * @param db FASTA database file path
	 * @param precursorTol Precursor tolerance in ppm
	 * @param fragmentTol Fragment tolerance in ppm
	 * @param enzyme Enzyme name
	 * @param missedCleavages Number of missed cleavages
	 * @param minPepLen Minimum peptide length
	 * @param maxPepLen Maximum peptide length
	 * @param fdr FDR threshold
	 * @param fixedMods Array of fixed modifications
	 * @param variableMods Array of variable modifications
	 */
	public void searchMzML(String[] mzMLFiles, String result, String db, double precursorTol, double fragmentTol,
			String enzyme, int missedCleavages, int minPepLen, int maxPepLen, double fdr, 
			String[] fixedMods, String[] variableMods) {
		
		File resultDir = new File(result);
		if (!resultDir.exists()) {
			resultDir.mkdirs();
		}

		File configFile = new File(resultDir, "sage_config.json");
		
		try (PrintWriter writer = new PrintWriter(configFile)) {
			writer.println("{");
			writer.println("  \"database\": {");
			writer.println("    \"fasta\": \"" + db.replace("\\", "\\\\") + "\",");
			writer.println("    \"decoy_tag\": \"rev_\",");
			writer.println("    \"generate_decoys\": true");
			writer.println("  },");
			writer.println("  \"precursor_tol\": {");
			writer.println("    \"ppm\": [0, " + precursorTol + "]");
			writer.println("  },");
			writer.println("  \"fragment_tol\": {");
			writer.println("    \"ppm\": [0, " + fragmentTol + "]");
			writer.println("  },");
			writer.println("  \"enzyme\": {");
			writer.println("    \"missed_cleavages\": " + missedCleavages + ",");
			writer.println("    \"min_len\": " + minPepLen + ",");
			writer.println("    \"max_len\": " + maxPepLen);
			writer.println("  },");
			
			// Add fixed modifications
			writer.println("  \"static_mods\": {");
			if (fixedMods != null && fixedMods.length > 0) {
				for (int i = 0; i < fixedMods.length; i++) {
					// Parse modification format: "Carbamidomethyl[C] 0" -> "C": 57.0214
					String mod = fixedMods[i];
					if (mod != null && !mod.trim().isEmpty()) {
						String[] parts = mod.split("\\[|\\]");
						if (parts.length >= 2) {
							String aminoAcid = parts[1]; // Extract amino acid from [C]
							String massShift = getModificationMass(parts[0]); // Get mass for modification name
							writer.print("    \"" + aminoAcid + "\": " + massShift);
							if (i < fixedMods.length - 1) writer.println(",");
							else writer.println();
						}
					}
				}
			}
			writer.println("  },");
			
			// Add variable modifications
			writer.println("  \"variable_mods\": {");
			if (variableMods != null && variableMods.length > 0) {
				for (int i = 0; i < variableMods.length; i++) {
					// Parse modification format: "Oxidation[M] 0" -> "M": 15.9949
					String mod = variableMods[i];
					if (mod != null && !mod.trim().isEmpty()) {
						String[] parts = mod.split("\\[|\\]");
						if (parts.length >= 2) {
							String aminoAcid = parts[1]; // Extract amino acid from [M]
							String massShift = getModificationMass(parts[0]); // Get mass for modification name
							writer.print("    \"" + aminoAcid + "\": " + massShift);
							if (i < variableMods.length - 1) writer.println(",");
							else writer.println();
						}
					}
				}
			}
			writer.println("  },");
			
			writer.println("  \"max_variable_mods\": 2,");
			writer.println("  \"output_directory\": \"" + result.replace("\\", "\\\\") + "\",");
			writer.print("  \"mzml_paths\": [");
			for (int i = 0; i < mzMLFiles.length; i++) {
				writer.print("\"" + mzMLFiles[i].replace("\\", "\\\\") + "\"");
				if (i < mzMLFiles.length - 1) writer.print(", ");
			}
			writer.println("],");
			writer.println("  \"write_checkpoint\": false,");
			writer.println("  \"fdr\": " + fdr);
			writer.println("}");
			
		} catch (IOException e) {
			LOGGER.error(taskName + ": error in writing Sage config file to " + configFile, e);
			System.out.println(format.format(new Date()) + "\t" + taskName + ": error in writing Sage config file to " + configFile);
			return;
		}
		
		addTask(configFile);
	}

	/**
	 * Search mzML files against a FASTA database using Sage with enhanced parameters
	 * 
	 * @param mzMLFiles Array of mzML file paths
	 * @param result Output directory for results
	 * @param db FASTA database file path
	 * @param precursorTolMin Minimum precursor tolerance in ppm
	 * @param precursorTolMax Maximum precursor tolerance in ppm
	 * @param fragmentTol Fragment tolerance in ppm
	 * @param enzyme Enzyme name
	 * @param missedCleavages Number of missed cleavages
	 * @param fixedMods Array of fixed modifications
	 * @param variableMods Array of variable modifications
	 * @param isQuant Whether to include quantification settings
	 */
	public void searchMzML(String[] mzMLFiles, String result, String db, double precursorTolMin, double precursorTolMax, double fragmentTol,
			String enzyme, int missedCleavages, String[] fixedMods, String[] variableMods, boolean isQuant) {
		
		File resultDir = new File(result);
		if (!resultDir.exists()) {
			resultDir.mkdirs();
		}

		File configFile = new File(resultDir, "sage_config.json");
		
		try (PrintWriter writer = new PrintWriter(configFile)) {
			writer.println("{");
			writer.println("  \"database\": {");
			writer.println("    \"fasta\": \"" + db.replace("\\", "\\\\") + "\",");
			writer.println("    \"decoy_tag\": \"rev_\",");
			writer.println("    \"generate_decoys\": true");
			writer.println("  },");
			writer.println("  \"precursor_tol\": {");
			writer.println("    \"ppm\": [" + precursorTolMin + ", " + precursorTolMax + "]");
			writer.println("  },");
			writer.println("  \"fragment_tol\": {");
			writer.println("    \"ppm\": [0, " + fragmentTol + "]");
			writer.println("  },");
			writer.println("  \"enzyme\": {");
			writer.println("    \"missed_cleavages\": " + missedCleavages);
			writer.println("  },");
			
			// Add fixed modifications
			writer.println("  \"static_mods\": {");
			if (fixedMods != null && fixedMods.length > 0) {
				for (int i = 0; i < fixedMods.length; i++) {
					// Parse modification format: "C:57.0214" -> "C": 57.0214
					String mod = fixedMods[i];
					if (mod != null && !mod.trim().isEmpty()) {
						String[] parts = mod.split(":");
						if (parts.length >= 2) {
							String aminoAcid = parts[0]; // Extract amino acid
							String massShift = parts[1]; // Extract mass value
							writer.print("    \"" + aminoAcid + "\": " + massShift);
							if (i < fixedMods.length - 1) writer.println(",");
							else writer.println();
						}
					}
				}
			}
			writer.println("  },");
			
			// Add variable modifications
			writer.println("  \"variable_mods\": {");
			if (variableMods != null && variableMods.length > 0) {
				for (int i = 0; i < variableMods.length; i++) {
					// Parse modification format: "Oxidation:M:15.9949" -> "M": 15.9949
					String mod = variableMods[i];
					if (mod != null && !mod.trim().isEmpty()) {
						String[] parts = mod.split(":");
						if (parts.length >= 3) {
							String aminoAcid = parts[1]; // Extract amino acid from middle part
							String massShift = parts[2]; // Extract mass value from last part
							writer.print("    \"" + aminoAcid + "\": " + massShift);
							if (i < variableMods.length - 1) writer.println(",");
							else writer.println();
						}
					}
				}
			}
			writer.println("  },");
			
			writer.println("  \"max_variable_mods\": 2,");
			writer.println("  \"output_directory\": \"" + result.replace("\\", "\\\\") + "\",");
			writer.print("  \"mzml_paths\": [");
			for (int i = 0; i < mzMLFiles.length; i++) {
				writer.print("\"" + mzMLFiles[i].replace("\\", "\\\\") + "\"");
				if (i < mzMLFiles.length - 1) writer.print(", ");
			}
			writer.println("]");
			
			// Add quantification settings if isQuant is true
			if (isQuant) {
				writer.println(",");
				writer.println("  \"quant\": {");
				writer.println("    \"lfq\": true,");
				writer.println("    \"lfq_settings\": {");
				writer.println("      \"peak_scoring\": \"Hybrid\",");
				writer.println("      \"integration\": \"Sum\",");
				writer.println("      \"spectral_angle\": 0.6,");
				writer.println("      \"ppm_tolerance\": 5.0");
				writer.println("    }");
				writer.println("  }");
			}
			
			writer.println("}");
			
		} catch (IOException e) {
			LOGGER.error(taskName + ": error in writing Sage config file to " + configFile, e);
			System.out.println(format.format(new Date()) + "\t" + taskName + ": error in writing Sage config file to " + configFile);
			return;
		}
		
		addTask(configFile);
	}

	/**
	 * Get the mass shift for a modification name
	 */
	private String getModificationMass(String modificationName) {
		switch (modificationName.trim()) {
			case "Carbamidomethyl":
				return "57.0214";
			case "Oxidation":
				return "15.9949";
			case "Acetyl":
				return "42.0106";
			case "Phospho":
				return "79.9663";
			case "Methyl":
				return "14.0157";
			case "Dimethyl":
				return "28.0314";
			case "Trimethyl":
				return "42.0471";
			case "Deamidated":
				return "0.9840";
			case "Pyro-glu":
				return "-17.0265";
			case "Pyro_glu":
				return "-17.0265";
			default:
				LOGGER.warn("Unknown modification: " + modificationName + ", using default mass 0");
				return "0";
		}
	}

	/**
	 * Set total thread count for parallel execution
	 */
	public void setTotalThread(int threadCount) {
		// Sage handles threading internally, but we can use this for task scheduling
		LOGGER.info("Set thread count to: " + threadCount);
	}

	/**
	 * Clean Sage result folder by removing temporary files
	 */
	public static void cleanSageResultFolder(File folder, boolean keepResults) {
		if (!folder.exists()) {
			return;
		}
		
		File[] files = folder.listFiles();
		if (files != null) {
			for (File file : files) {
				if (file.isFile()) {
					// Keep important result files, remove temporary ones
					String name = file.getName();
					if (!keepResults || (!name.endsWith(".sage.tsv") && !name.endsWith(".sage.parquet"))) {
						if (name.endsWith(".bat") || name.endsWith(".json") || name.endsWith(".log")) {
							file.delete();
						}
					}
				}
			}
		}
	}

	public File getSageFile() {
		return sageFile;
	}

	public void run(int waitHours) {
		run(1, waitHours);
	}

	public void run(int threadPool, int waitHours) {

		if (this.parList.size() == 0) {
			return;
		}

		ExecutorService executor = Executors.newFixedThreadPool(threadPool);
		for (int i = 0; i < parList.size(); i++) {

			String pari = parList.get(i);

			LOGGER.info(pari);

			executor.submit(new Runnable() {

				@Override
				public void run() {
					// TODO Auto-generated method stub

					try {

						ArrayList<String> commands = new ArrayList<String>();
						commands.add("cmd.exe");
						commands.add("/c");
						commands.add("start");
						commands.add("/min");

						String[] cs = pari.split("\t");
						for (String parString : cs) {
							commands.add(parString);
						}
						ProcessBuilder pb = new ProcessBuilder(commands);
						Process p = pb.start();
						BufferedInputStream in = new BufferedInputStream(p.getInputStream());
						BufferedReader inBr = new BufferedReader(new InputStreamReader(in));
						String line0 = "";
						String line1 = "";
						while ((line1 = inBr.readLine()) != null) {
							if (!line1.equals(line0)) {
								System.out.println(format.format(new Date()) + "\t" + line1);
								line0 = line1;
							}
						}

						if (p.waitFor() != 0) {
							if (p.exitValue() == 1)
								LOGGER.error("false");
						}
						inBr.close();
						in.close();

					} catch (IOException | InterruptedException e) {
						// TODO Auto-generated catch block
						LOGGER.error(e);
					}
				}
			});
		}

		try {

			executor.shutdown();

			boolean finish = executor.awaitTermination(waitHours, TimeUnit.HOURS);

			this.parList = new ArrayList<String>();

			if (finish) {
				System.out.println(format.format(new Date()) + "\t" + taskName + ": finished");
			} else {
				LOGGER.info(taskName
						+ ": task dosen't finish in a long time, please restart MetaLab after the Sage task finish");
				System.out.println(format.format(new Date()) + "\t" + taskName
						+ ": task dosen't finish in a long time, please restart MetaLab after the Sage task finish");
			}

		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			LOGGER.error(taskName + ": failed", e);
			System.out.println(format.format(new Date()) + "\t" + taskName + ": failed");
		}
	}

	public boolean runCMD(String cmd, JProgressBar bar) {

		try {

			int id = sageFile.getAbsolutePath().indexOf(":");
			if (id < 0) {
				if (id < 0) {
					LOGGER.error(taskName + ": cann't find the directory of sage.exe");
					System.out.println(
							format.format(new Date()) + "\t" + taskName + ": cann't find the directory of sage.exe");

					return false;
				}
			}

			PrintWriter writer = new PrintWriter(cmd);
			writer.println(sageFile.getAbsolutePath().substring(0, id) + ":");
			writer.println("cd " + sageFile.getParent());

			for (int i = 0; i < parList.size(); i++) {
				writer.println(sageFile.getName() + " " + parList.get(i));
			}

			writer.print("exit");
			writer.close();

		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			LOGGER.error(taskName + ": error in writing the executable command for sage.exe to " + cmd, e);
			System.out.println(format.format(new Date()) + "\t" + taskName
					+ ": error in writing the executable command for sage.exe to " + cmd);

			return false;
		}

		try {
			ArrayList<String> cmdList = new ArrayList<String>();
			cmdList.add("cmd.exe");
			cmdList.add("/c");
			cmdList.add("start");
			cmdList.add(cmd);

			ProcessBuilder pb = new ProcessBuilder(cmdList);
			Process p = pb.start();

			BufferedInputStream in = new BufferedInputStream(p.getInputStream());
			BufferedReader inBr = new BufferedReader(new InputStreamReader(in));
			String lineStr;
			while ((lineStr = inBr.readLine()) != null) {
				LOGGER.info(lineStr);
				System.out.println(format.format(new Date()) + "\t" + taskName + ": " + lineStr);
			}

			if (p.waitFor() != 0) {
				if (p.exitValue() == 1)
					LOGGER.error("false");

				bar.setValue(100);
			}
			inBr.close();
			in.close();
		} catch (Exception e) {
			LOGGER.error(taskName + ": error in running the executable command by sage.exe", e);
			System.out.println(format.format(new Date()) + "\t" + taskName
					+ ": error in running the executable command by sage.exe");

			return false;
		}

		return true;
	}

}
