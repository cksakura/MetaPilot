package uk.ac.quadram.metapilot.dbSearch.diann;

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

import uk.ac.quadram.metapilot.task.dia.DiaLibSearchPar;

public class DiaNNTask {

	private static String taskName = "DiaNN task";
	private static final Logger LOGGER = LogManager.getLogger(DiaNNTask.class);
	private SimpleDateFormat format = new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss");

	private File diaNNFile;
	private ArrayList<String> parList;
	private StringBuilder headsb;
	private StringBuilder tailsb;
	
	// Task state management
	private volatile boolean isRunning = false;
	private volatile boolean isCancelled = false;
	private volatile boolean isDone = false;
	private ExecutorService executor = null;

	public DiaNNTask(String diaNN) {
		this(new File(diaNN));
	}

	public DiaNNTask(File diaNN) {
		this.diaNNFile = diaNN;
		this.parList = new ArrayList<String>();
		this.headsb = new StringBuilder();
		this.headsb.append("@echo off\n");
		this.headsb.append("cd /d \"" + diaNNFile.getParent() + "\"\n");

		this.tailsb = new StringBuilder();
		this.tailsb.append("if errorlevel 1 (\n");
		this.tailsb.append("\techo Error: Failed to start the application.\n");
		this.tailsb.append(") else (\n");
		this.tailsb.append("\techo Application started successfully.\n");
		this.tailsb.append(")\n");
		this.tailsb.append("exit");
	}

	public File getDiaNNFile() {
		return diaNNFile;
	}

	/**
	 * convert raw files to dia files
	 * 
	 * @param par
	 * @param raws
	 */
	public void addTask(DiannParameter par, String[] raws) {
		addTask(par, raws, null);
	}
	
	/**
	 * convert raw files to dia files
	 * 
	 * @param par
	 * @param raws
	 * @param resultDir - directory where batch files should be created (if null, uses raw file directory)
	 */
	public void addTask(DiannParameter par, String[] raws, File resultDir) {
		String cmd = par.generateConvert();
		File targetDir = (resultDir != null) ? resultDir : (new File(raws[0])).getParentFile();
		File batFile = new File(targetDir, "convert.bat");
		File cfgFile = new File(targetDir, "convert.cfg");

		LOGGER.info("Creating DiaNN conversion files:");
		LOGGER.info("  Target directory: {}", targetDir.getAbsolutePath());
		LOGGER.info("  Batch file: {}", batFile.getAbsolutePath());
		LOGGER.info("  Config file: {}", cfgFile.getAbsolutePath());
		LOGGER.info("  DiaNN executable: {}", diaNNFile.getAbsolutePath());
		LOGGER.info("  Base command: {}", cmd);
		LOGGER.info("  Raw files to convert: {}", java.util.Arrays.toString(raws));

		try (PrintWriter writer = new PrintWriter(batFile)) {
			writer.println(headsb);
			/*
			 * writer.println("setlocal enabledelayedexpansion"); writer.println(); writer.
			 * println("REM Get all the file paths passed as parameters and store them in a variable"
			 * ); writer.println("set files="); writer.println(); writer.println(
			 * "REM Loop through each file path and append it to the variable with a space and a --f option"
			 * ); writer.println("for %%f in (%*) do (");
			 * writer.println("\tset files=!files! --f \"%%f\""); writer.println(")");
			 * writer.println();
			 * writer.println("REM Print the content of the start command");
			 * writer.println("echo !files!"); writer.println(); writer.println("start /B "
			 * + diaNNFile.getName() + " " + cmd + "!files!");
			 * 
			 * 
			 * writer.println(); writer.println("endlocal");
			 */
			writer.println("start /B " + diaNNFile.getName() + " --cfg " + cfgFile.getAbsolutePath());
			writer.println();
			writer.print("exit");
			writer.close();

			PrintWriter cfgWriter = new PrintWriter(cfgFile);
			cfgWriter.print(cmd);
			for (int i = 0; i < raws.length; i++) {
				cfgWriter.print(" --f \"" + raws[i] + "\"");
			}
			cfgWriter.close();
		} catch (IOException e) {
			LOGGER.error(taskName + ": error in writing the .bat file to " + batFile, e);
			System.out.println(
					format.format(new Date()) + "\t" + taskName + ": error in writing the .bat file to " + batFile);
		}
		/*
		 * StringBuilder sb = new StringBuilder(); sb.append(batFile.getAbsolutePath());
		 * for (String raw : raws) { sb.append("\t").append(raw); }
		 * this.parList.add(sb.toString());
		 */
		this.parList.add(batFile.getAbsolutePath());
	}

	/**
	 * search library
	 * 
	 * @param par
	 * @param raws
	 * @param out
	 * @param lib
	 */
	public void addTask(DiannParameter par, String[] raws, String out, String... libs) {
		String cmd = par.generateLibSearch(raws, out, libs);
		File batFile = new File(out + ".bat");
		File cfgFile = new File(out + ".cfg");
		try (PrintWriter writer = new PrintWriter(batFile)) {
			writer.println(headsb);
			writer.println("start /B " + diaNNFile.getName() + " --cfg " + cfgFile.getAbsolutePath());
			writer.print("exit");
			writer.close();

			PrintWriter cfgWriter = new PrintWriter(cfgFile);
			cfgWriter.print(cmd);
			for (int i = 0; i < raws.length; i++) {
				cfgWriter.print(" --f \"" + raws[i] + "\"");
			}
			cfgWriter.close();

		} catch (IOException e) {
			LOGGER.error(taskName + ": error in writing the .bat file to " + batFile, e);
			System.out.println(
					format.format(new Date()) + "\t" + taskName + ": error in writing the .bat file to " + batFile);
		}

		this.parList.add(batFile.getAbsolutePath());
	}

	/**
	 * search library
	 * 
	 * @param par
	 * @param raws
	 * @param out
	 * @param lib
	 */
	public void addTask(DiannParameter par, String[] raws, String out, boolean fastSearch, String... lib) {
		String cmd;
		if (fastSearch) {
			cmd = par.generateLibSearchFast(raws, out, lib);
		} else {
			cmd = par.generateLibSearch(raws, out, lib);
		}
		File batFile = new File(out + ".bat");
		File cfgFile = new File(out + ".cfg");

		try (PrintWriter writer = new PrintWriter(batFile)) {
			writer.println(headsb);
			/*
			 * writer.println("setlocal enabledelayedexpansion"); writer.println(); writer.
			 * println("REM Get all the file paths passed as parameters and store them in a variable"
			 * ); writer.println("set files="); writer.println(); writer.println(
			 * "REM Loop through each file path and append it to the variable with a space and a --f option"
			 * ); writer.println("for %%f in (%*) do (");
			 * writer.println("\tset files=!files! --f \"%%f\""); writer.println(")");
			 * writer.println();
			 * writer.println("REM Print the content of the start command");
			 * writer.println("echo !files!"); writer.println(); writer.println("start /B "
			 * + diaNNFile.getName() + " " + cmd + "!files!"); writer.println();
			 * writer.println("endlocal"); writer.println();
			 */
			writer.println("start /B " + diaNNFile.getName() + " --cfg " + cfgFile.getAbsolutePath());
			writer.println();
			writer.print("exit");
			writer.close();

			PrintWriter cfgWriter = new PrintWriter(cfgFile);
			cfgWriter.print(cmd);
			for (int i = 0; i < raws.length; i++) {
				cfgWriter.print(" --f \"" + raws[i] + "\"");
			}
			cfgWriter.close();
		} catch (IOException e) {
			LOGGER.error(taskName + ": error in writing the .bat file to " + batFile, e);
			System.out.println(format.format(new Date()) + "\t" + taskName + "batFile" + batFile);
		}
		/*
		 * StringBuilder sb = new StringBuilder(); sb.append(batFile.getAbsolutePath());
		 * for (String raw : raws) { sb.append("\t").append(raw); }
		 * this.parList.add(sb.toString());
		 */
		this.parList.add(batFile.getAbsolutePath());
	}

	/**
	 * search library
	 * 
	 * @param par
	 * @param raws
	 * @param out
	 * @param lib
	 */
	public void addTask(DiannParameter par, String[] raws, String out, String lib, boolean genSpecLib,
			boolean autoInfer) {
		String cmd = par.generateLibSearch(raws, out, lib, genSpecLib, autoInfer);
		File batFile = new File(out + ".bat");
		File cfgFile = new File(out + ".cfg");
		try (PrintWriter writer = new PrintWriter(batFile)) {
			writer.println(headsb);
			/*
			 * writer.println("setlocal enabledelayedexpansion"); writer.println(); writer.
			 * println("REM Get all the file paths passed as parameters and store them in a variable"
			 * ); writer.println("set files="); writer.println(); writer.println(
			 * "REM Loop through each file path and append it to the variable with a space and a --f option"
			 * ); writer.println("for %%f in (%*) do (");
			 * writer.println("\tset files=!files! --f \"%%f\""); writer.println(")");
			 * writer.println();
			 * writer.println("REM Print the content of the start command");
			 * writer.println("echo !files!"); writer.println(); writer.println("start /B "
			 * + diaNNFile.getName() + " " + cmd + "!files!"); writer.println();
			 * writer.println("endlocal"); writer.println();
			 */
			writer.println("start /B " + diaNNFile.getName() + " --cfg " + cfgFile.getAbsolutePath());
			writer.println();
			writer.print("exit");
			writer.close();

			PrintWriter cfgWriter = new PrintWriter(cfgFile);
			cfgWriter.print(cmd);
			for (int i = 0; i < raws.length; i++) {
				cfgWriter.print(" --f \"" + raws[i] + "\"");
			}
			cfgWriter.close();
		} catch (IOException e) {
			LOGGER.error(taskName + ": error in writing the .bat file to " + batFile, e);
			System.out.println(
					format.format(new Date()) + "\t" + taskName + ": error in writing the .bat file to " + batFile);
		}
		/*
		 * StringBuilder sb = new StringBuilder(); sb.append(batFile.getAbsolutePath());
		 * for (String raw : raws) { sb.append("\t").append(raw); }
		 * this.parList.add(sb.toString());
		 */
		this.parList.add(batFile.getAbsolutePath());
	}

	/**
	 * search library
	 * 
	 * @param par
	 * @param libSearchPar
	 * @param out
	 * @param lib
	 */
	public void addTask(DiannParameter par, DiaLibSearchPar libSearchPar, String out, String... lib) {
		String cmd = par.generateLibSearch(libSearchPar, out, lib);
		File batFile = new File(out + ".bat");
		File cfgFile = new File(out + ".cfg");
		try (PrintWriter writer = new PrintWriter(batFile)) {
			writer.println(headsb);
			/*
			 * writer.println("setlocal enabledelayedexpansion"); writer.println(); writer.
			 * println("REM Get all the file paths passed as parameters and store them in a variable"
			 * ); writer.println("set files="); writer.println(); writer.println(
			 * "REM Loop through each file path and append it to the variable with a space and a --f option"
			 * ); writer.println("for %%f in (%*) do (");
			 * writer.println("\tset files=!files! --f \"%%f\""); writer.println(")");
			 * writer.println();
			 * writer.println("REM Print the content of the start command");
			 * writer.println("echo !files!"); writer.println(); writer.println("start /B "
			 * + diaNNFile.getName() + " " + cmd + "!files!"); writer.println();
			 * writer.println("endlocal"); writer.println();
			 */
			writer.println("start /B " + diaNNFile.getName() + " --cfg " + cfgFile.getAbsolutePath());
			writer.println();
			writer.print("exit");
			writer.close();

			String[] raws = libSearchPar.getDiaFiles();
			PrintWriter cfgWriter = new PrintWriter(cfgFile);
			cfgWriter.print(cmd);
			for (int i = 0; i < raws.length; i++) {
				cfgWriter.print(" --f \"" + raws[i] + "\"");
			}
			cfgWriter.close();
		} catch (IOException e) {
			LOGGER.error(taskName + ": error in writing the .bat file to " + batFile, e);
			System.out.println(
					format.format(new Date()) + "\t" + taskName + ": error in writing the .bat file to " + batFile);
		}
		/*
		 * String[] raws = libSearchPar.getDiaFiles(); StringBuilder sb = new
		 * StringBuilder(); sb.append(batFile.getAbsolutePath()); for (String raw :
		 * raws) { sb.append("\t").append(raw); } this.parList.add(sb.toString());
		 */
		this.parList.add(batFile.getAbsolutePath());
	}

	/**
	 * search library
	 * 
	 * @param par
	 * @param libSearchPar
	 * @param out
	 * @param lib
	 */
	public void addFastTask(DiannParameter par, DiaLibSearchPar libSearchPar, String out, String... lib) {
		String cmd = par.generateLibSearchFast(libSearchPar, out, lib);
		File batFile = new File(out + ".bat");
		File cfgFile = new File(out + ".cfg");
		try (PrintWriter writer = new PrintWriter(batFile)) {
			writer.println(headsb);
			/*
			 * writer.println("setlocal enabledelayedexpansion"); writer.println(); writer.
			 * println("REM Get all the file paths passed as parameters and store them in a variable"
			 * ); writer.println("set files="); writer.println(); writer.println(
			 * "REM Loop through each file path and append it to the variable with a space and a --f option"
			 * ); writer.println("for %%f in (%*) do (");
			 * writer.println("\tset files=!files! --f \"%%f\""); writer.println(")");
			 * writer.println();
			 * writer.println("REM Print the content of the start command");
			 * writer.println("echo !files!"); writer.println(); writer.println("start /B "
			 * + diaNNFile.getName() + " " + cmd + "!files!"); writer.println();
			 * writer.println("endlocal"); writer.println();
			 */
			writer.println("start /B " + diaNNFile.getName() + " --cfg " + cfgFile.getAbsolutePath());
			writer.println();
			writer.print("exit");
			writer.close();

			String[] raws = libSearchPar.getDiaFiles();
			PrintWriter cfgWriter = new PrintWriter(cfgFile);
			cfgWriter.print(cmd);
			for (int i = 0; i < raws.length; i++) {
				cfgWriter.print(" --f \"" + raws[i] + "\"");
			}
			cfgWriter.close();
		} catch (IOException e) {
			LOGGER.error(taskName + ": error in writing the .bat file to " + batFile, e);
			System.out.println(
					format.format(new Date()) + "\t" + taskName + ": error in writing the .bat file to " + batFile);
		}
		/*
		 * String[] raws = libSearchPar.getDiaFiles(); StringBuilder sb = new
		 * StringBuilder(); sb.append(batFile.getAbsolutePath()); for (String raw :
		 * raws) { sb.append("\t").append(raw); } this.parList.add(sb.toString());
		 */
		this.parList.add(batFile.getAbsolutePath());
	}

	/**
	 * 
	 * @param par
	 * @param libSearchPar
	 * @param output
	 * @param fast
	 * @param genLib
	 * @param multiQvaule
	 * @param libs
	 */
	public void addTask(DiannParameter par, DiaLibSearchPar libSearchPar, String output, boolean fast, boolean genLib,
			boolean qvalueCutoff, String... libs) {
		String cmd = par.generateLibSearch(libSearchPar, output, fast, genLib, qvalueCutoff, libs);
		File batFile = new File(output + ".bat");
		File cfgFile = new File(output + ".cfg");
		try (PrintWriter writer = new PrintWriter(batFile)) {
			writer.println(headsb);
			/*
			 * writer.println("setlocal enabledelayedexpansion"); writer.println(); writer.
			 * println("REM Get all the file paths passed as parameters and store them in a variable"
			 * ); writer.println("set files="); writer.println(); writer.println(
			 * "REM Loop through each file path and append it to the variable with a space and a --f option"
			 * ); writer.println("for %%f in (%*) do (");
			 * writer.println("\tset files=!files! --f \"%%f\""); writer.println(")");
			 * writer.println();
			 * writer.println("REM Print the content of the start command");
			 * writer.println("echo !files!"); writer.println(); writer.println("start /B "
			 * + diaNNFile.getName() + " " + cmd + "!files!"); writer.println();
			 * writer.println("endlocal"); writer.println();
			 */
			writer.println("start /B " + diaNNFile.getName() + " --cfg " + cfgFile.getAbsolutePath());
			writer.println();
			writer.print("exit");
			writer.close();

			String[] raws = libSearchPar.getDiaFiles();
			PrintWriter cfgWriter = new PrintWriter(cfgFile);
			cfgWriter.print(cmd);
			for (int i = 0; i < raws.length; i++) {
				cfgWriter.print(" --f \"" + raws[i] + "\"");
			}
			cfgWriter.close();
		} catch (IOException e) {
			LOGGER.error(taskName + ": error in writing the .bat file to " + batFile, e);
			System.out.println(
					format.format(new Date()) + "\t" + taskName + ": error in writing the .bat file to " + batFile);
		}
		/*
		 * String[] raws = libSearchPar.getDiaFiles(); StringBuilder sb = new
		 * StringBuilder(); sb.append(batFile.getAbsolutePath()); for (String raw :
		 * raws) { sb.append("\t").append(raw); } this.parList.add(sb.toString());
		 */
		this.parList.add(batFile.getAbsolutePath());
	}

	/**
	 * search library
	 * 
	 * @param par
	 * @param libSearchPar
	 * @param out
	 * @param lib
	 */
	public void addFastLibTask(DiannParameter par, DiaLibSearchPar libSearchPar, String out, String... lib) {
		String cmd = par.generateLibSearchFastLib(libSearchPar, out, lib);
		File batFile = new File(out + ".bat");
		File cfgFile = new File(out + ".cfg");
		try (PrintWriter writer = new PrintWriter(batFile)) {
			writer.println(headsb);
			/*
			 * writer.println("setlocal enabledelayedexpansion"); writer.println(); writer.
			 * println("REM Get all the file paths passed as parameters and store them in a variable"
			 * ); writer.println("set files="); writer.println(); writer.println(
			 * "REM Loop through each file path and append it to the variable with a space and a --f option"
			 * ); writer.println("for %%f in (%*) do (");
			 * writer.println("\tset files=!files! --f \"%%f\""); writer.println(")");
			 * writer.println();
			 * writer.println("REM Print the content of the start command");
			 * writer.println("echo !files!"); writer.println(); writer.println("start /B "
			 * + diaNNFile.getName() + " " + cmd + "!files!"); writer.println();
			 * writer.println("endlocal"); writer.println();
			 */
			writer.println("start /B " + diaNNFile.getName() + " --cfg " + cfgFile.getAbsolutePath());
			writer.println();
			writer.print("exit");
			writer.close();

			String[] raws = libSearchPar.getDiaFiles();
			PrintWriter cfgWriter = new PrintWriter(cfgFile);
			cfgWriter.print(cmd);
			for (int i = 0; i < raws.length; i++) {
				cfgWriter.print(" --f \"" + raws[i] + "\"");
			}
			cfgWriter.close();
		} catch (IOException e) {
			LOGGER.error(taskName + ": error in writing the .bat file to " + batFile, e);
			System.out.println(
					format.format(new Date()) + "\t" + taskName + ": error in writing the .bat file to " + batFile);
		}
		/*
		 * String[] raws = libSearchPar.getDiaFiles(); StringBuilder sb = new
		 * StringBuilder(); sb.append(batFile.getAbsolutePath()); for (String raw :
		 * raws) { sb.append("\t").append(raw); } this.parList.add(sb.toString());
		 */
		this.parList.add(batFile.getAbsolutePath());
	}

	/**
	 * search fasta
	 * 
	 * @param par
	 * @param raws
	 * @param out
	 * @param fasta
	 * @param libOut
	 * @param fast
	 */
	public void addTask(DiannParameter par, String[] raws, String out, String fasta, String libOut, boolean fast) {
		String cmd = par.generateFastaSearch(raws, out, fasta, libOut, fast);
		File batFile = new File(out + ".bat");
		File cfgFile = new File(out + ".cfg");
		try (PrintWriter writer = new PrintWriter(batFile)) {
			writer.println(headsb);
			/*
			 * writer.println("setlocal enabledelayedexpansion"); writer.println(); writer.
			 * println("REM Get all the file paths passed as parameters and store them in a variable"
			 * ); writer.println("set files="); writer.println(); writer.println(
			 * "REM Loop through each file path and append it to the variable with a space and a --f option"
			 * ); writer.println("for %%f in (%*) do (");
			 * writer.println("\tset files=!files! --f \"%%f\""); writer.println(")");
			 * writer.println();
			 * writer.println("REM Print the content of the start command");
			 * writer.println("echo !files!"); writer.println(); writer.println("start /B "
			 * + diaNNFile.getName() + " " + cmd + "!files!"); writer.println();
			 * writer.println("endlocal"); writer.println();
			 */
			writer.println("start /B " + diaNNFile.getName() + " --cfg " + cfgFile.getAbsolutePath());
			writer.println();
			writer.print("exit");
			writer.close();

			PrintWriter cfgWriter = new PrintWriter(cfgFile);
			cfgWriter.print(cmd);
			for (int i = 0; i < raws.length; i++) {
				cfgWriter.print(" --f \"" + raws[i] + "\"");
			}
			cfgWriter.close();
		} catch (IOException e) {
			LOGGER.error(taskName + ": error in writing the .bat file to " + batFile, e);
			System.out.println(
					format.format(new Date()) + "\t" + taskName + ": error in writing the .bat file to " + batFile);
		}
		/*
		 * StringBuilder sb = new StringBuilder(); sb.append(batFile.getAbsolutePath());
		 * for (String raw : raws) { sb.append("\t").append(raw); }
		 * this.parList.add(sb.toString());
		 */
		this.parList.add(batFile.getAbsolutePath());
	}

	/**
	 * generate library from fasta
	 * 
	 * @param par
	 * @param raws
	 * @param out
	 * @param fasta
	 * @param libOut
	 * @param fast
	 */
	public void addTask(DiannParameter par, String fasta, String libOut, boolean fast) {
		String cmd = par.generateFasta(fasta, libOut, fast);
		File batFile = new File(libOut + ".bat");
		try (PrintWriter writer = new PrintWriter(batFile)) {
			writer.print(headsb);
			writer.println("start /B " + diaNNFile.getName() + " " + cmd);
			writer.print(tailsb);
			writer.close();
		} catch (IOException e) {
			LOGGER.error(taskName + ": error in writing the .bat file to " + batFile, e);
			System.out.println(format.format(new Date()) + "\t" + taskName + "batFile" + batFile);
		}

		this.parList.add(batFile.getAbsolutePath());
	}

	/**
	 * library from .tsv to .speclib
	 * 
	 * @param par
	 * @param tsv
	 */
	public void addTask(DiannParameter par, String tsv) {
		String cmd = par.generateTsv2Lib(tsv);
		File batFile = new File(tsv + ".bat");
		try (PrintWriter writer = new PrintWriter(batFile)) {
			writer.print(headsb);
			writer.println("start /B " + diaNNFile.getName() + " " + cmd);
			writer.print(tailsb);
			writer.close();
		} catch (IOException e) {
			LOGGER.error(taskName + ": error in writing the .bat file to " + batFile, e);
			System.out.println(format.format(new Date()) + "\t" + taskName + "batFile" + batFile);
		}

		this.parList.add(batFile.getAbsolutePath());
	}

	/**
	 * library from .tsv to .speclib, extracting protein names from the fasta files
	 * 
	 * @param par
	 * @param tsv
	 */
	public void addTask(DiannParameter par, String tsv, String[] fasta) {
		String cmd = par.generateTsv2Lib(tsv, fasta);
		File batFile = new File(tsv + ".bat");
		try (PrintWriter writer = new PrintWriter(batFile)) {
			writer.print(headsb);
			writer.println("start /B " + diaNNFile.getName() + " " + cmd);
			writer.print(tailsb);
			writer.close();
		} catch (IOException e) {
			LOGGER.error(taskName + ": error in writing the .bat file to " + batFile, e);
			System.out.println(format.format(new Date()) + "\t" + taskName + "batFile" + batFile);
		}

		this.parList.add(batFile.getAbsolutePath());
	}

	/**
	 * library from .speclib to .tsv
	 * 
	 * @param par
	 * @param tsv
	 */
	public void addTask(DiannParameter par, String lib, String tsv) {
		String cmd = par.generateLib2Tsv(lib, tsv);
		File batFile = new File(tsv + ".bat");
		try (PrintWriter writer = new PrintWriter(batFile)) {
			writer.print(headsb);
			writer.println("start /B " + diaNNFile.getName() + " " + cmd);
			writer.print(tailsb);
			writer.close();
		} catch (IOException e) {
			LOGGER.error(taskName + ": error in writing the .bat file to " + batFile, e);
			System.out.println(
					format.format(new Date()) + "\t" + taskName + ": error in writing the .bat file to " + batFile);
		}

		this.parList.add(batFile.getAbsolutePath());
	}

	public void run(int waitHours) {
		run(1, waitHours);
	}

	public void run(int threadPool, int waitHours) {
		// Reset task state
		isRunning = true;
		isCancelled = false;
		isDone = false;

		if (this.parList.size() == 0) {
			isDone = true;
			isRunning = false;
			return;
		}

		executor = Executors.newFixedThreadPool(threadPool);
		for (int i = 0; i < parList.size(); i++) {

			String pari = parList.get(i);

			LOGGER.info(pari);

			executor.submit(new Runnable() {

				@Override
				public void run() {
					// Check for cancellation before starting
					if (isCancelled) {
						LOGGER.info("DiaNN task cancelled before execution: {}", pari);
						return;
					}

					try {
						LOGGER.info("Executing DiaNN conversion: {}", pari);

						ArrayList<String> commands = new ArrayList<String>();
						commands.add("cmd.exe");
						commands.add("/c");
						commands.add("start");
						commands.add("/min");

						String[] cs = pari.split("\t");
						for (String parString : cs) {
							commands.add(parString);
						}
						
						LOGGER.info("Process command: {}", commands.toString());
						
						ProcessBuilder pb = new ProcessBuilder(commands);
						Process p = pb.start();
						
						// Capture both stdout and stderr
						BufferedInputStream in = new BufferedInputStream(p.getInputStream());
						BufferedInputStream err = new BufferedInputStream(p.getErrorStream());
						BufferedReader inBr = new BufferedReader(new InputStreamReader(in));
						BufferedReader errBr = new BufferedReader(new InputStreamReader(err));
						
						String line0 = "";
						String line1 = "";
						while ((line1 = inBr.readLine()) != null) {
							// Check for cancellation during execution
							if (isCancelled) {
								LOGGER.info("DiaNN task cancelled during execution: {}", pari);
								p.destroyForcibly();
								inBr.close();
								errBr.close();
								in.close();
								err.close();
								return;
							}
							
							if (!line1.equals(line0)) {
								LOGGER.info("DiaNN stdout: {}", line1);
								System.out.println(format.format(new Date()) + "\t" + line1);
								line0 = line1;
							}
						}
						
						// Read error stream
						String errorLine;
						while ((errorLine = errBr.readLine()) != null) {
							LOGGER.error("DiaNN stderr: {}", errorLine);
							System.err.println(format.format(new Date()) + "\tDiaNN ERROR: " + errorLine);
						}

						int exitCode = p.waitFor();
						LOGGER.info("DiaNN process completed with exit code: {}", exitCode);
						
						if (exitCode != 0) {
							LOGGER.error("DiaNN conversion failed with exit code: {}", exitCode);
							System.err.println(format.format(new Date()) + "\tDiaNN conversion failed with exit code: " + exitCode);
							System.err.println("ERROR: cannot read .raw files, please download and install .NET SDK 8.0.407 or later");
							System.err.println("https://dotnet.microsoft.com/en-us/download/dotnet/8.0");
						} else {
							LOGGER.info("DiaNN conversion completed successfully");
						}
						
						inBr.close();
						errBr.close();
						in.close();
						err.close();

					} catch (IOException | InterruptedException e) {
						LOGGER.error("Error executing DiaNN conversion: {}", e.getMessage(), e);
						System.err.println(format.format(new Date()) + "\tError executing DiaNN conversion: " + e.getMessage());
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
						+ ": task dosen't finish in a long time, please restart MetaLab after the DiaNN task finish");
				System.out.println(format.format(new Date()) + "\t" + taskName
						+ ": task dosen't finish in a long time, please restart MetaLab after the DiaNN task finish");
			}

		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			LOGGER.error(taskName + ": failed", e);
			System.out.println(format.format(new Date()) + "\t" + taskName + ": failed");
		} finally {
			isDone = true;
			isRunning = false;
		}
	}

	public boolean runCMD(String cmd, JProgressBar bar) {

		try {

			int id = diaNNFile.getAbsolutePath().indexOf(":");
			if (id < 0) {
				if (id < 0) {
					LOGGER.error(taskName + ": cann't find the directory of DiaNN.exe");
					System.out.println(
							format.format(new Date()) + "\t" + taskName + ": cann't find the directory of DiaNN.exe");

					return false;
				}
			}

			PrintWriter writer = new PrintWriter(cmd);
			writer.println(diaNNFile.getAbsolutePath().substring(0, id) + ":");
			writer.println("cd " + diaNNFile.getParent());

			for (int i = 0; i < parList.size(); i++) {
				writer.println(diaNNFile.getName() + " " + parList.get(i));
			}

			writer.print("exit");
			writer.close();

		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			LOGGER.error(taskName + ": error in writing the executable command for DiaNN.exe to " + cmd, e);
			System.out.println(format.format(new Date()) + "\t" + taskName
					+ ": error in writing the executable command for DiaNN.exe to " + cmd);

			return false;
		}

		try {
			int finishCount = 0;

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

				if (lineStr.startsWith("DIA-NN 1.8 (Data-Independent Acquisition by Neural Networks)")) {

					double percent = (double) finishCount / (double) parList.size() * 100.0;
					bar.setValue((int) percent);

					finishCount++;
				}
			}

			if (p.waitFor() != 0) {
				if (p.exitValue() == 1)
					LOGGER.error("false");

				bar.setValue(100);
			}
			inBr.close();
			in.close();
		} catch (Exception e) {
			LOGGER.error(taskName + ": error in running the executable command by DiaNN.exe", e);
			System.out.println(format.format(new Date()) + "\t" + taskName
					+ ": error in running the executable command by DiaNN.exe");

			return false;
		}

		return true;
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
	 * This will shutdown the executor service and terminate all running processes
	 */
	public void cancel() {
		if (isRunning && !isCancelled) {
			LOGGER.info(taskName + ": Cancellation requested");
			System.out.println(format.format(new Date()) + "\t" + taskName + ": Cancellation requested");
			
			isCancelled = true;
			
			// Shutdown the executor service if it's running
			if (executor != null && !executor.isShutdown()) {
				LOGGER.info(taskName + ": Shutting down executor service");
				System.out.println(format.format(new Date()) + "\t" + taskName + ": Shutting down executor service");
				
				executor.shutdownNow();
				
				try {
					// Wait a bit for the executor to terminate
					if (!executor.awaitTermination(5, TimeUnit.SECONDS)) {
						LOGGER.warn(taskName + ": Executor service did not terminate gracefully within 5 seconds");
						System.out.println(format.format(new Date()) + "\t" + taskName + ": Executor service did not terminate gracefully within 5 seconds");
					}
				} catch (InterruptedException e) {
					LOGGER.warn(taskName + ": Interrupted while waiting for executor termination", e);
					Thread.currentThread().interrupt();
				}
			}
			
			isDone = true;
			isRunning = false;
		}
	}
}
