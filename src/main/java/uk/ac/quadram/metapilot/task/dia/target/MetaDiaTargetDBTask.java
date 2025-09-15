package uk.ac.quadram.metapilot.task.dia.target;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.Date;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import javafx.scene.control.ProgressBar;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters.ToolType;
import uk.ac.quadram.metapilot.task.MetaAbstractTask;
import uk.ac.quadram.metapilot.task.dia.target.par.MetaParameterDiaTg;

public class MetaDiaTargetDBTask extends MetaAbstractTask {

	private String magDbDir;
	private String taxaFile;
	private String eggnogDir;

	public static final String[] FUNCTIONS = new String[] { "GOs", "EC", "KEGG_ko", "KEGG_Pathway", "KEGG_Module",
			"KEGG_Reaction", "KEGG_rclass", "BRITE", "KEGG_TC", "CAZy" };

	protected static final Logger LOGGER = LogManager.getLogger(MetaDiaTargetDBTask.class);
	protected static final String taskName = "Create taxa-function database";

	public MetaDiaTargetDBTask(MetaParameterDiaTg metaPar, MetaResourceParameters mrp, ProgressBar progressBar) {
		super(metaPar, mrp, progressBar);
		
		// Get MAG database directory from the parameter
		magDbDir = metaPar.getMagDbItem().getCurrentFile().getAbsolutePath();
		this.eggnogDir = magDbDir + "\\eggNOG";
		this.taxaFile = magDbDir + "\\genomes-all_metadata.tsv";
	}

	@Override
	protected String getTaskName() {
		return taskName;
	}

	@Override
	protected Logger getLogger() {
		return LOGGER;
	}

	@Override
	protected Boolean call() throws Exception {
		// Create the Python script for database generation
		File pythonFile = new File(magDbDir, "profuncSql.py");

		try (InputStream inputStream = MetaDiaTargetDBTask.class.getResourceAsStream("profuncSql.py");
				BufferedReader reader = new BufferedReader(new InputStreamReader(inputStream));
				PrintWriter writer = new PrintWriter(pythonFile)) {
			String line;
			while ((line = reader.readLine()) != null) {
				writer.println(line);
			}
			inputStream.close();
			reader.close();

			writer.println("if __name__ == \"__main__\":");
			writer.println("    folder_path = \"" + eggnogDir.replaceAll("\\\\", "\\\\\\\\") + "\"");
			writer.println("    db_path = \"" + magDbDir.replaceAll("\\\\", "\\\\\\\\") + "\\\\pro_func.db\"");
			writer.println("    func_names = [");
			StringBuilder sb = new StringBuilder();
			for (int i = 0; i < FUNCTIONS.length; i++) {
				sb.append("\"").append(FUNCTIONS[i]).append("\",");
			}
			sb.deleteCharAt(sb.length() - 1);
			writer.println("        " + sb);
			writer.println("    ]");

			writer.println("    for func_name in func_names:");
			writer.println("        process_tsv_files(folder_path, db_path, func_name)");
			writer.println();

			writer.println("    taxa_file = \"" + taxaFile.replaceAll("\\\\", "\\\\\\\\") + "\"");
			writer.println("    process_tsv_to_sqlite(taxa_file, db_path)");

			writer.close();
		} catch (IOException e) {
			System.err.println("An error occurred while reading the file: " + e.getMessage());
			return false;
		}

		try {
			File batFile = new File(pythonFile.getParent(), pythonFile.getName() + ".bat");
			PrintWriter batWriter = new PrintWriter(batFile);

			File pythonExeFile = new File(mrp.getToolPath(ToolType.PYTHON));

			StringBuilder sb = new StringBuilder();
			sb.append(pythonExeFile.getName()).append(" ");
			sb.append(pythonFile.getAbsolutePath());

			batWriter.println("@echo off");
			batWriter.println("cd /d \"" + pythonExeFile.getParent() + "\"");
			batWriter.println("start /B " + sb);

			batWriter.println("if errorlevel 1 (");
			batWriter.println("\techo Error: Failed to start the application.");
			batWriter.println(") else (");
			batWriter.println("\techo Application started successfully.");
			batWriter.println(")");
			batWriter.print("exit");
			batWriter.close();

			String[] args = { "cmd.exe", "/c", "start", batFile.getAbsolutePath() };
			ProcessBuilder pb = new ProcessBuilder(args);
			Process p = pb.start();
			BufferedInputStream in = new BufferedInputStream(p.getInputStream());
			BufferedReader inBr = new BufferedReader(new InputStreamReader(in));
			String lineStr = null;
			while ((lineStr = inBr.readLine()) != null) {
				System.out.println(lineStr);
			}

			if (p.waitFor() != 0) {
				if (p.exitValue() == 1)
					System.err.println("false");
			}
			inBr.close();
			in.close();
		} catch (Exception e) {
			LOGGER.error(taskName + ": error in creating the taxa-function SQL database", e);
			System.err.println(
					format.format(new Date()) + "\t" + taskName + ": error in creating the taxa-function SQL database");
			return false;
		}

		return true;
	}
}