package uk.ac.quadram.metapilot.task.dia;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.Date;

import org.apache.logging.log4j.Logger;

import javafx.scene.control.ProgressBar;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters.ToolType;
import uk.ac.quadram.metapilot.task.MetaAbstractTask;
import uk.ac.quadram.metapilot.task.dia.par.MetaParameterDiaTg;

public class MetaDiaTargetDBTask extends MetaAbstractTask {

	private String magDbDir;
	private String taxaFile;
	private String eggnogDir;

	public static final String[] FUNCTIONS = new String[] { "GOs", "EC", "KEGG_ko", "KEGG_Pathway", "KEGG_Module",
			"KEGG_Reaction", "KEGG_rclass", "BRITE", "KEGG_TC", "CAZy" };

	protected static final String taskName = "Create taxa-function database";

	public MetaDiaTargetDBTask(MetaParameterDiaTg metaPar, MetaResourceParameters mrp, ProgressBar progressBar,
			Logger logger) {
		super(metaPar, mrp, progressBar, logger);

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
	protected Boolean call() throws Exception {
		updateMessage(taskName + ": preparing Python database generation script");
		updateProgress(0, 100);

		// Create the Python script for database generation
		File pythonFile = new File(magDbDir, "profuncSql.py");

		try (InputStream inputStream = MetaDiaTargetDBTask.class.getResourceAsStream("/python_scripts/db/profuncSql.py");
				PrintWriter writer = new PrintWriter(pythonFile)) {
			if (inputStream == null) {
				logger.error(taskName + ": profuncSql.py not found in classpath at /python_scripts/db/profuncSql.py");
				System.err.println(
						format.format(new Date()) + "\t" + taskName + ": profuncSql.py not found in classpath");
				return false;
			}

			BufferedReader reader = new BufferedReader(new InputStreamReader(inputStream));
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
			File pythonExeFile = new File(mrp.getToolPath(ToolType.PYTHON));
			if (!pythonExeFile.exists()) {
				logger.error(taskName + ": Python executable was not found: " + pythonExeFile.getAbsolutePath());
				return false;
			}

			updateMessage(taskName + ": creating taxa-function SQL database");
			updateProgress(10, 100);

			ProcessBuilder pb = new ProcessBuilder(pythonExeFile.getAbsolutePath(), pythonFile.getAbsolutePath());
			pb.directory(new File(magDbDir));
			pb.redirectErrorStream(true);
			Process p = pb.start();
			try (BufferedReader inBr = new BufferedReader(new InputStreamReader(p.getInputStream()))) {
				String lineStr;
				while ((lineStr = inBr.readLine()) != null) {
					System.out.println(lineStr);
					logger.info(taskName + ": " + lineStr);
				}
			}

			int exitCode = p.waitFor();
			if (exitCode != 0) {
				logger.error(taskName + ": Python database generation failed, exit code " + exitCode);
				return false;
			}

			updateProgress(100, 100);
			updateMessage(taskName + ": database created successfully");
		} catch (Exception e) {
			logger.error(taskName + ": error in creating the taxa-function SQL database", e);
			System.err.println(
					format.format(new Date()) + "\t" + taskName + ": error in creating the taxa-function SQL database");
			return false;
		}

		return true;
	}
}
