package uk.ac.quadram.metapilot.core.deepdetect;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.util.function.Consumer;

import uk.ac.quadram.metapilot.dbSearch.deepDetect.DeepDetectParameter;

public final class DeepDetectWindowsRunner {

	private DeepDetectWindowsRunner() {
	}

	public static File createBatchFile(File deepDetectExe, DeepDetectParameter parameter, File batchFile)
			throws IOException {
		if (deepDetectExe == null || !deepDetectExe.exists()) {
			throw new IOException("DeepDetect executable not found: "
					+ (deepDetectExe == null ? "null" : deepDetectExe.getAbsolutePath()));
		}
		if (parameter == null || parameter.getInput() == null || parameter.getOutput() == null) {
			throw new IOException("DeepDetectParameter must have input and output set.");
		}
		if (batchFile == null) {
			File outputDir = parameter.getOutput().getParentFile();
			if (outputDir == null) {
				outputDir = parameter.getInput().getParentFile();
			}
			if (outputDir == null) {
				outputDir = new File(System.getProperty("java.io.tmpdir"));
			}
			String inputName = parameter.getInput().getName();
			int dot = inputName.lastIndexOf('.');
			String baseName = dot > 0 ? inputName.substring(0, dot) : inputName;
			batchFile = new File(outputDir, baseName + ".deepDetect.bat");
		}
		File parent = batchFile.getParentFile();
		if (parent != null) {
			parent.mkdirs();
		}
		try (PrintWriter writer = new PrintWriter(batchFile, StandardCharsets.UTF_8)) {
			writer.println("@echo off");
			writer.println("cd /d \"" + deepDetectExe.getParent() + "\"");
			writer.println("\"" + deepDetectExe.getAbsolutePath() + "\" " + parameter.toQuotedString());
		}
		return batchFile;
	}

	public static int runBatchFile(File batchFile, Consumer<String> outputHandler) throws IOException, InterruptedException {
		if (batchFile == null || !batchFile.exists()) {
			throw new IOException("DeepDetect batch file not found: "
					+ (batchFile == null ? "null" : batchFile.getAbsolutePath()));
		}
		Consumer<String> handler = outputHandler == null ? ignored -> {
		} : outputHandler;
		ProcessBuilder processBuilder = new ProcessBuilder("cmd.exe", "/c", batchFile.getAbsolutePath());
		processBuilder.redirectErrorStream(true);
		Process process = processBuilder.start();
		try (BufferedReader reader = new BufferedReader(
				new InputStreamReader(process.getInputStream(), StandardCharsets.UTF_8))) {
			String line;
			while ((line = reader.readLine()) != null) {
				handler.accept(line);
			}
		}
		return process.waitFor();
	}

	public static int run(File deepDetectExe, DeepDetectParameter parameter, File batchFile,
			Consumer<String> outputHandler) throws IOException, InterruptedException {
		File createdBatchFile = createBatchFile(deepDetectExe, parameter, batchFile);
		return runBatchFile(createdBatchFile, outputHandler);
	}

	public static void main(String[] args) throws IOException {
		DeepDetectParameter par = new DeepDetectParameter();
		par.setInput(new File("Z:\\PRIDE_dataset\\MetaPep\\MetaPep_peptide_sequence_database_for_deepdetect.fasta"));;
		par.setOutput(new File("Z:\\PRIDE_dataset\\MetaPep\\MetaPep_deepdetect.tsv"));
		DeepDetectWindowsRunner.createBatchFile(new File("C:\\Users\\daz25hak\\MetaProAI\\Resources\\DeepDetect"
				+ "\\DeepDetect.exe"), par, new File("Z:\\PRIDE_dataset\\MetaPep\\metapep.bat"));
	}
}
