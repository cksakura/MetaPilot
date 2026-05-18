package uk.ac.quadram.metapilot.dbSearch.deepDetect;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileFilter;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashSet;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;

import org.apache.commons.io.FileUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import javafx.application.Platform;
import javafx.scene.control.ProgressBar;
import javafx.scene.control.TextField;
import uk.ac.quadram.metapilot.core.deepdetect.DeepDetectPythonRunner;
import uk.ac.quadram.metapilot.core.deepdetect.DeepDetectWindowsRunner;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters.ToolType;
import uk.ac.quadram.metapilot.task.MetaAbstractTask;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;
import uk.ac.quadram.metapilot.task.par.MetaBaseParameter;

public class DeepDetectTask extends MetaAbstractTask {

	private static String taskName = "DeepDetect task";
	private static final Logger LOGGER = LogManager.getLogger(DeepDetectTask.class);
	private SimpleDateFormat format = new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss");

	private File deepDetectFile;
	private ArrayList<File> batFileList;
	private HashSet<String> genomeSet;

	private DeepDetectParameter par;
	private ExecutorService executor;
	private int threads;
	private TextField javafxTimeLeft;

	private MagDbItem magdb;

	public DeepDetectTask(MetaBaseParameter metaPar, MetaResourceParameters mrp, ProgressBar progressBar, Logger logger) {
		super(metaPar, mrp, progressBar, logger);
		this.deepDetectFile = new File(mrp.getToolPath(ToolType.DEEPDETECT));
		if (!deepDetectFile.exists()) {
			LOGGER.error(taskName + ": can't find DeepDetect.exe at " + deepDetectFile.getAbsolutePath());
			System.out.println(format.format(new Date()) + "\t" + taskName + ": can't find DeepDetect.exe");
		}
		this.batFileList = new ArrayList<File>();
	}

	public DeepDetectTask(MetaBaseParameter metaPar, MetaResourceParameters mrp, ProgressBar progressBar,
			HashSet<String> genomeSet, Logger logger) {
		this(metaPar, mrp, progressBar, logger);
		this.genomeSet = genomeSet;
	}

	// JavaFX constructor (par, magdb, threads, time display)
	public DeepDetectTask(MetaBaseParameter metaPar, MetaResourceParameters mrp, ProgressBar progressBar,
			DeepDetectParameter par, MagDbItem magdb, int threads, TextField javafxTimeLeft, Logger logger) {
		this(metaPar, mrp, progressBar, logger);
		this.par = par;
		this.magdb = magdb;
		this.threads = threads;
		this.executor = Executors.newFixedThreadPool(threads);
		this.javafxTimeLeft = javafxTimeLeft;
	}

	/**
	 * Predict all the genomes in a MAG database
	 */
	protected Boolean call() throws Exception {
		// TODO Auto-generated method stub

		if (javafxTimeLeft != null) {
			Platform.runLater(() -> {
				javafxTimeLeft.setText("Calculating...");
				System.out.println(format.format(new Date()) + "\t" + "Peptide detectability predection started");
			});
		}

		File dbfile = magdb.getOriginalDbFolder();
		File[] files = dbfile.listFiles();

		File predictedFile = new File(dbfile.getParentFile(), "predicted");
		if (!predictedFile.exists()) {
			predictedFile.mkdir();
		}

		for (int i = 0; i < files.length; i++) {
			String name = files[i].getName();
			if (name.endsWith(".faa") || name.endsWith(".fasta")) {
				name = name.substring(0, name.lastIndexOf("."));
				if (genomeSet == null || (genomeSet != null && genomeSet.contains(name))) {
					File predictedPepFile = new File(predictedFile, name + ".predicted.tsv");
					if (!predictedPepFile.exists()) {
						File inputToUse = files[i];
						try {
							inputToUse = ensureFastaHeadersHaveProteinName(files[i], predictedFile);
						} catch (IOException e) {
							LOGGER.error(taskName + ": error preparing FASTA for " + name, e);
							return false;
						}
						par.setInput(inputToUse);
						par.setOutput(predictedPepFile);

						File batFile = new File(predictedFile, name + ".deepDetect.bat");
						try {
							DeepDetectWindowsRunner.createBatchFile(deepDetectFile, par, batFile);
						} catch (IOException e) {
							LOGGER.error(taskName + ": error in writing the .bat file", e);
							return false;
						}
						this.batFileList.add(batFile);
					}
				}
			}
		}

		int totalGenomeCount = magdb.getSpeciesCount();

		AtomicInteger completedTasks = new AtomicInteger(magdb.getSpeciesCount() - batFileList.size());

		LOGGER.info(taskName + ": " + batFileList.size() + " files will be processed");
		System.out.println(format.format(new Date()) + "\t" + batFileList.size() + " files will be processed");

		for (int i = 0; i < batFileList.size(); i++) {
			File batFile = batFileList.get(i);

			executor.submit(new Runnable() {

				@Override
				public void run() {
					// TODO Auto-generated method stub

					try {

						long start = System.currentTimeMillis();

						int exitCode = DeepDetectWindowsRunner.runBatchFile(batFile,
								line -> Platform.runLater(() -> System.out.println(format.format(new Date()) + "\t" + line)));
						if (exitCode != 0) {
							LOGGER.error("DeepDetect exited with code " + exitCode);
						}

						int completed = completedTasks.incrementAndGet();
						double progress = (double) completed / (double) totalGenomeCount;
						// update the value of the progress bar
						progressBar.setProgress(progress);
						updateMessage(completed + "/" + totalGenomeCount);

						long end = System.currentTimeMillis();
						int runOneTask = (int) ((end - start) / 60000);
						int restTime = (int) ((double) (totalGenomeCount - completed) * runOneTask / (double) threads);
						if (javafxTimeLeft != null) {
							Platform.runLater(() -> {
								if (restTime == 0) {
									javafxTimeLeft.setText("About <1 minutes left");
								} else {
									javafxTimeLeft.setText("About " + restTime + " minutes left");
								}
							});
						}
						String batName = batFile.getName();

						LOGGER.info(
								taskName + ": predicting " + batName.substring(0, batName.length() - 4) + " finished");
						System.out.println(format.format(new Date()) + "\t" + "predicting "
								+ batName.substring(0, batName.length() - 4) + " finished");

					} catch (IOException | InterruptedException e) {
						// TODO Auto-generated catch block
						LOGGER.error(taskName + ": error in predicting peptide detectability", e);
						System.err.println(format.format(new Date()) + "\t" + taskName
								+ ": error in predicting peptide detectability");
					}
				}
			});
		}

		try {

			executor.shutdown();

			boolean finish = executor.awaitTermination(totalGenomeCount, TimeUnit.HOURS);

			cleanupTempFiles(predictedFile);
			batFileList = new ArrayList<File>();

			for (int i = 0; i < files.length; i++) {
				String name = files[i].getName();
				if (name.endsWith(".faa") || name.endsWith(".fasta")) {
					name = name.substring(0, name.lastIndexOf("."));

					File predictedPepFile = new File(predictedFile, name + ".predicted.tsv");
					if (!predictedPepFile.exists()) {

						File tempFastaFile = new File(predictedFile, name + ".faa");

						try {

							BufferedReader reader = new BufferedReader(new FileReader(files[i]));
							PrintWriter tempFastaWriter = new PrintWriter(tempFastaFile);
							StringBuilder sb = new StringBuilder();
							String line = null;
							boolean use = true;
							while ((line = reader.readLine()) != null) {
								if (line.startsWith(">")) {
									if (use && sb.length() > 0) {
										tempFastaWriter.print(sb);
									}

									use = true;
									sb = new StringBuilder();
									sb.append(line).append("\n");
								} else {
									if (line.contains("X")) {
										use = false;
									} else {
										sb.append(line).append("\n");
									}
								}
							}

							if (use && sb.length() > 0) {
								tempFastaWriter.print(sb);
							}

							reader.close();
							tempFastaWriter.close();

						} catch (IOException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
							LOGGER.error(taskName + ": error in writing the temp fasta file", e);
							return false;
						}

						File inputToUse = tempFastaFile;
						try {
							inputToUse = ensureFastaHeadersHaveProteinName(tempFastaFile, predictedFile);
						} catch (IOException e) {
							LOGGER.error(taskName + ": error preparing FASTA for " + name, e);
							return false;
						}
						par.setInput(inputToUse);
						par.setOutput(predictedPepFile);

						File batFile = new File(predictedFile, name + ".deepDetect.bat");
						try {
							DeepDetectWindowsRunner.createBatchFile(deepDetectFile, par, batFile);
						} catch (IOException e) {
							LOGGER.error(taskName + ": error in writing the .bat file", e);
							return false;
						}
						this.batFileList.add(batFile);
					}
				}
			}

			if (batFileList.size() > 0) {

				this.executor = Executors.newFixedThreadPool(threads);

				for (int i = 0; i < batFileList.size(); i++) {
					File batFile = batFileList.get(i);
					if (batFile.exists()) {
						executor.submit(new Runnable() {

							@Override
							public void run() {
								// TODO Auto-generated method stub

								try {

									long start = System.currentTimeMillis();

									int exitCode = DeepDetectWindowsRunner.runBatchFile(batFile,
											line -> Platform.runLater(
													() -> System.out.println(format.format(new Date()) + "\t" + line)));
									if (exitCode != 0) {
										LOGGER.error("DeepDetect exited with code " + exitCode);
									}

									int completed = completedTasks.incrementAndGet();
									double progress = (double) completed / (double) totalGenomeCount;
									// update the value of the progress bar
									progressBar.setProgress(progress);
									updateMessage(completed + "/" + totalGenomeCount);

									long end = System.currentTimeMillis();
									int runOneTask = (int) ((end - start) / 60000);
									int restTime = (int) ((double) (totalGenomeCount - completed) * runOneTask
											/ (double) threads);
									if (javafxTimeLeft != null) {
										Platform.runLater(() -> {
											if (restTime == 0) {
												javafxTimeLeft.setText("About <1 minutes left");
											} else {
												javafxTimeLeft.setText("About " + restTime + " minutes left");
											}
										});
									}

								} catch (IOException | InterruptedException e) {
									// TODO Auto-generated catch block
									LOGGER.error(taskName + ": error in predicting peptide detectability", e);
									System.err.println(format.format(new Date()) + "\t" + taskName
											+ ": error in predicting peptide detectability");
								}
							}
						});
					}
				}

				executor.shutdown();

				finish = executor.awaitTermination(totalGenomeCount, TimeUnit.HOURS);

				cleanupTempFiles(predictedFile);
			}

			progressBar.setProgress(1.0);
			updateMessage(totalGenomeCount + "/" + totalGenomeCount);

			if (javafxTimeLeft != null) {
				Platform.runLater(() -> javafxTimeLeft.setText("Task finished"));
			}

			this.batFileList = new ArrayList<File>();

			if (finish) {
				LOGGER.info(taskName + ": finished");
				System.out.println(format.format(new Date()) + "\t" + taskName + ": finished");
			} else {
				LOGGER.info(taskName
						+ ": task dosen't finish in a long time, please restart MetaLab after the DeepDetect task finish");
				System.out.println(format.format(new Date()) + "\t" + taskName
						+ ": task dosen't finish in a long time, please restart MetaLab after the DeepDetect task finish");

				return false;
			}

		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			LOGGER.error(taskName + ": failed", e);
			System.err.println(format.format(new Date()) + "\t" + taskName + ": failed");
		}

		int predictedFileCount = predictedFile.listFiles(new FileFilter() {

			@Override
			public boolean accept(File pathname) {
				// TODO Auto-generated method stub
				if (pathname.getName().endsWith(".predicted.tsv"))
					return true;
				return false;
			}
		}).length;

		if (predictedFileCount == magdb.getSpeciesCount()) {
			return true;
		} else {
			return false;
		}
	}

	public void forceStop() {
		if (executor != null) {
			this.executor.shutdown();
			Runtime run = Runtime.getRuntime();
			try {
				Process p = run.exec(new String[] { "taskkill", "/F", "/IM", "DeepDetect.exe" });
				BufferedInputStream in = new BufferedInputStream(p.getInputStream());
				BufferedReader inBr = new BufferedReader(new InputStreamReader(in));
				String lineStr;
				while ((lineStr = inBr.readLine()) != null) {
					if (lineStr.startsWith("SUCCESS")) {
						System.out.println(lineStr);
					}
				}
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
	}

	public void addTask(DeepDetectParameter par) {
		File inputFile = par.getInput();
		File outputFile = par.getOutput();
		String name = inputFile.getName();
		name = name.substring(0, name.lastIndexOf("."));

		File batFile = new File(outputFile.getParent(), name + ".deepDetect.bat");
		try {
			DeepDetectWindowsRunner.createBatchFile(deepDetectFile, par, batFile);
		} catch (IOException e) {
			LOGGER.error(taskName + ": error in writing the .bat file", e);
		}
		this.batFileList.add(batFile);
	}

	public void run(int waitHours) {
		run(1, waitHours);
	}

	public void run(int threadPool, int waitHours) {

		if (this.batFileList.size() == 0) {
			return;
		}

		ExecutorService executor = Executors.newFixedThreadPool(threadPool);
		for (int i = 0; i < batFileList.size(); i++) {
			File batFile = batFileList.get(i);

			executor.submit(new Runnable() {

				@Override
				public void run() {
					// TODO Auto-generated method stub

					try {
						String[] line0 = { "" };
						DeepDetectWindowsRunner.runBatchFile(batFile, line1 -> {
							if (!line1.equals(line0[0])) {
								System.out.println(format.format(new Date()) + "\t" + line1);
								line0[0] = line1;
							}
						});

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

			this.batFileList = new ArrayList<File>();

			if (finish) {
				LOGGER.info(taskName + ": finished");
				System.out.println(format.format(new Date()) + "\t" + taskName + ": finished");
			} else {
				LOGGER.info(taskName
						+ ": task dosen't finish in a long time, please restart MetaLab after the DeepDetect task finish");
				System.out.println(format.format(new Date()) + "\t" + taskName
						+ ": task dosen't finish in a long time, please restart MetaLab after the DeepDetect task finish");
			}

		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			LOGGER.error(taskName + ": failed", e);
			System.out.println(format.format(new Date()) + "\t" + taskName + ": failed");
		}
	}

	@Override
	protected String getTaskName() {
		// TODO Auto-generated method stub
		return taskName;
	}

	@Override
	protected Logger getLogger() {
		// TODO Auto-generated method stub
		return LOGGER;
	}

	/**
	 * Removes temporary .bat and _deepDetect_input.faa files from the predicted folder.
	 * Never touches original genome files or .predicted.tsv output files.
	 */
	private static void cleanupTempFiles(File predictedFolder) {
		if (predictedFolder == null || !predictedFolder.exists()) {
			return;
		}
		File[] tempFiles = predictedFolder.listFiles();
		if (tempFiles == null) {
			return;
		}
		for (File f : tempFiles) {
			String name = f.getName();
			if (name.endsWith(".bat") || name.endsWith("_deepDetect_input.faa")
					|| (name.endsWith(".faa") && !name.endsWith(".predicted.tsv"))) {
				FileUtils.deleteQuietly(f);
			}
		}
	}

	/**
	 * Prepares a FASTA file for DeepDetect by applying three rules:
	 * 1. Headers without a space get a duplicated accession as protein name
	 *    (DeepDetect requires ">accession proteinName" format)
	 * 2. Trailing '*' (stop codon) is stripped from protein sequences
	 * 3. Proteins containing non-standard amino acids (anything outside the 20
	 *    standard residues ACDEFGHIKLMNPQRSTVWY) are removed entirely
	 */
	public static File ensureFastaHeadersHaveProteinName(File inputFile, File outputDir) throws IOException {
		return DeepDetectPythonRunner.ensureFastaHeadersHaveProteinName(inputFile, outputDir);
	}

	/**
	 * Test method: runs DeepDetect once using only the exe file and
	 * DeepDetectParameter. Does not create a DeepDetectTask. Writes a .bat file
	 * and runs it via ExecutorService.
	 *
	 * @param deepDetectExe path to DeepDetect.exe (e.g. from MetaResourceParameters)
	 * @param par           parameter with input and output files set
	 * @param timeoutHours  max time to wait for the process
	 * @return true if the process exited with 0, false otherwise or on error
	 */
	public static boolean testDeepDetect(File deepDetectExe, DeepDetectParameter par, int timeoutHours) {
		if (deepDetectExe == null || !deepDetectExe.exists()) {
			LOGGER.error(taskName + ": testDeepDetect - exe file not found: " + (deepDetectExe != null ? deepDetectExe.getAbsolutePath() : "null"));
			return false;
		}
		if (par == null || par.getInput() == null || par.getOutput() == null) {
			LOGGER.error(taskName + ": testDeepDetect - DeepDetectParameter must have input and output set");
			return false;
		}
		File batDir = par.getOutput().getParentFile();
		if (batDir == null || !batDir.exists()) {
			batDir = par.getInput().getParentFile();
		}
		if (batDir == null || !batDir.exists()) {
			batDir = new File(System.getProperty("java.io.tmpdir"));
		}
		try {
			File modifiedInput = ensureFastaHeadersHaveProteinName(par.getInput(), batDir);
			if (modifiedInput != par.getInput()) {
				par.setInput(modifiedInput);
			}
		} catch (IOException e) {
			LOGGER.error(taskName + ": testDeepDetect - error preparing FASTA for DeepDetect", e);
			return false;
		}
		SimpleDateFormat fmt = new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss");
		File batFile = new File(batDir, "deepDetect_test.bat");
		try {
			DeepDetectWindowsRunner.createBatchFile(deepDetectExe, par, batFile);
		} catch (IOException e) {
			LOGGER.error(taskName + ": testDeepDetect - error writing .bat file", e);
			return false;
		}
		LOGGER.info(taskName + ": testDeepDetect - wrote " + batFile.getAbsolutePath());
		System.out.println(fmt.format(new Date()) + "\t" + taskName + ": testDeepDetect - running " + batFile.getAbsolutePath());

		ExecutorService executor = Executors.newFixedThreadPool(1);
		AtomicInteger exitCode = new AtomicInteger(-1);
		executor.submit(() -> {
			try {
				exitCode.set(DeepDetectWindowsRunner.runBatchFile(batFile,
						line -> System.out.println(fmt.format(new Date()) + "\t" + line)));
			} catch (IOException | InterruptedException e) {
				LOGGER.error(taskName + ": testDeepDetect - error running process", e);
			}
		});
		try {
			executor.shutdown();
			boolean finished = executor.awaitTermination(Math.max(1, timeoutHours), TimeUnit.HOURS);
			FileUtils.deleteQuietly(batFile);
			if (!finished) {
				LOGGER.warn(taskName + ": testDeepDetect - timed out after " + timeoutHours + " hour(s)");
				return false;
			}
			boolean ok = (exitCode.get() == 0);
			LOGGER.info(taskName + ": testDeepDetect - exit code " + exitCode.get() + ", success=" + ok);
			return ok;
		} catch (InterruptedException e) {
			LOGGER.error(taskName + ": testDeepDetect - interrupted", e);
			Thread.currentThread().interrupt();
			FileUtils.deleteQuietly(batFile);
			return false;
		}
	}
	
	public static void main(String[] args) {
		DeepDetectParameter par = new DeepDetectParameter();
		par.setInput(new File("Z:\\PRIDE_dataset\\PXD023217_campi\\GUT_DB2MG.faa"));
		par.setOutput(new File("Z:\\PRIDE_dataset\\PXD023217_campi\\GUT_DB2MG.predicted.tsv"));
		par.setMissCleavage(2);
		testDeepDetect(new File("C:\\Users\\daz25hak\\MetaProAI\\Resources\\DeepDetect\\DeepDetect.exe"), par, 2);
	}
}
