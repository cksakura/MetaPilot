package uk.ac.quadram.metapilot.task.dia;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import javafx.application.Platform;
// import javafx.collections.FXCollections; // Not used
import javafx.collections.ObservableList;
import javafx.scene.control.ProgressBar;
import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNParquetReader;
import uk.ac.quadram.metapilot.dbSearch.diann.DiannParameter;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.task.dia.par.MetaParameterDia;
import uk.ac.quadram.metapilot.task.dia.par.MetaParameterDiaRT;
import uk.ac.quadram.metapilot.task.common.gui.ConsolePanel;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.channels.FileChannel;
import java.nio.channels.FileLock;
import java.nio.file.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

public class MetaDiaMonitorTask extends MetaDiaTask {
	private final Path monitorFolderPath;
	private final Path destinationFolderPath;

	private File monitorFolder;
	private File destinationFolder;
	private ExecutorService executorService;
	private HashSet<String> fileNameSet;
	private ObservableList<TableRowData> resultTableItems;
	private ObservableList<String> rawListItems;
	private ConsolePanel consolePanel;
	
	private DiaLibSearchPar diaLibSearchPar;
	private DiannParameter diannPar;
	private File libFolderFile;
	private ArrayList<String> libPathList;
	private ArrayList<String> modelPathList;
	
	// Logger is inherited from MetaDiaTask -> MetaPilotMagTask -> MetaPilotTask as an instance field
	private static final String taskName = "Realtime DIA data analysis";
	
	/**
	 * Helper method to update table row status
	 */
	private void updateTableRowStatus(String fileName, int columnIndex, String status) {
		Platform.runLater(() -> {
			for (TableRowData row : resultTableItems) {
				if (row.matchesFileName(fileName)) {
					switch (columnIndex) {
						case 1: row.setRound1Status(status); break;
						case 2: row.setRound2Status(status); break;
						case 3: row.setRound3Status(status); break;
						case 4: row.setRound4Status(status); break;
						case 5: row.setRound5Status(status); break;
					}
					break;
				}
			}
		});
	}

	public MetaDiaMonitorTask(MetaParameterDiaRT metaPar, MetaResourceParameters mrp, ProgressBar bar1,
			ProgressBar bar2, File monitorFolder, File destinationFolder,
			ObservableList<TableRowData> resultTableItems, ObservableList<String> rawListItems, ConsolePanel consolePanel) {

		super(metaPar, mrp, bar1, bar2);
		this.monitorFolderPath = Paths.get(monitorFolder.getAbsolutePath());
		this.destinationFolderPath = Paths.get(destinationFolder.getAbsolutePath());

		this.monitorFolder = monitorFolder;
		this.destinationFolder = destinationFolder;

		this.resultTableItems = resultTableItems;
		this.rawListItems = rawListItems;
		this.consolePanel = consolePanel;
	}

	protected void initial() {
		this.magDb = ((MetaParameterDiaRT) metaPar).getMagDbItem();
		this.threadCount = metaPar.getThreadCount();

		this.diaLibSearchPar = new DiaLibSearchPar();
		this.diaLibSearchPar.setThreadCount(metaPar.getThreadCount());

		this.diannPar = new DiannParameter();
		this.diannPar.setThreads(metaPar.getThreadCount());

		setTaskCount();
		
		// Debug: Verify DiaNN task is initialized (should be initialized in parent constructor)
		System.out.println("DEBUG: MetaDiaMonitorTask.initial() - diaNNTask: " + (this.diaNNTask != null ? "initialized" : "NULL"));
		System.out.println("DEBUG: MetaDiaMonitorTask.initial() - diannPar: " + (this.diannPar != null ? "initialized" : "NULL"));
	}

	@SuppressWarnings("unchecked")
	@Override
	protected Boolean call() throws IOException {
		
		System.out.println("DEBUG: MetaDiaMonitorTask.call() started");
		System.out.println("DEBUG: Calling initial() method...");
		initial();
		System.out.println("DEBUG: initial() method completed");

		File magDbFile = this.magDb.getCurrentFile();
		this.libFolderFile = new File(magDbFile, "libraries");

		MetaParameterDiaRT diaPar = (MetaParameterDiaRT) metaPar;
		String[] usedLibs = diaPar.getLibrary();
		this.libPathList = new ArrayList<String>();

		for (int i = 0; i < usedLibs.length; i++) {
			if (!usedLibs[i].equals("host")) {
				libPathList.add(usedLibs[i]);
			}
		}

		File modelFolderFile = new File(magDbFile, "models");
		String[] usedModels = diaPar.getModels();
		this.modelPathList = new ArrayList<String>();
		for (int i = 0; i < usedModels.length; i++) {
			File modelFile = new File(new File(modelFolderFile, usedModels[i]), usedModels[i] + "_pep_rank.tsv");
			modelPathList.add(modelFile.getAbsolutePath());
		}

		executorService = Executors.newFixedThreadPool(maxTaskCount);

		this.fileNameSet = new HashSet<String>();

		List<String> diaFileList = new ArrayList<>();
		List<Future<String>> futures = new ArrayList<>();

		File[] files = monitorFolder.listFiles();
		if (files.length > 0) {
			logger.info(
					getTaskName() + ": " + files.length + " files already in the monitor folder, processing started");

			System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": " + files.length
					+ " files already in the monitor folder, processing started");
			
			// List the files that will be processed
			for (File file : files) {
				if (file.isFile()) {
					String fileName = file.getName().toLowerCase();
					if (fileName.endsWith(".raw") || fileName.endsWith(".d") || 
						fileName.endsWith(".dia") || fileName.endsWith(".mzml")) {
						System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": Found existing file: " + file.getName());
						logger.info(getTaskName() + ": Found existing file: " + file.getName());
					}
				}
			}

			for (int i = 0; i < files.length; i++) {
				final int index = i;
				System.out.println("DEBUG: Submitting file " + files[index].getName() + " to executor service");
				Future<String> future = executorService.submit(() -> process(files[index]));
				futures.add(future);
				System.out.println("DEBUG: File " + files[index].getName() + " submitted successfully");
			}
		}

		WatchService watchService = FileSystems.getDefault().newWatchService();
		monitorFolderPath.register(watchService, StandardWatchEventKinds.ENTRY_CREATE);

		while (!isCancelled()) {

			WatchKey key;
			try {
				key = watchService.take();
			} catch (InterruptedException e) {
				Thread.currentThread().interrupt();
				break;
			}

			for (WatchEvent<?> event : key.pollEvents()) {
				WatchEvent.Kind<?> kind = event.kind();
				if (kind == StandardWatchEventKinds.OVERFLOW) {
					continue;
				}

				WatchEvent<Path> ev = (WatchEvent<Path>) event;
				Path filename = ev.context();
				Path filePath = monitorFolderPath.resolve(filename);

				if (Files.isRegularFile(filePath)) {
					String fileName = filename.toString().toLowerCase();
					// Only process supported spectra file types
					if (fileName.endsWith(".raw") || fileName.endsWith(".d") || 
						fileName.endsWith(".dia") || fileName.endsWith(".mzml")) {

						logger.info(getTaskName() + ": found new file " + filePath);
						System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": New file detected: " + filename);

						Future<String> future = executorService.submit(() -> waitForFileCompletionAndProcess(filePath));
						futures.add(future);
					}
				}
			}

			boolean valid = key.reset();
			if (!valid) {
				break;
			}

			MetaParameterDiaRT parDiaRT = ((MetaParameterDiaRT) metaPar);
			boolean combine = parDiaRT.isCombineResult();

			if (parDiaRT.isStopNow()) {
				if (combine) {
					System.out.println(format.format(new Date()) + "\t" + getTaskName()
							+ ": the monitoring is stopped, the processing will stop after all the current tasks are finished. "
							+ "The taxonomic and functional analysis will be performed based on the combined result.");
					logger.info(getTaskName()
							+ ": the monitoring is stopped, the processing will stop after all the current tasks are finished. "
							+ "The taxonomic and functional analysis will be performed based on the combined result.");
				} else {
					System.out.println(format.format(new Date()) + "\t" + getTaskName()
							+ ": the monitoring is stopped, the processing will stop after all the current tasks are finished. "
							+ "The taxonomic and functional analysis will NOT be performed.");
					logger.info(getTaskName()
							+ ": the monitoring is stopped, the processing will stop after all the current tasks are finished. "
							+ "The taxonomic and functional analysis will NOT be performed.");
				}

				break;
			}
			if (parDiaRT.isStopComplete()) {
				int totalCount = parDiaRT.getTaskCount();
				if (combine) {
					System.out.println(format.format(new Date()) + "\t" + getTaskName()
							+ ": the total number of the files was set as " + totalCount
							+ ", the monitoring will stop after all " + totalCount
							+ " files are detected, the processing will stop after all the current tasks are finished."
							+ "The taxonomic and functional analysis will be performed based on the combined result.");
					logger.info(getTaskName() + ": the total number of the files was set as " + totalCount
							+ ", the monitoring will stop after all " + totalCount
							+ " files are detected, the processing will stop after all the current tasks are finished."
							+ "The taxonomic and functional analysis will be performed based on the combined result.");
				} else {

					System.out.println(format.format(new Date()) + "\t" + getTaskName()
							+ ": the total number of the files was set as " + totalCount
							+ ", the monitoring will stop after all " + totalCount
							+ " files are detected, the processing will stop after all the current tasks are finished."
							+ "The taxonomic and functional analysis will NOT be performed.");
					logger.info(getTaskName() + ": the total number of the files was set as " + totalCount
							+ ", the monitoring will stop after all " + totalCount
							+ " files are detected, the processing will stop after all the current tasks are finished."
							+ "The taxonomic and functional analysis will NOT be performed.");

				}

				if (this.fileNameSet.size() >= totalCount) {

					System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": The file count was set as "
							+ totalCount + " and " + this.fileNameSet.size() + " files have been detected."
							+ " The monitoring is stopped, the processing will stop after all the current tasks are finished. "
							+ "The taxonomic and functional analysis will NOT be performed.");
					logger.info(getTaskName() + ": The file count was set as " + totalCount + " and "
							+ this.fileNameSet.size() + " files have been detected."
							+ " The monitoring is stopped, the processing will stop after all the current tasks are finished. "
							+ "The taxonomic and functional analysis will NOT be performed.");

					break;
				}
			}
		}
		executorService.shutdown();

		for (Future<String> future : futures) {
			try {
				String result = future.get();
				if (result != null) {
					diaFileList.add(result);
				}
			} catch (InterruptedException | ExecutionException e) {
				e.printStackTrace();
			}
		}
		this.diaFiles = diaFileList.toArray(String[]::new);
		Arrays.sort(diaFiles);

		System.out.println(
				format.format(new Date()) + "\t" + getTaskName() + ": the processing of all the raw files finished");
		logger.info(getTaskName() + ": the processing of all the raw files finished");

		this.separateResultFolders = new File[fileNameSet.size()];
		int fildId = 0;
		for (String fName : fileNameSet) {
			this.separateResultFolders[fildId++] = new File(destinationFolder, fName);
		}
		this.resultFolderFile = new File(destinationFolder, "combined");

		cleanTempFiles();

		updateProgress(100, 100);

		return true;
	}
	
	private String process(File file) {
		try {
			System.out.println("DEBUG: process() method called for file: " + file.getName());

			String fileAbPath = file.getAbsolutePath();
			String fileName = file.getName();
			
			// Extract base filename for UI display (remove extensions)
			String baseFileName = fileName;
			if (fileName.endsWith(".raw")) {
				baseFileName = fileName.substring(0, fileName.length() - 4);
			} else if (fileName.endsWith(".d")) {
				baseFileName = fileName.substring(0, fileName.length() - 2);
			} else if (fileName.endsWith(".dia")) {
				if (fileName.endsWith(".raw.dia")) {
					baseFileName = fileName.substring(0, fileName.length() - 8);
				} else if (fileName.endsWith(".raw.d")) {
					baseFileName = fileName.substring(0, fileName.length() - 6);
				} else {
					baseFileName = fileName.substring(0, fileName.length() - 4);
				}
			}
			
			// Check if this base filename is already being processed to avoid duplicates
			synchronized (fileNameSet) {
				if (fileNameSet.contains(baseFileName)) {
					System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": Duplicate file detected, skipping: " + fileName + " (base: " + baseFileName + ")");
					logger.info(getTaskName() + ": Duplicate file detected, skipping: " + fileName + " (base: " + baseFileName + ")");
					return null;
				}
				// Add to set immediately to prevent other threads from processing the same base filename
				fileNameSet.add(baseFileName);
			}
			
			// IMMEDIATE UI UPDATE - Add file to UI as soon as processing starts
			final String finalBaseFileName = baseFileName;
			Platform.runLater(() -> {
				// Add to list if not already present
				if (!rawListItems.contains(finalBaseFileName)) {
					rawListItems.add(finalBaseFileName);
					System.out.println("DEBUG: IMMEDIATE UI UPDATE - Added file to list: " + finalBaseFileName);
				}
				
				// Add to table if not already present
				boolean existsInTable = resultTableItems.stream().anyMatch(row -> row.matchesFileName(finalBaseFileName));
				if (!existsInTable) {
					TableRowData newRow = new TableRowData(finalBaseFileName, "Processing...", "Waiting", "Waiting", "Waiting", "Waiting");
					resultTableItems.add(newRow);
					System.out.println("DEBUG: IMMEDIATE UI UPDATE - Added file to table: " + finalBaseFileName);
				}
			});

			logger.info(getTaskName() + ": processing raw files " + fileName + " started");
			System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": processing raw files " + fileName
					+ " started");

			// Use the DiaNN task initialized in the constructor, not create a new one
			System.out.println(
					"DEBUG: Using DiaNN task from constructor: " + (this.diaNNTask != null ? "available" : "NULL"));
			System.out.println(
					"DEBUG: Using diannPar from constructor: " + (this.diannPar != null ? "available" : "NULL"));

			// Create a new DiannParameter for this specific conversion task
			DiannParameter conversionDiannPar = new DiannParameter();
			conversionDiannPar.setThreads(threadCount / maxTaskCount + 1);

			File diaFile = new File(fileAbPath);
			if (fileName.endsWith(".raw") || fileName.endsWith(".RAW")) {

				diaFile = new File(file + ".dia");
				fileName = fileName.substring(0, fileName.length() - 4);
				// fileNameSet.add() already done earlier in the method

				if (!diaFile.exists()) {
					logger.info(getTaskName() + ": converting raw files " + fileName + " started");
					System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": converting raw files "
							+ fileName + " started");
					System.out.println("DEBUG: Input file: " + fileAbPath);
					System.out.println("DEBUG: Expected output: " + diaFile.getAbsolutePath());
					System.out.println("DEBUG: DiaNN task: " + (this.diaNNTask != null ? "initialized" : "NULL"));
					System.out
							.println("DEBUG: DiannParameter: " + (conversionDiannPar != null ? "initialized" : "NULL"));

					this.diaNNTask.addTask(conversionDiannPar, new String[] { fileAbPath });
					System.out.println("DEBUG: DiaNN task added, about to run...");
					this.diaNNTask.run(1);
					System.out.println("DEBUG: DiaNN task run completed");

					// Check if conversion was successful
					if (diaFile.exists()) {
						logger.info(getTaskName() + ": converting raw files " + fileName + " finished successfully");
						System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": converting raw files "
								+ fileName + " finished successfully");
						System.out.println("DEBUG: Output file created: " + diaFile.getAbsolutePath() + " (size: "
								+ diaFile.length() + " bytes)");
					} else {
						logger.error(getTaskName() + ": converting raw files " + fileName
								+ " FAILED - output file not created");
						System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": converting raw files "
								+ fileName + " FAILED - output file not created");
						System.out.println("DEBUG: Expected output file does not exist: " + diaFile.getAbsolutePath());
						System.out.println("DEBUG: Check DiaNN installation and .NET Runtime");
						System.out.println(
								"DEBUG: Download .NET Runtime 8.0.407 or later from: https://dotnet.microsoft.com/en-us/download/dotnet/8.0");
					}
				} else {
					logger.info(getTaskName() + ": .dia file already exists for " + fileName);
					System.out.println(format.format(new Date()) + "\t" + getTaskName()
							+ ": .dia file already exists for " + fileName);
					System.out.println("DEBUG: Using existing file: " + diaFile.getAbsolutePath() + " (size: "
							+ diaFile.length() + " bytes)");
				}

				// Continue processing the .dia file instead of returning null
				// The rest of the method will process the diaFile

			} else if (fileName.endsWith(".d") || fileName.endsWith(".D")) {
				fileName = fileName.substring(0, fileName.length() - 2);
				// fileNameSet.add() already done earlier in the method
			} else if (fileName.endsWith(".dia") || fileName.endsWith(".DIA")) {
				if (fileName.endsWith(".raw.dia") || fileName.endsWith(".RAW.DIA")) {
					fileName = fileName.substring(0, fileName.length() - 8);
				} else if (fileName.endsWith(".raw.d") || fileName.endsWith(".RAW.D")) {
					fileName = fileName.substring(0, fileName.length() - 6);
				} else {
					fileName = fileName.substring(0, fileName.length() - 4);
				}
				// fileNameSet.add() already done earlier in the method
			} else {
				// Only process supported spectra file types
				System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": Unsupported file type, skipping: " + fileName);
				logger.info(getTaskName() + ": Unsupported file type, skipping: " + fileName);
				return null;
			}

			final String finalFileName = baseFileName; // Use the base filename extracted earlier
			
			// UI updates already happened immediately at the start of processing

			// JavaFX models updated directly via Platform.runLater()
			System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": Added file to JavaFX models: "
					+ finalFileName);

			File resultFolderFile = new File(destinationFolderPath.toAbsolutePath().toString(), fileName);
			if (!resultFolderFile.exists()) {
				resultFolderFile.mkdir();
			}

			HashSet<String> initialLibPepSet = new HashSet<String>();
			String[] libPath = new String[libPathList.size()];
			for (int i = 0; i < libPathList.size(); i++) {
				String usedLibIString = libPathList.get(i);
				File libTsvFile = new File(new File(libFolderFile, usedLibIString), usedLibIString + ".parquet");
				File libFile = new File(new File(libFolderFile, usedLibIString), usedLibIString + ".speclib");
				libPath[i] = libFile.getAbsolutePath();

				initialLibPepSet.addAll(DiaNNParquetReader.getPeptideSet(libTsvFile.getAbsolutePath()));
			}

			File firstFile = new File(resultFolderFile, "firstSearch");
			if (!firstFile.exists()) {
				firstFile.mkdirs();
			}

			updateTableRowStatus(finalFileName, 1, "Working");
			logger.info(getTaskName() + ": Updated " + finalFileName + " Round 1 status to Working");

			// JavaFX models updated directly via Platform.runLater()
			System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": Updated " + finalFileName
					+ " Round 1 status to Working");

			File round1TsvFile = new File(firstFile, "firstSearch.parquet");
			if (!round1TsvFile.exists()) {

				logger.info(getTaskName() + ": search raw files " + fileName + " round 1 started");
				System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": search raw files " + fileName
						+ " round 1 started");

				int batchThreadCount = threadCount / maxTaskCount + 1;
				DiaLibSearchPar diaLibSearchPar = new DiaLibSearchPar();
				diaLibSearchPar.setDiaFiles(new String[] { diaFile.getAbsolutePath() });
				diaLibSearchPar.setThreadCount(batchThreadCount);
				diaLibSearchPar.setMbr(false);
				diaLibSearchPar.setqValue(0.01);

				DiannParameter diannPar = new DiannParameter();
				diannPar.setThreads(batchThreadCount);
				diannPar.setMatrices(false);
				diannPar.setRelaxed_prot_inf(true);
				diaNNTask.addTask(diannPar, diaLibSearchPar, round1TsvFile.getAbsolutePath(), libPath);
				diaNNTask.run(12);

				logger.info(getTaskName() + ": search raw files " + fileName + " round 1 finished");
				System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": converting raw files "
						+ fileName + " finished");
			} else {
				logger.info(
						getTaskName() + ": round 1 search result for " + fileName + " existed, go to the next step");
				System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": round 1 search result for "
						+ fileName + " existed, go to the next step");
			}

			if (round1TsvFile.exists()) {
				updateTableRowStatus(finalFileName, 1, "Done");
			} else {
				updateTableRowStatus(finalFileName, 1, "Failed");
				return null;
			}

			File modelFolderFile = new File(resultFolderFile, "model");
			if (!modelFolderFile.exists()) {
				modelFolderFile.mkdirs();
			}

			File rankedPepFile = new File(modelFolderFile, "model_pep_rank.tsv");
			if (!rankedPepFile.exists()) {

				logger.info(getTaskName() + ": creating model for " + fileName + " started");
				System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": creating model for " + fileName
						+ " started");

				HashMap<String, Double> pepIntenMap = new HashMap<String, Double>();
				try {
					pepIntenMap = DiaNNParquetReader.getPeptideMap(round1TsvFile.getAbsolutePath());
				} catch (IOException e) {
					logger.error(getTaskName() + ": error in reading parquet file " + round1TsvFile.getAbsolutePath(), e);
					System.out.println(format.format(new Date()) + "\t" + getTaskName()
							+ ": error in reading parquet file " + round1TsvFile.getAbsolutePath());
				}
				int batchThreadCount = threadCount / maxTaskCount + 1;
				DiaPyTorchModelTask diaModelTask = new DiaPyTorchModelTask((MetaParameterDia) metaPar, mrp, secondBar,
						pepIntenMap, modelFolderFile, batchThreadCount);
				diaModelTask.execute();

				updateTableRowStatus(finalFileName, 2, "Working");

				try {
					if (diaModelTask.get()) {
						logger.info(getTaskName() + ": creating model for " + fileName + " finished");
						System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": creating model for "
								+ fileName + " finished");

					} else {
						logger.info(getTaskName() + ": creating model for " + fileName + " failed");
						System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": creating model for "
								+ fileName + " failed");
					}

				} catch (ExecutionException e) {
					// TODO Auto-generated catch block
					logger.error(getTaskName() + ": creating model for " + fileName + " failed", e);
					System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": creating model for "
							+ fileName + " failed");
				}
			} else {
				logger.info(getTaskName() + ": model for " + fileName + " existed, go to the next step");
				System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": model for " + fileName
						+ " existed, go to the next step");
			}

			if (rankedPepFile.exists()) {
				updateTableRowStatus(finalFileName, 2, "Done");
			} else {
				updateTableRowStatus(finalFileName, 2, "Failed");
				return null;
			}

			ArrayList<String> modelListI = new ArrayList<String>(modelPathList);
			modelListI.add(rankedPepFile.getAbsolutePath());
			String[] modelPath = modelListI.toArray(String[]::new);

			File round1DbFile = new File(firstFile, "firstSearch.fasta");
			if (!round1DbFile.exists() || round1DbFile.length() == 0) {
				HashSet<String> pepSet1 = null;
				try {
					pepSet1 = DiaNNParquetReader.getPeptideSet(round1TsvFile.getAbsolutePath());
				} catch (Exception e) {
					logger.error(getTaskName() + ": error in reading " + round1TsvFile, e);
					System.err.println(
							format.format(new Date()) + "\t" + getTaskName() + ": error in reading " + round1TsvFile);
					return null;
				}
				int batchThreadCount = threadCount / maxTaskCount + 1;
				HashSet<String> genomeSet1 = getGenomeSet(firstFile, 1.0, pepSet1, 5, false, batchThreadCount);

				logger.info(getTaskName() + ": generating database for " + fileName + ", " + pepSet1.size()
						+ " peptide sequences were identified, " + genomeSet1.size() + " MAGs were used in this step");
				System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": generating database for "
						+ fileName + ", " + pepSet1.size() + " peptide sequences were identified, " + genomeSet1.size()
						+ " MAGs were used in this step");

				dbFromModelPeptide(modelPath, genomeSet1, pepSet1, initialLibPepSet, round1DbFile, 2500000);
			}

			File secondFile = new File(resultFolderFile, "secondSearch");
			if (!secondFile.exists()) {
				secondFile.mkdirs();
			}

			updateTableRowStatus(finalFileName, 3, "Working");

			File round2TsvFile = new File(secondFile, "secondSearch.parquet");
			if (!round2TsvFile.exists()) {

				logger.info(getTaskName() + ": search raw files " + fileName + " round 2 started");
				System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": search raw files " + fileName
						+ " round 2 started");

				File round2LibFile = new File(secondFile, "secondSearch.speclib");
				if (!round2LibFile.exists()) {
					File tempLibFile = new File(secondFile, "secondSearch.predicted.speclib");
					if (tempLibFile.exists()) {
						tempLibFile.renameTo(round2LibFile);
					} else {
						int batchThreadCount = threadCount / maxTaskCount + 1;
						diannPar.setThreads(batchThreadCount);
						diaNNTask.addTask(diannPar, round1DbFile.getAbsolutePath(), round2LibFile.getAbsolutePath(),
								false);
						diaNNTask.run(12);
					}
				}

				if (!round2LibFile.exists()) {
					File tempLibFile = new File(secondFile, "secondSearch.predicted.speclib");
					if (tempLibFile.exists()) {
						tempLibFile.renameTo(round2LibFile);
					} else {
						logger.error(getTaskName() + ": library file " + round2LibFile + " is not found");
						System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": library file "
								+ round2LibFile + " is not found");

					}
				}
				int batchThreadCount = threadCount / maxTaskCount + 1;
				DiaLibSearchPar diaLibSearchPar = new DiaLibSearchPar();
				diaLibSearchPar.setDiaFiles(new String[] { diaFile.getAbsolutePath() });
				diaLibSearchPar.setThreadCount(batchThreadCount);
				diaLibSearchPar.setMbr(false);
				diaLibSearchPar.setqValue(0.01);

				DiannParameter diannPar = new DiannParameter();
				diannPar.setThreads(batchThreadCount);
				diannPar.setMatrices(false);
				diannPar.setRelaxed_prot_inf(true);
				diaNNTask.addTask(diannPar, new String[] { diaFile.getAbsolutePath() }, round2TsvFile.getAbsolutePath(),
						round2LibFile.getAbsolutePath());
				diaNNTask.run(12);

				logger.info(getTaskName() + ": search raw files " + fileName + " round 2 finished");
				System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": search raw files " + fileName
						+ " round 2 finished");

			} else {
				logger.info(
						getTaskName() + ": round 2 search result for " + fileName + " existed, go to the next step");
				System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": round 2 search result for "
						+ fileName + " existed, go to the next step");
			}

			if (round2TsvFile.exists()) {
				updateTableRowStatus(finalFileName, 3, "Done");
			} else {
				updateTableRowStatus(finalFileName, 3, "Failed");
				return null;
			}

			File round2DbFile = new File(secondFile, "secondSearch.fasta");
			if (!round2DbFile.exists() || round2DbFile.length() == 0) {

				HashSet<String> pepSet2 = new HashSet<String>();
				try {
					pepSet2.addAll(DiaNNParquetReader.getPeptideSet(round1TsvFile.getAbsolutePath()));
				} catch (Exception e) {
					logger.error(getTaskName() + ": error in reading " + round1TsvFile, e);
					System.err.println(
							format.format(new Date()) + "\t" + getTaskName() + ": error in reading " + round1TsvFile);
					return null;
				}
				try {
					pepSet2.addAll(DiaNNParquetReader.getPeptideSet(round2TsvFile.getAbsolutePath()));
				} catch (Exception e) {
					logger.error(getTaskName() + ": error in reading " + round2TsvFile, e);
					System.err.println(
							format.format(new Date()) + "\t" + getTaskName() + ": error in reading " + round2TsvFile);
					return null;
				}

				int batchThreadCount = threadCount / maxTaskCount + 1;
				HashSet<String> genomeSet2 = getGenomeSet(secondFile, 1.0, pepSet2, 5, false, batchThreadCount);

				try (BufferedReader firstLibReader = new BufferedReader(
						new FileReader(new File(firstFile, "firstSearch.fasta")))) {

					String pline = null;
					while ((pline = firstLibReader.readLine()) != null) {
						if (!pline.startsWith(">")) {
							initialLibPepSet.add(pline);
						}
					}
					firstLibReader.close();

				} catch (Exception e) {
					logger.error(getTaskName() + ": error in reading fasta file from " + firstFile, e);
					System.err.println(format.format(new Date()) + "\t" + getTaskName()
							+ ": error in reading fasta file from " + firstFile);
					return null;
				}

				logger.info(getTaskName() + ": generating database for " + fileName + ", " + pepSet2.size()
						+ " peptide sequences were identified, " + genomeSet2.size()
						+ " MAGs were used in this step, the size of the excluded set is " + initialLibPepSet.size());
				System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": generating database for "
						+ fileName + ", " + pepSet2.size() + " peptide sequences were identified, " + genomeSet2.size()
						+ " MAGs were used in this step, the size of the excluded set is " + initialLibPepSet.size());

				dbFromModelPeptide(modelPath, genomeSet2, pepSet2, initialLibPepSet, round2DbFile, 2000000);
			}

			File thirdFile = new File(resultFolderFile, "thirdSearch");
			if (!thirdFile.exists()) {
				thirdFile.mkdirs();
			}

			updateTableRowStatus(finalFileName, 4, "Working");

			File round3TsvFile = new File(thirdFile, "thirdSearch.parquet");
			if (!round3TsvFile.exists()) {
				logger.info(getTaskName() + ": search raw files " + fileName + " round 3 started");
				System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": search raw files " + fileName
						+ " round 3 started");
				File round3LibFile = new File(thirdFile, "thirdSearch.speclib");
				if (!round3LibFile.exists()) {
					File tempLibFile = new File(thirdFile, "thirdSearch.predicted.speclib");
					if (tempLibFile.exists()) {
						tempLibFile.renameTo(round3LibFile);
					} else {
						int batchThreadCount = threadCount / maxTaskCount + 1;
						diannPar.setThreads(batchThreadCount);
						diaNNTask.addTask(diannPar, round2DbFile.getAbsolutePath(), round3LibFile.getAbsolutePath(),
								false);
						diaNNTask.run(12);
					}
				}
				if (!round3LibFile.exists()) {
					File tempLibFile = new File(thirdFile, "thirdSearch.predicted.speclib");
					if (tempLibFile.exists()) {
						tempLibFile.renameTo(round3LibFile);
					} else {
						logger.error(getTaskName() + ": library file " + round3LibFile + " is not found");
						System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": library file "
								+ round3LibFile + " is not found");
					}
				}
				int batchThreadCount = threadCount / maxTaskCount + 1;
				DiaLibSearchPar diaLibSearchPar = new DiaLibSearchPar();
				diaLibSearchPar.setThreadCount(batchThreadCount);
				diaLibSearchPar.setMbr(false);
				diaLibSearchPar.setqValue(0.01);

				DiannParameter diannPar = new DiannParameter();
				diannPar.setThreads(batchThreadCount);
				diannPar.setMatrices(false);
				diannPar.setRelaxed_prot_inf(true);
				diaNNTask.addTask(diannPar, new String[] { diaFile.getAbsolutePath() }, round3TsvFile.getAbsolutePath(),
						round3LibFile.getAbsolutePath());
				diaNNTask.run(12);

				logger.info(getTaskName() + ": search raw files " + fileName + " round 3 finished");
				System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": search raw files " + fileName
						+ " round 3 finished");
			} else {
				logger.info(
						getTaskName() + ": round 3 search result for " + fileName + " existed, go to the next step");
				System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": round 3 search result for "
						+ fileName + " existed, go to the next step");
			}

			if (round3TsvFile.exists()) {
				updateTableRowStatus(finalFileName, 4, "Done");
			} else {
				updateTableRowStatus(finalFileName, 4, "Failed");
				return null;
			}

			File round3DbFile = new File(thirdFile, "thirdSearch.fasta");
			if (!round3DbFile.exists() || round3DbFile.length() == 0) {

				HashSet<String> pepSet3 = new HashSet<String>();
				try {
					pepSet3.addAll(DiaNNParquetReader.getPeptideSet(round1TsvFile.getAbsolutePath()));
				} catch (Exception e) {
					logger.error(getTaskName() + ": error in reading " + round1TsvFile, e);
					System.err.println(
							format.format(new Date()) + "\t" + getTaskName() + ": error in reading " + round1TsvFile);
				}
				try {
					pepSet3.addAll(DiaNNParquetReader.getPeptideSet(round2TsvFile.getAbsolutePath()));
				} catch (Exception e) {
					logger.error(getTaskName() + ": error in reading " + round2TsvFile, e);
					System.err.println(
							format.format(new Date()) + "\t" + getTaskName() + ": error in reading " + round2TsvFile);
				}
				try {
					pepSet3.addAll(DiaNNParquetReader.getPeptideSet(round3TsvFile.getAbsolutePath()));
				} catch (Exception e) {
					logger.error(getTaskName() + ": error in reading " + round3TsvFile, e);
					System.err.println(
							format.format(new Date()) + "\t" + getTaskName() + ": error in reading " + round3TsvFile);
				}

				int batchThreadCount = threadCount / maxTaskCount + 1;
				HashSet<String> genomeSet3 = getGenomeSet(thirdFile, 1.0, pepSet3, 5, false, batchThreadCount);

				try (BufferedReader firstLibReader = new BufferedReader(
						new FileReader(new File(firstFile, "firstSearch.fasta")))) {

					String pline = null;
					while ((pline = firstLibReader.readLine()) != null) {
						if (!pline.startsWith(">")) {
							initialLibPepSet.add(pline);
						}
					}
					firstLibReader.close();

				} catch (Exception e) {
					logger.error(getTaskName() + ": error in reading fasta file from " + firstFile, e);
					System.err.println(format.format(new Date()) + "\t" + getTaskName()
							+ ": error in reading fasta file from " + firstFile);
					return null;
				}
				try (BufferedReader secondLibReader = new BufferedReader(
						new FileReader(new File(secondFile, "secondSearch.fasta")))) {

					String pline = null;
					while ((pline = secondLibReader.readLine()) != null) {
						if (!pline.startsWith(">")) {
							initialLibPepSet.add(pline);
						}
					}
					secondLibReader.close();

				} catch (Exception e) {
					logger.error(getTaskName() + ": error in reading fasta file from " + firstFile, e);
					System.err.println(format.format(new Date()) + "\t" + getTaskName()
							+ ": error in reading fasta file from " + firstFile);
					return null;
				}

				logger.info(getTaskName() + ": generating database for " + fileName + ", " + genomeSet3.size()
						+ " MAGs were used in this step, the size of the excluded set is " + initialLibPepSet.size());
				System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": generating database for "
						+ fileName + ", " + genomeSet3.size()
						+ " MAGs were used in this step, the size of the excluded set is " + initialLibPepSet.size());

				dbFromModelPeptide(modelPath, genomeSet3, pepSet3, initialLibPepSet, round3DbFile, 2500000);
			}

			updateTableRowStatus(finalFileName, 5, "Working");

			File lastTsvFile = new File(resultFolderFile, "lastSearch.parquet");
			if (!lastTsvFile.exists()) {

				logger.info(getTaskName() + ": search raw files " + fileName + " the last round started");
				System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": search raw files " + fileName
						+ " the last round started");

				File round4LibFile = new File(resultFolderFile, "lastSearch.speclib");
				if (!round4LibFile.exists()) {
					File tempLibFile = new File(resultFolderFile, "lastSearch.predicted.speclib");
					if (tempLibFile.exists()) {
						tempLibFile.renameTo(round4LibFile);
					} else {
						int batchThreadCount = threadCount / maxTaskCount + 1;
						diannPar.setThreads(batchThreadCount);
						diaNNTask.addTask(diannPar, round3DbFile.getAbsolutePath(), round4LibFile.getAbsolutePath(),
								false);
						diaNNTask.run(12);
					}
				}
				if (!round4LibFile.exists()) {
					File tempLibFile = new File(resultFolderFile, "lastSearch.predicted.speclib");
					if (tempLibFile.exists()) {
						tempLibFile.renameTo(round4LibFile);
					} else {
						logger.error(getTaskName() + ": library file " + round4LibFile + " is not found");
						System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": library file "
								+ round4LibFile + " is not found");
					}
				}
				int batchThreadCount = threadCount / maxTaskCount + 1;
				DiaLibSearchPar diaLibSearchPar = new DiaLibSearchPar();
				diaLibSearchPar.setThreadCount(batchThreadCount);
				diaLibSearchPar.setMbr(false);
				diaLibSearchPar.setqValue(0.01);

				DiannParameter diannPar = new DiannParameter();
				diannPar.setThreads(batchThreadCount);
				diannPar.setMatrices(false);
				diannPar.setRelaxed_prot_inf(true);
				diaNNTask.addTask(diannPar, new String[] { diaFile.getAbsolutePath() }, lastTsvFile.getAbsolutePath(),
						round4LibFile.getAbsolutePath());
				diaNNTask.run(12);

				logger.info(getTaskName() + ": search raw files " + fileName + " the last round finished");
				System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": search raw files " + fileName
						+ " the last round finished");
			} else {
				logger.info(getTaskName() + ": the last round search result for " + fileName
						+ " existed, go to the next step");
				System.out.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": the last round search result for " + fileName + " existed, go to the next step");
			}

			if (lastTsvFile.exists()) {
				updateTableRowStatus(finalFileName, 5, "Done");
			} else {
				updateTableRowStatus(finalFileName, 5, "Failed");
				return null;
			}

			int count = 0;
			File[] destFiles = destinationFolder.listFiles();
			for (int i = 0; i < destFiles.length; i++) {
				if (destFiles[i].isDirectory()) {
					File resultFileI = new File(destFiles[i], "lastSearch.parquet");
					if (resultFileI.exists() && resultFileI.length() > 0) {
						count++;
					}
				}
			}

			logger.info(getTaskName() + ": " + destFiles.length + " files are finished, the total number of file is "
					+ this.rawListItems.size());
			System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": " + destFiles.length
					+ " files are finished, the total number of file is " + this.rawListItems.size());

			if (count <= this.rawListItems.size()) {
				updateProgress((int) ((double) count / (double) this.rawListItems.size() * 100.0), 100);
			}

			return diaFile.getAbsolutePath();

		} catch (InterruptedException e) {
			e.printStackTrace();
			Thread.currentThread().interrupt();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		return null;
	}

	private String waitForFileCompletionAndProcess(Path filePath) {
		try {
			long previousSize = -1;
			while (true) {
				// Check the file size at regular intervals
				long currentSize = Files.size(filePath);
				if (currentSize == previousSize) {
					// Try to acquire an exclusive lock to see if any other process is still writing
					try (FileChannel fileChannel = FileChannel.open(filePath, StandardOpenOption.WRITE)) {
						FileLock lock = fileChannel.tryLock();
						if (lock != null) {
							lock.release(); // Release the lock immediately
							break; // File is not in use
						}
					} catch (Exception e) {
						// The file is still in use or another process is writing to it
						logger.info(getTaskName() + ": detect file " + filePath.toAbsolutePath().toString()
								+ ", waiting for the copy complete");
						System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": detect file "
								+ filePath.toAbsolutePath().toString() + ", waiting for the copy complete");
					}
				}
				previousSize = currentSize;
				Thread.sleep(1000); // Wait for 1 second before checking again
			}

			return process(new File(filePath.toAbsolutePath().toString()));

		} catch (IOException | InterruptedException e) {
			e.printStackTrace();
			Thread.currentThread().interrupt();
		}

		return null;
	}

	@Override
	protected void done() {
		try {
			if (!isCancelled()) {
				get(); // Trigger any exceptions thrown in doInBackground
				// Update GUI with completion status
				System.out.println("Monitoring completed. Total files processed: " + rawListItems.size());
			}
		} catch (Exception e) {
			e.printStackTrace();
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
		return logger;
	}
}
