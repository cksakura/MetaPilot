package uk.ac.quadram.metapilot.task.dia;

import org.apache.logging.log4j.Logger;

import javafx.application.Platform;
// import javafx.collections.FXCollections; // Not used
import javafx.collections.ObservableList;
import javafx.scene.control.ProgressBar;
import uk.ac.quadram.metapilot.core.dia.MagPredictedPeptideProteinMapper;
import uk.ac.quadram.metapilot.core.mag.MagGenomeSelectionService;
import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;
import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNParquetReader;
import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNPrecursor;
import uk.ac.quadram.metapilot.dbSearch.diann.DiannParameter;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.task.dia.par.MetaParameterDiaRT;
import uk.ac.quadram.metapilot.task.common.gui.ConsolePanel;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;
import uk.ac.quadram.metapilot.task.mag.MarkerGenomeSelector;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.channels.FileChannel;
import java.nio.channels.FileLock;
import java.nio.file.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
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
	
	// Logger is inherited from MetaDiaTask -> MetaPilotMagTask -> MetaPilotTask as an instance field
	private static final String taskName = "Realtime DIA data analysis";
	
	/**
	 * Helper method to update table row status
	 */
	private void updateTableStatus(String fileName, String status, String details) {
		Platform.runLater(() -> {
			for (TableRowData row : resultTableItems) {
				if (row.matchesFileName(fileName)) {
					row.setStatus(status);
					row.setDetails(details);
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

		executorService = Executors.newFixedThreadPool(maxTaskCount);

		this.fileNameSet = new HashSet<String>();

		List<String> diaFileList = new ArrayList<>();
		List<Future<String>> futures = new ArrayList<>();

		File[] files = monitorFolder.listFiles();
		if (files != null && files.length > 0) {
			logger.info(
					getTaskName() + ": " + files.length + " files already in the monitor folder, processing started");

			System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": " + files.length
					+ " files already in the monitor folder, processing started");
			
			// List the files that will be processed
			for (File file : files) {
				if (isSupportedInput(file)) {
					System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": Found existing file: " + file.getName());
					logger.info(getTaskName() + ": Found existing file: " + file.getName());
				}
			}

			for (int i = 0; i < files.length; i++) {
				if (!isSupportedInput(files[i])) {
					continue;
				}
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

				if (isSupportedInput(filePath)) {
					logger.info(getTaskName() + ": found new file " + filePath);
					System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": New file detected: " + filename);

					Future<String> future = executorService.submit(() -> waitForFileCompletionAndProcess(filePath));
					futures.add(future);
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
			String fileName = file.getName();
			String baseFileName = getBaseName(fileName);
			if (baseFileName == null || baseFileName.isEmpty()) {
				return null;
			}

			synchronized (fileNameSet) {
				if (fileNameSet.contains(baseFileName)) {
					logger.info(getTaskName() + ": Duplicate file detected, skipping {}", fileName);
					return null;
				}
				fileNameSet.add(baseFileName);
			}

			final String uiName = baseFileName;
			Platform.runLater(() -> {
				if (!rawListItems.contains(uiName)) {
					rawListItems.add(uiName);
				}
				boolean existsInTable = resultTableItems.stream().anyMatch(row -> row.matchesFileName(uiName));
				if (!existsInTable) {
					resultTableItems.add(new TableRowData(uiName, "Processing", "Queued"));
				}
			});
			updateTableStatus(baseFileName, "Processing", "Preparing input");

			File diaFile = toDiaInput(file);
			if (diaFile == null || !diaFile.exists()) {
				updateTableStatus(baseFileName, "Failed", "Input conversion failed");
				return null;
			}

			File perFileFolder = new File(destinationFolderPath.toFile(), baseFileName);
			if (!perFileFolder.exists()) {
				perFileFolder.mkdirs();
			}

			String rtMode = ((MetaParameterDiaRT) metaPar).getRtWorkflowMode();
			if (rtMode == null || rtMode.trim().isEmpty()) {
				rtMode = "SEED_DIRECT";
			}
			File finalParquet = new File(perFileFolder, baseFileName + ".parquet");

			if (rtMode.startsWith("MARKER_")) {
				updateTableStatus(baseFileName, "Processing", "Marker search");
				File markerDir = new File(perFileFolder, "marker");
				markerDir.mkdirs();
				File markerParquet = new File(markerDir, "marker.parquet");
				File markerLib = magDb != null ? magDb.resolveMarker4LibraryFile() : null;
				if (markerLib == null || !markerLib.exists()) {
					updateTableStatus(baseFileName, "Failed", "Marker4 library missing");
					return null;
				}
				if (!markerParquet.exists()) {
					if (!runNoQuantSearch(diaFile.getAbsolutePath(), markerParquet.getAbsolutePath(),
							new String[] { markerLib.getAbsolutePath() })) {
						updateTableStatus(baseFileName, "Failed", "Marker search failed");
						return null;
					}
				}
				updateTableStatus(baseFileName, "Processing", "Genome selection");
				String kneeMode = "MARKER_FOCUS".equals(rtMode) ? MarkerGenomeSelector.KNEE_MODE_FOCUS
						: MarkerGenomeSelector.KNEE_MODE_COMPREHENSIVE;
				HashSet<String> genomes = selectGenomesFromParquet(markerParquet, markerDir, kneeMode);
				if (genomes.isEmpty()) {
					updateTableStatus(baseFileName, "Failed", "No genomes selected");
					return null;
				}
				updateTableStatus(baseFileName, "Processing", "Building per-file library");
				File selectedFasta = new File(markerDir, baseFileName + ".fasta");
				if (!writeFastaFromGenomeSet(selectedFasta, genomes)) {
					updateTableStatus(baseFileName, "Failed", "Failed writing FASTA");
					return null;
				}
				File selectedLib = new File(markerDir, baseFileName + ".speclib");
				if (!buildSpeclibFromFasta(selectedFasta, selectedLib)) {
					updateTableStatus(baseFileName, "Failed", "Failed building speclib");
					return null;
				}
				updateTableStatus(baseFileName, "Processing", "Final identification search");
				if (!finalParquet.exists()) {
					if (!runNoQuantSearch(diaFile.getAbsolutePath(), finalParquet.getAbsolutePath(),
							new String[] { selectedLib.getAbsolutePath() })) {
						updateTableStatus(baseFileName, "Failed", "Final search failed");
						return null;
					}
				}
			} else {
				File[] seedLibs = resolveSeedLibraries();
				if (seedLibs.length == 0) {
					updateTableStatus(baseFileName, "Failed", "No seed libraries selected");
					return null;
				}
				File seedDir = new File(perFileFolder, "seed");
				seedDir.mkdirs();
				File seedParquet = new File(seedDir, "seed.parquet");
				updateTableStatus(baseFileName, "Processing", "Seed search");
				if (!seedParquet.exists()) {
					String[] seedLibPaths = new String[seedLibs.length];
					for (int i = 0; i < seedLibs.length; i++) {
						seedLibPaths[i] = seedLibs[i].getAbsolutePath();
					}
					if (!runNoQuantSearch(diaFile.getAbsolutePath(), seedParquet.getAbsolutePath(), seedLibPaths)) {
						updateTableStatus(baseFileName, "Failed", "Seed search failed");
						return null;
					}
				}

				if ("SEED_DIRECT".equals(rtMode)) {
					if (!finalParquet.exists()) {
						Files.copy(seedParquet.toPath(), finalParquet.toPath(), StandardCopyOption.REPLACE_EXISTING);
					}
				} else {
					updateTableStatus(baseFileName, "Processing", "Genome selection");
					String kneeMode = "SEED_FOCUS".equals(rtMode) ? MarkerGenomeSelector.KNEE_MODE_FOCUS
							: MarkerGenomeSelector.KNEE_MODE_COMPREHENSIVE;
					HashSet<String> genomes = selectGenomesFromParquet(seedParquet, seedDir, kneeMode);
					if (genomes.isEmpty()) {
						updateTableStatus(baseFileName, "Failed", "No genomes selected");
						return null;
					}
					updateTableStatus(baseFileName, "Processing", "Building per-file library");
					File selectedFasta = new File(seedDir, baseFileName + ".fasta");
					if (!writeFastaFromGenomeSet(selectedFasta, genomes)) {
						updateTableStatus(baseFileName, "Failed", "Failed writing FASTA");
						return null;
					}
					File selectedLib = new File(seedDir, baseFileName + ".speclib");
					if (!buildSpeclibFromFasta(selectedFasta, selectedLib)) {
						updateTableStatus(baseFileName, "Failed", "Failed building speclib");
						return null;
					}
					updateTableStatus(baseFileName, "Processing", "Final identification search");
					if (!finalParquet.exists()) {
						if (!runNoQuantSearch(diaFile.getAbsolutePath(), finalParquet.getAbsolutePath(),
								new String[] { selectedLib.getAbsolutePath() })) {
							updateTableStatus(baseFileName, "Failed", "Final search failed");
							return null;
						}
					}
				}
			}

			if (!finalParquet.exists()) {
				updateTableStatus(baseFileName, "Failed", "Final parquet not found");
				return null;
			}

			updateTableStatus(baseFileName, "Done", "Result: " + finalParquet.getName());
			updateProgress(getFinishedCount(), Math.max(1, rawListItems.size()));
			return diaFile.getAbsolutePath();
		} catch (Exception e) {
			logger.error(getTaskName() + ": failed processing {}", file.getAbsolutePath(), e);
			String base = getBaseName(file.getName());
			if (base != null) {
				updateTableStatus(base, "Failed", e.getMessage() == null ? "Unknown error" : e.getMessage());
			}
			return null;
		}
	}

	private boolean isSupportedInput(File file) {
		if (file == null || !file.exists()) {
			return false;
		}
		String lower = file.getName().toLowerCase();
		if (lower.endsWith(".d")) {
			return file.isDirectory();
		}
		return file.isFile() && (lower.endsWith(".raw") || lower.endsWith(".dia") || lower.endsWith(".mzml"));
	}

	private boolean isSupportedInput(Path filePath) {
		String lower = filePath.getFileName().toString().toLowerCase();
		if (!(lower.endsWith(".raw") || lower.endsWith(".d") || lower.endsWith(".dia") || lower.endsWith(".mzml"))) {
			return false;
		}
		return Files.exists(filePath);
	}

	private String getBaseName(String fileName) {
		if (fileName == null) {
			return null;
		}
		String lower = fileName.toLowerCase();
		if (lower.endsWith(".raw")) {
			return fileName.substring(0, fileName.length() - 4);
		}
		if (lower.endsWith(".dia")) {
			return fileName.substring(0, fileName.length() - 4);
		}
		if (lower.endsWith(".mzml")) {
			return fileName.substring(0, fileName.length() - 5);
		}
		if (lower.endsWith(".d")) {
			return fileName.substring(0, fileName.length() - 2);
		}
		return fileName;
	}

	private File toDiaInput(File inputFile) {
		String lower = inputFile.getName().toLowerCase();
		if (lower.endsWith(".raw")) {
			File diaFile = new File(inputFile.getAbsolutePath() + ".dia");
			if (!diaFile.exists()) {
				DiannParameter conversionPar = new DiannParameter();
				conversionPar.setThreads(Math.max(1, threadCount / Math.max(1, maxTaskCount)));
				diaNNTask.addTask(conversionPar, new String[] { inputFile.getAbsolutePath() });
				diaNNTask.run(1);
			}
			return diaFile.exists() ? diaFile : null;
		}
		return inputFile;
	}

	private boolean runNoQuantSearch(String diaPath, String outputParquet, String[] libs) {
		DiaLibSearchPar par = new DiaLibSearchPar();
		par.setDiaFiles(new String[] { diaPath });
		par.setThreadCount(Math.max(1, threadCount / Math.max(1, maxTaskCount)));
		par.setMbr(false);
		par.setqValue(0.01);

		DiannParameter dpar = new DiannParameter();
		dpar.setThreads(par.getThreadCount());
		dpar.setMatrices(false);
		dpar.setRelaxed_prot_inf(true);
		dpar.setReanalyse(true);
		dpar.setExportDecoy(false);

		diaNNTask.addNoQuantLibTask(dpar, par, outputParquet, libs);
		diaNNTask.run(12);
		return new File(outputParquet).exists();
	}

	private File[] resolveSeedLibraries() {
		ArrayList<File> files = new ArrayList<>();
		for (String libName : libPathList) {
			File direct = new File(libName);
			if (direct.exists()) {
				files.add(direct);
				continue;
			}
			if (magDb != null) {
				File resolved = magDb.resolveSpectralLibraryFile(libName);
				if (resolved != null && resolved.exists()) {
					files.add(resolved);
					continue;
				}
			}
			File libFile = new File(new File(libFolderFile, libName), libName + ".parquet");
			if (libFile.exists()) {
				files.add(libFile);
			}
		}
		return files.toArray(new File[0]);
	}

	private HashSet<String> selectGenomesFromParquet(File parquetFile, File outDir, String kneeMode) throws IOException {
		HashSet<String> genomes = new HashSet<>();
		if (!parquetFile.exists() || magDb == null || magDb.getCurrentFile() == null) {
			return genomes;
		}

		File markerFasta = new File(magDb.getCurrentFile(), "Marker4.fasta");
		File genomesTaxaTsv = new File(magDb.getCurrentFile(), "genomes-taxa.tsv");
		if (!markerFasta.exists() || !genomesTaxaTsv.exists()) {
			return genomes;
		}

		DiaNNParquetReader parquetReader = new DiaNNParquetReader(parquetFile.getAbsolutePath(), diannPar.getCut());
		DiaNNPrecursor[] precursors = parquetReader.getDiaNNPrecursors(0.01);
		if (precursors == null || precursors.length == 0) {
			return genomes;
		}
		File predictedFolder = new File(magDb.getCurrentFile(), "predicted");
		List<DiaNNPrecursor> withProtein = MagPredictedPeptideProteinMapper.mapPrecursorsToProteins(precursors, predictedFolder);
		if (withProtein == null || withProtein.isEmpty()) {
			return genomes;
		}

		MarkerGenomeSelector.MarkerPeptideRatioResult ratio = MarkerGenomeSelector.computeMarkerPeptideRatio(withProtein, markerFasta);
		File effectiveMarker = (ratio.ratio >= MarkerGenomeSelector.MARKER_PEPTIDE_RATIO_THRESHOLD) ? markerFasta : null;

		File evidenceFile = new File(outDir, "genome_evidence.tsv");
		File rankedFile = new File(outDir, "genome_ranked.tsv");
		File selectedFile = new File(outDir, "genome_selected.tsv");

		MagGenomeSelectionService.Request selectionRequest = new MagGenomeSelectionService.Request();
		selectionRequest.magDb = magDb;
		selectionRequest.genomeSelectorMode = MagDbItem.MARKER4;
		selectionRequest.genomeSelectionStrategy = kneeMode;
		selectionRequest.diaPrecursors = withProtein.toArray(new DiaNNPrecursor[0]);
		selectionRequest.markerFasta = effectiveMarker;
		selectionRequest.disableMarkerFasta = effectiveMarker == null;
		selectionRequest.genomesTaxaTsv = genomesTaxaTsv;
		selectionRequest.evidenceFile = evidenceFile;
		selectionRequest.rankedFile = rankedFile;
		selectionRequest.selectedFile = selectedFile;
		selectionRequest.maximumSelected = 300;
		try {
			return MagGenomeSelectionService.selectGenomes(selectionRequest, new ProgressCallback.NoOp());
		} catch (Exception e) {
			throw new IOException("Genome selection failed for " + parquetFile.getAbsolutePath(), e);
		}
	}

	private boolean writeFastaFromGenomeSet(File outFasta, HashSet<String> genomeSet) {
		File dbFolder = magDb != null ? magDb.getOriginalDbFolder() : null;
		if (dbFolder == null || !dbFolder.exists() || genomeSet == null || genomeSet.isEmpty()) {
			return false;
		}
		try (PrintWriter writer = new PrintWriter(outFasta)) {
			for (String genome : genomeSet) {
				File genomeDbFile = new File(dbFolder, genome + ".faa");
				if (!genomeDbFile.exists()) {
					genomeDbFile = new File(dbFolder, genome + ".fasta");
				}
				if (!genomeDbFile.exists()) {
					continue;
				}
				for (String line : Files.readAllLines(genomeDbFile.toPath())) {
					writer.println(line);
				}
			}
		} catch (IOException e) {
			return false;
		}
		return outFasta.exists() && outFasta.length() > 0;
	}

	private boolean buildSpeclibFromFasta(File fastaFile, File speclibFile) {
		if (speclibFile.exists() && speclibFile.length() > 0) {
			return true;
		}
		File predictedLib = new File(speclibFile.getParentFile(), speclibFile.getName().replace(".speclib", ".predicted.speclib"));
		if (predictedLib.exists()) {
			return predictedLib.renameTo(speclibFile);
		}
		DiannParameter buildPar = new DiannParameter();
		buildPar.setThreads(Math.max(1, threadCount));
		buildPar.setExportDecoy(true);
		diaNNTask.addTask(buildPar, fastaFile.getAbsolutePath(), speclibFile.getAbsolutePath(), false);
		diaNNTask.run(12);
		if (!speclibFile.exists() && predictedLib.exists()) {
			predictedLib.renameTo(speclibFile);
		}
		return speclibFile.exists() && speclibFile.length() > 0;
	}

	private int getFinishedCount() {
		int count = 0;
		File[] destFiles = destinationFolder.listFiles();
		if (destFiles == null) {
			return count;
		}
		for (File dir : destFiles) {
			if (!dir.isDirectory()) {
				continue;
			}
			File finalFile = new File(dir, dir.getName() + ".parquet");
			if (finalFile.exists() && finalFile.length() > 0) {
				count++;
			}
		}
		return count;
	}

	private String waitForFileCompletionAndProcess(Path filePath) {
		try {
			long previousSize = -1;
			int stableCount = 0;
			while (true) {
				long currentSize = Files.isDirectory(filePath) ? directorySize(filePath) : Files.size(filePath);
				if (currentSize == previousSize) {
					stableCount++;
					if (Files.isDirectory(filePath)) {
						if (stableCount >= 3) {
							break;
						}
					} else {
						try (FileChannel fileChannel = FileChannel.open(filePath, StandardOpenOption.WRITE)) {
							FileLock lock = fileChannel.tryLock();
							if (lock != null) {
								lock.release();
								break;
							}
						} catch (Exception e) {
							logger.info(getTaskName() + ": waiting for copy complete {}", filePath.toAbsolutePath());
						}
					}
				} else {
					stableCount = 0;
				}
				previousSize = currentSize;
				Thread.sleep(1000);
			}

			return process(new File(filePath.toAbsolutePath().toString()));

		} catch (IOException | InterruptedException e) {
			e.printStackTrace();
			Thread.currentThread().interrupt();
		}

		return null;
	}

	private long directorySize(Path dir) {
		final long[] size = new long[] { 0L };
		try {
			Files.walk(dir).forEach(p -> {
				try {
					if (Files.isRegularFile(p)) {
						size[0] += Files.size(p);
					}
				} catch (IOException ignored) {
				}
			});
		} catch (IOException ignored) {
		}
		return size[0];
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
