package uk.ac.quadram.metapilot.task.dia;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.concurrent.TimeUnit;

import org.apache.logging.log4j.Logger;

import javafx.scene.control.ProgressBar;

import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNParquetReader;
import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNPrecursor;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.task.dia.par.MetaParameterDiaFocus;
import uk.ac.quadram.metapilot.task.mag.MarkerGenomeSelector;

/**
 * DIA "Deep comprehensive" workflow task for the <b>seed strategy</b> only: reuses the same first step
 * as {@link MetaDiaSeedIterativeTask} (user library → result/seed/seed.parquet under base), then
 * rank-only genome selection, training, prediction, and iterative identification rounds.
 * No focus/ folder; seed output lives at base/seed/ (same folder name as seed comprehensive).
 */
public class MetaDiaDeepComprehensiveTask extends MetaDiaSeedIterativeTask {

	private static final String TASK_NAME = "DIA Deep Comprehensive analysis";
	private static final String BASE_FOLDER = "deep_comprehensive";
	/** Folder name for seed output (same as parent). */
	private static final String SEED_FOLDER = "seed";
	private static final String SEED_PARQUET = "seed.parquet";
	private static final String SEED_PR_MATRIX = "seed.pr_matrix.tsv";
	private static final String SEED_PG_MATRIX = "seed.pg_matrix.tsv";

	/** Passed from iterativeSearch to iden. */
	private List<String> deepComprehensiveSelectedGenomes;

	public MetaDiaDeepComprehensiveTask(MetaParameterDiaFocus focusParams, MetaResourceParameters mrp,
			ProgressBar bar1, ProgressBar bar2) {
		super(focusParams, mrp, bar1, bar2);
	}

	@Override
	protected String getTaskName() {
		return TASK_NAME;
	}

	@Override
	protected Logger getLogger() {
		return logger;
	}

	@Override
	protected File getIterativeResultFolder() {
		File base = metaPar.getDbSearchResultFile();
		return base != null ? new File(base, SEED_FOLDER) : new File(SEED_FOLDER);
	}

	/** Deep Comprehensive uses all user-selected libraries for the first search; no fallback to direct-only. */
	@Override
	protected boolean classifyLibraries() {
		String[] libPath = resolveLibraryPaths();
		if (libPath == null || libPath.length == 0) {
			logger.error(getTaskName() + ": No libraries resolved (MetaParameterDia and MAG required)");
			return false;
		}
		initialSearchLibPaths = libPath;
		fallbackDirectSearch = false;
		addOutLibForInitialSearch = false;
		standaloneSearchLibPaths = null;
		logger.info(getTaskName() + ": using all {} user-selected libraries for first search", libPath.length);
		return true;
	}

	/** First step: parent's runFirstSearchUserLibrary() → base/seed/seed.parquet. Then map peptides to MAG proteins, rank/training/prediction. */
	@Override
	protected boolean iterativeSearch() {
		if (!runFirstSearchUserLibrary()) {
			return false;
		}
		File base = metaPar.getDbSearchResultFile();
		File seedFolder = getIterativeResultFolder();
		File parquetFile = new File(seedFolder, SEED_PARQUET);
		if (!parquetFile.exists()) {
			deepComprehensiveSelectedGenomes = null;
			return true;
		}
		String cut = diannPar != null && diannPar.getCut() != null ? diannPar.getCut() : "KR";
		DiaNNParquetReader reader = new DiaNNParquetReader(parquetFile.getAbsolutePath(), cut);
		DiaNNPrecursor[] precursors = reader.getDiaNNPrecursors(0.01);
		if (precursors == null || precursors.length == 0) {
			deepComprehensiveSelectedGenomes = null;
			return true;
		}
		if (magDb == null) {
			buildTrainingDatasetFromPrecursors(base, precursors, null);
			deepComprehensiveSelectedGenomes = null;
			return true;
		}

		// Map peptides to MAG proteins (genome_proteinId) via predicted TSVs. setProteinInfoFromPredictedTsVs
		// modifies matched precursors in place (setProGroups to MAG IDs); returns only those that had a match.
		File predictedFolder = new File(magDb.getCurrentFile(), "predicted");
		List<DiaNNPrecursor> precursorsWithPro = null;
		if (predictedFolder.isDirectory()) {
			precursorsWithPro = MarkerGenomeSelector.setProteinInfoFromPredictedTsVs(precursors, predictedFolder);
			if (precursorsWithPro != null && !precursorsWithPro.isEmpty()) {
				logger.info(getTaskName() + ": mapped {} / {} precursors to MAG proteins", precursorsWithPro.size(), precursors.length);
			}
		}
		if (precursorsWithPro == null || precursorsWithPro.isEmpty()) {
			logger.warn(getTaskName() + ": no precursors mapped to MAG proteins (check predicted TSVs); skipping genome selection");
			buildTrainingDatasetFromPrecursors(base, precursors, null);
			deepComprehensiveSelectedGenomes = null;
			return true;
		}

		File markerFasta = new File(magDb.getCurrentFile(), "Marker4.fasta");
		File genomesTaxaTsv = new File(magDb.getCurrentFile(), "genomes-taxa.tsv");
		if (!markerFasta.exists() || !genomesTaxaTsv.exists()) {
			buildTrainingDatasetFromPrecursors(base, precursors, null);
			deepComprehensiveSelectedGenomes = null;
			return true;
		}
		updateMessage(getTaskName() + ": rank-only genome selection, training, prediction...");
		// Genome selection outputs in seed folder (same as seed comprehensive / seed iterative)
		File evidenceFile = new File(seedFolder, "genome_evidence.tsv");
		File rankedFile = new File(seedFolder, "genome_ranked.tsv");
		File selectedFile = new File(seedFolder, "genome_selected.tsv");
		MarkerGenomeSelector selector = new MarkerGenomeSelector();
		selector.setKneeSelectionMode(MarkerGenomeSelector.KNEE_MODE_COMPREHENSIVE);
		DiaNNPrecursor[] withProArray = precursorsWithPro.toArray(new DiaNNPrecursor[0]);
		selector.selectMarkerGenomesWithKnee(withProArray, markerFasta, genomesTaxaTsv,
				evidenceFile.getAbsolutePath(), rankedFile.getAbsolutePath(), selectedFile.getAbsolutePath(),
				300, true);
		List<String> selectedGenomes = selectGenomesUntilMarginalOne(rankedFile, withProArray);
		int nGenomes = selectedGenomes != null ? selectedGenomes.size() : 0;
		logger.info(getTaskName() + ": rank-only selection (marginal-one) produced {} genomes", nGenomes);
		// Use only mapped precursors so training/predict use MAG IDs for genomes-taxa.tsv and {genome}_eggNOG.tsv
		// (n_genomes, tax_sum_*, fun_sum_*). Unmapped precursors would have wrong keys and yield zeros.
		buildTrainingDatasetFromPrecursors(base, withProArray, selectedGenomes);
		runDeepComprehensiveTraining(base);
		buildPredictDatasetsForRuns(base, withProArray, selectedGenomes);
		runDeepComprehensivePrediction(base);
		deepComprehensiveSelectedGenomes = selectedGenomes;
		return true;
	}

	@Override
	protected boolean exportReport() {
		File base = metaPar.getDbSearchResultFile();
		if (base == null) {
			return true;
		}
		logger.info(getTaskName() + ": base folder {}", base.getAbsolutePath());
		logger.info(getTaskName() + ": seed folder {}", getIterativeResultFolder().getAbsolutePath());
		logger.info(getTaskName() + ": training dataset {}", DeepComprehensivePaths.getTrainingDatasetFile(base).getAbsolutePath());
		logger.info(getTaskName() + ": model dir {}", DeepComprehensivePaths.getModelDir(base).getAbsolutePath());
		return true;
	}

	@Override
	protected Boolean call() throws Exception {
		String resultPath = metaPar.getResult();
		if (resultPath == null || resultPath.trim().isEmpty()) {
			logger.error(getTaskName() + ": result path is null or empty");
			return false;
		}
		File base = new File(resultPath, BASE_FOLDER);
		if (!base.exists()) {
			base.mkdirs();
		}
		DeepComprehensivePaths.ensureTrainingFolders(base);
		if (quanExpNames != null) {
			for (int i = 0; i < quanExpNames.length; i++) {
				String runId = (quanExpNames != null && i < quanExpNames.length && quanExpNames[i] != null)
						? quanExpNames[i]
						: (rawFiles != null && i < rawFiles.length) ? DeepComprehensivePaths.runIdFromRawFilePath(rawFiles[i]) : "run";
				DeepComprehensivePaths.ensurePerRunFolder(base, runId);
			}
		}
		metaPar.setDbSearchResultFile(base);
		this.resultFolderFile = base;
		logger.info(getTaskName() + ": using base folder {}", base.getAbsolutePath());

		updateMessage(getTaskName());
		updateProgress(0, 100);
		logger.info(getTaskName() + ": started");
		System.out.println(new java.text.SimpleDateFormat("yyyy-MM-dd HH:mm:ss").format(new java.util.Date()) + "\t" + getTaskName() + ": started");

		updateProgress(5, 100);
		updateSecondaryProgress(ProgressBar.INDETERMINATE_PROGRESS);
		if (!copyReportTemplatesInBackground()) {
			logger.error(getTaskName() + ": Failed to copy report templates");
			return false;
		}

		logger.info("Step 1: Starting spectra format conversion...");
		updateMessage(getTaskName() + ": Converting spectra format...");
		updateProgress(10, 100);
		if (!spectraConvert()) {
			logger.error(getTaskName() + ": error in spectra format conversion");
			updateMessage(getTaskName() + ": Failed - Spectra conversion error");
			updateSecondaryProgress(0.0);
			return false;
		}

		logger.info("Step 2: Classifying libraries...");
		updateProgress(18, 100);
		updateMessage(getTaskName() + ": Checking library types...");
		if (!classifyLibraries()) {
			logger.error(getTaskName() + ": library classification failed");
			updateMessage(getTaskName() + ": Failed - Library classification error");
			updateSecondaryProgress(0.0);
			return false;
		}

		logger.info("Step 3: Validating user library...");
		updateProgress(20, 100);
		updateMessage(getTaskName() + ": Validating library...");
		if (!createSSLibrary()) {
			logger.error(getTaskName() + ": error in library validation");
			updateMessage(getTaskName() + ": Failed - Library validation error");
			updateSecondaryProgress(0.0);
			return false;
		}

		logger.info("Step 4: Iterative search (first search + rank/training/prediction)...");
		updateProgress(25, 100);
		updateMessage(getTaskName() + ": Performing iterative search...");
		updateSecondaryProgress(ProgressBar.INDETERMINATE_PROGRESS);
		if (!iterativeSearch()) {
			logger.error(getTaskName() + ": error in iterative search");
			updateMessage(getTaskName() + ": Failed - Iterative search error");
			updateSecondaryProgress(0.0);
			return false;
		}
		updateProgress(50, 100);

		logger.info("Step 5: Identification (iterative rounds)...");
		updateMessage(getTaskName() + ": Running identification...");
		updateProgress(55, 100);
		updateSecondaryProgress(ProgressBar.INDETERMINATE_PROGRESS);
		if (!iden()) {
			logger.error(getTaskName() + ": error in identification");
			updateMessage(getTaskName() + ": Failed - Identification error");
			updateSecondaryProgress(0.0);
			return false;
		}

		logger.info("Step 6: Export summary...");
		updateProgress(85, 100);
		updateMessage(getTaskName() + ": Generating reports...");
		updateSecondaryProgress(ProgressBar.INDETERMINATE_PROGRESS);
		if (!exportReport()) {
			logger.error(getTaskName() + ": error in report generation");
			updateMessage(getTaskName() + ": Failed - Report generation error");
			updateSecondaryProgress(0.0);
			return false;
		}

		updateProgress(100, 100);
		updateMessage(getTaskName() + ": Completed successfully!");
		updateSecondaryProgress(0.0);
		logger.info(getTaskName() + ": finished - all workflow steps completed successfully");
		System.out.println(new java.text.SimpleDateFormat("yyyy-MM-dd HH:mm:ss").format(new java.util.Date()) + "\t" + getTaskName() + ": finished successfully");
		return true;
	}

	@Override
	protected boolean iden() {
		File base = metaPar.getDbSearchResultFile();
		File seedFolder = getIterativeResultFolder();
		combinedFile = seedFolder;
		this.resultFolderFile = base;
		this.quan_pep_file = new File(seedFolder, SEED_PR_MATRIX);
		this.quan_pro_file = new File(seedFolder, SEED_PG_MATRIX);
		List<String> selected = deepComprehensiveSelectedGenomes != null ? deepComprehensiveSelectedGenomes : new ArrayList<>();
		runIterativeSearchRounds(base, selected, SEED_PARQUET);
		logger.info(getTaskName() + ": Identification (iterative rounds) completed");
		return true;
	}

	/**
	 * Read genome_ranked.tsv and select genomes in order until a row has Increment == 1.
	 * Uses the existing "Increment" column from the ranked file (no precursor-based recomputation).
	 * Ranked file format: header "Rank\tGenome\t...\tIncrement\t...", genome ID in col 1.
	 */
	private List<String> selectGenomesUntilMarginalOne(File rankedTsv, DiaNNPrecursor[] precursors) {
		List<String> result = new ArrayList<>();
		if (rankedTsv == null || !rankedTsv.exists()) {
			return result;
		}
		try (BufferedReader r = new BufferedReader(new FileReader(rankedTsv))) {
			String headerLine = r.readLine();
			if (headerLine == null) return result;
			String[] headers = headerLine.split("\t");
			int incrementCol = -1;
			for (int i = 0; i < headers.length; i++) {
				if ("Increment".equals(headers[i].trim())) {
					incrementCol = i;
					break;
				}
			}
			if (incrementCol < 0) {
				logger.warn(getTaskName() + ": no 'Increment' column in " + rankedTsv.getAbsolutePath());
				return result;
			}
			String line;
			while ((line = r.readLine()) != null) {
				line = line.trim();
				if (line.isEmpty() || line.startsWith("#")) continue;
				String[] parts = line.split("\t");
				if (parts.length <= 1 || parts[1].trim().isEmpty()) continue;
				String genome = parts[1].trim();
				result.add(genome);
				if (parts.length > incrementCol) {
					try {
						int inc = Integer.parseInt(parts[incrementCol].trim());
						if (inc == 1) {
							logger.info(getTaskName() + ": stopping at genome {} (Increment == 1)", genome);
							break;
						}
					} catch (NumberFormatException ignored) { }
				}
			}
		} catch (IOException e) {
			logger.warn(getTaskName() + ": could not read ranked file " + rankedTsv.getAbsolutePath(), e);
		}
		return result;
	}

	private void buildTrainingDatasetFromPrecursors(File base, DiaNNPrecursor[] precursors, List<String> selectedGenomes) {
		File outTsv = DeepComprehensivePaths.getTrainingDatasetFile(base);
		if (magDb == null) {
			logger.warn(getTaskName() + ": cannot build training dataset without MAG DB");
			return;
		}
		try {
			DeepComprehensiveTrainingDatasetBuilder builder = new DeepComprehensiveTrainingDatasetBuilder();
			int rows = builder.build(precursors, magDb, outTsv);
			logger.info(getTaskName() + ": wrote training dataset to {} ({} rows)", outTsv.getAbsolutePath(), rows);
		} catch (IOException e) {
			logger.error(getTaskName() + ": failed to write training dataset to " + outTsv.getAbsolutePath(), e);
		}
	}

	/**
	 * For each raw file (run), filter precursors by runString matching the run, then build
	 * predict_dataset.tsv in that run's folder using run-specific P_0 and candidate peptides from selected genomes.
	 */
	private void buildPredictDatasetsForRuns(File base, DiaNNPrecursor[] precursors, List<String> selectedGenomes) {
		if (magDb == null || precursors == null || precursors.length == 0) {
			return;
		}
		HashSet<String> genomeSet = new HashSet<>();
		if (selectedGenomes != null) {
			genomeSet.addAll(selectedGenomes);
		}
		if (genomeSet.isEmpty()) {
			logger.warn(getTaskName() + ": no selected genomes, skipping predict dataset build");
			return;
		}
		int numRuns = quanExpNames != null ? quanExpNames.length : (rawFiles != null ? rawFiles.length : 0);
		if (numRuns == 0) {
			logger.warn(getTaskName() + ": no runs (quanExpNames/rawFiles), skipping predict dataset build");
			return;
		}
		DeepComprehensivePredictDatasetBuilder builder = new DeepComprehensivePredictDatasetBuilder();
		List<Object[]> runPeptideCounts = new ArrayList<>();
		for (int i = 0; i < numRuns; i++) {
			String runId = (quanExpNames != null && i < quanExpNames.length && quanExpNames[i] != null)
					? quanExpNames[i]
					: (rawFiles != null && i < rawFiles.length) ? DeepComprehensivePaths.runIdFromRawFilePath(rawFiles[i]) : "run_" + i;
			String rawFileName = (rawFiles != null && i < rawFiles.length) ? new File(rawFiles[i]).getName() : null;
			List<DiaNNPrecursor> runList = new ArrayList<>();
			for (DiaNNPrecursor p : precursors) {
				if (p == null) continue;
				String rs = p.getRunString();
				if (rs == null) continue;
				if (rs.equals(runId)) { runList.add(p); continue; }
				if (rawFileName != null && rs.equals(rawFileName)) { runList.add(p); continue; }
				if (rs.endsWith(runId) || (rawFileName != null && rs.endsWith(rawFileName))) { runList.add(p); }
			}
			DiaNNPrecursor[] runPrecursors = runList.toArray(new DiaNNPrecursor[0]);
			File outTsv = DeepComprehensivePaths.getPredictDatasetFile(base, runId);
			if (outTsv.isFile() && outTsv.length() > 0) {
				int existingRows = countTsvDataRows(outTsv);
				logger.info(getTaskName() + ": predict_dataset.tsv already exists for run {} ({} peptides), skipping build", runId, existingRows);
				runPeptideCounts.add(new Object[] { runId, Integer.valueOf(existingRows) });
				continue;
			}
			try {
				int rows = builder.build(runPrecursors, genomeSet, magDb, outTsv);
				logger.info(getTaskName() + ": wrote predict dataset for run {} to {} ({} peptides)", runId, outTsv.getAbsolutePath(), rows);
				runPeptideCounts.add(new Object[] { runId, Integer.valueOf(rows) });
			} catch (IOException e) {
				logger.error(getTaskName() + ": failed to write predict dataset for run " + runId + " to " + outTsv.getAbsolutePath(), e);
			}
		}
		long totalPeptides = 0;
		for (Object[] pair : runPeptideCounts) {
			String runId = (String) pair[0];
			int count = ((Integer) pair[1]).intValue();
			totalPeptides += count;
			logger.info(getTaskName() + ": predict dataset run {}: {} peptides", runId, count);
		}
		logger.info(getTaskName() + ": predict dataset total peptides across runs: {}", totalPeptides);
	}

	/** Count data rows in a TSV (line count minus header). Returns 0 if file missing or empty. */
	private int countTsvDataRows(File tsvFile) {
		if (tsvFile == null || !tsvFile.isFile()) return 0;
		int count = 0;
		try (BufferedReader r = new BufferedReader(new FileReader(tsvFile))) {
			if (r.readLine() == null) return 0;
			while (r.readLine() != null) count++;
		} catch (IOException e) {
			return 0;
		}
		return count;
	}

	private static final String RESOURCE_TRAIN_SCRIPT = "/python_scripts/deep_comprehensive_train.py";
	private static final String RESOURCE_PREDICT_SCRIPT = "/python_scripts/deep_comprehensive_predict.py";

	/**
	 * Run the Deep Comprehensive training script to produce model.pt and scaler.npz in training/model/.
	 * Skips if Python is not configured or training_seed.tsv is missing/small.
	 */
	private void runDeepComprehensiveTraining(File base) {
		String python = mrp != null ? mrp.getToolPath(MetaResourceParameters.ToolType.PYTHON) : null;
		if (python == null || python.trim().isEmpty()) {
			logger.warn(getTaskName() + ": Python not configured, skipping model training");
			return;
		}
		File trainingDir = DeepComprehensivePaths.getTrainingDir(base);
		File trainingTsv = DeepComprehensivePaths.getTrainingDatasetFile(base);
		File modelDir = DeepComprehensivePaths.getModelDir(base);
		if (!trainingTsv.isFile() || trainingTsv.length() == 0) {
			logger.warn(getTaskName() + ": no training dataset at {}, skipping model training", trainingTsv.getAbsolutePath());
			return;
		}
		if (!extractScript(RESOURCE_TRAIN_SCRIPT, new File(trainingDir, "deep_comprehensive_train.py"))
				|| !extractScript(RESOURCE_PREDICT_SCRIPT, new File(trainingDir, "deep_comprehensive_predict.py"))) {
			logger.warn(getTaskName() + ": could not extract Python scripts, skipping model training");
			return;
		}
		List<String> cmd = new ArrayList<>();
		cmd.add(python);
		cmd.add(new File(trainingDir, "deep_comprehensive_train.py").getAbsolutePath());
		cmd.add("--input_tsv");
		cmd.add(trainingTsv.getAbsolutePath());
		cmd.add("--output_dir");
		cmd.add(modelDir.getAbsolutePath());
		// Prefer GPU if available (script falls back to CPU when cuda not present)
		cmd.add("--device");
		cmd.add("cuda");
		boolean ok = runPythonProcess(cmd, trainingDir, "training");
		if (ok) {
			logger.info(getTaskName() + ": model training completed");
		} else {
			logger.warn(getTaskName() + ": model training failed or timed out");
		}
	}

	/**
	 * Run the Deep Comprehensive prediction script for each run that has predict_dataset.tsv.
	 * Writes predictions.tsv per run. Skips if Python not configured or model not found.
	 */
	private void runDeepComprehensivePrediction(File base) {
		String python = mrp != null ? mrp.getToolPath(MetaResourceParameters.ToolType.PYTHON) : null;
		if (python == null || python.trim().isEmpty()) {
			logger.warn(getTaskName() + ": Python not configured, skipping prediction");
			return;
		}
		File modelDir = DeepComprehensivePaths.getModelDir(base);
		File modelPt = new File(modelDir, "model.pt");
		File scalerNpz = new File(modelDir, "scaler.npz");
		if (!modelPt.isFile() || !scalerNpz.isFile()) {
			logger.warn(getTaskName() + ": model.pt or scaler.npz not found in {}, skipping prediction", modelDir.getAbsolutePath());
			return;
		}
		File trainingDir = DeepComprehensivePaths.getTrainingDir(base);
		int numRuns = quanExpNames != null ? quanExpNames.length : (rawFiles != null ? rawFiles.length : 0);
		if (numRuns == 0) return;
		for (int i = 0; i < numRuns; i++) {
			String runId = (quanExpNames != null && i < quanExpNames.length && quanExpNames[i] != null)
					? quanExpNames[i]
					: (rawFiles != null && i < rawFiles.length) ? DeepComprehensivePaths.runIdFromRawFilePath(rawFiles[i]) : "run_" + i;
			File predictTsv = DeepComprehensivePaths.getPredictDatasetFile(base, runId);
			File outTsv = DeepComprehensivePaths.getPredictionsFile(base, runId);
			if (!predictTsv.isFile() || predictTsv.length() == 0) continue;
			List<String> cmd = new ArrayList<>();
			cmd.add(python);
			cmd.add(new File(trainingDir, "deep_comprehensive_predict.py").getAbsolutePath());
			cmd.add("--model_dir");
			cmd.add(modelDir.getAbsolutePath());
			cmd.add("--predict_tsv");
			cmd.add(predictTsv.getAbsolutePath());
			cmd.add("--output_tsv");
			cmd.add(outTsv.getAbsolutePath());
			cmd.add("--device");
			cmd.add("cuda");
			boolean ok = runPythonProcess(cmd, trainingDir, "predict-" + runId);
			if (ok) {
				logger.info(getTaskName() + ": predictions written for run {} to {}", runId, outTsv.getAbsolutePath());
			} else {
				logger.warn(getTaskName() + ": prediction failed for run {}", runId);
			}
		}
	}

	private boolean extractScript(String resourcePath, File targetFile) {
		try (InputStream in = MetaDiaDeepComprehensiveTask.class.getResourceAsStream(resourcePath)) {
			if (in == null) {
				logger.warn(getTaskName() + ": resource not found: {}", resourcePath);
				return false;
			}
			targetFile.getParentFile().mkdirs();
			try (BufferedReader reader = new BufferedReader(new InputStreamReader(in, StandardCharsets.UTF_8));
					PrintWriter writer = new PrintWriter(targetFile, StandardCharsets.UTF_8.name())) {
				String line;
				while ((line = reader.readLine()) != null) {
					writer.println(line);
				}
			}
			return true;
		} catch (IOException e) {
			logger.warn(getTaskName() + ": failed to extract " + resourcePath, e);
			return false;
		}
	}

	private boolean runPythonProcess(List<String> command, File workingDir, String label) {
		try {
			ProcessBuilder pb = new ProcessBuilder(command);
			pb.directory(workingDir);
			pb.redirectErrorStream(true);
			Process p = pb.start();
			try (BufferedReader in = new BufferedReader(new InputStreamReader(p.getInputStream(), StandardCharsets.UTF_8))) {
				String line;
				while ((line = in.readLine()) != null) {
					logger.info("{} [{}]: {}", getTaskName(), label, line);
				}
			}
			boolean finished = p.waitFor(30, TimeUnit.MINUTES);
			if (!finished) {
				p.destroyForcibly();
				logger.warn(getTaskName() + ": Python process timed out ({})", label);
				return false;
			}
			return p.exitValue() == 0;
		} catch (Exception e) {
			logger.warn(getTaskName() + ": Python process failed (" + label + ")", e);
			return false;
		}
	}

	/** Max iterative search rounds (build library from predictions + identified → search). */
	private static final int MAX_ITERATIVE_ROUNDS = 3;
	/** Top predicted peptides to add per round (in addition to all previously identified). */
	private static final int TOP_PREDICTED_PER_ROUND = 5000;

	/**
	 * Iterative rounds: build library from all peptides identified in the previous round plus top
	 * predicted peptides, search with that library, stop when no new peptides or max rounds.
	 * @param initialParquetName parquet file name in seed folder (e.g. "seed.parquet")
	 */
	private void runIterativeSearchRounds(File base, List<String> selectedGenomes, String initialParquetName) {
		if (magDb == null || selectedGenomes == null || selectedGenomes.isEmpty()) {
			return;
		}
		File seedFolder = getIterativeResultFolder();
		File prevParquet = new File(seedFolder, initialParquetName);
		if (!prevParquet.isFile()) {
			logger.warn(getTaskName() + ": no seed parquet at {}, skipping iterative rounds", prevParquet.getAbsolutePath());
			return;
		}
		HashSet<String> genomeSet = new HashSet<>(selectedGenomes);
		java.util.Map<String, Set<String>> candidatePeptideToProteins = loadCandidatePeptideToProteins(genomeSet);
		String cut = diannPar != null && diannPar.getCut() != null ? diannPar.getCut() : "KR";
		for (int round = 1; round <= MAX_ITERATIVE_ROUNDS; round++) {
			updateMessage(getTaskName() + ": iterative round " + round + " (build library + search)");
			java.util.Map<String, Set<String>> identifiedPepToProteins = getIdentifiedPeptidesAndProteins(prevParquet, cut, 0.01);
			Set<String> identifiedPeptides = new HashSet<>(identifiedPepToProteins.keySet());
			Set<String> toAdd = mergePredictionsAndTakeTop(base, identifiedPeptides, TOP_PREDICTED_PER_ROUND);
			Set<String> allProteins = new HashSet<>();
			for (Set<String> prots : identifiedPepToProteins.values()) {
				allProteins.addAll(prots);
			}
			for (String pep : toAdd) {
				Set<String> prots = candidatePeptideToProteins.get(pep);
				if (prots != null) allProteins.addAll(prots);
			}
			if (allProteins.isEmpty()) {
				logger.warn(getTaskName() + ": round {} no proteins to build library", round);
				break;
			}
			// Only use protein IDs that exist in the MAG DB (genome.faa/fasta). Drop IDs from user/library search (e.g. pro17511)
			Set<String> magProteins = DeepComprehensiveLibraryBuilder.filterProteinIdsToMagGenomes(allProteins, magDb);
			if (magProteins.size() < allProteins.size()) {
				logger.info(getTaskName() + ": round {} using {} protein IDs from MAG DB (dropped {} from other libraries)",
						round, magProteins.size(), allProteins.size() - magProteins.size());
			}
			if (magProteins.isEmpty()) {
				logger.warn(getTaskName() + ": round {} no MAG proteins to build library (all IDs were from other sources)", round);
				break;
			}
			File roundDir = new File(base, "round" + round);
			roundDir.mkdirs();
			File roundFasta = new File(roundDir, "combined.fasta");
			File roundSpeclib = new File(roundDir, "combined.speclib");
			try {
				DeepComprehensiveLibraryBuilder.writeFastaFromProteinIds(magProteins, magDb, roundFasta);
			} catch (IOException e) {
				logger.error(getTaskName() + ": failed to write round FASTA", e);
				break;
			}
			if (!roundFasta.isFile() || roundFasta.length() == 0) {
				logger.warn(getTaskName() + ": round {} FASTA empty", round);
				break;
			}
			diannPar.setThreads(threadCount);
			diaNNTask.addTask(diannPar, roundFasta.getAbsolutePath(), roundSpeclib.getAbsolutePath(), false);
			diaNNTask.run(12);
			File altLib = new File(roundDir, "combined.predicted.speclib");
			if (!roundSpeclib.isFile() && altLib.isFile()) {
				altLib.renameTo(roundSpeclib);
			}
			if (!roundSpeclib.isFile()) {
				logger.warn(getTaskName() + ": round {} speclib not created", round);
				break;
			}
			File roundParquet = new File(roundDir, "combined.parquet");
			diaNNTask.addTask(diannPar, diaLibSearchPar, roundParquet.getAbsolutePath(),
					new String[] { roundSpeclib.getAbsolutePath() });
			diaNNTask.run(diaFiles.length * 12);
			if (!roundParquet.isFile()) {
				logger.warn(getTaskName() + ": round {} search did not produce parquet", round);
				break;
			}
			DiaNNParquetReader r = new DiaNNParquetReader(roundParquet.getAbsolutePath(), cut);
			DiaNNPrecursor[] roundPrecursors = r.getDiaNNPrecursors(0.01);
			int prevCount = identifiedPeptides.size();
			Set<String> roundPeptides = new HashSet<>();
			if (roundPrecursors != null) {
				for (DiaNNPrecursor p : roundPrecursors) {
					if (p != null && p.getSeqString() != null) roundPeptides.add(p.getSeqString());
				}
			}
			int newCount = roundPeptides.size();
			logger.info(getTaskName() + ": round {} identified {} peptides (prev {})", round, newCount, prevCount);
			if (newCount <= prevCount && round > 1) {
				logger.info(getTaskName() + ": no new peptides, stopping iterative rounds");
				break;
			}
			prevParquet = roundParquet;
		}
	}

	private java.util.Map<String, Set<String>> getIdentifiedPeptidesAndProteins(File parquetFile, String cut, double qMax) {
		java.util.Map<String, Set<String>> out = new java.util.HashMap<>();
		try {
			DiaNNParquetReader reader = new DiaNNParquetReader(parquetFile.getAbsolutePath(), cut);
			DiaNNPrecursor[] precursors = reader.getDiaNNPrecursors(qMax);
			if (precursors == null) return out;
			for (DiaNNPrecursor p : precursors) {
				if (p == null || p.getIsDecoy()) continue;
				String pep = p.getSeqString();
				if (pep == null || pep.isEmpty()) continue;
				Set<String> prots = out.computeIfAbsent(pep, k -> new HashSet<>());
				if (p.getProGroups() != null) {
					for (String pro : p.getProGroups()) {
						if (pro != null && !pro.isEmpty()) prots.add(pro);
					}
				}
			}
		} catch (Exception e) {
			logger.warn(getTaskName() + ": failed to read parquet " + parquetFile.getAbsolutePath(), e);
		}
		return out;
	}

	private java.util.Map<String, Set<String>> loadCandidatePeptideToProteins(HashSet<String> selectedGenomes) {
		java.util.Map<String, Set<String>> out = new java.util.HashMap<>();
		File predFolder = magDb != null ? magDb.getPredictedPepFile() : null;
		if (predFolder == null || !predFolder.isDirectory()) return out;
		try {
			java.util.Map<String, DeepComprehensivePredictDatasetBuilder.CandidatePeptide> candidates =
					DeepComprehensivePredictDatasetBuilder.loadCandidatePeptidesFromPredicted(predFolder, selectedGenomes);
			for (java.util.Map.Entry<String, DeepComprehensivePredictDatasetBuilder.CandidatePeptide> e : candidates.entrySet()) {
				out.put(e.getKey(), new HashSet<>(e.getValue().proteins));
			}
		} catch (IOException e) {
			logger.warn(getTaskName() + ": failed to load candidate peptide→proteins", e);
		}
		return out;
	}

	/** Merge predictions from all runs, sort by prob_detect * (1 + pred_intensity), return top K peptides not in exclude. */
	private Set<String> mergePredictionsAndTakeTop(File base, Set<String> excludePeptides, int topK) {
		List<Object[]> withPep = new ArrayList<>();
		int numRuns = quanExpNames != null ? quanExpNames.length : (rawFiles != null ? rawFiles.length : 0);
		for (int i = 0; i < numRuns; i++) {
			String runId = (quanExpNames != null && i < quanExpNames.length && quanExpNames[i] != null)
					? quanExpNames[i]
					: (rawFiles != null && i < rawFiles.length) ? DeepComprehensivePaths.runIdFromRawFilePath(rawFiles[i]) : "run_" + i;
			File predFile = DeepComprehensivePaths.getPredictionsFile(base, runId);
			if (!predFile.isFile()) continue;
			try (BufferedReader r = new BufferedReader(new FileReader(predFile))) {
				String line = r.readLine();
				if (line == null || !line.contains("peptide")) continue;
				while ((line = r.readLine()) != null) {
					line = line.trim();
					if (line.isEmpty()) continue;
					String[] cs = line.split("\t");
					if (cs.length < 3) continue;
					String pep = cs[0].trim();
					if (excludePeptides != null && excludePeptides.contains(pep)) continue;
					double prob, intensity;
					try {
						prob = Double.parseDouble(cs[1].trim());
						intensity = Double.parseDouble(cs[2].trim());
					} catch (NumberFormatException ignored) { continue; }
					double score = prob * (1.0 + Math.max(0, intensity));
					withPep.add(new Object[] { Double.valueOf(score), pep });
				}
			} catch (IOException e) {
				logger.warn(getTaskName() + ": failed to read predictions " + predFile.getAbsolutePath(), e);
			}
		}
		withPep.sort((a, b) -> Double.compare(((Double) b[0]).doubleValue(), ((Double) a[0]).doubleValue()));
		Set<String> result = new HashSet<>();
		for (int i = 0; i < Math.min(topK, withPep.size()); i++) {
			result.add((String) withPep.get(i)[1]);
		}
		return result;
	}
}
