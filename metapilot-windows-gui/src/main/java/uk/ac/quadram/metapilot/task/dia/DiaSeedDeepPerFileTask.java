package uk.ac.quadram.metapilot.task.dia;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.logging.log4j.Logger;

import javafx.scene.control.ProgressBar;

import uk.ac.quadram.metapilot.core.dia.DiaDeepProteinMappingService;
import uk.ac.quadram.metapilot.core.dia.MagPredictedPeptideProteinMapper;
import uk.ac.quadram.metapilot.core.mag.MagGenomeSelectionService;
import uk.ac.quadram.metapilot.core.mag.Marker4ResourcePreparer;
import uk.ac.quadram.metapilot.core.runner.DiaNNNoQuantLibrarySearchService;
import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNParquetReader;
import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNPrecursor;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.task.dia.par.MetaParameterDiaMarker;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;

/**
 * DIA "Seed deep" workflow task for the <b>seed strategy</b>.
 * <p>
 * High-level flow:
 * <ol>
 * <li>Reuse seed first search from {@link DiaSeedCoreTask} to generate
 *     {@code seed.parquet} (or explicit per-file parquet outputs).</li>
 * <li>Map identified peptides to MAG proteins via MAG {@code predicted/} TSVs.</li>
 * <li>Run marker-style genome ranking to produce {@code genome_ranked.tsv}.</li>
 * <li>Run iterative rounds (shared round logic) to expand evidence.</li>
 * <li>Build combined FASTA/library and run final quantifying search.</li>
 * </ol>
 *
 * DIA-NN usage map in this class:
 * <ul>
 * <li>First search is in parent class ({@code runFirstSearchUserLibrary()}), not reimplemented here.</li>
 * <li>{@code addTask(... combined.fasta -> combined.speclib ...)}: build predicted library from FASTA.</li>
 * <li>{@code addNoQuantLibTask(... roundN.parquet ...)}: round-level non-quant identification.</li>
 * <li>{@code addTask(... combined.parquet ...)}: final quantifying search over DIA files.</li>
 * </ul>
 */
public class DiaSeedDeepPerFileTask extends DiaSeedCoreTask {

	private static final String TASK_NAME = "DIA Seed Deep analysis";
	private static final String BASE_FOLDER = "seed_deep";
	/** Folder name for seed output (same as parent). */
	private static final String SEED_FOLDER = "seed";
	private static final String SEED_PARQUET = "seed.parquet";

	/** Number of top genomes from genome_ranked.tsv per round. */
	private static final int DEFAULT_GENOMES_PER_ROUND = 5;
	/** Round stop trigger for peptide-based deep mode: stop after current genome if cumulative stripped peptides exceed this value. */
	private static final int PEPTIDE_ROUND_TRIGGER = 1_000_000;
	/** Maximum number of iterative genome rounds. */
	private static final int MAX_ITERATIVE_ROUNDS = 100;

	public DiaSeedDeepPerFileTask(MetaParameterDiaMarker markerParams, MetaResourceParameters mrp,
			ProgressBar bar1, ProgressBar bar2) {
		this(withMode(markerParams, DiaRefinementMode.DEEP_PER_FILE), mrp, bar1, bar2, true);
	}

	protected DiaSeedDeepPerFileTask(MetaParameterDiaMarker markerParams, MetaResourceParameters mrp, ProgressBar bar1,
			ProgressBar bar2, DiaRefinementMode refinementMode) {
		this(withMode(markerParams, refinementMode), mrp, bar1, bar2, true);
	}

	private DiaSeedDeepPerFileTask(MetaParameterDiaMarker markerParams, MetaResourceParameters mrp, ProgressBar bar1,
			ProgressBar bar2, boolean modeAlreadyApplied) {
		super(markerParams, mrp, bar1, bar2, null);
	}

	private static MetaParameterDiaMarker withMode(MetaParameterDiaMarker markerParams, DiaRefinementMode refinementMode) {
		if (markerParams != null && refinementMode != null) {
			markerParams.setMarkerRefinementMode(refinementMode);
		}
		return markerParams;
	}

	private int getGenomesPerRound() {
		if (metaPar instanceof MetaParameterDiaMarker) {
			int v = ((MetaParameterDiaMarker) metaPar).getGenomesPerRound();
			return v > 0 ? v : DEFAULT_GENOMES_PER_ROUND;
		}
		return DEFAULT_GENOMES_PER_ROUND;
	}

	@Override
	protected void initResultFolderIfNeeded() {
		// Seed deep manages its own folder structure (seed_deep) in call();
		// do not change dbSearchResultFile here.
		if (metaPar != null && metaPar.getDbSearchResultFile() != null) {
			this.resultFolderFile = metaPar.getDbSearchResultFile();
		}
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

	/**
	 * Seed-specific genome ranking stage.
	 * First-search itself is delegated to parent ({@link #runFirstSearchUserLibrary()}).
	 * The key seed-specific difference is peptide→MAG protein mapping using MAG {@code predicted/}
	 * before running marker-style genome ranking.
	 */
	@Override
	protected boolean iterativeSearch() {
		if (!runFirstSearchUserLibrary()) {
			return false;
		}
		File seedFolder = getIterativeResultFolder();
		File parquetFile = new File(seedFolder, SEED_PARQUET);
		String cut = diannPar != null && diannPar.getCut() != null ? diannPar.getCut() : "KR";
		DiaNNPrecursor[] precursors;
		if (parquetFile.exists()) {
			DiaNNParquetReader reader = new DiaNNParquetReader(parquetFile.getAbsolutePath(), cut);
			precursors = reader.getDiaNNPrecursors(0.01);
		} else {
			logger.warn(getTaskName()
					+ ": seed.parquet not created for deep seed first step; aggregating precursors from per_file_*.parquet");
			precursors = readPrecursorsFromPerFileParquets(seedFolder, 0.01);
		}
		if (precursors == null || precursors.length == 0) {
			logger.warn(getTaskName() + ": No precursors found in seed parquet; cannot rank genomes");
			return false;
		}
		if (magDb == null) {
			logger.error(getTaskName() + ": MAG database is required for deep seed genome selection");
			return false;
		}

		// Map peptides to MAG proteins (genome_proteinId) via predicted TSVs.
		File predictedFolder = new File(magDb.getCurrentFile(), "predicted");
		if (!predictedFolder.isDirectory()) {
			logger.error(getTaskName() + ": MAG predicted folder not found: " + predictedFolder.getAbsolutePath());
			return false;
		}
		List<DiaNNPrecursor> precursorsWithPro = MagPredictedPeptideProteinMapper.mapPrecursorsToProteins(precursors,
				predictedFolder);
		if (precursorsWithPro != null && !precursorsWithPro.isEmpty()) {
			logger.info(getTaskName() + ": mapped {} / {} precursors to MAG proteins", precursorsWithPro.size(),
					precursors.length);
		}
		if (precursorsWithPro == null || precursorsWithPro.isEmpty()) {
			logger.error(getTaskName() + ": no precursors mapped to MAG proteins (check predicted TSVs); cannot rank genomes");
			return false;
		}

		try {
			Marker4ResourcePreparer.ensureSelectionResources(magDb, createGuiProgressCallback());
		} catch (IOException e) {
			logger.error(getTaskName() + ": failed to prepare Marker4 resources", e);
			return false;
		}
		File markerFasta = new File(magDb.getCurrentFile(), "Marker4.fasta");
		File genomesTaxaTsv = new File(magDb.getCurrentFile(), "genomes-taxa.tsv");
		if (!markerFasta.exists() || !genomesTaxaTsv.exists()) {
			logger.error(getTaskName() + ": Marker4.fasta or genomes-taxa.tsv not found for genome selection");
			return false;
		}

		updateMessage(getTaskName() + ": genome ranking for deep seed workflow...");
		File evidenceFile = new File(seedFolder, "genome_evidence.tsv");
		File rankedFile = new File(seedFolder, "genome_ranked.tsv");
		File selectedFile = new File(seedFolder, "genome_selected.tsv");
		DiaNNPrecursor[] withProArray = precursorsWithPro.toArray(new DiaNNPrecursor[0]);
		MagGenomeSelectionService.Request selectionRequest = new MagGenomeSelectionService.Request();
		selectionRequest.magDb = magDb;
		selectionRequest.genomeSelectorMode = MagDbItem.MARKER4;
		selectionRequest.genomeSelectionStrategy = MagGenomeSelectionService.SELECTION_STRATEGY_COMPREHENSIVE;
		selectionRequest.diaPrecursors = withProArray;
		selectionRequest.markerFasta = markerFasta;
		selectionRequest.genomesTaxaTsv = genomesTaxaTsv;
		selectionRequest.evidenceFile = evidenceFile;
		selectionRequest.rankedFile = rankedFile;
		selectionRequest.selectedFile = selectedFile;
		selectionRequest.maximumSelected = 300;
		try {
			MagGenomeSelectionService.selectGenomes(selectionRequest, createGuiProgressCallback());
		} catch (Exception e) {
			logger.error(getTaskName() + ": genome selection failed", e);
			return false;
		}
		logger.info(getTaskName() + ": genome ranking for deep seed completed (genome_evidence.tsv, genome_ranked.tsv, genome_selected.tsv)");
		return true;
	}

	@Override
	protected boolean exportReport() {
		return super.exportReport();
	}

	@Override
	protected Boolean call() throws Exception {
		if (metaPar instanceof MetaParameterDiaMarker
				&& ((MetaParameterDiaMarker) metaPar).getMarkerRefinementMode() == DiaRefinementMode.DEEP_PER_FILE) {
			updateMessage("Seed Deep Per File is not supported yet.");
			logger.error("Seed Deep Per File is not supported yet.");
			return false;
		}

		String resultPath = metaPar.getResult();
		if (resultPath == null || resultPath.trim().isEmpty()) {
			logger.error(getTaskName() + ": result path is null or empty");
			return false;
		}
		File base = new File(resultPath, BASE_FOLDER);
		if (!base.exists()) {
			base.mkdirs();
		}
		// data/ and per-run <runId>/ folders are only for future deep model; not created in current workflow
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

		logger.info("Step 4: Iterative search (first search + genome rank and identification)...");
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
		updateSecondaryProgress(1.0);
		logger.info(getTaskName() + ": finished - all workflow steps completed successfully");
		System.out.println(new java.text.SimpleDateFormat("yyyy-MM-dd HH:mm:ss").format(new java.util.Date()) + "\t" + getTaskName() + ": finished successfully");
		return true;
	}

	/**
	 * Helper to aggregate precursors from per_file_*.parquet in the seed folder when
	 * seed.parquet could not be generated.
	 */
	private DiaNNPrecursor[] readPrecursorsFromPerFileParquets(File dir, Double qValue) {
		String localCut = diannPar != null && diannPar.getCut() != null ? diannPar.getCut() : "KR";
		return DeepIterativeCommonSupport.readPrecursorsFromPerFileParquets(dir, qValue, diaFiles, localCut, logger,
				getTaskName());
	}

	@Override
	protected boolean iden() {
		File base = metaPar.getDbSearchResultFile();
		if (base == null) {
			logger.error(getTaskName() + ": base result folder is null");
			return false;
		}
		File seedFolder = getIterativeResultFolder();
		// Iterative rounds produce base/combined/*. The final quantifying parquet is
		// created in this method after the shared round loop completes.
		File rankedFile = new File(seedFolder, "genome_ranked.tsv");
		if (!rankedFile.exists()) {
			logger.error(getTaskName() + ": genome_ranked.tsv not found for deep seed workflow: {}", rankedFile.getAbsolutePath());
			return false;
		}

		String cut = diannPar != null && diannPar.getCut() != null ? diannPar.getCut() : "KR";
		DeepIterativeCommonSupport.IterativeRoundLoopResult loop;
		MetaParameterDiaMarker markerPar = (metaPar instanceof MetaParameterDiaMarker) ? (MetaParameterDiaMarker) metaPar
				: null;
		boolean usePeptideBasedRound = markerPar != null && markerPar.isUsePeptideBasedRound();
		if (usePeptideBasedRound) {
			double detectabilityThreshold = markerPar != null ? markerPar.getPeptideDetectabilityThreshold() : 0.5;
			loop = runPeptideBasedIterativeRoundLoop(rankedFile, seedFolder, cut, detectabilityThreshold);
		} else {
			loop = DeepIterativeCommonSupport.runIterativeRoundLoop(
					rankedFile, getGenomesPerRound(), MAX_ITERATIVE_ROUNDS, "Deep seed round", cut,
					initialSearchSeparate, initialSearchSeparate && goodDiaFiles != null ? goodDiaFiles.length : 0,
					this::updateMessage,
					(round, genomes) -> runSeedRoundAppendPreviousPeptides(round, seedFolder, genomes),
					(inc) -> {
						StopRuleResult local = computeStopRoundFromIncrements(inc);
						return new DeepIterativeStopRule.StopRuleResult(local.stopRound, local.nRounds, local.kChord,
								local.kPiece, local.kKnee, local.k0, local.tMin);
					},
					logger, getTaskName());
		}

		if (loop.lastRound <= 0) {
			logger.error(getTaskName() + ": no iterative round completed");
			return false;
		}

		StopRuleResult lastStopRuleResult = loop.lastStopRuleResult == null ? null
				: new StopRuleResult(loop.lastStopRuleResult.stopRound, loop.lastStopRuleResult.nRounds,
						loop.lastStopRuleResult.kChord, loop.lastStopRuleResult.kPiece, loop.lastStopRuleResult.kKnee,
						loop.lastStopRuleResult.k0, loop.lastStopRuleResult.tMin);
		int lastRound = loop.lastRound;
		writeIterativeRoundsInfo(seedFolder, loop.roundsRecord, loop.peptidesInRoundRecord, loop.incrementsPerRound,
				loop.cumulativeRecord, loop.stopReason, lastStopRuleResult);

		if (usePeptideBasedRound) {
			File lastRoundDir = new File(seedFolder, "round" + lastRound);
			this.combinedFile = lastRoundDir;
			this.resultFolderFile = base;
			this.quan_pep_file = new File(lastRoundDir, "round" + lastRound + ".pr_matrix.tsv");
			this.quan_pro_file = new File(lastRoundDir, "round" + lastRound + ".pg_matrix.tsv");
			if (!this.quan_pep_file.exists()) {
				logger.error(getTaskName() + ": peptide-based deep seed last-round matrix missing: {}",
						this.quan_pep_file.getAbsolutePath());
				return false;
			}
			logger.info(getTaskName() + ": peptide-based deep seed completed at round {}. Using last-round outputs as final results.",
					Integer.valueOf(lastRound));
			return true;
		}

		File combinedDir = new File(base, "combined");
		if (!combinedDir.exists())
			combinedDir.mkdirs();
		File combinedParquet = new File(combinedDir, "combined.parquet");
		File combinedFasta = new File(combinedDir, "combined.fasta");
		File combinedSpeclib = new File(combinedDir, "combined.speclib");

		HashMap<String, HashSet<String>> pepProMap = buildPepToProteinsFromSeeds(seedFolder, lastRound);
		if (pepProMap.isEmpty()) {
			logger.error(getTaskName() + ": no peptides from round parquets for combined FASTA");
			return false;
		}
		if (!DeepIterativeCommonSupport.writeCombinedFastaFromPeptideProteinMap(combinedFasta, pepProMap, logger,
				getTaskName())) {
			return false;
		}
		if (!DeepIterativeCommonSupport.ensureCombinedSpeclib(combinedDir, combinedSpeclib, () -> {
				diannPar.setThreads(threadCount);
				diannPar.setExportDecoy(true);
				// DIA-NN call purpose: build post-round library for final quantifying search.
				diaNNTask.addTask(diannPar, combinedFasta.getAbsolutePath(), combinedSpeclib.getAbsolutePath(), false);
				diaNNTask.run(12);
			}, logger, getTaskName())) {
			logger.error(getTaskName() + ": failed to generate combined.speclib at " + combinedSpeclib.getAbsolutePath());
			return false;
		}
		String[] finalSearchFiles = initialSearchSeparate && goodDiaFiles != null ? goodDiaFiles : diaFiles;
		boolean combinedOk = DeepIterativeCommonSupport.ensureCombinedParquet(combinedParquet, () -> {
			diaLibSearchPar.setMbr(true);
			diaLibSearchPar.setDiaFiles(finalSearchFiles);
			diaLibSearchPar.setThreadCount(threadCount);
			diaLibSearchPar.setqValue(0.01);
			diannPar.setThreads(threadCount);
			diannPar.setMatrices(true);
			diannPar.setRelaxed_prot_inf(true);
			diannPar.setReanalyse(true);
			diannPar.setExportDecoy(true);
			diaNNTask.addTask(diannPar, diaLibSearchPar, combinedParquet.getAbsolutePath(),
					new String[] { combinedSpeclib.getAbsolutePath() });
			diaNNTask.run(finalSearchFiles.length * 12);
		});
		if (!combinedOk) {
			return false;
		}

		if (!combinedParquet.exists()) {
			logger.error(getTaskName() + ": combined deep seed search failed - combined.parquet not created");
			return false;
		}
		combinedFile = combinedDir;
		this.resultFolderFile = base;
		this.quan_pep_file = new File(combinedDir, "combined.pr_matrix.tsv");
		this.quan_pro_file = new File(combinedDir, "combined.pg_matrix.tsv");
		logger.info(getTaskName() + ": deep seed combined search completed: {} ({} bytes)",
				combinedParquet.getAbsolutePath(), combinedParquet.length());
		return true;
	}

	/**
	 * Writes iterative_rounds_info.tsv so users can see why the round search stopped.
	 * Columns: Round, Peptides_in_round, New_peptides, Cumulative_peptides, Stop_reason (on last row only).
	 * Stop is determined by relative gain: New_peptides/Cumulative_peptides < 1.5% (plateau on total discovery).
	 * If stopRule is non-null (e.g. last evaluation when rounds>=8), writes kChord, kPiece, kKnee, k0, tMin in comment lines.
	 */
	private void writeIterativeRoundsInfo(File resultDir, List<Integer> rounds, List<Integer> peptidesInRound,
			List<Integer> newPeptides, List<Integer> cumulative, String stopReason, StopRuleResult stopRule) {
		DeepIterativeStopRule.StopRuleResult shared = stopRule == null ? null
				: new DeepIterativeStopRule.StopRuleResult(stopRule.stopRound, stopRule.nRounds, stopRule.kChord,
						stopRule.kPiece, stopRule.kKnee, stopRule.k0, stopRule.tMin);
		DeepIterativeCommonSupport.writeIterativeRoundsInfo(resultDir, rounds, peptidesInRound, newPeptides, cumulative,
				stopReason, shared, CONSECUTIVE_ROUNDS_FOR_PLATEAU, logger, getTaskName());
	}

	private DeepIterativeCommonSupport.IterativeRoundLoopResult runPeptideBasedIterativeRoundLoop(File rankedFile,
			File seedFolder, String cut, double detectabilityThreshold) {
		final double earlyRelativeGainThreshold = 0.01;
		final int earlyConsecutiveRounds = 2;

		List<String> rankedGenomes = readAllGenomesFromRanked(rankedFile);
		if (rankedGenomes.isEmpty()) {
			return new DeepIterativeCommonSupport.IterativeRoundLoopResult(-1, 0, new ArrayList<>(), new ArrayList<>(),
					new ArrayList<>(), new ArrayList<>(), "empty_ranked_file", null);
		}

		HashMap<String, HashSet<String>> cumulativePepProMap = new HashMap<>();
		HashSet<String> seenIdentifiedPeptides = new HashSet<>();
		List<Integer> roundsRecord = new ArrayList<>();
		List<Integer> peptidesInRoundRecord = new ArrayList<>();
		List<Integer> incrementsPerRound = new ArrayList<>();
		List<Integer> cumulativeRecord = new ArrayList<>();
		String stopReason = "max_rounds_reached";
		DeepIterativeStopRule.StopRuleResult lastStopRuleResult = null;
		int cursor = 0;
		int lastRound = 0;
		int prevPeptideCount = 0;

		for (int round = 1; round <= MAX_ITERATIVE_ROUNDS; round++) {
			if (cursor >= rankedGenomes.size()) {
				stopReason = "end_of_ranked_list";
				break;
			}
			updateMessage(getTaskName() + ": peptide-based deep seed round " + round + " library and search");
			int beforeRoundSize = cumulativePepProMap.size();
			List<String> thisRoundGenomes = new ArrayList<>();
			while (cursor < rankedGenomes.size()) {
				String genome = rankedGenomes.get(cursor);
				thisRoundGenomes.add(genome);
				mergePredictedPeptidesForGenome(genome, detectabilityThreshold, cumulativePepProMap);
				cursor++;
				if (cumulativePepProMap.size() > PEPTIDE_ROUND_TRIGGER) {
					break;
				}
			}
			if (thisRoundGenomes.isEmpty()) {
				stopReason = "end_of_ranked_list";
				break;
			}

			File roundDir = new File(seedFolder, "round" + round);
			if (!roundDir.exists()) {
				roundDir.mkdir();
			}
			File combinedFasta = new File(roundDir, "combined.fasta");
			if (!writeFastaFromPeptideProteinMap(combinedFasta, cumulativePepProMap)) {
				logger.error(getTaskName() + ": round {}: failed to write peptide-based combined.fasta", round);
				return new DeepIterativeCommonSupport.IterativeRoundLoopResult(-1, prevPeptideCount, roundsRecord,
						peptidesInRoundRecord, incrementsPerRound, cumulativeRecord, "round_failed", lastStopRuleResult);
			}
			File combinedSpeclib = new File(roundDir, "combined.speclib");
			if (!DeepIterativeCommonSupport.ensureCombinedSpeclib(roundDir, combinedSpeclib, () -> {
				diannPar.setThreads(threadCount);
				diannPar.setExportDecoy(true);
				diaNNTask.addTask(diannPar, combinedFasta.getAbsolutePath(), combinedSpeclib.getAbsolutePath(), false);
				diaNNTask.run(12);
			}, logger, getTaskName())) {
				return new DeepIterativeCommonSupport.IterativeRoundLoopResult(-1, prevPeptideCount, roundsRecord,
						peptidesInRoundRecord, incrementsPerRound, cumulativeRecord, "round_failed", lastStopRuleResult);
			}

			File roundParquet = new File(roundDir, "round" + round + ".parquet");
			String[] searchFiles = initialSearchSeparate && goodDiaFiles != null ? goodDiaFiles : diaFiles;
			diaLibSearchPar.setMbr(false);
			diaLibSearchPar.setDiaFiles(searchFiles);
			diaLibSearchPar.setThreadCount(threadCount);
			diaLibSearchPar.setqValue(0.01);
			diannPar.setThreads(threadCount);
			diannPar.setMatrices(true);
			diannPar.setRelaxed_prot_inf(true);
			diannPar.setReanalyse(true);
			diannPar.setExportDecoy(true);
			diaNNTask.addTask(diannPar, diaLibSearchPar, roundParquet.getAbsolutePath(),
					new String[] { combinedSpeclib.getAbsolutePath() });
			diaNNTask.run(searchFiles.length * 12);
			if (!roundParquet.exists()) {
				logger.error(getTaskName() + ": round {}: search failed - {} not created", round,
						roundParquet.getAbsolutePath());
				return new DeepIterativeCommonSupport.IterativeRoundLoopResult(-1, prevPeptideCount, roundsRecord,
						peptidesInRoundRecord, incrementsPerRound, cumulativeRecord, "round_failed", lastStopRuleResult);
			}

			int peptideCount = DeepIterativeCommonSupport.countIdentifiedPeptides(roundParquet, cut, logger, getTaskName());
			int incThis = DeepIterativeCommonSupport.countNewPeptidesFromParquet(roundParquet, seenIdentifiedPeptides, cut,
					logger, getTaskName());
			mergeIdentifiedPeptidesFromParquet(roundParquet, cut, cumulativePepProMap);

			roundsRecord.add(Integer.valueOf(round));
			peptidesInRoundRecord.add(Integer.valueOf(peptideCount));
			incrementsPerRound.add(Integer.valueOf(incThis));
			cumulativeRecord.add(Integer.valueOf(seenIdentifiedPeptides.size()));
			lastRound = round;
			prevPeptideCount = peptideCount;

			logger.info(
					"{}: peptide-based round {} summary - genomes_added={}, cumulative_library_before={}, cumulative_library_after={}, identified={}, new={}, cumulative_identified={}",
					getTaskName(), Integer.valueOf(round), Integer.valueOf(thisRoundGenomes.size()),
					Integer.valueOf(beforeRoundSize), Integer.valueOf(cumulativePepProMap.size()), Integer.valueOf(peptideCount),
					Integer.valueOf(incThis), Integer.valueOf(seenIdentifiedPeptides.size()));

			if (round >= earlyConsecutiveRounds) {
				boolean allBelow = true;
				for (int k = 0; k < earlyConsecutiveRounds && allBelow; k++) {
					int idx = round - 1 - k;
					int cum = cumulativeRecord.get(idx);
					if (cum <= 0) {
						allBelow = false;
						break;
					}
					int incVal = incrementsPerRound.get(idx);
					double relativeGain = (double) incVal / (double) cum;
					if (relativeGain >= earlyRelativeGainThreshold) {
						allBelow = false;
					}
				}
				if (allBelow) {
					stopReason = "early_plateau_rule_1pct";
					break;
				}
			}
			if (round >= 8) {
				StopRuleResult local = computeStopRoundFromIncrements(incrementsPerRound);
				lastStopRuleResult = new DeepIterativeStopRule.StopRuleResult(local.stopRound, local.nRounds, local.kChord,
						local.kPiece, local.kKnee, local.k0, local.tMin);
				if (local.stopRound > 0 && local.stopRound == round) {
					stopReason = "knee_plateau_rule";
					break;
				}
			}
		}
		return new DeepIterativeCommonSupport.IterativeRoundLoopResult(lastRound, prevPeptideCount, roundsRecord,
				peptidesInRoundRecord, incrementsPerRound, cumulativeRecord, stopReason, lastStopRuleResult);
	}

	private List<String> readAllGenomesFromRanked(File rankedFile) {
		List<String> genomes = new ArrayList<>();
		try (BufferedReader r = new BufferedReader(new FileReader(rankedFile))) {
			if (r.readLine() == null) {
				return genomes;
			}
			String line;
			while ((line = r.readLine()) != null) {
				line = line.trim();
				if (line.isEmpty() || line.startsWith("#")) {
					continue;
				}
				String[] parts = line.split("\t");
				if (parts.length > 1 && !parts[1].trim().isEmpty()) {
					genomes.add(parts[1].trim());
				}
			}
		} catch (IOException e) {
			logger.warn(getTaskName() + ": error reading ranked file {}", rankedFile.getAbsolutePath(), e);
		}
		return genomes;
	}

	private void mergePredictedPeptidesForGenome(String genome, double detectabilityThreshold,
			HashMap<String, HashSet<String>> cumulativePepProMap) {
		File predictedFile = new File(new File(magDb.getCurrentFile(), "predicted"), genome + ".predicted.tsv");
		if (!predictedFile.isFile()) {
			logger.warn(getTaskName() + ": predicted peptide file missing for genome {}", genome);
			return;
		}
		int added = 0;
		try (BufferedReader reader = new BufferedReader(new FileReader(predictedFile))) {
			reader.readLine();
			String line;
			while ((line = reader.readLine()) != null) {
				String[] cs = line.split("\t");
				if (cs.length < 3) {
					continue;
				}
				double detectability;
				try {
					detectability = Double.parseDouble(cs[2]);
				} catch (NumberFormatException nfe) {
					continue;
				}
				if (detectability < detectabilityThreshold) {
					continue;
				}
				String protein = cs[0];
				String peptide = cs[1];
				if (peptide == null || peptide.isEmpty()) {
					continue;
				}
				HashSet<String> proSet = cumulativePepProMap.computeIfAbsent(peptide, k -> new HashSet<>());
				int before = proSet.size();
				proSet.add(protein);
				if (before == 0) {
					added++;
				}
			}
		} catch (IOException e) {
			logger.warn(getTaskName() + ": error reading predicted peptide file {}", predictedFile.getAbsolutePath(), e);
		}
		logger.info(getTaskName() + ": genome {} added {} predicted peptides (threshold={})",
				genome, Integer.valueOf(added), Double.valueOf(detectabilityThreshold));
	}

	private boolean writeFastaFromPeptideProteinMap(File outputFasta, HashMap<String, HashSet<String>> pepProMap) {
		return DeepIterativeCommonSupport.writeCombinedFastaFromPeptideProteinMap(outputFasta, pepProMap, logger,
				getTaskName());
	}

	private void mergeIdentifiedPeptidesFromParquet(File parquetFile, String cut,
			HashMap<String, HashSet<String>> cumulativePepProMap) {
		try {
			DiaNNParquetReader reader = new DiaNNParquetReader(parquetFile.getAbsolutePath(), cut);
			DiaNNPrecursor[] precursors = reader.getDiaNNPrecursors(0.01);
			if (precursors == null) {
				return;
			}
			for (DiaNNPrecursor p : precursors) {
				if (p == null || p.getIsDecoy()) {
					continue;
				}
				String seq = p.getSeqString();
				if (seq == null || seq.isEmpty()) {
					continue;
				}
				HashSet<String> proSet = cumulativePepProMap.computeIfAbsent(seq, k -> new HashSet<>());
				String[] pros = p.getProGroups();
				if (pros != null) {
					for (String pro : pros) {
						if (pro != null && !pro.isEmpty()) {
							proSet.add(pro);
						}
					}
				}
			}
		} catch (Exception e) {
			logger.warn(getTaskName() + ": failed to merge identified peptides from {}", parquetFile.getAbsolutePath(), e);
		}
	}

	// ===== Helper methods for deep seed iterative rounds (adapted from DiaMarkerDeepPerFileTask) =====

	/** Write combined FASTA from genome IDs using MAG DB original folder (genome.faa or genome.fasta). */
	private boolean writeCombinedFastaFromGenomes(List<String> genomes, File outputFasta) {
		File dbFolder = magDb.getOriginalDbFolder();
		if (dbFolder == null || !dbFolder.isDirectory()) {
			logger.error(getTaskName() + ": MAG original DB folder not available");
			return false;
		}
		try (PrintWriter w = new PrintWriter(outputFasta)) {
			for (String genome : genomes) {
				File genomeFaa = new File(dbFolder, genome + ".faa");
				if (!genomeFaa.exists()) {
					genomeFaa = new File(dbFolder, genome + ".fasta");
				}
				if (!genomeFaa.exists()) {
					logger.warn(getTaskName() + ": genome file not found: {}", genome);
					continue;
				}
				try (BufferedReader reader = new BufferedReader(new FileReader(genomeFaa))) {
					String line;
					while ((line = reader.readLine()) != null) {
						w.println(line);
					}
				}
			}
		} catch (IOException e) {
			logger.error(getTaskName() + ": error writing combined FASTA", e);
			return false;
		}
		return outputFasta.exists() && outputFasta.length() > 0;
	}

	private int runSeedRoundAppendPreviousPeptides(int round, File seedFolder, List<String> genomes) {
		String[] searchFiles = initialSearchSeparate && goodDiaFiles != null ? goodDiaFiles : diaFiles;
		return DeepIterativeCommonSupport.runAppendPreviousPeptidesRound(round, seedFolder, genomes, searchFiles, false,
				"deep seed", initialSearchSeparate,
				this::writeCombinedFastaFromGenomes,
				this::writeCombinedFastaFromGenomesAndPreviousPeptides,
				(combinedFasta, combinedSpeclib) -> {
					diannPar.setThreads(threadCount);
					diaNNTask.addTask(diannPar, combinedFasta.getAbsolutePath(), combinedSpeclib.getAbsolutePath(), false);
					diaNNTask.run(12);
				},
				(roundParquet, combinedSpeclibPath, files, reanalyse) -> {
					if (initialSearchSeparate) {
						File roundDir = roundParquet.getParentFile();
						int maxConcurrent = DeepIterativeCommonSupport.computeMaxConcurrentTasks();
						for (int i = 0; i < files.length; i++) {
							File pf = new File(roundDir, "per_file_" + i + ".parquet");
							if (pf.exists() && pf.length() > 0) {
								continue;
							}
							DiaNNNoQuantLibrarySearchService.configureNoQuantSearch(
									diannPar, diaLibSearchPar, new String[] { files[i] }, threadCount, reanalyse, false);
							diaNNTask.addNoQuantLibTask(diannPar, diaLibSearchPar, pf.getAbsolutePath(),
									combinedSpeclibPath);
						}
						diaNNTask.run(maxConcurrent, files.length * 12);
					} else {
						DiaNNNoQuantLibrarySearchService.configureNoQuantSearch(
								diannPar, diaLibSearchPar, files, threadCount, reanalyse, false);
						diaNNTask.addNoQuantLibTask(diannPar, diaLibSearchPar, roundParquet.getAbsolutePath(),
								combinedSpeclibPath);
						diaNNTask.run(files.length * 12);
					}
				},
				this::countIdentifiedPeptides, logger, getTaskName());
	}

	/**
	 * Write combined FASTA: genome proteins from the given genomes, then append one entry per unique identified peptide
	 * from previousRoundParquet (q≤0.01, non-decoy). Protein name = first non-decoy proGroup or a fallback.
	 */
	private boolean writeCombinedFastaFromGenomesAndPreviousPeptides(List<String> genomes, File previousRoundParquet,
			File outputFasta) {
		File dbFolder = magDb.getOriginalDbFolder();
		if (dbFolder == null || !dbFolder.isDirectory()) {
			logger.error(getTaskName() + ": MAG original DB folder not available");
			return false;
		}
		String cut = diannPar != null && diannPar.getCut() != null ? diannPar.getCut() : "KR";
		DiaNNPrecursor[] precursors;
		if (previousRoundParquet.exists()) {
			DiaNNParquetReader reader = new DiaNNParquetReader(previousRoundParquet.getAbsolutePath(), cut);
			precursors = reader.getDiaNNPrecursors(0.01);
		} else {
			File prevDir = previousRoundParquet.getParentFile();
			String[] filesForRead = initialSearchSeparate && goodDiaFiles != null ? goodDiaFiles : diaFiles;
			precursors = DeepIterativeCommonSupport.readPrecursorsFromPerFileParquets(
					prevDir, 0.01, filesForRead, cut, logger, getTaskName());
			logger.info(getTaskName() + ": read {} precursors from per-file parquets in {}",
					precursors != null ? precursors.length : 0, prevDir.getAbsolutePath());
		}
		Set<String> addedPeptideSequences = new HashSet<>();
		try (PrintWriter w = new PrintWriter(outputFasta)) {
			for (String genome : genomes) {
				File genomeFaa = new File(dbFolder, genome + ".faa");
				if (!genomeFaa.exists())
					genomeFaa = new File(dbFolder, genome + ".fasta");
				if (!genomeFaa.exists()) {
					logger.warn(getTaskName() + ": genome file not found: {}", genome);
					continue;
				}
				try (BufferedReader reader = new BufferedReader(new FileReader(genomeFaa))) {
					String line;
					while ((line = reader.readLine()) != null)
						w.println(line);
				}
			}
			if (precursors != null) {
				for (DiaNNPrecursor p : precursors) {
					if (p == null || p.getIsDecoy())
						continue;
					String seq = p.getSeqString();
					if (seq == null || seq.isEmpty())
						continue;
					if (addedPeptideSequences.contains(seq))
						continue;
					String[] proGroups = p.getProGroups();
					String proteinId = null;
					if (proGroups != null) {
						for (String pro : proGroups) {
							if (pro == null || pro.isEmpty())
								continue;
							if (pro.startsWith("REV_") || pro.startsWith("rev_"))
								continue;
							proteinId = pro;
							break;
						}
					}
					if (proteinId == null)
						proteinId = "pep_" + seq.length();
					w.println(">" + proteinId);
					w.println(seq);
					addedPeptideSequences.add(seq);
				}
			}
			logger.info(getTaskName() + ": appended {} unique peptides from previous round to combined.fasta",
					addedPeptideSequences.size());
		} catch (IOException e) {
			logger.error(getTaskName() + ": error writing combined FASTA (genomes + peptides)", e);
			return false;
		}
		return outputFasta.exists() && outputFasta.length() > 0;
	}

	/** Relative gain threshold: stop when New_peptides/Cumulative_peptides < this (e.g. 0.015 = 1.5%). */
	private static final double RELATIVE_GAIN_THRESHOLD = 0.015;
	/** Require this many consecutive rounds below threshold before stopping (trend-based, not single dip). */
	private static final int CONSECUTIVE_ROUNDS_FOR_PLATEAU = 2;

	/** Result of stop-rule evaluation: stop round (0 = don't stop) plus knee indices for iterative_rounds_info. */
	public static class StopRuleResult {
		public final int stopRound;
		public final int nRounds;
		public final int kChord;
		public final int kPiece;
		public final int kKnee;
		public final int k0;
		public final int tMin;

		public StopRuleResult(int stopRound, int nRounds, int kChord, int kPiece, int kKnee, int k0, int tMin) {
			this.stopRound = stopRound;
			this.nRounds = nRounds;
			this.kChord = kChord;
			this.kPiece = kPiece;
			this.kKnee = kKnee;
			this.k0 = k0;
			this.tMin = tMin;
		}
	}

	/**
	 * Iterative DIA stop rule: plateau on cumulative discovery by trend, not a single dip.
	 * Stop only when relative gain (New_peptides/Cumulative_peptides) is below 1.5% for CONSECUTIVE_ROUNDS_FOR_PLATEAU consecutive rounds,
	 * so we do not stop while the trend is still rising (e.g. one round at 4.7% but previous round at 5.3%).
	 * inc is 0-based list of length n (rounds 1..n). Returns result with stopRound in 1..n, or 0 if no plateau detected yet; knee indices set when n>=8.
	 */
	private StopRuleResult computeStopRoundFromIncrements(List<Integer> inc) {
		DeepIterativeStopRule.StopRuleResult shared = DeepIterativeStopRule.computeStopRoundFromIncrements(inc,
				RELATIVE_GAIN_THRESHOLD, CONSECUTIVE_ROUNDS_FOR_PLATEAU,
				DeepIterativeStopRule.AlgorithmProfile.DEEP_COMPREHENSIVE);
		return new StopRuleResult(shared.stopRound, shared.nRounds, shared.kChord, shared.kPiece, shared.kKnee,
				shared.k0, shared.tMin);
	}

	/**
	 * Build peptide→proteins map by scanning all round{N}/round{N}.parquet files (q≤0.01, non-decoy).
	 * Deep seed uses seed/roundN folders (next to genome_ranked.tsv).
	 */
	private HashMap<String, HashSet<String>> buildPepToProteinsFromSeeds(File seedFolder, int lastRound) {
		if (seedFolder == null || lastRound <= 0) {
			return new HashMap<>();
		}
		return DiaDeepProteinMappingService.buildPeptideProteinMap(
				DiaDeepProteinMappingService.Request.deepTogether(seedFolder, lastRound,
						initialSearchSeparate && goodDiaFiles != null ? goodDiaFiles : diaFiles,
						diannPar != null ? diannPar.getCut() : "KR"),
				logger, getTaskName());
	}

	private int countIdentifiedPeptides(File parquetFile) {
		String cut = diannPar != null && diannPar.getCut() != null ? diannPar.getCut() : "KR";
		return DeepIterativeCommonSupport.countIdentifiedPeptides(parquetFile, cut, logger, getTaskName());
	}

	/**
	 * Rebuild peptide→protein mapping from deep-seed round parquets in the report phase, so that
	 * each precursor's ProGroups reflect all proteins from iterative rounds (not just the single
	 * protein used in the peptide-only combined FASTA).
	 * Called from exportReport path only.
	 */
	@Override
	protected void enrichPrecursorsWithPeptideProteinMapping(DiaNNPrecursor[] precursors, File parquetFolder) {
		File seedFolder = getIterativeResultFolder();
		DiaDeepProteinMappingService.remapPrecursors(precursors,
				DiaDeepProteinMappingService.Request.deepTogether(seedFolder, 0,
						initialSearchSeparate && goodDiaFiles != null ? goodDiaFiles : diaFiles,
						diannPar != null ? diannPar.getCut() : "KR"),
				logger, getTaskName());
	}

	// (prediction-based candidate peptide helpers removed for the current non-deep implementation)

}
