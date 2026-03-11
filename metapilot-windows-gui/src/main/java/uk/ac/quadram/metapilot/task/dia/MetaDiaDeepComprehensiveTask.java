package uk.ac.quadram.metapilot.task.dia;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.logging.log4j.Logger;

import javafx.scene.control.ProgressBar;

import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNParquetReader;
import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNPrecursor;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.task.dia.par.MetaParameterDiaMarker;
import uk.ac.quadram.metapilot.task.mag.MarkerGenomeSelector;

/**
 * DIA "Seed deep" workflow task for the <b>seed strategy</b>.
 * <p>
 * High-level flow:
 * <ol>
 * <li>Reuse seed first search from {@link MetaDiaSeedIterativeTask} to generate
 *     {@code seed.parquet} (or per-file parquets on fallback).</li>
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
public class MetaDiaDeepComprehensiveTask extends MetaDiaSeedIterativeTask {

	private static final String TASK_NAME = "DIA Seed Deep analysis";
	private static final String BASE_FOLDER = "seed_deep";
	/** Folder name for seed output (same as parent). */
	private static final String SEED_FOLDER = "seed";
	private static final String SEED_PARQUET = "seed.parquet";

	/** Number of top genomes from genome_ranked.tsv per round. */
	private static final int GENOMES_PER_ROUND = 5;
	/** Maximum number of iterative genome rounds. */
	private static final int MAX_ITERATIVE_ROUNDS = 100;

	public MetaDiaDeepComprehensiveTask(MetaParameterDiaMarker markerParams, MetaResourceParameters mrp,
			ProgressBar bar1, ProgressBar bar2) {
		super(markerParams, mrp, bar1, bar2);
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
		List<DiaNNPrecursor> precursorsWithPro = MarkerGenomeSelector.setProteinInfoFromPredictedTsVs(precursors,
				predictedFolder);
		if (precursorsWithPro != null && !precursorsWithPro.isEmpty()) {
			logger.info(getTaskName() + ": mapped {} / {} precursors to MAG proteins", precursorsWithPro.size(),
					precursors.length);
		}
		if (precursorsWithPro == null || precursorsWithPro.isEmpty()) {
			logger.error(getTaskName() + ": no precursors mapped to MAG proteins (check predicted TSVs); cannot rank genomes");
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
		MarkerGenomeSelector selector = new MarkerGenomeSelector();
		selector.setKneeSelectionMode(MarkerGenomeSelector.KNEE_MODE_COMPREHENSIVE);
		DiaNNPrecursor[] withProArray = precursorsWithPro.toArray(new DiaNNPrecursor[0]);
		selector.selectMarkerGenomesWithKnee(withProArray, markerFasta, genomesTaxaTsv,
				evidenceFile.getAbsolutePath(), rankedFile.getAbsolutePath(), selectedFile.getAbsolutePath(),
				300, true);
		logger.info(getTaskName() + ": genome ranking for deep seed completed (genome_evidence.tsv, genome_ranked.tsv, genome_selected.tsv)");
		return true;
	}

	@Override
	protected boolean exportReport() {
		return super.exportReport();
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
		updateSecondaryProgress(0.0);
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
		DeepIterativeCommonSupport.IterativeRoundLoopResult loop = DeepIterativeCommonSupport.runIterativeRoundLoop(
				rankedFile, GENOMES_PER_ROUND, MAX_ITERATIVE_ROUNDS, "Deep seed round", cut, this::updateMessage,
				(round, genomes) -> runSeedRoundAppendPreviousPeptides(round, seedFolder, genomes),
				(inc) -> {
					StopRuleResult local = computeStopRoundFromIncrements(inc);
					return new DeepIterativeStopRule.StopRuleResult(local.stopRound, local.nRounds, local.kChord,
							local.kPiece, local.kKnee, local.k0, local.tMin);
				},
				logger, getTaskName());

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
		boolean combinedOk = DeepIterativeCommonSupport.ensureCombinedParquet(combinedParquet, () -> {
			diaLibSearchPar.setMbr(true);
			diaLibSearchPar.setDiaFiles(diaFiles);
			diaLibSearchPar.setThreadCount(threadCount);
			diaLibSearchPar.setqValue(0.01);
			diannPar.setThreads(threadCount);
			diannPar.setMatrices(true);
			diannPar.setRelaxed_prot_inf(true);
			diannPar.setReanalyse(true);
			diannPar.setExportDecoy(true);
			// DIA-NN call purpose: final quantifying search across all configured DIA files.
			diaNNTask.addTask(diannPar, diaLibSearchPar, combinedParquet.getAbsolutePath(),
					new String[] { combinedSpeclib.getAbsolutePath() });
			diaNNTask.run(diaFiles.length * 12);
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

	// ===== Helper methods for deep seed iterative rounds (adapted from MetaDiaDeepMarkerTask) =====

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
		return DeepIterativeCommonSupport.runAppendPreviousPeptidesRound(round, seedFolder, genomes, diaFiles, false,
				"deep seed", this::writeCombinedFastaFromGenomes,
				this::writeCombinedFastaFromGenomesAndPreviousPeptides,
				(combinedFasta, combinedSpeclib) -> {
					diannPar.setThreads(threadCount);
					// DIA-NN call purpose: round-level library build
					// (roundN/combined.fasta -> roundN/combined.speclib).
					diaNNTask.addTask(diannPar, combinedFasta.getAbsolutePath(), combinedSpeclib.getAbsolutePath(), false);
					diaNNTask.run(12);
				},
				(roundParquet, combinedSpeclibPath, files, reanalyse) -> {
					// DIA-NN call purpose: round-level non-quant identification using the
					// round-specific library; output is roundN/roundN.parquet.
					diaLibSearchPar.setMbr(false);
					diaLibSearchPar.setDiaFiles(files);
					diaLibSearchPar.setThreadCount(threadCount);
					diaLibSearchPar.setqValue(0.01);
					diannPar.setThreads(threadCount);
					diannPar.setMatrices(false);
					diannPar.setRelaxed_prot_inf(true);
					diannPar.setReanalyse(reanalyse);
					diannPar.setExportDecoy(false);
					diaNNTask.addNoQuantLibTask(diannPar, diaLibSearchPar, roundParquet.getAbsolutePath(),
							combinedSpeclibPath);
					diaNNTask.run(files.length * 12);
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
			DiaNNParquetReader reader = new DiaNNParquetReader(previousRoundParquet.getAbsolutePath(), cut);
			DiaNNPrecursor[] precursors = reader.getDiaNNPrecursors(0.01);
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
		HashMap<String, HashSet<String>> pepProMap = new HashMap<>();
		if (seedFolder == null || lastRound <= 0) {
			return pepProMap;
		}
		String cut = diannPar != null && diannPar.getCut() != null ? diannPar.getCut() : "KR";
		for (int round = 1; round <= lastRound; round++) {
			File roundDir = new File(seedFolder, "round" + round);
			File roundParquet = new File(roundDir, "round" + round + ".parquet");
			if (!roundParquet.isFile()) {
				continue;
			}
			try {
				DiaNNParquetReader reader = new DiaNNParquetReader(roundParquet.getAbsolutePath(), cut);
				DiaNNPrecursor[] precursors = reader.getDiaNNPrecursors(0.01);
				if (precursors == null) {
					continue;
				}
				for (DiaNNPrecursor p : precursors) {
					if (p == null || p.getIsDecoy()) {
						continue;
					}
					String pep = p.getSeqString();
					if (pep == null || pep.isEmpty()) {
						continue;
					}
					String[] pros = p.getProGroups();
					if (pros == null || pros.length == 0) {
						continue;
					}
					HashSet<String> set = pepProMap.computeIfAbsent(pep, k -> new HashSet<>());
					for (String pro : pros) {
						if (pro != null && !pro.isEmpty()) {
							set.add(pro);
						}
					}
				}
			} catch (Exception e) {
				logger.warn(getTaskName() + ": failed to read deep seed round parquet " + roundParquet.getAbsolutePath(), e);
			}
		}
		logger.info(getTaskName() + ": aggregated {} peptides from deep seed round parquets for combined FASTA", pepProMap.size());
		return pepProMap;
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
		if (seedFolder == null || !seedFolder.exists()) {
			return;
		}
		int lastRound = 0;
		for (int r = 1; r <= MAX_ITERATIVE_ROUNDS; r++) {
			File rd = new File(seedFolder, "round" + r);
			if (rd.isDirectory() && new File(rd, "round" + r + ".parquet").isFile()) {
				lastRound = r;
			}
		}
		if (lastRound <= 0) {
			return;
		}
		HashMap<String, HashSet<String>> pepProMap = buildPepToProteinsFromSeeds(seedFolder, lastRound);
		if (pepProMap.isEmpty() || precursors == null) {
			return;
		}
		// Count distinct proteins before remap
		HashSet<String> proteinsBefore = new HashSet<>();
		for (DiaNNPrecursor pp : precursors) {
			if (pp == null) continue;
			String[] pros = pp.getProGroups();
			if (pros == null) continue;
			for (String p : pros) {
				if (p != null && !p.isEmpty()) {
					proteinsBefore.add(p);
				}
			}
		}

		for (DiaNNPrecursor pp : precursors) {
			if (pp == null) {
				continue;
			}
			String seq = pp.getSeqString();
			if (seq == null || seq.isEmpty()) {
				continue;
			}
			HashSet<String> pros = pepProMap.get(seq);
			if (pros != null && !pros.isEmpty()) {
				pp.setProGroups(pros.toArray(new String[0]));
			}
		}
		// Count distinct proteins after remap
		HashSet<String> proteinsAfter = new HashSet<>();
		for (DiaNNPrecursor pp : precursors) {
			if (pp == null) continue;
			String[] pros = pp.getProGroups();
			if (pros == null) continue;
			for (String p : pros) {
				if (p != null && !p.isEmpty()) {
					proteinsAfter.add(p);
				}
			}
		}
		logger.info("{}: deep-seed remap from seed rounds – peptides with mapping: {}, proteins before remap: {}, proteins after remap: {}",
				getTaskName(), pepProMap.size(), proteinsBefore.size(), proteinsAfter.size());
	}

	// (prediction-based candidate peptide helpers removed for the current non-deep implementation)

}

