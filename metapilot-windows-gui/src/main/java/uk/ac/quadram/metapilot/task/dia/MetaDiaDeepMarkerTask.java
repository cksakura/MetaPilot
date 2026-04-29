package uk.ac.quadram.metapilot.task.dia;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.logging.log4j.Logger;

import javafx.scene.control.ProgressBar;

import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNParquetReader;
import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNPrecursor;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.core.prodb.FastaManager;
import uk.ac.quadram.metapilot.task.dia.par.MetaParameterDiaMarker;
import uk.ac.quadram.metapilot.task.mag.MarkerGenomeSelector;

/**
 * DIA Deep Marker workflow.
 * <p>
 * High-level flow:
 * <ol>
 * <li>Search marker library against DIA files to get marker evidence parquet.</li>
 * <li>Run marker-based genome ranking to produce {@code genome_ranked.tsv}.</li>
 * <li>Run iterative rounds (top-N genomes per round) where each round builds one
 *     temporary FASTA/speclib and performs non-quant identification.</li>
 * <li>Merge peptides from all rounds and run one final quantifying search to produce
 *     {@code combined.parquet}, {@code combined.pr_matrix.tsv}, {@code combined.pg_matrix.tsv}.</li>
 * </ol>
 *
 * DIA-NN usage map in this class:
 * <ul>
 * <li>{@code addNoQuantLibTask(... marker.parquet ...)}: marker-only identification for genome ranking input.</li>
 * <li>{@code addTask(... combined.fasta -> combined.speclib ...)}: build predicted library from FASTA.</li>
 * <li>{@code addNoQuantLibTask(... roundN.parquet ...)}: round-level non-quant identification.</li>
 * <li>{@code addTask(... combined.parquet ...)}: final quantifying search over all DIA files.</li>
 * </ul>
 */
public class MetaDiaDeepMarkerTask extends MetaDiaMarkerTask {

	private static final String TASK_NAME = "DIA Deep Marker analysis";
	/** Number of top genomes from genome_ranked.tsv per round. */
	private static final int GENOMES_PER_ROUND = 5;
	/** Stop iterative rounds when identification count does not increase (like Deep Comprehensive). */
	private static final int MAX_ITERATIVE_ROUNDS = 100;

	/** Last round index from iterative search; used in exportReport to rebuild peptide→protein mapping from round parquets. */
	private int lastRoundForExport;

	/** True when the initial marker search fell back to per-file search. */
	private boolean initialSearchSeparate;
	/** Subset of diaFiles that produced results. Set when initialSearchSeparate is true. */
	private String[] goodDiaFiles;

	public MetaDiaDeepMarkerTask(MetaParameterDiaMarker markerParams, MetaResourceParameters mrp,
			ProgressBar bar1, ProgressBar bar2) {
		super(markerParams, mrp, bar1, bar2);
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
	protected void initResultFolderIfNeeded() {
		if (metaPar == null) {
			return;
		}
		if (metaPar.getDbSearchResultFile() != null) {
			this.resultFolderFile = metaPar.getDbSearchResultFile();
			return;
		}
		if (!(metaPar instanceof MetaParameterDiaMarker)) {
			super.initResultFolderIfNeeded();
			return;
		}
		MetaParameterDiaMarker markerParams = (MetaParameterDiaMarker) metaPar;
		File base = new File(markerParams.getResult() != null ? markerParams.getResult() : ".");
		File dir = new File(base, "marker_deep");
		if (!dir.exists()) {
			dir.mkdirs();
			logger.info("Created result directory for DIA Deep Marker workflow: {}", dir.getAbsolutePath());
		}
		markerParams.setDbSearchResultFile(dir);
		this.resultFolderFile = dir;
	}

	/**
	 * Same as marker workflow: search marker library → marker.parquet, then run
	 * MarkerGenomeSelector to produce genome_evidence.tsv, genome_ranked.tsv, genome_selected.tsv.
	 * Stops after that; does not build combined.fasta or combined.speclib.
	 */
	@Override
	protected boolean iterativeSearch() {
		logger.info("Starting iterative search for DIA Deep Marker (search + genome selection only)...");
		System.out.println(new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss").format(new Date()) + "\t" + getTaskName()
				+ ": iterative search started");

		File markerLibFile = magDb.resolveMarker4LibraryFile();
		if (markerLibFile == null || !markerLibFile.exists()) {
			logger.error(getTaskName() + ": Marker library file not found for iterative search");
			System.err.println(new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss").format(new Date()) + "\t" + getTaskName()
					+ ": Marker library file not found");
			return false;
		}

		File markerFile = getIterativeResultFolder();
		if (!markerFile.exists()) {
			markerFile.mkdir();
		}

		File tempParquetFile = new File(markerFile, "marker.parquet");
		if (tempParquetFile.exists()) {
			logger.info("Marker library search result already exists in {}", tempParquetFile.getAbsolutePath());
		} else if (hasExistingPerFileParquets(markerFile)) {
			initialSearchSeparate = true;
			goodDiaFiles = collectGoodDiaFilesFromPerFileParquets(markerFile);
			int maxConcurrent = DeepIterativeCommonSupport.computeMaxConcurrentTasks();
			logger.info(getTaskName()
					+ ": per-file parquets already exist from previous run; {} good files; rounds will use per-file search (max concurrent: {})",
					goodDiaFiles.length, maxConcurrent);
		} else {
			diaLibSearchPar.setMbr(false);
			diaLibSearchPar.setDiaFiles(diaFiles);
			diaLibSearchPar.setThreadCount(threadCount);
			diaLibSearchPar.setqValue(0.01);
			diannPar.setThreads(threadCount);
			diannPar.setMatrices(false);
			diannPar.setRelaxed_prot_inf(true);
			diannPar.setReanalyse(true);
			diannPar.setExportDecoy(false);
			logger.info("Searching marker library against DIA files...");
			logger.info("  Marker library: {}", markerLibFile.getAbsolutePath());
			logger.info("  DIA files count: {}", diaFiles != null ? diaFiles.length : 0);
			diaNNTask.addNoQuantLibTask(diannPar, diaLibSearchPar, tempParquetFile.getAbsolutePath(),
					markerLibFile.getAbsolutePath());
			diaNNTask.run(diaFiles.length * 12);

			if (!tempParquetFile.exists() && diaFiles != null && diaFiles.length > 1) {
				List<String> damaged = runPerFileFallbackAndRetryCombinedNoQuant(
						tempParquetFile.getAbsolutePath(),
						new String[] { markerLibFile.getAbsolutePath() });
				if (!damaged.isEmpty()) {
					logger.warn(getTaskName() + ": The following raw files produced no results (may be damaged): {}",
							damaged);
					System.err.println(new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss").format(new Date()) + "\t"
							+ getTaskName() + ": The following raw files produced no results (may be damaged): "
							+ String.join(", ", damaged));
				}
				java.util.Set<String> damagedSet = new java.util.HashSet<>(damaged);
				java.util.ArrayList<String> goodList = new java.util.ArrayList<>();
				for (String f : diaFiles) {
					if (!damagedSet.contains(f)) {
						goodList.add(f);
					}
				}
				if (!goodList.isEmpty()) {
					initialSearchSeparate = true;
					goodDiaFiles = goodList.toArray(new String[0]);
					int maxConcurrent = DeepIterativeCommonSupport.computeMaxConcurrentTasks();
					logger.info(getTaskName()
							+ ": initial search was separate; {} good files, {} damaged; rounds will use per-file search (max concurrent: {})",
							goodDiaFiles.length, damaged.size(), maxConcurrent);
				}
			}
		}
		if (!initialSearchSeparate) {
			goodDiaFiles = diaFiles;
		}

		// Genome selection only (no combined.fasta / combined.speclib)
		logger.info("Reading PSMs from iterative search results...");
		DiaNNPrecursor[] precursors;
		if (tempParquetFile.exists()) {
			DiaNNParquetReader parquetReader = new DiaNNParquetReader(tempParquetFile.getAbsolutePath(),
					diannPar.getCut());
			precursors = parquetReader.getDiaNNPrecursors(0.01);
		} else {
			logger.warn(getTaskName()
					+ ": Marker library combined parquet not created; using per-file parquet results for genome selection");
			precursors = readPrecursorsFromPerFileParquets(markerFile, 0.01);
		}
		logger.info("Found {} precursors with Global.Q.Value < 0.01", precursors != null ? precursors.length : 0);

		if (precursors == null || precursors.length == 0) {
			logger.warn(getTaskName() + ": No precursors found in iterative search, cannot select genomes");
			return false;
		}

		MetaParameterDiaMarker markerParams = (metaPar instanceof MetaParameterDiaMarker) ? (MetaParameterDiaMarker) metaPar : null;
		// Map markerMode to MarkerGenomeSelector knee mode (FOCUS vs COMPREHENSIVE)
		String markerMode = markerParams != null ? markerParams.getMarkerMode() : null;
		String mode;
		if ("COMPREHENSIVE".equalsIgnoreCase(markerMode) || "DEEP".equalsIgnoreCase(markerMode)) {
			mode = MarkerGenomeSelector.KNEE_MODE_COMPREHENSIVE;
		} else {
			mode = MarkerGenomeSelector.KNEE_MODE_FOCUS;
		}
		logger.info("Selecting genomes using MarkerGenomeSelector ({} mode)...", mode);

		File currentFile = magDb.getCurrentFile();
		File markerFasta = new File(currentFile, "Marker4.fasta");
		File genomesTaxaTsv = new File(currentFile, "genomes-taxa.tsv");
		if (!markerFasta.exists() || !genomesTaxaTsv.exists()) {
			logger.error(getTaskName() + ": Marker FASTA or genomes-taxa.tsv not found for genome selection");
			System.err.println(new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss").format(new Date()) + "\t" + getTaskName()
					+ ": Marker FASTA or genomes-taxa.tsv not found");
			return false;
		}

		MarkerGenomeSelector selector = new MarkerGenomeSelector();
		selector.setKneeSelectionMode(mode);
		File evidenceFile = new File(markerFile, "genome_evidence.tsv");
		File rankedFile = new File(markerFile, "genome_ranked.tsv");
		File selectedFile = new File(markerFile, "genome_selected.tsv");
		int maximumSelected = 300;
		boolean enableAutoRep = true;
		selector.selectMarkerGenomesWithKnee(precursors, markerFasta, genomesTaxaTsv,
				evidenceFile.getAbsolutePath(), rankedFile.getAbsolutePath(), selectedFile.getAbsolutePath(),
				maximumSelected, enableAutoRep);

		logger.info(getTaskName() + ": genome selection completed (genome_evidence.tsv, genome_ranked.tsv, genome_selected.tsv); stopping before combined.fasta/speclib");
		System.out.println(new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss").format(new Date()) + "\t" + getTaskName()
				+ ": iterative search finished (genome selection only)");
		return true;
	}

	/**
	 * Reuse the same per-file parquet aggregation helper as the marker workflow.
	 */
	private DiaNNPrecursor[] readPrecursorsFromPerFileParquets(File dir, Double qValue) {
		String cut = diannPar != null && diannPar.getCut() != null ? diannPar.getCut() : "KR";
		return DeepIterativeCommonSupport.readPrecursorsFromPerFileParquets(dir, qValue, diaFiles, cut, logger,
				getTaskName());
	}

	/**
	 * Round controller after {@code genome_ranked.tsv} is available.
	 * Uses shared round-loop logic to:
	 * <ol>
	 * <li>take next genome window from rank list,</li>
	 * <li>run one non-quant identification round,</li>
	 * <li>evaluate stop rule based on incremental peptide gain,</li>
	 * <li>perform final combined quantifying search.</li>
	 * </ol>
	 */
	@Override
	protected boolean iden() {
		File markerFile = getIterativeResultFolder();
		File rankedFile = new File(markerFile, "genome_ranked.tsv");
		if (!rankedFile.exists()) {
			logger.error(getTaskName() + ": genome_ranked.tsv not found: {}", rankedFile.getAbsolutePath());
			return false;
		}

		String cut = diannPar != null && diannPar.getCut() != null ? diannPar.getCut() : "KR";
		DeepIterativeCommonSupport.IterativeRoundLoopResult loop = DeepIterativeCommonSupport.runIterativeRoundLoop(
				rankedFile, GENOMES_PER_ROUND, MAX_ITERATIVE_ROUNDS, "Round", cut,
				initialSearchSeparate, initialSearchSeparate ? goodDiaFiles.length : 0,
				this::updateMessage,
				(round, genomes) -> runRoundAppendPreviousPeptides(round, markerFile, genomes),
				(inc) -> {
					StopRuleResult local = computeStopRoundFromIncrements(inc);
					return new DeepIterativeStopRule.StopRuleResult(local.stopRound, local.nRounds, local.kChord,
							local.kPiece, local.kKnee, local.k0, local.tMin);
				},
				logger, getTaskName());

		if (loop.lastRound <= 0) {
			logger.error(getTaskName() + ": no round completed");
			return false;
		}

		StopRuleResult lastStopRuleResult = loop.lastStopRuleResult == null ? null
				: new StopRuleResult(loop.lastStopRuleResult.stopRound, loop.lastStopRuleResult.nRounds,
						loop.lastStopRuleResult.kChord, loop.lastStopRuleResult.kPiece, loop.lastStopRuleResult.kKnee,
						loop.lastStopRuleResult.k0, loop.lastStopRuleResult.tMin);
		int lastRound = loop.lastRound;
		int prevPeptideCount = loop.previousPeptideCount;
		writeIterativeRoundsInfo(markerFile, loop.roundsRecord, loop.peptidesInRoundRecord, loop.incrementsPerRound,
				loop.cumulativeRecord, loop.stopReason, lastStopRuleResult);

		logger.info(getTaskName() + ": iterative rounds finished at round {} ({} peptides); running combined search using combined FASTA built from all round parquets...", lastRound, prevPeptideCount);

		// Combined folder: build one combined FASTA from all round{N}.parquet, then combined.speclib, then search all raw files once.
		File combinedDir = new File(metaPar.getDbSearchResultFile(), "combined");
		if (!combinedDir.exists()) combinedDir.mkdirs();
		File combinedParquet = new File(combinedDir, "combined.parquet");
		File combinedPepMatrix = new File(combinedDir, "combined.pr_matrix.tsv");
		File combinedFasta = new File(combinedDir, "combined.fasta");
		File combinedSpeclib = new File(combinedDir, "combined.speclib");

		if (!combinedParquet.exists() || !combinedPepMatrix.exists()) {
			// 1) Aggregate peptide→all proteins from all round{N}.parquet
			HashMap<String, HashSet<String>> pepProMap = buildPepToProteinsFromRounds(markerFile, lastRound);
			if (pepProMap.isEmpty()) {
				logger.error(getTaskName() + ": no peptides collected from round parquets for combined FASTA");
				return false;
			}
			this.lastRoundForExport = lastRound;

			// 2) Write combined FASTA: one entry per unique peptide, using one of its proteins as the FASTA protein name.
			boolean hadCombinedFasta = combinedFasta.exists() && combinedFasta.length() > 0;
			if (!DeepIterativeCommonSupport.writeCombinedFastaFromPeptideProteinMap(combinedFasta, pepProMap, logger,
					getTaskName())) {
				return false;
			}
			if (!hadCombinedFasta) {
				logger.info(getTaskName() + ": combined FASTA written: {} ({} bytes)", combinedFasta.getAbsolutePath(), combinedFasta.length());
			} else {
				logger.info(getTaskName() + ": using existing combined FASTA at {}", combinedFasta.getAbsolutePath());
			}

			// DIA-NN call purpose: build library for the final quantifying pass
			// from round-aggregated FASTA (combined.fasta -> combined.speclib).
			if (!DeepIterativeCommonSupport.ensureCombinedSpeclib(combinedDir, combinedSpeclib, () -> {
					diannPar.setThreads(threadCount);
					diannPar.setExportDecoy(true);
					diaNNTask.addTask(diannPar, combinedFasta.getAbsolutePath(), combinedSpeclib.getAbsolutePath(), false);
					diaNNTask.run(12);
				}, logger, getTaskName())) {
				return false;
			}
			logger.info(getTaskName() + ": combined.speclib ready: {} ({} bytes)", combinedSpeclib.getAbsolutePath(), combinedSpeclib.length());

			String[] finalSearchFiles = initialSearchSeparate ? goodDiaFiles : diaFiles;
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
				logger.info(getTaskName() + ": combined search using combined.speclib {} → {} ({} files)",
						combinedSpeclib.getAbsolutePath(), combinedParquet.getAbsolutePath(), finalSearchFiles.length);
				diaNNTask.addTask(diannPar, diaLibSearchPar, combinedParquet.getAbsolutePath(),
						new String[] { combinedSpeclib.getAbsolutePath() });
				diaNNTask.run(finalSearchFiles.length * 12);
			});
			if (!combinedOk) {
				logger.error(getTaskName() + ": combined search failed - combined.parquet not created");
				return false;
			}
		} else {
			logger.info(getTaskName() + ": combined result already exists, skipping combined search");
			if (lastRoundForExport <= 0) {
				int gfc = goodDiaFiles != null ? goodDiaFiles.length : (diaFiles != null ? diaFiles.length : 0);
				for (int r = 1; r <= MAX_ITERATIVE_ROUNDS; r++) {
					File rd = new File(markerFile, "round" + r);
					if (rd.isDirectory() && (new File(rd, "round" + r + ".parquet").isFile()
							|| DeepIterativeCommonSupport.hasRoundPerFileParquets(rd, gfc))) {
						lastRoundForExport = r;
					}
				}
			}
		}

		combinedFile = combinedDir;
		this.quan_pep_file = new File(combinedDir, "combined.pr_matrix.tsv");
		this.quan_pro_file = new File(combinedDir, "combined.pg_matrix.tsv");
		logger.info(getTaskName() + ": combined search completed: {} ({} bytes); unified post-analysis will use combined folder", combinedParquet.getAbsolutePath(), combinedParquet.length());
		System.out.println(new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss").format(new Date()) + "\t" + getTaskName() + ": identification finished (combined)");
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

	/**
	 * Build peptide→proteins map by scanning all round{N}/round{N}.parquet files (q≤0.01, non-decoy).
	 */
	private HashMap<String, HashSet<String>> buildPepToProteinsFromRounds(File markerFile, int lastRound) {
		HashMap<String, HashSet<String>> pepProMap = new HashMap<>();
		if (markerFile == null || lastRound <= 0) {
			return pepProMap;
		}
		String cut = diannPar != null && diannPar.getCut() != null ? diannPar.getCut() : "KR";
		String[] filesForRead = initialSearchSeparate ? goodDiaFiles : diaFiles;
		for (int round = 1; round <= lastRound; round++) {
			File roundDir = new File(markerFile, "round" + round);
			DiaNNPrecursor[] precursors = null;
			File roundParquet = new File(roundDir, "round" + round + ".parquet");
			if (roundParquet.isFile()) {
				try {
					DiaNNParquetReader reader = new DiaNNParquetReader(roundParquet.getAbsolutePath(), cut);
					precursors = reader.getDiaNNPrecursors(0.01);
				} catch (Exception e) {
					logger.warn(getTaskName() + ": failed to read round parquet " + roundParquet.getAbsolutePath(), e);
				}
			} else {
				precursors = DeepIterativeCommonSupport.readPrecursorsFromPerFileParquets(
						roundDir, 0.01, filesForRead, cut, logger, getTaskName());
			}
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
		}
		logger.info(getTaskName() + ": aggregated {} peptides from round parquets for combined FASTA", pepProMap.size());
		return pepProMap;
	}

	/**
	 * Rebuild peptide→protein mapping from round parquets in the report phase, so that
	 * each precursor's ProGroups reflect all proteins from iterative rounds (not just
	 * the single protein from the combined FASTA). Called from exportReport path only.
	 */
	@Override
	protected void enrichPrecursorsWithPeptideProteinMapping(DiaNNPrecursor[] precursors, File parquetFolder) {
		File markerFile = getIterativeResultFolder();
		if (markerFile == null || !markerFile.exists()) {
			return;
		}
		int lastRound = lastRoundForExport;
		if (lastRound <= 0) {
			int gfc = goodDiaFiles != null ? goodDiaFiles.length : (diaFiles != null ? diaFiles.length : 0);
			for (int r = 1; r <= MAX_ITERATIVE_ROUNDS; r++) {
				File rd = new File(markerFile, "round" + r);
				if (rd.isDirectory() && (new File(rd, "round" + r + ".parquet").isFile()
						|| DeepIterativeCommonSupport.hasRoundPerFileParquets(rd, gfc))) {
					lastRound = r;
				}
			}
		}
		if (lastRound <= 0) {
			return;
		}
		HashMap<String, HashSet<String>> pepProMap = buildPepToProteinsFromRounds(markerFile, lastRound);
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
		logger.info("{}: deep-marker remap from round parquets – peptides with mapping: {}, proteins before remap: {}, proteins after remap: {}",
				getTaskName(), pepProMap.size(), proteinsBefore.size(), proteinsAfter.size());
	}

	/**
	 * Round strategy used by deep marker: one library per round.
	 * FASTA = [proteins from 5 genomes] + [identified peptides from previous round as single-protein entries].
	 * Returns identified peptide count or -1 on failure.
	 */
	private int runRoundAppendPreviousPeptides(int round, File markerFile, List<String> genomes) {
		String[] searchFiles = initialSearchSeparate ? goodDiaFiles : diaFiles;
		return DeepIterativeCommonSupport.runAppendPreviousPeptidesRound(round, markerFile, genomes, searchFiles, true,
				"append strategy", initialSearchSeparate,
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
							diaLibSearchPar.setMbr(false);
							diaLibSearchPar.setDiaFiles(new String[] { files[i] });
							diaLibSearchPar.setThreadCount(threadCount);
							diaLibSearchPar.setqValue(0.01);
							diannPar.setThreads(threadCount);
							diannPar.setMatrices(false);
							diannPar.setRelaxed_prot_inf(true);
							diannPar.setReanalyse(reanalyse);
							diannPar.setExportDecoy(false);
							diaNNTask.addNoQuantLibTask(diannPar, diaLibSearchPar, pf.getAbsolutePath(), combinedSpeclibPath);
						}
						diaNNTask.run(maxConcurrent, files.length * 12);
					} else {
						diaLibSearchPar.setMbr(false);
						diaLibSearchPar.setDiaFiles(files);
						diaLibSearchPar.setThreadCount(threadCount);
						diaLibSearchPar.setqValue(0.01);
						diannPar.setThreads(threadCount);
						diannPar.setMatrices(false);
						diannPar.setRelaxed_prot_inf(true);
						diannPar.setReanalyse(reanalyse);
						diannPar.setExportDecoy(false);
						diaNNTask.addNoQuantLibTask(diannPar, diaLibSearchPar, roundParquet.getAbsolutePath(), combinedSpeclibPath);
						diaNNTask.run(files.length * 12);
					}
				},
				this::countIdentifiedPeptides, logger, getTaskName());
	}

	/**
	 * Write combined FASTA: first genome proteins from the given genomes, then append one entry per unique identified peptide from previousRoundParquet (q≤0.01, non-decoy).
	 * Each peptide is written as ">proteinId\nsequence\n" using the first non-decoy protein from that peptide's proGroups so DIA-NN can link peptide to protein.
	 */
	private boolean writeCombinedFastaFromGenomesAndPreviousPeptides(List<String> genomes, File previousRoundParquet, File outputFasta) {
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
			String[] filesForRead = initialSearchSeparate ? goodDiaFiles : diaFiles;
			precursors = DeepIterativeCommonSupport.readPrecursorsFromPerFileParquets(
					prevDir, 0.01, filesForRead, cut, logger, getTaskName());
			logger.info(getTaskName() + ": read {} precursors from per-file parquets in {}",
					precursors != null ? precursors.length : 0, prevDir.getAbsolutePath());
		}
		Set<String> addedPeptideSequences = new HashSet<>();
		try (PrintWriter w = new PrintWriter(outputFasta)) {
			for (String genome : genomes) {
				File genomeFaa = new File(dbFolder, genome + ".faa");
				if (!genomeFaa.exists()) genomeFaa = new File(dbFolder, genome + ".fasta");
				if (!genomeFaa.exists()) {
					logger.warn(getTaskName() + ": genome file not found: {}", genome);
					continue;
				}
				try (BufferedReader reader = new BufferedReader(new FileReader(genomeFaa))) {
					String line;
					while ((line = reader.readLine()) != null) w.println(line);
				}
			}
			if (precursors != null) {
				for (DiaNNPrecursor p : precursors) {
					if (p == null || p.getIsDecoy()) continue;
					String seq = p.getSeqString();
					if (seq == null || seq.isEmpty()) continue;
					if (addedPeptideSequences.contains(seq)) continue;
					String[] proGroups = p.getProGroups();
					String proteinId = null;
					if (proGroups != null) {
						for (String pro : proGroups) {
							if (pro == null || pro.isEmpty()) continue;
							if (pro.startsWith(FastaManager.shuffle) || pro.startsWith("REV_") || pro.startsWith("rev_")) continue;
							proteinId = pro;
							break;
						}
					}
					if (proteinId == null) proteinId = "pep_" + seq.length();
					w.println(">" + proteinId);
					w.println(seq);
					addedPeptideSequences.add(seq);
				}
			}
			logger.info(getTaskName() + ": appended {} unique peptides from previous round to combined.fasta", addedPeptideSequences.size());
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
				DeepIterativeStopRule.AlgorithmProfile.DEEP_MARKER);
		return new StopRuleResult(shared.stopRound, shared.nRounds, shared.kChord, shared.kPiece, shared.kKnee,
				shared.k0, shared.tMin);
	}

	/** Count unique peptide sequences from parquet (q≤0.01, exclude decoys). */
	private int countIdentifiedPeptides(File parquetFile) {
		String cut = diannPar != null && diannPar.getCut() != null ? diannPar.getCut() : "KR";
		return DeepIterativeCommonSupport.countIdentifiedPeptides(parquetFile, cut, logger, getTaskName());
	}

	/**
	 * Check whether the marker folder already has per-file parquets from a previous separate search.
	 */
	private boolean hasExistingPerFileParquets(File markerDir) {
		if (diaFiles == null || diaFiles.length == 0 || markerDir == null) {
			return false;
		}
		for (int i = 0; i < diaFiles.length; i++) {
			File pf = new File(markerDir, "per_file_" + i + ".parquet");
			if (pf.exists() && pf.length() > 0) {
				return true;
			}
		}
		return false;
	}

	/**
	 * Collect the subset of diaFiles whose per-file parquets exist (non-damaged files).
	 */
	private String[] collectGoodDiaFilesFromPerFileParquets(File markerDir) {
		java.util.ArrayList<String> goodList = new java.util.ArrayList<>();
		if (diaFiles != null && markerDir != null) {
			for (int i = 0; i < diaFiles.length; i++) {
				File pf = new File(markerDir, "per_file_" + i + ".parquet");
				if (pf.exists() && pf.length() > 0) {
					goodList.add(diaFiles[i]);
				}
			}
		}
		return goodList.toArray(new String[0]);
	}

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
}
