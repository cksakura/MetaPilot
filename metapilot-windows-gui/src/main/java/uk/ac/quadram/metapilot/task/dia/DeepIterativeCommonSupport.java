package uk.ac.quadram.metapilot.task.dia;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.lang.management.ManagementFactory;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

import org.apache.logging.log4j.Logger;

import com.sun.management.OperatingSystemMXBean;

import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNParquetReader;
import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNPrecursor;

/**
 * Shared helper utilities for deep iterative DIA workflows.
 */
public final class DeepIterativeCommonSupport {

	@FunctionalInterface
	public interface TaskAction {
		void run();
	}

	@FunctionalInterface
	public interface WriteFastaFromGenomes {
		boolean write(List<String> genomes, File outputFasta);
	}

	@FunctionalInterface
	public interface WriteFastaFromGenomesAndPrevious {
		boolean write(List<String> genomes, File previousRoundParquet, File outputFasta);
	}

	@FunctionalInterface
	public interface BuildSpeclibAction {
		void run(File combinedFasta, File combinedSpeclib);
	}

	@FunctionalInterface
	public interface RunNoQuantSearchAction {
		void run(File roundParquet, String combinedSpeclibPath, String[] searchFiles, boolean reanalyse);
	}

	@FunctionalInterface
	public interface CountPeptidesAction {
		int count(File roundParquet);
	}

	@FunctionalInterface
	public interface UpdateMessageAction {
		void update(String message);
	}

	@FunctionalInterface
	public interface RunRoundAction {
		int run(int round, List<String> genomes);
	}

	@FunctionalInterface
	public interface ComputeStopRuleAction {
		DeepIterativeStopRule.StopRuleResult compute(List<Integer> increments);
	}

	public static final class IterativeRoundLoopResult {
		public final int lastRound;
		public final int previousPeptideCount;
		public final List<Integer> roundsRecord;
		public final List<Integer> peptidesInRoundRecord;
		public final List<Integer> incrementsPerRound;
		public final List<Integer> cumulativeRecord;
		public final String stopReason;
		public final DeepIterativeStopRule.StopRuleResult lastStopRuleResult;

		public IterativeRoundLoopResult(int lastRound, int previousPeptideCount, List<Integer> roundsRecord,
				List<Integer> peptidesInRoundRecord, List<Integer> incrementsPerRound, List<Integer> cumulativeRecord,
				String stopReason, DeepIterativeStopRule.StopRuleResult lastStopRuleResult) {
			this.lastRound = lastRound;
			this.previousPeptideCount = previousPeptideCount;
			this.roundsRecord = roundsRecord;
			this.peptidesInRoundRecord = peptidesInRoundRecord;
			this.incrementsPerRound = incrementsPerRound;
			this.cumulativeRecord = cumulativeRecord;
			this.stopReason = stopReason;
			this.lastStopRuleResult = lastStopRuleResult;
		}
	}

	private DeepIterativeCommonSupport() {
	}

	/**
	 * Max concurrent DIA-NN tasks based on system RAM: totalRAM_GB / 16.
	 * Each DIA-NN instance uses ~16 GB; this limits parallel tasks to avoid OOM.
	 */
	public static int computeMaxConcurrentTasks() {
		long totalMemoryGB = ((OperatingSystemMXBean) ManagementFactory.getOperatingSystemMXBean())
				.getTotalMemorySize() / 1024 / 1024 / 1024;
		int maxTasks = (int) (totalMemoryGB / 16);
		return Math.max(maxTasks, 1);
	}

	/**
	 * Count new (previously unseen) peptides from per-file parquets in a round directory.
	 * Mirrors {@link #countNewPeptidesFromParquet} but reads from per_file_*.parquet.
	 */
	public static int countNewPeptidesFromRoundPerFileParquets(File roundDir, int fileCount,
			HashSet<String> seenPeptides, String cut, Logger logger, String taskName) {
		if (roundDir == null || !roundDir.isDirectory()) {
			return 0;
		}
		String localCut = cut != null ? cut : "KR";
		int newCount = 0;
		for (int i = 0; i < fileCount; i++) {
			File pf = new File(roundDir, "per_file_" + i + ".parquet");
			if (!pf.exists() || pf.length() <= 0) {
				continue;
			}
			try {
				DiaNNParquetReader reader = new DiaNNParquetReader(pf.getAbsolutePath(), localCut);
				DiaNNPrecursor[] precursors = reader.getDiaNNPrecursors(0.01);
				if (precursors != null) {
					for (DiaNNPrecursor p : precursors) {
						if (p == null || p.getIsDecoy()) {
							continue;
						}
						String seq = p.getSeqString();
						if (seq == null || seq.isEmpty()) {
							continue;
						}
						if (!seenPeptides.contains(seq)) {
							seenPeptides.add(seq);
							newCount++;
						}
					}
				}
			} catch (Exception e) {
				logger.warn(taskName + ": failed to count new peptides from " + pf.getAbsolutePath(), e);
			}
		}
		return newCount;
	}

	/**
	 * Count unique identified peptides from per-file parquets in a round directory.
	 * Mirrors {@link #countIdentifiedPeptides} but reads from per_file_*.parquet.
	 */
	public static int countIdentifiedPeptidesFromRoundPerFileParquets(File roundDir, int fileCount,
			String cut, Logger logger, String taskName) {
		if (roundDir == null || !roundDir.isDirectory()) {
			return 0;
		}
		String localCut = cut != null ? cut : "KR";
		HashSet<String> peptides = new HashSet<>();
		for (int i = 0; i < fileCount; i++) {
			File pf = new File(roundDir, "per_file_" + i + ".parquet");
			if (!pf.exists() || pf.length() <= 0) {
				continue;
			}
			try {
				DiaNNParquetReader reader = new DiaNNParquetReader(pf.getAbsolutePath(), localCut);
				DiaNNPrecursor[] precursors = reader.getDiaNNPrecursors(0.01);
				if (precursors != null) {
					for (DiaNNPrecursor p : precursors) {
						if (p == null || p.getIsDecoy()) {
							continue;
						}
						String seq = p.getSeqString();
						if (seq != null && !seq.isEmpty()) {
							peptides.add(seq);
						}
					}
				}
			} catch (Exception e) {
				logger.warn(taskName + ": failed to read per-file parquet " + pf.getAbsolutePath(), e);
			}
		}
		return peptides.size();
	}

	/**
	 * Check whether a round directory contains at least one non-empty per-file parquet.
	 */
	public static boolean hasRoundPerFileParquets(File roundDir, int fileCount) {
		if (roundDir == null || !roundDir.isDirectory()) {
			return false;
		}
		for (int i = 0; i < fileCount; i++) {
			File pf = new File(roundDir, "per_file_" + i + ".parquet");
			if (pf.exists() && pf.length() > 0) {
				return true;
			}
		}
		return false;
	}

	public static DiaNNPrecursor[] readPrecursorsFromPerFileParquets(File dir, Double qValue, String[] diaFiles,
			String cut, Logger logger, String taskName) {
		if (dir == null || !dir.isDirectory()) {
			logger.warn(taskName + ": readPrecursorsFromPerFileParquets: invalid round directory: {}", dir);
			return null;
		}
		String localCut = cut != null ? cut : "KR";
		List<DiaNNPrecursor> all = new ArrayList<>();
		int totalFilesRead = 0;
		int totalFilesFailed = 0;
		int totalPrecursors = 0;
		int expectedFiles = (diaFiles != null ? diaFiles.length : -1);
		logger.info(taskName
				+ ": reading per-file parquets from {} (qValue={}, cut={}, expectedFiles={})",
				dir.getAbsolutePath(), qValue, localCut, expectedFiles >= 0 ? expectedFiles : "auto-detect");
		for (int i = 0; i < (diaFiles != null ? diaFiles.length : Integer.MAX_VALUE); i++) {
			File pf = new File(dir, "per_file_" + i + ".parquet");
			if (!pf.exists()) {
				logger.info(taskName + ": per_file_{}.parquet not found in {}, stop scanning further indices", i,
						dir.getAbsolutePath());
				break;
			}
			if (pf.length() <= 0) {
				logger.warn(taskName + ": skip empty parquet file: {}", pf.getAbsolutePath());
				continue;
			}
			try {
				DiaNNParquetReader reader = new DiaNNParquetReader(pf.getAbsolutePath(), localCut);
				DiaNNPrecursor[] prec = (qValue == null) ? reader.getDiaNNPrecursors() : reader.getDiaNNPrecursors(qValue);
				int fileCount = (prec == null ? 0 : prec.length);
				totalFilesRead++;
				totalPrecursors += fileCount;
				logger.info(taskName + ": per-file result {} -> {} precursors (qValue={})", pf.getName(), fileCount,
						qValue);
				if (prec != null) {
					for (DiaNNPrecursor p : prec) {
						if (p != null) {
							all.add(p);
						}
					}
				}
			} catch (Exception e) {
				totalFilesFailed++;
				logger.warn(taskName + ": failed to read per-file parquet " + pf.getAbsolutePath(), e);
			}
		}
		logger.info(taskName + ": per-file precursor read summary: filesRead={}, filesFailed={}, totalPrecursors={}",
				totalFilesRead, totalFilesFailed, totalPrecursors);
		return all.isEmpty() ? null : all.toArray(new DiaNNPrecursor[0]);
	}

	public static List<String> readGenomesFromRanked(File rankedTsv, int offset, int count, Logger logger, String taskName) {
		List<String> result = new ArrayList<>();
		try (BufferedReader r = new BufferedReader(new FileReader(rankedTsv))) {
			if (r.readLine() == null) {
				return result;
			}
			int skipped = 0;
			String line;
			while ((line = r.readLine()) != null && result.size() < count) {
				line = line.trim();
				if (line.isEmpty() || line.startsWith("#")) {
					continue;
				}
				String[] parts = line.split("\t");
				if (parts.length > 1 && !parts[1].trim().isEmpty()) {
					if (skipped < offset) {
						skipped++;
						continue;
					}
					result.add(parts[1].trim());
				}
			}
		} catch (IOException e) {
			logger.warn(taskName + ": error reading ranked file {}", rankedTsv.getAbsolutePath(), e);
		}
		return result;
	}

	public static int countNewPeptidesFromParquet(File parquetFile, HashSet<String> seenPeptides, String cut, Logger logger,
			String taskName) {
		if (parquetFile == null || !parquetFile.exists()) {
			logger.warn(taskName + ": countNewPeptidesFromParquet: parquet file not found: {}", parquetFile);
			return 0;
		}
		String localCut = cut != null ? cut : "KR";
		int newCount = 0;
		try {
			DiaNNParquetReader reader = new DiaNNParquetReader(parquetFile.getAbsolutePath(), localCut);
			DiaNNPrecursor[] precursors = reader.getDiaNNPrecursors(0.01);
			if (precursors != null) {
				for (DiaNNPrecursor p : precursors) {
					if (p == null || p.getIsDecoy()) {
						continue;
					}
					String seq = p.getSeqString();
					if (seq == null || seq.isEmpty()) {
						continue;
					}
					if (!seenPeptides.contains(seq)) {
						seenPeptides.add(seq);
						newCount++;
					}
				}
			}
		} catch (Exception e) {
			logger.warn(taskName + ": failed to count new peptides from {}", parquetFile.getAbsolutePath(), e);
		}
		return newCount;
	}

	public static int countIdentifiedPeptides(File parquetFile, String cut, Logger logger, String taskName) {
		if (parquetFile == null || !parquetFile.exists()) {
			return 0;
		}
		String localCut = cut != null ? cut : "KR";
		HashSet<String> peptides = new HashSet<>();
		try {
			DiaNNParquetReader reader = new DiaNNParquetReader(parquetFile.getAbsolutePath(), localCut);
			DiaNNPrecursor[] precursors = reader.getDiaNNPrecursors(0.01);
			if (precursors != null) {
				for (DiaNNPrecursor p : precursors) {
					if (p == null || p.getIsDecoy()) {
						continue;
					}
					String seq = p.getSeqString();
					if (seq == null || seq.isEmpty()) {
						continue;
					}
					peptides.add(seq);
				}
			}
		} catch (Exception e) {
			logger.warn(taskName + ": failed to read parquet " + parquetFile.getAbsolutePath(), e);
		}
		return peptides.size();
	}

	public static void writeIterativeRoundsInfo(File resultDir, List<Integer> rounds, List<Integer> peptidesInRound,
			List<Integer> newPeptides, List<Integer> cumulative, String stopReason,
			DeepIterativeStopRule.StopRuleResult stopRule, int consecutiveRoundsForPlateau, Logger logger, String taskName) {
		File out = new File(resultDir, "iterative_rounds_info.tsv");
		try (PrintWriter w = new PrintWriter(out)) {
			w.println("# Stop_rule: trend-based plateau - relative gain New_peptides/Cumulative_peptides < 1.5% for "
					+ consecutiveRoundsForPlateau + " consecutive rounds");
			if (stopRule != null && stopRule.nRounds >= 8) {
				w.println("# Stop_rule_params: kChord=" + stopRule.kChord + "\tkPiece=" + stopRule.kPiece + "\tkKnee="
						+ stopRule.kKnee + "\tk0=" + stopRule.k0 + "\ttMin=" + stopRule.tMin + "\tnRounds_evaluated="
						+ stopRule.nRounds + "\tstopRound=" + stopRule.stopRound);
			}
			w.println("Round\tPeptides_in_round\tNew_peptides\tCumulative_peptides\tStop_reason");
			for (int i = 0; i < rounds.size(); i++) {
				String reason = (i == rounds.size() - 1) ? stopReason : "";
				w.println(rounds.get(i) + "\t" + peptidesInRound.get(i) + "\t" + newPeptides.get(i) + "\t"
						+ cumulative.get(i) + "\t" + reason);
			}
			logger.info(taskName + ": wrote iterative rounds info to {}", out.getAbsolutePath());
		} catch (IOException e) {
			logger.warn(taskName + ": could not write iterative_rounds_info.tsv: {}", e.getMessage());
		}
	}

	public static boolean writeCombinedFastaFromPeptideProteinMap(File combinedFasta,
			HashMap<String, HashSet<String>> pepProMap, Logger logger, String taskName) {
		if (combinedFasta.exists() && combinedFasta.length() > 0) {
			return true;
		}
		try (PrintWriter fastaWriter = new PrintWriter(combinedFasta)) {
			int id = 0;
			for (String pep : pepProMap.keySet()) {
				if (pep == null || pep.isEmpty()) {
					continue;
				}
				HashSet<String> pros = pepProMap.get(pep);
				if (pros == null || pros.isEmpty()) {
					continue;
				}
				String proName = pros.iterator().next();
				if (proName == null || proName.isEmpty()) {
					proName = "pro" + id;
				}
				fastaWriter.println(">" + proName);
				fastaWriter.println(pep);
				id++;
			}
		} catch (IOException e) {
			logger.error(taskName + ": error writing combined FASTA " + combinedFasta.getAbsolutePath(), e);
			return false;
		}
		return true;
	}

	public static boolean ensureCombinedSpeclib(File combinedDir, File combinedSpeclib,
			TaskAction generateSpeclibAction, Logger logger, String taskName) {
		if (!combinedSpeclib.exists() || combinedSpeclib.length() == 0) {
			File tempLibFile = new File(combinedDir, "combined.predicted.speclib");
			if (tempLibFile.exists()) {
				tempLibFile.renameTo(combinedSpeclib);
			} else {
				generateSpeclibAction.run();
				if (tempLibFile.exists() && !combinedSpeclib.exists()) {
					tempLibFile.renameTo(combinedSpeclib);
				}
			}
		}
		if (!combinedSpeclib.exists() || combinedSpeclib.length() == 0) {
			logger.error(taskName + ": failed to generate combined.speclib at {}", combinedSpeclib.getAbsolutePath());
			return false;
		}
		return true;
	}

	public static boolean ensureCombinedParquet(File combinedParquet, TaskAction runCombinedSearchAction) {
		if (!combinedParquet.exists()) {
			runCombinedSearchAction.run();
		}
		return combinedParquet.exists();
	}

	public static int runAppendPreviousPeptidesRound(int round, File roundBaseDir, List<String> genomes, String[] searchFiles,
			boolean reanalyse, String roundLabel, WriteFastaFromGenomes writeFromGenomes,
			WriteFastaFromGenomesAndPrevious writeFromGenomesAndPrevious, BuildSpeclibAction buildSpeclibAction,
			RunNoQuantSearchAction runNoQuantSearchAction, CountPeptidesAction countPeptidesAction, Logger logger,
			String taskName) {
		return runAppendPreviousPeptidesRound(round, roundBaseDir, genomes, searchFiles, reanalyse, roundLabel,
				false, writeFromGenomes, writeFromGenomesAndPrevious, buildSpeclibAction,
				runNoQuantSearchAction, countPeptidesAction, logger, taskName);
	}

	/**
	 * Overload supporting separate (per-file) search mode.
	 * When {@code separateSearch} is true, per-file parquets are expected instead of a single roundN.parquet.
	 * The {@code writeFromGenomesAndPrevious} lambda should handle reading from per-file parquets when
	 * the passed parquet file does not exist.
	 */
	public static int runAppendPreviousPeptidesRound(int round, File roundBaseDir, List<String> genomes, String[] searchFiles,
			boolean reanalyse, String roundLabel, boolean separateSearch,
			WriteFastaFromGenomes writeFromGenomes,
			WriteFastaFromGenomesAndPrevious writeFromGenomesAndPrevious, BuildSpeclibAction buildSpeclibAction,
			RunNoQuantSearchAction runNoQuantSearchAction, CountPeptidesAction countPeptidesAction, Logger logger,
			String taskName) {
		File roundDir = new File(roundBaseDir, "round" + round);
		if (!roundDir.exists()) {
			roundDir.mkdir();
		}

		File combinedFasta = new File(roundDir, "combined.fasta");
		if (round == 1) {
			if (!writeFromGenomes.write(genomes, combinedFasta)) {
				logger.error(taskName + ": round {}: failed to write combined.fasta", round);
				return -1;
			}
			logger.info(taskName + ": round {} ({}): combined.fasta from {} genomes only", round, roundLabel,
					genomes.size());
		} else {
			File prevFolder = new File(roundBaseDir, "round" + (round - 1));
			File prevParquet = new File(prevFolder, "round" + (round - 1) + ".parquet");
			if (!prevParquet.exists() && !(separateSearch && hasRoundPerFileParquets(prevFolder, searchFiles.length))) {
				logger.error(taskName + ": round {}: previous round parquet not found: {}", round,
						prevParquet.getAbsolutePath());
				return -1;
			}
			if (!writeFromGenomesAndPrevious.write(genomes, prevParquet, combinedFasta)) {
				logger.error(taskName + ": round {}: failed to write combined.fasta (genomes + previous peptides)",
						round);
				return -1;
			}
			logger.info(taskName + ": round {} ({}): combined.fasta = {} genomes + peptides from previous round", round,
					roundLabel, genomes.size());
		}
		logger.info(taskName + ": round {}: combined.fasta: {} ({} bytes)", round, combinedFasta.getAbsolutePath(),
				combinedFasta.length());

		File combinedSpeclib = new File(roundDir, "combined.speclib");
		if (!ensureCombinedSpeclib(roundDir, combinedSpeclib, () -> buildSpeclibAction.run(combinedFasta, combinedSpeclib),
				logger, taskName)) {
			logger.error(taskName + ": round {}: failed to generate combined.speclib", round);
			return -1;
		}
		logger.info(taskName + ": round {}: combined.speclib: {} ({} bytes)", round, combinedSpeclib.getAbsolutePath(),
				combinedSpeclib.length());

		File roundParquet = new File(roundDir, "round" + round + ".parquet");
		if (!roundParquet.exists() && !(separateSearch && hasRoundPerFileParquets(roundDir, searchFiles.length))) {
			runNoQuantSearchAction.run(roundParquet, combinedSpeclib.getAbsolutePath(), searchFiles, reanalyse);
		}
		if (separateSearch) {
			if (!hasRoundPerFileParquets(roundDir, searchFiles.length)) {
				logger.error(taskName + ": round {}: separate search failed - no per-file parquets created in {}", round, roundDir.getAbsolutePath());
				return -1;
			}
			int count = countIdentifiedPeptidesFromRoundPerFileParquets(roundDir, searchFiles.length,
					null, logger, taskName);
			logger.info(taskName + ": round {} (separate): unique peptides (q<=0.01, non-decoy): {}", round, count);
			return count;
		}
		if (!roundParquet.exists()) {
			logger.error(taskName + ": round {}: search failed - {} not created", round, roundParquet.getName());
			return -1;
		}

		int count = countPeptidesAction.count(roundParquet);
		logger.info(taskName + ": round {}: unique peptides (q<=0.01, non-decoy): {}", round, count);
		return count;
	}

	public static IterativeRoundLoopResult runIterativeRoundLoop(File rankedFile, int genomesPerRound, int maxRounds,
			String roundHeaderLabel, String cut, UpdateMessageAction updateMessageAction,
			RunRoundAction runRoundAction, ComputeStopRuleAction computeStopRuleAction, Logger logger, String taskName) {
		return runIterativeRoundLoop(rankedFile, genomesPerRound, maxRounds, roundHeaderLabel, cut,
				false, 0, updateMessageAction, runRoundAction, computeStopRuleAction, logger, taskName);
	}

	/**
	 * Overload supporting separate (per-file) search mode.
	 * When {@code separateSearch} is true, new-peptide counting reads from per-file parquets
	 * (per_file_*.parquet in each roundN/ directory) instead of a single roundN.parquet.
	 *
	 * @param separateSearch  true if rounds produced per-file parquets instead of combined roundN.parquet
	 * @param goodFileCount   number of good (non-damaged) DIA files; used to enumerate per-file parquets
	 */
	public static IterativeRoundLoopResult runIterativeRoundLoop(File rankedFile, int genomesPerRound, int maxRounds,
			String roundHeaderLabel, String cut, boolean separateSearch, int goodFileCount,
			UpdateMessageAction updateMessageAction,
			RunRoundAction runRoundAction, ComputeStopRuleAction computeStopRuleAction, Logger logger, String taskName) {
		final double earlyRelativeGainThreshold = 0.01;
		final int earlyConsecutiveRounds = 2;

		int prevPeptideCount = 0;
		int lastRound = 0;
		List<Integer> incPerRound = new ArrayList<>();
		HashSet<String> seenPeptides = new HashSet<>();
		List<Integer> roundsRecord = new ArrayList<>();
		List<Integer> peptidesInRoundRecord = new ArrayList<>();
		List<Integer> cumulativeRecord = new ArrayList<>();
		String stopReason = "max_rounds_reached";
		DeepIterativeStopRule.StopRuleResult lastStopRuleResult = null;

		for (int round = 1; round <= maxRounds; round++) {
			updateMessageAction.update(taskName + ": round " + round + " library and search");
			logger.info(taskName + ": ========== {} {} ==========", roundHeaderLabel, round);

			int offset = (round - 1) * genomesPerRound;
			List<String> genomes = readGenomesFromRanked(rankedFile, offset, genomesPerRound, logger, taskName);
			if (genomes.isEmpty()) {
				logger.info(taskName + ": round {}: no genomes at rank {}-{} (end of ranked list), stopping", round,
						offset + 1, offset + genomesPerRound);
				stopReason = "end_of_ranked_list";
				break;
			}
			logger.info(taskName + ": round {}: rank range {}-{}, genomes: {}", round, offset + 1,
					offset + genomes.size(), genomes);

			int peptideCount = runRoundAction.run(round, genomes);
			if (peptideCount < 0) {
				return new IterativeRoundLoopResult(-1, prevPeptideCount, roundsRecord, peptidesInRoundRecord,
						incPerRound, cumulativeRecord, stopReason, lastStopRuleResult);
			}

			File roundDir = new File(rankedFile.getParentFile(), "round" + round);
			int incThis;
			if (separateSearch) {
				incThis = countNewPeptidesFromRoundPerFileParquets(roundDir, goodFileCount, seenPeptides, cut,
						logger, taskName);
			} else {
				File roundParquet = new File(roundDir, "round" + round + ".parquet");
				incThis = countNewPeptidesFromParquet(roundParquet, seenPeptides, cut, logger, taskName);
			}
			incPerRound.add(incThis);
			roundsRecord.add(round);
			peptidesInRoundRecord.add(peptideCount);
			cumulativeRecord.add(seenPeptides.size());

			logger.info(taskName + ": round {}: identified {} peptides (previous round: {}), new = {}, cumulative = {}",
					round, peptideCount, prevPeptideCount, incThis, seenPeptides.size());

			if (round >= earlyConsecutiveRounds) {
				boolean allBelow = true;
				for (int k = 0; k < earlyConsecutiveRounds && allBelow; k++) {
					int idx = round - 1 - k;
					int cumulative = cumulativeRecord.get(idx);
					if (cumulative <= 0) {
						allBelow = false;
						break;
					}
					int incVal = incPerRound.get(idx) != null ? incPerRound.get(idx).intValue() : 0;
					double relativeGain = (double) incVal / (double) cumulative;
					if (relativeGain >= earlyRelativeGainThreshold) {
						allBelow = false;
					}
				}
				if (allBelow) {
					logger.info(taskName + ": early stopping iterative rounds at round {} (early plateau rule: relative gain < 1.0% for {} consecutive rounds)",
							round, earlyConsecutiveRounds);
					stopReason = "early_plateau_rule_1pct";
					lastRound = round;
					prevPeptideCount = peptideCount;
					break;
				}
			}

			if (round >= 8) {
				DeepIterativeStopRule.StopRuleResult result = computeStopRuleAction.compute(incPerRound);
				lastStopRuleResult = result;
				logger.info(taskName + ": stopRound={} (rounds={}); kChord={}, kPiece={}, kKnee={}, k0={}, tMin={}",
						result.stopRound, round, result.kChord, result.kPiece, result.kKnee, result.k0, result.tMin);
				if (result.stopRound > 0 && result.stopRound == round) {
					logger.info(taskName + ": stopping iterative rounds at round {} (knee+plateau rule, relative gain < 1.5%)",
							round);
					stopReason = "knee_plateau_rule";
					lastRound = round;
					prevPeptideCount = peptideCount;
					break;
				}
			}
			lastRound = round;
			prevPeptideCount = peptideCount;
		}

		return new IterativeRoundLoopResult(lastRound, prevPeptideCount, roundsRecord, peptidesInRoundRecord, incPerRound,
				cumulativeRecord, stopReason, lastStopRuleResult);
	}
}
