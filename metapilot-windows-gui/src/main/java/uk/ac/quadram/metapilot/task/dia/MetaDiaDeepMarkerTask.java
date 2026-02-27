package uk.ac.quadram.metapilot.task.dia;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.Arrays;

import org.apache.logging.log4j.Logger;

import javafx.scene.control.ProgressBar;

import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNParquetReader;
import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNPrecursor;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.core.prodb.FastaManager;
import uk.ac.quadram.metapilot.task.dia.par.MetaParameterDiaFocus;
import uk.ac.quadram.metapilot.task.mag.MarkerGenomeSelector;

/**
 * DIA Deep Marker workflow: extends the marker strategy with a deep (two-head) model
 * for detectability and intensity prediction. Uses the same marker-based genome selection
 * as {@link MetaDiaMarkerTask}, with additional training and prediction steps.
 * <p>
 * Overrides {@link #iterativeSearch()} to run the same marker search and genome selection
 * (genome_evidence.tsv, genome_ranked.tsv, genome_selected.tsv) but stops there; does not
 * generate combined.fasta or combined.speclib.
 */
public class MetaDiaDeepMarkerTask extends MetaDiaMarkerTask {

	private static final String TASK_NAME = "DIA Deep Marker analysis";
	/** Number of top genomes from genome_ranked.tsv per round. */
	private static final int GENOMES_PER_ROUND = 5;
	/** Library parquet name in each round folder: round1-lib.parquet, round2-lib.parquet, ... (search result is round1.parquet, round2.parquet, ...). */
	private static final String ROUND_LIB_PARQUET_PREFIX = "round";
	private static final String ROUND_LIB_PARQUET_SUFFIX = "-lib.parquet";
	/** Stop iterative rounds when identification count does not increase (like Deep Comprehensive). */
	private static final int MAX_ITERATIVE_ROUNDS = 20;

	/**
	 * If true, use append-previous-peptides strategy: each round builds one FASTA = [proteins from 5 genomes] + [identified peptides from previous round as single-protein entries]; one library per round.
	 * If false, use runRound(): multiple libraries per round (round1-lib.parquet + round2/combined.speclib.parquet + ...).
	 * Toggle for testing which performs better.
	 */
	private static final boolean USE_APPEND_PREVIOUS_PEPTIDES = true;

	public MetaDiaDeepMarkerTask(MetaParameterDiaFocus focusParams, MetaResourceParameters mrp,
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

		File markerLibFile = magDb.getMarkerHighPepLibFile();
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
		if (!tempParquetFile.exists()) {
			diaLibSearchPar.setMbr(false);
			diaLibSearchPar.setDiaFiles(diaFiles);
			diaLibSearchPar.setThreadCount(threadCount);
			diaLibSearchPar.setqValue(0.01);
			diannPar.setThreads(threadCount);
			diannPar.setMatrices(false);
			diannPar.setRelaxed_prot_inf(true);
			diannPar.setReanalyse(true);
			diannPar.setExportDecoy(true);
			logger.info("Searching marker library against DIA files...");
			logger.info("  Marker library: {}", markerLibFile.getAbsolutePath());
			logger.info("  DIA files count: {}", diaFiles != null ? diaFiles.length : 0);
			diaNNTask.addTask(diannPar, diaLibSearchPar, tempParquetFile.getAbsolutePath(),
					new String[] { markerLibFile.getAbsolutePath() });
			diaNNTask.run(diaFiles.length * 12);
		} else {
			logger.info("Marker library search result already exists in {}", tempParquetFile.getAbsolutePath());
		}

		if (!tempParquetFile.exists()) {
			logger.error(getTaskName() + ": Iterative search failed - parquet file not created");
			System.err.println(new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss").format(new Date()) + "\t" + getTaskName()
					+ ": Iterative search failed - parquet file not created");
			return false;
		}

		// Genome selection only (no combined.fasta / combined.speclib)
		logger.info("Reading PSMs from iterative search results...");
		DiaNNParquetReader parquetReader = new DiaNNParquetReader(tempParquetFile.getAbsolutePath(), diannPar.getCut());
		DiaNNPrecursor[] precursors = parquetReader.getDiaNNPrecursors(0.01);
		logger.info("Found {} precursors with Global.Q.Value < 0.01", precursors != null ? precursors.length : 0);

		if (precursors == null || precursors.length == 0) {
			logger.warn(getTaskName() + ": No precursors found in iterative search, cannot select genomes");
			return false;
		}

		MetaParameterDiaFocus focusParams = (metaPar instanceof MetaParameterDiaFocus) ? (MetaParameterDiaFocus) metaPar : null;
		String strategy = focusParams != null ? focusParams.getGenomeSelectionStrategy() : null;
		String mode = (strategy != null && strategy.length() > 0) ? strategy : MarkerGenomeSelector.KNEE_MODE_FOCUS;
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
	 * Iterative rounds: each round adds next 5 genomes, builds fasta → speclib → (parquet for round>1),
	 * search; stop when identified peptide count does not increase (same as Deep Comprehensive).
	 * Round1 uses genLib so DIA-NN writes round1lib.parquet, which is renamed to round1-lib.parquet; round N>=2 writes roundN-lib.parquet.
	 */
	@Override
	protected boolean iden() {
		File markerFile = getIterativeResultFolder();
		File rankedFile = new File(markerFile, "genome_ranked.tsv");
		if (!rankedFile.exists()) {
			logger.error(getTaskName() + ": genome_ranked.tsv not found: {}", rankedFile.getAbsolutePath());
			return false;
		}

		List<File> libParquetPaths = new ArrayList<>();
		int prevPeptideCount = 0;
		File lastRoundDir = null;
		int lastRound = 0;
		// per-round incremental new peptide counts (inc[t]); used for knee+plateau stop rule
		List<Integer> incPerRound = new ArrayList<>();
		// set of peptides seen in all previous rounds (for computing inc[t] = new_vs_previous)
		HashSet<String> seenPeptides = new HashSet<>();

		for (int round = 1; round <= MAX_ITERATIVE_ROUNDS; round++) {
			updateMessage(getTaskName() + ": round " + round + " library and search");
			logger.info(getTaskName() + ": ========== Round {} ==========", round);

			int offset = (round - 1) * GENOMES_PER_ROUND;
			List<String> genomes = readGenomesFromRanked(rankedFile, offset, GENOMES_PER_ROUND);
			if (genomes.isEmpty()) {
				logger.info(getTaskName() + ": round {}: no genomes at rank {}-{} (end of ranked list), stopping", round, offset + 1, offset + GENOMES_PER_ROUND);
				break;
			}

			int rankStart = offset + 1;
			int rankEnd = offset + genomes.size();
			logger.info(getTaskName() + ": round {}: rank range {}-{}, genomes: {}", round, rankStart, rankEnd, genomes);

			int peptideCount = USE_APPEND_PREVIOUS_PEPTIDES
					? runRoundAppendPreviousPeptides(round, markerFile, genomes, libParquetPaths)
					: runRound(round, markerFile, rankedFile, genomes, libParquetPaths);
			if (peptideCount < 0) {
				return false;
			}

			lastRoundDir = new File(markerFile, "round" + round);
			lastRound = round;

			// Compute incremental new_vs_previous count for this round (inc[round]):
			File roundParquet = new File(lastRoundDir, "round" + round + ".parquet");
			int incThis = countNewPeptidesFromParquet(roundParquet, seenPeptides);
			incPerRound.add(incThis);
			int cumPeptides = seenPeptides.size();

			logger.info(getTaskName() + ": round {}: identified {} peptides in round (previous round: {}), new_vs_previous = {}, cumulative = {}",
					round, peptideCount, prevPeptideCount, incThis, cumPeptides);

			// Deep-comprehensive DIA stop rule (knee + local confirmation) applied online:
			// run rule only once we have enough rounds (>=8) and at least 5 points for local window.
			if (round >= 8) {
				int stopRound = computeStopRoundFromIncrements(incPerRound);
				logger.info(getTaskName() + ": current stopRound suggestion from knee+plateau rule = {} (rounds so far = {})",
						stopRound, round);
				if (stopRound == round) {
					logger.info(getTaskName() + ": stopping iterative rounds at round {} based on knee+local plateau rule", round);
					break;
				}
			}

			prevPeptideCount = peptideCount;
		}

		if (lastRoundDir == null) {
			logger.error(getTaskName() + ": no round completed");
			return false;
		}

		logger.info(getTaskName() + ": iterative rounds finished at round {} ({} peptides); running combined search with all round libraries...", lastRound, prevPeptideCount);

		// Combined folder: search all raw files with all {round}.lib.parquet libraries → combined.parquet (then unified post-analysis)
		File combinedDir = new File(metaPar.getDbSearchResultFile(), "combined");
		if (!combinedDir.exists()) combinedDir.mkdirs();
		File combinedParquet = new File(combinedDir, "combined.parquet");
		File combinedPepMatrix = new File(combinedDir, "combined.pr_matrix.tsv");

		if (libParquetPaths.isEmpty()) {
			logger.error(getTaskName() + ": no library parquets to run combined search");
			return false;
		}

		if (!combinedParquet.exists() || !combinedPepMatrix.exists()) {
			String[] allLibPaths = new String[libParquetPaths.size()];
			for (int i = 0; i < libParquetPaths.size(); i++) {
				allLibPaths[i] = libParquetPaths.get(i).getAbsolutePath();
			}
			logger.info(getTaskName() + ": combined search: {} libraries → {}", allLibPaths.length, combinedParquet.getAbsolutePath());
			for (int i = 0; i < allLibPaths.length; i++) {
				logger.info(getTaskName() + ":   combined lib {}: {}", i + 1, allLibPaths[i]);
			}

			diaLibSearchPar.setMbr(true);
			diaLibSearchPar.setDiaFiles(diaFiles);
			diaLibSearchPar.setThreadCount(threadCount);
			diaLibSearchPar.setqValue(0.01);
			diannPar.setThreads(threadCount);
			diannPar.setMatrices(true);
			diannPar.setRelaxed_prot_inf(true);
			diannPar.setReanalyse(false);
			diannPar.setExportDecoy(true);
			diaNNTask.addTask(diannPar, diaLibSearchPar, combinedParquet.getAbsolutePath(), allLibPaths);
			diaNNTask.run(diaFiles.length * 12);
		} else {
			logger.info(getTaskName() + ": combined result already exists, skipping combined search");
		}

		if (!combinedParquet.exists()) {
			logger.error(getTaskName() + ": combined search failed - combined.parquet not created");
			return false;
		}

		combinedFile = combinedDir;
		this.quan_pep_file = new File(combinedDir, "combined.pr_matrix.tsv");
		this.quan_pro_file = new File(combinedDir, "combined.pg_matrix.tsv");
		logger.info(getTaskName() + ": combined search completed: {} ({} bytes); unified post-analysis will use combined folder", combinedParquet.getAbsolutePath(), combinedParquet.length());
		System.out.println(new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss").format(new Date()) + "\t" + getTaskName() + ": identification finished (combined)");
		return true;
	}

	/**
	 * Run one round: genomes → fasta → speclib → (parquet if round>1) → search. Returns identified peptide count (q≤0.01, non-decoy) or -1 on failure.
	 */
	private int runRound(int round, File markerFile, File rankedFile, List<String> genomes,
			List<File> previousLibParquets) {
		String roundFolder = "round" + round;
		File roundDir = new File(markerFile, roundFolder);
		if (!roundDir.exists())
			roundDir.mkdir();

		// --- Combined FASTA ---
		File combinedFasta = new File(roundDir, "combined.fasta");
		if (!writeCombinedFastaFromGenomes(genomes, combinedFasta)) {
			logger.error(getTaskName() + ": round {}: failed to write combined.fasta", round);
			return -1;
		}
		logger.info(getTaskName() + ": round {}: combined.fasta: {} ({} bytes)", round, combinedFasta.getAbsolutePath(),
				combinedFasta.length());

		// --- Combined speclib ---
		File combinedSpeclib = new File(roundDir, "combined.speclib");
		if (!combinedSpeclib.exists()) {
			File altLib = new File(roundDir, "combined.predicted.speclib");
			if (altLib.exists()) {
				altLib.renameTo(combinedSpeclib);
			} else {
				diannPar.setThreads(threadCount);
				diaNNTask.addTask(diannPar, combinedFasta.getAbsolutePath(), combinedSpeclib.getAbsolutePath(), false);
				diaNNTask.run(12);
				if (!combinedSpeclib.exists() && altLib.exists()) {
					altLib.renameTo(combinedSpeclib);
				}
			}
		}
		if (!combinedSpeclib.exists()) {
			logger.error(getTaskName() + ": round {}: failed to generate combined.speclib", round);
			return -1;
		}
		logger.info(getTaskName() + ": round {}: combined.speclib: {} ({} bytes)", round,
				combinedSpeclib.getAbsolutePath(), combinedSpeclib.length());
		
		String[] libPaths = new String[] { combinedSpeclib.getAbsolutePath() };

/*
		File combinedParquetlib = new File(roundDir, "combined.speclib.parquet");
		if (!combinedParquetlib.exists()) {
			diannPar.setThreads(threadCount);
			diaNNTask.addSpeclib2ParquetTask(diannPar, combinedSpeclib.getAbsolutePath(),
					combinedParquetlib.getAbsolutePath());
			diaNNTask.run(12);
		}
		if (!combinedParquetlib.exists()) {
			logger.error(getTaskName() + ": round {}: failed to generate combined.speclib.parquet", round);
			return -1;
		}
		logger.info(getTaskName() + ": round {}: combined.speclib.parquet: {} ({} bytes)", round,
				combinedParquetlib.getAbsolutePath(), combinedParquetlib.length());

		// --- Libraries for search ---
		String[] libPaths;
		if (round == 1) {
			libPaths = new String[] { combinedParquetlib.getAbsolutePath() };
			logger.info(getTaskName() + ": round {}: search library: {}", round, combinedParquetlib.getAbsolutePath());
		} else {
			List<String> allLibPaths = new ArrayList<>();
			allLibPaths.add(combinedParquetlib.getAbsolutePath());
			for (File f : previousLibParquets) {
				if (f.exists())
					allLibPaths.add(f.getAbsolutePath());
			}
			libPaths = allLibPaths.toArray(new String[0]);
			logger.info(getTaskName() + ": round {}: search libraries ({}): previous round libs + round{}-lib.parquet",
					round, libPaths.length, round);
		}
*/
		// --- Search ---
		File roundParquet = new File(roundDir, "round" + round + ".parquet");
		File roundParquetLib = new File(roundDir, "round" + round + ".parquet.parquet");
		
		if (!roundParquet.exists() || !roundParquetLib.exists()) {
			diaLibSearchPar.setMbr(false);
			diaLibSearchPar.setDiaFiles(diaFiles);
			diaLibSearchPar.setThreadCount(threadCount);
			diaLibSearchPar.setqValue(0.01);
			diannPar.setThreads(threadCount);
			diannPar.setMatrices(false);
			diannPar.setRelaxed_prot_inf(true);
			diannPar.setReanalyse(true);
			diannPar.setExportDecoy(false);
//			diaNNTask.addTask(diannPar, diaLibSearchPar, roundParquet.getAbsolutePath(), false, false, true, libPaths);
			diaNNTask.addNoQuantLibTask(diannPar, diaLibSearchPar, roundParquet.getAbsolutePath(), libPaths);
			diaNNTask.run(diaFiles.length * 12);
		}

		if (!roundParquet.exists()) {
			logger.error(getTaskName() + ": round {}: search failed - {} not created", round, roundParquet.getName());
			return -1;
		}
		if (!roundParquetLib.exists()) {
			logger.error(getTaskName() + ": round {}: search failed - {} not created", round,
					roundParquetLib.getName());
			return -1;
		}
		logger.info(getTaskName() + ": round {}: result parquet: {} ({} bytes)", round, roundParquet.getAbsolutePath(),
				roundParquet.length());

		previousLibParquets.add(roundParquetLib);

		int count = countIdentifiedPeptides(roundParquet);
		logger.info(getTaskName() + ": round {}: unique peptides (q≤0.01, non-decoy): {}", round, count);
		return count;
	}

	/**
	 * Alternative round strategy: one library per round. FASTA = [proteins from 5 genomes] + [identified peptides from previous round as single-protein entries; DIA-NN accepts multiple entries with same protein name].
	 * Round 1: genomes only → speclib → search with genLib → round1-lib.parquet. Round 2+: genomes + peptides from round N-1 parquet → speclib → parquet → search with that single lib.
	 * Returns identified peptide count or -1 on failure.
	 */
	private int runRoundAppendPreviousPeptides(int round, File markerFile, List<String> genomes,
			List<File> libParquetPathsCollector) {
		String roundFolder = "round" + round;
		File roundDir = new File(markerFile, roundFolder);
		if (!roundDir.exists()) roundDir.mkdir();

		File combinedFasta = new File(roundDir, "combined.fasta");
		if (round == 1) {
			if (!writeCombinedFastaFromGenomes(genomes, combinedFasta)) {
				logger.error(getTaskName() + ": round {}: failed to write combined.fasta", round);
				return -1;
			}
			logger.info(getTaskName() + ": round {} (append strategy): combined.fasta from {} genomes only", round, genomes.size());
		} else {
			File prevFolder = new File(markerFile, "round" + (round - 1));
			File prevParquet = new File(prevFolder, "round" + (round - 1) + ".parquet");
			if (!prevParquet.exists()) {
				logger.error(getTaskName() + ": round {}: previous round parquet not found: {}", round, prevParquet.getAbsolutePath());
				return -1;
			}
			if (!writeCombinedFastaFromGenomesAndPreviousPeptides(genomes, prevParquet, combinedFasta)) {
				logger.error(getTaskName() + ": round {}: failed to write combined.fasta (genomes + previous peptides)", round);
				return -1;
			}
			logger.info(getTaskName() + ": round {} (append strategy): combined.fasta = {} genomes + peptides from previous round", round, genomes.size());
		}
		logger.info(getTaskName() + ": round {}: combined.fasta: {} ({} bytes)", round, combinedFasta.getAbsolutePath(), combinedFasta.length());

		File combinedSpeclib = new File(roundDir, "combined.speclib");
		if (!combinedSpeclib.exists()) {
			File altLib = new File(roundDir, "combined.predicted.speclib");
			if (altLib.exists()) altLib.renameTo(combinedSpeclib);
			else {
				diannPar.setThreads(threadCount);
				diaNNTask.addTask(diannPar, combinedFasta.getAbsolutePath(), combinedSpeclib.getAbsolutePath(), false);
				diaNNTask.run(12);
				if (!combinedSpeclib.exists() && altLib.exists()) altLib.renameTo(combinedSpeclib);
			}
		}
		if (!combinedSpeclib.exists()) {
			logger.error(getTaskName() + ": round {}: failed to generate combined.speclib", round);
			return -1;
		}
		logger.info(getTaskName() + ": round {}: combined.speclib: {} ({} bytes)", round, combinedSpeclib.getAbsolutePath(), combinedSpeclib.length());

		File roundParquet = new File(roundDir, "round" + round + ".parquet");

		if(!roundParquet.exists()) {
			diaLibSearchPar.setMbr(false);
			diaLibSearchPar.setDiaFiles(diaFiles);
			diaLibSearchPar.setThreadCount(threadCount);
			diaLibSearchPar.setqValue(0.01);
			diannPar.setThreads(threadCount);
			diannPar.setMatrices(false);
			diannPar.setRelaxed_prot_inf(true);
			diannPar.setReanalyse(true);
			diannPar.setExportDecoy(false);
//			diaNNTask.addTask(diannPar, diaLibSearchPar, roundParquet.getAbsolutePath(), false, true, true, combinedSpeclib.getAbsolutePath());
			diaNNTask.addNoQuantLibTask(diannPar, diaLibSearchPar, roundParquet.getAbsolutePath(), combinedSpeclib.getAbsolutePath());
			diaNNTask.run(diaFiles.length * 12);
		}
		

		if (!roundParquet.exists()) {
			logger.error(getTaskName() + ": round {}: search failed - {} not created", round, roundParquet.getName());
			return -1;
		}

		int count = countIdentifiedPeptides(roundParquet);
		logger.info(getTaskName() + ": round {}: unique peptides (q≤0.01, non-decoy): {}", round, count);
		return count;
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
			DiaNNParquetReader reader = new DiaNNParquetReader(previousRoundParquet.getAbsolutePath(), cut);
			DiaNNPrecursor[] precursors = reader.getDiaNNPrecursors(0.01);
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

	/**
	 * Count incremental new peptides for a round parquet (q≤0.01, non-decoy) versus all previously seen peptides.
	 * Updates seenPeptides with newly observed sequences.
	 */
	private int countNewPeptidesFromParquet(File parquetFile, HashSet<String> seenPeptides) {
		if (parquetFile == null || !parquetFile.exists()) {
			logger.warn(getTaskName() + ": countNewPeptidesFromParquet: parquet file not found: {}", parquetFile);
			return 0;
		}
		String cut = diannPar != null && diannPar.getCut() != null ? diannPar.getCut() : "KR";
		int newCount = 0;
		try {
			DiaNNParquetReader reader = new DiaNNParquetReader(parquetFile.getAbsolutePath(), cut);
			DiaNNPrecursor[] precursors = reader.getDiaNNPrecursors(0.01);
			if (precursors != null) {
				for (DiaNNPrecursor p : precursors) {
					if (p == null || p.getIsDecoy()) continue;
					String seq = p.getSeqString();
					if (seq == null || seq.isEmpty()) continue;
					if (!seenPeptides.contains(seq)) {
						seenPeptides.add(seq);
						newCount++;
					}
				}
			}
		} catch (Exception e) {
			logger.warn(getTaskName() + ": failed to count new peptides from {}", parquetFile.getAbsolutePath(), e);
		}
		return newCount;
	}

	/**
	 * Deep-comprehensive iterative DIA stop rule (knee + local confirmation) on per-round increments inc[t].
	 * inc is 0-based list of length n, representing inc[1..n]. Returns stopRound in 1..n.
	 */
	private int computeStopRoundFromIncrements(List<Integer> inc) {
		if (inc == null) return 0;
		int n = inc.size();
		if (n == 0) return 0;
		if (n < 8) {
			// For very short series, run all available rounds
			return n;
		}

		// Build cumulative curve cum[t] from inc[t]
		double[] cum = new double[n];
		cum[0] = inc.get(0);
		for (int i = 1; i < n; i++) {
			cum[i] = cum[i - 1] + inc.get(i);
		}

		// Step 1: global knee from cumulative curve (1-based indices)
		int kChord = computeKneeIndex1D(cum);
		int kPiece = computePiecewiseIndex1D(cum);
		int kKnee = computeKneedleIndex1D(cum);
		int k0 = medianOfThree(kChord, kPiece, kKnee);
		if (k0 < 1) k0 = 1;
		if (k0 > n) k0 = n;

		// Step 2: local confirmation of plateau on inc[t]
		int tMin = Math.max(8, k0);

		// Early peak scale from first up to 5 rounds
		int earlyEnd = Math.min(5, n);
		int peak = 0;
		for (int i = 0; i < earlyEnd; i++) {
			peak = Math.max(peak, inc.get(i));
		}
		double T = Math.max(500.0, 0.15 * (double) peak);

		int startT = Math.max(tMin, 5); // need at least 5 points for sum5
		for (int t = startT; t <= n; t++) {
			int i = t - 1; // 0-based
			int a = inc.get(i);
			int b = inc.get(i - 1);
			int c = inc.get(i - 2);
			int med3 = medianOfThree(a, b, c);

			int sum5 = a + b + c + inc.get(i - 3) + inc.get(i - 4);

			if (med3 < T && sum5 < 5.0 * T) {
				return t;
			}
		}
		return n;
	}

	// --- 1D knee helpers (duplicated from MarkerGenomeSelector for local use) ---

	private static int computeKneeIndex1D(double[] y) {
		if (y == null || y.length == 0) return 0;
		if (y.length <= 2) return y.length;

		int n = y.length;
		double x1 = 0.0, y1 = y[0];
		double x2 = 1.0, y2 = y[n - 1];
		double denom = Math.sqrt((y2 - y1) * (y2 - y1) + (x2 - x1) * (x2 - x1));
		if (denom == 0.0) return n;

		double maxDist = -1.0;
		int bestIndex = 1;
		for (int i = 0; i < n; i++) {
			double x = (n == 1) ? 0.0 : ((double) i / (double) (n - 1));
			double yi = y[i];
			double dist = Math.abs((y2 - y1) * x - (x2 - x1) * yi + x2 * y1 - y2 * x1) / denom;
			if (dist > maxDist) {
				maxDist = dist;
				bestIndex = i + 1; // 1-based
			}
		}
		return bestIndex;
	}

	private static int computeKneedleIndex1D(double[] y) {
		if (y == null || y.length == 0) return 0;
		if (y.length == 1) return 1;
		double min = Arrays.stream(y).min().orElse(0.0);
		double max = Arrays.stream(y).max().orElse(0.0);
		double range = max - min;
		double bestDiff = -Double.MAX_VALUE;
		int bestIdx = 1;
		for (int i = 0; i < y.length; i++) {
			double x = (y.length == 1) ? 0.0 : ((double) i / (double) (y.length - 1));
			double yn = (range == 0.0) ? 0.0 : ((y[i] - min) / range);
			double diff = yn - (1.0 - x);
			if (diff > bestDiff) {
				bestDiff = diff;
				bestIdx = i + 1;
			}
		}
		return bestIdx;
	}

	private static int computePiecewiseIndex1D(double[] y) {
		if (y == null || y.length == 0) return 0;
		if (y.length <= 3) return y.length;
		double bestSse = Double.MAX_VALUE;
		int bestK = 2;
		for (int k = 2; k <= y.length - 1; k++) {
			double sse1 = linearFitSse1D(y, 0, k - 1);
			double sse2 = linearFitSse1D(y, k - 1, y.length - 1);
			double sse = sse1 + sse2;
			if (sse < bestSse) {
				bestSse = sse;
				bestK = k;
			}
		}
		return bestK;
	}

	private static double linearFitSse1D(double[] y, int start, int end) {
		int n = end - start + 1;
		if (n <= 1) return 0.0;
		double sumX = 0.0, sumY = 0.0, sumXY = 0.0, sumXX = 0.0;
		for (int i = start; i <= end; i++) {
			double x = i - start;
			sumX += x;
			sumY += y[i];
			sumXY += x * y[i];
			sumXX += x * x;
		}
		double denom = n * sumXX - sumX * sumX;
		double b = denom == 0.0 ? 0.0 : (n * sumXY - sumX * sumY) / denom;
		double a = (sumY - b * sumX) / n;
		double sse = 0.0;
		for (int i = start; i <= end; i++) {
			double x = i - start;
			double fit = a + b * x;
			double err = y[i] - fit;
			sse += err * err;
		}
		return sse;
	}

	private static int medianOfThree(int a, int b, int c) {
		int[] v = new int[] { a, b, c };
		Arrays.sort(v);
		return v[1];
	}

	/** Count unique peptide sequences from parquet (q≤0.01, exclude decoys). */
	private int countIdentifiedPeptides(File parquetFile) {
		String cut = diannPar != null && diannPar.getCut() != null ? diannPar.getCut() : "KR";
		Set<String> peptides = new HashSet<>();
		try {
			DiaNNParquetReader reader = new DiaNNParquetReader(parquetFile.getAbsolutePath(), cut);
			DiaNNPrecursor[] precursors = reader.getDiaNNPrecursors(0.01);
			if (precursors != null) {
				for (DiaNNPrecursor p : precursors) {
					if (p == null || p.getIsDecoy()) continue;
					String seq = p.getSeqString();
					if (seq != null && !seq.isEmpty()) peptides.add(seq);
				}
			}
		} catch (Exception e) {
			logger.warn(getTaskName() + ": failed to count peptides from {}", parquetFile.getAbsolutePath(), e);
		}
		return peptides.size();
	}

	/** Read genome IDs from genome_ranked.tsv (column 1 = Genome), skip first `offset` data rows, then take up to `count`. */
	private List<String> readGenomesFromRanked(File rankedTsv, int offset, int count) {
		List<String> result = new ArrayList<>();
		try (BufferedReader r = new BufferedReader(new FileReader(rankedTsv))) {
			if (r.readLine() == null) return result;
			int skipped = 0;
			String line;
			while ((line = r.readLine()) != null && result.size() < count) {
				line = line.trim();
				if (line.isEmpty() || line.startsWith("#")) continue;
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
			logger.warn(getTaskName() + ": error reading ranked file {}", rankedTsv.getAbsolutePath(), e);
		}
		return result;
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
