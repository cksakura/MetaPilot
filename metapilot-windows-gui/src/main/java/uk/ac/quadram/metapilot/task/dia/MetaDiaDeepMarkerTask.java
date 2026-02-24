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

import org.apache.logging.log4j.Logger;

import javafx.scene.control.ProgressBar;

import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNParquetReader;
import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNPrecursor;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
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
	private static final int MAX_ITERATIVE_ROUNDS = 50;

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

			int peptideCount = runRound(round, markerFile, rankedFile, genomes, libParquetPaths);
			if (peptideCount < 0) {
				return false;
			}

			lastRoundDir = new File(markerFile, "round" + round);
			lastRound = round;
			logger.info(getTaskName() + ": round {}: identified {} peptides (previous round: {})", round, peptideCount, prevPeptideCount);
			if (round > 1 && peptideCount <= prevPeptideCount) {
				logger.info(getTaskName() + ": identification did not increase ({} <= {}), stopping iterative rounds (same as Deep Comprehensive)", peptideCount, prevPeptideCount);
				break;
			}
			prevPeptideCount = peptideCount;

			// For next round: each round folder has roundN-lib.parquet (round1-lib.parquet, round2-lib.parquet, ...)
			File roundLibFile = new File(lastRoundDir, ROUND_LIB_PARQUET_PREFIX + round + ROUND_LIB_PARQUET_SUFFIX);
			if (roundLibFile.exists()) {
				libParquetPaths.add(roundLibFile);
				logger.info(getTaskName() + ": round {}: library parquet for next round: {} ({} bytes)", round, roundLibFile.getAbsolutePath(), roundLibFile.length());
			}
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
			diannPar.setReanalyse(true);
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
	private int runRound(int round, File markerFile, File rankedFile, List<String> genomes, List<File> previousLibParquets) {
		String roundFolder = "round" + round;
		File roundDir = new File(markerFile, roundFolder);
		if (!roundDir.exists()) roundDir.mkdir();

		// --- Combined FASTA ---
		File combinedFasta = new File(roundDir, "combined.fasta");
		if (!writeCombinedFastaFromGenomes(genomes, combinedFasta)) {
			logger.error(getTaskName() + ": round {}: failed to write combined.fasta", round);
			return -1;
		}
		logger.info(getTaskName() + ": round {}: combined.fasta: {} ({} bytes)", round, combinedFasta.getAbsolutePath(), combinedFasta.length());

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
		logger.info(getTaskName() + ": round {}: combined.speclib: {} ({} bytes)", round, combinedSpeclib.getAbsolutePath(), combinedSpeclib.length());

		// --- Libraries for search ---
		String[] libPaths;
		boolean useGenLib = (round == 1);
		File roundLibParquet = new File(roundDir, ROUND_LIB_PARQUET_PREFIX + round + ROUND_LIB_PARQUET_SUFFIX);
		if (round == 1) {
			libPaths = new String[] { combinedSpeclib.getAbsolutePath() };
			logger.info(getTaskName() + ": round {}: search library: combined.speclib (genLib=true → round1-lib.parquet)", round);
		} else {
			if (!roundLibParquet.exists()) {
				diannPar.setThreads(threadCount);
				diaNNTask.addSpeclib2ParquetTask(diannPar, combinedSpeclib.getAbsolutePath(), roundLibParquet.getAbsolutePath());
				diaNNTask.run(12);
			}
			if (!roundLibParquet.exists()) {
				logger.error(getTaskName() + ": round {}: failed to create {}", round, roundLibParquet.getName());
				return -1;
			}
			logger.info(getTaskName() + ": round {}: {} ({} bytes)", round, roundLibParquet.getAbsolutePath(), roundLibParquet.length());

			List<String> allLibPaths = new ArrayList<>();
			for (File f : previousLibParquets) {
				if (f.exists()) allLibPaths.add(f.getAbsolutePath());
			}
			allLibPaths.add(roundLibParquet.getAbsolutePath());
			libPaths = allLibPaths.toArray(new String[0]);
			logger.info(getTaskName() + ": round {}: search libraries ({}): previous round libs + round{}-lib.parquet", round, libPaths.length, round);
		}

		// --- Search ---
		File roundParquet = new File(roundDir, "round" + round + ".parquet");
		diaLibSearchPar.setMbr(true);
		diaLibSearchPar.setDiaFiles(diaFiles);
		diaLibSearchPar.setThreadCount(threadCount);
		diaLibSearchPar.setqValue(0.01);
		diannPar.setThreads(threadCount);
		diannPar.setMatrices(true);
		diannPar.setRelaxed_prot_inf(true);
		diannPar.setReanalyse(true);
		diannPar.setExportDecoy(true);

		if (useGenLib) {
			diaNNTask.addTask(diannPar, diaLibSearchPar, roundParquet.getAbsolutePath(), false, true, true, libPaths[0]);
		} else {
			diaNNTask.addTask(diannPar, diaLibSearchPar, roundParquet.getAbsolutePath(), libPaths);
		}
		diaNNTask.run(diaFiles.length * 12);

		// DIA-NN writes round1lib.parquet when genLib=true; rename to round1-lib.parquet for consistent naming
		if (round == 1) {
			File diannLib = new File(roundDir, "round1lib.parquet");
			if (diannLib.exists() && !roundLibParquet.exists()) {
				if (!diannLib.renameTo(roundLibParquet)) {
					logger.warn(getTaskName() + ": round 1: could not rename {} to {}", diannLib.getName(), roundLibParquet.getName());
				} else {
					logger.info(getTaskName() + ": round 1: renamed {} → {}", diannLib.getName(), roundLibParquet.getName());
				}
			}
		}

		if (!roundParquet.exists()) {
			logger.error(getTaskName() + ": round {}: search failed - {} not created", round, roundParquet.getName());
			return -1;
		}
		logger.info(getTaskName() + ": round {}: result parquet: {} ({} bytes)", round, roundParquet.getAbsolutePath(), roundParquet.length());

		int count = countIdentifiedPeptides(roundParquet);
		logger.info(getTaskName() + ": round {}: unique peptides (q≤0.01, non-decoy): {}", round, count);
		return count;
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
