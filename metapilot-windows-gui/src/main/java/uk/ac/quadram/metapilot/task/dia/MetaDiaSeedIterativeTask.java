package uk.ac.quadram.metapilot.task.dia;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

import javafx.scene.control.ProgressBar;

import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNParquetReader;
import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNPrecursor;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.task.dia.par.MetaParameterDia;
import uk.ac.quadram.metapilot.task.dia.par.MetaParameterDiaMarker;
import uk.ac.quadram.metapilot.task.mag.MagDbItem.MagPepLib;
import uk.ac.quadram.metapilot.task.mag.MarkerGenomeSelector;

/**
 * DIA Seed iterative workflow: user-provided library → initial search → set protein info from MAG predicted TSVs
 * → MarkerGenomeSelector (Focus or Comprehensive) → combined FASTA/library → identification.
 * Precursors from the user library search do not have genome/protein info; protein mapping is done via
 * {@link MarkerGenomeSelector#setProteinInfoFromPredictedTsVs(DiaNNPrecursor[], File)} before genome selection.
 * Library type (microbiome vs standalone, from {@link MagPepLib#isMatchToMag()}):
 * <ul>
 *   <li>All standalone: workflow cannot run → fallback to direct search (iden only, no export).</li>
 *   <li>All microbiome: use full iterative workflow with all libraries.</li>
 *   <li>Mixed: use only microbiome libraries for the initial search and add --out-lib (library generation).</li>
 * </ul>
 * Report export uses the same logic as {@link MetaDiaMarkerTask} (genome/protein scores, peptides, etc.), with
 * combined FASTA taken from {@link #getIterativeResultFolder()} (seed_iterative).
 */
public class MetaDiaSeedIterativeTask extends MetaDiaMarkerTask {

	private static final String TASK_NAME = "DIA Seed iterative";
	private static final String SEED_ITERATIVE_FOLDER = "seed";
	private static final SimpleDateFormat format = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");

	protected File combinedFile;

	/** When true, iterative was skipped and direct iden was run (all standalone); no export. */
	protected boolean fallbackDirectSearch;
	/** Library paths to use for the initial DIA-NN search (all libs or microbiome-only when mixed). */
	protected String[] initialSearchLibPaths;
	/** When true, add --out-lib / library generation to the initial search (mixed case). */
	protected boolean addOutLibForInitialSearch;
	/** When mixed, standalone library paths for the final combined search (null otherwise). */
	protected String[] standaloneSearchLibPaths;

	/** True when the initial (first) search fell back to per-file search. */
	protected boolean initialSearchSeparate;
	/** Subset of diaFiles that produced results (excludes damaged files). Set when initialSearchSeparate is true. */
	protected String[] goodDiaFiles;

	public MetaDiaSeedIterativeTask(MetaParameterDiaMarker markerParams, MetaResourceParameters mrp, ProgressBar bar1,
			ProgressBar bar2) {
		super(markerParams, mrp, bar1, bar2);

		if (diannPar != null) {
			diannPar.setM_excision(markerParams.isMExcision());
			diannPar.setCarbami_C(markerParams.isCarbamiC());
			diannPar.setOxidation_M(markerParams.isOxidationM());
			diannPar.setAcety_proN(markerParams.isAcetyProN());
			diannPar.setPhospho_STY(markerParams.isPhosphoSTY());
			diannPar.setK_GG(markerParams.isKGG());
		}
		// Initialize result folder for seed iterative workflows if not already set:
		// seed_focus / seed_comprehensive, based on genome selection strategy.
		initResultFolderIfNeeded();

		logger.info("DIA Seed iterative workflow ({} mode) initialized",
				markerParams.getMarkerMode());
	}

	@Override
	protected String getTaskName() {
		return TASK_NAME;
	}

	@Override
	protected File getIterativeResultFolder() {
		return new File(metaPar.getDbSearchResultFile(), SEED_ITERATIVE_FOLDER);
	}

	@Override
	protected Boolean call() throws Exception {
		updateMessage(getTaskName());
		updateProgress(0, 100);
		logger.info(getTaskName() + ": started");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": started");

		updateProgress(5, 100);
		updateSecondaryProgress(ProgressBar.INDETERMINATE_PROGRESS);
		if (!copyReportTemplatesInBackground()) {
			logger.error(getTaskName() + ": Failed to copy report templates");
			return false;
		}

		logger.info("Step 1: Starting spectra format conversion...");
		updateMessage(getTaskName() + ": Converting spectra format...");
		updateProgress(10, 100);
		boolean convert = spectraConvert();
		if (!convert) {
			logger.error(getTaskName() + ": error in spectra format conversion");
			updateMessage(getTaskName() + ": Failed - Spectra conversion error");
			updateSecondaryProgress(0.0);
			return false;
		}

		logger.info("Step 2: Classifying libraries (microbiome vs standalone)...");
		updateProgress(18, 100);
		updateMessage(getTaskName() + ": Checking library types...");
		if (!classifyLibraries()) {
			logger.error(getTaskName() + ": library classification failed");
			updateMessage(getTaskName() + ": Failed - Library classification error");
			updateSecondaryProgress(0.0);
			return false;
		}

		if (fallbackDirectSearch) {
			logger.info("All libraries are standalone – fallback to direct search (identification only, no export).");
			updateProgress(25, 100);
			updateMessage(getTaskName() + ": Running direct identification...");
			updateSecondaryProgress(ProgressBar.INDETERMINATE_PROGRESS);
			boolean directOk = directIden();
			if (!directOk) {
				logger.error(getTaskName() + ": error in direct identification");
				updateMessage(getTaskName() + ": Failed - Identification error");
				updateSecondaryProgress(0.0);
				return false;
			}
			updateProgress(100, 100);
			updateMessage(getTaskName() + ": Direct identification completed (no report – standalone libraries).");
			updateSecondaryProgress(0.0);
			logger.info(getTaskName() + ": finished – direct search fallback, no export.");
			System.out.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": finished – direct search (standalone libraries, no export).");
			return true;
		}

		logger.info("Step 3: Validating user library...");
		updateProgress(20, 100);
		updateMessage(getTaskName() + ": Validating library...");
		boolean createOk = createSSLibrary();
		if (!createOk) {
			logger.error(getTaskName() + ": error in library validation");
			updateMessage(getTaskName() + ": Failed - Library validation error");
			updateSecondaryProgress(0.0);
			return false;
		}

		logger.info("Step 4: Starting iterative search (user library → protein mapping → genome selection)...");
		updateProgress(25, 100);
		updateMessage(getTaskName() + ": Performing iterative search...");
		updateSecondaryProgress(ProgressBar.INDETERMINATE_PROGRESS);
		boolean iterative = iterativeSearch();
		if (!iterative) {
			logger.error(getTaskName() + ": error in iterative search");
			updateMessage(getTaskName() + ": Failed - Iterative search error");
			updateSecondaryProgress(0.0);
			return false;
		}
		updateProgress(50, 100);

		logger.info("Step 5: Identification using combined library...");
		updateMessage(getTaskName() + ": Running identification...");
		updateProgress(55, 100);
		updateSecondaryProgress(ProgressBar.INDETERMINATE_PROGRESS);
		boolean iden = iden();
		if (!iden) {
			logger.error(getTaskName() + ": error in identification");
			updateMessage(getTaskName() + ": Failed - Identification error");
			updateSecondaryProgress(0.0);
			return false;
		}

		logger.info("Step 6: Starting report generation...");
		updateProgress(85, 100);
		updateMessage(getTaskName() + ": Generating reports...");
		updateSecondaryProgress(ProgressBar.INDETERMINATE_PROGRESS);
		boolean exportReport = exportReport();
		if (!exportReport) {
			logger.error(getTaskName() + ": error in report generation");
			updateMessage(getTaskName() + ": Failed - Report generation error");
			updateSecondaryProgress(0.0);
			return false;
		}

		updateProgress(100, 100);
		updateMessage(getTaskName() + ": Completed successfully!");
		updateSecondaryProgress(0.0);
		logger.info(getTaskName() + ": finished - all workflow steps completed successfully");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": finished successfully");
		return true;
	}

	/**
	 * Resolves selected library names/paths to absolute file paths. Used by {@link #classifyLibraries()}
	 * and by subclasses that need all lib paths (e.g. Deep Comprehensive).
	 * @return resolved paths, or null if metaPar/MAG/libraries are missing or invalid
	 */
	protected String[] resolveLibraryPaths() {
		if (metaPar == null || !(metaPar instanceof MetaParameterDia) || magDb == null) {
			return null;
		}
		String[] usedLibs = ((MetaParameterDia) metaPar).getLibrary();
		if (usedLibs == null || usedLibs.length == 0) {
			return null;
		}
		String[] libPath = new String[usedLibs.length];
		for (int i = 0; i < usedLibs.length; i++) {
			File directLibFile = new File(usedLibs[i]);
			if (directLibFile.exists()) {
				libPath[i] = directLibFile.getAbsolutePath();
			} else {
				File chosen = magDb.resolveSpectralLibraryFile(usedLibs[i]);
				if (chosen == null || !chosen.exists()) {
					logger.error(getTaskName() + ": no library file (.speclib/.tsv/.parquet) found for library '{}'", usedLibs[i]);
					return null;
				}
				libPath[i] = chosen.getAbsolutePath();
			}
		}
		return libPath;
	}

	/**
	 * Classifies selected libraries into microbiome vs standalone (MagPepLib.isMatchToMag).
	 * Sets {@link #fallbackDirectSearch}, {@link #initialSearchLibPaths}, {@link #addOutLibForInitialSearch}.
	 * @return true if classification succeeded and at least one library is present
	 */
	protected boolean classifyLibraries() {
		String[] libPath = resolveLibraryPaths();
		if (libPath == null || libPath.length == 0) {
			logger.error(getTaskName() + ": MetaParameterDia and MAG database required, or no libraries selected");
			return false;
		}
		String[] usedLibs = ((MetaParameterDia) metaPar).getLibrary();
		MagPepLib[] allLibs = magDb.getMagPepLibs();
		boolean[] matchToMag = new boolean[usedLibs.length];
		for (int i = 0; i < usedLibs.length; i++) {
			File directLibFile = new File(usedLibs[i]);
			if (directLibFile.exists()) {
				matchToMag[i] = false;
			} else {
				boolean micro = false;
				for (MagPepLib lib : allLibs) {
					if (lib.getLibName().equals(usedLibs[i])) {
						micro = lib.isMatchToMag();
						break;
					}
				}
				matchToMag[i] = micro;
			}
		}

		int microCount = 0;
		for (boolean b : matchToMag) {
			if (b) microCount++;
		}
		boolean allStandalone = (microCount == 0);
		boolean allMicrobiome = (microCount == usedLibs.length);
		boolean mixed = !allStandalone && !allMicrobiome;

		fallbackDirectSearch = allStandalone;
		addOutLibForInitialSearch = mixed;
		if (allStandalone) {
			initialSearchLibPaths = libPath;
			logger.info(getTaskName() + ": all libraries standalone – will fallback to direct search.");
		} else if (allMicrobiome) {
			initialSearchLibPaths = libPath;
			logger.info(getTaskName() + ": all libraries microbiome – full iterative workflow.");
		} else {
			List<String> microPaths = new ArrayList<>();
			List<String> standalonePaths = new ArrayList<>();
			for (int i = 0; i < usedLibs.length; i++) {
				if (matchToMag[i]) microPaths.add(libPath[i]);
				else standalonePaths.add(libPath[i]);
			}
			initialSearchLibPaths = microPaths.toArray(new String[0]);
			standaloneSearchLibPaths = standalonePaths.toArray(new String[0]);
			logger.info(getTaskName() + ": mixed libraries – using microbiome libraries only with --out-lib.");
		}
		return true;
	}

	/**
	 * Direct identification with all selected libraries (no iterative, no export).
	 * Used when all libraries are standalone. Writes results to combined/.
	 */
	private boolean directIden() {
		combinedFile = new File(metaPar.getDbSearchResultFile(), "combined");
		if (!combinedFile.exists()) {
			combinedFile.mkdir();
		}
		// Required for exportReport() / getExportQuanPepFile()
		this.quan_pep_file = new File(combinedFile, "combined.pr_matrix.tsv");
		this.quan_pro_file = new File(combinedFile, "combined.pg_matrix.tsv");
		File combineParquetFile = new File(combinedFile, "combined.parquet");

		diaLibSearchPar.setMbr(true);
		diaLibSearchPar.setDiaFiles(diaFiles);
		diaLibSearchPar.setThreadCount(threadCount);
		diaLibSearchPar.setqValue(0.01);
		diannPar.setThreads(threadCount);
		diannPar.setMatrices(true);
		diannPar.setRelaxed_prot_inf(true);
		diannPar.setReanalyse(true);
		diannPar.setExportDecoy(true);

		diaNNTask.addTask(diannPar, diaLibSearchPar, combineParquetFile.getAbsolutePath(), initialSearchLibPaths);
		diaNNTask.run(diaFiles.length * 12);

		if (!combineParquetFile.exists()) {
			logger.error(getTaskName() + ": Direct identification failed - parquet not created");
			return false;
		}
		logger.info("Direct identification completed.");
		return true;
	}

	@Override
	protected boolean createSSLibrary() {
		// No marker library to build; user library is already selected
		logger.info("User library mode - skipping sample-specific library creation");
		return true;
	}

	/**
	 * Runs the first search with user-selected (or microbiome-only) libraries and writes
	 * result to {@link #getIterativeResultFolder()}/seed.parquet when possible. If the
	 * combined search fails due to one or more damaged raw files, a per-file fallback is
	 * executed and per_file_*.parquet files are produced instead; subclasses can consume
	 * those directly (e.g. Deep Comprehensive).
	 * @return true if either seed.parquet exists or at least one non-empty per_file_*.parquet
	 *         exists under the iterative result folder
	 */
	protected boolean runFirstSearchUserLibrary() {
		File iterativeFolder = getIterativeResultFolder();
		if (!iterativeFolder.exists()) {
			iterativeFolder.mkdirs();
		}
		File tempParquetFile = new File(iterativeFolder, "seed.parquet");
		if (initialSearchLibPaths == null || initialSearchLibPaths.length == 0) {
			logger.error(getTaskName() + ": No library paths for first search (classification may have failed)");
			return false;
		}
		// Check if per-file parquets already exist from a previous separate search.
		boolean existingPerFile = false;
		if (!tempParquetFile.exists() && diaFiles != null && diaFiles.length > 0) {
			for (int i = 0; i < diaFiles.length; i++) {
				File pf = new File(iterativeFolder, "per_file_" + i + ".parquet");
				if (pf.exists() && pf.length() > 0) {
					existingPerFile = true;
					break;
				}
			}
		}

		if (tempParquetFile.exists()) {
			initialSearchSeparate = false;
			goodDiaFiles = diaFiles;
			return true;
		} else if (existingPerFile) {
			List<String> goodList = new ArrayList<>();
			for (int i = 0; i < diaFiles.length; i++) {
				File pf = new File(iterativeFolder, "per_file_" + i + ".parquet");
				if (pf.exists() && pf.length() > 0) {
					goodList.add(diaFiles[i]);
				}
			}
			initialSearchSeparate = true;
			goodDiaFiles = goodList.toArray(new String[0]);
			int maxConcurrent = DeepIterativeCommonSupport.computeMaxConcurrentTasks();
			logger.info(getTaskName()
					+ ": per-file parquets already exist from previous run; {} good files; rounds will use per-file search (max concurrent: {})",
					goodDiaFiles.length, maxConcurrent);
			return true;
		}

		// No existing results — run combined search.
		logger.info("Searching {} library against DIA files{}...",
				addOutLibForInitialSearch ? "microbiome" : "user",
				addOutLibForInitialSearch ? " (with --out-lib)" : "");
		diaLibSearchPar.setMbr(false);
		diaLibSearchPar.setDiaFiles(diaFiles);
		diaLibSearchPar.setThreadCount(threadCount);
		diaLibSearchPar.setqValue(0.01);
		diannPar.setThreads(threadCount);
		diannPar.setMatrices(false);
		diannPar.setRelaxed_prot_inf(true);
		diannPar.setReanalyse(true);
		diannPar.setExportDecoy(true);

		if (addOutLibForInitialSearch) {
			diannPar.setMatrices(false);
			diaNNTask.addTask(diannPar, diaLibSearchPar, tempParquetFile.getAbsolutePath(),
					false, true, true, initialSearchLibPaths);
		} else {
			diaNNTask.addNoQuantLibTask(diannPar, diaLibSearchPar, tempParquetFile.getAbsolutePath(),
					initialSearchLibPaths);
		}
		diaNNTask.run(diaFiles.length * 12);

		if (tempParquetFile.exists()) {
			initialSearchSeparate = false;
			goodDiaFiles = diaFiles;
			return true;
		}

		// Combined search failed — fall back to per-file search.
		if (diaFiles != null && diaFiles.length > 1 && initialSearchLibPaths != null) {
			List<String> damaged = addOutLibForInitialSearch
					? runPerFileFallbackAndRetryCombinedWithGenLib(tempParquetFile.getAbsolutePath(), initialSearchLibPaths)
					: runPerFileFallbackAndRetryCombinedNoQuant(tempParquetFile.getAbsolutePath(), initialSearchLibPaths);
			if (!damaged.isEmpty()) {
				logger.warn(getTaskName() + ": The following raw files produced no results (may be damaged): {}", damaged);
				System.err.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": The following raw files produced no results (may be damaged): " + String.join(", ", damaged));
			}
			java.util.Set<String> damagedSet = new java.util.HashSet<>(damaged);
			List<String> goodList = new ArrayList<>();
			for (String f : diaFiles) {
				if (!damagedSet.contains(f)) {
					goodList.add(f);
				}
			}
			if (!goodList.isEmpty()) {
				initialSearchSeparate = true;
				goodDiaFiles = goodList.toArray(new String[0]);
				logger.info(getTaskName() + ": initial search was separate; {} good files, {} damaged",
						goodDiaFiles.length, damaged.size());
				return true;
			}
		}

		logger.error(getTaskName() + ": First search failed - no seed.parquet and no per-file parquet outputs");
		return false;
	}

	@Override
	protected boolean iterativeSearch() {
		logger.info("Starting iterative search for DIA Seed (user library)...");

		MetaParameterDiaMarker markerParams = (metaPar instanceof MetaParameterDiaMarker) ? (MetaParameterDiaMarker) metaPar : null;
		if (markerParams == null || magDb == null) {
			logger.error(getTaskName() + ": MetaParameterDiaMarker and MAG database required for iterative seed workflow");
			return false;
		}

		if (!runFirstSearchUserLibrary()) {
			return false;
		}

		File iterativeFolder = getIterativeResultFolder();
		File tempParquetFile = new File(iterativeFolder, "seed.parquet");

		// Step 2: Read precursors and set protein info from MAG predicted TSVs
		logger.info("Reading precursors and setting protein information from MAG predicted TSVs...");
		DiaNNPrecursor[] precursors;
		if (tempParquetFile.exists()) {
			DiaNNParquetReader parquetReader = new DiaNNParquetReader(tempParquetFile.getAbsolutePath(),
					diannPar.getCut());
			precursors = parquetReader.getDiaNNPrecursors();
		} else {
			logger.warn(getTaskName()
					+ ": seed.parquet not found after fallback; aggregating precursors from per_file_*.parquet instead");
			precursors = readPrecursorsFromPerFileParquets(iterativeFolder, null);
		}
		if (precursors == null || precursors.length == 0) {
			logger.warn(getTaskName() + ": No precursors found in user library search");
			return false;
		}

		File predictedFolder = new File(magDb.getCurrentFile(), "predicted");
		if (!predictedFolder.isDirectory()) {
			logger.error(getTaskName() + ": MAG predicted folder not found: " + predictedFolder.getAbsolutePath());
			return false;
		}

		List<DiaNNPrecursor> precursorsWithPro = MarkerGenomeSelector.setProteinInfoFromPredictedTsVs(precursors,
				predictedFolder);
		if (precursorsWithPro == null || precursorsWithPro.isEmpty()) {
			logger.error(getTaskName() + ": No precursors could be mapped to MAG proteins; check predicted TSVs");
			return false;
		}
		logger.info("Mapped {} / {} precursors to MAG proteins", precursorsWithPro.size(), precursors.length);

		// Step 3: MarkerGenomeSelector (Focus or Comprehensive)
		String markerMode = markerParams.getMarkerMode();
		String mode;
		if ("COMPREHENSIVE".equalsIgnoreCase(markerMode) || "DEEP".equalsIgnoreCase(markerMode)) {
			mode = MarkerGenomeSelector.KNEE_MODE_COMPREHENSIVE;
		} else {
			mode = MarkerGenomeSelector.KNEE_MODE_FOCUS;
		}
		logger.info("Selecting genomes using MarkerGenomeSelector ({} mode)...", mode);

		File markerFasta = new File(magDb.getCurrentFile(), "Marker4.fasta");
		File genomesTaxaTsv = new File(magDb.getCurrentFile(), "genomes-taxa.tsv");
		if (!markerFasta.exists() || !genomesTaxaTsv.exists()) {
			logger.error(getTaskName() + ": Marker4.fasta or genomes-taxa.tsv not found for genome selection");
			return false;
		}

		// If marker peptides are < 5% of total identified peptides, use all peptides for genome selection (no privilege, no representative, greedy only)
		MarkerGenomeSelector.MarkerPeptideRatioResult ratioResult = MarkerGenomeSelector.computeMarkerPeptideRatio(precursorsWithPro, markerFasta);
		File effectiveMarkerFasta = markerFasta;
		if (ratioResult.ratio < MarkerGenomeSelector.MARKER_PEPTIDE_RATIO_THRESHOLD) {
			effectiveMarkerFasta = null;
			logger.info("Marker peptides ({}) are {}% of total identified peptides ({}); using all peptides for genome selection (greedy only, no PRIVILEGE/representative)",
					ratioResult.markerUniquePeptides, String.format("%.1f", ratioResult.ratio * 100), ratioResult.totalUniquePeptides);
		} else {
			logger.info("Marker peptide ratio: {} / {} = {}% (using marker-based genome selection)", ratioResult.markerUniquePeptides, ratioResult.totalUniquePeptides, String.format("%.1f", ratioResult.ratio * 100));
		}

		MarkerGenomeSelector selector = new MarkerGenomeSelector();
		selector.setKneeSelectionMode(mode);

		HashSet<String> genomeSet = new HashSet<>();
		File evidenceFile = new File(iterativeFolder, "genome_evidence.tsv");
		File rankedFile = new File(iterativeFolder, "genome_ranked.tsv");
		File selectedFile = new File(iterativeFolder, "genome_selected.tsv");
		selector.selectMarkerGenomesWithKnee(
				precursorsWithPro.toArray(new DiaNNPrecursor[0]),
				effectiveMarkerFasta, genomesTaxaTsv,
				evidenceFile.getAbsolutePath(), rankedFile.getAbsolutePath(), selectedFile.getAbsolutePath(),
				300, true);
		if (selectedFile.exists()) {
			try (BufferedReader reader = new BufferedReader(new FileReader(selectedFile))) {
				reader.readLine();
				String line;
				while ((line = reader.readLine()) != null) {
					line = line.trim();
					if (line.isEmpty() || line.startsWith("#")) break;
					String[] parts = line.split("\t");
					if (parts.length > 0 && !parts[0].trim().isEmpty()) {
						genomeSet.add(parts[0].trim());
					}
				}
			} catch (IOException e) {
				logger.error("Error reading selected genomes", e);
				return false;
			}
		}

		if (genomeSet.isEmpty()) {
			logger.warn(getTaskName() + ": No genomes selected");
			return false;
		}

		File combinedFastaFile = new File(iterativeFolder, "combined.fasta");
		boolean usePeptideBasedRound = markerParams.isUsePeptideBasedRound();
		if (usePeptideBasedRound) {
			double detectabilityThreshold = markerParams.getPeptideDetectabilityThreshold();
			if (!writeCombinedPeptideFastaFromSelectedGenomes(genomeSet, detectabilityThreshold, combinedFastaFile)) {
				logger.error(getTaskName() + ": failed to write peptide-based combined FASTA");
				return false;
			}
			logger.info(getTaskName() + ": using peptide-based iterative strategy for seed focus/comprehensive "
					+ "(threshold={})", detectabilityThreshold);
		} else {
			try (PrintWriter writer = new PrintWriter(combinedFastaFile)) {
				File dbFolder = magDb.getOriginalDbFolder();
				for (String genome : genomeSet) {
					File genomeDbFile = new File(dbFolder, genome + ".faa");
					if (!genomeDbFile.exists()) {
						genomeDbFile = new File(dbFolder, genome + ".fasta");
					}
					if (genomeDbFile.exists()) {
						try (BufferedReader reader = new BufferedReader(new FileReader(genomeDbFile))) {
							String line;
							while ((line = reader.readLine()) != null) {
								writer.println(line);
							}
						}
					} else {
						logger.warn("Genome file not found: {}", genome);
					}
				}
			} catch (IOException e) {
				logger.error("Error writing combined FASTA", e);
				return false;
			}
		}

		if (!combinedFastaFile.exists() || combinedFastaFile.length() == 0) {
			logger.error(getTaskName() + ": Failed to create combined FASTA");
			return false;
		}

		// Step 5: Generate combined library
		File combinedLibFile = new File(iterativeFolder, "combined.speclib");
		if (!combinedLibFile.exists()) {
			File altLibFile = new File(iterativeFolder, "combined.predicted.speclib");
			if (altLibFile.exists()) {
				altLibFile.renameTo(combinedLibFile);
			} else {
				diannPar.setThreads(threadCount);
				diaNNTask.addTask(diannPar, combinedFastaFile.getAbsolutePath(), combinedLibFile.getAbsolutePath(),
						false);
				diaNNTask.run(12);
				if (!combinedLibFile.exists() && altLibFile.exists()) {
					altLibFile.renameTo(combinedLibFile);
				}
			}
		}
		if (!combinedLibFile.exists()) {
			logger.error(getTaskName() + ": Failed to generate combined library");
			return false;
		}

		logger.info("Iterative search completed successfully");
		return true;
	}

	@Override
	protected boolean iden() {
		logger.info("Starting identification using combined library from iterative search...");
		updateMessage(getTaskName() + ": peptide and protein identification");

		combinedFile = new File(metaPar.getDbSearchResultFile(), "combined");
		if (!combinedFile.exists()) {
			combinedFile.mkdir();
		}
		// Required for exportReport() / getExportQuanPepFile()
		this.quan_pep_file = new File(combinedFile, "combined.pr_matrix.tsv");
		this.quan_pro_file = new File(combinedFile, "combined.pg_matrix.tsv");
		File combineParquetFile = new File(combinedFile, "combined.parquet");

		diaLibSearchPar.setMbr(true);
		diaLibSearchPar.setDiaFiles(diaFiles);
		diaLibSearchPar.setThreadCount(threadCount);
		diaLibSearchPar.setqValue(0.01);
		diannPar.setThreads(threadCount);
		diannPar.setMatrices(true);
		diannPar.setRelaxed_prot_inf(true);
		diannPar.setReanalyse(true);
		diannPar.setExportDecoy(true);

		if (addOutLibForInitialSearch) {
			// Mixed case: convert iterative microbiome combined.speclib to parquet, then
			// combine it with standalone parquet libraries for final search.
			File iterativeFolder = getIterativeResultFolder();
			File combinedLibFile = new File(iterativeFolder, "combined.speclib");
			if (!combinedLibFile.exists() || combinedLibFile.length() == 0) {
				logger.error(getTaskName() + ": Combined iterative microbiome library not found: "
						+ combinedLibFile.getAbsolutePath());
				return false;
			}
			File combinedLibParquetFile = new File(iterativeFolder, "combined.lib.parquet");
			if (!combinedLibParquetFile.exists() || combinedLibParquetFile.length() == 0) {
				logger.info(getTaskName()
						+ ": converting iterative combined.speclib to parquet library for mixed-library final search...");
				diannPar.setThreads(threadCount);
				diaNNTask.addSpeclib2ParquetTask(diannPar, combinedLibFile.getAbsolutePath(),
						combinedLibParquetFile.getAbsolutePath());
				diaNNTask.run(12);
			}
			if (!combinedLibParquetFile.exists() || combinedLibParquetFile.length() == 0) {
				logger.error(getTaskName() + ": failed to generate parquet library from combined.speclib: {}",
						combinedLibParquetFile.getAbsolutePath());
				return false;
			}
			if (!quan_pep_file.exists() || !combineParquetFile.exists()) {
				diannPar.setMatrices(true);
				List<String> combinedLibPaths = new ArrayList<>();
				combinedLibPaths.add(combinedLibParquetFile.getAbsolutePath());
				if (standaloneSearchLibPaths != null) {
					for (String p : standaloneSearchLibPaths) {
						if (p == null || !p.toLowerCase().endsWith(".parquet")) {
							logger.error(getTaskName()
									+ ": mixed-library mode accepts parquet libraries only; invalid standalone library: {}",
									p);
							return false;
						}
						combinedLibPaths.add(p);
					}
				}
				String[] combinedLibArr = combinedLibPaths.toArray(new String[0]);
				logger.info(
						"Searching with iterative microbiome parquet library + {} standalone parquet libraries -> combined.parquet",
						standaloneSearchLibPaths != null ? standaloneSearchLibPaths.length : 0);
				diaNNTask.addTask(diannPar, diaLibSearchPar, combineParquetFile.getAbsolutePath(), combinedLibArr);
				diaNNTask.run(diaFiles.length * 12);
			}
		} else {
			// All microbiome: use combined library from iterative search directly
			File iterativeFolder = getIterativeResultFolder();
			File combinedLibFile = new File(iterativeFolder, "combined.speclib");
			if (!combinedLibFile.exists()) {
				logger.error(getTaskName() + ": Combined library not found: " + combinedLibFile.getAbsolutePath());
				return false;
			}
			if (!quan_pep_file.exists() || !combineParquetFile.exists()) {
				diaNNTask.addTask(diannPar, diaLibSearchPar, combineParquetFile.getAbsolutePath(),
						new String[] { combinedLibFile.getAbsolutePath() });
				diaNNTask.run(diaFiles.length * 12);
			}
		}

		if (!quan_pep_file.exists() || !combineParquetFile.exists()) {
			logger.error(getTaskName() + ": Identification failed - result files not created");
			return false;
		}
		logger.info("Identification completed successfully");
		return true;
	}

	private boolean writeCombinedPeptideFastaFromSelectedGenomes(HashSet<String> genomeSet, double detectabilityThreshold,
			File outputFasta) {
		if (genomeSet == null || genomeSet.isEmpty()) {
			return false;
		}
		File predictedFolder = new File(magDb.getCurrentFile(), "predicted");
		if (!predictedFolder.isDirectory()) {
			logger.error(getTaskName() + ": predicted folder not found for peptide-based iterative FASTA: {}",
					predictedFolder.getAbsolutePath());
			return false;
		}
		HashMap<String, HashSet<String>> pepProMap = new HashMap<>();
		for (String genome : genomeSet) {
			File predictedFile = new File(predictedFolder, genome + ".predicted.tsv");
			if (!predictedFile.isFile()) {
				logger.warn(getTaskName() + ": predicted file missing for selected genome: {}", genome);
				continue;
			}
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
					pepProMap.computeIfAbsent(peptide, k -> new HashSet<>()).add(protein);
				}
			} catch (IOException e) {
				logger.warn(getTaskName() + ": error reading predicted file {}", predictedFile.getAbsolutePath(), e);
			}
		}
		if (pepProMap.isEmpty()) {
			logger.error(getTaskName() + ": no peptides collected from selected genomes under threshold={}",
					detectabilityThreshold);
			return false;
		}
		boolean ok = DeepIterativeCommonSupport.writeCombinedFastaFromPeptideProteinMap(outputFasta, pepProMap, logger,
				getTaskName());
		if (ok) {
			logger.info(getTaskName() + ": wrote peptide-based combined FASTA with {} stripped peptides", pepProMap.size());
		}
		return ok;
	}

	/**
	 * Aggregate DiaNN precursors from per_file_*.parquet files in the given directory.
	 *
	 * @param dir    directory containing per_file_*.parquet
	 * @param qValue optional q-value cutoff; if null, use default reader behaviour
	 * @return concatenated precursors or null if none
	 */
	private DiaNNPrecursor[] readPrecursorsFromPerFileParquets(File dir, Double qValue) {
		if (dir == null || !dir.isDirectory()) {
			return null;
		}
		String cut = diannPar != null && diannPar.getCut() != null ? diannPar.getCut() : "KR";
		List<DiaNNPrecursor> all = new ArrayList<>();
		for (int i = 0; i < (diaFiles != null ? diaFiles.length : Integer.MAX_VALUE); i++) {
			File pf = new File(dir, "per_file_" + i + ".parquet");
			if (!pf.exists()) {
				break;
			}
			if (pf.length() <= 0) {
				continue;
			}
			try {
				DiaNNParquetReader reader = new DiaNNParquetReader(pf.getAbsolutePath(), cut);
				DiaNNPrecursor[] prec = (qValue == null)
						? reader.getDiaNNPrecursors()
						: reader.getDiaNNPrecursors(qValue);
				if (prec != null) {
					for (DiaNNPrecursor p : prec) {
						if (p != null) {
							all.add(p);
						}
					}
				}
			} catch (Exception e) {
				logger.warn(getTaskName() + ": failed to read per-file parquet " + pf.getAbsolutePath(), e);
			}
		}
		return all.isEmpty() ? null : all.toArray(new DiaNNPrecursor[0]);
	}
}
