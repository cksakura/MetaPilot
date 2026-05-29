package uk.ac.quadram.metapilot.task.dia;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.ConcurrentHashMap;

import org.apache.logging.log4j.Logger;

import javafx.application.Platform;
import javafx.scene.control.ProgressBar;

import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.task.dia.par.MetaParameterDiaTg;
import uk.ac.quadram.metapilot.task.dia.gui.MetaPilotDiaTargetPanel.TargetLibraryInfo;
import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNPrecursor;
import uk.ac.quadram.metapilot.core.dia.MagPredictedPeptideProteinMapper;
import uk.ac.quadram.metapilot.core.prodb.FastaManager;

/**
 * MetaDiaTargetTask - A specialized DIA task for target library analysis
 * 
 * This class extends MetaDiaTask to provide specific functionality for DIA
 * Target workflow, which uses pre-generated target libraries instead of
 * creating sample-specific libraries.
 */
public class MetaDiaTargetTask extends MetaDiaTask {

	// Logger is inherited from MetaDiaTask -> MetaPilotMagTask -> MetaPilotTask as an instance field
	private static final String taskName = "DIA Target analysis";

	// DIA Target specific fields
	private String targetLibraryPath;
	private String[] selectedTaxa;
	private String[] selectedFunctions;
	private String taxonomicLevel;
	private String functionalLevel;
	private TargetLibraryInfo targetLibraryInfo;

	/**
	 * Constructor for MetaDiaTargetTask
	 * 
	 * @param targetParams The DIA Target specific parameters
	 * @param mrp          Resource parameters
	 * @param bar1         Main progress bar
	 * @param bar2         Secondary progress bar
	 */
	public MetaDiaTargetTask(MetaParameterDiaTg targetParams, MetaResourceParameters mrp, ProgressBar bar1,
			ProgressBar bar2) {
		// Call parent constructor with the base DIA parameters
		super(targetParams, mrp, bar1, bar2);

		this.targetLibraryPath = targetParams.getLibraryPath();
		this.selectedTaxa = targetParams.getSelectedTaxa();
		this.selectedFunctions = targetParams.getSelectedFunctions();
		this.taxonomicLevel = targetParams.getTaxonomicLevel();
		this.functionalLevel = targetParams.getFunctionalLevel();

		logger.info("Initializing MetaDiaTargetTask...");
		System.out.println(new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss").format(new Date()) + "\t" + getTaskName()
				+ ": initializing task");

		// Initialize result folder for DIA Target: use library path parent folder by default
		initResultFolderIfNeeded();

		logger.info("DIA Target workflow initialized with library search mode");
		logger.info("Target library path: {}", targetLibraryPath);
		logger.info("Selected taxa count: {}", selectedTaxa != null ? selectedTaxa.length : 0);
		logger.info("Selected functions count: {}", selectedFunctions != null ? selectedFunctions.length : 0);
	}

	@Override
	protected String getTaskName() {
		return taskName;
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
		if (!(metaPar instanceof MetaParameterDiaTg)) {
			super.initResultFolderIfNeeded();
			return;
		}
		MetaParameterDiaTg targetParams = (MetaParameterDiaTg) metaPar;
		File libFile = new File(targetParams.getLibraryPath());
		File dir = libFile.getParentFile();
		if (dir != null && !dir.exists()) {
			dir.mkdirs();
		}
		metaPar.setDbSearchResultFile(dir);
		this.resultFolderFile = dir;
	}

	/**
	 * Set target library information for peptide mapping
	 */
	public void setTargetLibraryInfo(TargetLibraryInfo libraryInfo) {
		this.targetLibraryInfo = libraryInfo;
		logger.info("Target library info set: {} genomes, {} functions, {} proteins", libraryInfo.getGenomes().size(),
				libraryInfo.getFunctions().size(), libraryInfo.getProteins().size());
	}

	@Override
	protected Boolean call() throws Exception {
		updateMessage(getTaskName());
		updateProgress(0, 100);

		logger.info(getTaskName() + ": started");
		System.out.println(
				new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss").format(new Date()) + "\t" + getTaskName() + ": started");

		// Log detailed parameter information
		logger.info("DIA Target workflow parameters:");
		logger.info("  Result directory: " + metaPar.getResult());
		logger.info("  Thread count: " + metaPar.getThreadCount());
		logger.info("  Target library path: " + targetLibraryPath);
		logger.info("  Selected taxa count: " + (selectedTaxa != null ? selectedTaxa.length : 0));
		logger.info("  Selected functions count: " + (selectedFunctions != null ? selectedFunctions.length : 0));

		if (metadata != null) {
			String[] rawFiles = metadata.getRawFiles();
			logger.info("  Raw files count: " + (rawFiles != null ? rawFiles.length : 0));
			if (rawFiles != null) {
				for (int i = 0; i < rawFiles.length; i++) {
					logger.info("    Raw file " + (i + 1) + ": " + rawFiles[i]);
				}
			}
		}

		// Step 0: Copy report templates in background
		updateProgress(5, 100);
		updateMessage(getTaskName() + ": Initializing analysis...");
		if (Platform.isFxApplicationThread()) {
			secondBar.setProgress(ProgressBar.INDETERMINATE_PROGRESS);
		} else {
			Platform.runLater(() -> secondBar.setProgress(ProgressBar.INDETERMINATE_PROGRESS));
		}

		if (!copyReportTemplatesInBackground()) {
			logger.error(getTaskName() + ": Failed to copy report templates");
			return false;
		}

		// Step 1: Spectra Conversion (5% - 25%)
		logger.info("Step 1: Starting spectra format conversion...");
		updateMessage(getTaskName() + ": Converting spectra format...");
		updateProgress(10, 100);
		if (Platform.isFxApplicationThread()) {
			secondBar.setProgress(ProgressBar.INDETERMINATE_PROGRESS);
		} else {
			Platform.runLater(() -> secondBar.setProgress(ProgressBar.INDETERMINATE_PROGRESS));
		}

		boolean convert = spectraConvert();
		if (!convert) {
			logger.error(getTaskName() + ": error in spectra format conversion");
			System.err.println(new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss").format(new Date()) + "\t" + getTaskName()
					+ ": error in spectra format conversion");
			updateMessage(getTaskName() + ": Failed - Spectra conversion error");
			if (Platform.isFxApplicationThread()) {
				secondBar.setProgress(0.0);
			} else {
				Platform.runLater(() -> secondBar.setProgress(0.0));
			}
			return false;
		}
		updateProgress(25, 100);
		logger.info("Step 1: Spectra format conversion completed successfully");

		// Step 2: Target library setup (25% - 35%)
		logger.info("Step 2: Using pre-generated target library...");
		updateMessage(getTaskName() + ": Using target library...");
		updateProgress(30, 100);
		if (Platform.isFxApplicationThread()) {
			secondBar.setProgress(ProgressBar.INDETERMINATE_PROGRESS);
		} else {
			Platform.runLater(() -> secondBar.setProgress(ProgressBar.INDETERMINATE_PROGRESS));
		}

		boolean useTargetLibrary = createSSLibrary();
		if (!useTargetLibrary) {
			logger.error(getTaskName() + ": error in target library setup");
			System.err.println(new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss").format(new Date()) + "\t" + getTaskName()
					+ ": error in target library setup");
			updateMessage(getTaskName() + ": Failed - Target library setup error");
			if (Platform.isFxApplicationThread()) {
				secondBar.setProgress(0.0);
			} else {
				Platform.runLater(() -> secondBar.setProgress(0.0));
			}
			return false;
		}
		updateProgress(35, 100);
		logger.info("Step 2: Target library setup completed successfully");

		// Step 3: Identification using target library (35% - 70%)
		logger.info("Step 3: Starting identification using target library...");
		updateMessage(getTaskName() + ": Running identification with target library...");
		updateProgress(40, 100);
		if (Platform.isFxApplicationThread()) {
			secondBar.setProgress(ProgressBar.INDETERMINATE_PROGRESS);
		} else {
			Platform.runLater(() -> secondBar.setProgress(ProgressBar.INDETERMINATE_PROGRESS));
		}

		boolean iden = iden();
		if (!iden) {
			logger.error(getTaskName() + ": error in identification");
			System.err.println(new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss").format(new Date()) + "\t" + getTaskName()
					+ ": error in identification");
			updateMessage(getTaskName() + ": Failed - Identification error");
			if (Platform.isFxApplicationThread()) {
				secondBar.setProgress(0.0);
			} else {
				Platform.runLater(() -> secondBar.setProgress(0.0));
			}
			return false;
		}
		updateProgress(70, 100);
		logger.info("Step 3: Identification completed successfully");

		// Step 4: Export Report (70% - 100%)
		logger.info("Step 4: Starting report generation...");
		updateMessage(getTaskName() + ": Generating reports...");
		updateProgress(75, 100);
		if (Platform.isFxApplicationThread()) {
			secondBar.setProgress(ProgressBar.INDETERMINATE_PROGRESS);
		} else {
			Platform.runLater(() -> secondBar.setProgress(ProgressBar.INDETERMINATE_PROGRESS));
		}

		boolean exportReport = exportReport();
		if (!exportReport) {
			logger.error(getTaskName() + ": error in report generation");
			System.err.println(new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss").format(new Date()) + "\t" + getTaskName()
					+ ": error in report generation");
			updateMessage(getTaskName() + ": Failed - Report generation error");
			if (Platform.isFxApplicationThread()) {
				secondBar.setProgress(0.0);
			} else {
				Platform.runLater(() -> secondBar.setProgress(0.0));
			}
			return false;
		}
		updateProgress(95, 100);
		logger.info("Step 4: Report generation completed successfully");

		// Update progress for completion
		updateProgress(100, 100);
		updateMessage(getTaskName() + ": Analysis completed successfully!");
		if (Platform.isFxApplicationThread()) {
			secondBar.setProgress(1.0);
		} else {
			Platform.runLater(() -> secondBar.setProgress(1.0));
		}

		logger.info(getTaskName() + ": finished - analysis completed successfully");
		System.out.println(new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss").format(new Date()) + "\t" + getTaskName()
				+ ": finished - analysis completed successfully");
		return true;
	}

	@Override
	protected boolean createSSLibrary() {
		logger.info("Starting target library setup for DIA Target workflow...");

		// For DIA Target, we use pre-generated target libraries instead of creating
		// sample-specific ones
		logger.info("DIA Target workflow: Using pre-generated target library");
		System.out.println(new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss").format(new Date()) + "\t" + getTaskName()
				+ ": using target library");

		// Validate target library path
		if (targetLibraryPath == null || targetLibraryPath.trim().isEmpty()) {
			logger.error(getTaskName() + ": Target library path is null or empty");
			System.err.println(new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss").format(new Date()) + "\t" + getTaskName()
					+ ": Target library path is null or empty");
			return false;
		}

		File targetLibraryFile = new File(targetLibraryPath);
		if (!targetLibraryFile.exists()) {
			logger.error(getTaskName() + ": Target library file does not exist: " + targetLibraryPath);
			System.err.println(new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss").format(new Date()) + "\t" + getTaskName()
					+ ": Target library file does not exist: " + targetLibraryPath);
			return false;
		}

		logger.info("Target library validated successfully: {}", targetLibraryPath);
		logger.info("Target library file size: {} bytes", targetLibraryFile.length());

		// Note: For DIA Target, we always use library search mode
		// The parent constructor handles setting the appropriate search mode

		logger.info("DIA Target library setup completed successfully");
		System.out.println(new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss").format(new Date()) + "\t" + getTaskName()
				+ ": target library setup finished");

		return true;
	}

	@Override
	protected boolean iden() {
		logger.info("Starting identification process using target library...");
		updateMessage(getTaskName() + ": peptide and protein identification with target library");
		System.out.println(new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss").format(new Date()) + "\t" + getTaskName()
				+ ": peptide and protein identification started");

		logger.info("Configuring DiaLibSearchPar for DIA Target:");
		logger.info("  DIA files count: {}", diaFiles != null ? diaFiles.length : 0);
		logger.info("  Thread count: {}", threadCount);
		logger.info("  Q-value threshold: 0.01");
		logger.info("  MBR enabled: true");
		logger.info("  Target library path: {}", targetLibraryPath);

		diaLibSearchPar.setMbr(true);
		diaLibSearchPar.setDiaFiles(diaFiles);
		diaLibSearchPar.setThreadCount(threadCount);
		diaLibSearchPar.setqValue(0.01);

		diannPar.setThreads(threadCount);
		diannPar.setMatrices(true);
		diannPar.setRelaxed_prot_inf(true);
		diannPar.setReanalyse(true);

		combinedFile = new File(metaPar.getDbSearchResultFile(), "combined");
		if (!combinedFile.exists()) {
			combinedFile.mkdir();
		}
		this.quan_pep_file = new File(combinedFile, "combined.pr_matrix.tsv");
		this.quan_pro_file = new File(combinedFile, "combined.pg_matrix.tsv");

		File combineParquetFile = new File(combinedFile, "combined.parquet");

		if (!quan_pep_file.exists() || !combineParquetFile.exists()) {
			logger.info(getTaskName() + ": target library search started");
			System.out.println(new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss").format(new Date()) + "\t" + getTaskName()
					+ ": target library search started");

			diannPar.setExportDecoy(true);

			// Use the target library path for DIA Target workflow
			logger.info("Using target library for identification: {}", targetLibraryPath);
			diaNNTask.addTask(diannPar, diaLibSearchPar, combineParquetFile.getAbsolutePath(),
					new String[] { targetLibraryPath });

			logger.info("Starting DiaNN identification with target library using {} threads", threadCount);
			diaNNTask.run(diaFiles.length * 12);

			logger.info("DiaNN identification with target library completed");
		} else {
			logger.info(getTaskName() + ": identification results already exist, skipping identification");
			System.out.println(new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss").format(new Date()) + "\t" + getTaskName()
					+ ": identification results already exist");
		}

		// Verify that the identification files were created
		if (!quan_pep_file.exists() || !combineParquetFile.exists()) {
			logger.error(getTaskName() + ": identification failed - result files not created");
			System.err.println(new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss").format(new Date()) + "\t" + getTaskName()
					+ ": identification failed - result files not created");
			return false;
		}

		logger.info("Target library identification completed successfully");
		logger.info("Peptide results file: {}", quan_pep_file.getAbsolutePath());
		logger.info("Protein results file: {}", quan_pro_file.getAbsolutePath());
		logger.info("Parquet results file: {}", combineParquetFile.getAbsolutePath());

		System.out.println(new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss").format(new Date()) + "\t" + getTaskName()
				+ ": identification finished");

		return true;
	}

	@Override
	protected boolean iterativeSearch() {
		// DIA Target workflow does not use iterative search - we use pre-generated
		// target libraries
		logger.info("DIA Target workflow: Skipping iterative search (using target library)");
		System.out.println(new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss").format(new Date()) + "\t" + getTaskName()
				+ ": skipping iterative search");
		return true;
	}

	@Override
	protected boolean exportReport() {
		return super.exportReport();
	}

	@Override
	protected DiaNNPrecursor[] mapPrecursorsToGenomes(DiaNNPrecursor[] precursors) {
		if (precursors == null || precursors.length == 0) {
			return precursors;
		}

		HashMap<String, String> targetProteinPeptideMap = readTargetLibraryProteinPeptides();
		HashSet<String> peptideSet = new HashSet<String>();
		HashMap<String, String> decoyPeptideTargetMap = new HashMap<String, String>();
		int decoyCount = 0;
		for (DiaNNPrecursor pp : precursors) {
			if (pp == null) {
				continue;
			}
			String sequence = pp.getSeqString();
			if (sequence == null || sequence.isEmpty()) {
				continue;
			}
			peptideSet.add(sequence);
			if (!pp.getIsDecoy()) {
				continue;
			}
			decoyCount++;
			String[] decoyPros = pp.getProGroups();
			if (decoyPros == null) {
				continue;
			}
			for (String decoyPro : decoyPros) {
				String targetSeq = targetProteinPeptideMap.get(decoyPro);
				if (targetSeq == null && decoyPro != null && decoyPro.startsWith(FastaManager.shuffle)) {
					targetSeq = targetProteinPeptideMap.get(decoyPro.substring(FastaManager.shuffle.length()));
				}
				if (targetSeq != null && targetSeq.length() == sequence.length()) {
					peptideSet.add(targetSeq);
					decoyPeptideTargetMap.put(sequence, targetSeq);
				}
			}
		}

		logger.info(getTaskName() + ": target mapping sequence count={}, decoy precursor count={}", peptideSet.size(),
				decoyCount);
		ConcurrentHashMap<String, HashMap<String, HashSet<String>>> matchedMap = getGenomeFromPep(peptideSet,
				threadCount);
		if (matchedMap == null || matchedMap.isEmpty()) {
			logger.warn(getTaskName() + ": no target peptides mapped to selected MAG proteins");
			return new DiaNNPrecursor[0];
		}

		HashMap<String, HashSet<String>> peptideProteinMap = new HashMap<String, HashSet<String>>();
		for (HashMap<String, HashSet<String>> proteinPeptideMap : matchedMap.values()) {
			for (String protein : proteinPeptideMap.keySet()) {
				for (String peptide : proteinPeptideMap.get(protein)) {
					peptideProteinMap.computeIfAbsent(peptide, k -> new HashSet<String>()).add(protein);
				}
			}
		}

		ArrayList<DiaNNPrecursor> mappedPrecursors = new ArrayList<DiaNNPrecursor>();
		HashSet<String> missingPeptides = new HashSet<String>();
		int mappedDecoys = 0;
		for (DiaNNPrecursor pp : precursors) {
			if (pp == null) {
				continue;
			}
			String sequence = pp.getSeqString();
			HashSet<String> proteins = peptideProteinMap.get(sequence);
			if (proteins != null && !proteins.isEmpty()) {
				pp.setProGroups(proteins.toArray(new String[0]));
				mappedPrecursors.add(pp);
				continue;
			}

			String targetSequence = decoyPeptideTargetMap.get(sequence);
			HashSet<String> targetProteins = targetSequence == null ? null : peptideProteinMap.get(targetSequence);
			if (targetProteins != null && !targetProteins.isEmpty()) {
				ArrayList<String> decoyProteins = new ArrayList<String>();
				for (String protein : targetProteins) {
					decoyProteins.add(FastaManager.shuffle + protein);
				}
				pp.setProGroups(decoyProteins.toArray(new String[0]));
				mappedPrecursors.add(pp);
				mappedDecoys++;
			} else if (sequence != null) {
				missingPeptides.add(sequence);
			}
		}

		logger.info(getTaskName()
				+ ": target-restricted precursor mapping finished, {} / {} precursors kept, {} decoys mapped, {} peptides missing",
				mappedPrecursors.size(), precursors.length, mappedDecoys, missingPeptides.size());
		System.out.println(new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss").format(new Date()) + "\t" + getTaskName()
				+ ": target-restricted precursor mapping finished, " + mappedPrecursors.size() + " / "
				+ precursors.length + " precursors kept");
		return mappedPrecursors.toArray(new DiaNNPrecursor[0]);
	}

	private HashMap<String, String> readTargetLibraryProteinPeptides() {
		HashMap<String, String> proteinPeptideMap = new HashMap<String, String>();
		if (targetLibraryPath == null || targetLibraryPath.isBlank()) {
			return proteinPeptideMap;
		}
		File targetLibraryFile = new File(targetLibraryPath);
		File libraryDir = targetLibraryFile.getParentFile();
		if (libraryDir == null) {
			return proteinPeptideMap;
		}
		String libraryName = libraryDir.getName();
		File targetFastaFile = new File(libraryDir, libraryName + ".fasta");
		if (!targetFastaFile.exists()) {
			logger.warn(getTaskName() + ": target library FASTA not found for decoy mapping: {}",
					targetFastaFile.getAbsolutePath());
			return proteinPeptideMap;
		}

		try (BufferedReader reader = new BufferedReader(new FileReader(targetFastaFile))) {
			String line;
			String protein = null;
			while ((line = reader.readLine()) != null) {
				if (line.startsWith(">")) {
					String[] parts = line.split("\\s+");
					protein = parts.length > 1 ? parts[1] : parts[0].substring(1);
				} else if (protein != null) {
					proteinPeptideMap.put(protein, line.trim());
					protein = null;
				}
			}
			logger.info(getTaskName() + ": read {} proteins from target library FASTA for decoy mapping",
					proteinPeptideMap.size());
		} catch (IOException e) {
			logger.warn(getTaskName() + ": error reading target library FASTA for decoy mapping: {}",
					targetFastaFile.getAbsolutePath(), e);
		}
		return proteinPeptideMap;
	}

	/**
	 * DIA Target specific implementation of getGenomeFromPep Maps identified
	 * peptides to genomes and proteins, but only for genomes that were included in
	 * the target library generation
	 */
	@Override
	protected ConcurrentHashMap<String, HashMap<String, HashSet<String>>> getGenomeFromPep(HashSet<String> allPepSet,
			int threadCount) {
		logger.info("Starting peptide-to-genome mapping for DIA Target workflow");

		if (targetLibraryInfo == null) {
			logger.error("Target library information not set - cannot perform peptide mapping");
			return new ConcurrentHashMap<>();
		}

		logger.info("Target library filtering info:");
		logger.info("  - Genome filtering: {} ({} genomes)",
				targetLibraryInfo.hasGenomeFiltering() ? "ENABLED" : "DISABLED", targetLibraryInfo.getGenomes().size());
		logger.info("  - Function filtering: {} ({} functions)",
				targetLibraryInfo.hasFunctionFiltering() ? "ENABLED" : "DISABLED",
				targetLibraryInfo.getFunctions().size());
		logger.info("  - Protein filtering: {} ({} proteins)",
				targetLibraryInfo.hasProteinFiltering() ? "ENABLED" : "DISABLED",
				targetLibraryInfo.getProteins().size());

		File currentFile = super.magDb.getCurrentFile();
		File predictPepFile = new File(currentFile, "predicted");
		HashSet<String> allowedGenomes = targetLibraryInfo.hasGenomeFiltering()
				? new HashSet<>(targetLibraryInfo.getGenomes())
				: null;
		HashSet<String> allowedProteins = targetLibraryInfo.hasProteinFiltering()
				? new HashSet<>(targetLibraryInfo.getProteins())
				: null;
		MagPredictedPeptideProteinMapper.Result result = MagPredictedPeptideProteinMapper.mapPeptides(allPepSet, predictPepFile,
				allowedGenomes, allowedProteins, threadCount);
		if (result == null) {
			return null;
		}

		allPepSet.clear();
		allPepSet.addAll(result.matchedPeptides);
		logger.info("DIA Target peptide mapping completed:");
		logger.info("  - Matched peptides: {}", allPepSet.size());
		logger.info("  - Mapped genomes: {}", result.genomeProteinPeptides.size());
		logger.info("  - Total proteins mapped: {}",
				result.genomeProteinPeptides.values().stream().mapToInt(HashMap::size).sum());
		System.out.println(new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss").format(new Date()) + "\t" + getTaskName()
				+ ": DIA Target peptide mapping finished, " + allPepSet.size() + " peptides from "
				+ result.genomeProteinPeptides.size() + " genomes were matched");

		return result.genomeProteinPeptides;
	}

	// Getters for DIA Target specific fields
	public String getTargetLibraryPath() {
		return targetLibraryPath;
	}

	public String[] getSelectedTaxa() {
		return selectedTaxa;
	}

	public String[] getSelectedFunctions() {
		return selectedFunctions;
	}

	public String getTaxonomicLevel() {
		return taxonomicLevel;
	}

	public String getFunctionalLevel() {
		return functionalLevel;
	}

}
