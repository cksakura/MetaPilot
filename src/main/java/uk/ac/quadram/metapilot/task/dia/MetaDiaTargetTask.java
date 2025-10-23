package uk.ac.quadram.metapilot.task.dia;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.math.BigInteger;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.sql.SQLException;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import org.apache.commons.lang3.math.NumberUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.dom4j.DocumentException;

import javafx.application.Platform;
import javafx.scene.control.ProgressBar;

import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters.ToolType;
import uk.ac.quadram.metapilot.task.dia.par.MetaParameterDiaTg;
import uk.ac.quadram.metapilot.task.io.MetaTreeHandler;
import uk.ac.quadram.metapilot.task.io.ReportManifestGenerator;
import uk.ac.quadram.metapilot.task.io.pep.MetaPeptide;
import uk.ac.quadram.metapilot.task.io.pro.MetaProtein;
import uk.ac.quadram.metapilot.task.io.pro.MetaProteinXMLReader2;
import uk.ac.quadram.metapilot.task.mag.HMGenomeHtmlWriter;
import uk.ac.quadram.metapilot.task.mag.MagFuncSearcher;
import uk.ac.quadram.metapilot.task.mag.MetaProteinAnnoMag;
import uk.ac.quadram.metapilot.task.dia.gui.MetaPilotDiaTargetPanel.TargetLibraryInfo;
import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNParquetReader;
import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNPrecursor;
import uk.ac.quadram.metapilot.dbSearch.diann.DiannPrMatrixReader;
import uk.ac.quadram.metapilot.quant.flashLFQ.FlashLfqQuanPeptide;
import uk.ac.quadram.metapilot.core.prodb.FastaManager;
import uk.ac.quadram.metapilot.core.taxonomy.TaxonomyRanks;
import uk.ac.quadram.metapilot.core.tools.FormatTool;
import uk.ac.quadram.metapilot.task.MetaPythonReportTask;
import uk.ac.quadram.metapilot.task.MetaPythonReportCopyTask;

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
			secondBar.setProgress(0.0);
		} else {
			Platform.runLater(() -> secondBar.setProgress(0.0));
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
		logger.info(getTaskName() + ": protein inference started for DIA Target workflow");
		System.out.println(new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss").format(new Date()) + "\t" + getTaskName()
				+ ": protein inference started");

		// Get the target library FASTA file
		File targetLibraryFile = new File(targetLibraryPath);
		File libraryDir = targetLibraryFile.getParentFile();
		String libraryName = libraryDir.getName();
		File targetFastaFile = new File(libraryDir, libraryName + ".fasta");

		logger.info("Using target library FASTA file: {}", targetFastaFile.getAbsolutePath());

		if (!targetFastaFile.exists()) {
			logger.error(getTaskName() + ": Target library FASTA file not found: " + targetFastaFile.getAbsolutePath());
			System.err.println(new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss").format(new Date()) + "\t" + getTaskName()
					+ ": Target library FASTA file not found: " + targetFastaFile.getAbsolutePath());
			return false;
		}

		HashMap<String, String> fastaProPepMap = new HashMap<String, String>();

		// Read the target library FASTA file
		try (BufferedReader reader = new BufferedReader(new FileReader(targetFastaFile))) {
			String line = null;
			while ((line = reader.readLine()) != null) {
				if (line.startsWith(">")) {
					String proName = line.split("\\s")[1];
					fastaProPepMap.put(proName, reader.readLine());
				}
			}
			logger.info("Read {} proteins from target library FASTA file", fastaProPepMap.size());
		} catch (IOException e) {
			logger.error(getTaskName() + ": error in reading the target library from " + targetFastaFile, e);
			System.err.println(new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss").format(new Date()) + "\t" + getTaskName()
					+ ": error in reading the target library from " + targetFastaFile);
			return false;
		}

		// Continue with the rest of the export report logic using the target library
		File parquetFile = new File(combinedFile, "combined.parquet");
		DiaNNParquetReader reader = new DiaNNParquetReader(parquetFile.getAbsolutePath(), this.diannPar.getCut());
		this.diaPp = reader.getDiaNNPrecursors(0.01);

		logger.info(getTaskName() + ": " + this.diaPp.length + " PSMs were obtained with Global.Q.Value<0.01");
		System.out.println(new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss").format(new Date()) + "\t" + getTaskName() + ": "
				+ this.diaPp.length + " PSMs were obtained with Global.Q.Value<0.01");

		int originalDecoyCount = 0;
		HashSet<String> pepSet = new HashSet<String>();
		HashMap<String, String> decoyPepMap = new HashMap<String, String>();
		if (diaPp != null && diaPp.length > 0) {
			for (int i = 0; i < diaPp.length; i++) {
				String sequence = diaPp[i].getSeqString();
				pepSet.add(sequence);
				if (diaPp[i].getIsDecoy()) {
					String[] decoyPros = diaPp[i].getProGroups();
					for (int j = 0; j < decoyPros.length; j++) {
						String targetSeq = fastaProPepMap.get(decoyPros[j]);
						if (targetSeq != null && targetSeq.length() == sequence.length()) {
							pepSet.add(targetSeq);
							decoyPepMap.put(sequence, targetSeq);
						}
					}
					originalDecoyCount++;
				}
			}
		}

		logger.info(getTaskName() + ": total sequence count: " + pepSet.size() + "; decoy precursor count: "
				+ originalDecoyCount);
		System.out
				.println(format.format(new Date()) + "\t" + getTaskName() + ": total sequence count: " + pepSet.size());

		ConcurrentHashMap<String, HashMap<String, HashSet<String>>> matchedMap = getGenomeFromPep(pepSet, threadCount);
		if (matchedMap == null) {
			return false;
		}

		HashMap<String, HashSet<String>> pepProMap = new HashMap<String, HashSet<String>>();
		for (String genome : matchedMap.keySet()) {
			HashMap<String, HashSet<String>> proPepMap = matchedMap.get(genome);
			for (String pro : proPepMap.keySet()) {
				HashSet<String> pset = proPepMap.get(pro);
				for (String pep : pset) {
					if (pepProMap.containsKey(pep)) {
						pepProMap.get(pep).add(pro);
					} else {
						HashSet<String> proSet = new HashSet<String>();
						proSet.add(pro);
						pepProMap.put(pep, proSet);
					}
				}
			}
		}

		logger.info(getTaskName() + ": matched sequence count: " + pepProMap.size());
		System.out.println(new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss").format(new Date()) + "\t" + getTaskName()
				+ ": matched sequence count: " + pepProMap.size());

		int decoyPPCount = 0;
		HashMap<String, HashSet<String>> genomeMap = new HashMap<String, HashSet<String>>();
		HashMap<String, HashSet<String>> proteinPepMap = new HashMap<String, HashSet<String>>();
		HashSet<String> tempPepSet = new HashSet<String>();
		HashSet<String> missSet = new HashSet<String>();
		for (DiaNNPrecursor pp : diaPp) {
			String sequence = pp.getSeqString();
			if (pepProMap.containsKey(sequence)) {
				for (String pro : pepProMap.get(sequence)) {
					String genome = pro.split("_")[0];
					if (genomeMap.containsKey(genome)) {
						genomeMap.get(genome).add(sequence);
					} else {
						HashSet<String> set = new HashSet<String>();
						set.add(sequence);
						genomeMap.put(genome, set);
					}
					if (proteinPepMap.containsKey(pro)) {
						proteinPepMap.get(pro).add(sequence);
					} else {
						HashSet<String> set = new HashSet<String>();
						set.add(sequence);
						proteinPepMap.put(pro, set);
					}
				}
				tempPepSet.add(sequence);
			} else {
				// DIA Target workflow: No host file processing, so skip hostPepProMap check
				if (decoyPepMap.containsKey(sequence)) {
					String targetSequence = decoyPepMap.get(sequence);
					if (pepProMap.containsKey(targetSequence)) {
						for (String pro : pepProMap.get(targetSequence)) {
							String genome = pro.split("_")[0];
							String decoyGenome = FastaManager.shuffle + genome;
							String decoyPro = FastaManager.shuffle + pro;

							if (genomeMap.containsKey(decoyGenome)) {
								genomeMap.get(decoyGenome).add(sequence);
							} else {
								HashSet<String> set = new HashSet<String>();
								set.add(sequence);
								genomeMap.put(decoyGenome, set);
							}

							if (proteinPepMap.containsKey(decoyPro)) {
								proteinPepMap.get(decoyPro).add(sequence);
							} else {
								HashSet<String> set = new HashSet<String>();
								set.add(sequence);
								proteinPepMap.put(decoyPro, set);
							}
						}
						tempPepSet.add(sequence);
						decoyPPCount++;
					}
				} else {
					missSet.add(sequence);
				}
			}
		}

		logger.info(getTaskName() + ": " + missSet.size()
				+ " peptide sequences were not found from the in-silico digested MAG database, "
				+ "maybe due to unspecific digestion, will be not be included in the following analysis");
		logger.info(getTaskName() + ": " + decoyPPCount
				+ " decoy precursor will be used for FDR control in the following steps.");
		System.out.println(new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss").format(new Date()) + "\t" + getTaskName() + ": "
				+ missSet.size() + " peptide sequences were not found from the in-silico digested MAG database, "
				+ "maybe due to unspecific digestion, will be not be included in the following analysis");

		ArrayList<DiaNNPrecursor> genomeFilteredList = new ArrayList<DiaNNPrecursor>();
		HashSet<String> temppepHashSet1 = new HashSet<String>();
		for (DiaNNPrecursor pp : diaPp) {
			String sequence = pp.getSeqString();
			// DIA Target workflow: No host file processing, so skip hostPepProMap check
			if (tempPepSet.contains(sequence)) {
				ArrayList<String> proList = new ArrayList<String>();
				if (pepProMap.containsKey(sequence)) {
					HashSet<String> pros = pepProMap.get(sequence);
					for (String pro : pros) {
						String genome = pro.split("_")[0];
						if (genomeMap.containsKey(genome) && proteinPepMap.containsKey(pro)) {
							proList.add(pro);
						}
					}
				} else {
					if (decoyPepMap.containsKey(sequence)) {
						String targetSeq = decoyPepMap.get(sequence);
						if (pepProMap.containsKey(targetSeq)) {
							HashSet<String> pros = pepProMap.get(targetSeq);
							for (String pro : pros) {
								String genome = pro.split("_")[0];
								if (genomeMap.containsKey(FastaManager.shuffle + genome)
										&& proteinPepMap.containsKey(FastaManager.shuffle + pro)) {
									proList.add(FastaManager.shuffle + pro);
								}
							}
						}
					}
				}

				if (proList.size() > 0) {
					pp.setProGroups(proList.toArray(String[]::new));
					genomeFilteredList.add(pp);
					temppepHashSet1.add(sequence);
				}
			}
		}

		DiaNNPrecursor[] filteredPPS = genomeFilteredList.toArray(DiaNNPrecursor[]::new);
		File genomeScoreFile = new File(combinedFile, "genomes_score.tsv");
		File proteinScoreFile = new File(combinedFile, "proteins_score.tsv");

		HashMap<String, double[]> genomeScoreMap = new HashMap<String, double[]>();
		HashMap<String, double[]> proteinScoreMap = new HashMap<String, double[]>();
		calGenomeProteinScore(filteredPPS, genomeScoreMap, proteinScoreMap, genomeScoreFile, proteinScoreFile);

		HashSet<String> tempPepHashSet = new HashSet<String>();
		ArrayList<DiaNNPrecursor> finalFilteredList = new ArrayList<DiaNNPrecursor>();
		for (DiaNNPrecursor pp : filteredPPS) {
			ArrayList<String> proList = new ArrayList<String>();
			String[] pros = pp.getProGroups();
			String razorProtein = "";
			double razorQvalue = 1.0;
			for (String pro : pros) {
				if (proteinScoreMap.containsKey(pro)) {
					proList.add(pro);
					if (proteinScoreMap.get(pro)[0] < razorQvalue) {
						razorQvalue = proteinScoreMap.get(pro)[0];
						razorProtein = pro;
					}
				}
			}

			if (proList.size() > 0) {
				tempPepHashSet.add(pp.getSeqString());
				pp.setProGroups(proList.toArray(String[]::new));
				pp.setRazorProtein(razorProtein);
				finalFilteredList.add(pp);
			}
		}

		this.diaPp = finalFilteredList.toArray(DiaNNPrecursor[]::new);
		DiannPrMatrixReader prReader = new DiannPrMatrixReader(this.quan_pep_file);
		FlashLfqQuanPeptide[] peps = prReader.getMetaPeptides(this.diaPp);
		this.quanExpNames = prReader.getIntensityTitle();

		HashMap<String, FlashLfqQuanPeptide> modseqPepMap = new HashMap<String, FlashLfqQuanPeptide>();
		for (FlashLfqQuanPeptide pep : peps) {
			modseqPepMap.put(pep.getModSeq(), pep);
		}

		logger.info(getTaskName() + ": protein inference finished, " + this.diaPp.length + " PSMs, "
				+ modseqPepMap.size() + " peptide, " + proteinScoreMap.size() + " proteins and " + genomeScoreMap.size()
				+ " genomes were kept");
		System.out.println(new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss").format(new Date()) + "\t" + getTaskName()
				+ ": protein inference finished, " + this.diaPp.length + " PSMs, " + modseqPepMap.size() + " peptide, "
				+ proteinScoreMap.size() + " proteins and " + genomeScoreMap.size() + " genomes were kept");

		logger.info(getTaskName() + ": exporting peptide report started");
		System.out.println(new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss").format(new Date()) + "\t" + getTaskName()
				+ ": exporting peptide report started");

		HashMap<String, double[]> genomeIntensityMap = new HashMap<String, double[]>();
		HashMap<String, double[]> proteinIntensityMap = new HashMap<String, double[]>();
		HashMap<String, Double> proteinTotalIntenMap = new HashMap<String, Double>();

		HashSet<String> pepSequenceSet = new HashSet<String>();
		HashMap<String, HashSet<String>> totalPepCountMap = new HashMap<String, HashSet<String>>();
		HashMap<String, HashSet<String>> razorPepCountMap = new HashMap<String, HashSet<String>>();
		HashMap<String, HashSet<String>> genomePepCountMap = new HashMap<String, HashSet<String>>();
		HashMap<String, HashSet<String>> genomeRazorPepCountMap = new HashMap<String, HashSet<String>>();
		ArrayList<FlashLfqQuanPeptide> peptideList = new ArrayList<FlashLfqQuanPeptide>();

		for (int i = 0; i < peps.length; i++) {
			String stripeSequence = peps[i].getSequence();
			String razorProtein = peps[i].getRazorProtein();

			pepSequenceSet.add(stripeSequence);
			peptideList.add(peps[i]);

			if (razorPepCountMap.containsKey(razorProtein)) {
				razorPepCountMap.get(razorProtein).add(stripeSequence);
			} else {
				HashSet<String> set = new HashSet<String>();
				set.add(stripeSequence);
				razorPepCountMap.put(razorProtein, set);
			}

			if (razorProtein.startsWith(magDb.getIdentifier())) {
				String genome = razorProtein.split("_")[0];
				if (razorPepCountMap.containsKey(genome)) {
					razorPepCountMap.get(genome).add(stripeSequence);
				} else {
					HashSet<String> set = new HashSet<String>();
					set.add(stripeSequence);
					razorPepCountMap.put(genome, set);
				}
			}
		}

		HashMap<String, Integer> psmCountMap = new HashMap<String, Integer>();

		int[] filePsmCount = new int[quanExpNames.length];
		HashSet<String>[] fileSequenceSets = new HashSet[quanExpNames.length];
		for (int i = 0; i < fileSequenceSets.length; i++) {
			fileSequenceSets[i] = new HashSet<String>();
		}
		
		for (int i = 0; i < peps.length; i++) {
			String[] proteins = peps[i].getProteins();
			int count = 0;
			for (int j = 0; j < proteins.length; j++) {
				if (razorPepCountMap.containsKey(proteins[j])) {
					count++;
				}
			}

			if (count > 0) {
				String modSequence = peps[i].getModSeq();
				int[] idenType = peps[i].getIdenType();
				for (int j = 0; j < fileSequenceSets.length; j++) {
					if (idenType[j] == 0) {
						fileSequenceSets[j].add(modSequence);
					}
				}

				String stripeSequence = peps[i].getSequence();
				double[] pepIntensity = new double[peps[i].getIntensity().length];
				for (int j = 0; j < pepIntensity.length; j++) {
					pepIntensity[j] = peps[i].getIntensity()[j] / (double) count;
				}

				int pepMs2Count = 0;
				int[] ms2Count = peps[i].getMs2Counts();
				for (int j = 0; j < ms2Count.length; j++) {
					pepMs2Count += ms2Count[j];
					filePsmCount[j] += ms2Count[j];
				}

				for (int j = 0; j < proteins.length; j++) {
					if (razorPepCountMap.containsKey(proteins[j])) {

						if (psmCountMap.containsKey(proteins[j])) {
							psmCountMap.put(proteins[j], psmCountMap.get(proteins[j]) + pepMs2Count);
						} else {
							psmCountMap.put(proteins[j], pepMs2Count);
						}

						double proteinTotalInten;
						double[] proteinIntensity;
						if (proteinIntensityMap.containsKey(proteins[j])) {
							proteinIntensity = proteinIntensityMap.get(proteins[j]);
							proteinTotalInten = proteinTotalIntenMap.get(proteins[j]);
						} else {
							proteinIntensity = new double[pepIntensity.length];
							proteinIntensityMap.put(proteins[j], proteinIntensity);
							proteinTotalInten = 0.0;
						}
						for (int k = 0; k < pepIntensity.length; k++) {
							proteinIntensity[k] += pepIntensity[k];
							proteinTotalInten += pepIntensity[k];
						}
						proteinTotalIntenMap.put(proteins[j], proteinTotalInten);

						if (totalPepCountMap.containsKey(proteins[j])) {
							totalPepCountMap.get(proteins[j]).add(stripeSequence);
						} else {
							HashSet<String> set = new HashSet<String>();
							set.add(stripeSequence);
							totalPepCountMap.put(proteins[j], set);
						}

						if (proteins[j].startsWith(magDb.getIdentifier())) {
							String genome = proteins[j].split("_")[0];
							if (genomeScoreMap.containsKey(genome)) {
								double[] genomeIntensity;
								if (genomeIntensityMap.containsKey(genome)) {
									genomeIntensity = genomeIntensityMap.get(genome);
								} else {
									genomeIntensity = new double[pepIntensity.length];
									genomeIntensityMap.put(genome, genomeIntensity);
								}
								for (int k = 0; k < pepIntensity.length; k++) {
									genomeIntensity[k] += pepIntensity[k];
								}

								if (psmCountMap.containsKey(genome)) {
									psmCountMap.put(genome, psmCountMap.get(genome) + pepMs2Count);
								} else {
									psmCountMap.put(genome, pepMs2Count);
								}

								if (totalPepCountMap.containsKey(genome)) {
									totalPepCountMap.get(genome).add(stripeSequence);
								} else {
									HashSet<String> set = new HashSet<String>();
									set.add(stripeSequence);
									totalPepCountMap.put(genome, set);
								}
							}
						}
					}
				}
			}
		}

		String[] proteins = proteinTotalIntenMap.keySet().toArray(String[]::new);
		Arrays.sort(proteins, new Comparator<String>() {

			@Override
			public int compare(String o1, String o2) {
				// TODO Auto-generated method stub
				if (proteinScoreMap.get(o1)[0] > proteinScoreMap.get(o2)[0]) {
					return 1;
				} else if (proteinScoreMap.get(o1)[0] < proteinScoreMap.get(o2)[0]) {
					return -1;
				} else {
					if (proteinTotalIntenMap.get(o1) > proteinTotalIntenMap.get(o2)) {
						return -1;
					} else if (proteinTotalIntenMap.get(o1) < proteinTotalIntenMap.get(o2)) {
						return 1;
					} else {
						return 0;
					}
				}
			}
		});

		HashMap<String, Integer> proteinIdMap = new HashMap<String, Integer>();
		for (int i = 0; i < proteins.length; i++) {
			proteinIdMap.put(proteins[i], i);
		}

		try (PrintWriter idWriter = new PrintWriter(new File(resultFolderFile, "id.tsv"))) {
			idWriter.println("Raw file\tUnique peptide count\tPSM count");
			for (int i = 0; i < quanExpNames.length; i++) {
				idWriter.println(quanExpNames[i] + "\t" + fileSequenceSets[i].size() + "\t" + filePsmCount[i]);
			}
			idWriter.close();
		} catch (IOException e) {
			System.out.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": error in writing identification information to " + new File(resultFolderFile, "id.tsv"));
			logger.error(getTaskName() + ": error in writing identification information to "
					+ new File(resultFolderFile, "id.tsv"), e);
		}

		updateProgress(85, 100);

		HashMap<String, ArrayList<Integer>> proPepIdMap = new HashMap<String, ArrayList<Integer>>();
		HashMap<String, ArrayList<Integer>> proPepRazorIdMap = new HashMap<String, ArrayList<Integer>>();
		try {
			this.final_pep_txt = new File(resultFolderFile, "final_peptides.tsv");

			PrintWriter writer = new PrintWriter(this.final_pep_txt);

			StringBuilder titlesb = new StringBuilder();
			titlesb.append("Sequence").append("\t");
			titlesb.append("Base sequence").append("\t");
			titlesb.append("Length").append("\t");
			titlesb.append("Missed cleavages").append("\t");
			titlesb.append("Mass").append("\t");
			titlesb.append("Proteins").append("\t");
			titlesb.append("Charges").append("\t");
			titlesb.append("Qvalue").append("\t");
			titlesb.append("PEP").append("\t");

			for (int i = 0; i < this.quanExpNames.length; i++) {
				titlesb.append("Identification type " + this.quanExpNames[i]).append("\t");
			}

			for (int i = 0; i < this.quanExpNames.length; i++) {
				titlesb.append("MS/MS count " + this.quanExpNames[i]).append("\t");
			}

			titlesb.append("Intensity").append("\t");

			for (int i = 0; i < quanExpNames.length; i++) {
				titlesb.append("Intensity " + this.quanExpNames[i]).append("\t");
			}

			titlesb.append("Reverse").append("\t");
			titlesb.append("Potential contaminant").append("\t");
			titlesb.append("id").append("\t");
			titlesb.append("Protein group IDs").append("\t");
			titlesb.append("MS/MS Count");

			writer.println(titlesb);

			int id = 0;
			for (FlashLfqQuanPeptide peptide : peptideList) {

				String fullSeq = peptide.getModSeq();
				String sequence = peptide.getSequence();

				StringBuilder sb = new StringBuilder();
				sb.append(fullSeq).append("\t");
				sb.append(sequence).append("\t");
				sb.append(sequence.length()).append("\t");
				sb.append(peptide.getMissCleave()).append("\t");
				sb.append(peptide.getMass()).append("\t");

				String[] pros = peptide.getProteins();

				for (String pro : pros) {
					if (proteinIdMap.containsKey(pro)) {
						sb.append(pro).append(";");
					}
				}
				sb.deleteCharAt(sb.length() - 1);
				sb.append("\t");

				int[] charge = peptide.getCharges();
				StringBuilder chargesb = new StringBuilder();
				for (int j = 0; j < charge.length; j++) {
					if (charge[j] == 1) {
						chargesb.append(j).append(";");
					}
				}
				if (chargesb.length() > 1) {
					sb.append(chargesb.subSequence(0, chargesb.length() - 1)).append("\t");
				} else {
					sb.append("2").append("\t");
				}

				sb.append(peptide.getScore()).append("\t");
				sb.append(peptide.getPEP()).append("\t");

				int[] idenType = peptide.getIdenType();
				for (int j = 0; j < idenType.length; j++) {
					sb.append(MetaPeptide.idenTypeStrings[idenType[j]]).append("\t");
				}

				int totalMsCount = 0;
				int[] ms2Count = peptide.getMs2Counts();
				for (int j = 0; j < ms2Count.length; j++) {
					totalMsCount += ms2Count[j];
					sb.append(ms2Count[j]).append("\t");
				}

				double[] intensity = peptide.getIntensity();
				double totalIntensity = 0;
				for (int j = 0; j < intensity.length; j++) {
					totalIntensity += intensity[j];
				}

				sb.append(BigInteger.valueOf((long) Math.floor(totalIntensity))).append("\t");

				for (int j = 0; j < intensity.length; j++) {
					sb.append(BigInteger.valueOf((long) Math.floor(intensity[j]))).append("\t");
				}

				sb.append("").append("\t");
				sb.append("").append("\t");

				sb.append(id).append("\t");

				for (String pro : pros) {
					if (proteinIdMap.containsKey(pro)) {
						sb.append(proteinIdMap.get(pro)).append(";");

						if (proPepIdMap.containsKey(pro)) {
							proPepIdMap.get(pro).add(id);
						} else {
							ArrayList<Integer> list = new ArrayList<Integer>();
							list.add(id);
							proPepIdMap.put(pro, list);
						}

						HashSet<String> rSet = razorPepCountMap.get(pro);
						if (proPepRazorIdMap.containsKey(pro)) {
							if (rSet.contains(sequence)) {
								proPepRazorIdMap.get(pro).add(1);
							} else {
								proPepRazorIdMap.get(pro).add(0);
							}
						} else {
							ArrayList<Integer> list = new ArrayList<Integer>();
							if (rSet.contains(sequence)) {
								list.add(1);
							} else {
								list.add(0);
							}
							proPepRazorIdMap.put(pro, list);
						}
					}
				}
				sb.deleteCharAt(sb.length() - 1);
				sb.append("\t");

				sb.append(totalMsCount);

				writer.println(sb);
				id++;
			}

			writer.close();

		} catch (Exception e) {
			// TODO Auto-generated catch block
			logger.error(getTaskName() + ": error in writing the peptides to  " + final_pep_txt, e);
			System.err.println(format.format(new Date()) + "\t" + getTaskName() + ": error in writing the peptides to  "
					+ final_pep_txt);
		}

		logger.info(getTaskName() + ": peptide report has been exported to " + final_pep_txt);
		System.out.println(
				format.format(new Date()) + "\t" + getTaskName() + ": peptide report has been exported to " + final_pep_txt);

		if (!final_pep_txt.exists() || final_pep_txt.length() == 0) {
			return false;
		}

		updateProgress(87, 100);

		logger.info(getTaskName() + ": exporting protein report started");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": exporting protein report started");

		ArrayList<MetaProtein> proteinList = new ArrayList<MetaProtein>();
		try {

			this.final_pro_txt = new File(resultFolderFile, "final_proteins.tsv");

			PrintWriter writer = new PrintWriter(this.final_pro_txt);

			StringBuilder titlesb = new StringBuilder();
			titlesb.append("Protein IDs").append("\t");
			titlesb.append("Peptide counts (all)").append("\t");
			titlesb.append("Peptide counts (razor)").append("\t");
			titlesb.append("Number of PSMs").append("\t");
			titlesb.append("Qvalue").append("\t");
			titlesb.append("Local FDR").append("\t");
			titlesb.append("Intensity").append("\t");
			for (int i = 0; i < quanExpNames.length; i++) {
				titlesb.append("Intensity " + quanExpNames[i]).append("\t");
			}

			titlesb.append("Reverse").append("\t");
			titlesb.append("Potential contaminant").append("\t");
			titlesb.append("id").append("\t");
			titlesb.append("Peptide IDs").append("\t");
			titlesb.append("Peptide is razor");

			writer.println(titlesb);

			DecimalFormat df = FormatTool.getDFE4();
			DecimalFormat df4 = FormatTool.getDF4();
			for (int i = 0; i < proteins.length; i++) {

				StringBuilder sb = new StringBuilder();
				sb.append(proteins[i]).append("\t");
				sb.append(totalPepCountMap.get(proteins[i]).size()).append("\t");
				sb.append(razorPepCountMap.get(proteins[i]).size()).append("\t");
				sb.append(psmCountMap.get(proteins[i])).append("\t");

				sb.append(df.format(proteinScoreMap.get(proteins[i])[0])).append("\t");
				sb.append(df4.format(proteinScoreMap.get(proteins[i])[1])).append("\t");

				double totalIntensity = proteinTotalIntenMap.get(proteins[i]);
				sb.append(df.format(totalIntensity)).append("\t");

				double[] intensity = proteinIntensityMap.get(proteins[i]);
				for (int j = 0; j < intensity.length; j++) {
					sb.append(df.format(intensity[j])).append("\t");
				}
				sb.append("\t");
				sb.append("\t");
				sb.append(i).append("\t");

				ArrayList<Integer> pepIdList = proPepIdMap.get(proteins[i]);
				ArrayList<Integer> pepRazorList = proPepRazorIdMap.get(proteins[i]);

				for (int j = 0; j < pepIdList.size(); j++) {
					sb.append(pepIdList.get(j)).append(";");
				}
				sb.deleteCharAt(sb.length() - 1);
				sb.append("\t");

				for (int j = 0; j < pepRazorList.size(); j++) {
					if (pepRazorList.get(j) == 1) {
						sb.append("true;");
					} else {
						sb.append("false;");
					}
				}
				sb.deleteCharAt(sb.length() - 1);
				writer.println(sb);

				MetaProtein mp = new MetaProtein(i, 0, proteins[i], totalPepCountMap.get(proteins[i]).size(),
						psmCountMap.get(proteins[i]), proteinScoreMap.get(proteins[i])[0], intensity);
				proteinList.add(mp);
			}

			writer.close();

		} catch (IOException e) {
			// TODO Auto-generated catch block
			logger.error(getTaskName() + ": error in writing final protein result to " + this.final_pro_txt.getName(), e);
			System.err.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": error in writing final protein result to " + this.final_pro_txt.getName());
		}

		logger.info(getTaskName() + ": protein report has been exported to " + final_pro_txt.getName());
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": protein report has been exported to "
				+ final_pro_txt.getName());

		updateProgress(89, 100);	

		logger.info(getTaskName() + ": functional annotation started");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": functional annotation started");

		File funcFile = new File(resultFolderFile, "functional_annotation");
		if (!funcFile.exists()) {
			funcFile.mkdir();
		}

		File funTsv = new File(funcFile, "functions.tsv");
		final_pro_xml_file = new File(funcFile, "functions.xml");

		HashMap<String, Integer> proTaxIdMap = new HashMap<String, Integer>();
		MetaProtein[] metapros = proteinList.toArray(new MetaProtein[proteinList.size()]);
		MetaProteinAnnoMag[] mpas = funcAnnoSqlite(metapros, quanExpNames);
		for (int i = 0; i < mpas.length; i++) {
			String name = mpas[i].getPro().getName();
			if (NumberUtils.isCreatable(mpas[i].getTax_id())) {
				proTaxIdMap.put(name, Integer.parseInt(mpas[i].getTax_id()));
			} else {
				proTaxIdMap.put(name, 131567);
			}
		}

		MetaProteinXMLReader2 funReader = null;
		try {
			funReader = new MetaProteinXMLReader2(final_pro_xml_file);
		} catch (DocumentException e) {
			// TODO Auto-generated catch block
			logger.error(getTaskName() + ": error in reading protein function information from "
					+ this.final_pro_xml_file.getName(), e);
			System.err.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": error in reading protein function information from " + this.final_pro_xml_file.getName());
			return false;
		}

		if (!funTsv.exists()) {
			System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": writing results to " + funTsv.getName()
					+ " started");
			try {
				funReader.exportTsv(funTsv);				
				// Export separate functional category reports
				File goReportTsv = new File(funcFile, "Gene_Ontology_id_report.tsv");
				File ecReportTsv = new File(funcFile, "EC_id_report.tsv");
				File keggReportTsv = new File(funcFile, "KEGG_ko_report.tsv");
				File briteReportTsv = new File(funcFile, "BRITE_report.tsv");
				File cogReportTsv = new File(funcFile, "COG_accession_report.tsv");
				File nogReportTsv = new File(funcFile, "NOG_accession_report.tsv");
				File cogCatTsv = new File(funcFile, "COG_category_report.tsv");
				File nogCatTsv = new File(funcFile, "NOG_category_report.tsv");
				
				funReader.exportGOReport(goReportTsv);
				funReader.exportECReport(ecReportTsv);
				funReader.exportKeggReport(keggReportTsv);
				funReader.exportBriteReport(briteReportTsv);
				funReader.exportCOGReport(cogReportTsv);
				funReader.exportNOGReport(nogReportTsv);
				funReader.exportCOGCategoryReport(cogCatTsv);
				funReader.exportNOGCategoryReport(nogCatTsv);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				logger.error(getTaskName() + ": error in writing protein function information to " + funTsv.getName(), e);
				System.err.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": error in writing protein function information to " + funTsv.getName());
			}
			System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": writing results to " + funTsv.getName()
					+ " finished");
		}

		File html = new File(funcFile, "functions.html");
		if (!html.exists()) {
			System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": writing results to "
					+ html.getName() + " started");
			try {
				funReader.exportHtml(html);

				// Copy the generated functions.html to the reports directory
				File reportFunctionsHtml = new File(metaPar.getReportDir(), "reports/functions.html");
				if (reportFunctionsHtml.getParentFile() != null && !reportFunctionsHtml.getParentFile().exists()) {
					reportFunctionsHtml.getParentFile().mkdirs();
				}
				Files.copy(html.toPath(), reportFunctionsHtml.toPath(), StandardCopyOption.REPLACE_EXISTING);
				logger.info("Copied functions.html to reports directory: " + reportFunctionsHtml.getAbsolutePath());

			} catch (IOException e) {
				// TODO Auto-generated catch block
				logger.error(getTaskName() + ": error in writing protein function information to " + html.getName(), e);
				System.err.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": error in writing protein function information to " + html.getName());
			}
			System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": writing results to "
					+ html.getName() + " finished");
		}
		
		// Copy functions.html to report directory for navigation access
		if (html.exists() && this.reportHtmlDir != null) {
			File reportFunctionsHtml = new File(this.reportHtmlDir, "functions.html");
			if (!reportFunctionsHtml.exists()) {
				try {
					java.nio.file.Files.copy(html.toPath(), reportFunctionsHtml.toPath(), 
						java.nio.file.StandardCopyOption.REPLACE_EXISTING);
					logger.info(getTaskName() + ": copied functions.html to report directory");
					System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": copied functions.html to report directory");
				} catch (IOException e) {
					logger.error(getTaskName() + ": error copying functions.html to report directory", e);
					System.err.println(format.format(new Date()) + "\t" + getTaskName() + ": error copying functions.html to report directory");
				}
			}
		}

		updateProgress(91, 100);

		// taxonomy analysis

		logger.info(getTaskName() + ": taxonomy analysis started");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": taxonomy analysis started");

		File taxFile = new File(resultFolderFile, "taxonomy_analysis");
		if (!taxFile.exists()) {
			taxFile.mkdir();
		}

		MagFuncSearcher searcher = null;
		try {
			searcher = new MagFuncSearcher(metaPar.getMicroDb(), mrp.getToolPath(ToolType.FUNCDB));

		} catch (NumberFormatException | SQLException | IOException e) {
			// TODO Auto-generated catch block
			logger.error(getTaskName() + ": error in genome analysis, task failed", e);
			System.err.println(format.format(new Date()) + "\t" + getTaskName() + ": error in genome analysis, task failed");

			return false;
		}

		String[] genomes = genomeIntensityMap.keySet().toArray(new String[genomeIntensityMap.size()]);
		HashMap<String, String[]> genomeTaxaMap = null;
		try {
			genomeTaxaMap = searcher.matchGenome(genomes);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			logger.error(getTaskName() + ": Error in searching the genome2taxa database", e);
			System.err.println(
					format.format(new Date()) + "\t" + getTaskName() + ": error in searching the genome2taxa database");
			return false;
		}

		File genomeHtmlFile = new File(taxFile, "Genome.html");
		try {
			HMGenomeHtmlWriter htmlWriter = new HMGenomeHtmlWriter(genomeHtmlFile);
			htmlWriter.write(genomeTaxaMap);
			htmlWriter.close();

			// Copy the generated Genome.html to the reports directory
			File reportGenomeHtml = new File(metaPar.getReportDir(), "reports/Genome.html");
			if (reportGenomeHtml.getParentFile() != null && !reportGenomeHtml.getParentFile().exists()) {
				reportGenomeHtml.getParentFile().mkdirs();
			}
			Files.copy(genomeHtmlFile.toPath(), reportGenomeHtml.toPath(), StandardCopyOption.REPLACE_EXISTING);
			logger.info("Copied Genome.html to reports directory: " + reportGenomeHtml.getAbsolutePath());
		} catch (IOException e) {
			// TODO Auto-generated catch block
			logger.error("Error in writing genome information to " + genomeHtmlFile, e);
			System.err.println(
					format.format(new Date()) + "\t" + "Error in writing genome information to " + genomeHtmlFile);
			return false;
		}
		
		// Copy Genome.html to report directory for navigation access
		if (genomeHtmlFile.exists() && this.reportHtmlDir != null) {
			File reportGenomeHtml = new File(this.reportHtmlDir, "Genome.html");
			if (!reportGenomeHtml.exists()) {
				try {
					java.nio.file.Files.copy(genomeHtmlFile.toPath(), reportGenomeHtml.toPath(), 
						java.nio.file.StandardCopyOption.REPLACE_EXISTING);
					logger.info(getTaskName() + ": copied Genome.html to report directory");
					System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": copied Genome.html to report directory");
				} catch (IOException e) {
					logger.error(getTaskName() + ": error copying Genome.html to report directory", e);
					System.err.println(format.format(new Date()) + "\t" + getTaskName() + ": error copying Genome.html to report directory");
				}
			}
		}

		TaxonomyRanks[] mainRanks = TaxonomyRanks.getMainRanks();
		double[] cellularIntensity = new double[quanExpNames.length];
		HashMap<String, String[]>[] taxRankMaps = new HashMap[mainRanks.length];
		HashMap<String, double[]>[] taxIntensityMaps = new HashMap[mainRanks.length];
		HashMap<String, Double>[] taxScoreMaps = new HashMap[mainRanks.length];
		for (int i = 0; i < taxRankMaps.length; i++) {
			taxRankMaps[i] = new HashMap<String, String[]>();
			taxIntensityMaps[i] = new HashMap<String, double[]>();
			taxScoreMaps[i] = new HashMap<String, Double>();
		}

		File genomeFile = new File(taxFile, "Genome.tsv");
		PrintWriter genomeWriter = null;
		try {
			genomeWriter = new PrintWriter(genomeFile);
			StringBuilder genomeTitlesb = new StringBuilder();
			genomeTitlesb.append("Genome\t");
			genomeTitlesb.append("p-value\t");
			genomeTitlesb.append("Local FDR\t");
			genomeTitlesb.append("PSM count\t");
			genomeTitlesb.append("Peptide count\t");
			genomeTitlesb.append("Razor peptide count\t");
			genomeTitlesb.append("GTDB_Name\t");
			genomeTitlesb.append("GTDB_Rank\t");

			for (int i = 0; i < mainRanks.length; i++) {
				genomeTitlesb.append("GTDB_").append(mainRanks[i].getName()).append("\t");
			}

			genomeTitlesb.append("NCBI_Name\t");
			genomeTitlesb.append("NCBI_Rank");

			for (int i = 0; i < mainRanks.length; i++) {
				genomeTitlesb.append("\t").append("NCBI_").append(mainRanks[i].getName());
			}

			for (String exp : expNames) {
				if (exp.startsWith("Intensity ")) {
					genomeTitlesb.append("\t").append(exp);
				} else {
					genomeTitlesb.append("\tIntensity ").append(exp);
				}
			}

			genomeWriter.println(genomeTitlesb);

			updateProgress(94, 100);

			for (int i = 0; i < genomes.length; i++) {
				String[] taxaContent = genomeTaxaMap.get(genomes[i]);
				if (genomeIntensityMap.containsKey(genomes[i])) {
					double[] genomeIntensity = genomeIntensityMap.get(genomes[i]);
					for (int j = 0; j < genomeIntensity.length; j++) {
						cellularIntensity[j] += genomeIntensity[j];
					}

					StringBuilder genomeTaxaSb = new StringBuilder();
					genomeTaxaSb.append(genomes[i]).append("\t");

					genomeTaxaSb.append(FormatTool.getDFE4().format(genomeScoreMap.get(genomes[i])[0])).append("\t");
					genomeTaxaSb.append(FormatTool.getDF4().format(genomeScoreMap.get(genomes[i])[1])).append("\t");
					
					if (psmCountMap.containsKey(genomes[i])) {
						genomeTaxaSb.append(psmCountMap.get(genomes[i])).append("\t");
					} else {
						genomeTaxaSb.append("0").append("\t");
					}
					if (totalPepCountMap.containsKey(genomes[i])) {
						genomeTaxaSb.append(totalPepCountMap.get(genomes[i]).size()).append("\t");
					} else {
						genomeTaxaSb.append("0").append("\t");
					}
					if (razorPepCountMap.containsKey(genomes[i])) {
						genomeTaxaSb.append(razorPepCountMap.get(genomes[i]).size());
					} else {
						genomeTaxaSb.append("0");
					}

					for (int j = 1; j < taxaContent.length; j++) {
						genomeTaxaSb.append("\t").append(taxaContent[j]);
					}

					for (int j = 0; j < genomeIntensity.length; j++) {
						genomeTaxaSb.append("\t").append(genomeIntensity[j]);
					}
					genomeWriter.println(genomeTaxaSb);

					for (int j = 0; j < taxRankMaps.length; j++) {
						String taxon = taxaContent[j + 3];
						if (taxon == null || taxon.length() == 0) {
							continue;
						}
						if (!taxRankMaps[j].containsKey(taxon)) {
							String[] lineage = new String[j + 1];
							System.arraycopy(taxaContent, 3, lineage, 0, lineage.length);
							taxRankMaps[j].put(taxon, lineage);
						}

						double[] taxIntensity;
						if (taxIntensityMaps[j].containsKey(taxon)) {
							taxIntensity = taxIntensityMaps[j].get(taxon);
						} else {
							taxIntensity = new double[genomeIntensity.length];
						}

						for (int k = 0; k < taxIntensity.length; k++) {
							taxIntensity[k] += genomeIntensity[k];
						}
						taxIntensityMaps[j].put(taxon, taxIntensity);

						if (taxScoreMaps[j].containsKey(taxon)) {
							taxScoreMaps[j].put(taxon, taxScoreMaps[j].get(taxon) * genomeScoreMap.get(genomes[i])[0]);
						} else {
							taxScoreMaps[j].put(taxon, genomeScoreMap.get(genomes[i])[0]);
						}
					}
				}
			}
			genomeWriter.close();

		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			logger.error("Error in writing genome information to " + genomeFile, e);
			System.err
					.println(format.format(new Date()) + "\t" + "Error in writing genome information to " + genomeFile);
			return false;
		}

		File taxaFile = new File(taxFile, "Taxa.tsv");
		File taxaReportFile = new File(taxFile, "Taxa_report.tsv");
		PrintWriter taxaWriter = null;
		PrintWriter taxaReportWriter = null;
		try {
			taxaWriter = new PrintWriter(taxaFile);
			taxaReportWriter = new PrintWriter(taxaReportFile);

			StringBuilder taxaTitlesb1 = new StringBuilder();
			taxaTitlesb1.append("Name\t");
			taxaTitlesb1.append("Score\t");
			taxaTitlesb1.append("Rank");

			for (int i = 0; i < mainRanks.length; i++) {
				taxaTitlesb1.append("\t").append(mainRanks[i].getName());
			}
			for (String exp : expNames) {
				if (exp.startsWith("Intensity ")) {
					taxaTitlesb1.append("\t").append(exp);
				} else {
					taxaTitlesb1.append("\tIntensity ").append(exp);
				}
			}
			StringBuilder taxaTitlesb2 = new StringBuilder();
			taxaTitlesb2.append("Name\t");
			taxaTitlesb2.append("Rank");

			for (int i = 0; i < mainRanks.length; i++) {
				taxaTitlesb2.append("\t").append(mainRanks[i].getName());
			}
			for (String exp : expNames) {
				if (exp.startsWith("Intensity ")) {
					taxaTitlesb2.append("\t").append(exp);
				} else {
					taxaTitlesb2.append("\tIntensity ").append(exp);
				}
			}
			taxaWriter.println(taxaTitlesb1);
			taxaReportWriter.println(taxaTitlesb2);

			for (int i = 0; i < taxRankMaps.length; i++) {
				for (String taxon : taxRankMaps[i].keySet()) {
					String[] taxaContent = taxRankMaps[i].get(taxon);
					double[] intensity = taxIntensityMaps[i].get(taxon);

					StringBuilder taxaSb = new StringBuilder();
					taxaSb.append(taxon).append("\t");
					
					double taxonScore = taxScoreMaps[i].get(taxon);
					if (taxonScore == 0) {
						taxaSb.append("1000").append("\t");
					} else {
						taxaSb.append(FormatTool.getDF4().format(-Math.log10(taxonScore))).append("\t");
					}

					taxaSb.append(mainRanks[i].getName());

					for (int j = 0; j < mainRanks.length; j++) {
						if (j < taxaContent.length) {
							taxaSb.append("\t").append(taxaContent[j]);
						} else {
							taxaSb.append("\t");
						}
					}

					for (int j = 0; j < intensity.length; j++) {
						taxaSb.append("\t").append(intensity[j]);
					}

					taxaWriter.println(taxaSb);
				}
			}

			taxaWriter.close();

			for (int i = 0; i < taxRankMaps.length; i++) {
				for (String taxon : taxRankMaps[i].keySet()) {
					String[] taxaContent = taxRankMaps[i].get(taxon);
					double[] intensity = taxIntensityMaps[i].get(taxon);

					StringBuilder taxaSb = new StringBuilder();
					taxaSb.append(taxon).append("\t");
					taxaSb.append(mainRanks[i].getName());

					for (int j = 0; j < mainRanks.length; j++) {
						if (j < taxaContent.length) {
							taxaSb.append("\t").append(taxaContent[j]);
						} else {
							taxaSb.append("\t");
						}
					}
					for (int j = 0; j < intensity.length; j++) {
						taxaSb.append("\t").append(intensity[j]);
					}

					taxaReportWriter.println(taxaSb);
				}
			}

			taxaReportWriter.close();

		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			logger.error("Error in writing taxa information to " + taxaFile, e);
			System.err.println(format.format(new Date()) + "\t" + "Error in writing taxa information to " + taxaFile);
		}

		logger.info(getTaskName() + ": taxonomy analysis finished");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": taxonomy analysis finished");

		updateProgress(95, 100);

		updateMessage(getTaskName() + ": exporting report");

		logger.info(getTaskName() + ": exporting report started");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": exporting report started");

		File reportDir = new File(resultFolderFile, "report");
		if (!reportDir.exists()) {
			reportDir.mkdirs();
		}

		// Copy report templates using MetaPythonReportCopyTask (same as DDA workflow)
		MetaPythonReportCopyTask reportCopyTask = new MetaPythonReportCopyTask(resultFolderFile);
		Thread copyThread = new Thread(reportCopyTask);
		copyThread.start();
		try {
			copyThread.join();
			if (!reportCopyTask.get()) {
				logger.error(getTaskName() + ": error in copying report template");
				System.err.println(format.format(new Date()) + "\t" + getTaskName() + ": error in copying report template");
				return false;
			}
		} catch (Exception e) {
			logger.error(getTaskName() + ": error in copying report template", e);
			System.err.println(format.format(new Date()) + "\t" + getTaskName() + ": error in copying report template: " + e.getMessage());
			return false;
		}

		metaPar.setReportDir(reportDir);

		this.reportHtmlDir = new File(reportDir + "\\reports");

		if (!reportHtmlDir.exists()) {
			reportHtmlDir.mkdirs();
		}
		
		File jsOutput = new File(metaPar.getReportDir() + "\\reports\\tree.js");
		if (!jsOutput.getParentFile().exists()) {
			jsOutput.getParentFile().mkdirs();
		}

		if (!jsOutput.exists() || jsOutput.length() == 0) {
			MetaTreeHandler.exportJS(taxRankMaps, taxIntensityMaps, cellularIntensity, expNames, jsOutput);
		}
		
		// Generate bubble chart report from tree.js data (optional - don't fail if this doesn't work)
		File bubbleChartOutput = new File(metaPar.getReportDir() + "\\reports\\bubble_chart_report.html");
		if (!bubbleChartOutput.exists() || bubbleChartOutput.length() == 0) {
			try {
				generateBubbleChartReport(jsOutput, bubbleChartOutput);
			} catch (Exception e) {
				logger.warn(getTaskName() + ": Bubble chart generation failed, continuing with other reports: " + e.getMessage());
				System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": Bubble chart generation failed, continuing with other reports");
			}
		}

		MetaPythonReportTask task = new MetaPythonReportTask(metaPar.getThreadCount(), mrp);
		task.setWorkflowType("DIA");
		task.initializeManifestGenerator(new File(metaPar.getReportDir(), "reports"));
		
		// Create initial manifest file to ensure index page can load
		ReportManifestGenerator manifestGen = new ReportManifestGenerator(new File(metaPar.getReportDir(), "reports"));
		manifestGen.createInitialManifest();
		
		this.summaryMetaFile = new File(this.resultFolderFile, "metainfo.tsv");
		try {
			metadata.exportMetadata(summaryMetaFile);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		this.exportSummary(task);

		if (!final_pep_txt.exists() || final_pep_txt.length() == 0) {
			return false;
		} else {
			File report_peptides_summary = new File(reportHtmlDir, "report_peptides_summary.html");
			if (!report_peptides_summary.exists() || report_peptides_summary.length() == 0) {
				if (summaryMetaFile.exists()) {
					task.addTask(MetaPythonReportTask.PEPTIDE, final_pep_txt, report_peptides_summary, summaryMetaFile);
				} else {
					task.addTask(MetaPythonReportTask.PEPTIDE, final_pep_txt, report_peptides_summary);
				}
			}
		}

		if (!final_pro_txt.exists() || final_pro_txt.length() == 0) {
			return false;
		} else {
			File report_proteinGroups_summary = new File(reportHtmlDir, "report_proteinGroups_summary.html");
			if (!report_proteinGroups_summary.exists() || report_proteinGroups_summary.length() == 0) {
				if (summaryMetaFile.exists()) {
					task.addTask(MetaPythonReportTask.PROTEIN, final_pro_txt, report_proteinGroups_summary, summaryMetaFile);
				} else {
					task.addTask(MetaPythonReportTask.PROTEIN, final_pro_txt, report_proteinGroups_summary);
				}
			}
		}

		if (!taxaReportFile.exists() || taxaReportFile.length() == 0) {
			return false;
		} else {
			File report_taxonomy_summary = new File(this.reportHtmlDir, "report_taxonomy_summary.html");
			if (!report_taxonomy_summary.exists() || report_taxonomy_summary.length() == 0) {

				if (summaryMetaFile != null && summaryMetaFile.exists()) {
					task.addTask(MetaPythonReportTask.TAXONOMY, taxaReportFile, report_taxonomy_summary, summaryMetaFile);
				} else {
					task.addTask(MetaPythonReportTask.TAXONOMY, taxaReportFile, report_taxonomy_summary);
				}
			}
		}

		if (!funcFile.exists() || funcFile.listFiles().length == 0) {
			return false;
		} else {
			File report_function_summary = new File(this.reportHtmlDir, "report_function_summary.html");
			if (!report_function_summary.exists() || report_function_summary.length() == 0) {

				if (summaryMetaFile != null && summaryMetaFile.exists()) {
					task.addTask(MetaPythonReportTask.FUNCTION, funcFile, report_function_summary, summaryMetaFile);
				} else {
					task.addTask(MetaPythonReportTask.FUNCTION, funcFile, report_function_summary);
				}
			}
		}

		task.run();

		try {
			boolean finish = task.get();
			if (finish) {

				logger.info(getTaskName() + ": report generation successes");
				System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": report generation successes");

				updateProgress(100, 100);
				return true;
			}

		} catch (InterruptedException | ExecutionException e) {
			// TODO Auto-generated catch block
			logger.error(getTaskName() + ": error in the report generation", e);
			System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": error in the report generation");
		}

		return true;
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

		ConcurrentHashMap<String, HashMap<String, HashSet<String>>> genomePepMap = new ConcurrentHashMap<String, HashMap<String, HashSet<String>>>();
		HashSet<String> matchedPepSet = new HashSet<String>();
		ExecutorService executor = Executors.newFixedThreadPool(threadCount);

		// Get the MAG database current file (contains predicted peptides)
		File currentFile = super.magDb.getCurrentFile();
		File predictPepFile = new File(currentFile, "predicted");

		if (!predictPepFile.exists()) {
			logger.error("Predicted peptide folder not found: {}", predictPepFile.getAbsolutePath());
			return genomePepMap;
		}

		File[] allFiles = predictPepFile.listFiles();
		if (allFiles == null || allFiles.length == 0) {
			logger.error("No predicted peptide files found in: {}", predictPepFile.getAbsolutePath());
			return genomePepMap;
		}

		// Filter files based on target library genome information
		File[] filesToProcess;
		if (targetLibraryInfo.hasGenomeFiltering()) {
			// Only process files for genomes that were included in target library
			filesToProcess = new File[0];
			int filteredCount = 0;
			for (File file : allFiles) {
				if (file.getName().endsWith(".predicted.tsv")) {
					String genomeName = file.getName().substring(0,
							file.getName().length() - ".predicted.tsv".length());
					if (targetLibraryInfo.getGenomes().contains(genomeName)) {
						// Resize array and add file
						File[] newFiles = new File[filesToProcess.length + 1];
						System.arraycopy(filesToProcess, 0, newFiles, 0, filesToProcess.length);
						newFiles[filesToProcess.length] = file;
						filesToProcess = newFiles;
						filteredCount++;
					}
				}
			}
			logger.info("Filtered to {} predicted peptide files (from {} total) based on target library genomes",
					filteredCount, allFiles.length);
		} else {
			// No genome filtering - process all files
			filesToProcess = allFiles;
			logger.info("No genome filtering applied - processing all {} predicted peptide files", allFiles.length);
		}

		for (File file : filesToProcess) {
			executor.execute(() -> {
				if (file.getName().endsWith(".predicted.tsv")) {
					try (BufferedReader reader = new BufferedReader(new FileReader(file))) {
						HashMap<String, HashSet<String>> proPepMap = new HashMap<String, HashSet<String>>();
						String pline = reader.readLine(); // Skip header

						while ((pline = reader.readLine()) != null) {
							String[] cs = pline.split("\t");
							if (cs.length >= 3 && allPepSet.contains(cs[1])) {
								// cs[0] = protein_name, cs[1] = peptide_sequence, cs[2] = detectability_score

								// Apply protein filtering if enabled
								if (targetLibraryInfo.hasProteinFiltering()
										&& !targetLibraryInfo.getProteins().contains(cs[0])) {
									continue; // Skip this protein if not in target library
								}

								matchedPepSet.add(cs[1]);

								if (proPepMap.containsKey(cs[0])) {
									proPepMap.get(cs[0]).add(cs[1]);
								} else {
									HashSet<String> pepSet = new HashSet<String>();
									pepSet.add(cs[1]);
									proPepMap.put(cs[0], pepSet);
								}
							}
						}

						if (proPepMap.size() > 0) {
							String preFileName = file.getName();
							String genome = preFileName.substring(0, preFileName.length() - ".predicted.tsv".length());
							genomePepMap.put(genome, proPepMap);

							logger.debug("Mapped {} proteins and {} peptides for genome: {}", proPepMap.size(),
									proPepMap.values().stream().mapToInt(HashSet::size).sum(), genome);
						}
					} catch (IOException e) {
						logger.error(getTaskName() + ": error in reading predicted peptide file " + file, e);
						System.err.println(new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss").format(new Date()) + "\t"
								+ getTaskName() + ": error in reading " + file);
					}
				}
			});
		}

		executor.shutdown();

		try {
			if (executor.awaitTermination(1, TimeUnit.HOURS)) {
				// Update the allPepSet to only include matched peptides
				allPepSet.clear();
				allPepSet.addAll(matchedPepSet);

				logger.info("DIA Target peptide mapping completed:");
				logger.info("  - Matched peptides: {}", allPepSet.size());
				logger.info("  - Mapped genomes: {}", genomePepMap.size());
				logger.info("  - Total proteins mapped: {}",
						genomePepMap.values().stream().mapToInt(HashMap::size).sum());

				System.out.println(new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss").format(new Date()) + "\t" + getTaskName()
						+ ": DIA Target peptide mapping finished, " + allPepSet.size() + " peptides from "
						+ genomePepMap.size() + " genomes were matched");

			} else {
				logger.error(getTaskName() + ": peptide mapping timed out after 1 hour");
				System.err.println(new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss").format(new Date()) + "\t" + getTaskName()
						+ ": peptide mapping timed out");
				return null;
			}
		} catch (InterruptedException e) {
			logger.error(getTaskName() + ": peptide mapping was interrupted", e);
			System.err.println(new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss").format(new Date()) + "\t" + getTaskName()
					+ ": peptide mapping was interrupted");
			return null;
		}

		return genomePepMap;
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
