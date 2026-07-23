package uk.ac.quadram.metapilot.task.mag;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.lang.management.ManagementFactory;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import uk.ac.quadram.metapilot.task.io.ReportManifestGenerator;
import uk.ac.quadram.metapilot.task.MetaPythonReportTask;
import java.sql.SQLException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutionException;

import org.apache.commons.math3.stat.inference.TTest;

import com.sun.management.OperatingSystemMXBean;

import uk.ac.quadram.metapilot.core.mag.DdaMagGenomeScoringService;
import uk.ac.quadram.metapilot.core.mag.DdaMagPostAnalysisService;
import uk.ac.quadram.metapilot.core.mag.DdaMagProteinInferenceService;
import uk.ac.quadram.metapilot.core.mag.MagPostAnalysisScoringService;
import uk.ac.quadram.metapilot.core.mag.MagResultLayoutUtil;
import uk.ac.quadram.metapilot.core.math.MathTool;
import uk.ac.quadram.metapilot.core.prodb.FastaManager;
import uk.ac.quadram.metapilot.core.runner.DdaMagReportExportService;
import uk.ac.quadram.metapilot.core.runner.MagFunctionalAnnotationService;
import uk.ac.quadram.metapilot.core.runner.MagTaxonomyAnalysisService;
import uk.ac.quadram.metapilot.core.taxonomy.TaxonomyRanks;
import uk.ac.quadram.metapilot.core.tools.FormatTool;
import uk.ac.quadram.metapilot.dbSearch.pfind.PFindResultReader;
import javafx.scene.control.ProgressBar;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters.ToolType;
import uk.ac.quadram.metapilot.task.MetaPilotTask;
import uk.ac.quadram.metapilot.task.MetaAbstractTask;
import uk.ac.quadram.metapilot.task.io.MetaTreeHandler;
import uk.ac.quadram.metapilot.task.io.pep.MetaPeptide;
import uk.ac.quadram.metapilot.task.io.pep.MetaPSM;
import uk.ac.quadram.metapilot.task.io.pro.MetaProtein;
import uk.ac.quadram.metapilot.task.par.MetaBaseParameter;
import uk.ac.quadram.metapilot.task.par.MetaConstants;
import uk.ac.quadram.metapilot.task.par.MetaData;

public abstract class MetaPilotMagTask extends MetaPilotTask {

	// Logger is inherited from MetaPilotTask as an instance field

	protected String[] rawFiles;
	protected String[] rawFileNames;
	protected String[] quanExpNames;
	protected File spectraDirFile;
	protected MagDbItem magDb;
	protected File[] separateResultFolders;
	protected int threadCount;
	protected double psmCoverThreshold = 0.95;

	protected File rawDir;
	protected File quan_pep_file;
	protected File quan_pro_file;
	protected File final_pep_txt;
	protected File final_pro_txt;
	protected HashMap<String, String> expNameMap;
	protected File summaryMetaFile;
	protected File combinedFile;

	protected int maxTaskCount;

	/** When non-null, DIA-style genome scoring was used: map genome -> double[]{S, p_emp, q_value} for taxonomy_analysis outputs. */
	protected HashMap<String, double[]> lastDiaStyleGenomeStats;

	/** When non-null, holds runinfo data for genome_rank_scores_runinfo.txt (DIA-style scoring only). */
	protected DdaMagGenomeScoringService.DdaGenomeScoreRuninfo lastDiaStyleRuninfo;

	public MetaPilotMagTask(MetaBaseParameter metaPar, MetaResourceParameters mrp, ProgressBar progressBar,
			ProgressBar secondBar) {
		super(metaPar, mrp, progressBar, secondBar);
		this.threadCount = metaPar.getThreadCount();

		setTaskCount();
		logger.info("Task count set to: {}", maxTaskCount);

		logger.info(getTaskName() + ": Constructor completed, report templates will be copied during execution");
	}

	protected void setTaskCount() {
		long totalMemorySize = ((OperatingSystemMXBean) ManagementFactory.getOperatingSystemMXBean())
				.getTotalMemorySize() / 1024 / 1024 / 1024;
		this.maxTaskCount = (int) (totalMemorySize / 20);
		if (this.maxTaskCount == 0) {
			this.maxTaskCount = 1;
		}

		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": the number of max memory is "
				+ totalMemorySize + " GB, the max number of tasks is set as " + maxTaskCount);
	}

	/**
	 * Whether this engine writes per-raw-file sample folders nested under
	 * {@code search_results/} (see {@link MagResultLayoutUtil}) instead of
	 * directly under the workflow result root.
	 *
	 * Defaults to {@code false} (legacy flat layout, one folder per raw file as
	 * a sibling of {@code combined/}, {@code taxonomy_analysis/}, etc.), kept
	 * only as the safe default for this base class. Both concrete subclasses
	 * that use it - {@code DDAMagPfindTask} and {@code DDAMagSageTask} - now
	 * override this to {@code true}. (FragPipe's GUI task, {@code
	 * DDAMagFragPipeTask}, also overrides it to {@code true} for consistency,
	 * though its real per-sample work happens in {@code
	 * MagFragPipeHeadlessRunner}, which resolves the same layout
	 * independently.)
	 */
	protected boolean useNestedSearchResultsLayout() {
		return false;
	}

	protected void createSeparateFolders() {
		MetaData metaData = metaPar.getMetadata();
		if (metaData == null) {
			logger.error("No metadata");
			throw new IllegalArgumentException("Metadata are required");
		}
		this.rawFiles = metaPar.getMetadata().getRawFiles();
		if (rawFiles == null || rawFiles.length == 0) {
			logger.error("No raw files found in metadata");
			throw new IllegalArgumentException("Raw files are required");
		}

		this.rawDir = (new File(rawFiles[0])).getParentFile();
		if (!rawDir.exists()) {
			logger.error("No raw files found in metadata");
			throw new IllegalArgumentException("Raw files are required");
		}

		logger.info("Initializing result directories in: {}", metaPar.getDbSearchResultFile().getAbsolutePath());

		this.quanExpNames = metaPar.getMetadata().getExpNames();
		this.combinedFile = new File(metaPar.getDbSearchResultFile(), "combined");
		if (!this.combinedFile.exists()) {
			this.combinedFile.mkdir();
		}

		File sampleFoldersRoot;
		if (useNestedSearchResultsLayout()) {
			sampleFoldersRoot = MagResultLayoutUtil.resolveSearchResultsRoot(metaPar.getDbSearchResultFile());
			logger.info("Using nested search_results layout, sample folders root: {}",
					sampleFoldersRoot.getAbsolutePath());
		} else {
			sampleFoldersRoot = metaPar.getDbSearchResultFile();
		}

		this.rawFileNames = new String[rawFiles.length];
		this.separateResultFolders = new File[rawFiles.length];
		this.expNameMap = new HashMap<String, String>();

		for (int i = 0; i < rawFiles.length; i++) {
			rawFileNames[i] = rawFiles[i].substring(rawFiles[i].lastIndexOf("\\") + 1, rawFiles[i].lastIndexOf("."));
			this.separateResultFolders[i] = new File(sampleFoldersRoot, rawFileNames[i]);
			if (!this.separateResultFolders[i].exists()) {
				this.separateResultFolders[i].mkdir();
				logger.info("  Created result directory: {}", this.separateResultFolders[i].getAbsolutePath());
			} else {
				logger.info("  Result directory already exists: {}", this.separateResultFolders[i].getAbsolutePath());
			}
			expNameMap.put(rawFileNames[i], quanExpNames[i]);

			File markerFile = new File(separateResultFolders[i], "marker");
			if (!markerFile.exists()) {
				markerFile.mkdirs();
				logger.info("  Created marker result directory: {}", markerFile.getAbsolutePath());
			} else {
				logger.info("  Marker result directory already exists: {}", markerFile.getAbsolutePath());
			}
		}

		this.spectraDirFile = metaPar.getSpectraFile();
	}

	/**
	 * Advanced genome refinement with fine-tuned parameters. Provides additional
	 * control over genome selection process.
	 * 
	 * USED BY: NONE - This advanced version is not used by any child class
	 * 
	 * @param genomePepMap  Map of genome names to their associated peptides
	 * @param totalPepCount Total number of unique peptides
	 * @param threshold     Coverage threshold for genome selection
	 * @param stepThreshold Step threshold for additional filtering
	 * @param maxium        Maximum number of genomes to select
	 * @param decoyFilter   Whether to apply decoy filtering
	 * @param out           Output file for genome statistics
	 * @return Set of selected genome names
	 * @deprecated
	 */
	protected HashSet<String> refineGenomes(HashMap<String, HashSet<String>> genomePepMap, int totalPepCount,
			double threshold, double stepThreshold, int maxium, boolean decoyFilter, File out) {

		String[] genomes = genomePepMap.keySet().toArray(new String[genomePepMap.size()]);
		Arrays.sort(genomes, (g1, g2) -> {
			int s1 = genomePepMap.get(g1).size();
			int s2 = genomePepMap.get(g2).size();
			return s2 - s1;
		});

		for (int iteratorCount = 0; iteratorCount < 5; iteratorCount++) {
			HashSet<String> totalSet = new HashSet<String>();
			HashMap<String, Integer> orderMap = new HashMap<String, Integer>();
			for (int i = 0; i < genomes.length; i++) {
				int currentSize = totalSet.size();
				totalSet.addAll(genomePepMap.get(genomes[i]));
				int increaseSize = totalSet.size() - currentSize;
				orderMap.put(genomes[i], increaseSize);
			}

			Arrays.sort(genomes, (g1, g2) -> {
				return orderMap.get(g2) - orderMap.get(g1);
			});
		}

		double[] percentage = new double[genomes.length];
		HashSet<String> genomeSet = new HashSet<String>();
		HashSet<String> totalSet = new HashSet<String>();
		int breakPoint1 = -1;
		int breakPoint2 = -1;

		try (PrintWriter writer = new PrintWriter(out)) {
			writer.println(
					"Rank\tGenome\tPeptide_count\tTotal_covered_peptide_count\tPercentage\tIncrement\tIncrement ratio");
			for (int i = 0; i < genomes.length; i++) {

				int currentSize = totalSet.size();
				totalSet.addAll(genomePepMap.get(genomes[i]));
				int add = totalSet.size() - currentSize;
				double additional = ((double) add / (double) totalPepCount);

				percentage[i] = ((double) totalSet.size() / (double) totalPepCount);

				writer.println(i + "\t" + genomes[i] + "\t" + genomePepMap.get(genomes[i]).size() + "\t"
						+ totalSet.size() + "\t" + percentage[i] + "\t" + add + "\t" + additional);

				if (add > 1) {
					breakPoint1 = i;
				}
				if (additional > stepThreshold) {
					breakPoint2 = i;
				}
			}
			writer.close();
		} catch (IOException e) {
			logger.error(getTaskName() + ": error in writing genomes to " + out, e);
			System.out
					.println(format.format(new Date()) + "\t" + getTaskName() + ": error in writing genomes to " + out);
		}

		int decoyCount = 0;
		for (int i = 0; i < genomes.length; i++) {
			genomeSet.add(genomes[i]);
			if (genomes[i].startsWith(FastaManager.shuffle)) {
				decoyCount++;
			}
			if (decoyFilter) {
				if (decoyCount >= 3) {
					if (percentage[i] >= threshold) {
						break;
					}
					if (i == breakPoint1 || i == breakPoint2) {
						break;
					}
					if (genomeSet.size() > maxium) {
						break;
					}
				}
			} else {
				if (percentage[i] >= threshold) {
					break;
				}
				if (threshold < 1.0) {
					if (i == breakPoint1 || i == breakPoint2) {
						break;
					}
					if (genomeSet.size() > maxium) {
						break;
					}
				}
			}
		}

		return genomeSet;
	}

	/**
	 * Refines genomes based on peptide intensity and FDR calculation. Uses default
	 * decoy symbol (PFindResultReader.REV). Reverser proteins are also included to
	 * calculate the FDR of genomes.
	 * 
	 * USED BY: - DDAMagPfindTask: Line 1330 (genome scoring for protein inference)
	 * - DDAMagSageTask: Line 1237 (genome scoring for protein inference)
	 * 
	 * @param quanPeps Array of quantified peptides
	 * @return Map of genome names to their FDR scores
	 * @deprecated Active DDA MAG workflows should pass PSMs and use
	 *             {@link #refineGenomes(MetaPeptide[], String, MetaPSM[])}. This
	 *             overload is retained only for legacy callers.
	 */
	@Deprecated
	protected HashMap<String, Double> refineGenomes(MetaPeptide[] quanPeps) {
		return refineGenomes(quanPeps, PFindResultReader.REV);
	}

	/**
	 * Refines genomes based on peptide intensity and FDR calculation with custom
	 * decoy symbol. Reverser proteins are also included to calculate the FDR of
	 * genomes.
	 * 
	 * USED BY: - DDAMagFragPipeTask: Line 1182 (genome scoring with
	 * FastaManager.shuffle decoy symbol)
	 * 
	 * @param quanPeps    Array of quantified peptides
	 * @param decoySymble Custom decoy symbol for FDR calculation
	 * @return Map of genome names to their FDR scores
	 * @deprecated Active DDA MAG workflows should pass PSMs and use
	 *             {@link #refineGenomes(MetaPeptide[], String, MetaPSM[])}. This
	 *             overload is retained only for legacy callers.
	 */
	@Deprecated
	protected HashMap<String, Double> refineGenomes(MetaPeptide[] quanPeps, String decoySymble) {
		return refineGenomes(quanPeps, decoySymble, null);
	}

	/**
	 * Refines genomes using DIA-style scoring when PSMs are available; otherwise
	 * uses the legacy DDA MAG intensity/decoy scoring path.
	 */
	protected HashMap<String, Double> refineGenomes(MetaPeptide[] quanPeps, String decoySymbol, MetaPSM[] psms) {
		String magIdentifier = this.magDb == null ? "" : this.magDb.getIdentifier();
		DdaMagGenomeScoringService.Result result = MagPostAnalysisScoringService.scoreDdaGenomes(quanPeps, psms,
				decoySymbol, magIdentifier);
		lastDiaStyleGenomeStats = result.diaStyleGenomeStats;
		lastDiaStyleRuninfo = result.diaStyleRuninfo;

		if (result.diaStyle) {
			logger.info(getTaskName() + ": DIA-style genome scoring: " + result.genomeScoreMap.size()
					+ " genomes (1 - q_value)" + (lastDiaStyleRuninfo != null && lastDiaStyleRuninfo.useIntensityOnly
							? ", intensity-only" : ""));
		} else {
			logger.info(getTaskName() + ": " + result.targetPeptideCount + " target peptides and "
					+ result.decoyPeptideCount + " decoy peptides were used to estimate the genome level FDR");
			logger.info(getTaskName() + ": genome before filtering: " + result.genomeCountBeforeFiltering);
			logger.info(getTaskName() + ": razor genome before filtering: "
					+ result.razorGenomeCountBeforeFiltering);
			logger.info(getTaskName() + ": decoy genome count: " + result.decoyGenomeCount);
			if (result.decoyGenomeCount == 0) {
				logger.info(getTaskName()
						+ ": no decoy genome were found, the p-Value for the genome identifications will be set as 0");
			} else {
				logger.info(getTaskName() + ": refined decoy genome:"
						+ result.refinedDecoyGenomeLogIntensities.length);
				logger.info(getTaskName() + ": decoy genome Log(Intensity):"
						+ Arrays.toString(result.refinedDecoyGenomeLogIntensities));
				if (result.refinedDecoyGenomeLogIntensities.length == 0) {
					logger.info(getTaskName() + ": the number of decoy genome was not sufficient for the t-test, "
							+ "the p-Value for the genome identifications will be set as 0");
				} else if (result.refinedDecoyGenomeLogIntensities.length < 6) {
					logger.info(getTaskName() + ": decoy genome used for scoring: 6");
				} else {
					logger.info(getTaskName() + ": decoy genome used for scoring: "
							+ result.refinedDecoyGenomeLogIntensities.length);
				}
			}
			logger.info(getTaskName() + ": genome after filtering: " + result.genomeCountAfterFiltering);
		}
		return result.genomeScoreMap;
	}
	/**
	 * Refines proteins using razor protein selection with default decoy symbol.
	 * Selects the best protein for each peptide based on genome scores and
	 * intensity.
	 * 
	 * USED BY: - DDAMagPfindTask: Line 1332 (razor protein selection for protein
	 * inference) - DDAMagSageTask: Line 1238 (razor protein selection for protein
	 * inference)
	 * 
	 * @param quanPeps       Array of quantified peptides
	 * @param genomeScoreMap Map of genome names to their FDR scores
	 * @return Array of selected razor proteins (one per peptide)
	 */
	protected String[] refineProteins(MetaPeptide[] quanPeps, HashMap<String, Double> genomeScoreMap) {
		return refineProteins(quanPeps, genomeScoreMap, PFindResultReader.REV);
	}

	/**
	 * Refines proteins using razor protein selection with custom decoy symbol.
	 * Selects the best protein for each peptide based on genome scores and
	 * intensity.
	 * 
	 * USED BY: - DDAMagFragPipeTask: Line 1184 (razor protein selection with
	 * FastaManager.shuffle decoy symbol)
	 * 
	 * @param quanPeps       Array of quantified peptides
	 * @param genomeScoreMap Map of genome names to their FDR scores
	 * @param decoySymble    Custom decoy symbol for protein filtering
	 * @return Array of selected razor proteins (one per peptide)
	 */
	protected String[] refineProteins(MetaPeptide[] quanPeps, HashMap<String, Double> genomeScoreMap,
			String decoySymble) {

		logger.info(getTaskName() + ": refineProteins called with decoySymble=[{}], quanPeps.length={}",
				decoySymble, quanPeps != null ? quanPeps.length : 0);

		String magIdentifier = this.magDb == null ? "" : this.magDb.getIdentifier();
		DdaMagProteinInferenceService.Result result = MagPostAnalysisScoringService.inferDdaProteins(quanPeps,
				genomeScoreMap, null, magIdentifier, decoySymble);
		logger.info(getTaskName() + ": refineProteins: peptides with only decoy proteins (excluded from razor)={}",
				result.peptidesWithOnlyDecoyProteins);
		return result.razorProteins;
	}

	protected boolean exportGenome(HashMap<String, double[]> genomeIntensityMap, HashMap<String, Integer> psmCountMap,
			HashMap<String, HashSet<String>> pepCountMap, HashMap<String, HashSet<String>> uniquePepCountMap,
			HashMap<String, ArrayList<Double>> genomeQvalueMap, File taxFile, String[] expNames,
			MetaPythonReportTask task) {

		MagFuncSearcher searcher = null;
		try {
			searcher = new MagFuncSearcher(metaPar.getMicroDb(), mrp.getToolPath(ToolType.FUNCDB), magDb);

		} catch (NumberFormatException | SQLException | IOException e) {
			// TODO Auto-generated catch block
			logger.error(getTaskName() + ": error in genome analysis, task failed", e);
			System.err.println(
					format.format(new Date()) + "\t" + getTaskName() + ": error in genome analysis, task failed");

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

			// Add to manifest
			File reportsDir = new File(metaPar.getReportDir(), "reports");
			ReportManifestGenerator manifestGen = new ReportManifestGenerator(reportsDir);
			manifestGen.addGenomeReport(reportGenomeHtml);

		} catch (IOException e) {
			// TODO Auto-generated catch block
			logger.error("Error in writing genome information to " + genomeHtmlFile, e);
			System.err.println(
					format.format(new Date()) + "\t" + "Error in writing genome information to " + genomeHtmlFile);
			return false;
		}

		File genomeFile = new File(taxFile, "Genome.tsv");
		PrintWriter genomeWriter = null;
		try {
			genomeWriter = new PrintWriter(genomeFile);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			logger.error("Error in writing genome information to " + genomeFile, e);
			System.err
					.println(format.format(new Date()) + "\t" + "Error in writing genome information to " + genomeFile);
			return false;
		}

		StringBuilder genomeTitlesb = new StringBuilder();
		genomeTitlesb.append("Genome\t");
		genomeTitlesb.append("Score\t");
		genomeTitlesb.append("PSM count\t");
		genomeTitlesb.append("Peptide count\t");
		genomeTitlesb.append("Unique peptide count\t");
		genomeTitlesb.append("GTDB_Name\t");
		genomeTitlesb.append("GTDB_Rank\t");
		TaxonomyRanks[] mainRanks = TaxonomyRanks.getMainRanks();
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

		double[] cellularIntensity = new double[expNames.length];
		HashMap<String, String[]>[] taxRankMaps = new HashMap[mainRanks.length];
		HashMap<String, double[]>[] taxIntensityMaps = new HashMap[mainRanks.length];
		HashMap<String, Double>[] taxScoreMaps = new HashMap[mainRanks.length];
		for (int i = 0; i < taxRankMaps.length; i++) {
			taxRankMaps[i] = new HashMap<String, String[]>();
			taxIntensityMaps[i] = new HashMap<String, double[]>();
			taxScoreMaps[i] = new HashMap<String, Double>();
		}

		Arrays.sort(genomes, new Comparator<String>() {
			@Override
			public int compare(String o1, String o2) {
				// TODO Auto-generated method stub
				return genomeQvalueMap.get(o2).size() - genomeQvalueMap.get(o1).size();
			}
		});

		double[] logGenomePepCount = new double[genomes.length];
		for (int i = 0; i < genomes.length; i++) {
			logGenomePepCount[i] = Math.log(genomeQvalueMap.get(genomes[i]).size());
		}
		double[] genomeScores = new double[genomes.length];
		double mean = MathTool.calculateArithmeticMean(logGenomePepCount);
		if (mean > 0) {
			TTest tTest = new TTest();
			double add = -Math.log10(tTest.tTest(logGenomePepCount[0], logGenomePepCount));
			genomeScores[0] = add * 2;
			for (int i = 1; i < logGenomePepCount.length; i++) {
				if (logGenomePepCount[i] > 0) {
					double pValue = tTest.tTest(logGenomePepCount[i], logGenomePepCount);
					genomeScores[i] = -Math.log10(pValue);
					if (logGenomePepCount[i] > mean) {
						genomeScores[i] = add + genomeScores[i];
					} else {
						genomeScores[i] = add - genomeScores[i];
					}

					if (genomeScores[i] < 0) {
						genomeScores[i] = 0;
					}
				}
			}
		}

		for (int i = 0; i < genomes.length; i++) {
			String[] taxaContent = genomeTaxaMap.get(genomes[i]);
			if (genomeIntensityMap.containsKey(genomes[i])) {
				double[] genomeIntensity = genomeIntensityMap.get(genomes[i]);
				for (int j = 0; j < genomeIntensity.length; j++) {
					cellularIntensity[j] += genomeIntensity[j];
				}

				StringBuilder genomeTaxaSb = new StringBuilder();
				genomeTaxaSb.append(genomes[i]).append("\t");

				genomeTaxaSb.append(FormatTool.getDF2().format(genomeScores[i])).append("\t");

				if (psmCountMap.containsKey(genomes[i])) {
					genomeTaxaSb.append(psmCountMap.get(genomes[i])).append("\t");
				} else {
					genomeTaxaSb.append("0").append("\t");
				}
				if (pepCountMap.containsKey(genomes[i])) {
					genomeTaxaSb.append(pepCountMap.get(genomes[i]).size()).append("\t");
				} else {
					genomeTaxaSb.append("0").append("\t");
				}
				if (uniquePepCountMap.containsKey(genomes[i])) {
					genomeTaxaSb.append(uniquePepCountMap.get(genomes[i]).size());
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
						taxScoreMaps[j].put(taxon, taxScoreMaps[j].get(taxon) + genomeScores[i]);
					} else {
						taxScoreMaps[j].put(taxon, genomeScores[i]);
					}
				}
			}
		}
		genomeWriter.close();

		File taxaFile = new File(taxFile, "Taxa.tsv");
		File taxaReportFile = new File(taxFile, "Taxa_report.tsv");
		PrintWriter taxaWriter = null;
		PrintWriter taxaReportWriter = null;
		try {
			taxaWriter = new PrintWriter(taxaFile);
			taxaReportWriter = new PrintWriter(taxaReportFile);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			logger.error("Error in writing taxa information to " + taxaFile, e);
			System.err.println(format.format(new Date()) + "\t" + "Error in writing taxa information to " + taxaFile);
			return false;
		}

		StringBuilder taxaTitlesb = new StringBuilder();
		taxaTitlesb.append("Name\t");
		taxaTitlesb.append("Rank");

		for (int i = 0; i < mainRanks.length; i++) {
			taxaTitlesb.append("\t").append(mainRanks[i].getName());
		}

		taxaTitlesb.append("\tScore");

		for (String exp : expNames) {
			if (exp.startsWith("Intensity ")) {
				taxaTitlesb.append("\t").append(exp);
			} else {
				taxaTitlesb.append("\tIntensity ").append(exp);
			}
		}
		taxaWriter.println(taxaTitlesb);
		taxaReportWriter.println(taxaTitlesb);

		double maxTaxonScore = 0;
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

				double taxonScore = taxScoreMaps[i].get(taxon);
				taxaSb.append("\t").append(FormatTool.getDF2().format(taxonScore));

				for (int j = 0; j < intensity.length; j++) {
					taxaSb.append("\t").append(intensity[j]);
				}

				taxaWriter.println(taxaSb);

				if (taxonScore > maxTaxonScore) {
					maxTaxonScore = taxonScore;
				}
			}
		}

		taxaWriter.close();

		for (int i = 0; i < taxRankMaps.length; i++) {
			for (String taxon : taxRankMaps[i].keySet()) {
				double taxonScore = taxScoreMaps[i].get(taxon);
				if (taxonScore < maxTaxonScore * 0.005) {
					continue;
				}
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

				taxaSb.append("\t").append(FormatTool.getDF2().format(taxonScore));

				for (int j = 0; j < intensity.length; j++) {
					taxaSb.append("\t").append(intensity[j]);
				}

				taxaReportWriter.println(taxaSb);
			}
		}

		taxaReportWriter.close();

		File jsOutput = new File(metaPar.getReportDir() + "\\reports\\tree.js");
		if (!jsOutput.getParentFile().exists()) {
			jsOutput.getParentFile().mkdirs();
		}

		if (!jsOutput.exists() || jsOutput.length() == 0) {
			MetaTreeHandler.exportJS(taxRankMaps, taxIntensityMaps, cellularIntensity, expNames, jsOutput);
		}

		// Generate bubble chart report from tree.js data (optional - don't fail if this
		// doesn't work)
		File bubbleChartOutput = new File(metaPar.getReportDir() + "\\reports\\bubble_chart_report.html");
		if (!bubbleChartOutput.exists() || bubbleChartOutput.length() == 0) {
			try {
				generateBubbleChartReport(jsOutput, bubbleChartOutput);
			} catch (Exception e) {
				logger.warn(getTaskName() + ": Bubble chart generation failed, continuing with other reports: "
						+ e.getMessage());
				System.out.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": Bubble chart generation failed, continuing with other reports");
			}
		}

		boolean finished = false;

		File report_taxonomy_summary = new File(metaPar.getReportDir(), "report_taxonomy_summary.html");
		if (!report_taxonomy_summary.exists() || report_taxonomy_summary.length() == 0) {

			if (summaryMetaFile != null && summaryMetaFile.exists()) {
				finished = task.addTask(MetaPythonReportTask.TAXONOMY, taxaReportFile, report_taxonomy_summary,
						summaryMetaFile);
			} else {
				finished = task.addTask(MetaPythonReportTask.TAXONOMY, taxaReportFile, report_taxonomy_summary);
			}
		} else {
			finished = true;
		}

		logger.info(getTaskName() + ": taxonomy analysis started");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": taxonomy analysis started");

		updateProgress(96, 100);

		return finished;
	}

	/**
	 * @deprecated MAG functional annotation is now handled by
	 *             {@link MagFunctionalAnnotationService}.
	 */
	@Deprecated
	protected MetaProteinAnnoMag[] funcAnnoSqlite(MetaProtein[] pros, String[] fileNames) {
		MagFunctionalAnnotationService.Request request = new MagFunctionalAnnotationService.Request();
		request.resultDir = metaPar.getDbSearchResultFile();
		request.magDb = magDb;
		request.microDbCatalogueId = metaPar.getMicroDb();
		request.functionDefinitionDb = new File(mrp.getToolPath(ToolType.FUNCDB));
		request.analysisType = MetaConstants.pFindMag;
		request.quantitativeMode = metaPar.getQuantitativeMode();
		request.intensityTitles = fileNames;
		request.metaProteins = pros;

		try {
			MagFunctionalAnnotationService.Result result = MagFunctionalAnnotationService.run(request);
			if (result.annotatedProteins.length > 0) {
				logger.info(getTaskName() + ": proteins with functional annotations have been exported to "
						+ result.annotatedProteinsTsv.getName());
				System.out.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": proteins with functional annotations have been exported to "
						+ result.annotatedProteinsTsv.getName());
			}
			return result.annotatedProteins;
		} catch (IOException | SQLException e) {
			logger.error(getTaskName() + ": error in proteins functional annotations, task failed");
			System.err.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": error in proteins functional annotations, task failed");
			return null;
		}
	}

	/**
	 * Generate identification summary report
	 */
	private boolean generateIdentificationSummary(HashMap<String, Integer> fileIdMap, int[] filePsmCount,
			HashSet<String>[] fileSequenceSets, String taskName) {
		try {
			DdaMagReportExportService.writeIdentificationSummary(combinedFile, metaPar.getDbSearchResultFile(),
					spectraDirFile, fileIdMap, filePsmCount, fileSequenceSets, expNameMap);
			return true;
		} catch (IOException e) {
			logger.error(taskName + ": error in writing identification information", e);
			System.err.println(
					format.format(new Date()) + "\t" + taskName + ": error in writing identification information");
			return false;
		}
	}

	/**
	 * Generate peptide report
	 */
	private boolean generatePeptideReport(ArrayList<MetaPeptide> peptideList, HashMap<String, MetaPSM> psmMap,
			HashMap<String, int[]> pepPsmCountMap, HashMap<String, int[]> pepChargeMap,
			HashMap<String, Integer> proteinIdMap, ArrayList<String> fileList, String quanMode, String[] intensityTitle,
			String taskName) {
		return generatePeptideReport(peptideList, psmMap, pepPsmCountMap, pepChargeMap, proteinIdMap, fileList,
				quanMode, intensityTitle, taskName, null);
	}

	private boolean generatePeptideReport(ArrayList<MetaPeptide> peptideList, HashMap<String, MetaPSM> psmMap,
			HashMap<String, int[]> pepPsmCountMap, HashMap<String, int[]> pepChargeMap,
			HashMap<String, Integer> proteinIdMap, ArrayList<String> fileList, String quanMode, String[] intensityTitle,
			String taskName, HashMap<String, String> razorProteinMap) {

		logger.info(taskName + ": exporting peptide report started");
		System.out.println(format.format(new Date()) + "\t" + taskName + ": exporting peptide report started");
		try {
			final_pep_txt = new File(this.metaPar.getDbSearchResultFile(), "final_peptides.tsv");
			DdaMagReportExportService.writePeptideReport(final_pep_txt, peptideList, psmMap, pepPsmCountMap,
					pepChargeMap, proteinIdMap, fileList, quanMode, intensityTitle, razorProteinMap);
		} catch (IOException e) {
			logger.error(taskName + ": error in writing the peptides to " + final_pep_txt, e);
			System.err.println(format.format(new Date()) + "\t" + taskName + ": error in writing the peptides to "
					+ final_pep_txt);
			return false;
		}
		logger.info(taskName + ": peptide report has been exported to " + final_pep_txt);
		System.out.println(
				format.format(new Date()) + "\t" + taskName + ": peptide report has been exported to " + final_pep_txt);
		return final_pep_txt.exists() && final_pep_txt.length() > 0;
	}

	private HashMap<String, String> buildRazorProteinMap(MetaPeptide[] peptides, String[] razorProteins) {
		HashMap<String, String> map = new HashMap<String, String>();
		if (peptides == null || razorProteins == null) {
			return map;
		}
		for (int i = 0; i < peptides.length && i < razorProteins.length; i++) {
			if (peptides[i] == null || razorProteins[i] == null) {
				continue;
			}
			if (peptides[i].getModSeq() != null) {
				map.put(peptides[i].getModSeq(), razorProteins[i]);
			}
			if (peptides[i].getSequence() != null) {
				map.put(peptides[i].getSequence(), razorProteins[i]);
			}
		}
		return map;
	}

	/**
	 * Generate protein report (backward-compatible: no peptide list, Peptide IDs/razor columns empty).
	 */
	protected boolean generateProteinReport(String[] proteins, HashMap<String, double[]> proteinIntensityMap,
			HashMap<String, Double> proteinTotalIntenMap, HashMap<String, HashSet<String>> totalPepCountMap,
			HashMap<String, HashSet<String>> razorPepCountMap, HashMap<String, Integer> psmCountMap,
			HashMap<String, Integer> proteinIdMap, String quanMode, String[] intensityTitle, String taskName) {
		return generateProteinReport(proteins, proteinIntensityMap, proteinTotalIntenMap, totalPepCountMap,
				razorPepCountMap, psmCountMap, proteinIdMap, quanMode, intensityTitle, taskName, null, null);
	}

	/**
	 * Generate protein report.
	 */
	private boolean generateProteinReport(String[] proteins, HashMap<String, double[]> proteinIntensityMap,
			HashMap<String, Double> proteinTotalIntenMap, HashMap<String, HashSet<String>> totalPepCountMap,
			HashMap<String, HashSet<String>> razorPepCountMap, HashMap<String, Integer> psmCountMap,
			HashMap<String, Integer> proteinIdMap, String quanMode, String[] intensityTitle, String taskName,
			ArrayList<MetaPeptide> peptideList) {
		return generateProteinReport(proteins, proteinIntensityMap, proteinTotalIntenMap, totalPepCountMap,
				razorPepCountMap, psmCountMap, proteinIdMap, quanMode, intensityTitle, taskName, peptideList, null);
	}

	private boolean generateProteinReport(String[] proteins, HashMap<String, double[]> proteinIntensityMap,
			HashMap<String, Double> proteinTotalIntenMap, HashMap<String, HashSet<String>> totalPepCountMap,
			HashMap<String, HashSet<String>> razorPepCountMap, HashMap<String, Integer> psmCountMap,
			HashMap<String, Integer> proteinIdMap, String quanMode, String[] intensityTitle, String taskName,
			ArrayList<MetaPeptide> peptideList, HashMap<String, double[]> proteinScoreMap) {

		logger.info(taskName + ": exporting protein report started");
		System.out.println(format.format(new Date()) + "\t" + taskName + ": exporting protein report started");
		try {
			final_pro_txt = new File(this.metaPar.getDbSearchResultFile(), "razor_proteins.tsv");
			DdaMagReportExportService.writeProteinReport(final_pro_txt, proteins, proteinIntensityMap,
					proteinTotalIntenMap, totalPepCountMap, razorPepCountMap, psmCountMap, quanMode, intensityTitle,
					peptideList, proteinScoreMap);
		} catch (IOException e) {
			logger.error(taskName + ": error in writing final protein result to " + final_pro_txt.getName(), e);
			System.err.println(format.format(new Date()) + "\t" + taskName
					+ ": error in writing final protein result to " + final_pro_txt.getName());
			return false;
		}
		logger.info(taskName + ": protein report has been exported to " + final_pro_txt.getName());
		System.out.println(format.format(new Date()) + "\t" + taskName + ": protein report has been exported to "
				+ final_pro_txt.getName());
		return true;
	}

	protected boolean generateAllProteinReport(String[] proteins, HashMap<String, double[]> proteinIntensityMap,
			HashMap<String, Double> proteinTotalIntenMap, HashMap<String, HashSet<String>> totalPepCountMap,
			HashMap<String, HashSet<String>> razorPepCountMap, HashMap<String, Integer> psmCountMap, String quanMode,
			String[] intensityTitle, String taskName, ArrayList<MetaPeptide> peptideList) {
		return generateAllProteinReport(proteins, proteinIntensityMap, proteinTotalIntenMap, totalPepCountMap,
				razorPepCountMap, psmCountMap, quanMode, intensityTitle, taskName, peptideList, null);
	}

	protected boolean generateAllProteinReport(String[] proteins, HashMap<String, double[]> proteinIntensityMap,
			HashMap<String, Double> proteinTotalIntenMap, HashMap<String, HashSet<String>> totalPepCountMap,
			HashMap<String, HashSet<String>> razorPepCountMap, HashMap<String, Integer> psmCountMap, String quanMode,
			String[] intensityTitle, String taskName, ArrayList<MetaPeptide> peptideList,
			HashMap<String, double[]> proteinScoreMap) {

		File allProteinFile = new File(this.metaPar.getDbSearchResultFile(), "all_proteins.tsv");
		logger.info(taskName + ": exporting all-protein report started");
		System.out.println(format.format(new Date()) + "\t" + taskName + ": exporting all-protein report started");
		try {
			DdaMagReportExportService.writeAllProteinReport(allProteinFile, proteins, proteinIntensityMap,
					proteinTotalIntenMap, totalPepCountMap, razorPepCountMap, psmCountMap, quanMode, intensityTitle,
					peptideList, proteinScoreMap);
		} catch (IOException e) {
			logger.error(taskName + ": error in writing all-protein result to " + allProteinFile.getName(), e);
			System.err.println(format.format(new Date()) + "\t" + taskName
					+ ": error in writing all-protein result to " + allProteinFile.getName());
			return false;
		}
		logger.info(taskName + ": all-protein report has been exported to " + allProteinFile.getName());
		System.out.println(format.format(new Date()) + "\t" + taskName
				+ ": all-protein report has been exported to " + allProteinFile.getName());
		return true;
	}

	/**
	 * Perform functional annotation.
	 */
	private boolean performFunctionalAnnotation(String[] proteins, HashMap<String, double[]> proteinIntensityMap,
			String[] intensityTitle, String taskName,
			HashMap<String, HashSet<String>> razorPepCountMap, HashMap<String, Integer> razorPsmCountMap) {

		logger.info(taskName + ": functional annotation started");
		System.out.println(format.format(new Date()) + "\t" + taskName + ": functional annotation started");
		MagFunctionalAnnotationService.Request request = new MagFunctionalAnnotationService.Request();
		request.resultDir = metaPar.getDbSearchResultFile();
		request.magDb = magDb;
		request.microDbCatalogueId = metaPar.getMicroDb();
		request.functionDefinitionDb = new File(mrp.getToolPath(ToolType.FUNCDB));
		request.analysisType = MetaConstants.pFindMag;
		request.quantitativeMode = metaPar.getQuantitativeMode();
		request.intensityTitles = intensityTitle;
		request.proteins = proteins;
		request.proteinIntensityMap = proteinIntensityMap;
		request.razorPepCountMap = razorPepCountMap;
		request.razorPsmCountMap = razorPsmCountMap;

		try {
			MagFunctionalAnnotationService.run(request);
		} catch (Exception e) {
			logger.error(taskName + ": error in proteins functional annotations, task failed", e);
			System.err.println(format.format(new Date()) + "\t" + taskName
					+ ": error in proteins functional annotations, task failed");
			return false;
		}

		logger.info(taskName + ": functional annotation finished");
		System.out.println(format.format(new Date()) + "\t" + taskName + ": functional annotation finished");
		return true;
	}

	/**
	 * Perform taxonomy analysis (DDA table schema; DIA workflows should call this with the same maps).
	 */
	protected boolean performTaxonomyAnalysis(HashMap<String, double[]> genomeIntensityMap,
			HashMap<String, Double> genomeScoreMap, HashMap<String, Integer> psmCountMap,
			HashMap<String, HashSet<String>> totalPepCountMap, HashMap<String, HashSet<String>> razorPepCountMap,
			HashMap<String, HashSet<String>> uniquePepCountMap,
			HashMap<String, HashSet<String>> totalPepWithILMap, HashMap<String, HashSet<String>> razorPepWithILMap,
			HashMap<String, HashSet<String>> uniquePepWithILMap,
			String[] intensityTitle, String taskName) {

		logger.info(taskName + ": taxonomy analysis started");
		System.out.println(format.format(new Date()) + "\t" + taskName + ": taxonomy analysis started");

		File taxFile = new File(metaPar.getDbSearchResultFile(), "taxonomy_analysis");
		MagTaxonomyAnalysisService.Request request = new MagTaxonomyAnalysisService.Request();
		request.taxonomyAnalysisDir = taxFile;
		request.magDb = magDb;
		request.generateFinalGenomeRanking = true;
		request.genomeIntensityMap = genomeIntensityMap;
		request.genomeScoreMap = genomeScoreMap;
		request.genomeStats = lastDiaStyleGenomeStats;
		request.genomeScoreRuninfo = lastDiaStyleRuninfo;
		request.psmCountMap = psmCountMap;
		request.totalPepCountMap = totalPepCountMap;
		request.razorPepCountMap = razorPepCountMap;
		request.uniquePepCountMap = uniquePepCountMap;
		request.totalPepWithILMap = totalPepWithILMap;
		request.razorPepWithILMap = razorPepWithILMap;
		request.uniquePepWithILMap = uniquePepWithILMap;
		request.intensityTitles = intensityTitle;

		MagTaxonomyAnalysisService.Result result;
		try {
			result = MagTaxonomyAnalysisService.run(request);
		} catch (IOException | RuntimeException e) {
			logger.error(taskName + ": error in taxonomy analysis, task failed", e);
			System.err.println(format.format(new Date()) + "\t" + taskName + ": error in taxonomy analysis, task failed");
			return false;
		}

		try {
			File reportGenomeHtml = new File(metaPar.getReportDir(), "reports/Genome.html");
			if (reportGenomeHtml.getParentFile() != null && !reportGenomeHtml.getParentFile().exists()) {
				reportGenomeHtml.getParentFile().mkdirs();
			}
			Files.copy(result.genomeHtml.toPath(), reportGenomeHtml.toPath(), StandardCopyOption.REPLACE_EXISTING);
			logger.info("Copied Genome.html to reports directory: " + reportGenomeHtml.getAbsolutePath());

			File reportsDir = new File(metaPar.getReportDir(), "reports");
			ReportManifestGenerator manifestGen = new ReportManifestGenerator(reportsDir);
			manifestGen.addGenomeReport(reportGenomeHtml);
		} catch (IOException e) {
			logger.error("Error copying genome report to reports directory", e);
			System.err.println(format.format(new Date()) + "\t" + "Error copying genome report to reports directory");
			return false;
		}

		if (result.rankedSelection != null && !result.rankedSelection.orderedGenomes.isEmpty()) {
			logger.info(
					"Genome list from ranked.tsv (GreedyGain export labels): focusKnee={}, comprehensiveKnee={}, {} genomes",
					result.rankedSelection.focusKnee, result.rankedSelection.comprehensiveKnee,
					Integer.valueOf(result.rankedSelection.orderedGenomes.size()));
		}

		logger.info(taskName + ": taxonomy analysis finished");
		System.out.println(format.format(new Date()) + "\t" + taskName + ": taxonomy analysis finished");
		return true;
	}

	/**
	 * Writes genome_rank_scores_runinfo.txt in taxonomy_analysis when DIA-style scoring was used (lastDiaStyleRuninfo set).
	 *
	 * @deprecated DDA MAG taxonomy analysis now writes this from
	 *             {@link MagTaxonomyAnalysisService}.
	 */
	@Deprecated
	private void writeGenomeRankScoresRuninfo(File taxFile, File rankScoresFile) {
		if (lastDiaStyleRuninfo == null) {
			return;
		}
		File runInfoFile = new File(taxFile, "genome_rank_scores_runinfo.txt");
		try (PrintWriter pw = new PrintWriter(runInfoFile)) {
			DdaMagGenomeScoringService.DdaGenomeScoreRuninfo r = lastDiaStyleRuninfo;
			pw.println("=================================================================================");
			pw.println("genome_rank_scores run info - diagnostic report (DDA MAG, DIA-style scoring)");
			pw.println("=================================================================================");
			pw.println();
			pw.println("--- Input / Output ---");
			pw.println("workflow\tDDA MAG post-analysis");
			pw.println("output_scores\t" + rankScoresFile.getAbsolutePath());
			pw.println();
			pw.println("--- Peptide counts ---");
			pw.println("n_peptides_target\t" + r.nTargetPeptides + "\t(after PSM collapse, unique peptides only m=1)");
			pw.println("n_peptides_decoy\t" + r.nDecoyPeptides);
			pw.println();
			pw.println("--- Peptide sharing (target) ---");
			pw.println("n_unique_peptides\t" + r.nUniquePeptides + "\t(m=1 genome)");
			pw.println("n_shared_peptides\t" + r.nSharedPeptides + "\t(m>1 genomes, skipped in unique-only mode)");
			pw.print("m_distribution\t");
			List<Integer> mKeys = new ArrayList<>(r.mDistribution.keySet());
			Collections.sort(mKeys);
			for (int i = 0; i < mKeys.size(); i++) {
				if (i > 0) pw.print("; ");
				pw.print("m=" + mKeys.get(i) + ":" + r.mDistribution.get(mKeys.get(i)));
			}
			pw.println();
			pw.println();
			pw.println("--- Decoy pool ---");
			pw.println("decoy_pool_size\t" + r.decoyPoolSize);
			pw.println("decoy_e_min\t" + String.format("%.4f", r.decoyPoolSize > 0 ? r.decoyEMin : 0));
			pw.println("decoy_e_max\t" + String.format("%.4f", r.decoyPoolSize > 0 ? r.decoyEMax : 0));
			pw.println("decoy_e_mean\t" + String.format("%.4f", r.decoyEMean));
			pw.println();
			pw.println("--- Score (e) distribution (e = PSM score capped to [0, cap]; if target e is all at cap, score does not distinguish genomes) ---");
			if (r.targetERaw != null && r.decoyERaw != null) {
				double[] targetEStats = percentileStats(r.targetERaw);
				double[] decoyEStats = percentileStats(r.decoyERaw);
				if (targetEStats != null) {
					pw.println("target_e\tmin=" + String.format("%.4f", targetEStats[0]) + "\tmax=" + String.format("%.4f", targetEStats[1])
							+ "\tmean=" + String.format("%.4f", targetEStats[2]) + "\tP50=" + String.format("%.4f", targetEStats[6])
							+ "\tP99=" + String.format("%.4f", targetEStats[9]));
				}
				if (decoyEStats != null) {
					pw.println("decoy_e\tmin=" + String.format("%.4f", decoyEStats[0]) + "\tmax=" + String.format("%.4f", decoyEStats[1])
							+ "\tmean=" + String.format("%.4f", decoyEStats[2]) + "\tP50=" + String.format("%.4f", decoyEStats[6])
							+ "\tP99=" + String.format("%.4f", decoyEStats[9]));
				}
				boolean eOverlap = targetEStats != null && decoyEStats != null
						&& targetEStats[0] <= decoyEStats[1] && decoyEStats[0] <= targetEStats[1];
				pw.println("e_overlap\t" + (eOverlap ? "YES (target and decoy e ranges intersect)" : "NO"));
				double capVal = r.cap;
				long targetAtCap = (r.targetERaw != null) ? r.targetERaw.stream().filter(x -> x >= capVal - 1e-6).count() : 0;
				int targetN = (r.targetERaw != null) ? r.targetERaw.size() : 0;
				double pctAtCap = targetN > 0 ? 100.0 * targetAtCap / targetN : 0;
				pw.println("target_peptides_at_e_cap\t" + targetAtCap + " / " + targetN + " (" + String.format("%.1f", pctAtCap) + "%)");
				pw.println("e_note\tIf most target e are at cap (" + String.format("%.1f", capVal) + "), score does not discriminate; genome S is driven mainly by intensity. Consider raising e_cap or using a score scale that spreads target above decoy.");
			}
			pw.println();
			pw.println("--- Scoring parameters ---");
			pw.println("peptide_mode\tunique only (m=1, shared skipped)");
			pw.println("scoring_mode\t" + (r.useIntensityOnly ? "intensity only (e=1, PSM score not used)" : "PSM score + MetaPeptide intensity (DDA)"));
			pw.println("intensity_factor\t0.5 + 0.5*norm(log2(1+totalIntensity))");
			pw.println("log2_I_range\t" + String.format("%.2f", r.minLog) + " to " + String.format("%.2f", r.maxLog));
			pw.println("top_k\t" + r.k);
			pw.println("resamples_B\t" + r.B);
			pw.println("e_cap\t" + String.format("%.1f", r.cap));
			pw.println();
			int nGenomes = r.sortedByScore.size();
			double minSid = Double.MAX_VALUE, maxSid = -Double.MAX_VALUE, sumSid = 0;
			HashMap<Double, Integer> sidCountMap = new HashMap<>();
			for (Double s : r.SMap.values()) {
				minSid = Math.min(minSid, s);
				maxSid = Math.max(maxSid, s);
				sumSid += s;
				sidCountMap.merge(s, 1, Integer::sum);
			}
			double meanSid = nGenomes > 0 ? sumSid / nGenomes : 0;
			int nGenomesQ01 = 0;
			for (int i = 0; i < r.qvals.length; i++) {
				if (r.qvals[i] < 0.01) nGenomesQ01++;
			}
			pw.println("--- Genome scores (S) summary ---");
			pw.println("n_genomes\t" + nGenomes);
			pw.println("S_min\t" + String.format("%.4f", nGenomes > 0 ? minSid : 0));
			pw.println("S_max\t" + String.format("%.4f", nGenomes > 0 ? maxSid : 0));
			pw.println("S_mean\t" + String.format("%.4f", meanSid));
			pw.println("n_genomes_q_lt_0.01\t" + nGenomesQ01);
			pw.print("S_ties\t");
			List<Double> sidKeys = new ArrayList<>(sidCountMap.keySet());
			Collections.sort(sidKeys, Comparator.reverseOrder());
			for (int i = 0; i < sidKeys.size(); i++) {
				if (i > 0) pw.print("; ");
				pw.print(String.format("%.4f", sidKeys.get(i)) + ":" + sidCountMap.get(sidKeys.get(i)));
			}
			pw.println();
			pw.println();
			pw.println("--- Per-genome details (by rank) ---");
			pw.println("Format: Rank, Genome, genome_score_S, n_peptides_total, n_top_k, genome_p_emp, genome_q_value");
			for (int i = 0; i < r.sortedByScore.size(); i++) {
				String g = r.sortedByScore.get(i);
				int M = r.genomeTotalPeptidesMap.containsKey(g) ? r.genomeTotalPeptidesMap.get(g) : 0;
				int ng = r.genomeNMap.containsKey(g) ? r.genomeNMap.get(g) : 0;
				pw.println((i + 1) + "\t" + g + "\t" + String.format("%.4f", r.SMap.get(g)) + "\t" + M + "\t" + ng
						+ "\t" + String.format("%.6f", r.pempArr[i]) + "\t" + String.format("%.6f", r.qvals[i]));
			}
			pw.println();
			pw.println("--- Top-k peptides per genome ---");
			for (int i = 0; i < r.sortedByScore.size(); i++) {
				String g = r.sortedByScore.get(i);
				pw.println("## Rank " + (i + 1) + " " + g + " S=" + String.format("%.4f", r.SMap.get(g)));
				List<DdaMagGenomeScoringService.DdaWeightedScore> topK = r.genomeTopKMap.get(g);
				if (topK != null) {
					for (int j = 0; j < topK.size(); j++) {
						DdaMagGenomeScoringService.DdaWeightedScore wps = topK.get(j);
						double mVal = (wps.w > 0) ? 1.0 / wps.w : 0;
						String seq = (wps.seq != null) ? wps.seq : "(n/a)";
						pw.println("  " + (j + 1) + "\t" + seq + "\te=" + String.format("%.4f", wps.e)
								+ "\tf(I)=" + String.format("%.4f", wps.intensityFactor)
								+ "\tw=" + String.format("%.4f", wps.w) + "\tm=" + (int) mVal
								+ "\tw*e*f=" + String.format("%.4f", wps.we));
					}
				}
				pw.println();
			}
		} catch (FileNotFoundException e) {
			logger.error("Error writing " + runInfoFile, e);
		}
	}

	/** Returns [min, max, mean, p1, p5, p25, p50, p75, p95, p99] or null if empty. */
	private static double[] percentileStats(ArrayList<Double> list) {
		if (list == null || list.isEmpty()) {
			return null;
		}
		double[] arr = new double[list.size()];
		for (int i = 0; i < arr.length; i++) {
			arr[i] = list.get(i);
		}
		Arrays.sort(arr);
		double min = arr[0];
		double max = arr[arr.length - 1];
		double sum = 0;
		for (double v : arr) {
			sum += v;
		}
		double mean = sum / arr.length;
		double p = (double) (arr.length - 1);
		double p1 = arr[Math.min((int) Math.round(0.01 * p), arr.length - 1)];
		double p5 = arr[Math.min((int) Math.round(0.05 * p), arr.length - 1)];
		double p25 = arr[Math.min((int) Math.round(0.25 * p), arr.length - 1)];
		double p50 = arr[Math.min((int) Math.round(0.50 * p), arr.length - 1)];
		double p75 = arr[Math.min((int) Math.round(0.75 * p), arr.length - 1)];
		double p95 = arr[Math.min((int) Math.round(0.95 * p), arr.length - 1)];
		double p99 = arr[Math.min((int) Math.round(0.99 * p), arr.length - 1)];
		return new double[] { min, max, mean, p1, p5, p25, p50, p75, p95, p99 };
	}

	/**
	 * Generate final reports using MetaPythonReportTask
	 */
	protected boolean generateFinalReports(String taskName) {

		updateMessage(taskName + ": exporting report");

		// Copy report templates in background
		if (!copyReportTemplatesInBackground()) {
			logger.error(getTaskName() + ": Failed to copy report templates");
			return false;
		}

		logger.info(taskName + ": exporting report started");
		System.out.println(format.format(new Date()) + "\t" + taskName + ": exporting report started");

		// Report templates are already copied in the main call() method
		// Just ensure the report directory exists and set the paths
		File reportDir = new File(this.metaPar.getDbSearchResultFile(), "report");
		if (!reportDir.exists()) {
			reportDir.mkdirs();
		}

		File reportHtmlDir = new File(reportDir + "\\reports");
		if (!reportHtmlDir.exists()) {
			reportHtmlDir.mkdirs();
		}

		// Export metadata
		File summaryMetaFile = new File(this.metaPar.getDbSearchResultFile(), "metainfo.tsv");
		try {
			metaPar.getMetadata().exportMetadata(summaryMetaFile);
		} catch (FileNotFoundException e) {
			logger.error(taskName + ": error in exporting metadata", e);
		}

		// Generate summary report
		MetaPythonReportTask reportTask = new MetaPythonReportTask(threadCount, mrp, logger);
		reportTask.setWorkflowType("MAG");
		reportTask.initializeManifestGenerator(reportHtmlDir);

		// Create initial manifest file to ensure index page can load
		ReportManifestGenerator manifestGen = new ReportManifestGenerator(reportHtmlDir);
		manifestGen.createInitialManifest();
		// Add individual report tasks
		boolean reportSuccess = true;

		// Add summary report
		File summaryReportFile = new File(this.metaPar.getDbSearchResultFile(), "final_summary.tsv");
		if (summaryReportFile.exists() && summaryReportFile.length() > 0) {
			File report_ID_summary = new File(reportHtmlDir, "report_ID_summary.html");
			if (summaryMetaFile != null && summaryMetaFile.exists()) {
				reportSuccess &= reportTask.addTask(MetaPythonReportTask.ID_SUMMARY, summaryReportFile,
						report_ID_summary, summaryMetaFile);
			} else {
				reportSuccess &= reportTask.addTask(MetaPythonReportTask.ID_SUMMARY, summaryReportFile,
						report_ID_summary);
			}
		}

		// Add peptide report
		if (final_pep_txt.exists() && final_pep_txt.length() > 0) {
			File report_peptides_summary = new File(reportHtmlDir, "report_peptides_summary.html");
			if (summaryMetaFile != null && summaryMetaFile.exists()) {
				reportSuccess &= reportTask.addTask(MetaPythonReportTask.PEPTIDE, final_pep_txt,
						report_peptides_summary, summaryMetaFile);
			} else {
				reportSuccess &= reportTask.addTask(MetaPythonReportTask.PEPTIDE, final_pep_txt,
						report_peptides_summary);
			}
		}

		// Add protein report
		if (final_pro_txt.exists() && final_pro_txt.length() > 0) {
			File report_proteinGroups_summary = new File(reportHtmlDir, "report_proteinGroups_summary.html");
			if (summaryMetaFile != null && summaryMetaFile.exists()) {
				reportSuccess &= reportTask.addTask(MetaPythonReportTask.PROTEIN, final_pro_txt,
						report_proteinGroups_summary, summaryMetaFile);
			} else {
				reportSuccess &= reportTask.addTask(MetaPythonReportTask.PROTEIN, final_pro_txt,
						report_proteinGroups_summary);
			}
		}

		// Add taxonomy report
		File taxaReportFile = new File(metaPar.getDbSearchResultFile(), "taxonomy_analysis\\Taxa_report.tsv");
		if (taxaReportFile.exists() && taxaReportFile.length() > 0) {
			File report_taxonomy_summary = new File(reportHtmlDir, "report_taxonomy_summary.html");
			if (summaryMetaFile != null && summaryMetaFile.exists()) {
				reportSuccess &= reportTask.addTask(MetaPythonReportTask.TAXONOMY, taxaReportFile,
						report_taxonomy_summary, summaryMetaFile);
			} else {
				reportSuccess &= reportTask.addTask(MetaPythonReportTask.TAXONOMY, taxaReportFile,
						report_taxonomy_summary);
			}
		}

		File genomeHtmlFile = new File(metaPar.getDbSearchResultFile(), "taxonomy_analysis\\Genome.html");
		if (genomeHtmlFile.exists()) {
			File reportGenomeHtml = new File(reportHtmlDir, "Genome.html");
			if (reportGenomeHtml.getParentFile() != null && !reportGenomeHtml.getParentFile().exists()) {
				reportGenomeHtml.getParentFile().mkdirs();
			}
			try {
				Files.copy(genomeHtmlFile.toPath(), reportGenomeHtml.toPath(), StandardCopyOption.REPLACE_EXISTING);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				logger.error("Copied Genome.html to reports directory failed", e);
			}
			logger.info("Copied Genome.html to reports directory: " + reportGenomeHtml.getAbsolutePath());
		} else {
			logger.info("Genome.html was not exist in " + genomeHtmlFile.getAbsolutePath());
		}

		// Add function report
		File funReportTsv = new File(metaPar.getDbSearchResultFile(), "functional_annotation");
		if (funReportTsv.exists()) {
			File report_function_summary = new File(reportHtmlDir, "report_function_summary.html");
			if (summaryMetaFile != null && summaryMetaFile.exists()) {
				reportSuccess &= reportTask.addTask(MetaPythonReportTask.FUNCTION, funReportTsv,
						report_function_summary, summaryMetaFile);
			} else {
				reportSuccess &= reportTask.addTask(MetaPythonReportTask.FUNCTION, funReportTsv,
						report_function_summary);
			}
		}

		File funcHtmlFile = new File(metaPar.getDbSearchResultFile(), "functional_annotation\\functions.html");
		if (funcHtmlFile.exists()) {
			File reportFunctionsHtml = new File(reportHtmlDir, "functions.html");
			try {
				Files.copy(funcHtmlFile.toPath(), reportFunctionsHtml.toPath(), StandardCopyOption.REPLACE_EXISTING);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				logger.error("Copied functions.html to reports directory failed", e);
			}
			logger.info("Copied functions.html to reports directory: " + reportFunctionsHtml.getAbsolutePath());
		} else {
			logger.info("functions.html was not exist in " + funcHtmlFile.getAbsolutePath());
		}

		if (!reportSuccess) {
			logger.warn(taskName + ": Some report tasks could not be added");
			System.out.println(format.format(new Date()) + "\t" + taskName + ": Some report tasks could not be added");
		}

		// Run the report generation task
		Thread reportThread = new Thread(reportTask);
		reportThread.start();
		try {
			reportThread.join();
		} catch (InterruptedException e) {
			logger.error(taskName + ": Report generation thread was interrupted", e);
			System.out.println(
					format.format(new Date()) + "\t" + taskName + ": Report generation thread was interrupted");
			Thread.currentThread().interrupt();
			return false;
		}

		try {
			boolean finish = reportTask.get();
			if (finish) {
				logger.info(taskName + ": report generation successes");
				System.out.println(format.format(new Date()) + "\t" + taskName + ": report generation successes");
				return true;
			}
		} catch (InterruptedException | ExecutionException e) {
			logger.error(taskName + ": error in the report generation", e);
			System.out.println(format.format(new Date()) + "\t" + taskName + ": error in the report generation");
			return false;
		}

		return true;
	}

	/**
	 * Generate bubble chart report from tree.js data using Python.
	 */
	private void generateBubbleChartReport(File treeJsFile, File outputFile) {
		try {
			logger.info(getTaskName() + ": Generating bubble chart report from " + treeJsFile.getName());
			System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": Generating bubble chart report");

			// Use MetaPythonReportTask to generate the bubble chart
			MetaPythonReportTask reportTask = new MetaPythonReportTask(1, mrp, logger);
			boolean success = reportTask.addTask(MetaPythonReportTask.BUBBLE_CHART, treeJsFile, outputFile);

			if (success) {
				// Wait for completion
				reportTask.run();
				boolean result = reportTask.get();

				if (result) {
					logger.info(
							getTaskName() + ": Bubble chart report generated successfully: " + outputFile.getName());
					System.out.println(format.format(new Date()) + "\t" + getTaskName()
							+ ": Bubble chart report generated successfully");

					// Add to manifest
					File reportsDir = new File(metaPar.getReportDir(), "reports");
					ReportManifestGenerator manifestGen = new ReportManifestGenerator(reportsDir);
					manifestGen.addReportEntry("Bubble Chart Report", "bubble_chart_report.html",
							"Interactive hierarchical bubble chart showing taxonomic abundance with treemap, sunburst, and distribution visualizations.",
							new String[] { "bubble-chart", "taxonomy", "hierarchical", "treemap" });
				} else {
					logger.error(getTaskName() + ": Failed to generate bubble chart report");
					System.err.println(format.format(new Date()) + "\t" + getTaskName()
							+ ": Failed to generate bubble chart report");
				}
			} else {
				logger.error(getTaskName() + ": Failed to add bubble chart task");
				System.err.println(
						format.format(new Date()) + "\t" + getTaskName() + ": Failed to add bubble chart task");
			}
		} catch (Exception e) {
			logger.error(getTaskName() + ": Error generating bubble chart report", e);
			System.err.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": Error generating bubble chart report: " + e.getMessage());
		}
	}

	/**
	 * Task for post-analysis processing in DDA MAG workflows. Handles protein
	 * inference, filtering, and report generation with proper progress tracking.
	 */
	public static class DDAMagPostAnalysisTask extends MetaAbstractTask {

		private MetaPeptide[] quanPeps;
		private MetaPSM[] psms;
		private HashMap<String, Integer> fileIdMap;
		private ArrayList<String> fileList;
		private String quanMode;
		private String[] intensityTitle;
		private String taskName;
		private String decoySymbol;
		private MetaPilotMagTask parentTask;

		public DDAMagPostAnalysisTask(MetaBaseParameter metaPar, MetaResourceParameters mrp, ProgressBar progressBar,
				MetaPeptide[] quanPeps, MetaPSM[] psms, HashMap<String, Integer> fileIdMap, ArrayList<String> fileList,
				String quanMode, String[] intensityTitle, String taskName, String decoySymbol,
				MetaPilotMagTask parentTask) {
			super(metaPar, mrp, progressBar, parentTask.getLogger());
			this.quanPeps = quanPeps;
			this.psms = psms;
			this.fileIdMap = fileIdMap;
			this.fileList = fileList;
			this.quanMode = quanMode;
			this.intensityTitle = intensityTitle;
			this.taskName = taskName;
			this.decoySymbol = decoySymbol;
			this.parentTask = parentTask;
		}

		@Override
		protected String getTaskName() {
			return taskName + " Post-Analysis";
		}

		@Override
		protected Boolean call() throws Exception {
			updateMessage("Starting post-analysis processing...");

			// Build PSM maps from the PSM array
			updateMessage("Building PSM maps...");
			updateProgress(0.05, 1.0);

			DdaMagPostAnalysisService.Result analysis = DdaMagPostAnalysisService.analyze(
					new DdaMagPostAnalysisService.Input(quanPeps, psms, fileIdMap, fileList, decoySymbol,
							parentTask.magDb.getIdentifier()));
			parentTask.lastDiaStyleGenomeStats = analysis.genomeDiaStats;
			parentTask.lastDiaStyleRuninfo = analysis.genomeScoring.diaStyleRuninfo;

			HashMap<String, MetaPSM> psmMap = analysis.psmMap;
			HashMap<String, int[]> pepPsmCountMap = analysis.pepPsmCountMap;
			HashMap<String, int[]> pepChargeMap = analysis.pepChargeMap;
			HashMap<String, Integer> psmCountMap = analysis.psmCountMap;
			int[] filePsmCount = analysis.filePsmCount;
			HashSet<String>[] fileSequenceSets = analysis.fileSequenceSets;
			HashMap<String, Double> genomeScoreMap = analysis.genomeScoreMap;
			HashMap<String, double[]> genomeIntensityMap = analysis.genomeIntensityMap;
			HashMap<String, double[]> proteinIntensityMap = analysis.proteinIntensityMap;
			HashMap<String, Double> proteinTotalIntenMap = analysis.proteinTotalIntenMap;
			HashMap<String, double[]> allProteinIntensityMap = analysis.allProteinIntensityMap;
			HashMap<String, Double> allProteinTotalIntenMap = analysis.allProteinTotalIntenMap;
			HashMap<String, HashSet<String>> allProteinPepCountMap = analysis.allProteinPepCountMap;
			HashMap<String, Integer> allProteinPsmCountMap = analysis.allProteinPsmCountMap;
			HashMap<String, double[]> allProteinScoreMap = analysis.allProteinScoreMap;
			HashMap<String, HashSet<String>> totalPepCountMap = analysis.totalPepCountMap;
			HashMap<String, HashSet<String>> razorPepCountMap = analysis.razorPepCountMap;
			HashMap<String, HashSet<String>> uniquePepCountMap = analysis.uniquePepCountMap;
			HashMap<String, HashSet<String>> totalPepWithILMap = analysis.totalPepWithILMap;
			HashMap<String, HashSet<String>> razorPepWithILMap = analysis.razorPepWithILMap;
			HashMap<String, HashSet<String>> uniquePepWithILMap = analysis.uniquePepWithILMap;
			ArrayList<MetaPeptide> peptideList = analysis.peptideList;
			String[] proteins = analysis.proteins;
			String[] allProteins = analysis.allProteins;
			HashMap<String, Integer> proteinIdMap = analysis.proteinIdMap;
			HashMap<String, Integer> razorPsmCountMap = analysis.razorPsmCountMap;
			HashMap<String, String> razorProteinMap = parentTask.buildRazorProteinMap(quanPeps, analysis.razorProteins);

			updateMessage("Refining genomes and proteins...");
			updateProgress(0.15, 1.0);
			parentTask.getLogger().info("{}: genome scoring={}, genomes={}", getTaskName(),
					analysis.genomeScoring.diaStyle ? "PSM-score/top-k/q-value" : "legacy fallback",
					genomeScoreMap.size());

			updateMessage("Sorting proteins and generating reports...");
			updateProgress(0.60, 1.0);
			// Generate identification summary
			updateMessage("Generating identification summary...");
			updateProgress(0.70, 1.0);
			if (!parentTask.generateIdentificationSummary(fileIdMap, filePsmCount, fileSequenceSets, taskName)) {
				return false;
			}

			// Generate peptide report
			updateMessage("Generating peptide report...");
			updateProgress(0.75, 1.0);
			if (!parentTask.generatePeptideReport(peptideList, psmMap, pepPsmCountMap, pepChargeMap, proteinIdMap,
					fileList, quanMode, intensityTitle, taskName, razorProteinMap)) {
				return false;
			}

			// Generate protein report
			updateMessage("Generating protein report...");
			updateProgress(0.80, 1.0);
			if (!parentTask.generateProteinReport(proteins, proteinIntensityMap, proteinTotalIntenMap, totalPepCountMap,
					razorPepCountMap, psmCountMap, proteinIdMap, quanMode, intensityTitle, taskName, peptideList,
					allProteinScoreMap)) {
				return false;
			}
			if (!parentTask.generateAllProteinReport(allProteins, allProteinIntensityMap, allProteinTotalIntenMap,
					allProteinPepCountMap, razorPepCountMap, allProteinPsmCountMap, quanMode, intensityTitle, taskName,
					peptideList, allProteinScoreMap)) {
				return false;
			}

			// Functional annotation (use razor peptide/PSM counts for functions.tsv)
			updateMessage("Performing functional annotation...");
			updateProgress(0.85, 1.0);
			if (!parentTask.performFunctionalAnnotation(proteins, proteinIntensityMap, intensityTitle, taskName,
					razorPepCountMap, razorPsmCountMap)) {
				return false;
			}

			// Taxonomy analysis
			updateMessage("Performing taxonomy analysis...");
			updateProgress(0.90, 1.0);
			if (!parentTask.performTaxonomyAnalysis(genomeIntensityMap, genomeScoreMap, psmCountMap, totalPepCountMap,
					razorPepCountMap, uniquePepCountMap, totalPepWithILMap, razorPepWithILMap, uniquePepWithILMap,
					intensityTitle, taskName)) {
				return false;
			}

			// Generate final reports
			updateMessage("Generating final HTML reports...");
			updateProgress(0.95, 1.0);
			if (!parentTask.generateFinalReports(taskName)) {
				return false;
			}

			// Complete
			updateMessage("Report generation completed!");
			updateProgress(1.0, 1.0);

			return true;
		}
	}

}
