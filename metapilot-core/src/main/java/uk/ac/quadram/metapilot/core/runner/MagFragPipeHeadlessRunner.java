package uk.ac.quadram.metapilot.core.runner;

import uk.ac.quadram.metapilot.core.fragpipe.FragPipeHeadlessRunner;
import uk.ac.quadram.metapilot.core.fragpipe.FragpipeWorkflowFile;
import uk.ac.quadram.metapilot.core.mag.MagGenomeSelectionService;
import uk.ac.quadram.metapilot.core.mag.MagResultLayoutUtil;
import uk.ac.quadram.metapilot.core.prodb.DecoyShuffler;
import uk.ac.quadram.metapilot.core.prodb.FastaManager;
import uk.ac.quadram.metapilot.core.report.MetaPilotHtmlReportService;
import uk.ac.quadram.metapilot.core.resources.MetaPilotResourceConfig;
import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;
import uk.ac.quadram.metapilot.dbSearch.fragpipe.FragpipePSMTsvReader;
import uk.ac.quadram.metapilot.task.dda.mag.par.DDAMagFragPipeParameters;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;
import uk.ac.quadram.metapilot.task.io.pep.MetaPSM;
import uk.ac.quadram.metapilot.task.par.MetaData;
import uk.ac.quadram.metapilot.task.par.MetaConstants;

import java.lang.management.ManagementFactory;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * Headless implementation of the MAG+FragPipe workflow (first CLI-ready slice).
 *
 * Implements: - marker search (Marker2/Marker4) per sample - genome selection
 * (MarkerGenomeSelector) and sample-specific FASTA creation - sample-specific
 * search per sample - combined search for quant outputs
 * (combined_modified_peptide.tsv / combined_protein.tsv)
 *
 * Uses shared core TSV/report services and skips GUI-only visualization steps.
 *
 * <p>
 * This is the single execution path for FragPipe MAG - the Windows GUI task
 * ({@code DDAMagFragPipeTask.callViaCoreRunner()}) delegates entirely to
 * {@link #run} rather than doing its own per-sample work, so this class is
 * where per-raw-file sample folders (marker search, sample-DB search) are
 * both written and resumed/skipped-if-already-done. Since this change, those
 * per-sample folders nest under {@code search_results/} (see
 * {@code MagResultLayoutUtil}) instead of sitting flat under {@code resultRoot}
 * - resolved once via {@code MagResultLayoutUtil.resolveSearchResultsRoot()}
 * near the top of {@code runInternal()}, before Step A's parallel execution
 * starts, so the one-time legacy-folder migration never races a parallel
 * worker thread. {@code combined/}, {@code taxonomy_analysis/},
 * {@code functional_annotation/}, and {@code report/} remain siblings of
 * {@code search_results/}, unaffected.
 */
public final class MagFragPipeHeadlessRunner {

	public static final String BASIC_SEARCH_WORKFLOW = "Basic-Search.workflow";
	public static final String LFQ_MBR_WORKFLOW = "LFQ-MBR.workflow";
	private static final String MARKER_SEARCH_DIR = "marker";
	private static final String SAMPLE_DB_SEARCH_DIR = "sample_db";
	private static final long BYTES_PER_GB = 1024L * 1024L * 1024L;
	private static final long MIN_MEMORY_PER_FRAGPIPE_TASK_BYTES = 60L * BYTES_PER_GB;
	private static final Object DECOY_FASTA_LOCK = new Object();

	// Coarse step-to-percent mapping for cb.onProgress(...), on a 0-100 scale
	// covering this entire run() call. Neither this class nor the lower-level
	// FragPipeHeadlessRunner process wrapper parses fragpipe.bat's own output
	// into a numeric fraction, so these step boundaries (plus per-sample
	// increments within Steps A/C, see runParallelSampleTasks) are the only
	// progress signal callers get.
	private static final double STEP_A_START_PERCENT = 5;
	private static final double STEP_A_END_PERCENT = 30;
	private static final double STEP_B_START_PERCENT = 30;
	private static final double STEP_B_END_PERCENT = 40;
	private static final double STEP_C_START_PERCENT = 40;
	private static final double STEP_C_END_PERCENT = 65;
	private static final double STEP_DEF_START_PERCENT = 65;
	private static final double STEP_DEF_END_PERCENT = 75;
	private static final double STEP_GHI_START_PERCENT = 75;
	private static final double STEP_GHI_END_PERCENT = 80;

	private MagFragPipeHeadlessRunner() {
	}

	public static final class Result {
		public final File resultRoot;
		public final File combinedDir;
		public final File combinedModifiedPeptideTsv;
		public final File combinedProteinTsv;

		Result(File resultRoot, File combinedDir, File combinedModifiedPeptideTsv, File combinedProteinTsv) {
			this.resultRoot = resultRoot;
			this.combinedDir = combinedDir;
			this.combinedModifiedPeptideTsv = combinedModifiedPeptideTsv;
			this.combinedProteinTsv = combinedProteinTsv;
		}
	}

	public static Result run(DDAMagFragPipeParameters params, File resourcesFragpipe, File workflowTemplate,
			File magdbConfigJson, ProgressCallback cb) throws Exception {
		return run(params, resourcesFragpipe, workflowTemplate, magdbConfigJson, null, cb);
	}

	public static Result run(DDAMagFragPipeParameters params, File resourcesFragpipe, File workflowTemplate,
			File magdbConfigJson, File funcDbFile, ProgressCallback cb) throws Exception {
		return run(params, resourcesFragpipe, workflowTemplate, magdbConfigJson, funcDbFile, cb, true, null);
	}

	public static Result run(DDAMagFragPipeParameters params, File resourcesFragpipe, File workflowTemplate,
			File magdbConfigJson, File funcDbFile, ProgressCallback cb, boolean generateHtmlReport,
			String reportPythonHint) throws Exception {
		return run(params, resourcesFragpipe, workflowTemplate, magdbConfigJson, funcDbFile, cb,
				generateHtmlReport, reportPythonHint, null);
	}

	public static Result run(DDAMagFragPipeParameters params, File resourcesFragpipe, File workflowTemplate,
			File magdbConfigJson, File funcDbFile, ProgressCallback cb, boolean generateHtmlReport,
			String reportPythonHint, MetaPilotResourceConfig reportResourceConfig) throws Exception {
		return runInternal(params, resourcesFragpipe, workflowTemplate, magdbConfigJson, funcDbFile, cb, true,
				generateHtmlReport, reportPythonHint, reportResourceConfig, false);
	}

	public static Result runWithInstalledDefaults(DDAMagFragPipeParameters params, File resourcesFragpipe,
			File magdbConfigJson, ProgressCallback cb) throws Exception {
		return runWithInstalledDefaults(params, resourcesFragpipe, magdbConfigJson, null, cb);
	}

	public static Result runWithInstalledDefaults(DDAMagFragPipeParameters params, File resourcesFragpipe,
			File magdbConfigJson, File funcDbFile, ProgressCallback cb) throws Exception {
		return runWithInstalledDefaults(params, resourcesFragpipe, magdbConfigJson, funcDbFile, cb, true, null);
	}

	public static Result runWithInstalledDefaults(DDAMagFragPipeParameters params, File resourcesFragpipe,
			File magdbConfigJson, File funcDbFile, ProgressCallback cb, boolean generateHtmlReport,
			String reportPythonHint) throws Exception {
		return runWithInstalledDefaults(params, resourcesFragpipe, magdbConfigJson, funcDbFile, cb,
				generateHtmlReport, reportPythonHint, null);
	}

	public static Result runWithInstalledDefaults(DDAMagFragPipeParameters params, File resourcesFragpipe,
			File magdbConfigJson, File funcDbFile, ProgressCallback cb, boolean generateHtmlReport,
			String reportPythonHint, MetaPilotResourceConfig reportResourceConfig) throws Exception {
		return runInternal(params, resourcesFragpipe, null, magdbConfigJson, funcDbFile, cb, true,
				generateHtmlReport, reportPythonHint, reportResourceConfig, false);
	}

	public static Result runWithoutPostAnalysis(DDAMagFragPipeParameters params, File resourcesFragpipe,
			File workflowTemplate, File magdbConfigJson, ProgressCallback cb) throws Exception {
		return runInternal(params, resourcesFragpipe, workflowTemplate, magdbConfigJson, null, cb, false,
				false, null, null, false);
	}

	/**
	 * Headless-only identification mode: runs Step A (marker search), Step B
	 * (genome selection + sample-specific FASTA build) and Step C
	 * (sample-specific search) and then stops. The combined FASTA build and
	 * combined quantification search (Step D/E/F, which produce
	 * combined_modified_peptide.tsv / combined_protein.tsv) and all downstream
	 * post-analysis export (Step G/H/I: final_peptides.tsv, taxonomy, functional
	 * annotation, HTML report - all of which depend on the combined quant
	 * output) are skipped entirely.
	 *
	 * This is a brand new entry point, not a flag on the existing {@code run}
	 * methods, so it cannot change behavior for any existing caller - in
	 * particular the Windows GUI task ({@code DDAMagFragPipeTask}), which only
	 * ever calls {@link #run} / {@link #runWithInstalledDefaults}, is
	 * unaffected. Only headless callers (CLI, HPC study runner) that explicitly
	 * opt in by calling this method will use identification-only mode.
	 */
	public static Result runIdentificationOnly(DDAMagFragPipeParameters params, File resourcesFragpipe,
			File workflowTemplate, File magdbConfigJson, ProgressCallback cb) throws Exception {
		return runInternal(params, resourcesFragpipe, workflowTemplate, magdbConfigJson, null, cb, false,
				false, null, null, true);
	}

	/** Identification-only mode using the FragPipe installation's own installed workflow templates. */
	public static Result runIdentificationOnlyWithInstalledDefaults(DDAMagFragPipeParameters params,
			File resourcesFragpipe, File magdbConfigJson, ProgressCallback cb) throws Exception {
		return runInternal(params, resourcesFragpipe, null, magdbConfigJson, null, cb, false, false, null,
				null, true);
	}

	private static Result runInternal(DDAMagFragPipeParameters params, File resourcesFragpipe, File workflowTemplate,
			File magdbConfigJson, File funcDbFile, ProgressCallback cb, boolean exportPostAnalysis,
			boolean generateHtmlReport, String reportPythonHint, MetaPilotResourceConfig reportResourceConfig,
			boolean identificationOnly)
			throws Exception {
		if (cb == null)
			cb = new ProgressCallback.NoOp();
		if (params == null)
			throw new IllegalArgumentException("params is null");

		MetaData md = params.getMetadata();
		if (md == null || md.getRawFiles() == null || md.getRawFiles().length == 0) {
			throw new IllegalArgumentException("No raw files in parameters metadata");
		}

		// Resolve MAG DB item from MAGDB_1_0.json using catalogueId + identifier + version
		MagDbItem mag = params.getMagDbItem();
		if (mag == null || mag.getCurrentFile() == null || !mag.getCurrentFile().exists()) {
			mag = resolveMagDb(params.getMagDbCatalogueId(), params.getMagDbIdentifier(),
					params.getMagDbVersion(), magdbConfigJson);
		}
		if (mag == null) {
			throw new IllegalArgumentException("MAG database not found in " + magdbConfigJson + " for catalogue="
					+ params.getMagDbCatalogueId() + " identifier=" + params.getMagDbIdentifier()
					+ " version=" + params.getMagDbVersion());
		}
		params.setMagDbItem(mag);

		cb.onMessage("Using MAG DB: " + mag.getCatalogueID() + " [" + mag.getIdentifier() + "] / "
				+ mag.getUsedVersion().getCatalogueVersion());
		cb.onMessage(
				"MAG DB folder: " + (mag.getCurrentFile() != null ? mag.getCurrentFile().getAbsolutePath() : "null"));

		// Prepare output folders
		File resultRoot = new File(params.getResult(), "fragpipe_MAG_search");
		resultRoot.mkdirs();
		cb.onMessage("Result root: " + resultRoot.getAbsolutePath());

		// Per-raw-file sample folders (marker search, sample-DB search) now nest
		// under search_results/ instead of sitting flat as siblings of combined/,
		// taxonomy_analysis/, functional_annotation/, report/ - matching pFind's
		// layout (see MagResultLayoutUtil). Resolved exactly once, here,
		// synchronously and before Step A's parallel execution begins, so the
		// one-time legacy-folder migration never races a parallel worker thread.
		File searchResultsRoot = MagResultLayoutUtil.resolveSearchResultsRoot(resultRoot);
		cb.onMessage("Search results root: " + searchResultsRoot.getAbsolutePath());

		String[] rawFiles = md.getRawFiles();
		String[] expNames = md.getExpNames();
		File spectraDirForSummary = null;
		if (MzmlSummaryFileWriter.allMzml(rawFiles)) {
			spectraDirForSummary = params.getSpectraFile();
			if (MzmlSummaryFileWriter.ensure(spectraDirForSummary, rawFiles)) {
				cb.onMessage("Prepared mzML summary files: " + spectraDirForSummary.getAbsolutePath());
			} else {
				cb.onMessage(
						"mzML summary files could not be prepared; final_summary.tsv will use identification fallback");
				spectraDirForSummary = null;
			}
		}
		int threads = Math.max(1, params.getThreadCount());
		int requestedParallelTasks = Math.max(0, params.getFragpipeParallelTasks());
		ParallelTaskPlan parallelPlan = determineFragpipeParallelTaskPlan(rawFiles.length, threads, requestedParallelTasks);
		int parallelTasks = parallelPlan.parallelTasks;
		int perTaskThreads = Math.max(1, threads / parallelTasks);
		if (requestedParallelTasks > 0 && parallelTasks > parallelPlan.memorySafeParallelTasks) {
			cb.onMessage("WARNING: requested FragPipe parallel tasks (" + requestedParallelTasks
					+ ") exceeds the memory-based safe estimate (" + parallelPlan.memorySafeParallelTasks + "). "
					+ "Each FragPipe task is assumed to need " + (MIN_MEMORY_PER_FRAGPIPE_TASK_BYTES / BYTES_PER_GB)
					+ " GB. This may cause the job to fail with an out-of-memory error.");
		}
		if (requestedParallelTasks > 0 && parallelTasks > threads) {
			cb.onMessage("WARNING: requested FragPipe parallel tasks (" + requestedParallelTasks
					+ ") exceeds totalThreads (" + threads + "); per-task threads will be reduced to 1.");
		}
		cb.onMessage("FragPipe parallel plan: totalThreads=" + threads + ", parallelTasks=" + parallelTasks
				+ ", perTaskThreads=" + perTaskThreads + ", memoryPerTaskGb="
				+ (MIN_MEMORY_PER_FRAGPIPE_TASK_BYTES / BYTES_PER_GB) + ", memorySafeParallelTasks="
				+ parallelPlan.memorySafeParallelTasks + ", requestedParallelTasks="
				+ (requestedParallelTasks > 0 ? requestedParallelTasks : "auto") + ", detectedMemoryGb="
				+ parallelPlan.detectedMemoryGb + ", memorySource=" + parallelPlan.memorySource);

		File basicWorkflowTemplate = resolveBasicWorkflowTemplate(resourcesFragpipe, workflowTemplate);
		// Identification-only runs never reach the combined/quant search, so don't
		// require the LFQ-MBR workflow template to be resolvable/installed for them.
		File combinedWorkflowTemplate = identificationOnly ? null
				: resolveCombinedWorkflowTemplate(resourcesFragpipe, workflowTemplate);

		// Step A: per-sample marker searches
		cb.onMessage("Step A: running marker searches (per sample) ...");
		runStepAParallel(rawFiles, expNames, params, mag, resourcesFragpipe, basicWorkflowTemplate, searchResultsRoot,
				parallelTasks, perTaskThreads, cb);

		// Step B: select genomes and build sample-specific FASTAs
		cb.onMessage("Step B: selecting genomes and building sample-specific FASTAs ...");
		for (int i = 0; i < rawFiles.length; i++) {
			File sampleDir = new File(searchResultsRoot, baseNameNoExt(rawFiles[i]));
			File markerDir = new File(sampleDir, MARKER_SEARCH_DIR);

			File fpOut = fragPipeResultFolder(markerDir, expNames, i);
			File psmTsv = new File(fpOut, "psm.tsv");
			if (!psmTsv.exists()) {
				throw new IOException("Missing FragPipe psm.tsv after marker search: " + psmTsv.getAbsolutePath());
			}

			boolean useMarker4 = MagDbItem.MARKER4.equalsIgnoreCase(params.getGenomeSelectorMode());
			String prefix = useMarker4 ? "marker_genome" : "genome";
			File evidence = new File(markerDir, prefix + "_evidence.tsv");
			File ranked = new File(markerDir, prefix + "_ranked.tsv");
			File selected = new File(markerDir, prefix + "_selected.tsv");

			HashSet<String> genomeSet = selectGenomes(mag, params, psmTsv, evidence, ranked, selected, cb);

			File sampleFasta = new File(markerDir, sampleDir.getName() + ".fasta");
			if (!sampleFasta.exists() || sampleFasta.length() == 0) {
				cb.onMessage("Writing sample-specific FASTA: " + sampleFasta.getAbsolutePath() + " (" + genomeSet.size()
						+ " genomes)");
				try (PrintWriter writer = new PrintWriter(sampleFasta)) {
					File originalDbFolder = mag.getOriginalDbFolder();
					if (originalDbFolder == null || !originalDbFolder.exists()) {
						throw new IOException("Original DB folder not found in MAG DB: " + originalDbFolder);
					}
					for (String genome : genomeSet) {
						File faa = new File(originalDbFolder, genome + ".faa");
						File fa = new File(originalDbFolder, genome + ".fasta");
						File src = faa.exists() ? faa : (fa.exists() ? fa : null);
						if (src == null)
							continue;
						appendFastaWithDecoys(src, writer);
					}
				}
			}
			if (!sampleFasta.exists() || sampleFasta.length() == 0) {
				throw new IOException("Sample-specific FASTA is empty: " + sampleFasta.getAbsolutePath());
			}

			double stepBPct = STEP_B_START_PERCENT
					+ (STEP_B_END_PERCENT - STEP_B_START_PERCENT) * (i + 1) / rawFiles.length;
			cb.onProgress(stepBPct, 100);
		}

		// Step C: sample-specific searches (per sample)
		cb.onMessage("Step C: running sample-specific searches (per sample) ...");
		runStepCParallel(rawFiles, expNames, params, resourcesFragpipe, basicWorkflowTemplate, searchResultsRoot,
				parallelTasks, perTaskThreads, cb);

		if (identificationOnly) {
			// Headless identification-only mode: stop here. Skip the combined FASTA
			// build and combined quantification search (Step D/E/F) and all
			// post-analysis export (Step G/H/I), since those all depend on the
			// combined quant output (combined_modified_peptide.tsv / combined_protein.tsv)
			// that this mode intentionally never produces.
			cb.onMessage(
					"Identification-only mode: skipping combined FASTA build, combined quantification search, and post-analysis export.");
			cb.onProgress(100, 100);
			cb.onMessage("MAG+FragPipe headless run completed (identification only).");
			cb.onComplete(resultRoot.getAbsolutePath());
			return new Result(resultRoot, null, null, null);
		}

		// Step D/E/F (CLI slice): build combined FASTA and run combined search for
		// quant outputs.
		cb.onProgress(STEP_DEF_START_PERCENT, 100);
		cb.onMessage("Step D/E/F: building combined FASTA and running combined search ...");
		File combinedDir = new File(resultRoot, "combined");
		combinedDir.mkdirs();

		File combinedFasta = new File(combinedDir, "combined.fasta");
		if (!combinedFasta.exists() || combinedFasta.length() == 0) {
			cb.onMessage("Building combined FASTA: " + combinedFasta.getAbsolutePath());
			LinkedHashMap<String, String> fastaMap = new LinkedHashMap<>();

			for (int i = 0; i < rawFiles.length; i++) {
				File sampleDir = new File(searchResultsRoot, baseNameNoExt(rawFiles[i]));
				File sampleSearchDir = new File(sampleDir, SAMPLE_DB_SEARCH_DIR);
				File fpOut = fragPipeResultFolder(sampleSearchDir, expNames, i);
				File proteinFas = new File(fpOut, "protein.fas");
				if (!proteinFas.exists()) {
					throw new IOException(
							"Missing FragPipe protein.fas for combined FASTA: " + proteinFas.getAbsolutePath());
				}
				fastaMap.putAll(readFastaSequences(proteinFas));
			}

			String cleavage = getCleavageResiduesOrDefault(FragpipeWorkflowFile.parse(basicWorkflowTemplate.toPath()));
			try (PrintWriter writer = new PrintWriter(combinedFasta)) {
				for (String proName : fastaMap.keySet()) {
					String seq = fastaMap.get(proName);
					writer.println(proName);
					writer.println(seq);

					String shuffleProName = ">shuf_" + proName.substring(1);
					String shuffle = DecoyShuffler.generateShuffledDecoy(seq, cleavage, 42);
					writer.println(shuffleProName);
					writer.println(shuffle);

					String revProName = ">rev_" + proName.substring(1);
					writer.println(revProName);
					writer.println(new StringBuilder(seq).reverse().toString());

					String revShufProName = ">rev_shuf_" + proName.substring(1);
					writer.println(revShufProName);
					writer.println(new StringBuilder(shuffle).reverse().toString());
				}

				// Optional host DB append
				if (params.isAppendHostDb() && params.getHostDB() != null && !params.getHostDB().isBlank()) {
					File hostFasta = resolveUserPath(params.getHostDB());
					if (hostFasta.exists()) {
						cb.onMessage("Appending host DB to combined FASTA: " + hostFasta.getAbsolutePath());
						appendHostFastaWithDecoys(hostFasta, writer, cleavage);
					} else {
						cb.onMessage("WARNING: host DB enabled but not found: " + hostFasta.getAbsolutePath());
					}
				}
			}
		}

		File combinedManifest = new File(combinedDir, "raws.fp-manifest");
		params.getMetadata().exportFragpipeManifest(combinedManifest);

		File combinedWorkflow = new File(combinedDir, "combined.workflow");
		configureWorkflowFromParams(FragpipeWorkflowFile.parse(combinedWorkflowTemplate.toPath()), params)
				.set("database.decoy-tag", "rev_").set("percolator.run-percolator", "true")
				.set("peptide-prophet.run-peptide-prophet", "false").set("protein-prophet.run-protein-prophet", "true")
				.setDatabase(combinedFasta.getAbsolutePath()).write(combinedWorkflow.toPath());

		File quanPep = new File(combinedDir, "combined_modified_peptide.tsv");
		File quanPro = new File(combinedDir, "combined_protein.tsv");
		cleanupMarkerSearchIntermediates(rawFiles, searchResultsRoot, cb);
		if (!quanPep.exists() || quanPep.length() == 0) {
			cb.onMessage("Running combined FragPipe search (quant outputs) ...");
			int exit = FragPipeHeadlessRunner.run(resourcesFragpipe, combinedWorkflow, combinedManifest, combinedDir,
					threads, cb);
			if (exit != 0) {
				throw new IOException("FragPipe combined search failed (exit=" + exit + ")");
			}
			cleanupCombinedProtXml(combinedDir, cb);
		} else {
			cb.onMessage("Combined peptide output already exists: " + quanPep.getName());
			cleanupCombinedProtXml(combinedDir, cb);
		}

		// Step G (headless): shared post-analysis and report exports.
		if (exportPostAnalysis) {
			cb.onProgress(STEP_GHI_START_PERCENT, 100);
			cb.onMessage("Step G: post-analysis + TSV export ...");
			FragpipePostAnalysisTsvExporter.Result post = FragpipePostAnalysisTsvExporter.exportFinalTsvs(mag, quanPep,
					resultRoot, spectraDirForSummary, "DDA MAG FragPipe", cb);
			cb.onProgress(STEP_GHI_START_PERCENT + (STEP_GHI_END_PERCENT - STEP_GHI_START_PERCENT) / 3, 100);

			cb.onMessage("Step H: taxonomy report export ...");
			exportTaxonomyTsvs(params, quanPep, mag, resultRoot, post, cb);
			cb.onProgress(STEP_GHI_START_PERCENT + (STEP_GHI_END_PERCENT - STEP_GHI_START_PERCENT) * 2 / 3, 100);

			if (funcDbFile != null && funcDbFile.isFile()) {
				cb.onMessage("Step I: functional annotation report export ...");
				MagFunctionalAnnotationService.Request request = new MagFunctionalAnnotationService.Request();
				request.resultDir = resultRoot;
				request.magDb = mag;
				request.microDbCatalogueId = params.getMagDbCatalogueId();
				request.functionDefinitionDb = funcDbFile;
				request.analysisType = MetaConstants.fragpipeMAG;
				request.quantitativeMode = params.getQuantitativeMode();
				request.intensityTitles = post.intensityTitles;
				request.metaProteins = post.metaProteins;
				MagFunctionalAnnotationService.run(request);
			} else {
				cb.onMessage("FragPipe functional annotation: skipped (func_def.db not configured)");
			}
			runHtmlReportIfRequested(resultRoot, params, generateHtmlReport, reportPythonHint,
					reportResourceConfig, cb);
			cb.onProgress(STEP_GHI_END_PERCENT, 100);
		}

		cb.onMessage("MAG+FragPipe headless run completed.");
		cb.onComplete(resultRoot.getAbsolutePath());
		return new Result(resultRoot, combinedDir, quanPep, quanPro);
	}

	private static void runHtmlReportIfRequested(File resultDir, DDAMagFragPipeParameters params,
			boolean generateHtmlReport, String reportPythonHint, MetaPilotResourceConfig reportResourceConfig,
			ProgressCallback cb) {
		if (!generateHtmlReport) {
			cb.onMessage("HTML report generation skipped (disabled)");
			return;
		}
		cb.onMessage("MAG+FragPipe: generating shared HTML report");
		MetaPilotHtmlReportService.Result report = MetaPilotHtmlReportService.runStandardReportWarningOnly(
				resultDir,
				"DDA_MAG_FRAGPIPE",
				params != null ? params.getThreadCount() : 1,
				reportPythonHint,
				reportResourceConfig,
				cb);
		if (!report.warnings.isEmpty()) {
			cb.onMessage("WARNING: HTML report generation completed with warnings: "
					+ String.join("; ", report.warnings));
		}
	}

	private static MagDbItem resolveMagDb(String catalogueId, String identifier, String version, File magdbConfigJson) {
		return uk.ac.quadram.metapilot.task.mag.MagDbItemResolver.resolve(
				catalogueId, identifier, version, magdbConfigJson);
	}

	private static void writeManifest(File manifest, String rawFile, String expName) throws IOException {
		manifest.getParentFile().mkdirs();
		try (PrintWriter pw = new PrintWriter(manifest)) {
			// file_path \t experiment_name \t fraction \t DDA
			pw.println(rawFile + "\t" + expName.replaceAll("-", "_") + "\t1\tDDA");
		}
	}

	private static String safeExpName(String[] expNames, int idx) {
		if (expNames == null || idx < 0 || idx >= expNames.length)
			return "sample_" + (idx + 1);
		String s = expNames[idx];
		if (s == null || s.isBlank())
			return "sample_" + (idx + 1);
		return s.trim();
	}

	private static String baseNameNoExt(String path) {
		if (path == null)
			return "sample";
		String name = new File(path).getName();
		int dot = name.lastIndexOf('.');
		return dot > 0 ? name.substring(0, dot) : name;
	}

	private static void runStepAParallel(String[] rawFiles, String[] expNames, DDAMagFragPipeParameters params,
			MagDbItem mag, File resourcesFragpipe, File basicWorkflowTemplate, File searchResultsRoot, int parallelTasks,
			int perTaskThreads, ProgressCallback cb) throws Exception {
		runParallelSampleTasks(rawFiles, parallelTasks, "FragPipe-StepA", STEP_A_START_PERCENT, STEP_A_END_PERCENT, cb,
				rawFile -> {
			File sampleDir = new File(searchResultsRoot, baseNameNoExt(rawFile));
			File markerDir = new File(sampleDir, MARKER_SEARCH_DIR);
			markerDir.mkdirs();

			int sampleIndex = findRawFileIndex(rawFiles, rawFile);
			File existingPsm = new File(fragPipeResultFolder(markerDir, expNames, sampleIndex), "psm.tsv");
			if (hasContent(existingPsm)) {
				cb.onMessage("Marker search already exists, skipping: " + existingPsm.getAbsolutePath());
				cleanupCombinedProtXml(markerDir, cb);
				return;
			}

			File manifest = new File(sampleDir, "raws.fp-manifest");
			writeManifest(manifest, rawFile, safeExpName(expNames, sampleIndex));

			String genomeSelectorMode = params.getGenomeSelectorMode();
			File hapFasta = mag.getHapFastaFile(genomeSelectorMode);
			File hapTdFasta = mag.getHapFastaTDFile(genomeSelectorMode);
			if (hapFasta == null || !hapFasta.exists()) {
				throw new IOException("HAP FASTA not found for selectorMode=" + genomeSelectorMode + ": " + hapFasta);
			}
			if (hapTdFasta == null) {
				throw new IOException("HAP TD FASTA path is null (MAG DB not configured correctly)");
			}
			ensureTargetDecoyFasta(hapFasta, hapTdFasta, cb);

			File wfOut = new File(markerDir, "marker.workflow");
			configureWorkflowFromParams(FragpipeWorkflowFile.parse(basicWorkflowTemplate.toPath()), params)
					.set("database.decoy-tag", "rev_").set("percolator.run-percolator", "true")
					.set("peptide-prophet.run-peptide-prophet", "false")
					.set("protein-prophet.run-protein-prophet", "true").setDatabase(hapTdFasta.getAbsolutePath())
					.write(wfOut.toPath());

			cb.onMessage("Marker search: " + rawFile + " (threads=" + perTaskThreads + ")");
			int exit = FragPipeHeadlessRunner.run(resourcesFragpipe, wfOut, manifest, markerDir, perTaskThreads, cb);
			File psmAfterRun = new File(fragPipeResultFolder(markerDir, expNames, sampleIndex), "psm.tsv");
			if (exit != 0) {
				if (hasContent(psmAfterRun)) {
					cb.onMessage("WARNING: FragPipe marker search returned exit=" + exit
							+ " but produced non-empty psm.tsv; continuing with marker evidence: "
							+ psmAfterRun.getAbsolutePath());
					cleanupCombinedProtXml(markerDir, cb);
					return;
				}
				throw new IOException("FragPipe marker search failed (exit=" + exit + ") for " + rawFile + "; workDir="
						+ markerDir.getAbsolutePath() + "; workflow=" + wfOut.getAbsolutePath() + "; manifest="
						+ manifest.getAbsolutePath());
			}
			cleanupCombinedProtXml(markerDir, cb);
		});
	}

	private static void ensureTargetDecoyFasta(File targetFasta, File targetDecoyFasta, ProgressCallback cb)
			throws IOException {
		if (targetFasta == null || !targetFasta.exists()) {
			throw new IOException("Target FASTA not found: " + targetFasta);
		}
		if (targetDecoyFasta == null) {
			throw new IOException("Target-decoy FASTA path is null");
		}
		synchronized (DECOY_FASTA_LOCK) {
			if (isValidTargetDecoyFasta(targetFasta, targetDecoyFasta)) {
				return;
			}
			File parent = targetDecoyFasta.getParentFile();
			if (parent != null && !parent.exists()) {
				parent.mkdirs();
			}
			File temp = new File(targetDecoyFasta.getParentFile(),
					targetDecoyFasta.getName() + ".tmp." + System.currentTimeMillis());
			cb.onMessage("Creating target-decoy FASTA: " + targetDecoyFasta.getAbsolutePath() + " from "
					+ targetFasta.getName());
			try {
				FastaManager.addDecoy(targetFasta.getAbsolutePath(), temp.getAbsolutePath());
				if (!isValidTargetDecoyFasta(targetFasta, temp)) {
					throw new IOException("Generated target-decoy FASTA failed validation: " + temp.getAbsolutePath());
				}
				moveReplacing(temp, targetDecoyFasta);
			} catch (IOException e) {
				Files.deleteIfExists(temp.toPath());
				throw e;
			}
		}
	}

	private static void moveReplacing(File source, File target) throws IOException {
		try {
			Files.move(source.toPath(), target.toPath(), StandardCopyOption.REPLACE_EXISTING,
					StandardCopyOption.ATOMIC_MOVE);
		} catch (IOException atomicMoveFailure) {
			Files.move(source.toPath(), target.toPath(), StandardCopyOption.REPLACE_EXISTING);
		}
	}

	private static boolean isValidTargetDecoyFasta(File targetFasta, File targetDecoyFasta) {
		if (targetDecoyFasta == null || !targetDecoyFasta.exists() || targetDecoyFasta.length() == 0) {
			return false;
		}
		if (targetDecoyFasta.lastModified() < targetFasta.lastModified()) {
			return false;
		}
		try {
			int targetHeaders = countFastaHeaders(targetFasta);
			TargetDecoyCounts tdCounts = countTargetDecoyHeaders(targetDecoyFasta);
			return targetHeaders > 0 && tdCounts.targetCount == targetHeaders && tdCounts.decoyCount == targetHeaders;
		} catch (IOException e) {
			return false;
		}
	}

	private static int countFastaHeaders(File fasta) throws IOException {
		int count = 0;
		try (BufferedReader reader = new BufferedReader(new FileReader(fasta))) {
			String line;
			while ((line = reader.readLine()) != null) {
				if (line.startsWith(">")) {
					count++;
				}
			}
		}
		return count;
	}

	private static TargetDecoyCounts countTargetDecoyHeaders(File fasta) throws IOException {
		int targetCount = 0;
		int decoyCount = 0;
		try (BufferedReader reader = new BufferedReader(new FileReader(fasta))) {
			String line;
			while ((line = reader.readLine()) != null) {
				if (!line.startsWith(">")) {
					continue;
				}
				if (line.startsWith(">rev_") || line.startsWith(">REV_")) {
					decoyCount++;
				} else {
					targetCount++;
				}
			}
		}
		return new TargetDecoyCounts(targetCount, decoyCount);
	}

	private static final class TargetDecoyCounts {
		final int targetCount;
		final int decoyCount;

		TargetDecoyCounts(int targetCount, int decoyCount) {
			this.targetCount = targetCount;
			this.decoyCount = decoyCount;
		}
	}

	private static void runStepCParallel(String[] rawFiles, String[] expNames, DDAMagFragPipeParameters params,
			File resourcesFragpipe, File basicWorkflowTemplate, File searchResultsRoot, int parallelTasks, int perTaskThreads,
			ProgressCallback cb) throws Exception {
		runParallelSampleTasks(rawFiles, parallelTasks, "FragPipe-StepC", STEP_C_START_PERCENT, STEP_C_END_PERCENT, cb,
				rawFile -> {
			File sampleDir = new File(searchResultsRoot, baseNameNoExt(rawFile));
			File markerDir = new File(sampleDir, MARKER_SEARCH_DIR);
			File sampleFasta = new File(markerDir, sampleDir.getName() + ".fasta");

			File sampleSearchDir = new File(sampleDir, SAMPLE_DB_SEARCH_DIR);
			sampleSearchDir.mkdirs();

			int sampleIndex = findRawFileIndex(rawFiles, rawFile);
			File existingProteinFas = new File(fragPipeResultFolder(sampleSearchDir, expNames, sampleIndex),
					"protein.fas");
			if (hasContent(existingProteinFas)) {
				cb.onMessage(
						"Sample-specific search already exists, skipping: " + existingProteinFas.getAbsolutePath());
				cleanupCombinedProtXml(sampleSearchDir, cb);
				return;
			}

			File manifest = new File(sampleDir, "raws.fp-manifest");
			writeManifest(manifest, rawFile, safeExpName(expNames, sampleIndex));

			File wfOut = new File(sampleSearchDir, "sample_db.workflow");
			configureWorkflowFromParams(FragpipeWorkflowFile.parse(basicWorkflowTemplate.toPath()), params)
					.set("database.decoy-tag", "rev_").set("percolator.run-percolator", "true")
					.set("peptide-prophet.run-peptide-prophet", "false")
					.set("protein-prophet.run-protein-prophet", "true").setDatabase(sampleFasta.getAbsolutePath())
					.write(wfOut.toPath());

			cb.onMessage("Sample-specific search: " + rawFile + " (threads=" + perTaskThreads + ")");
			int exit = FragPipeHeadlessRunner.run(resourcesFragpipe, wfOut, manifest, sampleSearchDir, perTaskThreads,
					cb);
			if (exit != 0) {
				throw new IOException("FragPipe sample-specific search failed (exit=" + exit + ") for " + rawFile);
			}
			cleanupCombinedProtXml(sampleSearchDir, cb);
		});
	}

	private static void cleanupMarkerSearchIntermediates(String[] rawFiles, File searchResultsRoot, ProgressCallback cb) {
		if (rawFiles == null || searchResultsRoot == null) {
			return;
		}
		for (String rawFile : rawFiles) {
			File markerDir = new File(new File(searchResultsRoot, baseNameNoExt(rawFile)), MARKER_SEARCH_DIR);
			deleteMatchingFiles(markerDir, cb, ".pepindex");

			File[] markerResultDirs = markerDir.listFiles(File::isDirectory);
			if (markerResultDirs == null) {
				continue;
			}
			for (File markerResultDir : markerResultDirs) {
				deleteMatchingFiles(markerResultDir, cb, ".xml");
			}
		}
	}

	private static void deleteMatchingFiles(File directory, ProgressCallback cb, String... suffixes) {
		if (directory == null || suffixes == null) {
			return;
		}
		File[] files = directory.listFiles(File::isFile);
		if (files == null) {
			return;
		}
		for (File file : files) {
			String fileName = file.getName().toLowerCase(java.util.Locale.ROOT);
			boolean matches = false;
			for (String suffix : suffixes) {
				String expectedSuffix = suffix.toLowerCase(java.util.Locale.ROOT);
				if (fileName.endsWith(expectedSuffix)) {
					matches = true;
					break;
				}
			}
			if (!matches) {
				continue;
			}
			try {
				if (Files.deleteIfExists(file.toPath())) {
					cb.onMessage("Deleted FragPipe marker intermediate: " + file.getAbsolutePath());
				}
			} catch (IOException e) {
				cb.onMessage("WARNING: could not delete FragPipe marker intermediate "
						+ file.getAbsolutePath() + ": " + e.getMessage());
			}
		}
	}

	private static void cleanupCombinedProtXml(File searchDir, ProgressCallback cb) {
		if (searchDir == null) {
			return;
		}
		File combinedProtXml = new File(searchDir, "combined.prot.xml");
		try {
			if (Files.deleteIfExists(combinedProtXml.toPath())) {
				cb.onMessage("Deleted FragPipe intermediate: " + combinedProtXml.getAbsolutePath());
			}
		} catch (IOException e) {
			cb.onMessage("WARNING: could not delete FragPipe intermediate "
					+ combinedProtXml.getAbsolutePath() + ": " + e.getMessage());
		}
	}

	/**
	 * Runs {@code task} once per raw file (in parallel if {@code parallelTasks
	 * > 1}), reporting coarse progress to {@code cb.onProgress(...)} as each
	 * sample completes - mapped linearly across {@code [stepStartPercent,
	 * stepEndPercent]} of the overall 0-100 run. This is the only progress
	 * signal Step A/C produce: neither this class nor the lower-level {@code
	 * FragPipeHeadlessRunner} process wrapper parses per-sample fragpipe.bat
	 * output into a numeric progress fraction, so without this the GUI
	 * progress bar previously never moved during the entire run (it only ever
	 * received {@code onMessage} text, never {@code onProgress}).
	 */
	private static void runParallelSampleTasks(String[] rawFiles, int parallelTasks, String threadNamePrefix,
			double stepStartPercent, double stepEndPercent, ProgressCallback cb, SampleTask task) throws Exception {
		int total = rawFiles.length;
		AtomicInteger completedCount = new AtomicInteger(0);
		SampleTask reportingTask = rawFile -> {
			task.run(rawFile);
			int completed = completedCount.incrementAndGet();
			double pct = stepStartPercent + (stepEndPercent - stepStartPercent) * completed / total;
			cb.onProgress(pct, 100);
		};

		if (parallelTasks <= 1 || rawFiles.length <= 1) {
			for (String rawFile : rawFiles) {
				reportingTask.run(rawFile);
			}
			return;
		}

		ExecutorService executor = Executors.newFixedThreadPool(parallelTasks, r -> {
			Thread t = new Thread(r, threadNamePrefix);
			t.setDaemon(true);
			return t;
		});
		List<Future<Void>> futures = new ArrayList<>();
		try {
			for (String rawFile : rawFiles) {
				futures.add(executor.submit(() -> {
					reportingTask.run(rawFile);
					return null;
				}));
			}
			for (Future<Void> future : futures) {
				try {
					future.get();
				} catch (ExecutionException e) {
					Throwable cause = e.getCause();
					if (cause instanceof Exception) {
						throw (Exception) cause;
					}
					throw new RuntimeException(cause);
				}
			}
		} finally {
			executor.shutdown();
			try {
				if (!executor.awaitTermination(5, TimeUnit.MINUTES)) {
					executor.shutdownNow();
				}
			} catch (InterruptedException e) {
				executor.shutdownNow();
				Thread.currentThread().interrupt();
				throw e;
			}
		}
	}

	private static ParallelTaskPlan determineFragpipeParallelTaskPlan(int sampleCount, int totalThreads,
			int requestedParallelTasks) {
		MemoryInfo memoryInfo = detectFragpipeMemoryLimit();
		int memoryLimitedTasks = 1;
		if (memoryInfo.bytes > 0) {
			memoryLimitedTasks = (int) Math.max(1L, memoryInfo.bytes / MIN_MEMORY_PER_FRAGPIPE_TASK_BYTES);
		}

		int sampleCap = Math.max(1, sampleCount);
		int threadCap = Math.max(1, totalThreads);
		int memorySafeTasks = Math.max(1, Math.min(sampleCap, Math.min(threadCap, memoryLimitedTasks)));
		int parallelTasks = memorySafeTasks;
		if (requestedParallelTasks > 0) {
			parallelTasks = Math.max(1, Math.min(sampleCap, requestedParallelTasks));
		}
		return new ParallelTaskPlan(parallelTasks, memorySafeTasks, memoryInfo.gb(), memoryInfo.source);
	}

	private static MemoryInfo detectFragpipeMemoryLimit() {
		long slurmMemPerNodeMb = readPositiveLongEnv("SLURM_MEM_PER_NODE");
		if (slurmMemPerNodeMb > 0) {
			return new MemoryInfo(slurmMemPerNodeMb * 1024L * 1024L, "SLURM_MEM_PER_NODE");
		}

		long slurmMemPerCpuMb = readPositiveLongEnv("SLURM_MEM_PER_CPU");
		long slurmCpus = readPositiveLongEnv("SLURM_CPUS_ON_NODE");
		if (slurmCpus <= 0) {
			slurmCpus = readPositiveLongEnv("SLURM_CPUS_PER_TASK");
		}
		if (slurmMemPerCpuMb > 0 && slurmCpus > 0) {
			return new MemoryInfo(slurmMemPerCpuMb * slurmCpus * 1024L * 1024L, "SLURM_MEM_PER_CPU");
		}

		Long cgroupV2 = readMemoryLimitFile("/sys/fs/cgroup/memory.max", true);
		if (cgroupV2 != null) {
			return new MemoryInfo(cgroupV2.longValue(), "cgroup-v2");
		}

		Long cgroupV1 = readMemoryLimitFile("/sys/fs/cgroup/memory/memory.limit_in_bytes", false);
		if (cgroupV1 != null) {
			return new MemoryInfo(cgroupV1.longValue(), "cgroup-v1");
		}

		java.lang.management.OperatingSystemMXBean osBean = ManagementFactory.getOperatingSystemMXBean();
		if (osBean instanceof com.sun.management.OperatingSystemMXBean) {
			return new MemoryInfo(((com.sun.management.OperatingSystemMXBean) osBean).getTotalMemorySize(),
					"physical-ram");
		}
		return new MemoryInfo(-1L, "unknown");
	}

	private static long readPositiveLongEnv(String name) {
		try {
			String raw = System.getenv(name);
			if (raw == null || raw.isBlank()) {
				return -1L;
			}
			long value = Long.parseLong(raw.trim());
			return value > 0 ? value : -1L;
		} catch (Exception e) {
			return -1L;
		}
	}

	private static Long readMemoryLimitFile(String path, boolean allowMax) {
		java.io.File file = new java.io.File(path);
		if (!file.isFile()) {
			return null;
		}
		try {
			String raw = java.nio.file.Files.readString(file.toPath()).trim();
			if (raw.isEmpty() || (allowMax && "max".equalsIgnoreCase(raw))) {
				return null;
			}
			long value = Long.parseLong(raw);
			if (value <= 0 || value >= Long.MAX_VALUE / 2) {
				return null;
			}
			return Long.valueOf(value);
		} catch (Exception e) {
			return null;
		}
	}

	private static final class ParallelTaskPlan {
		final int parallelTasks;
		final int memorySafeParallelTasks;
		final long detectedMemoryGb;
		final String memorySource;

		ParallelTaskPlan(int parallelTasks, int memorySafeParallelTasks, long detectedMemoryGb, String memorySource) {
			this.parallelTasks = parallelTasks;
			this.memorySafeParallelTasks = memorySafeParallelTasks;
			this.detectedMemoryGb = detectedMemoryGb;
			this.memorySource = memorySource;
		}
	}

	private static final class MemoryInfo {
		final long bytes;
		final String source;

		MemoryInfo(long bytes, String source) {
			this.bytes = bytes;
			this.source = source;
		}

		long gb() {
			return bytes > 0 ? bytes / BYTES_PER_GB : -1L;
		}
	}

	private static int findRawFileIndex(String[] rawFiles, String rawFile) {
		for (int i = 0; i < rawFiles.length; i++) {
			if (rawFile.equals(rawFiles[i])) {
				return i;
			}
		}
		return -1;
	}

	private static File fragPipeResultFolder(File workDir, String[] expNames, int sampleIndex) {
		String expFolder = safeExpName(expNames, sampleIndex).replaceAll("-", "_") + "_1";
		return new File(workDir, expFolder);
	}

	private static boolean hasContent(File file) {
		return file != null && file.exists() && file.length() > 0;
	}

	private static FragpipeWorkflowFile configureWorkflowFromParams(FragpipeWorkflowFile workflow,
			DDAMagFragPipeParameters params) {
		if (workflow == null || params == null) {
			return workflow;
		}

		String[] enzyme = getEnzymeInfo(params.getEnzyme());
		int missedCleavages = Math.max(0, params.getMissCleavages());
		workflow.set("msfragger.misc.fragger.enzyme-dropdown-1", enzyme[0])
				.set("msfragger.search_enzyme_name_1", enzyme[0]).set("msfragger.search_enzyme_cut_1", enzyme[1])
				.set("msfragger.search_enzyme_nocut_1", enzyme[2]).set("msfragger.num_enzyme_termini", "2")
				.set("msfragger.allowed_missed_cleavage_1", String.valueOf(missedCleavages))
				.set("msfragger.precursor_mass_lower", "-20.0").set("msfragger.precursor_mass_upper", "20.0")
				.set("msfragger.precursor_true_tolerance", "20.0").set("msfragger.precursor_mass_units", "1")
				.set("msfragger.precursor_true_units", "1").set("msfragger.fragment_mass_tolerance", "0.5")
				.set("msfragger.fragment_mass_units", "0");

		boolean ionMobility = params.getMs2ScanMode() != null
				&& params.getMs2ScanMode().toUpperCase().contains("IM-MS");
		workflow.set("workflow.input.data-type.im-ms", String.valueOf(ionMobility))
				.set("workflow.input.data-type.regular-ms", String.valueOf(!ionMobility));
		return workflow;
	}

	private static String[] getEnzymeInfo(String enzyme) {
		String normalized = enzyme == null ? "" : enzyme.trim().toLowerCase();
		switch (normalized) {
		case "trypsin":
			return new String[] { "stricttrypsin", "KR", "" };
		case "trypsin/p":
			return new String[] { "trypsin", "KR", "P" };
		case "chymotrypsin":
			return new String[] { "chymotrypsin", "FWY", "" };
		case "gluc":
			return new String[] { "gluc", "DE", "" };
		case "gluc-bicarb":
			return new String[] { "gluc-bicarb", "DE", "" };
		case "lysc":
			return new String[] { "lysc", "K", "" };
		case "lysc-p":
			return new String[] { "lysc-p", "K", "P" };
		case "argc":
			return new String[] { "argc", "R", "" };
		case "aspn":
			return new String[] { "aspn", "D", "" };
		case "elastase":
			return new String[] { "elastase", "AGILV", "" };
		case "proteinasek":
			return new String[] { "proteinasek", "AEFILMVWY", "" };
		case "thermolysin":
			return new String[] { "thermolysin", "AEFILMV", "" };
		default:
			return new String[] { "stricttrypsin", "KR", "" };
		}
	}

	@FunctionalInterface
	private interface SampleTask {
		void run(String rawFile) throws Exception;
	}

	private static HashSet<String> selectGenomes(MagDbItem mag, DDAMagFragPipeParameters params, File psmTsv,
			File evidenceFile, File rankedFile, File selectedFile, ProgressCallback cb) throws Exception {
		var psmList = new ArrayList<>(new FragpipePSMTsvReader(psmTsv).read());
		int totalPsms = psmList.size();
		psmList.removeIf(psm -> psm == null || psm.isDecoy());
		cb.onMessage("FragPipe marker selection PSMs: " + psmList.size() + "/" + totalPsms + " with Is Decoy=FALSE");

		MagGenomeSelectionService.Request request = new MagGenomeSelectionService.Request();
		request.magDb = mag;
		request.genomeSelectorMode = params.getGenomeSelectorMode();
		request.genomeSelectionStrategy = params.getGenomeSelectionStrategy();
		request.psms = new ArrayList<MetaPSM>(psmList);
		request.evidenceFile = evidenceFile;
		request.rankedFile = rankedFile;
		request.selectedFile = selectedFile;
		request.maximumSelected = MagGenomeSelectionService.DEFAULT_MAXIMUM_SELECTED;
		request.inferFromProteinNamesWhenEmpty = true;
		request.allowMarker4LibraryFallback = true;
		request.logDiagnostics = MagDbItem.MARKER4.equalsIgnoreCase(params.getGenomeSelectorMode());
		request.sourceDescription = psmTsv.getAbsolutePath();
		return MagGenomeSelectionService.selectGenomes(request, cb);
	}

	/**
	 * Append all sequences from {@code fastaFile} and create simple reverse decoys
	 * as >rev_ entries. (Matches the existing GUI behavior used for sample-specific
	 * DB construction.)
	 */
	private static void appendFastaWithDecoys(File fastaFile, PrintWriter writer) throws IOException {
		try (BufferedReader reader = new BufferedReader(new FileReader(fastaFile))) {
			String line;
			StringBuilder seq = null;
			String revHeader = null;

			while ((line = reader.readLine()) != null) {
				if (line.startsWith(">")) {
					if (revHeader != null && seq != null) {
						writer.println(revHeader);
						writer.println(seq.reverse());
					}
					seq = new StringBuilder();
					revHeader = ">rev_" + line.substring(1);
					writer.println(line);
				} else {
					if (seq != null)
						seq.append(line.trim());
					writer.println(line);
				}
			}
			if (revHeader != null && seq != null) {
				writer.println(revHeader);
				writer.println(seq.reverse());
			}
		}
	}

	private static String getCleavageResiduesOrDefault(FragpipeWorkflowFile wf) {
		try {
			String cut = wf.get("msfragger.search_enzyme_cut_1");
			if (cut != null && !cut.trim().isEmpty()) {
				return cut.trim();
			}
		} catch (Exception ignored) {
		}
		return "KR";
	}

	private static File resolveBasicWorkflowTemplate(File resourcesFragpipe, File explicitWorkflowTemplate)
			throws IOException {
		if (explicitWorkflowTemplate != null) {
			return explicitWorkflowTemplate;
		}
		return requireInstalledWorkflow(resourcesFragpipe, BASIC_SEARCH_WORKFLOW);
	}

	private static File resolveCombinedWorkflowTemplate(File resourcesFragpipe, File explicitWorkflowTemplate)
			throws IOException {
		if (explicitWorkflowTemplate != null) {
			return explicitWorkflowTemplate;
		}
		return requireInstalledWorkflow(resourcesFragpipe, LFQ_MBR_WORKFLOW);
	}

	public static File requireInstalledWorkflow(File resourcesFragpipe, String workflowName) throws IOException {
		File workflowsDir = resolveInstalledWorkflowsDir(resourcesFragpipe);
		File workflow = new File(workflowsDir, workflowName);
		if (!workflow.exists()) {
			throw new IOException("FragPipe workflow not found: " + workflow.getAbsolutePath());
		}
		if (!workflow.isFile()) {
			throw new IOException("FragPipe workflow is not a file: " + workflow.getAbsolutePath());
		}
		return workflow;
	}

	public static File resolveInstalledWorkflowsDir(File resourcesFragpipe) throws IOException {
		if (resourcesFragpipe == null) {
			throw new IOException("FragPipe executable path is null");
		}
		File fragpipeExe = resourcesFragpipe.getAbsoluteFile();
		File binDir = fragpipeExe.getParentFile();
		if (binDir == null) {
			throw new IOException("Cannot resolve FragPipe bin directory from: " + fragpipeExe.getAbsolutePath());
		}
		File fragpipeHome = binDir.getParentFile();
		if (fragpipeHome == null) {
			throw new IOException("Cannot resolve FragPipe home directory from: " + fragpipeExe.getAbsolutePath());
		}
		File workflowsDir = new File(fragpipeHome, "workflows");
		if (!workflowsDir.exists()) {
			throw new IOException("FragPipe workflows directory not found: " + workflowsDir.getAbsolutePath());
		}
		if (!workflowsDir.isDirectory()) {
			throw new IOException("FragPipe workflows path is not a directory: " + workflowsDir.getAbsolutePath());
		}
		return workflowsDir;
	}

	private static File resolveUserPath(String path) {
		if (path == null) {
			return null;
		}
		String trimmed = path.trim();
		if (trimmed.equals("~")) {
			String home = System.getProperty("user.home", "");
			return home == null || home.isBlank() ? new File(trimmed) : new File(home);
		}
		if (trimmed.startsWith("~/") || trimmed.startsWith("~\\")) {
			String home = System.getProperty("user.home", "");
			if (home == null || home.isBlank()) {
				return new File(trimmed);
			}
			return new File(home, trimmed.substring(2));
		}
		String home = System.getenv("HOME");
		if (home != null && !home.isBlank()) {
			if (trimmed.equals("$HOME")) {
				return new File(home);
			}
			if (trimmed.startsWith("$HOME/") || trimmed.startsWith("$HOME\\")) {
				return new File(home, trimmed.substring(6));
			}
		}
		return new File(trimmed);
	}

	private static void appendHostFastaWithDecoys(File hostFasta, PrintWriter writer, String cleavage)
			throws IOException {
		try (BufferedReader reader = new BufferedReader(new FileReader(hostFasta))) {
			String line;
			String header = null;
			StringBuilder seqSb = new StringBuilder();
			while ((line = reader.readLine()) != null) {
				line = line.trim();
				if (line.isEmpty())
					continue;
				if (line.startsWith(">")) {
					if (header != null) {
						String seq = seqSb.toString();
						writer.println(header);
						writer.println(seq);

						String shuffle = DecoyShuffler.generateShuffledDecoy(seq, cleavage, 42);
						writer.println(">shuf_" + header.substring(1));
						writer.println(shuffle);

						writer.println(">rev_" + header.substring(1));
						writer.println(new StringBuilder(seq).reverse().toString());

						writer.println(">rev_shuf_" + header.substring(1));
						writer.println(new StringBuilder(shuffle).reverse().toString());
					}
					header = line;
					seqSb.setLength(0);
				} else {
					seqSb.append(line);
				}
			}
			if (header != null) {
				String seq = seqSb.toString();
				writer.println(header);
				writer.println(seq);

				String shuffle = DecoyShuffler.generateShuffledDecoy(seq, cleavage, 42);
				writer.println(">shuf_" + header.substring(1));
				writer.println(shuffle);

				writer.println(">rev_" + header.substring(1));
				writer.println(new StringBuilder(seq).reverse().toString());

				writer.println(">rev_shuf_" + header.substring(1));
				writer.println(new StringBuilder(shuffle).reverse().toString());
			}
		}
	}

	private static LinkedHashMap<String, String> readFastaSequences(File fasta) throws IOException {
		LinkedHashMap<String, String> fastaMap = new LinkedHashMap<>();
		try (BufferedReader reader = new BufferedReader(new FileReader(fasta))) {
			String header = null;
			StringBuilder sequence = new StringBuilder();
			String line;
			while ((line = reader.readLine()) != null) {
				line = line.trim();
				if (line.isEmpty()) {
					continue;
				}
				if (line.startsWith(">")) {
					if (header != null) {
						fastaMap.put(header, sequence.toString());
					}
					int space = line.indexOf(' ');
					header = space > 0 ? line.substring(0, space) : line;
					sequence.setLength(0);
				} else if (header != null) {
					sequence.append(line);
				}
			}
			if (header != null) {
				fastaMap.put(header, sequence.toString());
			}
		}
		return fastaMap;
	}

	private static void exportTaxonomyTsvs(DDAMagFragPipeParameters params, File combinedModifiedPeptideTsv,
			MagDbItem mag, File resultRoot, FragpipePostAnalysisTsvExporter.Result post, ProgressCallback cb)
			throws Exception {
		if (cb == null)
			cb = new ProgressCallback.NoOp();
		if (mag == null)
			throw new IllegalArgumentException("mag is null");
		if (post == null)
			throw new IllegalArgumentException("post-analysis result is null");

		File taxDir = new File(resultRoot, "taxonomy_analysis");
		taxDir.mkdirs();

		File gtdbNcbi = new File("Resources", "gtdb_ncbi");
		if (!gtdbNcbi.exists()) {
			gtdbNcbi = null;
		}
		MagTaxonomyAnalysisService.Request request = new MagTaxonomyAnalysisService.Request();
		request.taxonomyAnalysisDir = taxDir;
		request.magDb = mag;
		request.gtdbNcbiMappingTsv = gtdbNcbi;
		request.generateFinalGenomeRanking = true;
		request.genomeIntensityMap = post.genomeIntensityMap;
		request.genomeScoreMap = post.genomeScoreMap;
		request.genomeStats = post.genomeDiaStats;
		request.psmCountMap = post.psmCountMap;
		request.totalPepCountMap = post.totalPepCountMap;
		request.razorPepCountMap = post.razorPepCountMap;
		request.uniquePepCountMap = post.uniquePepCountMap;
		request.totalPepWithILMap = post.totalPepWithILMap;
		request.razorPepWithILMap = post.razorPepWithILMap;
		request.uniquePepWithILMap = post.uniquePepWithILMap;
		request.intensityTitles = post.intensityTitles;
		request.progressCallback = cb;
		MagTaxonomyAnalysisService.run(request);
	}

}
