package uk.ac.quadram.metapilot.core.runner;

import uk.ac.quadram.metapilot.core.resources.MetaPilotResourceConfig;
import uk.ac.quadram.metapilot.core.report.MetaPilotHtmlReportService;
import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;
import uk.ac.quadram.metapilot.task.dda.mag.par.DDAMagPfindParameters;
import uk.ac.quadram.metapilot.task.mag.MagDbConfigIO;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;
import uk.ac.quadram.metapilot.task.par.IsobaricLabelSettings;
import uk.ac.quadram.metapilot.task.par.MetaConstants;

import java.io.File;
import java.io.IOException;

/**
 * Headless "resume from quantification output" runner for the DDA MAG pFind
 * workflow.
 *
 * Takes an already-produced peptide quantification + {@code psms.tsv} (the
 * files {@code DDAMagPfindTask} / {@link PfindMagSearchHeadlessRunner} write
 * under {@code <result>/pfind_MAG_search/combined/}) and runs the rest of the
 * pipeline headlessly: post-analysis TSV export, taxonomy analysis,
 * functional annotation, and the final HTML report - reusing exactly the
 * same core services the GUI pFind path and the FragPipe/DIA headless
 * runners already share. {@link PfindMagSearchHeadlessRunner} runs the
 * identification and quantification steps before this and then calls
 * {@link #postAnalysis}.
 *
 * The peptide table follows the quantitative mode in the parameters:
 * FlashLFQ's {@code QuantifiedPeptides.tsv} for label-free,
 * {@code IsobaricQuantifiedPeptides.tsv} (written by
 * {@code PFindIsobaricTask}) for isobaric labelling. Only that file and
 * {@code psms.tsv} are required to exist; no other intermediate result file
 * (marker/sample_db search output, ExperimentalDesign.tsv, etc.) is checked.
 */
public final class PfindMagHeadlessRunner {

	private PfindMagHeadlessRunner() {
	}

	public static final class Result {
		public final File resultRoot;
		public final File quantifiedPeptidesTsv;
		public final File psmsTsv;
		public final PfindMagPostAnalysisTsvExporter.Result post;

		Result(File resultRoot, File quantifiedPeptidesTsv, File psmsTsv,
				PfindMagPostAnalysisTsvExporter.Result post) {
			this.resultRoot = resultRoot;
			this.quantifiedPeptidesTsv = quantifiedPeptidesTsv;
			this.psmsTsv = psmsTsv;
			this.post = post;
		}
	}

	/**
	 * @param quantifiedPeptidesTsv override for QuantifiedPeptides.tsv; when null,
	 *                              defaults to
	 *                              {@code <result>/pfind_MAG_search/combined/QuantifiedPeptides.tsv}
	 *                              (the same location the Windows GUI writes it to).
	 * @param psmsTsv               override for psms.tsv; when null, defaults to
	 *                              {@code <result>/pfind_MAG_search/combined/psms.tsv}.
	 */
	public static Result run(DDAMagPfindParameters params, File magdbConfigJson, File quantifiedPeptidesTsv,
			File psmsTsv, File funcDbFile, ProgressCallback cb) throws Exception {
		return run(params, magdbConfigJson, quantifiedPeptidesTsv, psmsTsv, funcDbFile, cb, true, null, null);
	}

	public static Result run(DDAMagPfindParameters params, File magdbConfigJson, File quantifiedPeptidesTsv,
			File psmsTsv, File funcDbFile, ProgressCallback cb, boolean generateHtmlReport, String reportPythonHint,
			MetaPilotResourceConfig reportResourceConfig) throws Exception {
		return run(params, magdbConfigJson, quantifiedPeptidesTsv, psmsTsv, funcDbFile, cb, generateHtmlReport,
				reportPythonHint, reportResourceConfig, false);
	}

	/**
	 * Same as the 9-argument {@link #run}, with an additional opt-in
	 * {@code exportMetaAtlasParquet} flag (see {@code --export-metaatlas} in the CLI). This
	 * is a brand new overload, not a change to any existing signature, so it cannot affect
	 * any existing caller. Only CLI headless callers that explicitly opt in by calling this
	 * overload will produce the guide-aligned taxon_abundance.parquet /
	 * function_abundance.parquet exports.
	 */
	public static Result run(DDAMagPfindParameters params, File magdbConfigJson, File quantifiedPeptidesTsv,
			File psmsTsv, File funcDbFile, ProgressCallback cb, boolean generateHtmlReport, String reportPythonHint,
			MetaPilotResourceConfig reportResourceConfig, boolean exportMetaAtlasParquet) throws Exception {
		if (cb == null)
			cb = new ProgressCallback.NoOp();
		if (params == null)
			throw new IllegalArgumentException("params is null");

		MagDbItem mag = resolveMagDb(params.getMagDbCatalogueId(), params.getMagDbVersion(), magdbConfigJson);
		if (mag == null) {
			throw new IllegalArgumentException("MAG database not found in " + magdbConfigJson + " for catalogue="
					+ params.getMagDbCatalogueId() + " version=" + params.getMagDbVersion());
		}
		cb.onMessage("Using MAG DB: " + mag.getCatalogueID() + " / " + mag.getUsedVersion().getCatalogueVersion());

		File resultRoot = new File(params.getResult(), "pfind_MAG_search");
		resultRoot.mkdirs();
		cb.onMessage("Result root: " + resultRoot.getAbsolutePath());

		File combinedDir = new File(resultRoot, "combined");
		File quanPepFile = quantifiedPeptidesTsv != null ? quantifiedPeptidesTsv
				: defaultQuantifiedPeptides(combinedDir, params);
		File psmsFile = psmsTsv != null ? psmsTsv : new File(combinedDir, "psms.tsv");

		cb.onMessage("Resuming from existing quantification: " + quanPepFile.getAbsolutePath() + ", "
				+ psmsFile.getAbsolutePath());

		PfindMagPostAnalysisTsvExporter.Result post = postAnalysis(params, mag, resultRoot, quanPepFile, psmsFile,
				funcDbFile, cb, generateHtmlReport, reportPythonHint, reportResourceConfig, exportMetaAtlasParquet);

		cb.onMessage("MAG+pFind headless post-analysis completed.");
		cb.onComplete(resultRoot.getAbsolutePath());
		return new Result(resultRoot, quanPepFile, psmsFile, post);
	}

	/**
	 * The peptide table the parameters' quantitative mode produces in
	 * {@code combinedDir}: {@code IsobaricQuantifiedPeptides.tsv} for isobaric
	 * labelling, otherwise FlashLFQ's {@code QuantifiedPeptides.tsv}.
	 */
	static File defaultQuantifiedPeptides(File combinedDir, DDAMagPfindParameters params) {
		return isIsobaric(params) ? new File(combinedDir, "IsobaricQuantifiedPeptides.tsv")
				: new File(combinedDir, "QuantifiedPeptides.tsv");
	}

	static boolean isIsobaric(DDAMagPfindParameters params) {
		return MetaConstants.isobaricLabel
				.equals(IsobaricLabelSettings.normalizeQuantitativeMode(params.getQuantitativeMode()));
	}

	/**
	 * Post-analysis TSV export, taxonomy, functional annotation and the HTML
	 * report, from a peptide quantification and {@code psms.tsv}. Shared by
	 * {@link #run} and {@link PfindMagSearchHeadlessRunner}.
	 */
	static PfindMagPostAnalysisTsvExporter.Result postAnalysis(DDAMagPfindParameters params, MagDbItem mag,
			File resultRoot, File quanPepFile, File psmsFile, File funcDbFile, ProgressCallback cb,
			boolean generateHtmlReport, String reportPythonHint, MetaPilotResourceConfig reportResourceConfig,
			boolean exportMetaAtlasParquet) throws Exception {
		String quanMode = IsobaricLabelSettings.normalizeQuantitativeMode(params.getQuantitativeMode());
		if (quanMode == null || quanMode.isBlank()) {
			quanMode = MetaConstants.labelFree;
		}
		cb.onMessage("Step 1: post-analysis + TSV export (" + quanMode + ") ...");
		PfindMagPostAnalysisTsvExporter.Result post = PfindMagPostAnalysisTsvExporter.exportFinalTsvs(mag,
				quanPepFile, psmsFile, resultRoot, params.getMetadata(), "DDA MAG pFind", cb, quanMode,
				params.isCombineLabel(), new File(params.getResult(), "spectra"));

		cb.onMessage("Step 2: taxonomy report export ...");
		exportTaxonomyTsvs(mag, resultRoot, post, cb, exportMetaAtlasParquet);

		if (funcDbFile != null && funcDbFile.isFile()) {
			cb.onMessage("Step 3: functional annotation report export ...");
			MagFunctionalAnnotationService.Request request = new MagFunctionalAnnotationService.Request();
			request.resultDir = resultRoot;
			request.magDb = mag;
			request.microDbCatalogueId = params.getMagDbCatalogueId();
			request.functionDefinitionDb = funcDbFile;
			request.analysisType = MetaConstants.pFindMag;
			request.quantitativeMode = quanMode;
			request.intensityTitles = post.intensityTitles;
			request.metaProteins = post.metaProteins;
			request.exportMetaAtlasParquet = exportMetaAtlasParquet;
			MagFunctionalAnnotationService.run(request);
		} else {
			cb.onMessage("pFind functional annotation: skipped (func_def.db not configured)");
		}

		runHtmlReportIfRequested(resultRoot, params, generateHtmlReport, reportPythonHint, reportResourceConfig, cb);
		return post;
	}

	private static void exportTaxonomyTsvs(MagDbItem mag, File resultRoot, PfindMagPostAnalysisTsvExporter.Result post,
			ProgressCallback cb, boolean exportMetaAtlasParquet) throws IOException {
		File taxDir = new File(resultRoot, "taxonomy_analysis");
		taxDir.mkdirs();

		MagTaxonomyAnalysisService.Request request = new MagTaxonomyAnalysisService.Request();
		request.taxonomyAnalysisDir = taxDir;
		request.magDb = mag;
		request.generateFinalGenomeRanking = true;
		request.genomeIntensityMap = post.genomeIntensityMap;
		request.genomeUniqueIntensityMap = post.genomeUniqueIntensityMap;
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
		request.exportMetaAtlasParquet = exportMetaAtlasParquet;
		MagTaxonomyAnalysisService.run(request);
	}

	private static void runHtmlReportIfRequested(File resultDir, DDAMagPfindParameters params,
			boolean generateHtmlReport, String reportPythonHint, MetaPilotResourceConfig reportResourceConfig,
			ProgressCallback cb) {
		if (!generateHtmlReport) {
			cb.onMessage("HTML report generation skipped (disabled)");
			return;
		}
		cb.onMessage("MAG+pFind: generating shared HTML report");
		MetaPilotHtmlReportService.Result report = MetaPilotHtmlReportService.runStandardReportWarningOnly(resultDir,
				"DDA_MAG_PFIND", params != null ? params.getThreadCount() : 1, reportPythonHint, reportResourceConfig,
				cb, params != null ? params.getMetadata() : null);
		if (!report.warnings.isEmpty()) {
			cb.onMessage(
					"WARNING: HTML report generation completed with warnings: " + String.join("; ", report.warnings));
		}
	}

	static MagDbItem resolveMagDb(String catalogueId, String version, File magdbConfigJson) {
		if (catalogueId == null || catalogueId.isBlank())
			return null;
		if (magdbConfigJson == null || !magdbConfigJson.exists())
			return null;

		MagDbItem[][] all = MagDbConfigIO.parse(magdbConfigJson);
		if (all == null)
			return null;
		for (MagDbItem[] group : all) {
			if (group == null)
				continue;
			for (MagDbItem item : group) {
				if (item == null)
					continue;
				if (!catalogueId.equals(item.getCatalogueID()))
					continue;
				if (version != null && !version.isBlank()) {
					for (MagDbItem.MagVersion mv : item.getVersionList()) {
						if (mv != null && version.equals(mv.getCatalogueVersion())) {
							item.setUsedVersion(mv);
							return item;
						}
					}
				}
				MagDbItem.MagVersion used = item.getUsedVersion();
				if (used != null) {
					item.setUsedVersion(used);
					return item;
				}
			}
		}
		return null;
	}
}
