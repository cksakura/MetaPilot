package uk.ac.quadram.metapilot.core.runner;

import uk.ac.quadram.metapilot.core.resources.MetaPilotResourceConfig;
import uk.ac.quadram.metapilot.core.report.MetaPilotHtmlReportService;
import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;
import uk.ac.quadram.metapilot.task.dda.mag.par.DDAMagPfindParameters;
import uk.ac.quadram.metapilot.task.mag.MagDbConfigIO;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;
import uk.ac.quadram.metapilot.task.par.MetaConstants;

import java.io.File;
import java.io.IOException;

/**
 * Headless "resume from FlashLFQ output" runner for the DDA MAG pFind
 * workflow.
 *
 * pFind identification and FlashLFQ quantification themselves still require
 * the Windows GUI ({@code DDAMagPfindTask}) - there is no headless
 * replacement for those steps here. What this class does is take an
 * already-produced {@code QuantifiedPeptides.tsv} + {@code psms.tsv} (the
 * same files {@code DDAMagPfindTask.lfQuant()}/{@code combineResults()}
 * write under {@code <result>/pfind_MAG_search/combined/}) and run the rest
 * of the pipeline headlessly: post-analysis TSV export, taxonomy analysis,
 * functional annotation, and the final HTML report - reusing exactly the
 * same core services the GUI pFind path and the FragPipe/DIA headless
 * runners already share.
 *
 * Only {@code QuantifiedPeptides.tsv} and {@code psms.tsv} are required to
 * exist; no other intermediate result file (marker/sample_db search output,
 * ExperimentalDesign.tsv, etc.) is checked.
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
				: new File(combinedDir, "QuantifiedPeptides.tsv");
		File psmsFile = psmsTsv != null ? psmsTsv : new File(combinedDir, "psms.tsv");

		cb.onMessage("Resuming from existing quantification: " + quanPepFile.getAbsolutePath() + ", "
				+ psmsFile.getAbsolutePath());

		cb.onMessage("Step 1: post-analysis + TSV export ...");
		PfindMagPostAnalysisTsvExporter.Result post = PfindMagPostAnalysisTsvExporter.exportFinalTsvs(mag,
				quanPepFile, psmsFile, resultRoot, params.getMetadata(), "DDA MAG pFind", cb);

		cb.onMessage("Step 2: taxonomy report export ...");
		exportTaxonomyTsvs(mag, resultRoot, post, cb);

		if (funcDbFile != null && funcDbFile.isFile()) {
			cb.onMessage("Step 3: functional annotation report export ...");
			MagFunctionalAnnotationService.Request request = new MagFunctionalAnnotationService.Request();
			request.resultDir = resultRoot;
			request.magDb = mag;
			request.microDbCatalogueId = params.getMagDbCatalogueId();
			request.functionDefinitionDb = funcDbFile;
			request.analysisType = MetaConstants.pFindMag;
			request.quantitativeMode = MetaConstants.labelFree;
			request.intensityTitles = post.intensityTitles;
			request.metaProteins = post.metaProteins;
			MagFunctionalAnnotationService.run(request);
		} else {
			cb.onMessage("pFind functional annotation: skipped (func_def.db not configured)");
		}

		runHtmlReportIfRequested(resultRoot, params, generateHtmlReport, reportPythonHint, reportResourceConfig, cb);

		cb.onMessage("MAG+pFind headless post-analysis completed.");
		cb.onComplete(resultRoot.getAbsolutePath());
		return new Result(resultRoot, quanPepFile, psmsFile, post);
	}

	private static void exportTaxonomyTsvs(MagDbItem mag, File resultRoot, PfindMagPostAnalysisTsvExporter.Result post,
			ProgressCallback cb) throws IOException {
		File taxDir = new File(resultRoot, "taxonomy_analysis");
		taxDir.mkdirs();

		MagTaxonomyAnalysisService.Request request = new MagTaxonomyAnalysisService.Request();
		request.taxonomyAnalysisDir = taxDir;
		request.magDb = mag;
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
				cb);
		if (!report.warnings.isEmpty()) {
			cb.onMessage(
					"WARNING: HTML report generation completed with warnings: " + String.join("; ", report.warnings));
		}
	}

	private static MagDbItem resolveMagDb(String catalogueId, String version, File magdbConfigJson) {
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
