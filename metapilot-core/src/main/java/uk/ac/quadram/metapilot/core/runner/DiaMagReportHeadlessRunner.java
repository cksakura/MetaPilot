package uk.ac.quadram.metapilot.core.runner;

import uk.ac.quadram.metapilot.core.dia.DiaDeepProteinMappingService;
import uk.ac.quadram.metapilot.core.report.MetaPilotHtmlReportService;
import uk.ac.quadram.metapilot.core.resources.MetaPilotResourceConfig;
import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;
import uk.ac.quadram.metapilot.task.dia.DiaRefinementMode;
import uk.ac.quadram.metapilot.task.dia.par.MetaParameterDiaMarker;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Headless "resume from DIA-NN search output" runner for the DIA MAG (DIA
 * Marker) workflow.
 *
 * The DIA-NN search itself (Marker4 no-quant selection pass + combined-library
 * quant pass) still requires DIA-NN and normally runs on the user's own
 * machine, via the Windows GUI or {@link DiaMarkerHeadlessRunner}. This class
 * takes an already-produced {@code combined.parquet} +
 * {@code combined.pr_matrix.tsv} - the same files {@link DiaMarkerHeadlessRunner}
 * writes under {@code <result>/marker_focus|marker_comprehensive/combined/}
 * (or {@code <result>/marker_deep/combined/} for Deep Together and Deep Per
 * File) - and runs only the memory-heavy remainder of the pipeline headlessly: post-analysis
 * TSV export (protein/genome inference + scoring), taxonomy analysis,
 * functional annotation, and the final HTML report.
 *
 * Unlike {@link DiaMarkerHeadlessRunner}, no DIA-NN executable, Marker4/combined
 * FASTA or library build, or raw DIA spectra files are required or read here -
 * only {@code combined.parquet}, {@code combined.pr_matrix.tsv} (both already
 * self-describing: run/sample names are read from the pr_matrix column
 * headers), and the MAG database. This mirrors {@link PfindMagHeadlessRunner}
 * for the DDA pFind MAG workflow.
 *
 * The {@code result} field inside the DIA Marker params JSON is exported by
 * the Windows GUI and is a local path on that machine, unrelated to where this
 * runs. When {@code resultDirOverride} is supplied, it is used verbatim as the
 * output directory and {@code params.getResult()} is never read - so nothing
 * in {@code params}, {@code magDbConfigJson} (once it is a fresh HPC-local
 * registry), or the two combined result files needs to reference the local
 * machine at all.
 *
 * Deep Together / Deep Per File caveat - THIS ONE IS NOT OPTIONAL: for those two
 * marker modes, {@code combined.fasta} (and therefore {@code combined.speclib}
 * and {@code combined.parquet}) is built with exactly one FASTA record per
 * unique peptide across the whole run ({@code DeepIterativeCommonSupport
 * .writeCombinedFastaFromPeptideProteinMap}), tagged with a single arbitrarily
 * chosen protein ID out of a {@code HashSet} ({@code pros.iterator().next()}).
 * DIA-NN therefore cannot report more than that one protein per peptide in
 * {@code combined.parquet} - a peptide genuinely shared across several genomes
 * shows only whichever genome happened to win iteration order. Razor-protein
 * assignment, unique-vs-shared peptide classification, and genome confidence
 * scoring all depend on knowing every genome a peptide matches, so scoring
 * directly off {@code combined.parquet}'s own protein groups for these two
 * modes would silently produce wrong genome/protein inference, not merely a
 * less-complete one.
 *
 * The GUI/{@code dia-marker} pipeline avoids this by re-mapping each peptide
 * back to every protein it matched across the iterative marker rounds (see
 * {@link DiaDeepProteinMappingService}) using the {@code round1.parquet,
 * round2.parquet, ...} files under {@code <result>/marker} (Deep Together) or
 * {@code <result>/marker/<rawname>/round*} (Deep Per File) - those round
 * parquets are built from real full-length genome protein sequences and do
 * carry the complete peptide-to-protein membership.
 *
 * Rather than requiring that whole round-parquet folder tree to be copied to
 * HPC, the recommended path is to run {@link DiaMagPeptideProteinMapExportRunner}
 * once, locally, where the round parquets already sit; it reads them and writes
 * one compact peptide-to-proteins TSV. Only that small file needs to travel to
 * HPC ({@code peptideProteinMapOverride} here). The round-parquet folder itself
 * is still accepted as a fallback ({@code roundRootOverride}) for anyone who'd
 * rather copy it directly. For Deep modes, one of the two is required unless
 * {@code allowIncompleteDeepProteinMapping} is explicitly set - this runner
 * fails loudly rather than silently scoring off the reduced
 * single-protein-per-peptide combined.parquet.
 */
public final class DiaMagReportHeadlessRunner {

	private DiaMagReportHeadlessRunner() {
	}

	public static final class Result {
		public final File resultDir;
		public final File combinedParquet;
		public final File prMatrixTsv;
		public final DiaMagPostAnalysisTsvExporter.Result post;

		Result(File resultDir, File combinedParquet, File prMatrixTsv, DiaMagPostAnalysisTsvExporter.Result post) {
			this.resultDir = resultDir;
			this.combinedParquet = combinedParquet;
			this.prMatrixTsv = prMatrixTsv;
			this.post = post;
		}
	}

	public static Result run(MetaParameterDiaMarker params, File magDbConfigJson, File resultDirOverride,
			File combinedParquetOverride, File combinedPrMatrixOverride, File peptideProteinMapOverride,
			File roundRootOverride, File funcDbFile, ProgressCallback cb) throws Exception {
		return run(params, magDbConfigJson, resultDirOverride, combinedParquetOverride, combinedPrMatrixOverride,
				peptideProteinMapOverride, roundRootOverride, funcDbFile, cb, true, null, null, false, false);
	}

	/**
	 * @param resultDirOverride          override for the output directory; when null,
	 *                                    falls back to deriving it from
	 *                                    {@code params.getResult()} (a path on the
	 *                                    machine the GUI ran on - only usable when
	 *                                    this runs on that same machine/filesystem).
	 *                                    When set, this exact directory is used (no
	 *                                    marker_focus/marker_comprehensive/marker_deep
	 *                                    subfolder is appended) and params.getResult()
	 *                                    is never read.
	 * @param combinedParquetOverride    override for combined.parquet; when null,
	 *                                    defaults to {@code <result-dir>/combined/combined.parquet}
	 *                                    (the same relative location the Windows
	 *                                    GUI / {@code dia-marker} write it to).
	 * @param combinedPrMatrixOverride    override for combined.pr_matrix.tsv; when
	 *                                    null, same default rule.
	 * @param peptideProteinMapOverride   Deep Together / Deep Per File only, and the
	 *                                    recommended input for them: a compact TSV
	 *                                    produced once, locally, by
	 *                                    {@link DiaMagPeptideProteinMapExportRunner}.
	 *                                    When set, this is used directly and
	 *                                    {@code roundRootOverride} / the round-parquet
	 *                                    folder are never read. Takes priority over
	 *                                    {@code roundRootOverride} when both are set.
	 * @param roundRootOverride          Deep Together / Deep Per File only: fallback
	 *                                    override for the {@code marker} round-parquet
	 *                                    folder itself (see class docs) when you'd
	 *                                    rather copy it directly instead of using
	 *                                    {@code peptideProteinMapOverride}. When null,
	 *                                    defaults to {@code <result-dir>/marker}.
	 *                                    Ignored for Core/Focus/Comprehensive.
	 * @param exportMetaAtlasParquet      opt-in MetaAtlas-aligned parquet export
	 *                                    (see {@code --export-metaatlas} in the
	 *                                    CLI).
	 * @param allowIncompleteDeepProteinMapping  Deep Together / Deep Per File only: when true,
	 *                                    proceeds even if neither peptideProteinMapOverride
	 *                                    nor the round-parquet folder is available, scoring
	 *                                    directly off combined.parquet's own protein groups.
	 *                                    Per class docs, that collapses every peptide to a
	 *                                    single arbitrary protein and will misclassify
	 *                                    peptides shared across genomes - only for debugging /
	 *                                    when you understand and accept that limitation.
	 */
	public static Result run(MetaParameterDiaMarker params, File magDbConfigJson, File resultDirOverride,
			File combinedParquetOverride, File combinedPrMatrixOverride, File peptideProteinMapOverride,
			File roundRootOverride, File funcDbFile, ProgressCallback cb, boolean generateHtmlReport,
			String reportPythonHint, MetaPilotResourceConfig reportResourceConfig, boolean exportMetaAtlasParquet,
			boolean allowIncompleteDeepProteinMapping) throws Exception {
		if (cb == null)
			cb = new ProgressCallback.NoOp();
		if (params == null)
			throw new IllegalArgumentException("params is null");

		DiaRefinementMode markerMode = params.getMarkerRefinementMode();
		if (markerMode == DiaRefinementMode.DIRECT) {
			throw new IOException("diann-mag-report does not support DIRECT mode. DIRECT is only valid for DIA Seed workflows.");
		}

		MagDbItem magDb = DiaMarkerHeadlessRunner.resolveMagDb(params, magDbConfigJson);
		if (magDb == null) {
			throw new IOException("MAG database not found in " + magDbConfigJson + " for micro_db="
					+ params.getMicroDb() + " identifier=" + params.getMagDbIdentifier());
		}
		cb.onMessage("Using MAG DB: " + magDb.getCatalogueID()
				+ (magDb.getUsedVersion() != null ? " / " + magDb.getUsedVersion().getCatalogueVersion() : ""));

		File resultDir;
		File defaultCombinedDir;
		if (resultDirOverride != null) {
			resultDir = resultDirOverride;
			defaultCombinedDir = new File(resultDir, "combined");
		} else if (markerMode == DiaRefinementMode.DEEP_TOGETHER || markerMode == DiaRefinementMode.DEEP_PER_FILE) {
			File studyRoot = new File(params.getResult());
			File markerDeepDir = new File(studyRoot, "marker_deep");
			resultDir = markerDeepDir;
			defaultCombinedDir = new File(markerDeepDir, "combined");
		} else {
			File studyRoot = new File(params.getResult());
			String workFolderName = markerMode.usesComprehensiveSelection() ? "marker_comprehensive" : "marker_focus";
			resultDir = new File(studyRoot, workFolderName);
			defaultCombinedDir = new File(resultDir, "combined");
		}
		if (!resultDir.exists())
			resultDir.mkdirs();
		cb.onMessage("DIA MAG report: result folder -> " + resultDir.getAbsolutePath());

		File combinedParquet = combinedParquetOverride != null ? combinedParquetOverride
				: new File(defaultCombinedDir, "combined.parquet");
		File prMatrix = combinedPrMatrixOverride != null ? combinedPrMatrixOverride
				: new File(defaultCombinedDir, "combined.pr_matrix.tsv");

		if (!combinedParquet.exists()) {
			throw new IOException("combined.parquet not found: " + combinedParquet.getAbsolutePath());
		}
		if (!prMatrix.exists()) {
			throw new IOException("combined.pr_matrix.tsv not found: " + prMatrix.getAbsolutePath());
		}
		cb.onMessage("Resuming from existing DIA-NN search output: " + combinedParquet.getAbsolutePath() + ", "
				+ prMatrix.getAbsolutePath());

		File gtdbNcbi = null;
		File maybeGtdbNcbi = new File("Resources", "gtdb_ncbi");
		if (maybeGtdbNcbi.exists())
			gtdbNcbi = maybeGtdbNcbi;

		DiaDeepProteinMappingService.Request deepProteinMappingRequest = buildDeepProteinMappingRequest(markerMode,
				resultDir, peptideProteinMapOverride, roundRootOverride, params, cb, allowIncompleteDeepProteinMapping);

		cb.onMessage("Step 1: post-analysis (protein/genome inference) + TSV export ...");
		DiaMagPostAnalysisTsvExporter.Result post = DiaMagPostAnalysisTsvExporter.export(resultDir, magDb,
				combinedParquet, prMatrix, gtdbNcbi, params.getMarkerMode(), funcDbFile, params.getMicroDb(),
				params.getQuantitativeMode(), cb, deepProteinMappingRequest, exportMetaAtlasParquet);

		runHtmlReportIfRequested(resultDir, params, generateHtmlReport, reportPythonHint, reportResourceConfig, cb);

		cb.onMessage("DIA MAG headless post-analysis completed.");
		cb.onComplete(resultDir.getAbsolutePath());
		return new Result(resultDir, combinedParquet, prMatrix, post);
	}

	/**
	 * Resolution order for Deep Together / Deep Per File: (1) peptideProteinMapOverride, if
	 * set - read directly, no round parquet files touched; (2) the round-parquet folder
	 * (roundRootOverride or its default location under resultDir); (3) neither found - hard
	 * failure unless allowIncompleteDeepProteinMapping.
	 */
	private static DiaDeepProteinMappingService.Request buildDeepProteinMappingRequest(DiaRefinementMode markerMode,
			File resultDir, File peptideProteinMapOverride, File roundRootOverride, MetaParameterDiaMarker params,
			ProgressCallback cb, boolean allowIncompleteDeepProteinMapping) throws IOException {
		if (markerMode != DiaRefinementMode.DEEP_TOGETHER && markerMode != DiaRefinementMode.DEEP_PER_FILE) {
			return null;
		}

		if (peptideProteinMapOverride != null) {
			cb.onMessage("DIA MAG report: reading precomputed peptide-protein map from "
					+ peptideProteinMapOverride.getAbsolutePath());
			HashMap<String, HashSet<String>> map = DiaDeepProteinMappingService.readMapTsv(peptideProteinMapOverride);
			cb.onMessage("DIA MAG report: loaded " + map.size() + " peptide-protein mappings");
			return DiaDeepProteinMappingService.Request.prebuiltMap(map);
		}

		String[] diaFiles = params.getMetadata() != null ? params.getMetadata().getRawFiles() : null;
		if (markerMode == DiaRefinementMode.DEEP_TOGETHER) {
			File roundRoot = roundRootOverride != null ? roundRootOverride : new File(resultDir, "marker");
			if (!roundRoot.isDirectory()) {
				failOrWarnMissingRoundRoot(roundRoot, allowIncompleteDeepProteinMapping, cb);
				return null;
			}
			cb.onMessage("DIA MAG report: remapping deep peptides using round parquets under "
					+ roundRoot.getAbsolutePath());
			return DiaDeepProteinMappingService.Request.deepTogether(roundRoot, 0, diaFiles, "KR");
		}
		// DEEP_PER_FILE: Request.deepPerFile expects the parent of "marker" (i.e. marker_deep/), not
		// marker/ itself - see DiaDeepProteinMappingService.buildPeptideProteinMap perFileLayout branch.
		File roundRoot = roundRootOverride != null ? roundRootOverride : resultDir;
		File markerRoot = new File(roundRoot, "marker");
		if (!markerRoot.isDirectory()) {
			failOrWarnMissingRoundRoot(markerRoot, allowIncompleteDeepProteinMapping, cb);
			return null;
		}
		cb.onMessage("DIA MAG report: remapping deep peptides using per-file round parquets under "
				+ markerRoot.getAbsolutePath());
		return DiaDeepProteinMappingService.Request.deepPerFile(roundRoot, diaFiles, "KR");
	}

	private static void failOrWarnMissingRoundRoot(File missingDir, boolean allowIncompleteDeepProteinMapping,
			ProgressCallback cb) throws IOException {
		String message = "no peptide-protein mapping input found for this Deep-mode run: neither "
				+ "--peptide-protein-map was given, nor does a round-parquet folder exist at "
				+ missingDir.getAbsolutePath() + ". For Deep Together / Deep Per File, combined.parquet alone tags "
				+ "every peptide with a single arbitrary protein (see class docs) - scoring off it directly would "
				+ "silently misclassify peptides shared across genomes. Run diann-mag-peptide-map locally (where "
				+ "the round parquets already exist) and pass its output via --peptide-protein-map, or copy the "
				+ "'marker' folder and use --round-root, or pass --allow-incomplete-deep-protein-mapping to proceed "
				+ "anyway and accept that limitation.";
		if (!allowIncompleteDeepProteinMapping) {
			throw new IOException("diann-mag-report: " + message);
		}
		cb.onMessage("WARNING: diann-mag-report: " + message + " Proceeding anyway because "
				+ "--allow-incomplete-deep-protein-mapping was set - genome/protein inference for this run will be "
				+ "less accurate than the original local result.");
	}

	private static void runHtmlReportIfRequested(File resultDir, MetaParameterDiaMarker params,
			boolean generateHtmlReport, String reportPythonHint, MetaPilotResourceConfig reportResourceConfig,
			ProgressCallback cb) {
		if (!generateHtmlReport) {
			cb.onMessage("HTML report generation skipped (disabled)");
			return;
		}
		cb.onMessage("DIA MAG report: generating shared HTML report");
		MetaPilotHtmlReportService.Result report = MetaPilotHtmlReportService.runStandardReportWarningOnly(resultDir,
				"DIA_MAG", params != null ? params.getThreadCount() : 1, reportPythonHint, reportResourceConfig, cb);
		if (!report.warnings.isEmpty()) {
			cb.onMessage(
					"WARNING: HTML report generation completed with warnings: " + String.join("; ", report.warnings));
		}
	}
}
