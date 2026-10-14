package uk.ac.quadram.metapilot.core.runner;

import uk.ac.quadram.metapilot.core.mag.DdaMagGenomeScoringService;
import uk.ac.quadram.metapilot.core.mag.GenomeEvidenceTest;
import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Shared taxonomy analysis for MAG workflows.
 * <p>
 * Use this service from GUI/task code that already has MAG genome-level abundance
 * maps and peptide support maps. It resolves the MAG taxonomy TSV, optionally
 * writes a fresh final genome ranking, carries over ranked-selection labels and
 * genome ordering, and delegates the actual TSV/HTML writing to
 * {@link MagTaxonomyTsvExporter}.
 * <p>
 * Current intended callers:
 * <ul>
 * <li>DDA MAG GUI tasks via {@code MetaPilotMagTask.performTaxonomyAnalysis}</li>
 * <li>DIA Marker/Seed GUI tasks via {@code MetaDiaTask.exportReport}</li>
 * <li>DDA MAG headless runners after post-analysis abundance is available</li>
 * </ul>
 * Workflows that do not produce MAG genome abundance maps, such as DIA Target
 * and standalone Seed Direct runs, should not call this service.
 */
public final class MagTaxonomyAnalysisService {

	private MagTaxonomyAnalysisService() {
	}

	public static final class Request {
		public File taxonomyAnalysisDir;
		public MagDbItem magDb;
		public File gtdbNcbiMappingTsv;
		/**
		 * Optional manual report-ranking TSV. Ignored when generateFinalGenomeRanking is true.
		 */
		public File rankedTsv;
		public boolean generateFinalGenomeRanking;
		public HashMap<String, double[]> genomeIntensityMap;
		public HashMap<String, Double> genomeScoreMap;
		public HashMap<String, double[]> genomeStats;
		public DdaMagGenomeScoringService.DdaGenomeScoreRuninfo genomeScoreRuninfo;
		/**
		 * The open-search genome test behind {@link #genomeStats}; when set, genome_scores.tsv and
		 * its run-info are written instead of genome_rank_scores.tsv.
		 */
		public GenomeEvidenceTest.Result genomeEvidence;
		public HashMap<String, Integer> psmCountMap;
		public HashMap<String, HashSet<String>> totalPepCountMap;
		public HashMap<String, HashSet<String>> razorPepCountMap;
		public HashMap<String, HashSet<String>> uniquePepCountMap;
		public HashMap<String, HashSet<String>> totalPepWithILMap;
		public HashMap<String, HashSet<String>> razorPepWithILMap;
		public HashMap<String, HashSet<String>> uniquePepWithILMap;
		public String[] intensityTitles;
		public ProgressCallback progressCallback;
		/**
		 * Opt-in: also write the MetaAtlas-aligned taxon_abundance.parquet (see
		 * {@link MagTaxonomyTsvExporter#exportMetaAtlasTaxonAbundance}). Defaults to false,
		 * so existing callers (in particular GUI tasks) get no behavior change unless they
		 * explicitly opt in. Intended to be set from a CLI flag (e.g. --export-metaatlas),
		 * not from GUI code.
		 */
		public boolean exportMetaAtlasParquet;
	}

	public static final class Result {
		public File taxonomyAnalysisDir;
		public File genomeHtml;
		public File genomeTsv;
		public File taxaTsv;
		public File taxaReportTsv;
		/** genome_rank_scores.tsv, or genome_scores.tsv for the open-search genome test. */
		public File genomeRankScoresTsv;
		public File genomeRankScoresRuninfo;
		public RankedTsvLabelParser.RankedGenomeSelectionResult rankedSelection;
	}

	public static Result run(Request request) throws IOException {
		if (request == null) throw new IllegalArgumentException("request is null");
		if (request.magDb == null) throw new IllegalArgumentException("magDb is null");
		if (request.taxonomyAnalysisDir == null) throw new IllegalArgumentException("taxonomyAnalysisDir is null");
		ProgressCallback cb = request.progressCallback != null ? request.progressCallback : new ProgressCallback.NoOp();
		request.taxonomyAnalysisDir.mkdirs();

		File rankedTsv = request.generateFinalGenomeRanking ? null : request.rankedTsv;
		if (request.generateFinalGenomeRanking) {
			MagFinalGenomeRankingService.Request rankingRequest = new MagFinalGenomeRankingService.Request();
			rankingRequest.taxonomyAnalysisDir = request.taxonomyAnalysisDir;
			rankingRequest.genomeIntensityMap = request.genomeIntensityMap;
			rankingRequest.genomePsmCountMap = request.psmCountMap;
			rankingRequest.totalPepCountMap = request.totalPepCountMap;
			rankingRequest.razorPepCountMap = request.razorPepCountMap;
			rankingRequest.uniquePepCountMap = request.uniquePepCountMap;
			MagFinalGenomeRankingService.Result finalRanking = MagFinalGenomeRankingService.write(rankingRequest);
			rankedTsv = finalRanking.rankedTsv;
			cb.onMessage("Final genome ranking written: " + rankedTsv.getAbsolutePath());
		}

		RankedTsvLabelParser.RankedGenomeSelectionResult rankedSelection = new RankedTsvLabelParser.RankedGenomeSelectionResult();
		if (rankedTsv != null && rankedTsv.isFile()) {
			if (request.genomeIntensityMap != null) {
				rankedSelection = RankedTsvLabelParser.buildExportSelectionFromRankedTsv(rankedTsv,
						request.genomeIntensityMap.keySet());
			}
		}

		MagTaxonomyTsvExporter.exportGenomeAndTaxaReports(request.taxonomyAnalysisDir,
				request.magDb.getOriginalTaxaFile(), request.gtdbNcbiMappingTsv, request.genomeIntensityMap,
				request.genomeScoreMap, request.genomeStats, rankedSelection.orderedGenomes, request.psmCountMap,
				rankedSelection.genomeSelectionLabel, request.totalPepCountMap, request.razorPepCountMap,
				request.uniquePepCountMap, request.totalPepWithILMap, request.razorPepWithILMap,
				request.uniquePepWithILMap, request.intensityTitles, request.genomeScoreRuninfo,
				request.genomeEvidence, cb);

		if (request.exportMetaAtlasParquet) {
			try {
				MagTaxonomyTsvExporter.exportMetaAtlasTaxonAbundance(request.taxonomyAnalysisDir,
						request.magDb.getOriginalTaxaFile(), request.genomeIntensityMap, request.genomeStats,
						request.genomeScoreMap, request.intensityTitles, cb);
			} catch (Exception e) {
				cb.onMessage("Skipping MetaAtlas taxon_abundance.parquet: " + e.getMessage());
			}
		}

		Result result = new Result();
		result.taxonomyAnalysisDir = request.taxonomyAnalysisDir;
		result.genomeHtml = new File(request.taxonomyAnalysisDir, "Genome.html");
		result.genomeTsv = new File(request.taxonomyAnalysisDir, "Genome.tsv");
		result.taxaTsv = new File(request.taxonomyAnalysisDir, "Taxa.tsv");
		result.taxaReportTsv = new File(request.taxonomyAnalysisDir, "Taxa_report.tsv");
		boolean evidence = request.genomeEvidence != null;
		result.genomeRankScoresTsv = new File(request.taxonomyAnalysisDir,
				evidence ? MagTaxonomyTsvExporter.GENOME_SCORES_TSV : "genome_rank_scores.tsv");
		result.genomeRankScoresRuninfo = new File(request.taxonomyAnalysisDir,
				evidence ? MagTaxonomyTsvExporter.GENOME_SCORES_RUNINFO : "genome_rank_scores_runinfo.txt");
		result.rankedSelection = rankedSelection;
		return result;
	}
}
