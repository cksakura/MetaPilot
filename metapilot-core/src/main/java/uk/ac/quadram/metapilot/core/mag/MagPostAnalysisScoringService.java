package uk.ac.quadram.metapilot.core.mag;

import uk.ac.quadram.metapilot.core.dia.DiaMagScoringService;
import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNPrecursor;
import uk.ac.quadram.metapilot.task.io.pep.MetaPSM;
import uk.ac.quadram.metapilot.task.io.pep.MetaPeptide;

import java.util.HashMap;
import java.util.Set;

/**
 * Shared orchestration point for final MAG scoring and protein inference.
 */
public final class MagPostAnalysisScoringService {

	private MagPostAnalysisScoringService() {
	}

	public static final class DdaResult {
		public final DdaMagGenomeScoringService.Result genomeScoring;
		public final HashMap<String, Double> confidentGenomeScoreMap;
		public final DdaMagProteinInferenceService.Result proteinInference;
		public final HashMap<String, double[]> proteinScoreMap;

		private DdaResult(DdaMagGenomeScoringService.Result genomeScoring,
				HashMap<String, Double> confidentGenomeScoreMap,
				DdaMagProteinInferenceService.Result proteinInference,
				HashMap<String, double[]> proteinScoreMap) {
			this.genomeScoring = genomeScoring;
			this.confidentGenomeScoreMap = confidentGenomeScoreMap;
			this.proteinInference = proteinInference;
			this.proteinScoreMap = proteinScoreMap;
		}
	}

	public static final class DiaResult {
		public final DiaNNPrecursor[] filteredPrecursors;
		public final HashMap<String, double[]> genomeScoreMap;
		public final HashMap<String, double[]> proteinScoreMap;

		private DiaResult(DiaNNPrecursor[] filteredPrecursors, HashMap<String, double[]> genomeScoreMap,
				HashMap<String, double[]> proteinScoreMap) {
			this.filteredPrecursors = filteredPrecursors;
			this.genomeScoreMap = genomeScoreMap;
			this.proteinScoreMap = proteinScoreMap;
		}
	}

	public static DdaResult scoreDda(MetaPeptide[] peptides, MetaPSM[] psms, String decoySymbol,
			String magIdentifier) {
		return scoreDda(peptides, psms, peptides, psms, decoySymbol, magIdentifier);
	}

	/**
	 * @param genomeScoringPeptides {@code peptides} plus any decoy peptides that
	 *                              only the genome null should see
	 * @param genomeScoringPsms     {@code psms} plus the matching decoy PSMs
	 */
	public static DdaResult scoreDda(MetaPeptide[] peptides, MetaPSM[] psms, MetaPeptide[] genomeScoringPeptides,
			MetaPSM[] genomeScoringPsms, String decoySymbol, String magIdentifier) {
		DdaMagGenomeScoringService.Result genomeScoring = scoreDdaGenomes(genomeScoringPeptides, genomeScoringPsms,
				decoySymbol, magIdentifier);
		HashMap<String, Double> confidentGenomeScoreMap = filterConfidentGenomeScores(genomeScoring.genomeScoreMap,
				genomeScoring.diaStyleGenomeStats);
		DdaMagProteinInferenceService.Result proteinInference = inferDdaProteins(peptides, confidentGenomeScoreMap,
				genomeScoring.diaStyleGenomeStats, magIdentifier, decoySymbol);
		HashMap<String, double[]> proteinScoreMap = DdaMagGenomeScoringService.computeProteinScoreMap(peptides, psms,
				decoySymbol);
		return new DdaResult(genomeScoring, confidentGenomeScoreMap, proteinInference, proteinScoreMap);
	}

	/**
	 * Razor assignment and protein statistics from a {@link GenomeEvidenceTest}
	 * result: its genome stats keep the {@code {score, p, q}} layout, so the
	 * confident genomes are those it accepted.
	 */
	public static DdaResult scoreDda(MetaPeptide[] peptides, GenomeEvidenceTest.Result evidence, String decoySymbol,
			String magIdentifier) {
		DdaMagGenomeScoringService.Result genomeScoring = DdaMagGenomeScoringService.Result
				.fromGenomeEvidence(evidence.genomeScoreMap, evidence.genomeStats);
		HashMap<String, Double> confidentGenomeScoreMap = filterConfidentGenomeScores(evidence.genomeScoreMap,
				evidence.genomeStats);
		// With the protein filter (identification v2), razor assignment considers passing proteins only.
		DdaMagProteinInferenceService.Result proteinInference = DdaMagProteinInferenceService.refineProteins(peptides,
				confidentGenomeScoreMap, evidence.genomeStats, magIdentifier, decoySymbol,
				DdaMagProteinInferenceService.DEFAULT_GENOME_Q_VALUE_THRESHOLD, evidence.passingProteins);
		return new DdaResult(genomeScoring, confidentGenomeScoreMap, proteinInference, evidence.proteinStats);
	}

	public static DdaMagGenomeScoringService.Result scoreDdaGenomes(MetaPeptide[] peptides, MetaPSM[] psms,
			String decoySymbol, String magIdentifier) {
		return DdaMagGenomeScoringService.refineGenomes(peptides, decoySymbol, psms, magIdentifier);
	}

	public static DdaMagProteinInferenceService.Result inferDdaProteins(MetaPeptide[] peptides,
			HashMap<String, Double> genomeScoreMap, HashMap<String, double[]> genomeStats, String magIdentifier,
			String decoySymbol) {
		return DdaMagProteinInferenceService.refineProteins(peptides, genomeScoreMap, genomeStats, magIdentifier,
				decoySymbol);
	}

	public static DiaResult scoreDia(DiaNNPrecursor[] precursors) {
		DiaMagScoringService.FilteredScoringResult scoring = DiaMagScoringService.filterAndScorePrecursors(precursors);
		return new DiaResult(scoring.filteredPrecursors, scoring.genomeScoreMap, scoring.proteinScoreMap);
	}

	public static HashMap<String, Double> filterConfidentGenomeScores(HashMap<String, Double> genomeScoreMap,
			HashMap<String, double[]> genomeStats) {
		HashMap<String, Double> out = new HashMap<>();
		if (genomeScoreMap == null) {
			return out;
		}
		Set<String> confident = DdaMagProteinInferenceService.confidentGenomes(genomeScoreMap, genomeStats,
				DdaMagProteinInferenceService.DEFAULT_GENOME_Q_VALUE_THRESHOLD);
		for (String genome : confident) {
			if (genomeScoreMap.containsKey(genome)) {
				out.put(genome, genomeScoreMap.get(genome));
			}
		}
		return out;
	}
}
