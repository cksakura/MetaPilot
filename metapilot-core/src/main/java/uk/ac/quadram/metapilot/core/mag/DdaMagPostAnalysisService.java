package uk.ac.quadram.metapilot.core.mag;

import uk.ac.quadram.metapilot.task.io.pep.MetaPSM;
import uk.ac.quadram.metapilot.task.io.pep.MetaPeptide;
import uk.ac.quadram.metapilot.task.io.pro.MetaProtein;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Shared, platform-independent post-analysis for DDA MAG workflows.
 */
public final class DdaMagPostAnalysisService {

	private DdaMagPostAnalysisService() {
	}

	private static <T> T[] concat(T[] first, T[] second, java.util.function.IntFunction<T[]> newArray) {
		if (second == null || second.length == 0) {
			return first;
		}
		T[] out = newArray.apply(first.length + second.length);
		System.arraycopy(first, 0, out, 0, first.length);
		System.arraycopy(second, 0, out, first.length, second.length);
		return out;
	}

	@FunctionalInterface
	public interface FileIndexResolver {
		int resolve(String fileName, Map<String, Integer> fileIdMap, List<String> fileList);
	}

	public static final class Input {
		public final MetaPeptide[] peptides;
		public final MetaPSM[] psms;
		public final HashMap<String, Integer> fileIdMap;
		public final ArrayList<String> fileList;
		public final String decoySymbol;
		public final String magIdentifier;
		public final FileIndexResolver fileIndexResolver;
		public final boolean preferPeptideMs2Counts;
		/**
		 * Decoy peptides and PSMs seen only by genome scoring, as its null
		 * background; they are never counted, reported or used for protein
		 * inference. For callers whose {@link #psms} hold targets only. Null when
		 * {@link #psms} already carry the decoys.
		 */
		public final MetaPeptide[] genomeBackgroundPeptides;
		public final MetaPSM[] genomeBackgroundPsms;
		/**
		 * Genome and protein statistics already computed by
		 * {@link GenomeEvidenceTest}; when set they replace the resampling genome
		 * test and protein scores, and the genome background is not used.
		 */
		public final GenomeEvidenceTest.Result genomeEvidence;

		public Input(MetaPeptide[] peptides, MetaPSM[] psms, HashMap<String, Integer> fileIdMap,
				ArrayList<String> fileList, String decoySymbol, String magIdentifier) {
			this(peptides, psms, fileIdMap, fileList, decoySymbol, magIdentifier, null, false);
		}

		public Input(MetaPeptide[] peptides, MetaPSM[] psms, HashMap<String, Integer> fileIdMap,
				ArrayList<String> fileList, String decoySymbol, String magIdentifier,
				FileIndexResolver fileIndexResolver, boolean preferPeptideMs2Counts) {
			this(peptides, psms, fileIdMap, fileList, decoySymbol, magIdentifier, fileIndexResolver,
					preferPeptideMs2Counts, null, null);
		}

		public Input(MetaPeptide[] peptides, MetaPSM[] psms, HashMap<String, Integer> fileIdMap,
				ArrayList<String> fileList, String decoySymbol, String magIdentifier,
				FileIndexResolver fileIndexResolver, boolean preferPeptideMs2Counts,
				MetaPeptide[] genomeBackgroundPeptides, MetaPSM[] genomeBackgroundPsms) {
			this(peptides, psms, fileIdMap, fileList, decoySymbol, magIdentifier, fileIndexResolver,
					preferPeptideMs2Counts, genomeBackgroundPeptides, genomeBackgroundPsms, null);
		}

		public Input(MetaPeptide[] peptides, MetaPSM[] psms, HashMap<String, Integer> fileIdMap,
				ArrayList<String> fileList, String decoySymbol, String magIdentifier,
				FileIndexResolver fileIndexResolver, boolean preferPeptideMs2Counts,
				MetaPeptide[] genomeBackgroundPeptides, MetaPSM[] genomeBackgroundPsms,
				GenomeEvidenceTest.Result genomeEvidence) {
			this.genomeEvidence = genomeEvidence;
			this.genomeBackgroundPeptides = genomeBackgroundPeptides;
			this.genomeBackgroundPsms = genomeBackgroundPsms;
			this.peptides = peptides;
			this.psms = psms;
			this.fileIdMap = fileIdMap;
			this.fileList = fileList;
			this.decoySymbol = decoySymbol;
			this.magIdentifier = magIdentifier;
			this.fileIndexResolver = fileIndexResolver;
			this.preferPeptideMs2Counts = preferPeptideMs2Counts;
		}
	}

	public static final class Result {
		public MagResultFinalisationService.Selection finalSelection;
		public final HashMap<String, MetaPSM> psmMap;
		public final HashMap<String, int[]> pepPsmCountMap;
		public final HashMap<String, int[]> pepChargeMap;
		public final HashMap<String, double[]> pepExpectMap;
		public final HashMap<String, Integer> psmCountMap;
		public final int[] filePsmCount;
		public final HashSet<String>[] fileSequenceSets;
		public final DdaMagGenomeScoringService.Result genomeScoring;
		public final DdaMagProteinInferenceService.Result proteinRefinement;
		public final HashMap<String, Double> genomeScoreMap;
		public final HashMap<String, double[]> genomeDiaStats;
		public final String[] razorProteins;
		public HashMap<String, double[]> genomeIntensityMap;
		/** Intensity summed from peptides that map to exactly one confident genome. */
		public HashMap<String, double[]> genomeUniqueIntensityMap = new HashMap<>();
		public final HashMap<String, double[]> proteinIntensityMap;
		public final HashMap<String, Double> proteinTotalIntenMap;
		public final HashMap<String, double[]> allProteinIntensityMap;
		public final HashMap<String, Double> allProteinTotalIntenMap;
		public final HashMap<String, HashSet<String>> allProteinPepCountMap;
		public final HashMap<String, Integer> allProteinPsmCountMap;
		public final HashMap<String, double[]> allProteinScoreMap;
		public final HashMap<String, HashSet<String>> totalPepCountMap;
		public final HashMap<String, HashSet<String>> razorPepCountMap;
		public final HashMap<String, HashSet<String>> uniquePepCountMap;
		public final HashMap<String, HashSet<String>> totalPepWithILMap;
		public final HashMap<String, HashSet<String>> razorPepWithILMap;
		public final HashMap<String, HashSet<String>> uniquePepWithILMap;
		public final ArrayList<MetaPeptide> peptideList;
		public final String[] proteins;
		public final String[] allProteins;
		public final HashMap<String, Integer> proteinIdMap;
		public final HashMap<String, Integer> razorPsmCountMap;
		public final MetaProtein[] metaProteins;
		/** Diagnostic count; excluded evidence never enters final peptideList. */
		public final int peptidesWithoutConfidentRazorProtein;

		private Result(HashMap<String, MetaPSM> psmMap, HashMap<String, int[]> pepPsmCountMap,
				HashMap<String, int[]> pepChargeMap, HashMap<String, double[]> pepExpectMap,
				HashMap<String, Integer> psmCountMap, int[] filePsmCount, HashSet<String>[] fileSequenceSets,
				DdaMagGenomeScoringService.Result genomeScoring, HashMap<String, Double> genomeScoreMap,
				DdaMagProteinInferenceService.Result proteinRefinement, HashMap<String, double[]> genomeIntensityMap,
				HashMap<String, double[]> proteinIntensityMap, HashMap<String, Double> proteinTotalIntenMap,
				HashMap<String, double[]> allProteinIntensityMap, HashMap<String, Double> allProteinTotalIntenMap,
				HashMap<String, HashSet<String>> allProteinPepCountMap,
				HashMap<String, Integer> allProteinPsmCountMap,
				HashMap<String, double[]> allProteinScoreMap,
				HashMap<String, HashSet<String>> totalPepCountMap,
				HashMap<String, HashSet<String>> razorPepCountMap,
				HashMap<String, HashSet<String>> uniquePepCountMap,
				HashMap<String, HashSet<String>> totalPepWithILMap,
				HashMap<String, HashSet<String>> razorPepWithILMap,
				HashMap<String, HashSet<String>> uniquePepWithILMap, ArrayList<MetaPeptide> peptideList,
				String[] proteins, String[] allProteins, HashMap<String, Integer> proteinIdMap,
				HashMap<String, Integer> razorPsmCountMap, MetaProtein[] metaProteins,
				int peptidesWithoutConfidentRazorProtein) {
			this.psmMap = psmMap;
			this.pepPsmCountMap = pepPsmCountMap;
			this.pepChargeMap = pepChargeMap;
			this.pepExpectMap = pepExpectMap;
			this.psmCountMap = psmCountMap;
			this.filePsmCount = filePsmCount;
			this.fileSequenceSets = fileSequenceSets;
			this.genomeScoring = genomeScoring;
			this.proteinRefinement = proteinRefinement;
			this.genomeScoreMap = genomeScoreMap;
			this.genomeDiaStats = genomeScoring.diaStyleGenomeStats;
			this.razorProteins = proteinRefinement.razorProteins;
			this.genomeIntensityMap = genomeIntensityMap;
			this.proteinIntensityMap = proteinIntensityMap;
			this.proteinTotalIntenMap = proteinTotalIntenMap;
			this.allProteinIntensityMap = allProteinIntensityMap;
			this.allProteinTotalIntenMap = allProteinTotalIntenMap;
			this.allProteinPepCountMap = allProteinPepCountMap;
			this.allProteinPsmCountMap = allProteinPsmCountMap;
			this.allProteinScoreMap = allProteinScoreMap;
			this.totalPepCountMap = totalPepCountMap;
			this.razorPepCountMap = razorPepCountMap;
			this.uniquePepCountMap = uniquePepCountMap;
			this.totalPepWithILMap = totalPepWithILMap;
			this.razorPepWithILMap = razorPepWithILMap;
			this.uniquePepWithILMap = uniquePepWithILMap;
			this.peptideList = peptideList;
			this.proteins = proteins;
			this.allProteins = allProteins;
			this.proteinIdMap = proteinIdMap;
			this.razorPsmCountMap = razorPsmCountMap;
			this.metaProteins = metaProteins;
			this.peptidesWithoutConfidentRazorProtein = peptidesWithoutConfidentRazorProtein;
		}
	}

	public static Result analyze(Input input) {
		if (input == null) {
			throw new IllegalArgumentException("input is null");
		}
		MetaPeptide[] peptides = input.peptides != null ? input.peptides : new MetaPeptide[0];
		MetaPSM[] psms = input.psms != null ? input.psms : new MetaPSM[0];
		HashMap<String, Integer> fileIdMap = input.fileIdMap != null ? input.fileIdMap : new HashMap<>();
		ArrayList<String> fileList = input.fileList != null ? input.fileList : new ArrayList<>();
		String decoySymbol = input.decoySymbol != null ? input.decoySymbol : "";
		String magIdentifier = input.magIdentifier != null ? input.magIdentifier : "";
		int fileCount = Math.max(fileIdMap.size(), fileList.size());

		HashMap<String, MetaPSM> psmMap = new HashMap<>();
		HashMap<String, int[]> pepPsmCountMap = new HashMap<>();
		HashMap<String, int[]> pepChargeMap = new HashMap<>();
		HashMap<String, double[]> pepExpectMap = new HashMap<>();
		HashMap<String, Integer> psmCountMap = new HashMap<>();

		int[] filePsmCount = new int[fileCount];
		@SuppressWarnings("unchecked")
		HashSet<Integer>[] fileScanSets = new HashSet[fileCount];
		@SuppressWarnings("unchecked")
		HashSet<String>[] fileSequenceSets = new HashSet[fileCount];
		for (int i = 0; i < fileCount; i++) {
			fileScanSets[i] = new HashSet<>();
			fileSequenceSets[i] = new HashSet<>();
		}

		for (MetaPSM psm : psms) {
			if (psm == null) {
				continue;
			}
			String modseq = psm.getModSeq();
			if (modseq == null) {
				continue;
			}
			double pep = psm.getPEP();
			int fileId = resolveFileIndex(psm.getFileName(), fileIdMap, fileList, input.fileIndexResolver);

			if (isValidFileId(fileId, fileCount)) {
				int[] psmCount = pepPsmCountMap.computeIfAbsent(modseq, k -> new int[fileCount]);
				psmCount[fileId]++;
				// Decoy PSMs stay in the protein/genome counts (needed for genome FDR) but are not
				// identifications, so they are left out of the per-file summary counts.
				if (!mapsOnlyToDecoys(psm.getProteins(), decoySymbol)) {
					filePsmCount[fileId]++;
					fileScanSets[fileId].add(psm.getScan());
					String sequence = psm.getBaseSequence();
					if (sequence == null) {
						sequence = psm.getSequence();
					}
					if (sequence != null) {
						fileSequenceSets[fileId].add(sequence);
					}
				}

				double[] expect = pepExpectMap.computeIfAbsent(modseq, k -> {
					double[] values = new double[fileCount];
					Arrays.fill(values, 1.0);
					return values;
				});
				if (pep < expect[fileId]) {
					expect[fileId] = pep;
				}
			} else {
				pepPsmCountMap.computeIfAbsent(modseq, k -> new int[fileCount]);
			}

			int[] charges = pepChargeMap.computeIfAbsent(modseq, k -> new int[7]);
			int charge = psm.getCharge();
			if (charge >= 0 && charge < charges.length) {
				charges[charge] = 1;
			}

			if (!psmMap.containsKey(modseq) || pep < psmMap.get(modseq).getPEP()) {
				psmMap.put(modseq, psm);
			}

			String[] proteins = psm.getProteins();
			for (String protein : proteins != null ? proteins : new String[0]) {
				if (protein == null) {
					continue;
				}
				psmCountMap.put(protein, psmCountMap.getOrDefault(protein, 0) + 1);
				if (isMagProtein(protein, magIdentifier)) {
					String genome = genomeOf(protein);
					psmCountMap.put(genome, psmCountMap.getOrDefault(genome, 0) + 1);
				}
			}
		}

		for (int i = 0; i < fileCount; i++) {
			if (!(fileScanSets[i].size() == 1 && filePsmCount[i] > 1)) {
				filePsmCount[i] = fileScanSets[i].size();
			}
		}
		if (input.preferPeptideMs2Counts && hasPeptideMs2Counts(peptides)) {
			pepPsmCountMap.clear();
			Arrays.fill(filePsmCount, 0);
			for (HashSet<String> fileSequenceSet : fileSequenceSets) {
				fileSequenceSet.clear();
			}
			applyPeptideMs2Counts(peptides, pepPsmCountMap, filePsmCount, fileSequenceSets, fileCount, decoySymbol);
		} else if (input.preferPeptideMs2Counts && sum(filePsmCount) == 0) {
			applyPeptideMs2Counts(peptides, pepPsmCountMap, filePsmCount, fileSequenceSets, fileCount, decoySymbol);
		}

		for (MetaPeptide peptide : peptides) {
			if (peptide == null) {
				continue;
			}
			String modseq = peptide.getModSeq();
			if (modseq != null && pepExpectMap.containsKey(modseq)) {
				peptide.setIndividualPEP(pepExpectMap.get(modseq));
			} else {
				double[] expect = new double[fileCount];
				Arrays.fill(expect, 1.0);
				peptide.setIndividualPEP(expect);
			}
		}

//		enrichPeptideProteinsFromPsms(peptides, psms);

		MagPostAnalysisScoringService.DdaResult scoring = input.genomeEvidence != null
				? MagPostAnalysisScoringService.scoreDda(peptides, input.genomeEvidence, decoySymbol, magIdentifier)
				: MagPostAnalysisScoringService.scoreDda(peptides, psms,
						concat(peptides, input.genomeBackgroundPeptides, MetaPeptide[]::new),
						concat(psms, input.genomeBackgroundPsms, MetaPSM[]::new), decoySymbol, magIdentifier);
		Result result = buildResult(peptides, decoySymbol, magIdentifier, psmMap, pepPsmCountMap, pepChargeMap, pepExpectMap,
				psmCountMap, filePsmCount, fileSequenceSets, scoring.genomeScoring, scoring.confidentGenomeScoreMap,
				scoring.proteinInference, scoring.proteinScoreMap);
        if (!input.preferPeptideMs2Counts) {
            Set<String> retainedKeys = new HashSet<>();
            for (MetaPeptide p : result.peptideList) retainedKeys.add(p.getModSeq());
            Arrays.fill(result.filePsmCount, 0);
            for (int i = 0; i < fileCount; i++) { fileScanSets[i].clear(); result.fileSequenceSets[i].clear(); }
            int[] rawCounts = new int[fileCount];
            for (MetaPSM psm : psms) {
                if (psm == null || !retainedKeys.contains(psm.getModSeq()) || mapsOnlyToDecoys(psm.getProteins(), decoySymbol)) continue;
                int file = resolveFileIndex(psm.getFileName(), fileIdMap, fileList, input.fileIndexResolver);
                if (!isValidFileId(file, fileCount)) continue;
                rawCounts[file]++;
                fileScanSets[file].add(psm.getScan());
                String seq = psm.getBaseSequence() != null ? psm.getBaseSequence() : psm.getSequence();
                if (seq != null) result.fileSequenceSets[file].add(seq);
            }
            for (int i = 0; i < fileCount; i++) result.filePsmCount[i] = fileScanSets[i].size() == 1 && rawCounts[i] > 1
                    ? rawCounts[i] : fileScanSets[i].size();
        }
        return result;

	}

	/**
	 * Post-analysis for peptide input that carries no PSMs, scores or decoys (a plain
	 * peptide intensity table mapped to MAG proteins). Genome scoring is skipped: the
	 * caller supplies the accepted genome set, every accepted genome is treated as
	 * confident, and razor assignment plus genome/protein intensity roll-up then run
	 * through the same {@link DdaMagProteinInferenceService#refineProteins} and
	 * {@code buildResult} path as {@link #analyze}. PSM-derived maps in the result are
	 * empty.
	 *
	 * @param peptides         peptides with {@code getProteins()} already set to their
	 *                         mapped MAG proteins and {@code getIntensity()} of length
	 *                         {@code fileCount}
	 * @param acceptedGenomes  genomes allowed to receive razor peptides and intensity
	 */
	public static Result analyzeWithAcceptedGenomes(MetaPeptide[] peptides, int fileCount, String magIdentifier,
			Set<String> acceptedGenomes) {
		MetaPeptide[] peps = peptides != null ? peptides : new MetaPeptide[0];
		String magId = magIdentifier != null ? magIdentifier : "";
		HashMap<String, Double> genomeScoreMap = new HashMap<>();
		if (acceptedGenomes != null) {
			for (String genome : acceptedGenomes) {
				if (genome != null) {
					genomeScoreMap.put(genome, 1.0);
				}
			}
		}

		@SuppressWarnings("unchecked")
		HashSet<String>[] fileSequenceSets = new HashSet[fileCount];
		for (int i = 0; i < fileCount; i++) {
			fileSequenceSets[i] = new HashSet<>();
		}
		for (MetaPeptide peptide : peps) {
			if (peptide == null || peptide.getSequence() == null || peptide.getIntensity() == null) {
				continue;
			}
			double[] intensity = peptide.getIntensity();
			for (int i = 0; i < Math.min(fileCount, intensity.length); i++) {
				if (intensity[i] > 0) {
					fileSequenceSets[i].add(peptide.getSequence());
				}
			}
		}

		// decoySymbol == null: the input has no decoys, and both services treat a null
		// decoy symbol as "nothing is a decoy" (an empty string would match everything
		// in DdaMagProteinInferenceService.isDecoyProtein).
		DdaMagProteinInferenceService.Result refinement = DdaMagProteinInferenceService.refineProteins(peps,
				genomeScoreMap, null, magId, null);
		DdaMagGenomeScoringService.Result genomeScoring = DdaMagGenomeScoringService.Result
				.acceptedGenomesOnly(genomeScoreMap);
		return buildResult(peps, null, magId, new HashMap<>(), new HashMap<>(), new HashMap<>(), new HashMap<>(),
				new HashMap<>(), new int[fileCount], fileSequenceSets, genomeScoring, genomeScoreMap, refinement,
				new HashMap<>());
	}

	/**
	 * Not used currently
	 * @param peptides
	 * @param psms
	 */
	private static void enrichPeptideProteinsFromPsms(MetaPeptide[] peptides, MetaPSM[] psms) {
		if (peptides == null || psms == null || psms.length == 0) {
			return;
		}
		HashMap<String, LinkedHashSet<String>> proteinsByKey = new HashMap<>();
		for (MetaPSM psm : psms) {
			if (psm == null || psm.getProteins() == null || psm.getProteins().length == 0) {
				continue;
			}
			addProteinKeys(proteinsByKey, psm.getModSeq(), psm.getProteins());
			addProteinKeys(proteinsByKey, psm.getBaseSequence(), psm.getProteins());
			addProteinKeys(proteinsByKey, psm.getSequence(), psm.getProteins());
		}

		for (MetaPeptide peptide : peptides) {
			if (peptide == null) {
				continue;
			}
			LinkedHashSet<String> merged = new LinkedHashSet<>();
			addProteins(merged, peptide.getProteins());
			addProteins(merged, proteinsByKey.get(peptide.getModSeq()));
			addProteins(merged, proteinsByKey.get(peptide.getSequence()));
			if (!merged.isEmpty()) {
				peptide.setProteins(merged.toArray(new String[0]));
			}
		}
	}

	private static void addProteinKeys(HashMap<String, LinkedHashSet<String>> proteinsByKey, String key,
			String[] proteins) {
		if (key == null || key.isBlank()) {
			return;
		}
		addProteins(proteinsByKey.computeIfAbsent(key, k -> new LinkedHashSet<>()), proteins);
	}

	private static void addProteins(LinkedHashSet<String> target, String[] proteins) {
		if (target == null || proteins == null) {
			return;
		}
		for (String protein : proteins) {
			if (protein != null && !protein.trim().isEmpty()) {
				target.add(protein.trim());
			}
		}
	}

	private static void addProteins(LinkedHashSet<String> target, Set<String> proteins) {
		if (target == null || proteins == null) {
			return;
		}
		for (String protein : proteins) {
			if (protein != null && !protein.trim().isEmpty()) {
				target.add(protein.trim());
			}
		}
	}

	private static Result buildResult(MetaPeptide[] peptides, String decoySymbol, String magIdentifier,
			HashMap<String, MetaPSM> psmMap, HashMap<String, int[]> pepPsmCountMap,
			HashMap<String, int[]> pepChargeMap, HashMap<String, double[]> pepExpectMap,
			HashMap<String, Integer> psmCountMap, int[] filePsmCount, HashSet<String>[] fileSequenceSets,
			DdaMagGenomeScoringService.Result genomeScoring, HashMap<String, Double> confidentGenomeScoreMap,
			DdaMagProteinInferenceService.Result proteinRefinement, HashMap<String, double[]> allProteinScoreMap) {
        Map<MetaPeptide, String[]> originalMappings = new java.util.IdentityHashMap<>();
        for (MetaPeptide p : peptides) {
            if (p == null) continue;
            originalMappings.put(p, p.getProteins());
            p.setProteins(MagResultFinalisationService.retainedProteins(p.getProteins(), magIdentifier,
                    decoySymbol, confidentGenomeScoreMap.keySet()));
        }
        DdaMagProteinInferenceService.Result gatedRefinement = DdaMagProteinInferenceService.refineProteins(
                peptides, confidentGenomeScoreMap, null, magIdentifier, decoySymbol);
        Result provisional = buildProvisionalResult(peptides, decoySymbol, magIdentifier, psmMap,
                pepPsmCountMap, pepChargeMap, pepExpectMap, psmCountMap, filePsmCount, fileSequenceSets,
                genomeScoring, confidentGenomeScoreMap, gatedRefinement, allProteinScoreMap);
        for (Map.Entry<MetaPeptide, String[]> e : originalMappings.entrySet()) e.getKey().setProteins(e.getValue());
        MagResultFinalisationService.Selection selection = MagResultFinalisationService.select(
                provisional.genomeIntensityMap, genomeScoring.diaStyleGenomeStats,
                confidentGenomeScoreMap.keySet(), psmCountMap, provisional.totalPepCountMap,
                provisional.razorPepCountMap, provisional.uniquePepCountMap);
        MetaPeptide[] retained = MagResultFinalisationService.prune(peptides, magIdentifier, decoySymbol, selection);
        HashMap<String, Double> frozenScores = new HashMap<>(confidentGenomeScoreMap);
        frozenScores.keySet().retainAll(selection.genomes);
        DdaMagProteinInferenceService.Result refined = DdaMagProteinInferenceService.refineProteins(
                retained, frozenScores, null, magIdentifier, decoySymbol);
        Arrays.fill(filePsmCount, 0);
        for (HashSet<String> seqs : fileSequenceSets) seqs.clear();
        psmCountMap.clear();
        for (MetaPeptide p : retained) {
            if (p == null) continue;
            int[] counts = pepPsmCountMap.get(p.getModSeq());
            double[] intensity = p.getIntensity();
            for (int i = 0; i < filePsmCount.length; i++) {
                if (counts != null && i < counts.length) filePsmCount[i] += counts[i];
                // "Peptide Sequences Identified" per file: MS/MS identifications only. A peptide quantified
                // in a file by match-between-runs is not identified there. Intensity stands in only when
                // the input carries no MS2 counts at all (peptide-table workflows).
                boolean identified = counts != null
                        ? i < counts.length && counts[i] > 0
                        : intensity != null && i < intensity.length && intensity[i] > 0;
                if (identified) fileSequenceSets[i].add(p.getSequence());
            }
            int n = counts == null ? 0 : Arrays.stream(counts).sum();
            Set<String> gs = new HashSet<>();
            for (String pro : p.getProteins()) {
                psmCountMap.merge(pro, n, Integer::sum);
                if (isMagProtein(pro, magIdentifier)) gs.add(genomeOf(pro));
            }
            for (String g : gs) psmCountMap.merge(g, n, Integer::sum);
        }
        Result result = buildProvisionalResult(retained, decoySymbol, magIdentifier, psmMap,
                pepPsmCountMap, pepChargeMap, pepExpectMap, psmCountMap, filePsmCount, fileSequenceSets,
                genomeScoring, frozenScores, refined, allProteinScoreMap);
        result.finalSelection = selection;
        result.genomeIntensityMap = new MagResultFinalisationService.FinalGenomeIntensityMap(result.genomeIntensityMap, selection);
        result.genomeUniqueIntensityMap = new MagResultFinalisationService.FinalGenomeIntensityMap(
                result.genomeUniqueIntensityMap, selection);
        return result;
    }

	private static Result buildProvisionalResult(MetaPeptide[] peptides, String decoySymbol, String magIdentifier,
			HashMap<String, MetaPSM> psmMap, HashMap<String, int[]> pepPsmCountMap,
			HashMap<String, int[]> pepChargeMap, HashMap<String, double[]> pepExpectMap,
			HashMap<String, Integer> psmCountMap, int[] filePsmCount, HashSet<String>[] fileSequenceSets,
			DdaMagGenomeScoringService.Result genomeScoring, HashMap<String, Double> confidentGenomeScoreMap,
			DdaMagProteinInferenceService.Result proteinRefinement, HashMap<String, double[]> allProteinScoreMap) {
		HashMap<String, Double> genomeScoreMap = confidentGenomeScoreMap;
		String[] razorProteins = proteinRefinement.razorProteins;
		HashMap<String, double[]> genomeIntensityMap = new HashMap<>();
		HashMap<String, double[]> genomeUniqueIntensityMap = new HashMap<>();
		HashMap<String, double[]> proteinIntensityMap = new HashMap<>();
		HashMap<String, Double> proteinTotalIntenMap = new HashMap<>();
		HashMap<String, double[]> allProteinIntensityMap = new HashMap<>();
		HashMap<String, Double> allProteinTotalIntenMap = new HashMap<>();
		HashMap<String, HashSet<String>> allProteinPepCountMap = new HashMap<>();
		HashMap<String, Integer> allProteinPsmCountMap = new HashMap<>();
		HashMap<String, HashSet<String>> totalPepCountMap = new HashMap<>();
		HashMap<String, HashSet<String>> razorPepCountMap = new HashMap<>();
		HashMap<String, HashSet<String>> uniquePepCountMap = new HashMap<>();
		HashMap<String, HashSet<String>> totalPepWithILMap = new HashMap<>();
		HashMap<String, HashSet<String>> razorPepWithILMap = new HashMap<>();
		HashMap<String, HashSet<String>> uniquePepWithILMap = new HashMap<>();
		ArrayList<MetaPeptide> peptideList = new ArrayList<>();

		// Final tables require a retained razor protein; excluded evidence is exported separately.
		int peptidesWithoutConfidentRazorProtein = 0;
		for (int i = 0; i < peptides.length; i++) {
			if (peptides[i] == null) {
				continue;
			}
			String sequence = peptides[i].getSequence();
			if (sequence == null) {
				continue;
			}
			if (mapsOnlyToDecoys(peptides[i], decoySymbol)) {
				continue;
			}
			if (i >= razorProteins.length || razorProteins[i] == null) {
				peptidesWithoutConfidentRazorProtein++;
				continue;
			}
			peptideList.add(peptides[i]);
			boolean hasIorL = hasIorL(sequence);
			razorPepCountMap.computeIfAbsent(razorProteins[i], k -> new HashSet<>()).add(sequence);
			if (isMagProtein(razorProteins[i], magIdentifier)) {
				String genome = genomeOf(razorProteins[i]);
				razorPepCountMap.computeIfAbsent(genome, k -> new HashSet<>()).add(sequence);
				if (hasIorL) {
					razorPepWithILMap.computeIfAbsent(genome, k -> new HashSet<>()).add(sequence);
				}
			}
		}

		for (MetaPeptide peptide : peptides) {
			if (peptide == null || peptide.getProteins() == null) {
				continue;
			}
			String[] proteins = peptide.getProteins();
			int count = 0;
			for (String protein : proteins) {
				if (protein != null && razorPepCountMap.containsKey(protein)) {
					count++;
				}
			}
			if (count <= 0) {
				continue;
			}
			String sequence = peptide.getSequence();
			if (sequence == null) {
				continue;
			}
			addAllProteinEvidence(peptide, decoySymbol, allProteinIntensityMap, allProteinTotalIntenMap,
					allProteinPepCountMap, allProteinPsmCountMap, pepPsmCountMap);
            // Support counts follow every retained association, including proteins
            // represented only in all_proteins.tsv. Allocation still follows the
            // established razor-bearing-protein policy below.
            Set<String> mappedGenomes = new HashSet<>();
            boolean hasNonMag = false;
            for (String protein : proteins) {
                if (protein == null || isDecoyProtein(protein, decoySymbol)) continue;
                totalPepCountMap.computeIfAbsent(protein, k -> new HashSet<>()).add(sequence);
                if (isMagProtein(protein, magIdentifier) && genomeScoreMap.containsKey(genomeOf(protein))) {
                    String genome = genomeOf(protein);
                    mappedGenomes.add(genome);
                    totalPepCountMap.computeIfAbsent(genome, k -> new HashSet<>()).add(sequence);
                    if (hasIorL(sequence)) totalPepWithILMap.computeIfAbsent(genome, k -> new HashSet<>()).add(sequence);
                } else hasNonMag = true;
            }
            if (!hasNonMag && mappedGenomes.size() == 1) {
                String genome = mappedGenomes.iterator().next();
                uniquePepCountMap.computeIfAbsent(genome, k -> new HashSet<>()).add(sequence);
                if (hasIorL(sequence)) uniquePepWithILMap.computeIfAbsent(genome, k -> new HashSet<>()).add(sequence);
            }
			boolean hasIorL = hasIorL(sequence);
			String uniqueGenome = null;
			boolean uniqueToOneGenome = true;
			for (String protein : proteins) {
				if (isMagProtein(protein, magIdentifier)) {
					String genome = genomeOf(protein);
					if (uniqueGenome == null) {
						uniqueGenome = genome;
					} else if (!uniqueGenome.equals(genome)) {
						uniqueToOneGenome = false;
						break;
					}
				} else {
					uniqueToOneGenome = false;
					break;
				}
			}

			double[] intensity = peptide.getIntensity();
			if (intensity == null) {
				intensity = new double[filePsmCount.length];
			}
			if (!hasNonMag && mappedGenomes.size() == 1) {
				addIntensity(genomeUniqueIntensityMap, mappedGenomes.iterator().next(), intensity);
			}
			double[] sharedIntensity = new double[intensity.length];
			for (int i = 0; i < sharedIntensity.length; i++) {
				sharedIntensity[i] = intensity[i] / (double) count;
			}

			for (String protein : proteins) {
				if (protein == null || !razorPepCountMap.containsKey(protein) || isDecoyProtein(protein, decoySymbol)) {
					continue;
				}
				if (isMagProtein(protein, magIdentifier)) {
					String genome = genomeOf(protein);
					if (!genomeScoreMap.containsKey(genome)) {
						continue;
					}
					addIntensity(genomeIntensityMap, genome, sharedIntensity);
					addProteinIntensity(proteinIntensityMap, proteinTotalIntenMap, protein, sharedIntensity);
					totalPepCountMap.computeIfAbsent(genome, k -> new HashSet<>()).add(sequence);
					totalPepCountMap.computeIfAbsent(protein, k -> new HashSet<>()).add(sequence);
					if (hasIorL) {
						totalPepWithILMap.computeIfAbsent(genome, k -> new HashSet<>()).add(sequence);
					}
					if (uniqueToOneGenome && uniqueGenome != null && uniqueGenome.equals(genome)) {
						uniquePepCountMap.computeIfAbsent(genome, k -> new HashSet<>()).add(sequence);
						if (hasIorL) {
							uniquePepWithILMap.computeIfAbsent(genome, k -> new HashSet<>()).add(sequence);
						}
					}
				} else {
					addProteinIntensity(proteinIntensityMap, proteinTotalIntenMap, protein, sharedIntensity);
					totalPepCountMap.computeIfAbsent(protein, k -> new HashSet<>()).add(sequence);
				}
			}
		}

		String[] proteins = proteinTotalIntenMap.keySet().toArray(new String[0]);
		Arrays.sort(proteins,
				Comparator.comparingDouble((String p) -> proteinTotalIntenMap.getOrDefault(p, 0.0)).reversed());
		String[] allProteins = allProteinTotalIntenMap.keySet().toArray(new String[0]);
		Arrays.sort(allProteins,
				Comparator.comparingDouble((String p) -> allProteinTotalIntenMap.getOrDefault(p, 0.0)).reversed());
		HashMap<String, Integer> proteinIdMap = new HashMap<>();
		for (int i = 0; i < proteins.length; i++) {
			proteinIdMap.put(proteins[i], i);
		}
		HashMap<String, Integer> razorPsmCountMap = buildRazorPsmCountMap(peptides, razorProteins, pepPsmCountMap);
		MetaProtein[] metaProteins = buildMetaProteins(proteins, proteinIntensityMap, razorPsmCountMap,
				razorPepCountMap);

		Result result = new Result(psmMap, pepPsmCountMap, pepChargeMap, pepExpectMap, psmCountMap, filePsmCount,
				fileSequenceSets, genomeScoring, genomeScoreMap, proteinRefinement, genomeIntensityMap,
				proteinIntensityMap, proteinTotalIntenMap, allProteinIntensityMap, allProteinTotalIntenMap,
				allProteinPepCountMap, allProteinPsmCountMap, allProteinScoreMap, totalPepCountMap,
				razorPepCountMap, uniquePepCountMap, totalPepWithILMap, razorPepWithILMap, uniquePepWithILMap,
				peptideList, proteins, allProteins,
				proteinIdMap, razorPsmCountMap, metaProteins, peptidesWithoutConfidentRazorProtein);
		result.genomeUniqueIntensityMap = genomeUniqueIntensityMap;
		return result;
	}

	private static void addAllProteinEvidence(MetaPeptide peptide, String decoySymbol,
			HashMap<String, double[]> allProteinIntensityMap, HashMap<String, Double> allProteinTotalIntenMap,
			HashMap<String, HashSet<String>> allProteinPepCountMap, HashMap<String, Integer> allProteinPsmCountMap,
			HashMap<String, int[]> pepPsmCountMap) {
		String[] proteins = peptide.getProteins();
		if (proteins == null || proteins.length == 0) {
			return;
		}
		ArrayList<String> targets = new ArrayList<>();
		for (String protein : proteins) {
			if (protein != null && !isDecoyProtein(protein, decoySymbol)) {
				targets.add(protein);
			}
		}
		if (targets.isEmpty()) {
			return;
		}
		double[] intensity = peptide.getIntensity();
		if (intensity == null) {
			intensity = new double[0];
		}
		double[] sharedIntensity = new double[intensity.length];
		for (int i = 0; i < intensity.length; i++) {
			sharedIntensity[i] = intensity[i] / (double) targets.size();
		}
		String sequence = peptide.getSequence();
		int psmCount = sum(pepPsmCountMap.get(peptide.getModSeq()));
		for (String protein : targets) {
			addProteinIntensity(allProteinIntensityMap, allProteinTotalIntenMap, protein, sharedIntensity);
			if (sequence != null) {
				allProteinPepCountMap.computeIfAbsent(protein, k -> new HashSet<>()).add(sequence);
			}
			allProteinPsmCountMap.put(protein, allProteinPsmCountMap.getOrDefault(protein, 0) + psmCount);
		}
	}

	private static void applyPeptideMs2Counts(MetaPeptide[] peptides, HashMap<String, int[]> pepPsmCountMap,
			int[] filePsmCount, HashSet<String>[] fileSequenceSets, int fileCount, String decoySymbol) {
		for (MetaPeptide peptide : peptides) {
			if (peptide == null || peptide.getModSeq() == null) {
				continue;
			}
			int[] counts = pepPsmCountMap.computeIfAbsent(peptide.getModSeq(), k -> new int[fileCount]);
			int[] ms2Counts = peptide.getMs2Counts();
			if (ms2Counts == null) {
				continue;
			}
			boolean decoy = mapsOnlyToDecoys(peptide, decoySymbol);
			for (int i = 0; i < Math.min(ms2Counts.length, counts.length); i++) {
				counts[i] = ms2Counts[i];
				if (decoy) {
					continue;
				}
				filePsmCount[i] += ms2Counts[i];
				if (ms2Counts[i] > 0 && peptide.getSequence() != null) {
					fileSequenceSets[i].add(peptide.getSequence());
				}
			}
		}
	}

	private static boolean hasPeptideMs2Counts(MetaPeptide[] peptides) {
		for (MetaPeptide peptide : peptides) {
			if (peptide == null || peptide.getMs2Counts() == null) {
				continue;
			}
			for (int count : peptide.getMs2Counts()) {
				if (count > 0) {
					return true;
				}
			}
		}
		return false;
	}

	private static HashMap<String, Integer> buildRazorPsmCountMap(MetaPeptide[] peptides, String[] razorProteins,
			HashMap<String, int[]> pepPsmCountMap) {
		HashMap<String, Integer> out = new HashMap<>();
		for (int i = 0; i < peptides.length; i++) {
			if (i >= razorProteins.length || peptides[i] == null || razorProteins[i] == null) {
				continue;
			}
			int total = 0;
			int[] counts = pepPsmCountMap.get(peptides[i].getModSeq());
			if (counts != null) {
				for (int count : counts) {
					total += count;
				}
			}
			out.merge(razorProteins[i], total, Integer::sum);
		}
		return out;
	}

	private static MetaProtein[] buildMetaProteins(String[] proteins, HashMap<String, double[]> proteinIntensityMap,
			HashMap<String, Integer> razorPsmCountMap, HashMap<String, HashSet<String>> razorPepCountMap) {
		ArrayList<MetaProtein> list = new ArrayList<>();
		for (int i = 0; i < proteins.length; i++) {
			String protein = proteins[i];
			double[] intensity = proteinIntensityMap.get(protein);
			if (intensity == null) {
				intensity = new double[0];
			}
			int pepCount = razorPepCountMap.containsKey(protein) ? razorPepCountMap.get(protein).size() : 0;
			int psmCount = razorPsmCountMap.getOrDefault(protein, 0);
			list.add(new MetaProtein(i, 0, protein, pepCount, psmCount, 0, intensity));
		}
		return list.toArray(new MetaProtein[0]);
	}

	private static void addIntensity(HashMap<String, double[]> map, String key, double[] values) {
		double[] target = map.computeIfAbsent(key, k -> new double[values.length]);
		for (int i = 0; i < values.length; i++) {
			target[i] += values[i];
		}
	}

	private static void addProteinIntensity(HashMap<String, double[]> proteinIntensityMap,
			HashMap<String, Double> proteinTotalIntenMap, String protein, double[] values) {
		double[] target = proteinIntensityMap.computeIfAbsent(protein, k -> new double[values.length]);
		double total = proteinTotalIntenMap.getOrDefault(protein, 0.0);
		for (int i = 0; i < values.length; i++) {
			target[i] += values[i];
			total += values[i];
		}
		proteinTotalIntenMap.put(protein, total);
	}

	private static int resolveFileIndex(String fileName, Map<String, Integer> fileIdMap, List<String> fileList,
			FileIndexResolver resolver) {
		if (resolver != null) {
			return resolver.resolve(fileName, fileIdMap, fileList);
		}
		if (fileName == null || fileName.isBlank()) {
			return -1;
		}
		String key = fileName.trim();
		Integer index = fileIdMap.get(key);
		if (index != null) {
			return index;
		}
		String base = baseNameNoExt(key);
		for (Map.Entry<String, Integer> entry : fileIdMap.entrySet()) {
			String candidate = entry.getKey();
			if (candidate == null) {
				continue;
			}
			if (base.equals(baseNameNoExt(candidate.trim()))) {
				return entry.getValue();
			}
		}
		for (int i = 0; i < fileList.size(); i++) {
			String candidate = fileList.get(i);
			if (candidate == null) {
				continue;
			}
			if (key.equals(candidate.trim()) || base.equals(baseNameNoExt(candidate.trim()))) {
				return i;
			}
		}
		return -1;
	}

	private static String baseNameNoExt(String path) {
		String name = new java.io.File(path).getName();
		int dot = name.lastIndexOf('.');
		return dot > 0 ? name.substring(0, dot) : name;
	}

	private static boolean isValidFileId(int fileId, int fileCount) {
		return fileId >= 0 && fileId < fileCount;
	}

	private static int sum(int[] values) {
		int total = 0;
		if (values == null) {
			return total;
		}
		for (int value : values) {
			total += value;
		}
		return total;
	}

	private static boolean isDecoyProtein(String protein, String decoySymbol) {
		return protein == null || (decoySymbol != null && !decoySymbol.isEmpty() && protein.startsWith(decoySymbol));
	}

	private static boolean mapsOnlyToDecoys(MetaPeptide peptide, String decoySymbol) {
		return mapsOnlyToDecoys(peptide.getProteins(), decoySymbol);
	}

	private static boolean mapsOnlyToDecoys(String[] proteins, String decoySymbol) {
		if (decoySymbol == null || decoySymbol.isEmpty()) {
			return false;
		}
		if (proteins == null || proteins.length == 0) {
			return false;
		}
		for (String protein : proteins) {
			if (!isDecoyProtein(protein, decoySymbol)) {
				return false;
			}
		}
		return true;
	}

	private static boolean isMagProtein(String protein, String magIdentifier) {
		return protein != null && magIdentifier != null && !magIdentifier.isEmpty() && protein.startsWith(magIdentifier);
	}

	private static String genomeOf(String protein) {
		int idx = protein.indexOf("_");
		return idx > 0 ? protein.substring(0, idx) : protein;
	}

	private static boolean hasIorL(String sequence) {
		return sequence != null && (sequence.contains("I") || sequence.contains("L"));
	}
}
