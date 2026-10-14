package uk.ac.quadram.metapilot.core.mag;

import org.apache.commons.math3.stat.inference.TTest;
import uk.ac.quadram.metapilot.core.math.MathTool;
import uk.ac.quadram.metapilot.task.io.pep.MetaPSM;
import uk.ac.quadram.metapilot.task.io.pep.MetaPeptide;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * DDA MAG genome/protein evidence scoring.
 */
public final class DdaMagGenomeScoringService {

	public static final int IDX_S = 0;
	public static final int IDX_P_EMP = 1;
	public static final int IDX_Q_VALUE = 2;

	private static final boolean DDA_GENOME_SCORE_INTENSITY_ONLY = false;

	private DdaMagGenomeScoringService() {
	}

	/**
	 * True when every protein of a PSM is a decoy. A target sequence that also
	 * occurs in a reversed protein (short or near-palindromic tryptic peptides
	 * such as ARIPLAK) is target evidence with a target-level score; counting it
	 * as a decoy puts target scores into the resampling null. On PXD050781 that
	 * contamination (553 such peptides, 68 at target-level weight) moved genome
	 * p_emp from 0.004 to 0.18 and dropped 18 genomes from the MSFragger + Sage run.
	 */
	static boolean isDecoyOnly(String[] proteins, String decoyPrefix) {
		if (proteins == null || decoyPrefix == null || decoyPrefix.isEmpty()) {
			return false;
		}
		boolean anyDecoy = false;
		for (String protein : proteins) {
			if (protein == null || protein.isEmpty()) {
				continue;
			}
			if (!protein.startsWith(decoyPrefix)) {
				return false;
			}
			anyDecoy = true;
		}
		return anyDecoy;
	}

	/** Scoring unit: the base peptide with I and L merged, so charge states and modified forms count once. */
	static String peptideKey(MetaPSM psm) {
		String seq = psm.getBaseSequence();
		if (seq == null || seq.isEmpty()) {
			seq = stripModifications(psm.getModSeq());
		}
		return seq.replace('I', 'L');
	}

	private static String stripModifications(String modseq) {
		StringBuilder sb = new StringBuilder(modseq.length());
		int depth = 0;
		for (int i = 0; i < modseq.length(); i++) {
			char c = modseq.charAt(i);
			if (c == '[' || c == '(') {
				depth++;
			} else if (c == ']' || c == ')') {
				depth = Math.max(0, depth - 1);
			} else if (depth == 0 && c >= 'A' && c <= 'Z') {
				sb.append(c);
			}
		}
		return sb.toString();
	}

	private static String uniqueDecoyGenome(Set<String> proteins, String revIdentifier, String magIdentifier) {
		String genome = null;
		for (String pro : proteins) {
			if (!pro.startsWith(revIdentifier)) {
				continue;
			}
			String g = (magIdentifier + pro.substring(revIdentifier.length())).split("_")[0];
			if (genome == null) {
				genome = g;
			} else if (!genome.equals(g)) {
				return null;
			}
		}
		return genome;
	}

	public static final class Result {
		public final HashMap<String, Double> genomeScoreMap;
		public final HashMap<String, double[]> diaStyleGenomeStats;
		public final DdaGenomeScoreRuninfo diaStyleRuninfo;
		public final int targetPeptideCount;
		public final int decoyPeptideCount;
		public final int genomeCountBeforeFiltering;
		public final int razorGenomeCountBeforeFiltering;
		public final int decoyGenomeCount;
		public final int genomeCountAfterFiltering;
		public final double[] refinedDecoyGenomeLogIntensities;
		public final boolean diaStyle;

		private Result(HashMap<String, Double> genomeScoreMap, HashMap<String, double[]> diaStyleGenomeStats,
				DdaGenomeScoreRuninfo diaStyleRuninfo, int targetPeptideCount, int decoyPeptideCount,
				int genomeCountBeforeFiltering, int razorGenomeCountBeforeFiltering, int decoyGenomeCount,
				int genomeCountAfterFiltering, double[] refinedDecoyGenomeLogIntensities, boolean diaStyle) {
			this.genomeScoreMap = genomeScoreMap;
			this.diaStyleGenomeStats = diaStyleGenomeStats;
			this.diaStyleRuninfo = diaStyleRuninfo;
			this.targetPeptideCount = targetPeptideCount;
			this.decoyPeptideCount = decoyPeptideCount;
			this.genomeCountBeforeFiltering = genomeCountBeforeFiltering;
			this.razorGenomeCountBeforeFiltering = razorGenomeCountBeforeFiltering;
			this.decoyGenomeCount = decoyGenomeCount;
			this.genomeCountAfterFiltering = genomeCountAfterFiltering;
			this.refinedDecoyGenomeLogIntensities = refinedDecoyGenomeLogIntensities;
			this.diaStyle = diaStyle;
		}

		/**
		 * Result for input on which no genome scoring can run (no PSMs, scores or
		 * decoys - e.g. a plain peptide intensity table). The caller has already
		 * decided the accepted genome set; every entry of the map is accepted and no
		 * statistics (S / p_emp / q_value) are attached.
		 */
		public static Result acceptedGenomesOnly(HashMap<String, Double> acceptedGenomeScoreMap) {
			HashMap<String, Double> map = acceptedGenomeScoreMap != null ? acceptedGenomeScoreMap : new HashMap<>();
			return new Result(map, null, null, 0, 0, 0, 0, 0, map.size(), new double[0], false);
		}

		/** Result carrying {@link GenomeEvidenceTest} statistics ({@code {kept peptides, p, q}} per genome). */
		public static Result fromGenomeEvidence(HashMap<String, Double> genomeScoreMap,
				HashMap<String, double[]> genomeStats) {
			int accepted = DdaMagProteinInferenceService.confidentGenomes(genomeScoreMap, genomeStats,
					DdaMagProteinInferenceService.DEFAULT_GENOME_Q_VALUE_THRESHOLD).size();
			return new Result(genomeScoreMap, genomeStats, null, 0, 0, genomeStats.size(), 0, 0, accepted,
					new double[0], true);
		}
	}

	public static Result refineGenomes(MetaPeptide[] quanPeps, String decoySymbol, MetaPSM[] psms,
			String magIdentifier) {
		if (psms != null && psms.length > 0) {
			return computeGenomeScoreMapDiaStyle(quanPeps, psms, decoySymbol, magIdentifier);
		}
		return refineGenomesLegacy(quanPeps, decoySymbol, magIdentifier);
	}

	public static HashMap<String, double[]> computeProteinScoreMap(MetaPeptide[] quanPeps, MetaPSM[] psms,
			String decoySymbol) {
		HashMap<String, double[]> proteinScoreMap = new HashMap<>();
		if (quanPeps == null || psms == null || psms.length == 0) {
			return proteinScoreMap;
		}
		String decoyPrefix = decoySymbol != null ? decoySymbol : "";

		HashMap<String, MetaPeptide> modSeqToPeptide = new HashMap<>();
		for (MetaPeptide peptide : quanPeps) {
			if (peptide != null && peptide.getModSeq() != null) {
				modSeqToPeptide.put(peptide.getModSeq(), peptide);
			}
		}

		HashMap<String, DdaPeptideEvidence> peptideMap = new HashMap<>();
		for (MetaPSM psm : psms) {
			if (psm == null) {
				continue;
			}
			String modseq = psm.getModSeq();
			if (modseq == null || modseq.isEmpty()) {
				continue;
			}
			String[] proteins = psm.getProteins();
			boolean isDecoy = isDecoyOnly(proteins, decoyPrefix);
			String key = modseq + "|" + psm.getCharge() + "|" + isDecoy;
			DdaPeptideEvidence ev = peptideMap.get(key);
			if (ev == null) {
				ev = new DdaPeptideEvidence();
				ev.bestScore = psm.getScore();
				ev.isDecoy = isDecoy;
				MetaPeptide peptide = modSeqToPeptide.get(modseq);
				if (peptide != null && peptide.getIntensity() != null) {
					ev.totalIntensity = MathTool.getTotal(peptide.getIntensity());
				}
				peptideMap.put(key, ev);
			} else {
				ev.bestScore = Math.max(ev.bestScore, psm.getScore());
			}
			for (String protein : proteins != null ? proteins : new String[0]) {
				if (protein != null && !protein.isEmpty()) {
					ev.proteins.add(protein);
				}
			}
		}

		double minLog = Double.MAX_VALUE;
		double maxLog = -Double.MAX_VALUE;
		for (DdaPeptideEvidence ev : peptideMap.values()) {
			double logI = Math.log(1.0 + Math.max(0.0, ev.totalIntensity)) / Math.log(2.0);
			minLog = Math.min(minLog, logI);
			maxLog = Math.max(maxLog, logI);
		}
		if (peptideMap.isEmpty()) {
			return proteinScoreMap;
		}
		double logRange = Math.max(1e-10, maxLog - minLog);

		final double cap = 50.0;
		ArrayList<Double> decoyEPool = new ArrayList<>();
		HashMap<String, ArrayList<Double>> proteinPeptides = new HashMap<>();
		for (DdaPeptideEvidence ev : peptideMap.values()) {
			double e = Math.min(cap, Math.max(0.0, ev.bestScore));
			double logI = Math.log(1.0 + Math.max(0.0, ev.totalIntensity)) / Math.log(2.0);
			double norm = (logI - minLog) / logRange;
			double wei = e * (0.5 + 0.5 * norm);

			if (ev.isDecoy) {
				decoyEPool.add(wei);
				continue;
			}
			for (String protein : ev.proteins) {
				if (protein == null || (!decoyPrefix.isEmpty() && protein.startsWith(decoyPrefix))) {
					continue;
				}
				proteinPeptides.computeIfAbsent(protein, k -> new ArrayList<>()).add(wei);
			}
		}

		HashMap<String, Double> SMap = new HashMap<>();
		HashMap<String, Integer> proteinNMap = new HashMap<>();
		for (Map.Entry<String, ArrayList<Double>> entry : proteinPeptides.entrySet()) {
			double S = 0.0;
			for (double wei : entry.getValue()) {
				S += wei;
			}
			SMap.put(entry.getKey(), S);
			proteinNMap.put(entry.getKey(), entry.getValue().size());
		}

		if (decoyEPool.isEmpty()) {
			decoyEPool.add(0.0);
		}
		double[] decoyArr = new double[decoyEPool.size()];
		for (int i = 0; i < decoyArr.length; i++) {
			decoyArr[i] = decoyEPool.get(i);
		}
		Random rnd = new Random(42);
		int B = 10000;
		HashMap<String, Double> pempMap = new HashMap<>();
		for (String protein : SMap.keySet()) {
			double Sp = SMap.get(protein);
			int M = proteinNMap.get(protein);
			int count = 0;
			for (int b = 0; b < B; b++) {
				double Sb = 0.0;
				for (int i = 0; i < M; i++) {
					Sb += decoyArr[rnd.nextInt(decoyArr.length)];
				}
				if (Sb >= Sp) {
					count++;
				}
			}
			pempMap.put(protein, (1.0 + count) / (1.0 + B));
		}

		HashMap<String, Double> qValueMap = computeBenjaminiHochbergQValues(pempMap);
		for (String protein : pempMap.keySet()) {
			double pemp = pempMap.get(protein);
			proteinScoreMap.put(protein, new double[] { SMap.get(protein), pemp, qValueMap.get(protein) });
		}
		return proteinScoreMap;
	}

	/**
	 * @deprecated Legacy DDA intensity/t-test fallback. Active DDA MAG workflows
	 *             should pass PSMs and use the PSM-score/top-k path.
	 */
	@Deprecated
	private static Result refineGenomesLegacy(MetaPeptide[] quanPeps, String decoySymbol, String magIdentifier) {
		HashMap<String, Double> empty = new HashMap<>();
		if (quanPeps == null || magIdentifier == null || magIdentifier.isEmpty()) {
			return new Result(empty, null, null, 0, 0, 0, 0, 0, 0, new double[0], false);
		}

		String revIdentifier = decoySymbol + magIdentifier;
		HashMap<String, Double> genomeIntensityScoreMap = new HashMap<>();
		double[] peptideIntensityScore = new double[quanPeps.length];
		int[] pepCount = new int[2];

		for (int i = 0; i < quanPeps.length; i++) {
			MetaPeptide quanPep = quanPeps[i];
			if (quanPep == null) {
				continue;
			}
			double[] pepIntensity = quanPep.getIntensity();
			double totalIntensity = MathTool.getTotal(pepIntensity);
			int[] idenType = quanPep.getIdenType();
			int idenCount = 0;
			for (int j = 0; idenType != null && j < idenType.length; j++) {
				if (idenType[j] == MetaPeptide.IDEN_MSMS) {
					idenCount++;
				}
			}

			if (totalIntensity == 0 || idenCount == 0) {
				continue;
			}

			double[] indiPEP = quanPep.getIndividualPEP();
			for (int j = 0; pepIntensity != null && j < pepIntensity.length; j++) {
				if (idenType[j] == MetaPeptide.IDEN_MSMS) {
					if (indiPEP[j] == 0) {
						peptideIntensityScore[i] += pepIntensity[j] * 1.0E7;
					} else {
						peptideIntensityScore[i] += pepIntensity[j] / indiPEP[j];
					}
				}
			}

			String[] proteins = quanPep.getProteins();
			boolean isDecoy = false;
			for (int j = 0; proteins != null && j < proteins.length; j++) {
				if (proteins[j].startsWith(magIdentifier)) {
					String genome = proteins[j].split("_")[0];
					genomeIntensityScoreMap.put(genome,
							genomeIntensityScoreMap.getOrDefault(genome, 0.0) + peptideIntensityScore[i]);
				} else if (proteins[j].startsWith(revIdentifier)) {
					isDecoy = true;
					String pro = magIdentifier + proteins[j].substring(revIdentifier.length());
					String genome = decoySymbol + pro.split("_")[0];
					genomeIntensityScoreMap.put(genome,
							genomeIntensityScoreMap.getOrDefault(genome, 0.0) + peptideIntensityScore[i]);
				}
			}
			if (isDecoy) {
				pepCount[1]++;
			} else {
				pepCount[0]++;
			}
		}

		HashMap<String, Double> razorGenomeIntenMap = new HashMap<>();
		for (int i = 0; i < quanPeps.length; i++) {
			if (peptideIntensityScore[i] > 0 && quanPeps[i] != null) {
				String[] proteins = quanPeps[i].getProteins();
				String razorGenome = "";
				double razorGenomeIntensity = 0;
				for (int j = 0; proteins != null && j < proteins.length; j++) {
					if (proteins[j].startsWith(magIdentifier)) {
						String genome = proteins[j].split("_")[0];
						if (genomeIntensityScoreMap.containsKey(genome)) {
							double genomeIntensity = genomeIntensityScoreMap.get(genome);
							if (genomeIntensity > razorGenomeIntensity) {
								razorGenome = genome;
								razorGenomeIntensity = genomeIntensity;
							}
						}
					} else if (proteins[j].startsWith(revIdentifier)) {
						String pro = magIdentifier + proteins[j].substring(revIdentifier.length());
						String genome = decoySymbol + pro.split("_")[0];
						if (genomeIntensityScoreMap.containsKey(genome)) {
							double genomeIntensity = genomeIntensityScoreMap.get(genome);
							if (genomeIntensity > razorGenomeIntensity) {
								razorGenome = genome;
								razorGenomeIntensity = genomeIntensity;
							}
						}
					}
				}

				if (razorGenome.startsWith(magIdentifier) || razorGenome.startsWith(revIdentifier)) {
					razorGenomeIntenMap.put(razorGenome, 0.0);
				}
			}
		}

		for (int i = 0; i < quanPeps.length; i++) {
			if (peptideIntensityScore[i] > 0 && quanPeps[i] != null) {
				String[] proteins = quanPeps[i].getProteins();
				int count = 0;
				for (int j = 0; proteins != null && j < proteins.length; j++) {
					if (proteins[j].startsWith(magIdentifier)) {
						String genome = proteins[j].split("_")[0];
						if (razorGenomeIntenMap.containsKey(genome)) {
							count++;
						}
					} else if (proteins[j].startsWith(revIdentifier)) {
						String pro = magIdentifier + proteins[j].substring(revIdentifier.length());
						String genome = decoySymbol + pro.split("_")[0];
						if (razorGenomeIntenMap.containsKey(genome)) {
							count++;
						}
					}
				}

				if (count > 0) {
					for (int j = 0; proteins != null && j < proteins.length; j++) {
						if (proteins[j].startsWith(magIdentifier)) {
							String genome = proteins[j].split("_")[0];
							if (razorGenomeIntenMap.containsKey(genome)) {
								razorGenomeIntenMap.put(genome,
										razorGenomeIntenMap.get(genome) + peptideIntensityScore[i]);
							}
						} else if (proteins[j].startsWith(revIdentifier)) {
							String pro = magIdentifier + proteins[j].substring(revIdentifier.length());
							String genome = decoySymbol + pro.split("_")[0];
							if (razorGenomeIntenMap.containsKey(genome)) {
								razorGenomeIntenMap.put(genome,
										razorGenomeIntenMap.get(genome) + peptideIntensityScore[i]);
							}
						}
					}
				}
			}
		}

		List<Map.Entry<String, Double>> sortedRazorGenomeList = razorGenomeIntenMap.entrySet().stream()
				.sorted(Map.Entry.comparingByValue(Collections.reverseOrder())).collect(Collectors.toList());
		String[] razorGenomes = sortedRazorGenomeList.stream().map(Map.Entry::getKey).toArray(String[]::new);

		ArrayList<Double> revLogIntenList = new ArrayList<>();
		for (String razorGenome : razorGenomes) {
			if (razorGenome.startsWith(revIdentifier)) {
				double intensity = razorGenomeIntenMap.get(razorGenome);
				if (intensity > 0) {
					revLogIntenList.add(Math.log10(intensity));
				}
			}
		}

		HashMap<String, Double> genomeScoreMap = new HashMap<>();
		double[] revLogIntensities = new double[0];
		if (revLogIntenList.size() == 0) {
			for (String razorGenome : razorGenomes) {
				double intensity = razorGenomeIntenMap.get(razorGenome);
				if (intensity == 0) {
					break;
				}
				genomeScoreMap.put(razorGenome, intensity);
			}
		} else {
			TTest tTest = new TTest();
			revLogIntensities = new double[revLogIntenList.size()];
			for (int i = 2; i < revLogIntenList.size(); i++) {
				revLogIntensities = new double[i];
				for (int j = 0; j < revLogIntensities.length; j++) {
					revLogIntensities[j] = revLogIntenList.get(j);
				}

				double pValue = tTest.tTest(revLogIntensities[0], revLogIntensities);
				if (pValue < 0.01) {
					break;
				}
			}

			if (revLogIntensities.length == 0) {
				for (String razorGenome : razorGenomes) {
					double intensity = razorGenomeIntenMap.get(razorGenome);
					if (razorGenome.startsWith(revIdentifier) || intensity == 0) {
						break;
					}
					genomeScoreMap.put(razorGenome, 0.0);
				}
			} else if (revLogIntensities.length < 6) {
				double[] extended = new double[6];
				for (int i = 0; i < extended.length; i++) {
					extended[i] = revLogIntensities[i % revLogIntensities.length];
				}
				for (String razorGenome : razorGenomes) {
					double intensity = razorGenomeIntenMap.get(razorGenome);
					if (razorGenome.startsWith(revIdentifier) || intensity == 0) {
						break;
					}
					double pvalue = tTest.tTest(Math.log10(intensity), extended);
					if (pvalue < 0.01) {
						genomeScoreMap.put(razorGenome, pvalue);
					}
				}
			} else {
				for (String razorGenome : razorGenomes) {
					double intensity = razorGenomeIntenMap.get(razorGenome);
					if (razorGenome.startsWith(revIdentifier) || intensity == 0) {
						break;
					}
					double pvalue = tTest.tTest(Math.log10(intensity), revLogIntensities);
					if (pvalue < 0.01) {
						genomeScoreMap.put(razorGenome, pvalue);
					}
				}
			}
		}

		return new Result(genomeScoreMap, null, null, pepCount[0], pepCount[1], genomeIntensityScoreMap.size(),
				razorGenomeIntenMap.size(), revLogIntenList.size(), genomeScoreMap.size(), revLogIntensities, false);
	}

	private static Result computeGenomeScoreMapDiaStyle(MetaPeptide[] quanPeps, MetaPSM[] psms, String decoySymbol,
			String magIdentifier) {
		HashMap<String, Double> result = new HashMap<>();
		if (quanPeps == null || psms == null || magIdentifier == null || magIdentifier.isEmpty()) {
			return new Result(result, null, null, 0, 0, 0, 0, 0, 0, new double[0], true);
		}
		String revIdentifier = decoySymbol + magIdentifier;

		HashMap<String, MetaPeptide> modSeqToPeptide = new HashMap<>();
		for (MetaPeptide p : quanPeps) {
			if (p != null && p.getModSeq() != null) {
				modSeqToPeptide.put(p.getModSeq(), p);
			}
		}

		// Intensity is per modified form; everything else is per peptide unit
		// (peptideKey), which takes its strongest form's weight.
		HashMap<String, DdaPeptideEvidence> formMap = new HashMap<>();
		HashMap<String, DdaPeptideEvidence> unitMap = new HashMap<>();
		for (MetaPSM psm : psms) {
			if (psm == null) continue;
			String modseq = psm.getModSeq();
			if (modseq == null || modseq.isEmpty()) continue;
			double score = psm.getScore();
			String[] pros = psm.getProteins();
			boolean isDecoy = isDecoyOnly(pros, revIdentifier);
			String unitKey = peptideKey(psm) + "|" + isDecoy;
			String formKey = modseq + "|" + isDecoy;
			DdaPeptideEvidence form = formMap.get(formKey);
			if (form == null) {
				form = new DdaPeptideEvidence();
				form.bestScore = score;
				form.isDecoy = isDecoy;
				form.unitKey = unitKey;
				MetaPeptide pep = modSeqToPeptide.get(modseq);
				if (pep != null && pep.getIntensity() != null) {
					form.totalIntensity = MathTool.getTotal(pep.getIntensity());
				}
				formMap.put(formKey, form);
			} else {
				form.bestScore = Math.max(form.bestScore, score);
			}
			DdaPeptideEvidence unit = unitMap.get(unitKey);
			if (unit == null) {
				unit = new DdaPeptideEvidence();
				unit.bestScore = score;
				unit.isDecoy = isDecoy;
				unitMap.put(unitKey, unit);
			} else {
				unit.bestScore = Math.max(unit.bestScore, score);
			}
			if (pros != null) {
				for (String pro : pros) {
					if (pro != null && !pro.isEmpty()) unit.proteins.add(pro);
				}
			}
		}

		double minLog = Double.MAX_VALUE, maxLog = -Double.MAX_VALUE;
		for (DdaPeptideEvidence ev : formMap.values()) {
			double logI = Math.log(1.0 + Math.max(0, ev.totalIntensity)) / Math.log(2.0);
			minLog = Math.min(minLog, logI);
			maxLog = Math.max(maxLog, logI);
		}
		double logRange = Math.max(1e-10, maxLog - minLog);
		final double cap = 50.0;
		for (DdaPeptideEvidence form : formMap.values()) {
			double eRaw = Math.min(cap, Math.max(0, form.bestScore));
			double e = DDA_GENOME_SCORE_INTENSITY_ONLY ? 1.0 : eRaw;
			double logI = Math.log(1.0 + Math.max(0, form.totalIntensity)) / Math.log(2.0);
			double intensityFactor = 0.5 + 0.5 * (logI - minLog) / logRange;
			DdaPeptideEvidence unit = unitMap.get(form.unitKey);
			if (e * intensityFactor > unit.weight) {
				unit.weight = e * intensityFactor;
				unit.eRaw = eRaw;
				unit.e = e;
				unit.intensityFactor = intensityFactor;
			}
		}
		ArrayList<Double> decoyEPool = new ArrayList<>();
		HashMap<String, DdaGenomePeptideEntry> genomePeptides = new HashMap<>();
		int nTargetPeptides = 0;
		int nDecoyPeptides = 0;
		int nUniquePeptides = 0;
		int nSharedPeptides = 0;
		HashMap<Integer, Integer> mDistribution = new HashMap<>();
		ArrayList<Double> targetERaw = new ArrayList<>();
		ArrayList<Double> decoyERaw = new ArrayList<>();

		HashMap<String, ArrayList<Double>> decoyGenomeWeights = new HashMap<>();
		for (Map.Entry<String, DdaPeptideEvidence> ent : unitMap.entrySet()) {
			DdaPeptideEvidence ev = ent.getValue();
			double eRaw = ev.eRaw;
			double e = ev.e;
			double intensityFactor = ev.intensityFactor;
			String peptide = ent.getKey().split("\\|", 2)[0];

			if (ev.isDecoy) {
				nDecoyPeptides++;
				decoyERaw.add(eRaw);
				decoyEPool.add(ev.weight);
				String decoyGenome = uniqueDecoyGenome(ev.proteins, revIdentifier, magIdentifier);
				if (decoyGenome != null) {
					decoyGenomeWeights.computeIfAbsent(decoyGenome, g -> new ArrayList<>()).add(ev.weight);
				}
				continue;
			}
			targetERaw.add(eRaw);

			HashSet<String> genomes = new HashSet<>();
			for (String pro : ev.proteins) {
				if (pro.startsWith(revIdentifier)) continue;
				if (pro.startsWith(magIdentifier)) {
					String[] parts = pro.split("_");
					if (parts.length > 0 && !parts[0].isEmpty()) genomes.add(parts[0]);
				}
			}
			if (genomes.isEmpty()) continue;
			int m = genomes.size();
			boolean uniquePeptidesOnly = true;
			if (m != 1 && uniquePeptidesOnly) {
				nSharedPeptides++;
				mDistribution.merge(m, 1, Integer::sum);
				continue;
			}
			nTargetPeptides++;
			nUniquePeptides++;
			mDistribution.merge(1, 1, Integer::sum);
			double w = (m == 1) ? 1.0 : (1.0 / m);
			double wei = w * e * intensityFactor;

			for (String g : genomes) {
				DdaWeightedScore wps = new DdaWeightedScore();
				wps.seq = peptide;
				wps.e = e;
				wps.w = w;
				wps.intensityFactor = intensityFactor;
				wps.we = wei;
				DdaGenomePeptideEntry ge = genomePeptides.get(g);
				if (ge == null) {
					ge = new DdaGenomePeptideEntry();
					genomePeptides.put(g, ge);
				}
				ge.weightedScores.add(wps);
			}
		}

		int k = 5;
		Random rnd = new Random(42);
		int B = 10000;
		HashMap<String, Double> SMap = new HashMap<>();
		HashMap<String, Integer> genomeNMap = new HashMap<>();
		HashMap<String, Integer> genomeTotalPeptidesMap = new HashMap<>();
		for (String g : genomePeptides.keySet()) {
			List<DdaWeightedScore> scores = genomePeptides.get(g).weightedScores;
			genomeTotalPeptidesMap.put(g, scores.size());
			scores.sort((a, b) -> Double.compare(b.we, a.we));
			int n = Math.min(k, scores.size());
			double S = 0;
			for (int i = 0; i < n; i++) S += scores.get(i).we;
			SMap.put(g, S);
			genomeNMap.put(g, n);
		}

		if (decoyEPool.isEmpty()) decoyEPool.add(0.0);
		double[] decoyArr = new double[decoyEPool.size()];
		for (int i = 0; i < decoyArr.length; i++) decoyArr[i] = decoyEPool.get(i);
		int decoyLen = decoyArr.length;
		int maxM = 0;
		for (Integer m : genomeTotalPeptidesMap.values()) maxM = Math.max(maxM, m);
		double[] tmp = new double[Math.max(maxM, 1)];
		HashMap<String, Double> pempMap = new HashMap<>();
		for (String g : SMap.keySet()) {
			double Sg = SMap.get(g);
			int M = genomeTotalPeptidesMap.get(g);
			int ng = genomeNMap.get(g);
			int count = 0;
			for (int b = 0; b < B; b++) {
				for (int i = 0; i < M; i++) tmp[i] = decoyArr[rnd.nextInt(decoyLen)];
				Arrays.sort(tmp, 0, M);
				double Sb = 0;
				for (int i = M - 1; i >= M - ng && i >= 0; i--) Sb += tmp[i];
				if (Sb >= Sg) count++;
			}
			double pemp = (1.0 + count) / (1.0 + B);
			pempMap.put(g, pemp);
		}

		// Reversed proteomes through the same statistic, for calibration only.
		HashMap<String, Double> decoyGenomeScorePMap = new HashMap<>();
		HashMap<String, Integer> decoyGenomeScoreNMap = new HashMap<>();
		for (Map.Entry<String, ArrayList<Double>> ent : decoyGenomeWeights.entrySet()) {
			ArrayList<Double> weights = ent.getValue();
			weights.sort(Collections.reverseOrder());
			int M = weights.size();
			int ng = Math.min(k, M);
			double Sg = 0;
			for (int i = 0; i < ng; i++) Sg += weights.get(i);
			double[] draws = new double[M];
			int count = 0;
			for (int b = 0; b < B; b++) {
				for (int i = 0; i < M; i++) draws[i] = decoyArr[rnd.nextInt(decoyLen)];
				Arrays.sort(draws);
				double Sb = 0;
				for (int i = M - 1; i >= M - ng; i--) Sb += draws[i];
				if (Sb >= Sg) count++;
			}
			decoyGenomeScorePMap.put(ent.getKey(), (1.0 + count) / (1.0 + B));
			decoyGenomeScoreNMap.put(ent.getKey(), M);
		}

		String[] genomes = pempMap.keySet().toArray(new String[0]);
		Arrays.sort(genomes, Comparator.comparingDouble(pempMap::get));
		int nGenomes = genomes.length;
		HashMap<String, Double> qValueMap = computeBenjaminiHochbergQValues(pempMap);
		HashMap<String, double[]> statsForTaxonomy = new HashMap<>();
		for (int i = 0; i < nGenomes; i++) {
			String g = genomes[i];
			double pemp = pempMap.get(g);
			double q = qValueMap.get(g);
			result.put(g, 1.0 - q);
			double S = SMap.get(g);
			statsForTaxonomy.put(g, new double[] { S, pemp, q });
		}

		List<String> sortedByScore = new ArrayList<>(result.keySet());
		sortedByScore.sort(Comparator.comparingDouble(g -> -SMap.get(g)));
		double[] pempArr = new double[sortedByScore.size()];
		double[] qvals = new double[sortedByScore.size()];
		for (int i = 0; i < sortedByScore.size(); i++) {
			String g = sortedByScore.get(i);
			pempArr[i] = pempMap.get(g);
			qvals[i] = statsForTaxonomy.get(g)[IDX_Q_VALUE];
		}
		int decoyPoolSize = decoyEPool.size();
		double decoyEMin = Double.MAX_VALUE, decoyEMax = -Double.MAX_VALUE, decoyESum = 0;
		for (Double de : decoyEPool) {
			decoyEMin = Math.min(decoyEMin, de);
			decoyEMax = Math.max(decoyEMax, de);
			decoyESum += de;
		}
		double decoyEMean = decoyPoolSize > 0 ? decoyESum / decoyPoolSize : 0;
		HashMap<String, List<DdaWeightedScore>> genomeTopKMap = new HashMap<>();
		for (String g : genomePeptides.keySet()) {
			List<DdaWeightedScore> scores = genomePeptides.get(g).weightedScores;
			int n = Math.min(k, scores.size());
			genomeTopKMap.put(g, new ArrayList<>(scores.subList(0, n)));
		}

		DdaGenomeScoreRuninfo runinfo = new DdaGenomeScoreRuninfo();
		runinfo.decoyPoolSize = decoyPoolSize;
		runinfo.nDecoyPeptides = nDecoyPeptides;
		runinfo.nTargetPeptides = nTargetPeptides;
		runinfo.nUniquePeptides = nUniquePeptides;
		runinfo.nSharedPeptides = nSharedPeptides;
		runinfo.mDistribution = mDistribution;
		runinfo.decoyEMin = decoyEMin;
		runinfo.decoyEMax = decoyEMax;
		runinfo.decoyEMean = decoyEMean;
		runinfo.genomeTotalPeptidesMap = genomeTotalPeptidesMap;
		runinfo.genomeNMap = genomeNMap;
		runinfo.sortedByScore = sortedByScore;
		runinfo.pempArr = pempArr;
		runinfo.qvals = qvals;
		runinfo.genomeTopKMap = genomeTopKMap;
		runinfo.SMap = SMap;
		runinfo.minLog = minLog;
		runinfo.maxLog = maxLog;
		runinfo.logRange = logRange;
		runinfo.k = k;
		runinfo.B = B;
		runinfo.cap = cap;
		runinfo.targetERaw = targetERaw;
		runinfo.decoyERaw = decoyERaw;
		runinfo.useIntensityOnly = DDA_GENOME_SCORE_INTENSITY_ONLY;
		runinfo.decoyGenomeScorePMap = decoyGenomeScorePMap;
		runinfo.decoyGenomeScoreNMap = decoyGenomeScoreNMap;

		return new Result(result, statsForTaxonomy, runinfo, nTargetPeptides, nDecoyPeptides, 0, 0, decoyPoolSize,
				result.size(), new double[0], true);
	}

	private static HashMap<String, Double> computeBenjaminiHochbergQValues(HashMap<String, Double> pValueMap) {
		HashMap<String, Double> qValueMap = new HashMap<>();
		if (pValueMap == null || pValueMap.isEmpty()) {
			return qValueMap;
		}
		String[] ids = pValueMap.keySet().toArray(new String[0]);
		Arrays.sort(ids, Comparator.comparingDouble(pValueMap::get));
		int n = ids.length;
		double minAdjusted = 1.0;
		for (int i = n - 1; i >= 0; i--) {
			double p = pValueMap.get(ids[i]);
			double adjusted = Math.min(1.0, p * n / (i + 1));
			minAdjusted = Math.min(minAdjusted, adjusted);
			qValueMap.put(ids[i], minAdjusted);
		}
		return qValueMap;
	}

	public static final class DdaGenomeScoreRuninfo {
		public int decoyPoolSize;
		public int nDecoyPeptides;
		public int nTargetPeptides;
		public int nUniquePeptides;
		public int nSharedPeptides;
		public HashMap<Integer, Integer> mDistribution;
		public double decoyEMin, decoyEMax, decoyEMean;
		public HashMap<String, Integer> genomeTotalPeptidesMap;
		public HashMap<String, Integer> genomeNMap;
		public List<String> sortedByScore;
		public double[] pempArr;
		public double[] qvals;
		public HashMap<String, List<DdaWeightedScore>> genomeTopKMap;
		public HashMap<String, Double> SMap;
		public double minLog, maxLog, logRange;
		public int k;
		public int B;
		public double cap;
		public ArrayList<Double> targetERaw;
		public ArrayList<Double> decoyERaw;
		public boolean useIntensityOnly;
		/** Reversed proteomes through the genome statistic: p-value and unique decoy peptide count. */
		public HashMap<String, Double> decoyGenomeScorePMap = new HashMap<>();
		public HashMap<String, Integer> decoyGenomeScoreNMap = new HashMap<>();
	}

	private static final class DdaPeptideEvidence {
		double bestScore;
		double totalIntensity;
		boolean isDecoy;
		HashSet<String> proteins = new HashSet<>();
		/** Genome scoring only: the peptide unit a modified form belongs to, and the unit's weight. */
		String unitKey;
		double weight = -1;
		double eRaw;
		double e;
		double intensityFactor;
	}

	public static final class DdaWeightedScore {
		public double we;
		public String seq;
		public double e;
		public double w;
		public double intensityFactor = 1.0;
	}

	private static final class DdaGenomePeptideEntry {
		ArrayList<DdaWeightedScore> weightedScores = new ArrayList<>();
	}
}
