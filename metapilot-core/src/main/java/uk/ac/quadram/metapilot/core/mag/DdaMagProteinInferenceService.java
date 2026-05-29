package uk.ac.quadram.metapilot.core.mag;

import uk.ac.quadram.metapilot.core.math.MathTool;
import uk.ac.quadram.metapilot.task.io.pep.MetaPeptide;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * DDA MAG razor protein inference from confident genomes.
 */
public final class DdaMagProteinInferenceService {

	public static final double DEFAULT_GENOME_Q_VALUE_THRESHOLD = 0.01;

	private DdaMagProteinInferenceService() {
	}

	public static final class Result {
		public final String[] razorProteins;
		public final HashMap<String, Double> proteinIntensityMap;
		public final int peptidesWithOnlyDecoyProteins;

		private Result(String[] razorProteins, HashMap<String, Double> proteinIntensityMap,
				int peptidesWithOnlyDecoyProteins) {
			this.razorProteins = razorProteins;
			this.proteinIntensityMap = proteinIntensityMap;
			this.peptidesWithOnlyDecoyProteins = peptidesWithOnlyDecoyProteins;
		}
	}

	public static Result refineProteins(MetaPeptide[] quanPeps, Map<String, Double> genomeScoreMap,
			String magIdentifier, String decoySymbol) {
		return refineProteins(quanPeps, genomeScoreMap, null, magIdentifier, decoySymbol,
				DEFAULT_GENOME_Q_VALUE_THRESHOLD);
	}

	public static Result refineProteins(MetaPeptide[] quanPeps, Map<String, Double> genomeScoreMap,
			Map<String, double[]> genomeStats, String magIdentifier, String decoySymbol) {
		return refineProteins(quanPeps, genomeScoreMap, genomeStats, magIdentifier, decoySymbol,
				DEFAULT_GENOME_Q_VALUE_THRESHOLD);
	}

	public static Result refineProteins(MetaPeptide[] quanPeps, Map<String, Double> genomeScoreMap,
			Map<String, double[]> genomeStats, String magIdentifier, String decoySymbol, double qValueThreshold) {
		if (quanPeps == null) {
			return new Result(new String[0], new HashMap<>(), 0);
		}
		if (genomeScoreMap == null) {
			genomeScoreMap = new HashMap<>();
		}
		if (magIdentifier == null) {
			magIdentifier = "";
		}

		Set<String> confidentGenomes = confidentGenomes(genomeScoreMap, genomeStats, qValueThreshold);
		HashMap<String, Double> proteinIntensityMap = new HashMap<>();
		HashMap<String, Double> genomeIntensityMap = new HashMap<>();
		for (MetaPeptide quanPep : quanPeps) {
			if (quanPep == null) {
				continue;
			}
			double[] pepIntensity = quanPep.getIntensity();
			String[] proteins = quanPep.getProteins();
			if (proteins == null) {
				continue;
			}

			int count = 0;
			for (String protein : proteins) {
				if (isDecoyProtein(protein, decoySymbol)) {
					continue;
				}
				if (isMagProtein(protein, magIdentifier)) {
					String genome = genomeOf(protein);
					if (confidentGenomes.contains(genome)) {
						count++;
					}
				} else {
					count++;
				}
			}

			if (count > 0) {
				double intensity = MathTool.getTotal(pepIntensity) / (double) count;
				for (String protein : proteins) {
					if (isDecoyProtein(protein, decoySymbol)) {
						continue;
					}
					if (isMagProtein(protein, magIdentifier)) {
						String genome = genomeOf(protein);
						if (confidentGenomes.contains(genome)) {
							addIntensity(proteinIntensityMap, protein, intensity);
							addIntensity(genomeIntensityMap, genome, intensity);
						}
					} else {
						addIntensity(proteinIntensityMap, protein, intensity);
					}
				}
			}
		}

		String[] razorProteins = new String[quanPeps.length];
		int peptidesWithOnlyDecoyProteins = countPeptidesWithOnlyDecoyProteins(quanPeps, decoySymbol);
		for (int i = 0; i < quanPeps.length; i++) {
			MetaPeptide quanPep = quanPeps[i];
			if (quanPep == null || quanPep.getProteins() == null) {
				continue;
			}
			String[] proteins = quanPep.getProteins();
			double proteinIntensity = 0;
			double genomeIntensity = 0;

			for (String protein : proteins) {
				if (isDecoyProtein(protein, decoySymbol)) {
					continue;
				}
				if (isMagProtein(protein, magIdentifier)) {
					String genome = genomeOf(protein);
					if (confidentGenomes.contains(genome)) {
						Double currentProteinIntensity = proteinIntensityMap.get(protein);
						if (currentProteinIntensity == null) {
							continue;
						}
						double currentGenomeIntensity = genomeIntensityMap.getOrDefault(genome, 0.0);
						if (isBetterRazorCandidate(currentProteinIntensity, currentGenomeIntensity, protein,
								proteinIntensity, genomeIntensity, razorProteins[i])) {
							razorProteins[i] = protein;
							proteinIntensity = currentProteinIntensity;
							genomeIntensity = currentGenomeIntensity;
						}
					}
				} else {
					Double currentProteinIntensity = proteinIntensityMap.get(protein);
					if (currentProteinIntensity != null && isBetterRazorCandidate(currentProteinIntensity, 0.0, protein,
							proteinIntensity, genomeIntensity, razorProteins[i])) {
						proteinIntensity = currentProteinIntensity;
						genomeIntensity = 0.0;
						razorProteins[i] = protein;
					}
				}
			}
		}

		return new Result(razorProteins, proteinIntensityMap, peptidesWithOnlyDecoyProteins);
	}

	private static void addIntensity(HashMap<String, Double> proteinIntensityMap, String protein, double intensity) {
		proteinIntensityMap.put(protein, proteinIntensityMap.getOrDefault(protein, 0.0) + intensity);
	}

	private static boolean isBetterRazorCandidate(double candidateProteinIntensity, double candidateGenomeIntensity,
			String candidateProtein, double currentProteinIntensity, double currentGenomeIntensity,
			String currentProtein) {
		if (candidateProteinIntensity > currentProteinIntensity) {
			return true;
		}
		if (candidateProteinIntensity < currentProteinIntensity) {
			return false;
		}
		if (candidateGenomeIntensity > currentGenomeIntensity) {
			return true;
		}
		if (candidateGenomeIntensity < currentGenomeIntensity) {
			return false;
		}
		return currentProtein == null || candidateProtein.compareTo(currentProtein) < 0;
	}

	public static Set<String> confidentGenomes(Map<String, Double> genomeScoreMap, Map<String, double[]> genomeStats,
			double qValueThreshold) {
		HashSet<String> confident = new HashSet<>();
		if (genomeStats != null && !genomeStats.isEmpty()) {
			for (Map.Entry<String, double[]> entry : genomeStats.entrySet()) {
				double[] stats = entry.getValue();
				if (stats != null && stats.length > DdaMagGenomeScoringService.IDX_Q_VALUE
						&& stats[DdaMagGenomeScoringService.IDX_Q_VALUE] <= qValueThreshold) {
					confident.add(entry.getKey());
				}
			}
			return confident;
		}
		if (genomeScoreMap != null) {
			confident.addAll(genomeScoreMap.keySet());
		}
		return confident;
	}

	private static int countPeptidesWithOnlyDecoyProteins(MetaPeptide[] quanPeps, String decoySymbol) {
		int count = 0;
		for (MetaPeptide quanPep : quanPeps) {
			if (quanPep == null || quanPep.getProteins() == null) {
				continue;
			}
			String[] proteins = quanPep.getProteins();
			boolean hasTarget = false;
			for (String protein : proteins) {
				if (!isDecoyProtein(protein, decoySymbol)) {
					hasTarget = true;
					break;
				}
			}
			if (!hasTarget && proteins.length > 0) {
				count++;
			}
		}
		return count;
	}

	private static boolean isDecoyProtein(String protein, String decoySymbol) {
		return protein == null || (decoySymbol != null && protein.startsWith(decoySymbol));
	}

	private static boolean isMagProtein(String protein, String magIdentifier) {
		return protein != null && magIdentifier != null && !magIdentifier.isEmpty() && protein.startsWith(magIdentifier);
	}

	private static String genomeOf(String protein) {
		int idx = protein.indexOf("_");
		return idx > 0 ? protein.substring(0, idx) : protein;
	}
}
