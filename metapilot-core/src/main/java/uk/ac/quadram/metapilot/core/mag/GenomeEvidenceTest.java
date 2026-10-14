package uk.ac.quadram.metapilot.core.mag;

import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.apache.commons.math3.special.Gamma;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collections;
import java.util.Comparator;
import java.util.EnumMap;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

/**
 * Genome acceptance by expected false counts, with near-neighbour exclusion and
 * a picked protein FDR within the accepted genomes
 * ({@code docs/developer/genome-catalogue-resolution-plan.md}, version 4, steps
 * 1, 2, 3 and 5).
 * <ol>
 * <li><b>Peptide FDR per class.</b> Unit: base peptide with I = L, scored by its
 * best PSM. Separate target-decoy walks for unmodified, known-shift and
 * unexplained-shift units; q &lt;= {@link Settings#peptideFdr}. Unmodified and
 * known shift are evidence; unexplained shift never is.</li>
 * <li><b>Expected false count.</b> A peptide is unique to genome g in run r when
 * every protein of r's database containing it (search digest) belongs to g; it
 * counts for g only when every run in which it was accepted calls it unique to
 * g. The decoy side uses the same rule with the reversed proteomes. {@code n_g}
 * = target evidence peptides unique to g; {@code F} = decoy evidence peptides
 * unique to one reversed genome; {@code pi_g = (w_g + 1) / (W + G)}, where
 * {@code w_g} counts {@code rev_g} among the {@code W} = top {@code K * F}
 * unique decoy peptides of <em>any</em> class by score (an evidence-class
 * window is too small to stay in the high-score range and reproduces the
 * all-score bias) and {@code G} is the number of candidate genomes. The added
 * one keeps a genome with no window decoy (typically one searched in few runs)
 * from passing on a single peptide. {@code p_g = P(Poisson(F * pi_g) >= n_g)}.</li>
 * <li><b>Sequential exclusion.</b> Genomes in ascending {@code log p}; before g
 * is tested, its peptides within {@link Settings#maxSubstitutions} of a peptide
 * of an accepted genome's full proteome are removed, and the same on the decoy
 * side (which scales {@code lambda_g}). A genome joins the accepted set when its
 * p-value is under the running Benjamini-Hochberg threshold; reported q-values
 * are Benjamini-Hochberg over the post-exclusion p-values. Every reversed genome
 * goes through the same calculation as a calibration row.</li>
 * <li>(Step 4, catalogue resolution, is a separate later change.)</li>
 * <li><b>Picked protein FDR.</b> Each evidence peptide goes to the protein of an
 * accepted genome holding the most evidence peptides (ties by accession); a
 * protein scores its best peptide; target {@code P} competes with
 * {@code rev_P}, built the same way from the decoy evidence peptides.</li>
 * </ol>
 * Peptide uniqueness and neighbourhoods come from digesting the runs'
 * databases, never from the protein column of a PSM table, which names one
 * protein per hit.
 */
public final class GenomeEvidenceTest {

	public enum PeptideClass {
		UNMODIFIED, KNOWN_SHIFT, UNEXPLAINED,
		/**
		 * Closed-stage PSMs of the closed-then-open identification. Evidence; its decoys count in
		 * {@code F} but stay out of the {@code pi_g} window, which keeps the open search's score
		 * scale.
		 */
		CLOSED;

		boolean isEvidence() {
			return this != UNEXPLAINED;
		}
	}

	public static final class Settings {
		public double peptideFdr = 0.01;
		/** Must equal the threshold razor assignment applies, so accepted genomes are the confident ones. */
		public double genomeFdr = DdaMagProteinInferenceService.DEFAULT_GENOME_Q_VALUE_THRESHOLD;
		/** Reported only (count at this q in the run-info); proteins are not filtered on it. */
		public double proteinFdr = 0.01;
		/** The {@code pi_g} window holds the top {@code windowK * F} unique decoy peptides. */
		public int windowK = 30;
		public int maxSubstitutions = 1;
		/** Digest; must match the search that produced the PSMs. */
		public int missedCleavages = 2;
		public int minLength = 7;
		public int maxLength = 50;
		public String decoyPrefix = "rev_";
	}

	/** Streams (accession, sequence) for every protein of the union of the runs' databases, once each. */
	public interface ProteinSource {
		void forEach(ProteinConsumer consumer) throws IOException;
	}

	@FunctionalInterface
	public interface ProteinConsumer {
		void accept(String accession, String sequence);
	}

	/** Spectrum-level hits, aggregated on arrival into peptide units. */
	public static final class Hits {
		private final Map<String, Integer> runIndex = new LinkedHashMap<>();
		private final EnumMap<PeptideClass, HashMap<String, Unit>> targets = new EnumMap<>(PeptideClass.class);
		private final EnumMap<PeptideClass, HashMap<String, Unit>> decoys = new EnumMap<>(PeptideClass.class);
		private final HashMap<String, Unit> windowDecoys = new HashMap<>();
		private final int[][] psmCounts = new int[PeptideClass.values().length][2];

		public Hits() {
			for (PeptideClass c : PeptideClass.values()) {
				targets.put(c, new HashMap<>());
				decoys.put(c, new HashMap<>());
			}
		}

		/**
		 * @param basePeptide unmodified sequence; I is merged into L here
		 * @param score       higher is better
		 */
		public void add(String run, String basePeptide, boolean target, PeptideClass cls, double score) {
			if (run == null || basePeptide == null || basePeptide.isEmpty() || cls == null || !Double.isFinite(score)) {
				return;
			}
			int r = runIndex.computeIfAbsent(run, k -> runIndex.size());
			String key = basePeptide.replace('I', 'L');
			(target ? targets : decoys).get(cls).computeIfAbsent(key, k -> new Unit()).add(r, score);
			if (!target && cls != PeptideClass.CLOSED) {
				windowDecoys.computeIfAbsent(key, k -> new Unit()).add(r, score);
			}
			psmCounts[cls.ordinal()][target ? 0 : 1]++;
		}

		/**
		 * An open-search decoy winner on a spectrum the closed stage identified: in no class, but
		 * part of the open decoys the {@code pi_g} window is drawn from.
		 */
		public void addWindowDecoy(String run, String basePeptide, double score) {
			if (run == null || basePeptide == null || basePeptide.isEmpty() || !Double.isFinite(score)) {
				return;
			}
			int r = runIndex.computeIfAbsent(run, k -> runIndex.size());
			windowDecoys.computeIfAbsent(basePeptide.replace('I', 'L'), k -> new Unit()).add(r, score);
		}

		public Set<String> runs() {
			return Collections.unmodifiableSet(runIndex.keySet());
		}
	}

	private static final class Unit {
		double best = Double.NEGATIVE_INFINITY;
		final BitSet runs = new BitSet();

		void add(int run, double score) {
			runs.set(run);
			if (score > best) {
				best = score;
			}
		}
	}

	/** One target or reversed genome. */
	public static final class GenomeRow {
		public String genome;
		public boolean decoy;
		public int runsSearched;
		/** Evidence peptides unique to the genome before exclusion. */
		public int uniquePeptides;
		public int excludedNearNeighbour;
		public int keptPeptides;
		/** Accepted decoy peptides unique to the reversed counterpart (target rows). */
		public int decoyUniquePeptides;
		public double expectedFalse;
		public double logP;
		public double qValue;
		public boolean accepted;
		/** Position in the test order, from 1; 0 for reversed genomes. */
		public int acceptanceOrder;
		boolean runningAccepted;
	}

	public static final class Calibration {
		public final int[][] psms = new int[PeptideClass.values().length][2];
		public final int[][] acceptedUnits = new int[PeptideClass.values().length][2];
		public final TreeMap<String, Integer> targetUniqueness = new TreeMap<>();
		public final TreeMap<String, Integer> decoyUniqueness = new TreeMap<>();
		public int runs;
		public int candidateGenomes;
		public int decoyEvidencePeptides;
		public int windowK;
		public int windowSize;
		public int windowAvailable;
		public double windowEdgeScore = Double.NaN;
		public double chiSquare = Double.NaN;
		public int chiSquareDf;
		public double chiSquareP = Double.NaN;
		/**
		 * The same test on the open-stage evidence decoys alone; set when a closed class exists.
		 * {@code F} then mixes closed-stage decoys with a window of open-search decoys, so only a
		 * failure of this reading points at {@code pi_g}.
		 */
		public double chiSquareOpen = Double.NaN;
		public int chiSquareOpenDf;
		public double chiSquareOpenP = Double.NaN;
		public int openDecoyEvidencePeptides;
		public int maxSubstitutions;
		public int targetGenomesAccepted;
		public int decoyGenomesAccepted;
		public int runningVersusFinalMismatches;
		public int targetProteinsScored;
		public int decoyProteinsScored;
		public int targetProteinsAtFdr;
		public double genomeFdr;
		public double peptideFdr;
		public double proteinFdr;
	}

	public static final class Result {
		/** Target genomes in test order, then reversed genomes by name. */
		public final List<GenomeRow> rows;
		/** Target genome -> {kept peptides, p-value, q-value}; the layout {@link DdaMagProteinInferenceService#confidentGenomes} reads. */
		public final HashMap<String, double[]> genomeStats;
		/** Target genome -> 1 - q. */
		public final HashMap<String, Double> genomeScoreMap;
		/** Target protein -> {score, score of its reversed counterpart or NaN, picked q-value}. */
		public final HashMap<String, double[]> proteinStats;
		public final Calibration calibration;

		private Result(List<GenomeRow> rows, HashMap<String, double[]> genomeStats,
				HashMap<String, Double> genomeScoreMap, HashMap<String, double[]> proteinStats,
				Calibration calibration) {
			this.rows = rows;
			this.genomeStats = genomeStats;
			this.genomeScoreMap = genomeScoreMap;
			this.proteinStats = proteinStats;
			this.calibration = calibration;
		}
	}

	private static final String NOT_IN_DIGEST = "not_in_digest";
	private static final String NOT_IN_RUN_DB = "not_in_run_database";
	private static final String INCONSISTENT = "unique_to_different_genomes_in_different_runs";
	private static final String SHARED = "shared";
	private static final String UNIQUE = "unique";

	private GenomeEvidenceTest() {
	}

	/**
	 * @param runGenomes genomes in each run's database, keyed by the run names
	 *                   used in {@link Hits#add}
	 * @param proteins   the union of the runs' databases, targets and decoys
	 */
	public static Result run(Hits hits, Map<String, ? extends Set<String>> runGenomes, ProteinSource proteins,
			Settings settings) throws IOException {
		Settings s = settings != null ? settings : new Settings();
		Calibration cal = new Calibration();
		cal.peptideFdr = s.peptideFdr;
		cal.genomeFdr = s.genomeFdr;
		cal.proteinFdr = s.proteinFdr;
		cal.windowK = s.windowK;
		cal.maxSubstitutions = s.maxSubstitutions;
		for (int c = 0; c < hits.psmCounts.length; c++) {
			cal.psms[c] = hits.psmCounts[c].clone();
		}

		List<String> genomes = new ArrayList<>();
		for (Set<String> set : runGenomes.values()) {
			genomes.addAll(set);
		}
		genomes = new ArrayList<>(new java.util.TreeSet<>(genomes));
		HashMap<String, Integer> genomeIndex = new HashMap<>();
		for (String g : genomes) {
			genomeIndex.put(g, genomeIndex.size());
		}
		BitSet[] runDb = new BitSet[hits.runIndex.size()];
		for (Map.Entry<String, Integer> run : hits.runIndex.entrySet()) {
			Set<String> db = runGenomes.get(run.getKey());
			if (db == null) {
				throw new IllegalArgumentException("No database for run " + run.getKey()
						+ "; genome uniqueness is decided against each run's own database");
			}
			BitSet bits = new BitSet();
			for (String g : db) {
				bits.set(genomeIndex.get(g));
			}
			runDb[run.getValue()] = bits;
		}
		cal.runs = runDb.length;
		cal.candidateGenomes = genomes.size();

		// Step 1
		HashMap<String, Evidence> targetEvidence = new HashMap<>();
		HashMap<String, Evidence> decoyEvidence = new HashMap<>();
		java.util.HashSet<String> openDecoyEvidence = new java.util.HashSet<>();
		boolean hasClosedClass = false;
		for (PeptideClass c : PeptideClass.values()) {
			HashMap<String, Evidence> classDecoys = c.isEvidence() ? new HashMap<>() : null;
			int[] accepted = peptideFdr(hits.targets.get(c), hits.decoys.get(c), s.peptideFdr,
					c.isEvidence() ? targetEvidence : null, classDecoys);
			cal.acceptedUnits[c.ordinal()] = accepted;
			if (classDecoys != null) {
				for (Map.Entry<String, Evidence> e : classDecoys.entrySet()) {
					Evidence ev = decoyEvidence.computeIfAbsent(e.getKey(), x -> new Evidence());
					ev.runs.or(e.getValue().runs);
					ev.score = Math.max(ev.score, e.getValue().score);
				}
				if (c != PeptideClass.CLOSED) {
					openDecoyEvidence.addAll(classDecoys.keySet());
				}
			}
			hasClosedClass |= c == PeptideClass.CLOSED
					&& !(hits.targets.get(c).isEmpty() && hits.decoys.get(c).isEmpty());
		}

		// Step 2: exact membership in the digest
		Membership targetMembership = new Membership(targetEvidence.keySet());
		HashMap<String, Unit> windowUnits = hits.windowDecoys;
		java.util.HashSet<String> decoyKeys = new java.util.HashSet<>(decoyEvidence.keySet());
		decoyKeys.addAll(windowUnits.keySet());
		Membership decoyMembership = new Membership(decoyKeys);
		List<String> proteinNames = new ArrayList<>();
		HashMap<String, Integer> proteinIndex = new HashMap<>();
		Digest digest = new Digest(s);
		proteins.forEach((accession, sequence) -> {
			boolean decoy = accession.startsWith(s.decoyPrefix);
			Integer g = genomeIndex.get(genomeOf(decoy ? accession.substring(s.decoyPrefix.length()) : accession));
			if (g == null) {
				return;
			}
			Membership m = decoy ? decoyMembership : targetMembership;
			HashMap<String, Evidence> evidence = decoy ? decoyEvidence : targetEvidence;
			String seq = sequence.replace('I', 'L');
			int[] protein = { -1 };
			digest.forEach(seq, (a, b) -> {
				String key = seq.substring(a, b);
				BitSet bits = m.genomes.get(key);
				if (bits == null) {
					return;
				}
				bits.set(g);
				Evidence ev = evidence.get(key);
				if (ev != null) {
					if (protein[0] < 0) {
						protein[0] = proteinIndex.computeIfAbsent(accession, k -> {
							proteinNames.add(k);
							return proteinNames.size() - 1;
						});
					}
					ev.addProtein(protein[0]);
				}
			});
		});

		HashMap<Integer, List<String>> unique = new HashMap<>();
		HashMap<String, Integer> ownTarget = new HashMap<>();
		for (Map.Entry<String, Evidence> e : targetEvidence.entrySet()) {
			Object status = uniqueGenome(e.getValue().runs, targetMembership.genomes.get(e.getKey()), runDb);
			count(cal.targetUniqueness, status);
			if (status instanceof Integer) {
				unique.computeIfAbsent((Integer) status, k -> new ArrayList<>()).add(e.getKey());
				ownTarget.put(e.getKey(), (Integer) status);
			}
		}
		HashMap<Integer, List<String>> decoyUnique = new HashMap<>();
		HashMap<String, Integer> ownDecoy = new HashMap<>();
		for (Map.Entry<String, Evidence> e : decoyEvidence.entrySet()) {
			Object status = uniqueGenome(e.getValue().runs, decoyMembership.genomes.get(e.getKey()), runDb);
			count(cal.decoyUniqueness, status);
			if (status instanceof Integer) {
				decoyUnique.computeIfAbsent((Integer) status, k -> new ArrayList<>()).add(e.getKey());
				ownDecoy.put(e.getKey(), (Integer) status);
			}
		}
		int F = 0;
		for (List<String> keys : decoyUnique.values()) {
			F += keys.size();
		}
		cal.decoyEvidencePeptides = F;

		List<String> window = new ArrayList<>();
		HashMap<String, Integer> windowGenome = new HashMap<>();
		for (Map.Entry<String, Unit> e : windowUnits.entrySet()) {
			Object status = uniqueGenome(e.getValue().runs, decoyMembership.genomes.get(e.getKey()), runDb);
			if (status instanceof Integer) {
				window.add(e.getKey());
				windowGenome.put(e.getKey(), (Integer) status);
			}
		}
		window.sort(Comparator.comparingDouble((String k) -> -windowUnits.get(k).best)
				.thenComparing(Comparator.naturalOrder()));
		cal.windowAvailable = window.size();
		long windowSize = Math.min((long) window.size(), Math.max(1L, (long) s.windowK * F));
		window = new ArrayList<>(window.subList(0, (int) windowSize));
		cal.windowSize = window.size();
		if (!window.isEmpty()) {
			cal.windowEdgeScore = windowUnits.get(window.get(window.size() - 1)).best;
		}
		int G = genomes.size();
		int[] windowCount = new int[G];
		HashMap<Integer, List<String>> windowKeys = new HashMap<>();
		for (String k : window) {
			int g = windowGenome.get(k);
			windowCount[g]++;
			windowKeys.computeIfAbsent(g, x -> new ArrayList<>()).add(k);
			ownDecoy.putIfAbsent(k, g);
		}
		double lambdaUnit = Math.max(F, 1) / (double) (window.size() + G);
		double[] lambda0 = new double[G];
		for (int g = 0; g < G; g++) {
			lambda0[g] = lambdaUnit * (windowCount[g] + 1);
		}
		double[] chi = chiSquare(decoyUnique, windowCount, window.size(), F);
		cal.chiSquare = chi[0];
		cal.chiSquareDf = (int) chi[1];
		cal.chiSquareP = chi[2];
		if (hasClosedClass) {
			HashMap<Integer, List<String>> openDecoyUnique = new HashMap<>();
			int openF = 0;
			for (Map.Entry<Integer, List<String>> e : decoyUnique.entrySet()) {
				for (String k : e.getValue()) {
					if (openDecoyEvidence.contains(k)) {
						openDecoyUnique.computeIfAbsent(e.getKey(), x -> new ArrayList<>()).add(k);
						openF++;
					}
				}
			}
			cal.openDecoyEvidencePeptides = openF;
			double[] open = chiSquare(openDecoyUnique, windowCount, window.size(), openF);
			cal.chiSquareOpen = open[0];
			cal.chiSquareOpenDf = (int) open[1];
			cal.chiSquareOpenP = open[2];
		}

		// Step 3 neighbourhoods
		List<String> targetQueries = new ArrayList<>(ownTarget.keySet());
		List<String> decoyQueries = new ArrayList<>(ownDecoy.keySet());
		BitSet[] nearTarget = new BitSet[targetQueries.size()];
		BitSet[] nearDecoy = new BitSet[decoyQueries.size()];
		Neighbourhood targetNear = new Neighbourhood(targetQueries, nearTarget, s.maxSubstitutions);
		Neighbourhood decoyNear = new Neighbourhood(decoyQueries, nearDecoy, s.maxSubstitutions);
		proteins.forEach((accession, sequence) -> {
			boolean decoy = accession.startsWith(s.decoyPrefix);
			Integer g = genomeIndex.get(genomeOf(decoy ? accession.substring(s.decoyPrefix.length()) : accession));
			if (g == null) {
				return;
			}
			Neighbourhood n = decoy ? decoyNear : targetNear;
			String seq = sequence.replace('I', 'L');
			digest.forEach(seq, (a, b) -> n.probe(seq, a, b, g));
		});
		HashMap<String, BitSet> nearT = new HashMap<>();
		for (int i = 0; i < targetQueries.size(); i++) {
			BitSet bits = nearTarget[i] != null ? nearTarget[i] : new BitSet();
			bits.clear(ownTarget.get(targetQueries.get(i)));
			nearT.put(targetQueries.get(i), bits);
		}
		HashMap<String, BitSet> nearD = new HashMap<>();
		for (int i = 0; i < decoyQueries.size(); i++) {
			BitSet bits = nearDecoy[i] != null ? nearDecoy[i] : new BitSet();
			bits.clear(ownDecoy.get(decoyQueries.get(i)));
			nearD.put(decoyQueries.get(i), bits);
		}

		// Steps 2-3: sequential test
		double[] logP0 = new double[G];
		Integer[] order = new Integer[G];
		for (int g = 0; g < G; g++) {
			logP0[g] = logPoissonSf(sizeOf(unique.get(g)), lambda0[g]);
			order[g] = g;
		}
		List<String> names = genomes;
		Arrays.sort(order, Comparator.comparingDouble((Integer g) -> logP0[g])
				.thenComparingDouble(g -> -(sizeOf(unique.get(g)) - lambda0[g]) / Math.sqrt(lambda0[g]))
				.thenComparing(names::get));
		BitSet accepted = new BitSet();
		GenomeRow[] targetRows = new GenomeRow[G];
		for (int i = 0; i < G; i++) {
			int g = order[i];
			int kept = surviving(unique.get(g), nearT, accepted);
			double lambda = scaledLambda(lambdaUnit, windowKeys.get(g), nearD, accepted);
			GenomeRow row = new GenomeRow();
			row.genome = names.get(g);
			row.runsSearched = runsSearched(runDb, g);
			row.uniquePeptides = sizeOf(unique.get(g));
			row.keptPeptides = kept;
			row.excludedNearNeighbour = row.uniquePeptides - kept;
			row.decoyUniquePeptides = sizeOf(decoyUnique.get(g));
			row.expectedFalse = lambda;
			row.logP = logPoissonSf(kept, lambda);
			row.acceptanceOrder = i + 1;
			row.runningAccepted = row.logP <= Math.log(s.genomeFdr * (i + 1) / G);
			if (row.runningAccepted) {
				accepted.set(g);
			}
			targetRows[i] = row;
		}
		benjaminiHochberg(targetRows, s.genomeFdr);
		GenomeRow[] decoyRows = new GenomeRow[G];
		for (int g = 0; g < G; g++) {
			BitSet exclude = (BitSet) accepted.clone();
			exclude.clear(g);
			GenomeRow row = new GenomeRow();
			row.genome = s.decoyPrefix + names.get(g);
			row.decoy = true;
			row.runsSearched = runsSearched(runDb, g);
			row.uniquePeptides = sizeOf(decoyUnique.get(g));
			row.keptPeptides = surviving(decoyUnique.get(g), nearD, exclude);
			row.excludedNearNeighbour = row.uniquePeptides - row.keptPeptides;
			row.decoyUniquePeptides = row.uniquePeptides;
			row.expectedFalse = scaledLambda(lambdaUnit, windowKeys.get(g), nearD, exclude);
			row.logP = logPoissonSf(row.keptPeptides, row.expectedFalse);
			decoyRows[g] = row;
		}
		benjaminiHochberg(decoyRows, s.genomeFdr);

		List<GenomeRow> rows = new ArrayList<>(Arrays.asList(targetRows));
		rows.addAll(Arrays.asList(decoyRows));
		HashMap<String, double[]> genomeStats = new HashMap<>();
		HashMap<String, Double> genomeScoreMap = new HashMap<>();
		BitSet finalAccepted = new BitSet();
		for (GenomeRow row : targetRows) {
			genomeStats.put(row.genome, new double[] { row.keptPeptides, Math.exp(row.logP), row.qValue });
			genomeScoreMap.put(row.genome, 1.0 - row.qValue);
			if (row.accepted) {
				cal.targetGenomesAccepted++;
				finalAccepted.set(genomeIndex.get(row.genome));
			}
			if (row.accepted != row.runningAccepted) {
				cal.runningVersusFinalMismatches++;
			}
		}
		for (GenomeRow row : decoyRows) {
			if (row.accepted) {
				cal.decoyGenomesAccepted++;
			}
		}

		// Step 5
		HashMap<String, Double> targetProteinScores = proteinScores(targetEvidence, proteinNames, finalAccepted,
				genomeIndex, s.decoyPrefix);
		HashMap<String, Double> decoyProteinScores = proteinScores(decoyEvidence, proteinNames, finalAccepted,
				genomeIndex, s.decoyPrefix);
		cal.targetProteinsScored = targetProteinScores.size();
		cal.decoyProteinsScored = decoyProteinScores.size();
		HashMap<String, double[]> proteinStats = pickedProteinFdr(targetProteinScores, decoyProteinScores,
				s.decoyPrefix);
		for (double[] st : proteinStats.values()) {
			if (st[2] <= s.proteinFdr) {
				cal.targetProteinsAtFdr++;
			}
		}
		return new Result(rows, genomeStats, genomeScoreMap, proteinStats, cal);
	}

	// ---------------------------------------------------------------- step 1

	private static final class Evidence {
		double score = Double.NEGATIVE_INFINITY;
		final BitSet runs = new BitSet();
		int[] proteins = new int[0];
		int proteinCount;

		void addProtein(int p) {
			for (int i = 0; i < proteinCount; i++) {
				if (proteins[i] == p) {
					return;
				}
			}
			if (proteinCount == proteins.length) {
				proteins = Arrays.copyOf(proteins, Math.max(2, proteinCount * 2));
			}
			proteins[proteinCount++] = p;
		}
	}

	/** @return {accepted targets, accepted decoys} */
	private static int[] peptideFdr(HashMap<String, Unit> targets, HashMap<String, Unit> decoys, double fdr,
			HashMap<String, Evidence> targetEvidence, HashMap<String, Evidence> decoyEvidence) {
		int n = targets.size() + decoys.size();
		String[] keys = new String[n];
		Unit[] units = new Unit[n];
		boolean[] isTarget = new boolean[n];
		int i = 0;
		for (Map.Entry<String, Unit> e : targets.entrySet()) {
			keys[i] = e.getKey();
			units[i] = e.getValue();
			isTarget[i++] = true;
		}
		for (Map.Entry<String, Unit> e : decoys.entrySet()) {
			keys[i] = e.getKey();
			units[i++] = e.getValue();
		}
		Integer[] idx = new Integer[n];
		for (int j = 0; j < n; j++) {
			idx[j] = j;
		}
		Arrays.sort(idx, Comparator.comparingDouble((Integer j) -> -units[j].best)
				.thenComparing(j -> isTarget[j]).thenComparing(j -> keys[j]));
		double[] q = new double[n];
		int t = 0;
		int d = 0;
		for (int j = 0; j < n; j++) {
			if (isTarget[idx[j]]) {
				t++;
			} else {
				d++;
			}
			q[j] = d / (double) Math.max(1, t);
		}
		for (int j = n - 2; j >= 0; j--) {
			q[j] = Math.min(q[j], q[j + 1]);
		}
		int[] accepted = new int[2];
		for (int j = 0; j < n && q[j] <= fdr; j++) {
			int k = idx[j];
			accepted[isTarget[k] ? 0 : 1]++;
			HashMap<String, Evidence> out = isTarget[k] ? targetEvidence : decoyEvidence;
			if (out != null) {
				Evidence ev = out.computeIfAbsent(keys[k], x -> new Evidence());
				ev.runs.or(units[k].runs);
				ev.score = Math.max(ev.score, units[k].best);
			}
		}
		return accepted;
	}

	// ---------------------------------------------------------------- step 2

	private static final class Membership {
		final HashMap<String, BitSet> genomes = new HashMap<>();

		Membership(java.util.Collection<String> keys) {
			for (String k : keys) {
				genomes.put(k, new BitSet());
			}
		}
	}

	/** Genome index, or a status string when the peptide is not unique to one genome. */
	private static Object uniqueGenome(BitSet runs, BitSet genomes, BitSet[] runDb) {
		if (genomes == null || genomes.isEmpty()) {
			return NOT_IN_DIGEST;
		}
		BitSet common = null;
		for (int r = runs.nextSetBit(0); r >= 0; r = runs.nextSetBit(r + 1)) {
			BitSet inRun = (BitSet) genomes.clone();
			inRun.and(runDb[r]);
			if (common == null) {
				common = inRun;
			} else if (!common.equals(inRun)) {
				return INCONSISTENT;
			}
		}
		if (common == null || common.isEmpty()) {
			return NOT_IN_RUN_DB;
		}
		return common.cardinality() == 1 ? (Object) common.nextSetBit(0) : SHARED;
	}

	private static void count(TreeMap<String, Integer> counts, Object status) {
		counts.merge(status instanceof Integer ? UNIQUE : (String) status, 1, Integer::sum);
	}

	/** @return {chi-square, df, p}; NaN when there is too little to test */
	private static double[] chiSquare(HashMap<Integer, List<String>> decoyUnique, int[] windowCount,
			int windowSize, int F) {
		double[] none = { Double.NaN, 0, Double.NaN };
		if (windowSize == 0 || F == 0) {
			return none;
		}
		double observedSum = 0;
		double expectedSum = 0;
		List<double[]> kept = new ArrayList<>();
		for (int g = 0; g < windowCount.length; g++) {
			double e = F * (double) windowCount[g] / windowSize;
			if (e >= 5) {
				double o = sizeOf(decoyUnique.get(g));
				kept.add(new double[] { o, e });
				observedSum += o;
				expectedSum += e;
			}
		}
		if (kept.size() < 2) {
			return none;
		}
		double chi = 0;
		for (double[] oe : kept) {
			double e = oe[1] * observedSum / expectedSum;
			chi += (oe[0] - e) * (oe[0] - e) / e;
		}
		int df = kept.size() - 1;
		return new double[] { chi, df, 1.0 - new ChiSquaredDistribution(df).cumulativeProbability(chi) };
	}

	// ---------------------------------------------------------------- step 3

	/**
	 * Peptides to test, indexed by {@code maxSubstitutions + 1} segments: a
	 * peptide within that many substitutions of a query matches it exactly on at
	 * least one segment (pigeonhole). Candidates are verified by Hamming distance.
	 */
	private static final class Neighbourhood {
		final List<String> queries;
		final BitSet[] near;
		final int d;
		final HashMap<Long, int[]> index = new HashMap<>();
		final BitSet lengths = new BitSet();

		Neighbourhood(List<String> queries, BitSet[] near, int maxSubstitutions) {
			this.queries = queries;
			this.near = near;
			this.d = maxSubstitutions;
			HashMap<Long, List<Integer>> build = new HashMap<>();
			for (int q = 0; q < queries.size(); q++) {
				String p = queries.get(q);
				lengths.set(p.length());
				int[] bounds = bounds(p.length());
				for (int s = 0; s + 1 < bounds.length; s++) {
					build.computeIfAbsent(key(p, 0, p.length(), s, bounds), k -> new ArrayList<>()).add(q);
				}
			}
			for (Map.Entry<Long, List<Integer>> e : build.entrySet()) {
				index.put(e.getKey(), e.getValue().stream().mapToInt(Integer::intValue).toArray());
			}
		}

		int[] bounds(int length) {
			int[] b = new int[d + 2];
			for (int i = 0; i < b.length; i++) {
				b[i] = (int) Math.round(i * length / (double) (d + 1));
			}
			return b;
		}

		long key(String seq, int from, int length, int segment, int[] bounds) {
			long h = 0xcbf29ce484222325L;
			for (int i = from + bounds[segment]; i < from + bounds[segment + 1]; i++) {
				h ^= seq.charAt(i);
				h *= 0x100000001b3L;
			}
			return h * 31 + length * 64L + segment;
		}

		void probe(String seq, int from, int to, int genome) {
			int length = to - from;
			if (!lengths.get(length)) {
				return;
			}
			int[] bounds = bounds(length);
			for (int s = 0; s + 1 < bounds.length; s++) {
				int[] candidates = index.get(key(seq, from, length, s, bounds));
				if (candidates == null) {
					continue;
				}
				for (int q : candidates) {
					String query = queries.get(q);
					if (query.length() == length && withinSubstitutions(seq, from, query, d)) {
						if (near[q] == null) {
							near[q] = new BitSet();
						}
						near[q].set(genome);
					}
				}
			}
		}
	}

	static boolean withinSubstitutions(String seq, int from, String query, int d) {
		int mismatches = 0;
		for (int i = 0; i < query.length(); i++) {
			if (seq.charAt(from + i) != query.charAt(i) && ++mismatches > d) {
				return false;
			}
		}
		return true;
	}

	private static int surviving(List<String> keys, HashMap<String, BitSet> near, BitSet exclude) {
		if (keys == null) {
			return 0;
		}
		int n = 0;
		for (String k : keys) {
			BitSet bits = near.get(k);
			if (bits == null || !bits.intersects(exclude)) {
				n++;
			}
		}
		return n;
	}

	/** {@code F * (surviving window decoys of rev_g + 1) / (window + G)}. */
	private static double scaledLambda(double lambdaUnit, List<String> windowKeys, HashMap<String, BitSet> nearDecoy,
			BitSet exclude) {
		return lambdaUnit * (surviving(windowKeys, nearDecoy, exclude) + 1);
	}

	private static int runsSearched(BitSet[] runDb, int g) {
		int n = 0;
		for (BitSet db : runDb) {
			if (db.get(g)) {
				n++;
			}
		}
		return n;
	}

	/** Benjamini-Hochberg over {@code exp(logP)}, ordered in log space so underflowed p-values keep their order. */
	private static void benjaminiHochberg(GenomeRow[] rows, double fdr) {
		int m = rows.length;
		Integer[] idx = new Integer[m];
		for (int i = 0; i < m; i++) {
			idx[i] = i;
		}
		Arrays.sort(idx, Comparator.comparingDouble(i -> rows[i].logP));
		double running = 1.0;
		for (int j = m - 1; j >= 0; j--) {
			GenomeRow row = rows[idx[j]];
			running = Math.min(running, Math.exp(row.logP + Math.log(m / (double) (j + 1))));
			row.qValue = running;
		}
		for (GenomeRow row : rows) {
			row.accepted = row.qValue <= fdr;
		}
	}

	/** {@code log P(Poisson(lambda) >= n)}, finite where the plain tail underflows. */
	static double logPoissonSf(int n, double lambda) {
		if (n <= 0) {
			return 0.0;
		}
		double p = Gamma.regularizedGammaP(n, lambda);
		if (p > 1e-300) {
			return Math.log(p);
		}
		double total = 1.0;
		double term = 1.0;
		for (int j = 1; term > 1e-17 * total; j++) {
			term *= lambda / (n + j);
			total += term;
		}
		return -lambda + n * Math.log(lambda) - Gamma.logGamma(n + 1.0) + Math.log(total);
	}

	// ---------------------------------------------------------------- step 5

	private static HashMap<String, Double> proteinScores(HashMap<String, Evidence> evidence,
			List<String> proteinNames, BitSet acceptedGenomes, HashMap<String, Integer> genomeIndex,
			String decoyPrefix) {
		HashMap<Integer, Integer> peptidesPerProtein = new HashMap<>();
		List<Evidence> usable = new ArrayList<>();
		for (Evidence ev : evidence.values()) {
			boolean any = false;
			for (int i = 0; i < ev.proteinCount; i++) {
				if (inAccepted(proteinNames.get(ev.proteins[i]), acceptedGenomes, genomeIndex, decoyPrefix)) {
					peptidesPerProtein.merge(ev.proteins[i], 1, Integer::sum);
					any = true;
				}
			}
			if (any) {
				usable.add(ev);
			}
		}
		HashMap<String, Double> scores = new HashMap<>();
		for (Evidence ev : usable) {
			int best = -1;
			int bestCount = -1;
			for (int i = 0; i < ev.proteinCount; i++) {
				int p = ev.proteins[i];
				Integer c = peptidesPerProtein.get(p);
				if (c == null) {
					continue;
				}
				if (c > bestCount || (c == bestCount && proteinNames.get(p).compareTo(proteinNames.get(best)) < 0)) {
					best = p;
					bestCount = c;
				}
			}
			scores.merge(proteinNames.get(best), ev.score, Math::max);
		}
		return scores;
	}

	private static boolean inAccepted(String accession, BitSet acceptedGenomes, HashMap<String, Integer> genomeIndex,
			String decoyPrefix) {
		String target = accession.startsWith(decoyPrefix) ? accession.substring(decoyPrefix.length()) : accession;
		Integer g = genomeIndex.get(genomeOf(target));
		return g != null && acceptedGenomes.get(g);
	}

	/** Target {@code P} against {@code rev_P}; ties go to the decoy. A target that loses gets q = 1. */
	static HashMap<String, double[]> pickedProteinFdr(Map<String, Double> targetScores, Map<String, Double> decoyScores,
			String decoyPrefix) {
		java.util.TreeSet<String> bases = new java.util.TreeSet<>(targetScores.keySet());
		for (String decoy : decoyScores.keySet()) {
			bases.add(decoy.substring(decoyPrefix.length()));
		}
		List<Object[]> winners = new ArrayList<>();
		for (String base : bases) {
			Double t = targetScores.get(base);
			Double d = decoyScores.get(decoyPrefix + base);
			boolean targetWins = t != null && (d == null || t > d);
			winners.add(new Object[] { base, targetWins, targetWins ? t : d });
		}
		winners.sort(Comparator.comparingDouble((Object[] w) -> -(Double) w[2])
				.thenComparing(w -> (Boolean) w[1]).thenComparing(w -> (String) w[0]));
		double[] q = new double[winners.size()];
		int t = 0;
		int d = 0;
		for (int i = 0; i < winners.size(); i++) {
			if ((Boolean) winners.get(i)[1]) {
				t++;
			} else {
				d++;
			}
			q[i] = d / (double) Math.max(1, t);
		}
		for (int i = q.length - 2; i >= 0; i--) {
			q[i] = Math.min(q[i], q[i + 1]);
		}
		HashMap<String, double[]> out = new HashMap<>();
		for (int i = 0; i < winners.size(); i++) {
			String base = (String) winners.get(i)[0];
			Double target = targetScores.get(base);
			if (target == null) {
				continue;
			}
			Double decoy = decoyScores.get(decoyPrefix + base);
			out.put(base, new double[] { target, decoy != null ? decoy : Double.NaN,
					(Boolean) winners.get(i)[1] ? Math.min(1.0, q[i]) : 1.0 });
		}
		return out;
	}

	// ---------------------------------------------------------------- shared

	/** Strict trypsin (after K/R, no proline rule), protein N-terminal Met optionally clipped. */
	private static final class Digest {
		final Settings s;

		Digest(Settings s) {
			this.s = s;
		}

		interface Range {
			void accept(int from, int to);
		}

		void forEach(String seq, Range out) {
			int len = seq.length();
			int[] sites = new int[len + 2];
			int n = 0;
			sites[n++] = 0;
			for (int i = 0; i < len - 1; i++) {
				char c = seq.charAt(i);
				if (c == 'K' || c == 'R') {
					sites[n++] = i + 1;
				}
			}
			sites[n++] = len;
			boolean met = len > 0 && seq.charAt(0) == 'M';
			for (int a = 0; a < n - 1; a++) {
				for (int b = a + 1; b <= Math.min(a + 1 + s.missedCleavages, n - 1); b++) {
					emit(sites[a], sites[b], out);
					if (a == 0 && met) {
						emit(1, sites[b], out);
					}
				}
			}
		}

		private void emit(int from, int to, Range out) {
			int l = to - from;
			if (l >= s.minLength && l <= s.maxLength) {
				out.accept(from, to);
			}
		}
	}

	static String genomeOf(String accession) {
		int idx = accession.indexOf('_');
		return idx > 0 ? accession.substring(0, idx) : accession;
	}

	private static int sizeOf(List<String> list) {
		return list == null ? 0 : list.size();
	}
}
