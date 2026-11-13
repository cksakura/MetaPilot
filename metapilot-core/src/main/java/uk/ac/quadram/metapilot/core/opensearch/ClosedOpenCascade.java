package uk.ac.quadram.metapilot.core.opensearch;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.EnumMap;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.function.Consumer;

/**
 * Closed-then-open identification ({@code docs/developer/closed-then-open-identification-plan.md}):
 * the closed stage's validation, the identification classes and the peptide-level filter.
 * <ol>
 * <li><b>Closed stage.</b> The closed search's PSMs, scored by a classifier trained on them
 * alone; per-spectrum target-decoy q, accepted at q &lt;= 1%.</li>
 * <li><b>Open stage.</b> {@link ModificationValidator#validateRemaining} on the open winners of
 * the spectra the closed stage did not accept.</li>
 * <li><b>Peptides.</b> Per class (closed, open unmodified, open accepted shift), base peptide
 * (I = L) at its best score, walked to q &lt;= 1%. A true peptide is often accepted in more than
 * one class and counts once in their union while decoy peptides rarely repeat, so the union can
 * exceed 1%; every class threshold is then lowered by one factor until it does not.</li>
 * </ol>
 * Each stage controls its spectrum-level FDR within its own spectra, so the union of the two is
 * at most 1% (Kertesz-Farkas, Keich and Noble, J Proteome Res 2015).
 */
public final class ClosedOpenCascade {

	/** 13C - 12C. */
	public static final double ISOTOPE_SPACING = 1.0033548;
	/** Number of variable modifications on a closed PSM; an external feature of the closed classifier. */
	public static final String VARIABLE_MODIFICATIONS_FEATURE = "closed_variable_modifications";
	/**
	 * Closed-stage classifier features: search scores, mass error after the isotope offset, charge,
	 * length, cleavage and the variable modifications. No delta-mass or peak feature: a closed PSM
	 * has neither.
	 */
	public static final String[] CLOSED_FEATURES = { "charge", "precursor_mass", "peptide_length", "num_tol_term",
			"num_missed_cleavages", "num_matched_ions", "tot_num_ions", "matched_ion_fraction", "hyperscore",
			"nextscore", "hyperscore_minus_nextscore", "neg_log10_expect", "normalized_mass_diff_ppm",
			"peak_isotope_offset", "log_same_peptide_count", "log_same_protein_count",
			VARIABLE_MODIFICATIONS_FEATURE };
	/** Lowering step of the peptide-union check. */
	static final double UNION_STEP = 0.9;

	public enum IdClass {
		CLOSED("closed"), OPEN_UNMODIFIED("open_unmodified"), OPEN_ACCEPTED_SHIFT("open_accepted_shift");

		public final String label;

		IdClass(String label) {
			this.label = label;
		}
	}

	private ClosedOpenCascade() {
	}

	/** File and scan; the unit both stages share. */
	public static String spectrumKey(OpenSearchPsm psm) {
		return psm.getFileName() + "\u0000" + psm.getScan();
	}

	/** Upper case with I merged into L, as the genome test counts peptides. */
	public static String basePeptide(String peptide) {
		return peptide == null ? "" : peptide.toUpperCase(Locale.ROOT).replace('I', 'L');
	}

	/**
	 * Isotope offset ({@code 0..2}, the closed search's {@code isotope_error}) and the mass error
	 * left after it, in {@code normalizedMassDiff}; the variable-modification count as an external
	 * feature. A closed PSM keeps {@code peakId = -1}.
	 */
	public static void prepareClosed(List<OpenSearchPsm> closed) {
		for (OpenSearchPsm psm : closed) {
			double massDiff = psm.getMassDiff();
			int k = Double.isFinite(massDiff) ? (int) Math.round(massDiff / ISOTOPE_SPACING) : 0;
			k = Math.max(0, Math.min(2, k));
			psm.setPeakIsotopeOffset(k);
			psm.setNormalizedMassDiff(Double.isFinite(massDiff) ? massDiff - k * ISOTOPE_SPACING : Double.NaN);
			Map<String, Double> external = new HashMap<>(psm.getExternalFeatures());
			external.put(VARIABLE_MODIFICATIONS_FEATURE, (double) psm.getClosedModificationCount());
			psm.setExternalFeatures(external);
		}
	}

	/** The closed stage's outcome. */
	public static final class ClosedStage {
		/** One PSM per spectrum, scored. */
		public final List<OpenSearchPsm> winners;
		/** Winners at q &lt;= the threshold, targets and decoys. */
		public final List<OpenSearchPsm> accepted;
		/** Accepted winners by {@link #spectrumKey}. */
		public final Map<String, OpenSearchPsm> acceptedBySpectrum;
		public final int acceptedTargets;
		public final int acceptedDecoys;

		ClosedStage(List<OpenSearchPsm> winners, List<OpenSearchPsm> accepted) {
			this.winners = winners;
			this.accepted = accepted;
			this.acceptedBySpectrum = new HashMap<>();
			int t = 0;
			for (OpenSearchPsm psm : accepted) {
				acceptedBySpectrum.put(spectrumKey(psm), psm);
				if (psm.isTarget()) {
					t++;
				}
			}
			this.acceptedTargets = t;
			this.acceptedDecoys = accepted.size() - t;
		}
	}

	/**
	 * Scores the closed PSMs with a classifier trained on them alone and accepts the spectra at
	 * {@code fdr}. Keeps the best-ranked PSM per spectrum; the search reports one.
	 */
	public static ClosedStage validateClosed(List<OpenSearchPsm> closed, SemiSupervisedPsmClassifier.Config cfg,
			double fdr, Consumer<String> progress) throws Exception {
		return validateClosed(closed, cfg, fdr, progress, false, new String[0]);
	}

	/**
	 * @param rerank        score every closed candidate (all reported ranks, with {@code hit_rank} as a
	 *                      feature) and take each spectrum's best-scoring one as its winner; otherwise
	 *                      only the rank-1 PSM is scored
	 * @param extraFeatures external features added to {@link #CLOSED_FEATURES} (the predicted-spectrum
	 *                      block)
	 */
	public static ClosedStage validateClosed(List<OpenSearchPsm> closed, SemiSupervisedPsmClassifier.Config cfg,
			double fdr, Consumer<String> progress, boolean rerank, String[] extraFeatures) throws Exception {
		Map<String, OpenSearchPsm> bySpectrum = new LinkedHashMap<>();
		for (OpenSearchPsm psm : closed) {
			OpenSearchPsm held = bySpectrum.get(spectrumKey(psm));
			if (held == null || psm.getHitRank() < held.getHitRank()) {
				bySpectrum.put(spectrumKey(psm), psm);
			}
		}
		List<OpenSearchPsm> rankOne = new ArrayList<>(bySpectrum.values());
		List<OpenSearchPsm> pool = rerank ? closed : rankOne;
		prepareClosed(pool);
		// Counts over the rank-1 PSMs, as without re-ranking; a lower-ranked candidate gets its
		// peptide's and protein's rank-1 counts, so the counts do not grow with the report depth.
		OpenSearchFeatures.computeAggregateCounts(rankOne);
		if (rerank) {
			Map<String, Integer> peptideCount = new HashMap<>();
			Map<String, Integer> proteinCount = new HashMap<>();
			for (OpenSearchPsm psm : rankOne) {
				if (psm.getPeptide() != null) {
					peptideCount.put(psm.getPeptide(), psm.getSamePeptideCount());
				}
				if (psm.getProtein() != null) {
					proteinCount.put(psm.getProtein(), psm.getSameProteinCount());
				}
			}
			for (OpenSearchPsm psm : pool) {
				psm.setSamePeptideCount(psm.getPeptide() == null ? 0 : peptideCount.getOrDefault(psm.getPeptide(), 0));
				psm.setSameProteinCount(psm.getProtein() == null ? 0 : proteinCount.getOrDefault(psm.getProtein(), 0));
			}
		}
		List<String> names = new ArrayList<>(Arrays.asList(CLOSED_FEATURES));
		if (rerank) {
			names.add("hit_rank");
		}
		names.addAll(Arrays.asList(extraFeatures));
		if (progress != null) {
			progress.accept(String.format(Locale.ROOT, "Closed stage: scoring %,d candidates of %,d spectra on %d features%s",
					pool.size(), rankOne.size(), names.size(), rerank ? " (all reported ranks, winner after rescoring)" : ""));
		}
		SemiSupervisedPsmClassifier.Config used = cfg;
		if (rerank) {
			// A spectrum's candidates compete on the logit; the clamped logit ties confident ones.
			used = copyOf(cfg);
			used.unclampedLogits = true;
		}
		SemiSupervisedPsmClassifier.Result scored = SemiSupervisedPsmClassifier.classify(pool, used,
				names.toArray(new String[0]), progress);
		for (int i = 0; i < pool.size(); i++) {
			OpenSearchPsm psm = pool.get(i);
			psm.setClassScore(scored.scores[i]);
			psm.setClassLogit(scored.logits[i]);
			psm.setEngineSupport("MSFRAGGER");
		}
		List<OpenSearchPsm> winners = rankOne;
		if (rerank) {
			Map<String, OpenSearchPsm> best = new LinkedHashMap<>();
			for (OpenSearchPsm psm : rankOne) {
				best.put(spectrumKey(psm), psm);
			}
			int overruled = 0;
			for (OpenSearchPsm psm : pool) {
				String key = spectrumKey(psm);
				OpenSearchPsm held = best.get(key);
				double a = psm.getComparableScore();
				double b = held.getComparableScore();
				if (a > b || (a == b && psm.getHitRank() < held.getHitRank())) {
					best.put(key, psm);
				}
			}
			winners = new ArrayList<>(best.values());
			for (OpenSearchPsm psm : winners) {
				if (psm.getHitRank() != 1) {
					overruled++;
				}
			}
			if (progress != null) {
				progress.accept(String.format(Locale.ROOT,
						"Closed stage: the rescored winner is not the closed rank 1 for %,d of %,d spectra (%.1f%%)",
						overruled, winners.size(), 100.0 * overruled / Math.max(1, winners.size())));
			}
		}
		double[] scores = new double[winners.size()];
		boolean[] targets = new boolean[winners.size()];
		for (int i = 0; i < winners.size(); i++) {
			OpenSearchPsm psm = winners.get(i);
			scores[i] = psm.getComparableScore();
			targets[i] = psm.isTarget();
		}
		double[] q = TargetDecoyFdr.qValues(scores, targets);
		List<OpenSearchPsm> accepted = new ArrayList<>();
		for (int i = 0; i < winners.size(); i++) {
			winners.get(i).setQValue(q[i]);
			if (q[i] <= fdr) {
				accepted.add(winners.get(i));
			}
		}
		return new ClosedStage(winners, accepted);
	}

	private static SemiSupervisedPsmClassifier.Config copyOf(SemiSupervisedPsmClassifier.Config cfg) {
		SemiSupervisedPsmClassifier.Config c = new SemiSupervisedPsmClassifier.Config();
		c.folds = cfg.folds;
		c.maxIterations = cfg.maxIterations;
		c.trainQValue = cfg.trainQValue;
		c.positiveFilter = cfg.positiveFilter;
		c.convergenceTolerance = cfg.convergenceTolerance;
		c.minPositives = cfg.minPositives;
		c.stratifyPositiveSeeding = cfg.stratifyPositiveSeeding;
		c.negativeSource = cfg.negativeSource;
		c.maxTrainingInstances = cfg.maxTrainingInstances;
		c.modelFactory = cfg.modelFactory;
		c.unclampedLogits = cfg.unclampedLogits;
		return c;
	}

	/** Peptide-level filter and the union check. */
	public static final class PeptideFilter {
		/** Accepted PSMs whose peptide passed in their class. */
		public final List<OpenSearchPsm> kept = new ArrayList<>();
		/** Per class {target, decoy} peptides at the nominal threshold. */
		public final EnumMap<IdClass, int[]> peptidesBefore = new EnumMap<>(IdClass.class);
		/** Per class {target, decoy} peptides at the final threshold. */
		public final EnumMap<IdClass, int[]> peptidesAfter = new EnumMap<>(IdClass.class);
		/** Per class {target, decoy} spectrum-accepted PSMs in and kept. */
		public final EnumMap<IdClass, int[]> psmsIn = new EnumMap<>(IdClass.class);
		public final EnumMap<IdClass, int[]> psmsKept = new EnumMap<>(IdClass.class);
		public double fdr;
		/** Every class threshold is {@code fdr * factor}. */
		public double factor = 1.0;
		public int unionTargetsBefore;
		public int unionDecoysBefore;
		public int unionTargetsAfter;
		public int unionDecoysAfter;

		public double unionFdrBefore() {
			return unionDecoysBefore / (double) Math.max(1, unionTargetsBefore);
		}

		public double unionFdrAfter() {
			return unionDecoysAfter / (double) Math.max(1, unionTargetsAfter);
		}
	}

	/**
	 * @param pools    per class every scored winner (targets and decoys), the walk's population
	 * @param accepted per class the spectrum-accepted PSMs, a subset of {@code pools}
	 */
	public static PeptideFilter filterPeptides(Map<IdClass, List<OpenSearchPsm>> pools,
			Map<IdClass, List<OpenSearchPsm>> accepted, double fdr) {
		PeptideFilter out = new PeptideFilter();
		out.fdr = fdr;
		EnumMap<IdClass, Map<String, Double>> unitQ = new EnumMap<>(IdClass.class);
		for (IdClass c : IdClass.values()) {
			List<OpenSearchPsm> pool = pools.get(c);
			unitQ.put(c, peptideQValues(pool == null ? new ArrayList<>() : pool));
		}

		double factor = 1.0;
		int[] union = union(unitQ, fdr * factor);
		out.unionTargetsBefore = union[0];
		out.unionDecoysBefore = union[1];
		for (IdClass c : IdClass.values()) {
			out.peptidesBefore.put(c, count(unitQ.get(c), fdr));
		}
		for (int step = 0; union[1] > fdr * Math.max(1, union[0]) && step < 200; step++) {
			factor *= UNION_STEP;
			union = union(unitQ, fdr * factor);
		}
		out.factor = factor;
		out.unionTargetsAfter = union[0];
		out.unionDecoysAfter = union[1];

		double threshold = fdr * factor;
		for (IdClass c : IdClass.values()) {
			out.peptidesAfter.put(c, count(unitQ.get(c), threshold));
			int[] in = new int[2];
			int[] kept = new int[2];
			List<OpenSearchPsm> list = accepted.get(c);
			if (list != null) {
				Map<String, Double> q = unitQ.get(c);
				for (OpenSearchPsm psm : list) {
					int side = psm.isTarget() ? 0 : 1;
					in[side]++;
					Double unit = q.get(unitKey(psm));
					if (unit != null && unit <= threshold) {
						kept[side]++;
						out.kept.add(psm);
					}
				}
			}
			out.psmsIn.put(c, in);
			out.psmsKept.put(c, kept);
		}
		return out;
	}

	private static String unitKey(OpenSearchPsm psm) {
		return (psm.isTarget() ? "T:" : "D:") + basePeptide(psm.getPeptide());
	}

	/** Peptide unit key -> monotonised q, units scored by their best PSM; ties rank decoys first. */
	static Map<String, Double> peptideQValues(List<OpenSearchPsm> pool) {
		Map<String, Double> best = new HashMap<>();
		for (OpenSearchPsm psm : pool) {
			double score = psm.getComparableScore();
			if (psm.getPeptide() == null || !Double.isFinite(score)) {
				continue;
			}
			best.merge(unitKey(psm), score, Math::max);
		}
		List<String> keys = new ArrayList<>(best.keySet());
		keys.sort(Comparator.comparingDouble((String k) -> -best.get(k))
				.thenComparing(k -> k.startsWith("T:")).thenComparing(Comparator.naturalOrder()));
		double[] q = new double[keys.size()];
		int t = 0;
		int d = 0;
		for (int i = 0; i < keys.size(); i++) {
			if (keys.get(i).startsWith("T:")) {
				t++;
			} else {
				d++;
			}
			q[i] = d / (double) Math.max(1, t);
		}
		for (int i = q.length - 2; i >= 0; i--) {
			q[i] = Math.min(q[i], q[i + 1]);
		}
		Map<String, Double> out = new HashMap<>();
		for (int i = 0; i < keys.size(); i++) {
			out.put(keys.get(i), q[i]);
		}
		return out;
	}

	private static int[] count(Map<String, Double> unitQ, double threshold) {
		int[] n = new int[2];
		for (Map.Entry<String, Double> e : unitQ.entrySet()) {
			if (e.getValue() <= threshold) {
				n[e.getKey().startsWith("T:") ? 0 : 1]++;
			}
		}
		return n;
	}

	/** {targets, decoys} among the distinct peptide units accepted in any class. */
	private static int[] union(EnumMap<IdClass, Map<String, Double>> unitQ, double threshold) {
		Set<String> keys = new HashSet<>();
		for (Map<String, Double> q : unitQ.values()) {
			for (Map.Entry<String, Double> e : q.entrySet()) {
				if (e.getValue() <= threshold) {
					keys.add(e.getKey());
				}
			}
		}
		int[] n = new int[2];
		for (String k : keys) {
			n[k.startsWith("T:") ? 0 : 1]++;
		}
		return n;
	}

	/** Lines {@code name<TAB>value} for the run summary. */
	public static List<String> describe(PeptideFilter filter) {
		List<String> lines = new ArrayList<>();
		lines.add(String.format(Locale.ROOT, "peptide_fdr\t%.4f", filter.fdr));
		lines.add(String.format(Locale.ROOT, "peptide_union_targets_before\t%d", filter.unionTargetsBefore));
		lines.add(String.format(Locale.ROOT, "peptide_union_decoys_before\t%d", filter.unionDecoysBefore));
		lines.add(String.format(Locale.ROOT, "peptide_union_fdr_before\t%.5f", filter.unionFdrBefore()));
		lines.add(String.format(Locale.ROOT, "peptide_threshold_factor\t%.5f", filter.factor));
		lines.add(String.format(Locale.ROOT, "peptide_union_targets_after\t%d", filter.unionTargetsAfter));
		lines.add(String.format(Locale.ROOT, "peptide_union_decoys_after\t%d", filter.unionDecoysAfter));
		lines.add(String.format(Locale.ROOT, "peptide_union_fdr_after\t%.5f", filter.unionFdrAfter()));
		for (IdClass c : IdClass.values()) {
			lines.add(c.label + "_peptides_before\t" + Arrays.toString(filter.peptidesBefore.get(c)));
			lines.add(c.label + "_peptides_after\t" + Arrays.toString(filter.peptidesAfter.get(c)));
			lines.add(c.label + "_psms_spectrum_accepted\t" + Arrays.toString(filter.psmsIn.get(c)));
			lines.add(c.label + "_psms_kept\t" + Arrays.toString(filter.psmsKept.get(c)));
		}
		return lines;
	}
}
