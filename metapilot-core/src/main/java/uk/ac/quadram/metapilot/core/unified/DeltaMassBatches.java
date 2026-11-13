package uk.ac.quadram.metapilot.core.unified;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.function.Consumer;

import uk.ac.quadram.metapilot.core.opensearch.DeltaMassPeakFitter;
import uk.ac.quadram.metapilot.core.opensearch.ModificationMassDatabase;
import uk.ac.quadram.metapilot.core.opensearch.SubstitutionMassTable;

/**
 * Delta-mass peaks as FDR batches: detected and filtered from mass alone, mapped to known
 * modifications only afterwards.
 * <p>
 * <b>Detection</b> is the existing KDE + Gaussian fit run on the rank-1 open hits of targets and
 * decoys <em>pooled</em>: every hit is passed as a "target", so no target-excess density is
 * formed and no label is read. Isotope families are folded by the fitter (F3).
 * <p>
 * <b>Filters</b>, none of which reads a label, a model score or a database:
 * <ul>
 * <li>F2 - present in at least {@link Config#minFiles} runs;</li>
 * <li>F4 - the zero family is the exact explanation, not a batch;</li>
 * <li>F5 - the delta sits on the peptide: the share of the peak's rank-1 hits whose best placement
 * gains at least {@link Config#minShiftedIons} shifted fragment ions is above the share among hits
 * at off-peak deltas (one-sided binomial, Bonferroni over the peaks tested);</li>
 * <li>F6 - at least {@link Config#minBatchHits} rank-1 hits; a smaller peak joins the nearest
 * retained peak within {@link Config#mergeSigmas} of its sigma, else its hits are unexplained.</li>
 * </ul>
 * <b>Mapping</b> of the retained peaks (bundled modification set, single and pairwise, isotope
 * variants, single amino-acid substitutions) is annotation: it changes no batch.
 */
public final class DeltaMassBatches {

	static final double ISOTOPE = 1.0033548;

	public static final class Config {
		/** F1: rank-1 hit density inside +-3 sigma over the flanks' density, targets and decoys pooled. */
		public double minEnrichment = 3.0;
		/** F1: width of each flank beyond +-3 sigma (Da). */
		public double flankDa = 0.05;
		public int minFiles = 2;
		public int minShiftedIons = 2;
		public double fragmentSupportAlpha = 0.01;
		/** CONGA's minimum group size, 2K with K = 40. */
		public int minBatchHits = 80;
		public double mergeSigmas = 3.0;
		public double mappingToleranceDa = 0.005;
		/** Off-peak deltas used as the fragment-support background must be at least this far from 0 (Da). */
		public double backgroundMinAbsDelta = 1.0;
	}

	/** One family-root peak and what the filters decided. */
	public static final class Candidate {
		public final int id;
		public final double centroid;
		public final double sigma;
		public int rank1Hits;
		public int files;
		public double enrichment = Double.NaN;
		public boolean zero;
		public int supportTested;
		public int supportSuccesses;
		public double supportP = Double.NaN;
		public boolean passFiles;
		public boolean passSupport;
		public boolean passSize;
		/** Retained batch this peak's hits belong to (itself, or the peak it merged into), or -1. */
		public int batch = -1;
		public String verdict = "";
		public String mapping = "";

		Candidate(int id, double centroid, double sigma) {
			this.id = id;
			this.centroid = centroid;
			this.sigma = sigma;
		}
	}

	public final DeltaMassPeakFitter.Result fit;
	public final Map<Integer, Candidate> candidates = new HashMap<>();
	public double backgroundSupportRate = Double.NaN;
	/** Raw rank-1 delta masses, sorted: the F1 density contrast. */
	double[] rank1Deltas = new double[0];

	/** F1: hits per Da inside +-3 sigma over hits per Da in the two flanks. */
	double enrichment(Candidate c, double flankDa) {
		double half = 3 * Math.max(c.sigma, 0.001);
		double inside = count(c.centroid - half, c.centroid + half);
		double flank = count(c.centroid - half - flankDa, c.centroid - half)
				+ count(c.centroid + half, c.centroid + half + flankDa);
		double insideDensity = inside / (2 * half);
		double flankDensity = (flank + 1) / (2 * flankDa);
		return insideDensity / flankDensity;
	}

	private int count(double lo, double hi) {
		int a = java.util.Arrays.binarySearch(rank1Deltas, lo);
		int b = java.util.Arrays.binarySearch(rank1Deltas, hi);
		a = a < 0 ? -a - 1 : a;
		b = b < 0 ? -b - 1 : b;
		return Math.max(0, b - a);
	}

	private DeltaMassBatches(DeltaMassPeakFitter.Result fit) {
		this.fit = fit;
	}

	/**
	 * Detects peaks on the open rows' delta masses.
	 *
	 * @param rank1   rank-1 rows, the only ones used for detection
	 * @param weight  per row, 1 / (rank-1 rows of its spectrum): one vote per spectrum
	 */
	public static DeltaMassBatches detect(double[] delta, boolean[] rank1, double[] weight, String[] files,
			DeltaMassPeakFitter.Config fitterConfig, Consumer<String> log) {
		boolean[] pooled = new boolean[delta.length];
		java.util.Arrays.fill(pooled, true);
		fitterConfig.detectionRows = rank1;
		DeltaMassPeakFitter.Result fit = DeltaMassPeakFitter.fit(delta, rank1, weight, pooled, files, fitterConfig, log);
		DeltaMassBatches out = new DeltaMassBatches(fit);
		Map<Integer, Set<String>> fileSets = new HashMap<>();
		for (DeltaMassPeakFitter.Peak p : fit.peaks) {
			if (p.familyPeakId() == p.id) {
				Candidate c = new Candidate(p.id, p.centroid, Math.max(p.sigma, 1e-4));
				c.zero = p.id == fit.zeroPeakId;
				out.candidates.put(p.id, c);
				fileSets.put(p.id, new HashSet<>());
			}
		}
		for (int i = 0; i < delta.length; i++) {
			if (!rank1[i]) {
				continue;
			}
			Candidate c = out.candidates.get(fit.familyAssignment[i]);
			if (c != null) {
				c.rank1Hits++;
				fileSets.get(c.id).add(files[i]);
			}
		}
		for (Candidate c : out.candidates.values()) {
			c.files = fileSets.get(c.id).size();
		}
		double[] sorted = new double[delta.length];
		int m = 0;
		for (int i = 0; i < delta.length; i++) {
			if (rank1[i] && Double.isFinite(delta[i])) {
				sorted[m++] = delta[i];
			}
		}
		out.rank1Deltas = java.util.Arrays.copyOf(sorted, m);
		java.util.Arrays.sort(out.rank1Deltas);
		log.accept(String.format(Locale.ROOT, "Delta-mass peaks (label-blind): %d family roots detected, zero peak %s",
				out.candidates.size(), fit.zeroPeakId >= 0 ? "found" : "not found"));
		return out;
	}

	/** Family root of open row i, or -1 when it sits on no detected peak. */
	public int familyOf(int row) {
		return fit.familyAssignment[row];
	}

	public boolean isZero(int row) {
		return fit.zeroPeakId >= 0 && fit.familyAssignment[row] == fit.zeroPeakId;
	}

	/**
	 * Applies F2, F5 and F6 and maps the retained peaks.
	 *
	 * @param gained per rank-1 open row with an annotation, shifted ions gained at the best site
	 *               (Integer.MIN_VALUE when not annotated)
	 */
	public void filter(boolean[] rank1, int[] gained, double[] normalizedDelta, Config cfg, Consumer<String> log) {
		// Background: off-peak rank-1 hits far from zero.
		long bgTested = 0;
		long bgSuccess = 0;
		Map<Integer, int[]> support = new HashMap<>();
		for (int i = 0; i < rank1.length; i++) {
			if (!rank1[i] || gained[i] == Integer.MIN_VALUE) {
				continue;
			}
			int fam = fit.familyAssignment[i];
			boolean success = gained[i] >= cfg.minShiftedIons;
			if (fam < 0) {
				if (Math.abs(normalizedDelta[i]) >= cfg.backgroundMinAbsDelta) {
					bgTested++;
					if (success) {
						bgSuccess++;
					}
				}
			} else {
				int[] s = support.computeIfAbsent(fam, k -> new int[2]);
				s[0]++;
				if (success) {
					s[1]++;
				}
			}
		}
		backgroundSupportRate = bgTested == 0 ? 0.0 : bgSuccess / (double) bgTested;
		int tested = 0;
		for (Candidate c : candidates.values()) {
			if (!c.zero) {
				tested++;
			}
		}
		double alpha = cfg.fragmentSupportAlpha / Math.max(1, tested);
		List<Candidate> retained = new ArrayList<>();
		List<Candidate> small = new ArrayList<>();
		for (Candidate c : candidates.values()) {
			c.enrichment = enrichment(c, cfg.flankDa);
			if (c.zero) {
				c.verdict = "zero family: exact explanation";
				continue;
			}
			int[] s = support.getOrDefault(c.id, new int[2]);
			c.supportTested = s[0];
			c.supportSuccesses = s[1];
			c.supportP = binomialUpperTail(s[1], s[0], backgroundSupportRate);
			c.passFiles = c.files >= cfg.minFiles;
			c.passSupport = s[0] > 0 && c.supportP <= alpha;
			c.passSize = c.rank1Hits >= cfg.minBatchHits;
			if (c.enrichment < cfg.minEnrichment) {
				c.verdict = String.format(Locale.ROOT, "F1: density %.2fx the flanks", c.enrichment);
			} else if (!c.passFiles) {
				c.verdict = "F2: seen in " + c.files + " run(s)";
			} else if (!c.passSupport) {
				c.verdict = "F5: no fragment support";
			} else if (c.passSize) {
				c.batch = c.id;
				c.verdict = "retained";
				retained.add(c);
			} else {
				small.add(c);
			}
		}
		for (Candidate c : small) {
			Candidate host = null;
			for (Candidate r : retained) {
				double d = Math.abs(c.centroid - r.centroid) / r.sigma;
				if (d <= cfg.mergeSigmas && (host == null || d < Math.abs(c.centroid - host.centroid) / host.sigma)) {
					host = r;
				}
			}
			if (host != null) {
				c.batch = host.id;
				c.verdict = "F6: merged into peak " + host.id;
			} else {
				c.verdict = "F6: " + c.rank1Hits + " rank-1 hits, below " + cfg.minBatchHits;
			}
		}
		for (Candidate r : retained) {
			r.mapping = map(r.centroid, cfg.mappingToleranceDa);
		}
		log.accept(String.format(Locale.ROOT,
				"Delta-mass batches: %d retained of %d non-zero peaks (fragment-support background %.4f, alpha %.2g); "
						+ "%d small peaks merged",
				retained.size(), tested, backgroundSupportRate, alpha,
				small.stream().filter(c -> c.batch >= 0).count()));
	}

	/** Retained batch of a family root, or -1. */
	public int batchOf(int familyRoot) {
		Candidate c = candidates.get(familyRoot);
		return c == null ? -1 : c.batch;
	}

	public Candidate candidate(int id) {
		return candidates.get(id);
	}

	/** Annotation of a retained peak: known modifications, pairs, isotope variants, substitutions. */
	static String map(double centroid, double tol) {
		List<String> out = new ArrayList<>();
		for (ModificationMassDatabase.Match m : ModificationMassDatabase.matching(centroid, tol)) {
			out.add(m.entry.title);
			if (out.size() >= 3) {
				break;
			}
		}
		for (int k = 1; k <= 2 && out.isEmpty(); k++) {
			for (ModificationMassDatabase.Match m : ModificationMassDatabase.matching(centroid - k * ISOTOPE, tol)) {
				out.add(m.entry.title + " +" + k + " isotope");
				break;
			}
		}
		for (SubstitutionMassTable.Substitution s : SubstitutionMassTable.matching(centroid, tol)) {
			out.add(s.toString());
			if (out.size() >= 5) {
				break;
			}
		}
		if (out.isEmpty()) {
			String pair = pair(centroid, tol);
			if (pair != null) {
				out.add(pair);
			}
		}
		return String.join("; ", out);
	}

	private static String pair(double centroid, double tol) {
		List<ModificationMassDatabase.Entry> all = ModificationMassDatabase.all();
		String best = null;
		double bestErr = tol;
		for (int i = 0; i < all.size(); i++) {
			double rest = centroid - all.get(i).monoMass;
			for (ModificationMassDatabase.Match m : ModificationMassDatabase.matching(rest, tol)) {
				if (Math.abs(m.massErrorDa) <= bestErr) {
					bestErr = Math.abs(m.massErrorDa);
					best = all.get(i).title + " + " + m.entry.title;
				}
				break;
			}
		}
		return best;
	}

	/** P(X >= k) for X ~ Binomial(n, p), in log space. */
	static double binomialUpperTail(int k, int n, double p) {
		if (k <= 0) {
			return 1.0;
		}
		if (p <= 0) {
			return 0.0;
		}
		if (p >= 1) {
			return 1.0;
		}
		double logP = Math.log(p);
		double logQ = Math.log1p(-p);
		double maxTerm = Double.NEGATIVE_INFINITY;
		double[] terms = new double[n - k + 1];
		for (int i = k; i <= n; i++) {
			double t = lnChoose(n, i) + i * logP + (n - i) * logQ;
			terms[i - k] = t;
			maxTerm = Math.max(maxTerm, t);
		}
		double sum = 0;
		for (double t : terms) {
			sum += Math.exp(t - maxTerm);
		}
		return Math.min(1.0, Math.exp(maxTerm + Math.log(sum)));
	}

	private static double lnChoose(int n, int k) {
		return lnGamma(n + 1) - lnGamma(k + 1) - lnGamma(n - k + 1);
	}

	/** Lanczos approximation. */
	private static double lnGamma(double x) {
		double[] c = { 76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155,
				0.1208650973866179e-2, -0.5395239384953e-5 };
		double y = x;
		double tmp = x + 5.5;
		tmp -= (x + 0.5) * Math.log(tmp);
		double ser = 1.000000000190015;
		for (double cc : c) {
			ser += cc / ++y;
		}
		return -tmp + Math.log(2.5066282746310005 * ser / x);
	}

	/** Peak report: every family root, each filter's verdict, the retained batches' mapping. */
	public void writeReport(File out) throws IOException {
		out.getParentFile().mkdirs();
		try (PrintWriter w = new PrintWriter(out, "UTF-8")) {
			w.println("# Delta-mass peaks detected label-blind (targets and decoys pooled); filters F2, F5, F6 read no "
					+ "label; mapping is annotation of retained peaks only.");
			w.println("# fragment-support background (off-peak rank-1 hits gaining shifted ions): "
					+ String.format(Locale.ROOT, "%.5f", backgroundSupportRate));
			w.println("peak_id\tcentroid_da\tsigma_da\trank1_hits\truns\tenrichment\tsupport_tested\tsupport_successes\t"
					+ "support_p\tbatch\tverdict\tmapping");
			List<Candidate> sorted = new ArrayList<>(candidates.values());
			sorted.sort((a, b) -> Double.compare(a.centroid, b.centroid));
			for (Candidate c : sorted) {
				w.println(c.id + "\t" + String.format(Locale.ROOT, "%.5f\t%.5f", c.centroid, c.sigma) + "\t" + c.rank1Hits
						+ "\t" + c.files + "\t"
						+ (Double.isNaN(c.enrichment) ? "" : String.format(Locale.ROOT, "%.2f", c.enrichment)) + "\t"
						+ c.supportTested + "\t" + c.supportSuccesses + "\t"
						+ (Double.isNaN(c.supportP) ? "" : String.format(Locale.ROOT, "%.3g", c.supportP)) + "\t"
						+ c.batch + "\t" + c.verdict + "\t" + c.mapping);
			}
		}
	}
}
