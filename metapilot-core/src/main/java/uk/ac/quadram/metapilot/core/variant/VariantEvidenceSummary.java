package uk.ac.quadram.metapilot.core.variant;

import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

/**
 * {@code variant_pangenome_summary.tsv}: per substitution type, how often the
 * variant peptides occur in UHGP-100 against their matched controls.
 * <p>
 * Percentages are over <em>unique</em> peptides, not candidates: a peptide seen
 * in fifty PSMs is one sequence, found or not, and counting it fifty times
 * would let one abundant peptide decide a group. Controls are grouped with the
 * candidate they were drawn for. {@code enrichment = pct / control_pct}; NA
 * when the control percentage is zero or there are no controls. The Fisher
 * p-value is one-sided (variants found more often than controls) and treats
 * every unique peptide as independent; its q-value is Benjamini-Hochberg over
 * the rows of one {@code level}.
 * <p>
 * Rows, by {@code level}: {@code substitution} per substitution type (for open
 * search also per {@code peak_substitution_shaped}, pooled over peaks);
 * {@code substitution_peak}, open search only, per delta-mass peak;
 * {@code all}, one per (engine, ms2_decidable); and one {@code base} row: the
 * identified base peptides of all candidates, which come from databases built
 * from UHGG genomes and so should be found close to 100% of the time. A low
 * base row means the lookup, not the biology, is off.
 * <p>
 * The open-search MS2 form ranking does not split the rows: its verdicts
 * follow the substitution type more than the biology, so a verdict-filtered
 * subset would mix composition into the enrichment.
 */
final class VariantEvidenceSummary {

	static final String[] COLUMNS = { "engine", "level", "substitution", "ms2_decidable", "peak_id",
			"peak_substitution_shaped", "n_candidates", "n_unique_variant_peptides", "n_unique_in_pangenome_all",
			"pct_in_pangenome_all", "n_unique_in_pangenome_detected", "pct_in_pangenome_detected", "control_n",
			"control_n_in_pangenome_all", "control_pct_in_pangenome_all", "control_n_in_pangenome_detected",
			"control_pct_in_pangenome_detected", "enrichment_all", "enrichment_detected", "fisher_p_all",
			"fisher_p_detected", "fisher_q_all", "fisher_q_detected" };

	static final String LEVEL_SUBSTITUTION = "substitution";
	static final String LEVEL_PEAK = "substitution_peak";
	static final String LEVEL_ALL = "all";
	static final String LEVEL_BASE = "base";

	private static final String NA = "NA";

	static final class Row {
		String engine;
		String level;
		String substitution;
		String ms2Decidable;
		String peakId = NA;
		String peakShaped = NA;
		int nCandidates;
		final Set<String> variants = new HashSet<>();
		final Set<String> controls = new HashSet<>();
		/** {@code [variants found ALL, DETECTED, controls found ALL, DETECTED]}. */
		int[] found;
		double pAll = Double.NaN;
		double pDetected = Double.NaN;
		double qAll = Double.NaN;
		double qDetected = Double.NaN;
		String[] values;
	}

	private VariantEvidenceSummary() {
	}

	static List<Row> summarise(List<VariantCandidate> candidates, List<MatchedControlGenerator.Control> controls,
			PangenomeVariantCheck.Evidence ev, boolean byPeak) {
		Map<String, Row> pooled = new TreeMap<>();
		Map<String, Row> peaks = new TreeMap<>();
		Map<String, Row> totals = new TreeMap<>();
		Row base = new Row();
		Map<VariantCandidate, List<Row>> rowsOf = new java.util.IdentityHashMap<>();
		for (VariantCandidate c : candidates) {
			String engine = PangenomeVariantCheck.engineLabel(c.engine);
			String decidable = String.valueOf(c.ms2Decidable);
			List<Row> rows = new ArrayList<>(3);
			rows.add(pooled(pooled, c, engine, decidable, byPeak));
			if (byPeak && c.peakId >= 0) {
				rows.add(peak(peaks, c, engine, decidable));
			}
			rows.add(total(totals, engine, decidable));
			for (Row r : rows) {
				r.nCandidates++;
				r.variants.add(PeptideAutomaton.normalise(c.variantPeptide));
			}
			rowsOf.put(c, rows);
			base.engine = engine;
			base.nCandidates++;
			base.variants.add(PeptideAutomaton.normalise(c.basePeptide));
		}
		for (MatchedControlGenerator.Control k : controls) {
			for (Row r : rowsOf.getOrDefault(k.parent, List.of())) {
				r.controls.add(PeptideAutomaton.normalise(k.peptide));
			}
		}
		List<Row> out = new ArrayList<>();
		out.addAll(pooled.values());
		out.addAll(peaks.values());
		out.addAll(totals.values());
		if (base.nCandidates > 0) {
			base.level = LEVEL_BASE;
			base.substitution = "base";
			base.ms2Decidable = NA;
			out.add(base);
		}
		for (Row r : out) {
			r.found = found(r, ev);
			if (!LEVEL_BASE.equals(r.level) && !r.controls.isEmpty()) {
				if (ev.reportsAll()) {
					r.pAll = fisherGreater(r.found[0], r.variants.size(), r.found[2], r.controls.size());
				}
				if (ev.reportsDetected()) {
					r.pDetected = fisherGreater(r.found[1], r.variants.size(), r.found[3], r.controls.size());
				}
			}
		}
		for (String level : List.of(LEVEL_SUBSTITUTION, LEVEL_PEAK, LEVEL_ALL)) {
			List<Row> family = new ArrayList<>();
			for (Row r : out) {
				if (level.equals(r.level)) {
					family.add(r);
				}
			}
			double[] qAll = benjaminiHochberg(family.stream().mapToDouble(r -> r.pAll).toArray());
			double[] qDetected = benjaminiHochberg(family.stream().mapToDouble(r -> r.pDetected).toArray());
			for (int i = 0; i < family.size(); i++) {
				family.get(i).qAll = qAll[i];
				family.get(i).qDetected = qDetected[i];
			}
		}
		for (Row r : out) {
			r.values = values(r, ev);
		}
		return out;
	}

	private static Row pooled(Map<String, Row> groups, VariantCandidate c, String engine, String decidable,
			boolean byPeak) {
		String shaped = byPeak ? PangenomeVariantCheck.bool(c.peakSubstitutionShaped) : NA;
		String key = String.join("\u0001", engine, c.substitutionLabel(), decidable, shaped);
		return groups.computeIfAbsent(key, k -> {
			Row r = new Row();
			r.engine = engine;
			r.level = LEVEL_SUBSTITUTION;
			r.substitution = c.substitutionLabel();
			r.ms2Decidable = decidable;
			r.peakId = byPeak ? "ALL" : NA;
			r.peakShaped = shaped;
			return r;
		});
	}

	private static Row peak(Map<String, Row> groups, VariantCandidate c, String engine, String decidable) {
		String peak = String.valueOf(c.peakId);
		String shaped = PangenomeVariantCheck.bool(c.peakSubstitutionShaped);
		String key = String.join("\u0001", engine, c.substitutionLabel(), decidable,
				String.format(Locale.ROOT, "%09d", c.peakId), shaped);
		return groups.computeIfAbsent(key, k -> {
			Row r = new Row();
			r.engine = engine;
			r.level = LEVEL_PEAK;
			r.substitution = c.substitutionLabel();
			r.ms2Decidable = decidable;
			r.peakId = peak;
			r.peakShaped = shaped;
			return r;
		});
	}

	private static Row total(Map<String, Row> totals, String engine, String decidable) {
		return totals.computeIfAbsent(engine + "\u0001" + decidable, k -> {
			Row r = new Row();
			r.engine = engine;
			r.level = LEVEL_ALL;
			r.substitution = "ALL";
			r.ms2Decidable = decidable;
			return r;
		});
	}

	private static String[] values(Row r, PangenomeVariantCheck.Evidence ev) {
		boolean all = ev.reportsAll();
		boolean det = ev.reportsDetected();
		boolean baseRow = LEVEL_BASE.equals(r.level);
		int nv = r.variants.size();
		int nc = r.controls.size();
		boolean controlled = !baseRow && nc > 0;
		return new String[] { r.engine, r.level, r.substitution, r.ms2Decidable, r.peakId, r.peakShaped,
				String.valueOf(r.nCandidates), String.valueOf(nv), all ? String.valueOf(r.found[0]) : NA,
				all ? pct(r.found[0], nv) : NA, det ? String.valueOf(r.found[1]) : NA, det ? pct(r.found[1], nv) : NA,
				baseRow ? NA : String.valueOf(nc), all && controlled ? String.valueOf(r.found[2]) : NA,
				all && controlled ? pct(r.found[2], nc) : NA, det && controlled ? String.valueOf(r.found[3]) : NA,
				det && controlled ? pct(r.found[3], nc) : NA,
				all && controlled ? ratio(r.found[0], nv, r.found[2], nc) : NA,
				det && controlled ? ratio(r.found[1], nv, r.found[3], nc) : NA, p(r.pAll), p(r.pDetected),
				p(r.qAll), p(r.qDetected) };
	}

	private static int[] found(Row r, PangenomeVariantCheck.Evidence ev) {
		int[] v = count(r.variants, ev);
		int[] k = count(r.controls, ev);
		return new int[] { v[0], v[1], k[0], k[1] };
	}

	/** {@code [found ALL, found DETECTED]}. */
	private static int[] count(Set<String> peptides, PangenomeVariantCheck.Evidence ev) {
		int all = 0;
		int detected = 0;
		for (String p : peptides) {
			int q = ev.index(p);
			if (ev.inAll(q)) {
				all++;
			}
			if (Boolean.TRUE.equals(ev.inDetected(q))) {
				detected++;
			}
		}
		return new int[] { all, detected };
	}

	private static String pct(int k, int n) {
		return n == 0 ? NA : String.format(Locale.ROOT, "%.2f", 100.0 * k / n);
	}

	private static String ratio(int k, int n, int controlK, int controlN) {
		if (n == 0 || controlN == 0 || controlK == 0) {
			return NA;
		}
		return String.format(Locale.ROOT, "%.3f", ((double) k / n) / ((double) controlK / controlN));
	}

	private static String p(double p) {
		return Double.isNaN(p) ? NA : String.format(Locale.ROOT, "%.3e", p);
	}

	/**
	 * One-sided Fisher exact test: the probability, with the margins fixed, of
	 * at least {@code k} of the {@code n} variants found when {@code controlK}
	 * of the {@code controlN} controls were. NaN when either group is empty.
	 */
	static double fisherGreater(int k, int n, int controlK, int controlN) {
		if (n == 0 || controlN == 0) {
			return Double.NaN;
		}
		int total = n + controlN;
		int found = k + controlK;
		double[] lf = new double[total + 1];
		for (int i = 2; i <= total; i++) {
			lf[i] = lf[i - 1] + Math.log(i);
		}
		double logDenominator = lf[total] - lf[n] - lf[controlN];
		int from = Math.max(k, n - (total - found));
		int to = Math.min(n, found);
		double p = 0;
		for (int x = from; x <= to; x++) {
			p += Math.exp(lf[found] - lf[x] - lf[found - x] + lf[total - found] - lf[n - x]
					- lf[total - found - n + x] - logDenominator);
		}
		return Math.min(1.0, p);
	}

	/** Benjamini-Hochberg q-values; NaN entries are left out of the family and stay NaN. */
	static double[] benjaminiHochberg(double[] p) {
		double[] q = new double[p.length];
		java.util.Arrays.fill(q, Double.NaN);
		List<Integer> order = new ArrayList<>();
		for (int i = 0; i < p.length; i++) {
			if (!Double.isNaN(p[i])) {
				order.add(i);
			}
		}
		order.sort((a, b) -> Double.compare(p[a], p[b]));
		int m = order.size();
		double running = 1.0;
		for (int rank = m; rank >= 1; rank--) {
			int i = order.get(rank - 1);
			running = Math.min(running, p[i] * m / rank);
			q[i] = running;
		}
		return q;
	}

	static void write(Writer w, List<Row> rows) throws IOException {
		w.write(String.join("\t", COLUMNS));
		w.write('\n');
		for (Row r : rows) {
			w.write(String.join("\t", r.values));
			w.write('\n');
		}
	}
}
