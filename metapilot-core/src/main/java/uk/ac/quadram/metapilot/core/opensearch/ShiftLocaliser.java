package uk.ac.quadram.metapilot.core.opensearch;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;

import org.apache.commons.math3.special.Gamma;

import uk.ac.quadram.metapilot.core.variant.MzmlScanReader;
import uk.ac.quadram.metapilot.task.io.qpx.UnimodCatalog;

/**
 * Fragment-supported localisation of an open-search mass shift, and the explanation at the localised
 * site ({@code docs/developer/shift-localisation-explanation-design.md}, sections 0 and 2).
 * <p>
 * Two kinds of hypothesis per PSM: the shift is not on the backbone (H0), or it sits on residue k
 * (Hk), shifting every b ion with number at least k + 1 and every y ion with number at least n - k.
 * Ions are b/y, numbers 1..n-1, charge 1..min(2, z). A hypothesis scores
 * {@code S = -log10 P(X >= matched)}, X binomial over its ions with the mean chance match from the
 * local peak density (±50 Th) and the fragment tolerance; this is length- and density-calibrated, which
 * the search's hyperscore gain is not. An ion whose shifted position the tolerance cannot tell from the
 * unshifted ion or its 13C peaks ({@link #unresolvable}) is evidence for no hypothesis and left out of all.
 * <ul>
 * <li><b>Support</b> - S of the best site's shifted ions minus S of the unshifted ions they replace:
 * the fragment evidence that the shift is on the backbone at all.</li>
 * <li><b>Site probability</b> - P(k) proportional to 10^S(Hk) over every residue.</li>
 * <li><b>Explanation</b> - at the best site only: a substitution of that residue with the right mass
 * ({@link SubstitutionMassTable}), or a UNIMOD entry listing that residue or terminus
 * ({@link UnimodCatalog}); UNIMOD's own amino acid substitution entries are left out, since the
 * substitution table already says so ({@link Result#explanationClass()}).</li>
 * </ul>
 * A PSM is {@link State#LOCALISED_EXPLAINED} only with support, a confident site and an explanation
 * there; {@link State#LOCALISED_UNEXPLAINED} with support and a confident site but nothing known there;
 * {@link State#AMBIGUOUS} otherwise, with nothing claimed beyond its mass. Measured on the
 * pseudo-reversed entrapment runs (design doc, sections 6 and 7): confident sites are right 84-85% of
 * the time against substitution and extension truth on fecal; exported as their own FDR group, the
 * localised-explained PSMs add 1.5% peptides at matched entrapment FDP with every level calibrated.
 * Nothing here reads a target or decoy label.
 */
public final class ShiftLocaliser {

	public enum State {
		/** The fragments place the shift on one residue, where a substitution or UNIMOD entry explains it. */
		LOCALISED_EXPLAINED,
		/** The fragments place the shift, but nothing known explains it there. */
		LOCALISED_UNEXPLAINED,
		/** The fragments carry the shift, but no residue is confident. */
		AMBIGUOUS,
		/** The fragments do not carry the shift: it is not on the backbone. */
		UNSUPPORTED,
		/**
		 * An accepted peak's shift that the fragments do not place, given to a residue or terminus of the
		 * peptide where UNIMOD allows a modification of the peak's mass ({@link ShiftLocaliser#assign}),
		 * as a closed search places a variable modification. The site is the table's, not the spectrum's.
		 */
		ASSIGNED
	}

	/** The kind of explanation a localised shift has at its site. */
	public enum ExplanationClass {
		/** Only UNIMOD modifications fit the residue. */
		MODIFICATION,
		/** Only an amino acid substitution fits the residue. */
		SUBSTITUTION,
		/** Both fit: the fragments are the same, so the spectrum cannot decide. */
		MODIFICATION_OR_SUBSTITUTION
	}

	/** Prefix of a substitution explanation label, {@code subst X->Y}. */
	public static final String SUBSTITUTION_LABEL = "subst ";

	/** Thresholds measured in the design doc, section 6 (stable for support 1-3, probability 0.75-0.9). */
	public static final double DEFAULT_MIN_SUPPORT = 2.0;
	public static final double DEFAULT_MIN_SITE_PROBABILITY = 0.75;
	/** Accepted-shift check: an explained residue (one where a UNIMOD entry explains the shift) within this many log10 units of the best site is consistent. */
	public static final double EXPLAINED_CONSISTENCY_LOG10 = 1.0;
	/** Explanation mass tolerance, as the census and {@link SubstitutionMassTable#matching}. */
	public static final double EXPLANATION_TOLERANCE_DA = 0.005;
	static final double DENSITY_HALF_WINDOW = 50.0;
	static final double PROTON = 1.007276466812;
	static final double H2O = 18.010564684;

	/** One PSM's localisation. {@code site} is 0-based, -1 when not scored. */
	public static final class Result {
		public final double support;
		public final int site;
		public final double siteProbability;
		public final State state;
		/** Candidate explanations at the site ("subst X->Y", "unimod Title"); empty unless explained. */
		public final List<String> explanations;
		/**
		 * The theoretical mass of the first explanation (a substitution's residue difference, or the
		 * UNIMOD entry's mass); NaN unless explained. Exported in place of the observed shift, so every
		 * PSM of one peptide and explanation reports the same mass.
		 */
		public final double explainedMass;

		Result(double support, int site, double siteProbability, State state, List<String> explanations) {
			this(support, site, siteProbability, state, explanations, Double.NaN);
		}

		public Result(double support, int site, double siteProbability, State state, List<String> explanations,
				double explainedMass) {
			this.support = support;
			this.site = site;
			this.siteProbability = siteProbability;
			this.state = state;
			this.explanations = explanations;
			this.explainedMass = explainedMass;
		}

		/** Whether an explanation and its mass are claimed: localised and explained, or assigned from UNIMOD. */
		public boolean isExplained() {
			return state == State.LOCALISED_EXPLAINED || state == State.ASSIGNED;
		}

		/**
		 * What kind of explanation the site has; null unless explained. A substitution and a
		 * modification on one residue move the same fragments by the same mass, so with both the
		 * spectrum cannot say which it is.
		 */
		public ExplanationClass explanationClass() {
			if (!isExplained()) {
				return null;
			}
			boolean substitution = false;
			boolean modification = false;
			for (String label : explanations) {
				if (label.startsWith(SUBSTITUTION_LABEL)) {
					substitution = true;
				} else {
					modification = true;
				}
			}
			return substitution && modification ? ExplanationClass.MODIFICATION_OR_SUBSTITUTION
					: substitution ? ExplanationClass.SUBSTITUTION : ExplanationClass.MODIFICATION;
		}

		/** The residue the first substitution explanation puts at the site; 0 without one. */
		public char substitutedResidue() {
			for (String label : explanations) {
				if (label.startsWith(SUBSTITUTION_LABEL) && label.length() == SUBSTITUTION_LABEL.length() + 4) {
					return label.charAt(SUBSTITUTION_LABEL.length() + 3);
				}
			}
			return 0;
		}

		/** {@link #explanationClass()} as a table cell: lower case, empty unless explained. */
		public String explanationClassText() {
			ExplanationClass c = explanationClass();
			return c == null ? "" : c.name().toLowerCase(Locale.ROOT);
		}

		static final Result NOT_SCORED = new Result(Double.NaN, -1, Double.NaN, State.AMBIGUOUS, List.of());

		public String explanationText() {
			return String.join("|", explanations);
		}
	}

	private final double tolerancePpm;
	private final double minSupport;
	private final double minSiteProbability;
	private final UnimodCatalog unimod;

	public ShiftLocaliser(double tolerancePpm) {
		this(tolerancePpm, DEFAULT_MIN_SUPPORT, DEFAULT_MIN_SITE_PROBABILITY, UnimodCatalog.get());
	}

	ShiftLocaliser(double tolerancePpm, double minSupport, double minSiteProbability, UnimodCatalog unimod) {
		this.tolerancePpm = tolerancePpm;
		this.minSupport = minSupport;
		this.minSiteProbability = minSiteProbability;
		this.unimod = unimod;
	}

	/** Site scores of one PSM: S(Hk) per residue, the best site, and the fragment support there. */
	private static final class Fit {
		double[] siteScores;
		int best;
		double support;

		/** P(k) at a residue, from the site scores. */
		double probability(int k) {
			double norm = 0;
			for (double sk : siteScores) {
				norm += Math.pow(10, Math.max(-30, sk - siteScores[best]));
			}
			return Math.pow(10, Math.max(-30, siteScores[k] - siteScores[best])) / norm;
		}
	}

	/** Null when the PSM cannot be scored (no spectrum, an unknown residue, a non-finite shift). */
	private Fit fit(String peptide, int charge, double delta, MzmlScanReader.Peaks peaks) {
		if (peptide == null || peaks == null || peaks.mz.length == 0 || !Double.isFinite(delta)) {
			return null;
		}
		int n = peptide.length();
		if (n < 2) {
			return null;
		}
		if (n < 2 || peaks == null || peaks.mz.length == 0 || !Double.isFinite(delta)) {
			return null;
		}
		double[] prefix = new double[n + 1];
		for (int i = 0; i < n; i++) {
			double m = SubstitutionMassTable.residueMass(peptide.charAt(i));
			if (m <= 0) {
				return null;
			}
			prefix[i + 1] = prefix[i] + m;
		}
		int maxCharge = Math.max(1, Math.min(2, charge));
		int all = 2 * (n - 1) * maxCharge;
		boolean[] isB = new boolean[all];
		int[] number = new int[all];
		boolean[] matchedPlain = new boolean[all];
		boolean[] matchedShifted = new boolean[all];
		double[] chancePlain = new double[all];
		double[] chanceShifted = new double[all];
		int t = 0;
		for (int j = 1; j <= n - 1; j++) {
			double b = prefix[j];
			double y = prefix[n] - prefix[n - j] + H2O;
			for (int c = 1; c <= maxCharge; c++) {
				for (int s = 0; s < 2; s++) {
					double neutral = s == 0 ? b : y;
					double mz = (neutral + c * PROTON) / c;
					double mzShifted = (neutral + delta + c * PROTON) / c;
					if (unresolvable(neutral, delta, c, mzShifted)) {
						// The fragments cannot tell this ion shifted from unshifted: its shifted position is
						// within tolerance of the unshifted ion itself or of one of its 13C peaks, and peak
						// lists often hold a large fragment's 13C peak without its monoisotopic one. It is
						// evidence for no hypothesis, so it is left out of all of them.
						continue;
					}
					isB[t] = s == 0;
					number[t] = j;
					matchedPlain[t] = matched(peaks, mz);
					chancePlain[t] = chance(peaks, mz);
					if (mzShifted > 0) {
						matchedShifted[t] = matched(peaks, mzShifted);
						chanceShifted[t] = chance(peaks, mzShifted);
					}
					t++;
				}
			}
		}
		int ions = t;
		if (ions == 0) {
			return null;
		}
		double[] sk = new double[n];
		int best = 0;
		for (int k = 0; k < n; k++) {
			int hits = 0;
			double p = 0;
			for (int i = 0; i < ions; i++) {
				boolean shifted = carriesSite(isB[i], number[i], k, n);
				hits += (shifted ? matchedShifted[i] : matchedPlain[i]) ? 1 : 0;
				p += shifted ? chanceShifted[i] : chancePlain[i];
			}
			sk[k] = significance(hits, ions, p / ions);
			if (sk[k] > sk[best]) {
				best = k;
			}
		}
		int m = 0;
		int hitsShifted = 0;
		int hitsPlain = 0;
		double pShifted = 0;
		double pPlain = 0;
		for (int i = 0; i < ions; i++) {
			if (carriesSite(isB[i], number[i], best, n)) {
				m++;
				hitsShifted += matchedShifted[i] ? 1 : 0;
				hitsPlain += matchedPlain[i] ? 1 : 0;
				pShifted += chanceShifted[i];
				pPlain += chancePlain[i];
			}
		}
		double support = m == 0 ? 0.0
				: significance(hitsShifted, m, pShifted / m) - significance(hitsPlain, m, pPlain / m);
		Fit f = new Fit();
		f.siteScores = sk;
		f.best = best;
		f.support = support;
		return f;
	}

	/** Isotope peaks of an unshifted ion checked against its shifted position: +0..+3 x 13C. */
	static final int UNRESOLVABLE_ISOTOPES = 3;
	static final double C13 = 1.0033548;

	/**
	 * Whether an ion's shifted position ({@code mzShifted}, the ion with {@code delta} added) lies within
	 * the fragment tolerance of the unshifted ion or of one of its 13C peaks (+1, +2, +3 x 1.00335 Da),
	 * so that no spectrum at this tolerance can tell the two apart. Depends on the masses and the
	 * tolerance only, never on the spectrum. At 20 ppm, deamidation (+0.98402, 19.3 mDa below the +1
	 * isotope step) is unresolvable for singly charged ions above m/z 967 and doubly charged ones above
	 * m/z 483. On 12mix, accepted unmodified PSMs scored at +0.984 were supported in 8.6% of cases
	 * without this rule (216 of 244 decisive ions sat on a 13C peak whose monoisotopic peak was absent)
	 * and 2.9% with it; closed 0.984 PSMs went from 81% to 42%.
	 */
	boolean unresolvable(double neutral, double delta, int charge, double mzShifted) {
		double tol = mzShifted * tolerancePpm * 1e-6;
		for (int k = 0; k <= UNRESOLVABLE_ISOTOPES; k++) {
			double isotope = (neutral + k * C13 + charge * PROTON) / charge;
			if (Math.abs(mzShifted - isotope) <= tol) {
				return true;
			}
		}
		return false;
	}

	/** Localises shift {@code delta} on {@code peptide} (charge {@code charge}) in {@code peaks}. */
	public Result localise(String peptide, int charge, double delta, MzmlScanReader.Peaks peaks) {
		Fit f = fit(peptide, charge, delta, peaks);
		if (f == null) {
			return Result.NOT_SCORED;
		}
		double siteProbability = f.probability(f.best);
		if (f.support < minSupport) {
			return new Result(f.support, f.best, siteProbability, State.UNSUPPORTED, List.of());
		}
		if (siteProbability < minSiteProbability) {
			return new Result(f.support, f.best, siteProbability, State.AMBIGUOUS, List.of());
		}
		List<String> explanations = explanationsAt(peptide, f.best, delta);
		if (explanations.isEmpty()) {
			return new Result(f.support, f.best, siteProbability, State.LOCALISED_UNEXPLAINED, explanations);
		}
		return new Result(f.support, f.best, siteProbability, State.LOCALISED_EXPLAINED, explanations,
				massOf(explanations.get(0), delta));
	}

	/**
	 * Checks a shift on an accepted peak (design doc, section 8, where it is called the allow-list
	 * check): {@link State#UNSUPPORTED} without fragment support; {@link State#LOCALISED_EXPLAINED}
	 * at the best explained residue (one where a UNIMOD entry explains the shift) when the fragments
	 * do not contradict it (its score within {@link #EXPLAINED_CONSISTENCY_LOG10} of the best site);
	 * {@link State#LOCALISED_UNEXPLAINED} at the best site when they contradict every explained residue.
	 * Sites the fragments cannot tell apart (residues 1 and 2 differ only by b1) share probability, so
	 * the test is consistency, not a site-probability threshold.
	 */
	public Result checkAcceptedShift(String peptide, int charge, double delta, MzmlScanReader.Peaks peaks) {
		Fit f = fit(peptide, charge, delta, peaks);
		if (f == null) {
			return Result.NOT_SCORED;
		}
		if (f.support < minSupport) {
			return new Result(f.support, f.best, f.probability(f.best), State.UNSUPPORTED, List.of());
		}
		int bestExplained = -1;
		List<String> explanations = List.of();
		for (int k = 0; k < peptide.length(); k++) {
			List<String> e = explanationsAt(peptide, k, delta);
			if (!e.isEmpty() && (bestExplained < 0 || f.siteScores[k] > f.siteScores[bestExplained])) {
				bestExplained = k;
				explanations = e;
			}
		}
		if (bestExplained < 0 || f.siteScores[bestExplained] < f.siteScores[f.best] - EXPLAINED_CONSISTENCY_LOG10) {
			return new Result(f.support, f.best, f.probability(f.best), State.LOCALISED_UNEXPLAINED, List.of());
		}
		return new Result(f.support, bestExplained, f.probability(bestExplained), State.LOCALISED_EXPLAINED,
				explanations, massOf(explanations.get(0), delta));
	}

	/**
	 * The modification an accepted peak's mass is on this peptide by UNIMOD alone, for a PSM whose
	 * fragments do not place the shift: the first residue or terminus where an entry of that mass is
	 * listed as a common specificity, else the first where one is listed at all (rarely used
	 * specificities count; isotopic labels and amino acid substitutions do not), as
	 * {@link State#ASSIGNED}. Common first, because UNIMOD lists oxidation as rare on most residues and
	 * the first residue of a peptide would otherwise take it from the methionine. Null when the peptide
	 * has no such site: the mass is then not accounted for and the PSM is no identification.
	 * <p>
	 * Why: such a PSM was exported as the unmodified peptide, so FlashLFQ looked for a mass the
	 * fragmented precursor did not have. On the round-5 entrapment runs (2026-11-10) 425 of 6,219
	 * (faecal) and 8,628 of 51,251 (12mix) such PSMs had a UNIMOD site, with 0 and 8 entrapment PSMs;
	 * the rest were almost all in the +1.95-1.98, +4.0-4.1 and +4.9-5.0 Da bands, which no UNIMOD
	 * entry has, and held 0.4% and 0.7% of the peptides.
	 *
	 * @param mass    the peak's mass (its centroid; a single PSM's own delta carries its mass error)
	 * @param support the PSM's fragment support, kept for the record
	 */
	public Result assign(String peptide, double mass, double support) {
		if (peptide == null || !Double.isFinite(mass)) {
			return null;
		}
		int any = -1;
		int common = -1;
		for (int k = 0; k < peptide.length() && common < 0; k++) {
			for (UnimodCatalog.Site site : sitesAt(peptide, k)) {
				for (UnimodCatalog.Entry e : unimod.allByMass(mass, site, EXPLANATION_TOLERANCE_DA)) {
					if (e.isAminoAcidSubstitutionOnly()) {
						continue;
					}
					any = any < 0 ? k : any;
					common = common < 0 && e.listsCommonly(site) ? k : common;
				}
			}
		}
		int at = common >= 0 ? common : any;
		if (at < 0) {
			return null;
		}
		List<String> modifications = new ArrayList<>();
		for (String label : explanationsAt(peptide, at, mass)) {
			if (!label.startsWith(SUBSTITUTION_LABEL)) {
				modifications.add(label);
			}
		}
		return new Result(support, at, Double.NaN, State.ASSIGNED, modifications, massOf(modifications.get(0), mass));
	}

	/** The UNIMOD sites residue k of a peptide can be: the residue, and a terminus when it is one. */
	private static List<UnimodCatalog.Site> sitesAt(String peptide, int k) {
		char residue = Character.toUpperCase(peptide.charAt(k));
		List<UnimodCatalog.Site> sites = new ArrayList<>();
		sites.add(UnimodCatalog.Site.residue(residue));
		if (k == 0) {
			sites.add(UnimodCatalog.Site.nTermResidue(residue));
			sites.add(UnimodCatalog.Site.nTerm());
		}
		if (k == peptide.length() - 1) {
			sites.add(UnimodCatalog.Site.cTermResidue(residue));
			sites.add(UnimodCatalog.Site.cTerm());
		}
		return sites;
	}

	/**
	 * Fragment support of shift {@code delta} at its best site (as {@link Result#support}); NaN when the
	 * PSM cannot be scored. The per-peak fragment test ({@link PeakFragmentTest}) scores unmodified PSMs
	 * with this at a peak's mass to learn how often support reaches {@link #minSupport()} by chance.
	 */
	public double supportAt(String peptide, int charge, double delta, MzmlScanReader.Peaks peaks) {
		Fit f = fit(peptide, charge, delta, peaks);
		return f == null ? Double.NaN : f.support;
	}

	/** The support a PSM needs before its shift counts as carried by the fragments. */
	public double minSupport() {
		return minSupport;
	}

	/**
	 * {@code r} as {@link State#UNSUPPORTED}: the peak's shift is not on the fragments
	 * ({@link PeakFragmentTest}), so nothing is localised or explained, whatever this one PSM showed.
	 */
	static Result asUnsupported(Result r) {
		if (r == null || r.state == State.UNSUPPORTED) {
			return r;
		}
		return new Result(r.support, r.site, r.siteProbability, State.UNSUPPORTED, List.of());
	}

	/** Theoretical mass of an explanation label; {@code delta} when the label cannot be resolved. */
	double massOf(String label, double delta) {
		if (label.startsWith("subst ") && label.length() == "subst ".length() + 4) {
			char from = label.charAt("subst ".length());
			char to = label.charAt("subst ".length() + 3);
			return SubstitutionMassTable.residueMass(to) - SubstitutionMassTable.residueMass(from);
		}
		if (label.startsWith("unimod ")) {
			UnimodCatalog.Entry e = unimod.byTitle(label.substring("unimod ".length()));
			if (e != null) {
				return e.monoMass;
			}
		}
		return delta;
	}

	/** Whether ion b_number (or y_number) contains residue k of an n-residue peptide. */
	static boolean carriesSite(boolean bIon, int number, int k, int n) {
		return bIon ? number >= k + 1 : number >= n - k;
	}

	/** Substitutions of residue k with mass {@code delta}, then UNIMOD entries listing that residue or terminus. */
	List<String> explanationsAt(String peptide, int k, double delta) {
		char residue = Character.toUpperCase(peptide.charAt(k));
		int n = peptide.length();
		List<String> out = new ArrayList<>();
		for (SubstitutionMassTable.Substitution s : SubstitutionMassTable.matching(delta, EXPLANATION_TOLERANCE_DA)) {
			if (s.from == residue || s.from == 'L' && residue == 'I' || s.from == 'I' && residue == 'L') {
				String label = "subst " + residue + "->" + s.to;
				if (!out.contains(label)) {
					out.add(label);
				}
			}
		}
		for (UnimodCatalog.Site site : sitesAt(peptide, k)) {
			for (UnimodCatalog.Entry e : unimod.allByMass(delta, site, EXPLANATION_TOLERANCE_DA)) {
				if (e.isAminoAcidSubstitutionOnly()) {
					// A sequence difference under a UNIMOD name: the substitution table above is its
					// explanation, and listing it again would make every substitution look like a tie.
					continue;
				}
				String label = "unimod " + e.title;
				if (!out.contains(label)) {
					out.add(label);
				}
			}
		}
		return out;
	}

	boolean matched(MzmlScanReader.Peaks peaks, double mz) {
		double tol = mz * tolerancePpm * 1e-6;
		return peaks.maxIntensity(mz - tol, mz + tol) > 0f;
	}

	double chance(MzmlScanReader.Peaks peaks, double mz) {
		int lo = lowerBound(peaks.mz, mz - DENSITY_HALF_WINDOW);
		int hi = lowerBound(peaks.mz, mz + DENSITY_HALF_WINDOW);
		return Math.min(1.0, (hi - lo) / (2 * DENSITY_HALF_WINDOW) * 2 * mz * tolerancePpm * 1e-6);
	}

	private static int lowerBound(double[] a, double v) {
		int i = Arrays.binarySearch(a, v);
		return i < 0 ? -i - 1 : i;
	}

	/** -log10 P(X >= k), X ~ Binomial(n, p). */
	static double significance(int k, int n, double p) {
		if (n == 0 || k == 0 || p >= 1) {
			return 0.0;
		}
		if (p <= 0) {
			return 30.0;
		}
		double lp = Math.log(p);
		double lq = Math.log1p(-p);
		double[] terms = new double[n - k + 1];
		double max = Double.NEGATIVE_INFINITY;
		for (int x = k; x <= n; x++) {
			terms[x - k] = Gamma.logGamma(n + 1.0) - Gamma.logGamma(x + 1.0) - Gamma.logGamma(n - x + 1.0) + x * lp
					+ (n - x) * lq;
			max = Math.max(max, terms[x - k]);
		}
		double sum = 0;
		for (double v : terms) {
			sum += Math.exp(v - max);
		}
		return -(max + Math.log(sum)) / Math.log(10);
	}

	@Override
	public String toString() {
		return String.format(Locale.ROOT, "fragment support >= %.1f, site probability >= %.2f, %.0f ppm", minSupport,
				minSiteProbability, tolerancePpm);
	}
}
