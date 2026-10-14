package uk.ac.quadram.metapilot.core.variant;

import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Locale;
import java.util.Set;

import uk.ac.quadram.metapilot.core.opensearch.SubstitutionMassTable;

/**
 * Ranks the forms that could explain one open-search PSM, within its own
 * spectrum (Kai's "proteoform ranking", agreed 2026-09-14).
 * <p>
 * The forms for a candidate X->Y at site k with delta d on assigned peptide P:
 * <ul>
 * <li><b>substitution</b> - P with Y at k, no delta: this candidate.</li>
 * <li><b>other substitution</b> - every other single substitution at any site
 * whose mass matches d ({@link SubstitutionMassTable#matching}).</li>
 * <li><b>modification</b> - P with d placed on residue j, for every j. A delta
 * on a residue no catalogued modification targets is still a valid open-search
 * explanation.</li>
 * </ul>
 * A modification at k and the substitution at k give identical b/y m/z; only
 * fragment intensities separate them. Each form is therefore scored with
 * predicted intensities: a substituted peptide with its own prediction, a
 * modification form with the assigned backbone's prediction and d added to
 * the fragments that contain its site.
 * <p>
 * The verdict compares only forms that change site k. Forms at other sites
 * differ in fragment m/z as well, and the best of about one placement per
 * residue would beat the single substitution form by selection alone; they
 * are reported separately as the "elsewhere" score.
 * <p>
 * The score is the unweighted spectral-entropy similarity (Li et al., 2021;
 * MSBooster's {@code unweighted_spectral_entropy}) between the predicted
 * intensities and the most intense observed peak within tolerance of each
 * predicted fragment. It is compared only between forms of the same spectrum,
 * never against an absolute threshold: the baseline similarity moves with
 * length and charge.
 */
final class VariantFormRanker {

	static final double PROTON = 1.007276466812;
	static final double H2O = 18.010564684;
	private static final double LN4 = Math.log(4.0);

	enum Kind {
		SUBSTITUTION, OTHER_SUBSTITUTION, MODIFICATION
	}

	/** One explanation of the spectrum. */
	static final class Form {
		final Kind kind;
		/** Sequence the fragments are predicted for. */
		final String backbone;
		/** 0-based residue this form changes, substituted or modified. */
		final int site;
		/** 0-based residue carrying {@link #shift}; -1 for none. */
		final int shiftSite;
		final double shift;
		final String label;

		Form(Kind kind, String backbone, int site, int shiftSite, double shift, String label) {
			this.kind = kind;
			this.backbone = backbone;
			this.site = site;
			this.shiftSite = shiftSite;
			this.shift = shift;
			this.label = label;
		}
	}

	interface Predictions {
		/** Null when the backbone was not predicted. */
		DiannFragmentPredictor.Prediction get(String backbone, int charge);
	}

	private VariantFormRanker() {
	}

	/**
	 * Every form for this candidate; the candidate's own substitution first.
	 * Empty when the peptide holds a residue without a mass.
	 */
	static List<Form> forms(VariantCandidate c, double substitutionToleranceDa) {
		List<Form> out = new ArrayList<>();
		String base = c.basePeptide;
		if (!standard(base) || !standard(c.variantPeptide)) {
			return out;
		}
		out.add(new Form(Kind.SUBSTITUTION, c.variantPeptide, c.position - 1, -1, 0.0,
				"substitution " + c.fromAa + c.position + c.toAa));
		Set<String> seen = new LinkedHashSet<>();
		seen.add(PeptideAutomaton.normalise(c.variantPeptide));
		List<SubstitutionMassTable.Substitution> matching = SubstitutionMassTable.matching(c.deltaMass,
				substitutionToleranceDa);
		for (int j = 0; j < base.length(); j++) {
			char residue = base.charAt(j);
			for (SubstitutionMassTable.Substitution s : matching) {
				if (s.from != residue) {
					continue;
				}
				StringBuilder v = new StringBuilder(base);
				v.setCharAt(j, s.to);
				String variant = v.toString();
				if (seen.add(PeptideAutomaton.normalise(variant))) {
					out.add(new Form(Kind.OTHER_SUBSTITUTION, variant, j, -1, 0.0,
							"substitution " + s.from + (j + 1) + s.to));
				}
			}
		}
		for (int j = 0; j < base.length(); j++) {
			out.add(new Form(Kind.MODIFICATION, base, j, j, c.deltaMass,
					String.format(Locale.ROOT, "modification %+.4f@%c%d", c.deltaMass, base.charAt(j), j + 1)));
		}
		return out;
	}

	private static boolean standard(String peptide) {
		if (peptide == null || peptide.isEmpty()) {
			return false;
		}
		for (int i = 0; i < peptide.length(); i++) {
			if ("ACDEFGHIKLMNPQRSTVWY".indexOf(peptide.charAt(i)) < 0) {
				return false;
			}
		}
		return true;
	}

	/** Fragment m/z, with {@code shift} on the fragments that contain {@code shiftSite}. */
	static double fragmentMz(String backbone, int ionNumber, boolean yIon, int charge, int shiftSite, double shift) {
		int n = backbone.length();
		int from = yIon ? n - ionNumber : 0;
		int to = yIon ? n : ionNumber;
		double mass = yIon ? H2O : 0.0;
		for (int i = from; i < to; i++) {
			mass += SubstitutionMassTable.residueMass(backbone.charAt(i));
		}
		if (shiftSite >= from && shiftSite < to) {
			mass += shift;
		}
		return (mass + charge * PROTON) / charge;
	}

	/** Similarity of a form to the spectrum; NaN when there is nothing to compare. */
	static double score(Form f, DiannFragmentPredictor.Prediction p, MzmlScanReader.Peaks spectrum,
			int precursorCharge, double tolerancePpm) {
		if (p == null || spectrum == null) {
			return Double.NaN;
		}
		int n = f.backbone.length();
		double[] predicted = new double[p.size()];
		double[] observed = new double[p.size()];
		int used = 0;
		for (int i = 0; i < p.size(); i++) {
			int ion = p.ionNumber[i];
			if (p.intensity[i] <= 0 || ion < 1 || ion >= n || p.charge[i] > Math.max(1, precursorCharge)) {
				continue;
			}
			double mz = fragmentMz(f.backbone, ion, p.yIon[i], p.charge[i], f.shiftSite, f.shift);
			double tol = mz * tolerancePpm * 1e-6;
			predicted[used] = p.intensity[i];
			observed[used] = spectrum.maxIntensity(mz - tol, mz + tol);
			used++;
		}
		if (used == 0) {
			return Double.NaN;
		}
		return entropySimilarity(predicted, observed, used);
	}

	/**
	 * Unweighted spectral-entropy similarity over aligned intensity vectors:
	 * {@code 1 - (2 H(P+Q) - H(P) - H(Q)) / ln 4} with each vector summing to 1.
	 * 1 for identical shapes, 0 for disjoint ones or an empty observation.
	 */
	static double entropySimilarity(double[] p, double[] q, int n) {
		double sp = 0;
		double sq = 0;
		for (int i = 0; i < n; i++) {
			sp += p[i];
			sq += q[i];
		}
		if (sp <= 0 || sq <= 0) {
			return 0.0;
		}
		double hp = 0;
		double hq = 0;
		double hm = 0;
		for (int i = 0; i < n; i++) {
			double a = p[i] / sp;
			double b = q[i] / sq;
			double m = (a + b) / 2;
			hp -= xlogx(a);
			hq -= xlogx(b);
			hm -= xlogx(m);
		}
		double sim = 1.0 - (2 * hm - hp - hq) / LN4;
		return Math.max(0.0, Math.min(1.0, sim));
	}

	private static double xlogx(double x) {
		return x > 0 ? x * Math.log(x) : 0.0;
	}

	/**
	 * Scores every form and sets the candidate's MS2 fields. The verdict
	 * compares the forms at the candidate's own site - its substitution, the
	 * modification there and any other same-mass substitution there: if the top
	 * two are within {@code margin}, AMBIGUOUS; otherwise the top one. The best
	 * form at any other site only sets the elsewhere score and flag.
	 */
	static void rank(VariantCandidate c, List<Form> forms, Predictions predictions, MzmlScanReader.Peaks spectrum,
			double tolerancePpm, double margin) {
		int site = c.position - 1;
		double own = Double.NaN;
		double other = Double.NaN;
		double mod = Double.NaN;
		double elsewhere = Double.NaN;
		String bestLabel = "";
		double best = Double.NEGATIVE_INFINITY;
		for (Form f : forms) {
			double s = score(f, predictions.get(f.backbone, c.charge), spectrum, c.charge, tolerancePpm);
			if (Double.isNaN(s)) {
				continue;
			}
			if (f.kind == Kind.SUBSTITUTION) {
				own = s;
			} else if (f.site != site) {
				elsewhere = max(elsewhere, s);
			} else if (f.kind == Kind.OTHER_SUBSTITUTION) {
				other = max(other, s);
			} else {
				mod = max(mod, s);
			}
			if (s > best) {
				best = s;
				bestLabel = f.label;
			}
		}
		c.ms2ScoreSubstitution = own;
		c.ms2ScoreOtherSubstitution = other;
		c.ms2ScoreModification = mod;
		c.ms2ScoreElsewhere = elsewhere;
		c.ms2BestForm = bestLabel;
		if (Double.isNaN(own) || Double.isNaN(mod)) {
			c.ms2Verdict = VariantCandidate.Ms2Verdict.NOT_RANKED;
			c.ms2LocalisesElsewhere = null;
			return;
		}
		double sameSite = Math.max(own, Math.max(mod, Double.isNaN(other) ? Double.NEGATIVE_INFINITY : other));
		c.ms2LocalisesElsewhere = Double.isNaN(elsewhere) ? null : elsewhere - sameSite > margin;
		double[] scores = { own, Double.isNaN(other) ? Double.NEGATIVE_INFINITY : other, mod };
		int top = 0;
		for (int i = 1; i < 3; i++) {
			if (scores[i] > scores[top]) {
				top = i;
			}
		}
		double runnerUp = Double.NEGATIVE_INFINITY;
		for (int i = 0; i < 3; i++) {
			if (i != top) {
				runnerUp = Math.max(runnerUp, scores[i]);
			}
		}
		if (scores[top] - runnerUp <= margin) {
			c.ms2Verdict = VariantCandidate.Ms2Verdict.AMBIGUOUS;
		} else if (top == 0) {
			c.ms2Verdict = VariantCandidate.Ms2Verdict.SUBSTITUTION;
		} else if (top == 1) {
			c.ms2Verdict = VariantCandidate.Ms2Verdict.OTHER_SUBSTITUTION;
		} else {
			c.ms2Verdict = VariantCandidate.Ms2Verdict.MODIFICATION;
		}
	}

	private static double max(double current, double s) {
		return Double.isNaN(current) ? s : Math.max(current, s);
	}
}
