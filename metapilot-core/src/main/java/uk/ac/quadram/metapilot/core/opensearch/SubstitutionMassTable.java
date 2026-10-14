package uk.ac.quadram.metapilot.core.opensearch;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.TreeMap;


/**
 * Every single amino-acid substitution, indexed by the mass difference it
 * produces.
 *
 * <h2>Why this exists</h2>
 * An open search assigns peptide P with a delta mass when the spectrum actually
 * came from P', a peptide differing from P by one residue. The engine has no way
 * to prefer P' unless P' is in the database and outscores the shifted
 * assignment, and in a MAG metaproteomics search it frequently is not: genome
 * selection picks one strain, the peptide came from a relative, and the two
 * differ by a single residue.
 * <p>
 * Measured on 12mix (2026-09-14), the two peaks at +/-14.0157 Da - the mass of
 * CH2, i.e. G/A, S/T, V/L, V/I, N/Q, D/E - were the two worst accepted peaks in
 * the entire run: 0.94% and 0.96% decoy FDR after the per-peak local-FDR walk,
 * against 18.35% and 17.85% <em>true</em> error measured by entrapment. The
 * third worst was +28.0312 (A-&gt;V, C-&gt;M, i.e. two CH2) at 14.10%. Everything
 * else accepted sat at or below 4.56%.
 * <p>
 * Target-decoy cannot see this error class at all, because the wrong peptide is
 * itself a target sequence from a real protein: a reversed or shuffled decoy
 * never generates that competitor. Neither does the delta-mass fragment-gain
 * test in {@link DeltaMassEvidenceGate}, because a substituted residue genuinely
 * shifts part of the ion series - the evidence is real, the assignment is wrong.
 * This table is the only handle on it.
 *
 * <h2>Residue masses, and why they are not taken from AminoAcidProperty</h2>
 * The obvious source is {@link AminoAcidProperty#AMINOACID_FORMULAR}, and the
 * first version of this class used it. It is wrong here. That table encodes one
 * particular sample-preparation chemistry - its own comment says "C is the
 * alkylate form and K is the guanidylate form" - so its K is guanidylated
 * lysine at 170.11676, not lysine at 128.09496. Built from it, this table
 * silently produced {@code R->K} and {@code K->R} as +/-14.01565 Da
 * substitutions, because guanidyl-K and R happen to differ by exactly CH2. Every
 * PSM on the CH2 peaks would then have been offered an R/K candidate that cannot
 * occur in a tryptic search, and a candidate at a K or R is precisely where a
 * tryptic peptide's cleavage site is.
 * <p>
 * So the masses are stated here explicitly, for underivatised residues, and the
 * fixed cysteine modification is added separately and named. That makes this
 * table's chemistry a property of the search it describes rather than an
 * inherited assumption: the workflow's searches set
 * {@code add_C_cysteine = 57.02146} and do not guanidylate, so C carries the fixed
 * mod and K does not.
 *
 * <h2>What is deliberately excluded</h2>
 * <ul>
 * <li><b>I and L</b> - identical residue mass. The substitution is real and
 * common but produces a delta of exactly 0, so it is invisible here and is not a
 * delta-mass artifact. Pairs involving only I/L interchange are dropped.</li>
 * <li><b>B, J, O, U, X, Z</b> - not residues; their formula rows are all
 * zero.</li>
 * </ul>
 *
 * <h2>Collisions are expected and must not be resolved here</h2>
 * Several substitution masses coincide exactly with real modifications:
 * N-&gt;D is 0.98401 Da against deamidation at 0.98402; A-&gt;S and F-&gt;Y are
 * 15.99491 against oxidation; A-&gt;Q and G-&gt;N are 57.02146 against
 * carbamidomethyl; G-&gt;P is 40.03130; A-&gt;P is 26.01565. On the same run
 * those four peaks carried 22,186 / 4,790 / 1,913 / 6,490 PSMs at 2.80% / 1.68%
 * / 1.17% / 1.76% entrapment - i.e. they are overwhelmingly the modification,
 * not the substitution. A mass match is therefore a <em>flag</em>, never a
 * verdict; {@link SubstitutionReassigner} resolves it with the localization
 * evidence, which is what separates "the delta sits on an M" (oxidation) from
 * "the delta sits on an A" (A-&gt;S).
 */
public final class SubstitutionMassTable {

	/** One ordered substitution: the assigned residue became the true one. */
	public static final class Substitution {
		/** Residue as the search assigned it. */
		public final char from;
		/** Residue the spectrum actually carried. */
		public final char to;
		/** {@code mass(to) - mass(from)}; the delta the open search reports. */
		public final double massDelta;

		Substitution(char from, char to, double massDelta) {
			this.from = from;
			this.to = to;
			this.massDelta = massDelta;
		}

		@Override
		public String toString() {
			return String.format(Locale.ROOT, "%c->%c(%+.5f)", from, to, massDelta);
		}
	}

	/** The 20 residues. I and L are both present; their equal mass is handled below. */
	private static final String RESIDUES = "ACDEFGHIKLMNPQRSTVWY";

	/**
	 * Fixed carbamidomethyl on cysteine, matching {@code add_C_cysteine} in this
	 * workflow's {@code fragger_open.params}. Stated rather than assumed: a search
	 * without alkylation would need this at zero, and a substitution to or from C
	 * would then be a different mass.
	 */
	public static final double FIXED_CYSTEINE_MOD = 57.02146;

	/** Underivatised monoisotopic residue masses. See the class docs. */
	private static final double[] RESIDUE_MASS = new double[128];
	/** Ordered by mass delta, so a tolerance lookup is a sub-map. */
	private static final TreeMap<Double, List<Substitution>> BY_MASS = new TreeMap<>();

	static {
		RESIDUE_MASS['G'] = 57.02146;
		RESIDUE_MASS['A'] = 71.03711;
		RESIDUE_MASS['S'] = 87.03203;
		RESIDUE_MASS['P'] = 97.05276;
		RESIDUE_MASS['V'] = 99.06841;
		RESIDUE_MASS['T'] = 101.04768;
		RESIDUE_MASS['C'] = 103.00919 + FIXED_CYSTEINE_MOD;
		RESIDUE_MASS['L'] = 113.08406;
		RESIDUE_MASS['I'] = 113.08406;
		RESIDUE_MASS['N'] = 114.04293;
		RESIDUE_MASS['D'] = 115.02694;
		RESIDUE_MASS['Q'] = 128.05858;
		RESIDUE_MASS['K'] = 128.09496;
		RESIDUE_MASS['E'] = 129.04259;
		RESIDUE_MASS['M'] = 131.04049;
		RESIDUE_MASS['H'] = 137.05891;
		RESIDUE_MASS['F'] = 147.06841;
		RESIDUE_MASS['R'] = 156.10111;
		RESIDUE_MASS['Y'] = 163.06333;
		RESIDUE_MASS['W'] = 186.07931;
		for (int i = 0; i < RESIDUES.length(); i++) {
			for (int j = 0; j < RESIDUES.length(); j++) {
				if (i == j) {
					continue;
				}
				char from = RESIDUES.charAt(i);
				char to = RESIDUES.charAt(j);
				double delta = RESIDUE_MASS[to] - RESIDUE_MASS[from];
				// I/L and any other exact-mass pair: no delta, nothing to explain.
				if (Math.abs(delta) < 1e-6) {
					continue;
				}
				// A substitution that creates or destroys a tryptic cleavage site
				// does not yield the same peptide with one residue changed - it
				// yields a different digest entirely, with different termini. Those
				// are real events, but they are not what this table models, and
				// proposing one would generate a candidate the enzyme could never
				// have produced.
				if (from == 'K' || from == 'R' || to == 'K' || to == 'R') {
					continue;
				}
				BY_MASS.computeIfAbsent(round(delta), k -> new ArrayList<>())
						.add(new Substitution(from, to, delta));
			}
		}
	}

	private SubstitutionMassTable() {
	}

	/** Monoisotopic residue mass; C carries the fixed carbamidomethyl. */
	public static double residueMass(char aa) {
		char u = Character.toUpperCase(aa);
		return u < RESIDUE_MASS.length ? RESIDUE_MASS[u] : 0.0;
	}

	/**
	 * Every substitution whose mass difference is within {@code toleranceDa} of
	 * {@code delta}, or an empty list.
	 * <p>
	 * The tolerance should be tight. Fitted peaks in this workflow have sigma
	 * around 2 mDa, so 5 mDa admits the peak and little else; opening it to
	 * 20 mDa starts merging distinct substitutions and pulling in unrelated
	 * modifications.
	 */
	public static List<Substitution> matching(double delta, double toleranceDa) {
		List<Substitution> out = new ArrayList<>();
		for (Map.Entry<Double, List<Substitution>> e : BY_MASS
				.subMap(round(delta - toleranceDa), true, round(delta + toleranceDa), true).entrySet()) {
			out.addAll(e.getValue());
		}
		return out;
	}

	/** True when any substitution explains this delta. */
	public static boolean isSubstitutionMass(double delta, double toleranceDa) {
		return !BY_MASS.subMap(round(delta - toleranceDa), true, round(delta + toleranceDa), true).isEmpty();
	}

	/** Distinct substitution masses in the table; diagnostics only. */
	public static List<Double> masses() {
		return Collections.unmodifiableList(new ArrayList<>(BY_MASS.keySet()));
	}

	private static double round(double v) {
		return Math.round(v * 1e6) / 1e6;
	}
}
