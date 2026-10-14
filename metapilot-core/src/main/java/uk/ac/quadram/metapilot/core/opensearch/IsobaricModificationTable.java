package uk.ac.quadram.metapilot.core.opensearch;

import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import java.util.Set;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Which single-residue substitutions have exactly the mass of a modification on
 * the <em>same residue</em>, and so cannot be told apart from it by backbone
 * fragments.
 * <p>
 * Same residue is the whole point. A substitution whose mass matches a
 * modification that sits on a different residue is separable by localization:
 * on 12mix A-&gt;S and oxidation were, because oxidation sits on M. A
 * substitution that acts on the residue the modification acts on is not: the
 * shifted fragment series is identical either way.
 * <p>
 * Built from the bundled modification reference ({@link ModificationMassDatabase})
 * crossed with each entry's stated specificity, plus Unimod:34 methylation,
 * which the bundled MaxQuant set only lists on K, R and E. That explicit entry
 * is what makes S-&gt;T, D-&gt;E and N-&gt;Q isobaric with methylation of the same
 * residue; G-&gt;A and V-&gt;I/L are not, because Unimod does not list methylation
 * on G or V.
 * <p>
 * Shared by {@link SequenceVariantCaller} (its {@code isobaric_modifications}
 * column) and the pangenome variant check (its {@code ms2_decidable} flag), so
 * the two can never disagree about which substitutions are undecidable.
 */
public final class IsobaricModificationTable {

	private static final double UNIMOD_METHYL_DA = 14.01565;
	private static final String UNIMOD_METHYL_SITES = "CDEHIKLNQRST";

	/** The 20 residues, I and L both present; their equal mass is skipped below. */
	private static final String RESIDUES = "ACDEFGHIKLMNPQRSTVWY";

	private static final Pattern SITE_SPEC = Pattern.compile("\\(([^()]*)\\)\\s*$");

	private IsobaricModificationTable() {
	}

	/**
	 * Modifications on {@code residue} whose mass is within {@code toleranceDa} of
	 * {@code substitutionMassDa}, sorted; empty when the substitution is
	 * decidable.
	 *
	 * @param residue      the residue the substitution acts on (the assigned one)
	 * @param peptideNTerm whether the residue is the peptide's first, which admits
	 *                     peptide N-terminal modifications
	 */
	public static Set<String> sameResidue(char residue, double substitutionMassDa, boolean peptideNTerm,
			double toleranceDa) {
		Set<String> out = new TreeSet<>();
		for (ModificationMassDatabase.Match m : ModificationMassDatabase.matching(substitutionMassDa, toleranceDa)) {
			if (appliesTo(m.entry.title, residue, peptideNTerm)) {
				out.add(m.entry.title);
			}
		}
		if (Math.abs(substitutionMassDa - UNIMOD_METHYL_DA) <= toleranceDa && UNIMOD_METHYL_SITES.indexOf(residue) >= 0) {
			out.add("Methyl (" + residue + ") [Unimod:34]");
		}
		return out;
	}

	/**
	 * {@code mass(to) - mass(from)} for underivatised residues, with
	 * {@code fixedCysteineModDa} added to C. An MSFragger open search carries
	 * carbamidomethyl as a fixed mod ({@link SubstitutionMassTable#FIXED_CYSTEINE_MOD});
	 * a pFind open search writes no fixed mods, so C is unmodified there.
	 */
	public static double substitutionMass(char from, char to, double fixedCysteineModDa) {
		return residueMass(to, fixedCysteineModDa) - residueMass(from, fixedCysteineModDa);
	}

	private static double residueMass(char aa, double fixedCysteineModDa) {
		char u = Character.toUpperCase(aa);
		double m = SubstitutionMassTable.residueMass(u);
		return u == 'C' ? m - SubstitutionMassTable.FIXED_CYSTEINE_MOD + fixedCysteineModDa : m;
	}

	/**
	 * Every undecidable ordered substitution, one line each, for logging at the
	 * start of a run. K and R are left out (they change the digest and are never
	 * candidates) and so is I&lt;-&gt;L (no mass difference). Substitutions on a
	 * peptide's first residue can additionally collide with N-terminal
	 * modifications; those are listed with an {@code [N-term]} tag.
	 */
	public static List<String> describe(double toleranceDa, double fixedCysteineModDa) {
		List<String> out = new ArrayList<>();
		for (int i = 0; i < RESIDUES.length(); i++) {
			for (int j = 0; j < RESIDUES.length(); j++) {
				char from = RESIDUES.charAt(i);
				char to = RESIDUES.charAt(j);
				if (from == to || from == 'K' || from == 'R' || to == 'K' || to == 'R'
						|| (isIL(from) && isIL(to))) {
					continue;
				}
				double delta = substitutionMass(from, to, fixedCysteineModDa);
				Set<String> any = sameResidue(from, delta, false, toleranceDa);
				Set<String> nTermOnly = sameResidue(from, delta, true, toleranceDa);
				nTermOnly.removeAll(any);
				if (any.isEmpty() && nTermOnly.isEmpty()) {
					continue;
				}
				StringBuilder sb = new StringBuilder(
						String.format(Locale.ROOT, "%c->%c %+.5f Da: ", from, to, delta));
				sb.append(String.join("; ", any));
				if (!nTermOnly.isEmpty()) {
					sb.append(any.isEmpty() ? "" : "; ").append("[N-term] ").append(String.join("; ", nTermOnly));
				}
				out.add(sb.toString());
			}
		}
		return out;
	}

	/** MaxQuant titles carry their specificity, e.g. {@code Methyl (KR)}, {@code Acetyl (Protein N-term)}. */
	static boolean appliesTo(String title, char residue, boolean peptideNTerm) {
		if (title == null) {
			return false;
		}
		Matcher m = SITE_SPEC.matcher(title);
		if (!m.find()) {
			return false;
		}
		String spec = m.group(1).trim();
		if (spec.contains("N-term")) {
			return peptideNTerm && !spec.startsWith("Protein");
		}
		return spec.matches("[A-Z]+") && spec.indexOf(residue) >= 0;
	}

	private static boolean isIL(char c) {
		return c == 'I' || c == 'L';
	}
}
