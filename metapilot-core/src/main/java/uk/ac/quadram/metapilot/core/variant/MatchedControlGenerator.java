package uk.ac.quadram.metapilot.core.variant;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.TreeSet;

/**
 * Matched random-substitution controls: the background a candidate's pangenome
 * hit rate is compared against.
 * <p>
 * For each real candidate X-&gt;Y, pick a random identified base peptide of the
 * same run and engine that contains X, pick a random X in it, and apply the
 * same X-&gt;Y. The control has the substitution type and the peptide background
 * of the real candidates but no spectrum behind it, so how often it turns up in
 * UHGP-100 is what chance alone produces for that substitution. I and L count
 * as one residue when looking for X, as everywhere else in the check.
 * <p>
 * Reproducible: the pool is sorted, the candidates are visited in the order
 * given, and every choice comes from one {@link Random} seeded with
 * {@code pangenome_control_seed}.
 */
final class MatchedControlGenerator {

	/** One control peptide, tied to the candidate it was drawn for. */
	static final class Control {
		final VariantCandidate parent;
		final String basePeptide;
		/** 1-based. */
		final int position;
		final String peptide;

		Control(VariantCandidate parent, String basePeptide, int position, String peptide) {
			this.parent = parent;
			this.basePeptide = basePeptide;
			this.position = position;
			this.peptide = peptide;
		}
	}

	private MatchedControlGenerator() {
	}

	/**
	 * @param candidates   in a deterministic order
	 * @param basePeptides identified base peptides of the same run and engine
	 * @return up to {@code perCandidate} controls per candidate; none for a
	 *         candidate whose from-residue no identified peptide contains
	 */
	static List<Control> generate(List<VariantCandidate> candidates, Collection<String> basePeptides, long seed,
			int perCandidate) {
		List<Control> out = new ArrayList<>();
		if (perCandidate <= 0 || candidates.isEmpty()) {
			return out;
		}
		List<String> pool = new ArrayList<>(new TreeSet<>(basePeptides));
		Map<Character, List<String>> byResidue = new HashMap<>();
		Random rnd = new Random(seed);
		for (VariantCandidate c : candidates) {
			char from = fold(c.fromAa);
			List<String> containing = byResidue.computeIfAbsent(from, r -> containing(pool, r));
			if (containing.isEmpty()) {
				continue;
			}
			for (int k = 0; k < perCandidate; k++) {
				String base = containing.get(rnd.nextInt(containing.size()));
				List<Integer> sites = sites(base, from);
				int site = sites.get(rnd.nextInt(sites.size()));
				StringBuilder sb = new StringBuilder(base);
				sb.setCharAt(site, c.toAa);
				out.add(new Control(c, base, site + 1, sb.toString()));
			}
		}
		return out;
	}

	private static List<String> containing(List<String> pool, char foldedResidue) {
		List<String> out = new ArrayList<>();
		for (String p : pool) {
			if (!sites(p, foldedResidue).isEmpty()) {
				out.add(p);
			}
		}
		return out;
	}

	private static List<Integer> sites(String peptide, char foldedResidue) {
		List<Integer> out = new ArrayList<>();
		for (int i = 0; i < peptide.length(); i++) {
			if (fold(peptide.charAt(i)) == foldedResidue) {
				out.add(i);
			}
		}
		return out;
	}

	private static char fold(char aa) {
		char u = Character.toUpperCase(aa);
		return u == 'I' ? 'L' : u;
	}
}
