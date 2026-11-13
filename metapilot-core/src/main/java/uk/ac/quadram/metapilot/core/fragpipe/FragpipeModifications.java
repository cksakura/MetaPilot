package uk.ac.quadram.metapilot.core.fragpipe;

import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

/**
 * Writes MetaPilot's fixed and variable modifications into a FragPipe workflow's
 * {@code msfragger.table.fix-mods} and {@code msfragger.table.var-mods}.
 * <p>
 * Modifications are given as strings: fixed {@code site:mass} (for example
 * {@code C:57.02146} or {@code protein n-term:42.0106}) and variable
 * {@code name:sites:mass} (for example {@code Oxidation:M:15.9949}). Sites are
 * residue letters or one of {@code peptide n-term}, {@code peptide c-term},
 * {@code protein n-term} and {@code protein c-term}; MSFragger's own site
 * syntax ({@code n^}, {@code [^}, {@code nQnC}, ...) is passed through.
 */
public final class FragpipeModifications {

	/** Carbamidomethyl on C. */
	public static final List<String> DEFAULT_FIXED = List.of("C:57.02146");

	/** Oxidation on M and acetylation at the protein N-terminus. */
	public static final List<String> DEFAULT_VARIABLE = List.of("Oxidation:M:15.9949",
			"Acetylation:protein n-term:42.0106");

	/** FragPipe's variable modification table has this many rows. */
	static final int VARIABLE_SLOTS = 16;

	/** FragPipe's fixed modification table rows, in its order. */
	private static final String[] FIXED_SLOTS = { "C-Term Peptide", "N-Term Peptide", "C-Term Protein",
			"N-Term Protein", "G (glycine)", "A (alanine)", "S (serine)", "P (proline)", "V (valine)",
			"T (threonine)", "C (cysteine)", "L (leucine)", "I (isoleucine)", "N (asparagine)",
			"D (aspartic acid)", "Q (glutamine)", "K (lysine)", "E (glutamic acid)", "M (methionine)",
			"H (histidine)", "F (phenylalanine)", "R (arginine)", "Y (tyrosine)", "W (tryptophan)", "B ", "J", "O",
			"U", "X", "Z" };

	private static final Map<String, String> TERMINUS_FIXED_SLOT = Map.of("peptide c-term", "C-Term Peptide",
			"peptide n-term", "N-Term Peptide", "protein c-term", "C-Term Protein", "protein n-term",
			"N-Term Protein");

	private static final Map<String, String> TERMINUS_SITE_CODE = Map.of("peptide n-term", "n^", "peptide c-term",
			"c^", "protein n-term", "[^", "protein c-term", "]^");

	private FragpipeModifications() {
	}

	/**
	 * Replaces both modification tables of {@code workflow}.
	 *
	 * @throws IllegalArgumentException for a site FragPipe has no slot for, a
	 *                                  malformed entry, or more variable
	 *                                  modifications than FragPipe's table holds
	 */
	public static FragpipeWorkflowFile apply(FragpipeWorkflowFile workflow, List<String> fixed,
			List<String> variable) {
		return workflow.set("msfragger.table.fix-mods", fixedTable(fixed))
				.set("msfragger.table.var-mods", variableTable(variable));
	}

	static String fixedTable(List<String> fixed) {
		Map<String, Double> massBySlot = new LinkedHashMap<>();
		for (String slot : FIXED_SLOTS) {
			massBySlot.put(slot, 0.0);
		}
		for (String entry : fixed == null ? List.<String>of() : fixed) {
			int colon = entry.lastIndexOf(':');
			if (colon < 1) {
				throw new IllegalArgumentException("Fixed modification is not 'site:mass': " + entry);
			}
			String slot = fixedSlot(entry.substring(0, colon).trim());
			massBySlot.put(slot, massBySlot.get(slot) + parseMass(entry, entry.substring(colon + 1)));
		}
		List<String> rows = new ArrayList<>();
		for (Map.Entry<String, Double> e : massBySlot.entrySet()) {
			rows.add(formatMass(e.getValue()) + "," + e.getKey() + ",true,-1");
		}
		return String.join("; ", rows);
	}

	static String variableTable(List<String> variable) {
		List<String> entries = variable == null ? List.of() : variable;
		if (entries.size() > VARIABLE_SLOTS) {
			throw new IllegalArgumentException("FragPipe accepts at most " + VARIABLE_SLOTS
					+ " variable modifications, got " + entries.size());
		}
		List<String> rows = new ArrayList<>();
		for (String entry : entries) {
			int massColon = entry.lastIndexOf(':');
			int siteColon = massColon < 0 ? -1 : entry.lastIndexOf(':', massColon - 1);
			if (siteColon < 0) {
				throw new IllegalArgumentException("Variable modification is not 'name:sites:mass': " + entry);
			}
			String sites = siteCode(entry.substring(siteColon + 1, massColon).trim(), entry);
			double mass = parseMass(entry, entry.substring(massColon + 1));
			int maxPerPeptide = sites.contains("^") ? 1 : 3;
			rows.add(formatMass(mass) + "," + sites + ",true," + maxPerPeptide);
		}
		for (int slot = rows.size() + 1; slot <= VARIABLE_SLOTS; slot++) {
			rows.add("0.0,site_" + slot + ",false,1");
		}
		return String.join("; ", rows);
	}

	private static String fixedSlot(String site) {
		String terminus = TERMINUS_FIXED_SLOT.get(site.toLowerCase(Locale.ROOT));
		if (terminus != null) {
			return terminus;
		}
		if (site.length() == 1) {
			String letter = site.toUpperCase(Locale.ROOT);
			for (String slot : FIXED_SLOTS) {
				if (slot.trim().equals(letter) || slot.startsWith(letter + " (")) {
					return slot;
				}
			}
		}
		throw new IllegalArgumentException("FragPipe has no fixed modification slot for site '" + site + "'");
	}

	private static String siteCode(String sites, String entry) {
		String terminus = TERMINUS_SITE_CODE.get(sites.toLowerCase(Locale.ROOT));
		if (terminus != null) {
			return terminus;
		}
		String code = sites.replace(" ", "");
		if (code.isEmpty() || !code.matches("[A-Za-z\\[\\]^*]+")) {
			throw new IllegalArgumentException("Variable modification has no usable sites: " + entry);
		}
		// lower-case n and c are MSFragger's terminus markers, so case is kept
		return code;
	}

	private static double parseMass(String entry, String mass) {
		try {
			return Double.parseDouble(mass.trim());
		} catch (NumberFormatException e) {
			throw new IllegalArgumentException("Modification mass is not a number: " + entry, e);
		}
	}

	private static String formatMass(double mass) {
		return mass == 0.0 ? "0.0" : BigDecimal.valueOf(mass).stripTrailingZeros().toPlainString();
	}
}
