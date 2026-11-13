package uk.ac.quadram.metapilot.dbSearch.pfind;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

/**
 * pFind's enzyme vocabulary: the names the pFind workflows offer and the
 * definition each one writes to the {@code enzyme=} line of {@code pFind.cfg}.
 * <p>
 * A definition is {@code name sites blocker terminus}, where blocker {@code _}
 * means no residue blocks cleavage. The names are pFind's own and the
 * {@code _P} suffix reads the opposite of the convention elsewhere in
 * MetaPilot: pFind's {@code Trypsin} ({@code KR _ C}) cleaves before proline,
 * and {@code Trypsin_P} ({@code KR P C}) does not.
 */
public final class PFindEnzymes {

	/**
	 * The default, and the enzyme every pFind search used before the selected
	 * enzyme was written to the configuration: pFind's own template carries
	 * {@code enzyme=Trypsin KR _ C}.
	 */
	public static final String DEFAULT_ENZYME = "Trypsin";

	/** Explains the naming where the enzyme is chosen. */
	public static final String TOOLTIP = "pFind enzyme names: a name ending in \"_P\" does not cleave before proline. "
			+ "\"Trypsin\" cleaves after K/R also before proline (called Trypsin/P elsewhere); "
			+ "\"Trypsin_P\" is conventional trypsin.";

	private static final Map<String, String> DEFINITIONS = new LinkedHashMap<>();

	static {
		add("Trypsin KR _ C");
		add("Trypsin_P KR P C");
		add("Asp-N D _ N");
		add("Arg-C R P C");
		add("Chymotrypsin FYWML _ C");
		add("Chymotrypsin_P FYWML P C");
		add("CNBr M _ C");
		add("FormicAcid D _ C");
		add("Lys-K_P K P C");
		add("Lys-K K _ C");
		add("PepsinA-FL FL _ C");
		add("TrypCNBr-KRM KRM P C");
		add("TrypChymo-FYWLKR FYWLKR P C");
		add("V8-BDEZ_P BDEZ P C");
		add("V8-EZ_P EZ P C");
	}

	private static void add(String definition) {
		DEFINITIONS.put(definition.substring(0, definition.indexOf(' ')), definition);
	}

	private PFindEnzymes() {
	}

	/** Enzyme names in display order, the default first. */
	public static List<String> names() {
		return List.copyOf(DEFINITIONS.keySet());
	}

	/**
	 * The {@code pFind.cfg} definition for a pFind enzyme name, matched
	 * case-insensitively.
	 *
	 * @throws IllegalArgumentException if the name is not a pFind enzyme; the
	 *                                  search must not fall back to another digest
	 */
	public static String definition(String name) {
		if (name != null) {
			String wanted = name.trim().toLowerCase(Locale.ROOT);
			for (Map.Entry<String, String> e : DEFINITIONS.entrySet()) {
				if (e.getKey().toLowerCase(Locale.ROOT).equals(wanted)) {
					return e.getValue();
				}
			}
		}
		throw new IllegalArgumentException(
				"Unknown pFind enzyme '" + name + "'. pFind enzymes are: " + String.join(", ", names()));
	}
}
