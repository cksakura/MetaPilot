package uk.ac.quadram.metapilot.core.opensearch;

import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.dom4j.Document;
import org.dom4j.DocumentException;
import org.dom4j.Element;
import org.dom4j.io.SAXReader;

import uk.ac.quadram.metapilot.core.aminoacid.AminoAcidProperty;

/**
 * Chemical-modification mass reference for explaining an open search's
 * accepted delta-mass peaks, built from the bundled {@code modifications.xml}
 * (586 {@code <modification>} entries, MaxQuant/Unimod-derived).
 * <p>
 * This is the local set Kai Cheng chose (2026-09-22) over fetching the full
 * ~2,500-entry Unimod database for this stage: "Use the bundled 586-entry set
 * in current stage." No parser for this file existed anywhere in the codebase
 * before this class - {@link uk.ac.quadram.metapilot.core.mod.UnimodReader}
 * parses a different, larger schema ({@code elements}/{@code modifications}/
 * {@code amino_acids}/{@code mod_bricks} with the mass stated directly on each
 * {@code <delta>}) that no bundled resource in this repository actually
 * supplies.
 * <p>
 * Each entry here states only an elemental composition (e.g.
 * {@code "C(2) H(2) O"}), never a mass attribute, so the mass is computed from
 * standard monoisotopic element masses - the same constants
 * {@link AminoAcidProperty} and {@link SubstitutionMassTable} already use, so
 * this table cannot silently disagree with the rest of the workflow about what
 * an atom weighs.
 *
 * <h2>Element vocabulary</h2>
 * Checked against every {@code composition} attribute the bundled file
 * contains (not just the top-level 586 - including the {@code diagnostic} and
 * {@code neutralloss} sub-elements, a superset): the only tokens used are
 * {@code C, H, N, O, S, P, Na}, plus MaxQuant's heavy-isotope labels
 * {@code Cx, Hx, Nx, Ox} (13C, 2H, 15N, 18O respectively - e.g. {@code Pro5}'s
 * {@code "Cx(5) C(-5)"} for {@code Proline-13C5}). An unrecognised element
 * token would mean the file grew a new one outside this checked set; that
 * entry is skipped rather than silently mis-massed, and the skip count is
 * logged.
 */
public final class ModificationMassDatabase {

	private static final Logger LOGGER = LogManager.getLogger(ModificationMassDatabase.class);

	/** Classpath location of the bundled reference (root of {@code src/main/resources}). */
	private static final String RESOURCE = "/modifications.xml";

	/** One {@code <modification>} entry: its name and monoisotopic mass. */
	public static final class Entry {
		/** Already site-qualified by MaxQuant's own convention, e.g. "Acetyl (K)". */
		public final String title;
		public final String description;
		public final String composition;
		public final double monoMass;

		Entry(String title, String description, String composition, double monoMass) {
			this.title = title;
			this.description = description == null ? "" : description;
			this.composition = composition;
			this.monoMass = monoMass;
		}
	}

	/** One candidate match: the entry, and how far its mass sits from the peak. */
	public static final class Match {
		public final Entry entry;
		/** {@code delta - entry.monoMass}; the residual the peak leaves unexplained. */
		public final double massErrorDa;

		Match(Entry entry, double massErrorDa) {
			this.entry = entry;
			this.massErrorDa = massErrorDa;
		}
	}

	private static final Map<String, Double> ELEMENT_MASS = new HashMap<>();
	static {
		ELEMENT_MASS.put("C", AminoAcidProperty.MONOW_C);
		ELEMENT_MASS.put("H", AminoAcidProperty.MONOW_H);
		ELEMENT_MASS.put("N", AminoAcidProperty.MONOW_N);
		ELEMENT_MASS.put("O", AminoAcidProperty.MONOW_O);
		ELEMENT_MASS.put("S", AminoAcidProperty.MONOW_S32);
		ELEMENT_MASS.put("P", AminoAcidProperty.MONOW_P);
		// Not in AminoAcidProperty - stated here explicitly. Monoisotopic 23Na.
		ELEMENT_MASS.put("Na", 22.98976928);
		// MaxQuant's heavy-isotope element names: 13C, 2H (deuterium), 15N, 18O.
		ELEMENT_MASS.put("Cx", AminoAcidProperty.MONOW_C13);
		ELEMENT_MASS.put("Hx", AminoAcidProperty.MONOW_D);
		ELEMENT_MASS.put("Nx", AminoAcidProperty.MONOW_N15);
		ELEMENT_MASS.put("Ox", AminoAcidProperty.MONOW_O18);
	}

	/** e.g. "C(2)", "O", "Cx(-5)" - a bare symbol with no parenthesised count means count 1. */
	private static final Pattern TOKEN = Pattern.compile("([A-Za-z]+)(?:\\((-?\\d+)\\))?");

	/** Ordered by mass, so a tolerance lookup is a sub-map - same idiom as {@link SubstitutionMassTable#BY_MASS}. */
	private static final TreeMap<Double, List<Entry>> BY_MASS = new TreeMap<>();
	private static final List<Entry> ALL = new ArrayList<>();
	private static boolean loaded = false;

	private ModificationMassDatabase() {
	}

	private static synchronized void ensureLoaded() {
		if (loaded) {
			return;
		}
		// Set first: a parse failure must not be retried on every peak lookup for
		// the rest of the run - fall back to reporting no candidates instead.
		loaded = true;
		try (InputStream in = ModificationMassDatabase.class.getResourceAsStream(RESOURCE)) {
			if (in == null) {
				LOGGER.warn("Modification reference {} not found on the classpath; potential-modification "
						+ "matching will report no candidates for every peak.", RESOURCE);
				return;
			}
			parse(in);
		} catch (DocumentException | IOException e) {
			LOGGER.error("Failed to parse modification reference " + RESOURCE
					+ "; potential-modification matching will report no candidates for every peak.", e);
		}
	}

	@SuppressWarnings("unchecked")
	private static void parse(InputStream in) throws DocumentException {
		SAXReader reader = new SAXReader();
		Document document = reader.read(in);
		Element root = document.getRootElement();

		int parsed = 0;
		int skipped = 0;
		List<Element> mods = root.elements("modification");
		for (Element mod : mods) {
			String title = mod.attributeValue("title");
			String description = mod.attributeValue("description");
			String composition = mod.attributeValue("composition");
			Double mass = title == null || composition == null ? null : massOf(composition);
			if (mass == null) {
				skipped++;
				continue;
			}
			Entry entry = new Entry(title, description, composition, mass);
			ALL.add(entry);
			BY_MASS.computeIfAbsent(round(mass), k -> new ArrayList<>()).add(entry);
			parsed++;
		}
		LOGGER.info("Modification reference: {} of {} entries in {} parsed ({} skipped - unrecognised "
				+ "composition element or missing attribute)", parsed, parsed + skipped, RESOURCE, skipped);
	}

	/**
	 * Monoisotopic mass from a composition string like {@code "C(2) H(2) O"} or
	 * {@code "Cx(5) C(-5)"}, or null if any element token is not in
	 * {@link #ELEMENT_MASS} or the string has no recognisable token at all.
	 */
	static Double massOf(String composition) {
		double mass = 0.0;
		boolean any = false;
		Matcher m = TOKEN.matcher(composition);
		while (m.find()) {
			any = true;
			Double elementMass = ELEMENT_MASS.get(m.group(1));
			if (elementMass == null) {
				return null;
			}
			int count = m.group(2) == null ? 1 : Integer.parseInt(m.group(2));
			mass += elementMass * count;
		}
		return any ? mass : null;
	}

	/**
	 * Every reference entry whose mass is within {@code toleranceDa} of
	 * {@code delta}, nearest first.
	 * <p>
	 * The tolerance should be tight, matching
	 * {@link SubstitutionMassTable#matching}'s reasoning: fitted peaks in this
	 * workflow have sigma around 2 mDa, so 5 mDa admits the peak and little else.
	 */
	public static List<Match> matching(double delta, double toleranceDa) {
		ensureLoaded();
		List<Match> out = new ArrayList<>();
		for (Map.Entry<Double, List<Entry>> e : BY_MASS
				.subMap(round(delta - toleranceDa), true, round(delta + toleranceDa), true).entrySet()) {
			for (Entry entry : e.getValue()) {
				out.add(new Match(entry, delta - entry.monoMass));
			}
		}
		out.sort((a, b) -> Double.compare(Math.abs(a.massErrorDa), Math.abs(b.massErrorDa)));
		return out;
	}

	/** Every parsed entry; diagnostics only. */
	public static List<Entry> all() {
		ensureLoaded();
		return Collections.unmodifiableList(ALL);
	}

	private static double round(double v) {
		return Math.round(v * 1e6) / 1e6;
	}
}
