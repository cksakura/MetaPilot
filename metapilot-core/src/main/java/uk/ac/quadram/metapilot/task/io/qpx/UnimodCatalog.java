package uk.ac.quadram.metapilot.task.io.qpx;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * UNIMOD lookup for the QPX {@code modifications} array and ProForma accessions.
 * <p>
 * Backed by {@code /unimod_compact.tsv} on the classpath, a flat extract of the
 * {@code unimod.xml} that already ships in this repository under
 * {@code dbSearch/pia/} (that XML sits in the Java source tree, which Maven does not
 * package, so it is not reachable at run time). Columns: record_id, title,
 * mono_mass, and {@code site@position} specificities with a trailing {@code !}
 * marking UNIMOD's "hidden" (rarely used) ones.
 * <p>
 * Two lookups, matching what the engines report:
 * <ul>
 * <li>{@link #byTitle} for engines that name modifications (pFind, MaxQuant);</li>
 * <li>{@link #byMass} for engines that only give a mass (FragPipe, Sage, the
 * open-search fixed modifications). A mass maps to a UNIMOD entry only when exactly
 * one entry both matches the mass and lists the site; ties are broken in favour of
 * non-hidden specificities, and anything still ambiguous stays mass-only rather
 * than being given a guessed name.</li>
 * </ul>
 */
public final class UnimodCatalog {

	private static final Logger LOGGER = LogManager.getLogger(UnimodCatalog.class);

	/** Default mass tolerance: engines print 4 decimals (57.0215 vs UNIMOD 57.021464). */
	public static final double DEFAULT_MASS_TOLERANCE = 0.001;

	public static final String RESOURCE = "/unimod_compact.tsv";

	/** A UNIMOD record. */
	public static final class Entry {
		public final int recordId;
		public final String title;
		public final double monoMass;
		/** site@position -> hidden */
		final Map<String, Boolean> specificities;

		Entry(int recordId, String title, double monoMass, Map<String, Boolean> specificities) {
			this.recordId = recordId;
			this.title = title;
			this.monoMass = monoMass;
			this.specificities = specificities;
		}

		public String accession() {
			return "UNIMOD:" + recordId;
		}

		/**
		 * @return null if this entry does not list the site, otherwise whether that
		 *         specificity is hidden
		 */
		Boolean hiddenFor(Site site) {
			Boolean hidden = null;
			for (String key : site.keys()) {
				Boolean h = specificities.get(key);
				if (h != null && (hidden == null || !h)) {
					hidden = h;
				}
			}
			return hidden;
		}
	}

	/** Where a modification sits, in UNIMOD's site/position vocabulary. */
	public static final class Site {
		final char residue; // 0 for a bare terminus
		final boolean nTerm;
		final boolean cTerm;

		private Site(char residue, boolean nTerm, boolean cTerm) {
			this.residue = residue;
			this.nTerm = nTerm;
			this.cTerm = cTerm;
		}

		public static Site residue(char aa) {
			return new Site(Character.toUpperCase(aa), false, false);
		}

		/** A residue that is also the peptide N-terminal one (e.g. pyro-Glu on Q). */
		public static Site nTermResidue(char aa) {
			return new Site(Character.toUpperCase(aa), true, false);
		}

		public static Site cTermResidue(char aa) {
			return new Site(Character.toUpperCase(aa), false, true);
		}

		public static Site nTerm() {
			return new Site((char) 0, true, false);
		}

		public static Site cTerm() {
			return new Site((char) 0, false, true);
		}

		List<String> keys() {
			List<String> keys = new ArrayList<>(6);
			if (residue != 0) {
				keys.add(residue + "@Anywhere");
				if (nTerm) {
					keys.add(residue + "@Any N-term");
					keys.add(residue + "@Protein N-term");
				}
				if (cTerm) {
					keys.add(residue + "@Any C-term");
					keys.add(residue + "@Protein C-term");
				}
			} else if (nTerm) {
				keys.add("N-term@Any N-term");
				keys.add("N-term@Protein N-term");
			} else if (cTerm) {
				keys.add("C-term@Any C-term");
				keys.add("C-term@Protein C-term");
			}
			return keys;
		}
	}

	private static volatile UnimodCatalog shared;

	private final List<Entry> entries;
	private final Map<String, Entry> byExactTitle = new HashMap<>();
	private final Map<String, Entry> byLowerTitle = new HashMap<>();

	private UnimodCatalog(List<Entry> entries) {
		this.entries = Collections.unmodifiableList(entries);
		for (Entry e : entries) {
			byExactTitle.putIfAbsent(e.title, e);
			byLowerTitle.putIfAbsent(e.title.toLowerCase(Locale.ROOT), e);
		}
	}

	/**
	 * The catalog loaded from the classpath. If the resource is missing the catalog is
	 * empty (every modification is then written mass-only) and a warning is logged
	 * once - the PSM export must never fail because of it.
	 */
	public static UnimodCatalog get() {
		UnimodCatalog c = shared;
		if (c == null) {
			synchronized (UnimodCatalog.class) {
				c = shared;
				if (c == null) {
					c = load();
					shared = c;
				}
			}
		}
		return c;
	}

	private static UnimodCatalog load() {
		List<Entry> list = new ArrayList<>();
		try (InputStream in = UnimodCatalog.class.getResourceAsStream(RESOURCE)) {
			if (in == null) {
				LOGGER.warn("{} not found on the classpath; QPX modifications will carry masses only", RESOURCE);
				return new UnimodCatalog(list);
			}
			try (BufferedReader br = new BufferedReader(new InputStreamReader(in, StandardCharsets.UTF_8))) {
				String line;
				boolean header = true;
				while ((line = br.readLine()) != null) {
					if (line.startsWith("#") || line.isEmpty()) {
						continue;
					}
					if (header) {
						header = false;
						continue;
					}
					String[] cs = line.split("\t", -1);
					if (cs.length < 4) {
						continue;
					}
					Map<String, Boolean> spec = new HashMap<>();
					for (String s : cs[3].split(";")) {
						if (s.isEmpty()) {
							continue;
						}
						boolean hidden = s.endsWith("!");
						spec.put(hidden ? s.substring(0, s.length() - 1) : s, hidden);
					}
					list.add(new Entry(Integer.parseInt(cs[0]), cs[1], Double.parseDouble(cs[2]), spec));
				}
			}
		} catch (IOException | RuntimeException e) {
			LOGGER.warn("Could not read {}; QPX modifications will carry masses only: {}", RESOURCE, e.getMessage());
			list.clear();
		}
		return new UnimodCatalog(list);
	}

	public int size() {
		return entries.size();
	}

	/** Exact UNIMOD title, then case-insensitive; null if neither matches. */
	public Entry byTitle(String title) {
		if (title == null) {
			return null;
		}
		String t = title.trim();
		Entry e = byExactTitle.get(t);
		return e != null ? e : byLowerTitle.get(t.toLowerCase(Locale.ROOT));
	}

	/**
	 * The single UNIMOD entry whose mass is within {@code tolerance} of {@code mass}
	 * and which lists {@code site}; null when there is none or the match is ambiguous.
	 */
	public Entry byMass(double mass, Site site, double tolerance) {
		if (Double.isNaN(mass) || site == null) {
			return null;
		}
		Entry visible = null;
		int visibleCount = 0;
		Entry hidden = null;
		int hiddenCount = 0;
		for (Entry e : entries) {
			if (Math.abs(e.monoMass - mass) > tolerance) {
				continue;
			}
			Boolean h = e.hiddenFor(site);
			if (h == null) {
				continue;
			}
			if (h) {
				hidden = e;
				hiddenCount++;
			} else {
				visible = e;
				visibleCount++;
			}
		}
		if (visibleCount == 1) {
			return visible;
		}
		if (visibleCount == 0 && hiddenCount == 1) {
			return hidden;
		}
		return null;
	}

	public Entry byMass(double mass, Site site) {
		return byMass(mass, site, DEFAULT_MASS_TOLERANCE);
	}
}
