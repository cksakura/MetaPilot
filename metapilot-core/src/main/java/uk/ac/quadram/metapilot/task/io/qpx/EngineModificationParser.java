package uk.ac.quadram.metapilot.task.io.qpx;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Turns each search engine's own modification notation into QPX modifications.
 * One method per notation actually reaching {@code MetaPSM} in MetaPilot:
 * <ul>
 * <li>{@link #fromMassBrackets} - Sage's {@code peptide} column and the open-search
 * canonical sequence: {@code [+42.0106]-PEPT[+79.9663]IDE-[-0.9840]};</li>
 * <li>{@link #fromFragPipeAssigned} - FragPipe's "Assigned Modifications":
 * {@code N-term(42.0106), 3C(57.0215), 5M(15.9949)};</li>
 * <li>{@link #fromPFind} - the pFind full sequence {@code PFindPSM} builds:
 * {@code LG(Gly->Pro)HDMF}, name after the residue it sits on, a leading name for
 * the N-terminus;</li>
 * <li>{@link #fromMaxQuant} - MaxQuant's "Modified sequence":
 * {@code _(Acetyl (Protein N-term))M(Oxidation (M))PEPTIDE_}.</li>
 * </ul>
 * Mass-only notations are mapped to UNIMOD by mass and site
 * ({@link UnimodCatalog#byMass}); named notations by title.
 */
public final class EngineModificationParser {

	private static final Pattern FRAGPIPE_RESIDUE = Pattern.compile("(\\d+)([A-Za-z])\\((-?[0-9.]+)\\)");
	private static final Pattern FRAGPIPE_TERMINUS = Pattern.compile("([NC])-term\\((-?[0-9.]+)\\)",
			Pattern.CASE_INSENSITIVE);

	/** Legacy MaxQuant two-letter abbreviations (pre-1.5 "Modified sequence"). */
	private static final Map<String, String> MAXQUANT_ABBREVIATIONS;
	static {
		Map<String, String> m = new HashMap<>();
		m.put("ox", "Oxidation");
		m.put("ac", "Acetyl");
		m.put("ph", "Phospho");
		m.put("de", "Deamidated");
		m.put("gl", "Gln->pyro-Glu");
		m.put("cm", "Carbamidomethyl");
		MAXQUANT_ABBREVIATIONS = Collections.unmodifiableMap(m);
	}

	private static volatile Map<String, Double> pfindMasses;

	private EngineModificationParser() {
	}

	/** Base sequence plus the modifications parsed from an engine string. */
	public static final class Parsed {
		public final String sequence;
		public final List<QpxModification> modifications;

		Parsed(String sequence, List<QpxModification> modifications) {
			this.sequence = sequence;
			this.modifications = modifications;
		}
	}

	// ------------------------------------------------------------------ mass brackets

	/**
	 * {@code [+42.0106]-PEPT[+79.9663]IDE-[-0.9840]}: a bracket after a residue sits on
	 * that residue, a bracket followed by {@code -} at the start is N-terminal, one
	 * after {@code -} at the end is C-terminal. Non-numeric bracket content is ignored
	 * (callers handle their own non-ProForma tags before calling this).
	 */
	public static Parsed fromMassBrackets(String modSeq) {
		List<QpxModification> mods = new ArrayList<>();
		if (modSeq == null) {
			return new Parsed("", mods);
		}
		StringBuilder seq = new StringBuilder();
		List<int[]> pendingPositions = new ArrayList<>();
		List<Double> pendingMasses = new ArrayList<>();
		int i = 0;
		int n = modSeq.length();
		boolean afterCTermDash = false;
		while (i < n) {
			char c = modSeq.charAt(i);
			if (c == '[') {
				int close = modSeq.indexOf(']', i);
				if (close < 0) {
					break;
				}
				double mass = parseSigned(modSeq.substring(i + 1, close));
				int next = close + 1;
				int position;
				if (seq.length() == 0) {
					position = 0; // N-terminal, with or without the trailing '-'
					if (next < n && modSeq.charAt(next) == '-') {
						next++;
					}
				} else if (afterCTermDash) {
					position = -1; // resolved to length+1 once the sequence is complete
				} else {
					position = seq.length();
				}
				if (!Double.isNaN(mass)) {
					pendingPositions.add(new int[] { position });
					pendingMasses.add(mass);
				}
				i = next;
				continue;
			}
			if (c == '-') {
				afterCTermDash = true;
				i++;
				continue;
			}
			if (Character.isLetter(c)) {
				seq.append(Character.toUpperCase(c));
			}
			i++;
		}
		String sequence = seq.toString();
		for (int k = 0; k < pendingMasses.size(); k++) {
			int p = pendingPositions.get(k)[0];
			if (p < 0) {
				p = sequence.length() + 1;
			}
			mods.add(massMod(pendingMasses.get(k), sequence, p));
		}
		return new Parsed(sequence, mods);
	}

	// ------------------------------------------------------------------ FragPipe

	/**
	 * FragPipe "Assigned Modifications" against its "Peptide" column. Signed masses
	 * are accepted (e.g. {@code 1Q(-17.0265)}).
	 */
	public static List<QpxModification> fromFragPipeAssigned(String sequence, String assigned) {
		List<QpxModification> mods = new ArrayList<>();
		if (assigned == null || assigned.trim().isEmpty() || sequence == null) {
			return mods;
		}
		for (String part : assigned.split(",")) {
			String t = part.trim();
			if (t.isEmpty()) {
				continue;
			}
			Matcher term = FRAGPIPE_TERMINUS.matcher(t);
			if (term.matches()) {
				boolean n = Character.toUpperCase(term.group(1).charAt(0)) == 'N';
				double mass = parseSigned(term.group(2));
				mods.add(massMod(mass, sequence, n ? 0 : sequence.length() + 1));
				continue;
			}
			Matcher res = FRAGPIPE_RESIDUE.matcher(t);
			if (res.matches()) {
				int pos = Integer.parseInt(res.group(1));
				mods.add(massMod(parseSigned(res.group(3)), sequence, pos));
			}
		}
		return mods;
	}

	// ------------------------------------------------------------------ pFind

	/**
	 * The pFind full sequence built by {@code PFindPSM.buildModSeq}: each
	 * {@code (name)} follows the residue it modifies; one before the first residue is
	 * N-terminal. Names may themselves contain parentheses
	 * ({@code Label:13C(6)15N(2)}), so brackets are matched by depth.
	 */
	public static Parsed fromPFind(String fullSequence) {
		List<QpxModification> mods = new ArrayList<>();
		if (fullSequence == null) {
			return new Parsed("", mods);
		}
		RawNamed raw = fromPFindLike(fullSequence);
		UnimodCatalog catalog = UnimodCatalog.get();
		for (int k = 0; k < raw.names.size(); k++) {
			String name = pFindName(raw.names.get(k));
			int pos = raw.positions.get(k);
			UnimodCatalog.Entry e = catalog.byTitle(name);
			if (e == null && name.indexOf('_') > 0) {
				// pFind spells UNIMOD's "Delta:H(2)C(2)" as "Delta_H(2)C(2)"
				e = catalog.byTitle(name.replaceFirst("_", ":"));
				if (e == null) {
					int last = name.lastIndexOf('_'); // "BEMAD_ST_2H(6)" -> "BEMAD_ST:2H(6)"
					e = catalog.byTitle(name.substring(0, last) + ":" + name.substring(last + 1));
				}
			}
			if (pos == raw.sequence.length() && e != null && !listsSite(e, raw.sequence, pos) && listsCTerm(e)) {
				pos = raw.sequence.length() + 1; // a C-terminal modification written after the last residue
			}
			double mass = e != null ? e.monoMass : pFindMass(name);
			mods.add(QpxModification.at(e, mass, name, raw.sequence, pos));
		}
		return new Parsed(raw.sequence, mods);
	}

	// ------------------------------------------------------------------ MaxQuant

	/** MaxQuant "Modified sequence", with or without the surrounding underscores. */
	public static Parsed fromMaxQuant(String modifiedSequence) {
		List<QpxModification> mods = new ArrayList<>();
		if (modifiedSequence == null) {
			return new Parsed("", mods);
		}
		String s = modifiedSequence.trim();
		if (s.startsWith("_")) {
			s = s.substring(1);
		}
		if (s.endsWith("_")) {
			s = s.substring(0, s.length() - 1);
		}
		RawNamed raw = fromPFindLike(s);
		UnimodCatalog catalog = UnimodCatalog.get();
		for (int k = 0; k < raw.names.size(); k++) {
			String title = maxQuantTitle(raw.names.get(k));
			int pos = raw.positions.get(k);
			UnimodCatalog.Entry e = catalog.byTitle(title);
			if (pos == raw.sequence.length() && e != null && !listsSite(e, raw.sequence, pos) && listsCTerm(e)) {
				pos = raw.sequence.length() + 1;
			}
			mods.add(QpxModification.at(e, Double.NaN, title, raw.sequence, pos));
		}
		return new Parsed(raw.sequence, mods);
	}

	/** "Oxidation (M)" -> "Oxidation"; "Acetyl (Protein N-term)" -> "Acetyl"; "ox" -> "Oxidation". */
	static String maxQuantTitle(String name) {
		String t = name.trim();
		String abbreviation = MAXQUANT_ABBREVIATIONS.get(t);
		if (abbreviation != null) {
			return abbreviation;
		}
		int paren = t.lastIndexOf(" (");
		if (paren > 0 && t.endsWith(")")) {
			t = t.substring(0, paren);
		}
		return t;
	}

	private static final class RawNamed {
		final String sequence;
		final List<String> names;
		final List<Integer> positions;

		RawNamed(String sequence, List<String> names, List<Integer> positions) {
			this.sequence = sequence;
			this.names = names;
			this.positions = positions;
		}
	}

	/**
	 * Repairs the name of a legacy aliased pFind token. Before PFindPSM.modificationName
	 * existed, {@code Methyl[S](Ser->Thr[S])} was written as {@code (Methyl[S](Ser->Thr)};
	 * psms.tsv files from those runs are still read on resume.
	 */
	static String pFindName(String name) {
		if (count(name, '(') > count(name, ')')) {
			int k = name.indexOf("](");
			if (k > 0) {
				String primary = name.substring(0, k + 1);
				int open = primary.lastIndexOf('[');
				return open > 0 ? primary.substring(0, open) : primary;
			}
		}
		return name;
	}

	private static int legacyTokenEnd(String s, int open) {
		int k = s.indexOf("](", open);
		if (k > 0) {
			int depth = 1;
			for (int j = k + 2; j < s.length(); j++) {
				char d = s.charAt(j);
				if (d == '(') {
					depth++;
				} else if (d == ')') {
					depth--;
					if (depth == 0) {
						return j;
					}
				}
			}
		}
		int close = s.indexOf(')', open);
		return close < 0 ? s.length() : close;
	}

	private static int count(String s, char c) {
		int n = 0;
		for (int i = 0; i < s.length(); i++) {
			if (s.charAt(i) == c) {
				n++;
			}
		}
		return n;
	}

	private static RawNamed fromPFindLike(String s) {
		StringBuilder seq = new StringBuilder();
		List<String> names = new ArrayList<>();
		List<Integer> positions = new ArrayList<>();
		int i = 0;
		int n = s.length();
		while (i < n) {
			char c = s.charAt(i);
			if (c == '(') {
				int depth = 0;
				int j = i;
				for (; j < n; j++) {
					char d = s.charAt(j);
					if (d == '(') {
						depth++;
					} else if (d == ')') {
						depth--;
						if (depth == 0) {
							break;
						}
					}
				}
				if (j >= n) {
					// Unbalanced: a legacy aliased pFind token, written as
					// "(" + NAME[SITE]( + alias-without-its-[SITE]) + ")" - e.g.
					// (Dimethyl[K](Ethyl[K]/Delta_H(4)C(2)) - so it ends where the alias
					// group, opened at "](", closes; fall back to the first ')'.
					j = legacyTokenEnd(s, i);
				}
				names.add(s.substring(i + 1, Math.min(j, n)));
				positions.add(seq.length());
				i = j + 1;
				continue;
			}
			if (Character.isLetter(c)) {
				seq.append(Character.toUpperCase(c));
			}
			i++;
		}
		return new RawNamed(seq.toString(), names, positions);
	}

	// ------------------------------------------------------------------ helpers

	/** A mass-only modification at a QPX position, mapped to UNIMOD by mass and site if unambiguous. */
	public static QpxModification massMod(double mass, String sequence, int position) {
		UnimodCatalog.Entry e = UnimodCatalog.get().byMass(mass, siteFor(sequence, position));
		return QpxModification.at(e, mass, sequence, position);
	}

	static UnimodCatalog.Site siteFor(String sequence, int position) {
		int n = sequence == null ? 0 : sequence.length();
		if (position <= 0) {
			return UnimodCatalog.Site.nTerm();
		}
		if (position > n) {
			return UnimodCatalog.Site.cTerm();
		}
		char aa = sequence.charAt(position - 1);
		if (position == 1) {
			return UnimodCatalog.Site.nTermResidue(aa);
		}
		if (position == n) {
			return UnimodCatalog.Site.cTermResidue(aa);
		}
		return UnimodCatalog.Site.residue(aa);
	}

	private static boolean listsSite(UnimodCatalog.Entry e, String sequence, int position) {
		return position >= 1 && position <= sequence.length() && e.hiddenFor(siteFor(sequence, position)) != null;
	}

	private static boolean listsCTerm(UnimodCatalog.Entry e) {
		return e.hiddenFor(UnimodCatalog.Site.cTerm()) != null;
	}

	static double parseSigned(String s) {
		if (s == null) {
			return Double.NaN;
		}
		try {
			return Double.parseDouble(s.trim());
		} catch (NumberFormatException e) {
			return Double.NaN;
		}
	}

	/**
	 * Mass of a pFind modification name that is not a UNIMOD title, from pFind's own
	 * {@code modification_default.ini} when that resource is on the classpath (the
	 * Windows GUI bundles it); NaN otherwise.
	 */
	static double pFindMass(String name) {
		Map<String, Double> m = pfindMasses;
		if (m == null) {
			synchronized (EngineModificationParser.class) {
				m = pfindMasses;
				if (m == null) {
					m = loadPFindMasses();
					pfindMasses = m;
				}
			}
		}
		Double mass = m.get(name);
		return mass == null ? Double.NaN : mass;
	}

	private static Map<String, Double> loadPFindMasses() {
		Map<String, Double> map = new HashMap<>();
		try (InputStream in = EngineModificationParser.class.getResourceAsStream("/modification_default.ini")) {
			if (in == null) {
				return map;
			}
			try (BufferedReader br = new BufferedReader(new InputStreamReader(in, StandardCharsets.UTF_8))) {
				String line;
				while ((line = br.readLine()) != null) {
					// e.g. "Oxidation[M]=M NORMAL 15.994915 15.994915 0 O(1)"
					int eq = line.indexOf('=');
					if (eq <= 0 || line.startsWith("name") || line.startsWith("@")) {
						continue;
					}
					String key = line.substring(0, eq);
					int bracket = key.lastIndexOf('[');
					String name = bracket > 0 ? key.substring(0, bracket) : key;
					String[] cs = line.substring(eq + 1).trim().split("\\s+");
					if (cs.length >= 3) {
						try {
							map.putIfAbsent(name, Double.parseDouble(cs[2]));
						} catch (NumberFormatException ignored) {
							// not a modification line
						}
					}
				}
			}
		} catch (Exception e) {
			map.clear();
		}
		return map;
	}
}
