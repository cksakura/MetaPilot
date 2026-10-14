package uk.ac.quadram.metapilot.core.variant;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.function.Consumer;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import uk.ac.quadram.metapilot.core.opensearch.IsobaricModificationTable;
import uk.ac.quadram.metapilot.dbSearch.pfind.PFindPSM;
import uk.ac.quadram.metapilot.dbSearch.pfind.PFindPsmReader;

/**
 * Front end for pFind open search: substitution candidates from pFind's own
 * calls.
 * <p>
 * pFind writes an amino-acid substitution as an ordinary entry of the
 * {@code Modification} column, {@code pos,Xxx->Yyy[X]}, with {@code pos} 1-based
 * on the residue (0 is the peptide N-terminus) and {@code Sequence} still the
 * database peptide. Those calls are taken as pFind made them - same PSMs, same
 * FDR - and only translated into {@link VariantCandidate}s:
 * <ul>
 * <li>An entry like {@code Methyl[S](Ser->Thr[S])} is pFind saying the two
 * readings have the same mass on the same residue. It becomes an
 * undecidable candidate, with pFind's alias listed as the isobaric
 * modification.</li>
 * <li>The residue the token names must be the one at that position (I and L
 * interchangeable, {@code Xle} = either). A mismatch means the column was not
 * read the way pFind meant it; the PSM is skipped and counted, not
 * guessed.</li>
 * <li>A substitution to or from K or R changes the tryptic digest, so the
 * variant peptide would not be the peptide the protein yields; such PSMs are
 * skipped. I&lt;-&gt;L sites carry no information and are dropped.</li>
 * <li>Several substitutions in one PSM give one candidate per site, all with
 * the variant peptide that has every substitution applied.</li>
 * </ul>
 * Mass deltas use an unmodified cysteine: in open-search mode MetaPilot writes
 * no fixed modifications into the pFind configuration
 * ({@code PFindConfig.applyModifications}). MS2 decidability comes from the
 * shared {@link IsobaricModificationTable} with {@code undecidable_tolerance_da}.
 */
public final class PFindVariantCandidateSource {

	/** Candidates and the control pool of one pFind run, with what was skipped. */
	public static final class Result {
		public final List<VariantCandidate> candidates = new ArrayList<>();
		/** Target base sequences of every accepted PSM; the matched-control pool. */
		public final Set<String> basePeptides = new TreeSet<>();
		public int psms;
		public int substitutionPsms;
		public int skippedKR;
		public int skippedResidueMismatch;
		public int skippedUnknownResidue;
		public int skippedNonMicrobial;
		final Map<String, Integer> mismatchExamples = new TreeMap<>();

		public String describe() {
			return String.format(Locale.ROOT,
					"%,d target PSMs, %,d with substitutions -> %,d candidates; skipped: %,d K/R, %,d residue "
							+ "mismatch%s, %,d unknown residue, %,d without a microbial protein",
					psms, substitutionPsms, candidates.size(), skippedKR, skippedResidueMismatch,
					mismatchExamples.isEmpty() ? "" : " " + mismatchExamples, skippedUnknownResidue,
					skippedNonMicrobial);
		}
	}

	private static final Pattern PLAIN = Pattern.compile("^([A-Z][a-z]{2})->([A-Z][a-z]{2})\\[([^\\]]+)\\]$");
	/** {@code Methyl[S](Ser->Thr[S])}: a modification aliased to a substitution. */
	private static final Pattern ALIAS_TAIL = Pattern
			.compile("^(.+)\\(([A-Z][a-z]{2})->([A-Z][a-z]{2})\\[([^\\]]+)\\]\\)$");
	/** {@code Ser->Thr[S](Methyl[S])}: the same, written the other way round. */
	private static final Pattern ALIAS_HEAD = Pattern
			.compile("^([A-Z][a-z]{2})->([A-Z][a-z]{2})\\[([^\\]]+)\\]\\((.+)\\)$");

	private static final Map<String, Character> THREE_LETTER = new HashMap<>();
	static {
		String[][] codes = { { "Ala", "A" }, { "Arg", "R" }, { "Asn", "N" }, { "Asp", "D" }, { "Cys", "C" },
				{ "Gln", "Q" }, { "Glu", "E" }, { "Gly", "G" }, { "His", "H" }, { "Ile", "I" }, { "Leu", "L" },
				{ "Lys", "K" }, { "Met", "M" }, { "Phe", "F" }, { "Pro", "P" }, { "Ser", "S" }, { "Thr", "T" },
				{ "Trp", "W" }, { "Tyr", "Y" }, { "Val", "V" }, { "Xle", "L" } };
		for (String[] c : codes) {
			THREE_LETTER.put(c[0], c[1].charAt(0));
		}
	}

	/** Cysteine as pFind's open search sees it. */
	public static final double FIXED_CYSTEINE_MOD_DA = 0.0;

	private PFindVariantCandidateSource() {
	}

	/**
	 * @param spectraFiles          each sample's {@code sample_db/pFind-Filtered.spectra}
	 * @param microbialProteinPrefix when non-empty, PSMs none of whose proteins
	 *                              start with it (host-only PSMs) are left out of
	 *                              both candidates and control pool: UHGP-100
	 *                              holds gut-microbial proteins only
	 */
	public static Result read(List<File> spectraFiles, String microbialProteinPrefix, double undecidableToleranceDa,
			Consumer<String> log) {
		Result result = new Result();
		for (File f : spectraFiles) {
			if (f == null || !f.isFile() || f.length() == 0) {
				if (log != null) {
					log.accept("Pangenome variant check: pFind result not found, skipped: " + f);
				}
				continue;
			}
			for (PFindPSM psm : new PFindPsmReader(f).getPsms()) {
				String seq = psm.getBaseSequence();
				if (!psm.isTarget() || seq == null || seq.contains("X")) {
					continue;
				}
				if (microbialProteinPrefix != null && !microbialProteinPrefix.isEmpty()
						&& !hasPrefix(psm.getProteins(), microbialProteinPrefix)) {
					result.skippedNonMicrobial++;
					continue;
				}
				result.psms++;
				result.basePeptides.add(seq);
				List<VariantCandidate> fromPsm = fromPsm(psm.getFileName(), psm.getScan(), psm.getCharge(), seq,
						psm.getMod(), undecidableToleranceDa, result);
				if (!fromPsm.isEmpty()) {
					double mz = psm.getCharge() > 0 && psm.getExpMH() > 0
							? (psm.getExpMH() - OpenSearchVariantCandidateSource.PROTON) / psm.getCharge()
									+ OpenSearchVariantCandidateSource.PROTON
							: Double.NaN;
					for (VariantCandidate c : fromPsm) {
						c.qValue = psm.getqValue();
						c.precursorMz = mz;
					}
					String[] proteins = psm.getProteins();
					if (proteins != null) {
						for (VariantCandidate c : fromPsm) {
							for (String protein : proteins) {
								if (protein != null && !protein.isEmpty()) {
									c.assignedProteins.add(protein);
								}
							}
						}
					}
					result.substitutionPsms++;
					result.candidates.addAll(fromPsm);
				}
			}
		}
		return result;
	}

	private static boolean hasPrefix(String[] proteins, String prefix) {
		if (proteins == null) {
			return false;
		}
		for (String p : proteins) {
			if (p != null && p.startsWith(prefix)) {
				return true;
			}
		}
		return false;
	}

	/** One parsed substitution entry. */
	static final class Token {
		final int position;
		final char from;
		final char to;
		/** pFind's isobaric alias, or null for a plain substitution. */
		final String alias;
		final String raw;

		Token(int position, char from, char to, String alias, String raw) {
			this.position = position;
			this.from = from;
			this.to = to;
			this.alias = alias;
			this.raw = raw;
		}
	}

	/**
	 * Substitution entries of a {@code Modification} value, in column order.
	 * Entries that are not substitutions are ignored; an entry naming a residue
	 * outside the 20 (plus {@code Xle}) is returned with {@code from} or
	 * {@code to} = 0.
	 */
	static List<Token> parse(String modification) {
		List<Token> out = new ArrayList<>();
		if (modification == null || modification.isBlank()) {
			return out;
		}
		for (String entry : modification.split(";")) {
			String e = entry.trim();
			int comma = e.indexOf(',');
			if (comma <= 0) {
				continue;
			}
			int pos;
			try {
				pos = Integer.parseInt(e.substring(0, comma).trim());
			} catch (NumberFormatException ex) {
				continue;
			}
			String name = e.substring(comma + 1).trim();
			if (!name.contains("->")) {
				continue;
			}
			Matcher m = PLAIN.matcher(name);
			if (m.matches()) {
				out.add(new Token(pos, residue(m.group(1)), residue(m.group(2)), null, e));
				continue;
			}
			m = ALIAS_TAIL.matcher(name);
			if (m.matches()) {
				out.add(new Token(pos, residue(m.group(2)), residue(m.group(3)), m.group(1), e));
				continue;
			}
			m = ALIAS_HEAD.matcher(name);
			if (m.matches()) {
				out.add(new Token(pos, residue(m.group(1)), residue(m.group(2)), m.group(4), e));
			}
		}
		return out;
	}

	private static char residue(String threeLetter) {
		Character c = THREE_LETTER.get(threeLetter);
		return c == null ? 0 : c;
	}

	/** Candidates for one PSM; empty when it has no usable substitution. */
	static List<VariantCandidate> fromPsm(String file, int scan, int charge, String sequence, String modification,
			double undecidableToleranceDa, Result stats) {
		List<Token> tokens = parse(modification);
		List<VariantCandidate> out = new ArrayList<>();
		if (tokens.isEmpty()) {
			return out;
		}
		Map<Integer, Token> byPosition = new TreeMap<>();
		for (Token t : tokens) {
			if (t.from == 0 || t.to == 0) {
				stats.skippedUnknownResidue++;
				return out;
			}
			if (isKR(t.from) || isKR(t.to)) {
				stats.skippedKR++;
				return out;
			}
			if (t.position < 1 || t.position > sequence.length()
					|| !sameResidue(sequence.charAt(t.position - 1), t.from)) {
				stats.skippedResidueMismatch++;
				if (stats.mismatchExamples.size() < 5) {
					stats.mismatchExamples.merge(t.raw + " on " + sequence, 1, Integer::sum);
				}
				return out;
			}
			if (isIL(t.from) && isIL(t.to)) {
				continue;
			}
			byPosition.putIfAbsent(t.position, t);
		}
		if (byPosition.isEmpty()) {
			return out;
		}
		StringBuilder variant = new StringBuilder(sequence);
		for (Token t : byPosition.values()) {
			variant.setCharAt(t.position - 1, t.to);
		}
		String variantPeptide = variant.toString();
		for (Token t : byPosition.values()) {
			VariantCandidate c = new VariantCandidate();
			c.engine = VariantCandidate.Engine.PFIND;
			c.file = file;
			c.scan = scan;
			c.psmId = charge;
			c.charge = charge;
			c.basePeptide = sequence;
			c.position = t.position;
			c.fromAa = sequence.charAt(t.position - 1);
			c.toAa = t.to;
			c.variantPeptide = variantPeptide;
			c.deltaMass = IsobaricModificationTable.substitutionMass(c.fromAa, c.toAa, FIXED_CYSTEINE_MOD_DA);
			c.isobaricModifications.addAll(IsobaricModificationTable.sameResidue(c.fromAa, c.deltaMass,
					t.position == 1, undecidableToleranceDa));
			if (t.alias != null) {
				c.isobaricModifications.add("pFind alias " + t.alias);
			}
			c.ms2Decidable = c.isobaricModifications.isEmpty();
			c.frontEndEvidence = t.raw;
			out.add(c);
		}
		return out;
	}

	/** The searched databases: each sample's DBReducer output, or its marker database when DBReducer was off. */
	public static List<File> searchedDatabases(File[] sampleFolders, String sampleDbDir, String markerDir) {
		Set<File> out = new LinkedHashSet<>();
		if (sampleFolders == null) {
			return new ArrayList<>(out);
		}
		for (File sample : sampleFolders) {
			File reduced = new File(new File(new File(sample, sampleDbDir), "DBReducer"), "DBReducer.fasta");
			File marker = new File(new File(sample, markerDir), sample.getName() + ".fasta");
			if (reduced.isFile() && reduced.length() > 0) {
				out.add(reduced);
			} else if (marker.isFile() && marker.length() > 0) {
				out.add(marker);
			}
		}
		return new ArrayList<>(out);
	}

	private static boolean sameResidue(char a, char b) {
		char x = Character.toUpperCase(a);
		char y = Character.toUpperCase(b);
		return x == y || (isIL(x) && isIL(y));
	}

	private static boolean isKR(char c) {
		return c == 'K' || c == 'R';
	}

	private static boolean isIL(char c) {
		return c == 'I' || c == 'L';
	}
}
