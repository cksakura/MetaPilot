package uk.ac.quadram.metapilot.core.dia;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.List;
import java.util.SplittableRandom;

import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNPrecursor;

/**
 * Paired shuffled entrapment proteins for the entrapment Deep mode ({@code dia-marker-deep-entrapment}).
 * <p>
 * Each protein is cut after every K/R into segments; in each segment the first and last residues stay in place and
 * the standard residues in between are permuted. So the shuffled protein digests exactly like its source with
 * Trypsin and Trypsin/P (K/R sites, the residue after each site and the protein N-terminus are unchanged), and every
 * peptide maps 1:1, by position, to a source peptide of the same length, composition and mass. A shuffled segment
 * equal to the source (I = L) is reshuffled (up to {@link #MAX_TRIES}). The same rule as
 * {@code scripts/05_reliability/build_shadow_catalogue.py}; one RNG per protein and segment, seeded from the seed and
 * the IDs, so the result is reproducible.
 * <p>
 * IDs: MGYG000000001_00005 -> MGYG800000001_00005 (digit {@link #DEFAULT_ID_DIGIT}); versioned genomes keep their
 * version (MGYG000000761.1_00014 -> MGYG800000761.1_00014). Peptides of shuffled proteins that also occur as target
 * peptides are classified at analysis time.
 */
public final class EntrapmentShuffle {

	public static final int DEFAULT_ID_DIGIT = 8;
	public static final int MAX_TRIES = 20;
	private static final String STANDARD = "ACDEFGHILMNPQSTVWY";

	private EntrapmentShuffle() {
	}

	/** Entrapment ID of a genome or protein ID ({@code MGYG0...} -> {@code MGYG<digit>...}). */
	public static String entrapmentId(String id, int digit) {
		int us = id.indexOf('_');
		String genome = us < 0 ? id : id.substring(0, us);
		String rest = us < 0 ? "" : id.substring(us);
		String base = genome;
		String version = "";
		int dot = genome.indexOf('.');
		if (dot >= 0) {
			base = genome.substring(0, dot);
			version = genome.substring(dot);
		}
		if (!(base.startsWith("MGYG") && base.length() == 13 && base.charAt(4) == '0'
				&& base.substring(4).chars().allMatch(Character::isDigit))) {
			throw new IllegalArgumentException("unexpected genome ID " + genome + " (need MGYG0xxxxxxxx[.N])");
		}
		return "MGYG" + digit + base.substring(5) + version + rest;
	}

	public static boolean isEntrapment(String proteinOrGenomeId, int digit) {
		return proteinOrGenomeId != null && proteinOrGenomeId.startsWith("MGYG" + digit);
	}

	/** Target precursor: at least one non-decoy protein that is not an entrapment protein. */
	public static boolean hasTargetProtein(DiaNNPrecursor p, int digit) {
		String[] pros = p == null ? null : p.getProGroups();
		if (pros == null) return false;
		for (String pro : pros) {
			if (pro == null || pro.isEmpty() || isDecoyId(pro)) continue;
			if (!isEntrapment(pro, digit)) return true;
		}
		return false;
	}

	/** Entrapment precursor: all its non-decoy proteins are entrapment proteins. */
	public static boolean onlyEntrapmentProteins(DiaNNPrecursor p, int digit) {
		String[] pros = p == null ? null : p.getProGroups();
		if (pros == null) return false;
		boolean any = false;
		for (String pro : pros) {
			if (pro == null || pro.isEmpty() || isDecoyId(pro)) continue;
			if (!isEntrapment(pro, digit)) return false;
			any = true;
		}
		return any;
	}

	public static boolean isDecoyId(String id) {
		return id.startsWith("REV_") || id.startsWith("rev_") || id.startsWith("DECOY_") || id.startsWith("decoy_");
	}

	/** Cleavage segments: [start, end) after every K/R (no proline rule). */
	static List<int[]> segments(String seq) {
		List<int[]> out = new ArrayList<>();
		int start = 0;
		for (int i = 0; i < seq.length(); i++) {
			char a = seq.charAt(i);
			if ((a == 'K' || a == 'R') && i + 1 < seq.length()) {
				out.add(new int[] { start, i + 1 });
				start = i + 1;
			}
		}
		if (start < seq.length()) out.add(new int[] { start, seq.length() });
		return out;
	}

	/** Shuffled copy of {@code sequence} (see class comment). */
	public static String shuffle(String sequence, String entrapmentProteinId, long seed) {
		String seq = sequence.toUpperCase();
		char[] out = seq.toCharArray();
		int k = 0;
		for (int[] seg : segments(seq)) {
			int s = seg[0], e = seg[1];
			List<Integer> movable = new ArrayList<>();
			for (int i = s + 1; i < e - 1; i++) {
				if (STANDARD.indexOf(seq.charAt(i)) >= 0) movable.add(i);
			}
			boolean distinct = false;
			for (int i = 1; i < movable.size() && !distinct; i++) {
				distinct = seq.charAt(movable.get(i)) != seq.charAt(movable.get(0));
			}
			if (movable.size() >= 2 && distinct) {
				String srcKey = seq.substring(s, e).replace('I', 'L');
				SplittableRandom rng = new SplittableRandom(hash64(seed + "|" + entrapmentProteinId + "|" + k));
				char[] cand = null;
				for (int attempt = 0; attempt < MAX_TRIES; attempt++) {
					char[] res = new char[movable.size()];
					for (int i = 0; i < res.length; i++) res[i] = seq.charAt(movable.get(i));
					for (int i = res.length - 1; i > 0; i--) {
						int j = rng.nextInt(i + 1);
						char t = res[i];
						res[i] = res[j];
						res[j] = t;
					}
					cand = seq.substring(s, e).toCharArray();
					for (int i = 0; i < res.length; i++) cand[movable.get(i) - s] = res[i];
					if (!new String(cand).replace('I', 'L').equals(srcKey)) break;
				}
				System.arraycopy(cand, 0, out, s, cand.length);
			}
			k++;
		}
		return new String(out);
	}

	/**
	 * Appends the shuffled copy of every protein of {@code genomeFaa} to {@code w}, with entrapment IDs (header
	 * descriptions kept). Returns the number of proteins written.
	 */
	public static int appendShuffledGenome(File genomeFaa, PrintWriter w, int digit, long seed) throws IOException {
		int n = 0;
		try (BufferedReader r = new BufferedReader(new FileReader(genomeFaa, StandardCharsets.UTF_8))) {
			String header = null;
			StringBuilder seq = new StringBuilder();
			String line;
			while ((line = r.readLine()) != null) {
				if (line.startsWith(">")) {
					if (header != null) n += writeShuffled(header, seq.toString(), w, digit, seed);
					header = line.substring(1).trim();
					seq.setLength(0);
				} else {
					seq.append(line.trim());
				}
			}
			if (header != null) n += writeShuffled(header, seq.toString(), w, digit, seed);
		}
		return n;
	}

	private static int writeShuffled(String header, String seq, PrintWriter w, int digit, long seed) {
		if (seq.isEmpty()) return 0;
		int sp = header.indexOf(' ');
		String id = sp < 0 ? header : header.substring(0, sp);
		String desc = sp < 0 ? "" : header.substring(sp);
		String eid = entrapmentId(id, digit);
		w.println(">" + eid + desc);
		w.println(shuffle(seq, eid, seed));
		return 1;
	}

	/** 64-bit FNV-1a of the UTF-8 bytes (stable across JVMs). */
	static long hash64(String s) {
		long h = 0xcbf29ce484222325L;
		for (byte b : s.getBytes(StandardCharsets.UTF_8)) {
			h ^= (b & 0xff);
			h *= 0x100000001b3L;
		}
		return h;
	}
}
