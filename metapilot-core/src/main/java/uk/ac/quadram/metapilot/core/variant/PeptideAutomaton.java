package uk.ac.quadram.metapilot.core.variant;

import java.util.Arrays;
import java.util.List;

/**
 * Aho-Corasick automaton over a set of peptides, compiled to a complete DFA so
 * scanning costs one table lookup per residue.
 * <p>
 * I and L are one symbol, in the queries and in the scanned proteins: a mass
 * spectrometer cannot tell them apart, so a database hit on either reading
 * counts. Any byte outside the 20 residues (X, B, Z, U, O, '*') resets the
 * automaton, so no query can match across it.
 * <p>
 * Read-only after construction and safe to share between scanning threads.
 * Memory is {@code 20 x 4} bytes per trie state, i.e. bounded by the total
 * length of the query set - never by the size of what is scanned.
 */
final class PeptideAutomaton {

	/** 19 residue symbols (I folded into L) plus one reset symbol. */
	static final int ALPHABET = 20;
	private static final int RESET = ALPHABET - 1;
	private static final byte[] CODE = new byte[256];
	static {
		Arrays.fill(CODE, (byte) RESET);
		String residues = "ACDEFGHKLMNPQRSTVWY";
		for (int i = 0; i < residues.length(); i++) {
			char c = residues.charAt(i);
			CODE[c] = (byte) i;
			CODE[Character.toLowerCase(c)] = (byte) i;
		}
		CODE['I'] = CODE['L'];
		CODE['i'] = CODE['L'];
	}

	/** Receives one match; offsets are relative to the scanned range, end inclusive. */
	interface MatchSink {
		void match(int query, int start, int end);
	}

	private final int[] next;
	private final int[] output;
	private final int[] dictLink;
	private final int[] queryLength;
	private final int stateCount;

	/**
	 * @param queries distinct peptides; a query containing a non-residue letter
	 *                can never match and is kept only so indices stay aligned
	 */
	PeptideAutomaton(List<String> queries) {
		queryLength = new int[queries.size()];
		int capacity = 1;
		for (String q : queries) {
			capacity += q.length();
		}
		int[] trie = new int[capacity * ALPHABET];
		Arrays.fill(trie, -1);
		int[] out = new int[capacity];
		Arrays.fill(out, -1);
		int states = 1;
		for (int qi = 0; qi < queries.size(); qi++) {
			String q = queries.get(qi);
			queryLength[qi] = q.length();
			if (q.isEmpty() || !allResidues(q)) {
				continue;
			}
			int s = 0;
			for (int i = 0; i < q.length(); i++) {
				int sym = CODE[q.charAt(i) & 0xff];
				int t = trie[s * ALPHABET + sym];
				if (t < 0) {
					t = states++;
					trie[s * ALPHABET + sym] = t;
				}
				s = t;
			}
			if (out[s] < 0) {
				out[s] = qi;
			}
		}
		stateCount = states;
		next = Arrays.copyOf(trie, states * ALPHABET);
		output = Arrays.copyOf(out, states);
		dictLink = new int[states];
		Arrays.fill(dictLink, -1);
		int[] fail = new int[states];
		int[] queue = new int[states];
		int head = 0;
		int tail = 0;
		for (int sym = 0; sym < ALPHABET; sym++) {
			int child = next[sym];
			if (sym == RESET || child < 0) {
				next[sym] = 0;
			} else {
				fail[child] = 0;
				queue[tail++] = child;
			}
		}
		while (head < tail) {
			int u = queue[head++];
			for (int sym = 0; sym < ALPHABET; sym++) {
				int v = next[u * ALPHABET + sym];
				if (sym == RESET) {
					next[u * ALPHABET + sym] = 0;
				} else if (v < 0) {
					next[u * ALPHABET + sym] = next[fail[u] * ALPHABET + sym];
				} else {
					int f = next[fail[u] * ALPHABET + sym];
					fail[v] = f;
					dictLink[v] = output[f] >= 0 ? f : dictLink[f];
					queue[tail++] = v;
				}
			}
		}
	}

	int stateCount() {
		return stateCount;
	}

	int queryLength(int query) {
		return queryLength[query];
	}

	/** Reports every query occurrence in {@code seq[from, to)}. */
	void scan(byte[] seq, int from, int to, MatchSink sink) {
		int s = 0;
		for (int i = from; i < to; i++) {
			s = next[s * ALPHABET + CODE[seq[i] & 0xff]];
			int t = output[s] >= 0 ? s : dictLink[s];
			while (t >= 0) {
				int q = output[t];
				sink.match(q, i - from - queryLength[q] + 1, i - from);
				t = dictLink[t];
			}
		}
	}

	/** I and L collapsed to L; everything upper-cased. */
	static String normalise(String peptide) {
		return peptide.toUpperCase(java.util.Locale.ROOT).replace('I', 'L');
	}

	private static boolean allResidues(String q) {
		for (int i = 0; i < q.length(); i++) {
			char c = q.charAt(i);
			if (c > 255 || CODE[c] == RESET) {
				return false;
			}
		}
		return true;
	}
}
