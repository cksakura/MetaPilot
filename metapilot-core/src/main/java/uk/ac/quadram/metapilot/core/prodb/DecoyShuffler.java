package uk.ac.quadram.metapilot.core.prodb;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class DecoyShuffler {

	/**
	 * Generate a shuffled decoy sequence from a target sequence, such that
	 * digestion at the given cleavage residues produces peptides with the same
	 * length distribution as the original.
	 *
	 * Strategy: digest -> shuffle each peptide (preserve terminal K/R) -> join.
	 *
	 * @param seq              target protein sequence (uppercase amino acids)
	 * @param cleavageResidues residues to cleave AFTER (e.g. "KR")
	 * @param seed             random seed for reproducibility
	 * @return shuffled decoy sequence
	 */
	public static String generateShuffledDecoy(String seq, String cleavageResidues, long seed) {
		Random rnd = new Random(seed);
		List<String> peptides = digest(seq, cleavageResidues);
		StringBuilder decoy = new StringBuilder(seq.length());

		for (String pep : peptides) {
			String shuffled = shufflePeptide(pep, cleavageResidues, rnd);
			decoy.append(shuffled);
		}
		return decoy.toString();
	}

	/**
	 * Simple in-silico digestion: cleave AFTER any residue in cleavageResidues.
	 * Does NOT implement special rules like "no cleavage before P".
	 */
	private static List<String> digest(String seq, String cleavageResidues) {
		List<String> peptides = new ArrayList<>();
		StringBuilder current = new StringBuilder();

		for (int i = 0; i < seq.length(); i++) {
			char aa = seq.charAt(i);
			current.append(aa);

			if (cleavageResidues.indexOf(aa) >= 0) {
				// Cleavage after this residue
				peptides.add(current.toString());
				current.setLength(0);
			}
		}

		// Trailing fragment after last cleavage site
		if (current.length() > 0) {
			peptides.add(current.toString());
		}
		return peptides;
	}

	/**
	 * Shuffle a peptide while preserving the terminal cleavage residue (if it is in
	 * cleavageResidues). Uses Fisher–Yates shuffle on the internal residues. Tries
	 * a few times to avoid returning the original.
	 */
	private static String shufflePeptide(String peptide, String cleavageResidues, Random rnd) {

		int len = peptide.length();
		if (len <= 1) {
			return peptide; // nothing to shuffle
		}

		char last = peptide.charAt(len - 1);
		boolean hasTerminalCleavage = cleavageResidues.indexOf(last) >= 0;

		// Part we will shuffle
		int shuffleEnd = hasTerminalCleavage ? len - 1 : len;
		if (shuffleEnd <= 1) {
			// Only 1 residue to shuffle -> cannot change anything
			return peptide;
		}

		char[] original = peptide.substring(0, shuffleEnd).toCharArray();
		char[] shuffled = original.clone();

		// Try a few times to ensure shuffled != original
		final int maxTries = 10;
		for (int attempt = 0; attempt < maxTries; attempt++) {
			fisherYatesShuffle(shuffled, rnd);
			if (!arrayEquals(original, shuffled)) {
				break;
			}
		}

		StringBuilder sb = new StringBuilder(len);
		sb.append(shuffled, 0, shuffleEnd);

		if (hasTerminalCleavage) {
			sb.append(last);
		} else if (len > shuffleEnd) {
			// Append whatever tail exists (shouldn't happen in this design,
			// but kept for completeness)
			sb.append(peptide.substring(shuffleEnd));
		}

		return sb.toString();
	}

	private static void fisherYatesShuffle(char[] array, Random rnd) {
		for (int i = array.length - 1; i > 0; i--) {
			int j = rnd.nextInt(i + 1);
			char tmp = array[i];
			array[i] = array[j];
			array[j] = tmp;
		}
	}

	private static boolean arrayEquals(char[] a, char[] b) {
		if (a.length != b.length)
			return false;
		for (int i = 0; i < a.length; i++) {
			if (a[i] != b[i])
				return false;
		}
		return true;
	}

	// Simple demo
	public static void main(String[] args) {
		String seq = "MKWVTFISLLFLFSSAYSRGVFRRDTHKSEIAHRFKDLGE";
		String cleavage = "KR";
		String decoy = generateShuffledDecoy(seq, cleavage, 42L);
		System.out.println("Target: " + seq);
		System.out.println("Decoy : " + decoy);
	}
}
