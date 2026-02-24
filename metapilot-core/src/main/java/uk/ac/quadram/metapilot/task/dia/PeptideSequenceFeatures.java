package uk.ac.quadram.metapilot.task.dia;

import uk.ac.quadram.metapilot.core.aminoacid.Aminoacids;

/**
 * Peptide sequence-derived features for the Deep Comprehensive training dataset:
 * GRAVY (hydrophobicity), basic/acidic/aromatic counts, optional net charge.
 * All values are computable from the stripped sequence only.
 */
public final class PeptideSequenceFeatures {

	private PeptideSequenceFeatures() {}

	private static final double[] HYDRO = Aminoacids.HydroScore;
	private static final int HYDRO_OFFSET = 'A';

	/** Basic (K, R, H), acidic (D, E), aromatic (F, Y, W). */
	public static final String BASIC = "KRH";
	public static final String ACIDIC = "DE";
	public static final String AROMATIC = "FYW";

	/**
	 * GRAVY = mean Kyte–Doolittle hydropathicity over residues.
	 * Unknown residues are treated as 0.
	 */
	public static double gravy(String sequence) {
		if (sequence == null || sequence.isEmpty()) {
			return 0.0;
		}
		String s = sequence.toUpperCase();
		double sum = 0.0;
		int n = 0;
		for (int i = 0; i < s.length(); i++) {
			char c = s.charAt(i);
			int idx = c - HYDRO_OFFSET;
			if (idx >= 0 && idx < HYDRO.length) {
				sum += HYDRO[idx];
				n++;
			}
		}
		return n > 0 ? sum / n : 0.0;
	}

	/** Count of basic residues (K, R, H). */
	public static int basicCount(String sequence) {
		return countChars(sequence, BASIC);
	}

	/** Count of acidic residues (D, E). */
	public static int acidicCount(String sequence) {
		return countChars(sequence, ACIDIC);
	}

	/** Count of aromatic residues (F, Y, W). */
	public static int aromaticCount(String sequence) {
		return countChars(sequence, AROMATIC);
	}

	/** Net charge proxy at neutral pH: (K+R+H) - (D+E). */
	public static int netCharge(String sequence) {
		return basicCount(sequence) - acidicCount(sequence);
	}

	/**
	 * Tryptic missed cleavages from sequence: count of K/R not at the C-terminus.
	 * Used when precursor miss-cleavage is not available (e.g. prediction from predicted TSV).
	 */
	public static int missCleavageFromSequence(String sequence) {
		if (sequence == null || sequence.length() <= 1) return 0;
		String s = sequence.toUpperCase();
		int count = 0;
		for (int i = 0; i < s.length() - 1; i++) {
			char c = s.charAt(i);
			if (c == 'K' || c == 'R') count++;
		}
		return count;
	}

	private static int countChars(String sequence, String set) {
		if (sequence == null || set == null) {
			return 0;
		}
		String s = sequence.toUpperCase();
		int count = 0;
		for (int i = 0; i < s.length(); i++) {
			if (set.indexOf(s.charAt(i)) >= 0) {
				count++;
			}
		}
		return count;
	}
}
