package uk.ac.quadram.metapilot.core.opensearch;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * One search candidate - a (spectrum, peptide) pair - carrying whichever
 * engines proposed it.
 * <p>
 * This is the unit the combined classifier scores, replacing the per-engine
 * row. Two engines proposing the same peptide for the same spectrum produce one
 * candidate with both blocks populated, not two rows competing with each other.
 *
 * <h2>Why missing blocks are not a leak</h2>
 * Most candidates are proposed by one engine only, so most instances have one
 * block absent, and absence marks which engine found it. That was the stated
 * reason for scoring per engine instead. Measured on the 2026-09-10 combined
 * run the objection does not hold:
 * <pre>
 *   MSFragger only   3,768,634 candidates   42.36% decoy
 *   Sage only        4,072,018 candidates   47.29% decoy
 *   both engines     2,162,486 candidates   15.14% decoy
 * </pre>
 * Agreement is the strongest single discriminant available - a candidate both
 * engines propose is roughly three times cleaner - while the two single-engine
 * groups sit five points apart, so engine identity separates almost nothing. A
 * model given presence indicators learns agreement, which is evidence, and
 * cannot meaningfully learn "which engine", which is not.
 * <p>
 * Absence is represented by an explicit indicator plus genuine missing values,
 * never by imputation: a mean-filled block would assert a score the engine
 * never reported.
 */
public final class CombinedPsmCandidate {

	/**
	 * Key separator. A tab cannot occur in a file name, a scan number or a
	 * peptide sequence read from pepXML or a Sage TSV - the TSV reader would
	 * already have split on it - so keys cannot collide by concatenation.
	 */
	private static final char KEY_SEPARATOR = '\t';

	private final String fileName;
	private final int scan;
	private final String peptide;
	private final boolean target;

	private OpenSearchPsm msfragger;
	private OpenSearchPsm sage;

	CombinedPsmCandidate(String fileName, int scan, String peptide, boolean target) {
		this.fileName = fileName;
		this.scan = scan;
		this.peptide = peptide;
		this.target = target;
	}

	/**
	 * The row carrying shared, engine-independent properties: charge, precursor
	 * mass, peptide geometry and the delta-mass peak context. MSFragger is
	 * preferred purely for determinism when both are present; what these getters
	 * read are properties of the candidate, not of the search that found it.
	 */
	public OpenSearchPsm primary() {
		return msfragger != null ? msfragger : sage;
	}

	public OpenSearchPsm msfragger() {
		return msfragger;
	}

	public OpenSearchPsm sage() {
		return sage;
	}

	public boolean foundByMsfragger() {
		return msfragger != null;
	}

	public boolean foundBySage() {
		return sage != null;
	}

	public boolean foundByBoth() {
		return msfragger != null && sage != null;
	}

	public String getFileName() {
		return fileName;
	}

	public int getScan() {
		return scan;
	}

	public String getPeptide() {
		return peptide;
	}

	public boolean isTarget() {
		return target;
	}

	/** Spectrum key, used for fold assignment and the spectrum-level collapse. */
	public String spectrumKey() {
		return fileName + KEY_SEPARATOR + scan;
	}

	/**
	 * Merges engine rows into candidates, preserving order of first appearance so
	 * downstream indices are reproducible across runs.
	 * <p>
	 * A PSM with no peptide or no file name cannot be keyed; it is counted into
	 * {@code dropped} rather than merged into a bucket it does not belong to, so
	 * the caller can report the shortfall instead of meeting it later as an
	 * unexplained missing count.
	 */
	public static List<CombinedPsmCandidate> merge(List<OpenSearchPsm> psms, int[] dropped) {
		Map<String, CombinedPsmCandidate> byKey = new LinkedHashMap<>(Math.max(16, psms.size() / 2));
		int skipped = 0;
		for (OpenSearchPsm psm : psms) {
			if (psm.getPeptide() == null || psm.getFileName() == null) {
				skipped++;
				continue;
			}
			String key = psm.getFileName() + KEY_SEPARATOR + psm.getScan() + KEY_SEPARATOR + psm.getPeptide();
			CombinedPsmCandidate candidate = byKey.get(key);
			if (candidate == null) {
				candidate = new CombinedPsmCandidate(psm.getFileName(), psm.getScan(), psm.getPeptide(),
						psm.isTarget());
				byKey.put(key, candidate);
			}
			if (psm.getEngine() == OpenSearchPsm.SearchEngine.SAGE) {
				// An engine can report one sequence twice when two modification
				// sites collapse to the same peptide string; keep its better rank.
				if (candidate.sage == null || psm.getHitRank() < candidate.sage.getHitRank()) {
					candidate.sage = psm;
				}
			} else {
				if (candidate.msfragger == null || psm.getHitRank() < candidate.msfragger.getHitRank()) {
					candidate.msfragger = psm;
				}
			}
		}
		if (dropped != null && dropped.length > 0) {
			dropped[0] = skipped;
		}
		return new ArrayList<>(byKey.values());
	}
}
