package uk.ac.quadram.metapilot.core.opensearch;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * One (spectrum, peptide) hypothesis of the joint closed+open identification, carrying every search
 * that proposed it: the closed MSFragger search, the open MSFragger search and, when it ran, Sage's
 * open search.
 * <p>
 * The closed and open searches of a sample run against the same refined database, so a peptide
 * both propose for one spectrum is one hypothesis with two bodies of evidence, not two rows
 * competing for the spectrum. The same peptide on the same spectrum also has the same precursor
 * mass, so the open delta is the closed interpretation (variable modifications plus isotope
 * offset) and the two never disagree on what the peptide is.
 * <p>
 * The FDR {@link Group} is a property of the candidate's evidence, decided before scoring and in
 * the same way for targets and decoys, which is what makes a target-decoy walk inside a group
 * valid.
 */
public final class JointCandidate {

	/** FDR group: which kind of evidence identifies the spectrum. Order is priority, highest first. */
	public enum Group {
		/** Proposed by the closed search: the peptide explains the precursor at exact mass. */
		EXACT,
		/** Open search only, on the zero family (unmodified, isotope errors folded in). */
		OPEN_UNMODIFIED,
		/** Open search only, on a delta-mass peak accepted by the evidence rule ({@link ModificationValidator#filterPeaks}). */
		OPEN_ACCEPTED_SHIFT,
		/**
		 * Open search only, any other delta that the fragments place on one residue, that a UNIMOD
		 * entry explains there, and that no substitution could produce there ({@link ShiftLocaliser},
		 * identification v2). Assigned to spectrum winners after the model has scored, so the model
		 * never sees it.
		 */
		OPEN_LOCALISED,
		/**
		 * As {@link #OPEN_LOCALISED}, on a residue that this delta could turn into another amino acid.
		 * The fragments cannot tell a modification from a sequence difference there, and on 12mix 17.7%
		 * of such PSMs lay outside the true genomes against 4.0% otherwise, unseen by decoys. Walked
		 * and reported, and read by the variant outputs, but not an identification. An accepted-shift
		 * winner explained on such a residue is moved here too (2026-11-10): one rule on and off a peak.
		 */
		OPEN_SUBSTITUTION_SITE,
		/** Open search only, any other delta: walked and reported, not an identification by default. */
		OPEN_UNEXPLAINED
	}

	private final String fileName;
	private final int scan;
	private final String peptide;
	private final boolean target;

	OpenSearchPsm closed;
	OpenSearchPsm fragger;
	OpenSearchPsm sage;
	Group group;
	/**
	 * Precursor queries (identification v2): the open rows explain a precursor Crystal-C reassigned
	 * away from the closed search's (another charge, or a mass not within 0.02 Da of the closed one
	 * after whole isotope spacings). Such a candidate is its own query and never merges with a closed
	 * row.
	 */
	boolean reassigned;
	/** Shift localisation of an unexplained-shift spectrum winner (identification v2); null when not run. */
	ShiftLocaliser.Result localisation;

	JointCandidate(String fileName, int scan, String peptide, boolean target) {
		this.fileName = fileName;
		this.scan = scan;
		this.peptide = peptide;
		this.target = target;
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

	public Group getGroup() {
		return group;
	}

	public OpenSearchPsm closed() {
		return closed;
	}

	public OpenSearchPsm fragger() {
		return fragger;
	}

	public OpenSearchPsm sage() {
		return sage;
	}

	public boolean foundByClosed() {
		return closed != null;
	}

	public boolean foundByOpen() {
		return fragger != null || sage != null;
	}

	/** The open row carrying the delta-mass peak context: MSFragger's when present. */
	public OpenSearchPsm openPrimary() {
		return fragger != null ? fragger : sage;
	}

	/**
	 * The row reported for this candidate. The closed row when the closed search proposed it: its
	 * interpretation is exact (variable modifications, isotope offset), which is what the
	 * quantification and the reports need.
	 */
	public OpenSearchPsm representative() {
		return closed != null ? closed : openPrimary();
	}

	/** Which searches proposed this candidate, e.g. {@code CLOSED+MSFRAGGER}. */
	public String support() {
		StringBuilder sb = new StringBuilder();
		if (closed != null) {
			sb.append("CLOSED");
		}
		if (fragger != null) {
			sb.append(sb.length() > 0 ? "+" : "").append("MSFRAGGER");
		}
		if (sage != null) {
			sb.append(sb.length() > 0 ? "+" : "").append("SAGE");
		}
		return sb.toString();
	}

	/**
	 * The fold and collapse key. Same string as {@link SemiSupervisedPsmClassifier}'s fold key, so a
	 * spectrum lands in the same fold whichever classifier scores it.
	 */
	public String spectrumKey() {
		return fileName + ":" + scan;
	}

	public ShiftLocaliser.Result getLocalisation() {
		return localisation;
	}

	public boolean isReassigned() {
		return reassigned;
	}

	/**
	 * The competition key: one winner per query. The spectrum, or the spectrum's reassigned precursor
	 * (identification v2). Folds stay keyed on {@link #spectrumKey}, so both queries of a scan share a
	 * fold.
	 */
	public String queryKey() {
		return reassigned ? spectrumKey() + ":R" : spectrumKey();
	}

	/**
	 * Merges open rows (either engine, any rank) and closed rows into candidates, in order of first
	 * appearance. A row without a peptide or file name is counted into {@code dropped[0]}; a row whose
	 * target/decoy label differs from the candidate's (first row's) label into {@code dropped[1]}
	 * when the array has room.
	 */
	public static List<JointCandidate> merge(List<OpenSearchPsm> open, List<OpenSearchPsm> closed, int[] dropped) {
		return merge(open, closed, dropped, false);
	}

	/** Mass tolerance of "the same precursor" after whole isotope spacings, Da. */
	static final double SAME_PRECURSOR_DA = 0.02;
	/** Isotope spacings ({@code |k|}) within which a precursor mass is still the same precursor. */
	static final int SAME_PRECURSOR_MAX_ISOTOPES = 3;

	/**
	 * As {@link #merge(List, List, int[])}; with {@code byPrecursor}, an open row is merged only with
	 * rows explaining the same precursor ({@link #reassigned}).
	 *
	 * @throws IllegalStateException when a scan's open rows carry more than one reassigned precursor
	 *                               (a query key holds one)
	 */
	public static List<JointCandidate> merge(List<OpenSearchPsm> open, List<OpenSearchPsm> closed, int[] dropped,
			boolean byPrecursor) {
		int size = (open == null ? 0 : open.size()) + (closed == null ? 0 : closed.size());
		Map<String, JointCandidate> byKey = new LinkedHashMap<>(Math.max(16, size));
		Map<String, double[]> closedPrecursor = new java.util.HashMap<>();
		Map<String, String> reassignedPrecursor = new java.util.HashMap<>();
		if (byPrecursor && closed != null) {
			for (OpenSearchPsm psm : closed) {
				if (psm.getFileName() != null) {
					closedPrecursor.putIfAbsent(psm.getFileName() + ":" + psm.getScan(),
							new double[] { psm.getPrecursorNeutralMass(), psm.getCharge() });
				}
			}
		}
		int skipped = 0;
		int labelConflicts = 0;
		int multiplePrecursors = 0;
		for (List<OpenSearchPsm> rows : java.util.Arrays.asList(open, closed)) {
			if (rows == null) {
				continue;
			}
			for (OpenSearchPsm psm : rows) {
				if (psm.getPeptide() == null || psm.getFileName() == null) {
					skipped++;
					continue;
				}
				String scanKey = psm.getFileName() + ":" + psm.getScan();
				boolean reassigned = byPrecursor && !psm.isClosedStage()
						&& !samePrecursor(closedPrecursor.get(scanKey), psm);
				if (reassigned) {
					String precursor = String.format(java.util.Locale.ROOT, "%.3f/%d", psm.getPrecursorNeutralMass(),
							psm.getCharge());
					String held = reassignedPrecursor.putIfAbsent(scanKey, precursor);
					if (held != null && !held.equals(precursor)) {
						multiplePrecursors++;
					}
				}
				String key = psm.getFileName() + '\t' + psm.getScan() + '\t' + psm.getPeptide() + (reassigned ? "\tR" : "");
				JointCandidate c = byKey.get(key);
				if (c == null) {
					c = new JointCandidate(psm.getFileName(), psm.getScan(), psm.getPeptide(), psm.isTarget());
					c.reassigned = reassigned;
					byKey.put(key, c);
				} else if (c.target != psm.isTarget()) {
					labelConflicts++;
				}
				if (psm.isClosedStage()) {
					if (c.closed == null || psm.getHitRank() < c.closed.getHitRank()) {
						c.closed = psm;
					}
				} else if (psm.getEngine() == OpenSearchPsm.SearchEngine.SAGE) {
					if (c.sage == null || psm.getHitRank() < c.sage.getHitRank()) {
						c.sage = psm;
					}
				} else if (c.fragger == null || psm.getHitRank() < c.fragger.getHitRank()) {
					c.fragger = psm;
				}
			}
		}
		if (multiplePrecursors > 0) {
			throw new IllegalStateException(multiplePrecursors + " open rows carry a second reassigned precursor for "
					+ "their scan; a precursor query holds one reassigned precursor per scan");
		}
		if (dropped != null && dropped.length > 0) {
			dropped[0] = skipped;
		}
		if (dropped != null && dropped.length > 1) {
			dropped[1] = labelConflicts;
		}
		return new ArrayList<>(byKey.values());
	}

	/**
	 * Whether an open row explains the closed search's precursor ({mass, charge}; null when the scan has
	 * no closed row, which makes the open precursor the only one): same charge, and neutral masses
	 * within {@link #SAME_PRECURSOR_DA} after a whole number of isotope spacings.
	 */
	static boolean samePrecursor(double[] closedPrecursor, OpenSearchPsm open) {
		if (closedPrecursor == null) {
			return true;
		}
		double om = open.getPrecursorNeutralMass();
		if (!Double.isFinite(om) || !Double.isFinite(closedPrecursor[0])) {
			return true;
		}
		if (open.getCharge() != (int) closedPrecursor[1]) {
			return false;
		}
		double dm = om - closedPrecursor[0];
		long k = Math.round(dm / ClosedOpenCascade.ISOTOPE_SPACING);
		return Math.abs(k) <= SAME_PRECURSOR_MAX_ISOTOPES
				&& Math.abs(dm - k * ClosedOpenCascade.ISOTOPE_SPACING) <= SAME_PRECURSOR_DA;
	}
}
