package uk.ac.quadram.metapilot.task.io.pep;

import java.util.Collections;
import java.util.List;
import java.util.Map;

import uk.ac.quadram.metapilot.task.io.qpx.ProForma;
import uk.ac.quadram.metapilot.task.io.qpx.QpxModification;
import uk.ac.quadram.metapilot.task.io.qpx.QpxScore;

/**
 * Abstract base class for PSM (Peptide Spectrum Match) data. Provides the common
 * interface every DDA workflow's PSMs reach the post-analysis through, and - via the
 * "normalised accessors" below - the one contract the QPX PSM export
 * ({@code QpxPsmParquetWriter}) reads, so that it needs no knowledge of any engine.
 * <p>
 * <b>Legacy accessors</b> ({@link #getExpMass()}, {@link #getCalMass()},
 * {@link #getRt()}, {@link #getMassDiff()}, {@link #getPEP()}) keep their historical,
 * engine-specific meaning because the peptide/protein reports depend on them. In
 * particular {@link #getPEP()} returns whatever each engine put in the {@code pep}
 * field - a q-value for every engine MetaPilot currently reads - and the peptide
 * report's "PEP" column and best-PSM selection use it as such.
 * <p>
 * <b>Normalised accessors</b> have one meaning across engines, and return NaN /
 * null / empty when an engine cannot supply the value (never a stand-in):
 * <ul>
 * <li>{@link #getRtSeconds()} - seconds;</li>
 * <li>{@link #getObservedNeutralMass()} / {@link #getTheoreticalNeutralMass()} -
 * neutral monoisotopic masses, the theoretical one including every modification in
 * {@link #getProForma()};</li>
 * <li>{@link #getPosteriorErrorProbability()} - a real PEP, never a q-value;</li>
 * <li>{@link #getPsmQValue()} - the PSM-level q-value;</li>
 * <li>{@link #isDecoy()} - the search engine's own target/decoy decision (the
 * population its PSM FDR was estimated against). MetaPilot's shuffled genome-FDR
 * entries are targets here and are reported separately via
 * {@link #isGenomeFdrDecoy()}.</li>
 * </ul>
 */
public abstract class MetaPSM {

	protected String fileName;
	protected int scan;
	protected String baseSequence;
	protected String fullSequence;
	protected int charge;
	protected String[] proteins;
	protected double score;
	protected double pep; // legacy: a q-value for every engine - see class javadoc

	/**
	 * Constructor for MetaPSM
	 */
	public MetaPSM(String fileName, int scan, String baseSequence, String fullSequence, int charge, String[] proteins,
			double score, double pep) {
		this.fileName = fileName;
		this.scan = scan;
		this.baseSequence = baseSequence;
		this.fullSequence = fullSequence;
		this.charge = charge;
		this.proteins = proteins;
		this.score = score;
		this.pep = pep;
	}

	// Abstract methods that subclasses must implement for workflow-specific data
	public abstract double getExpMass();

	public abstract double getRt();

	public abstract double getCalMass();

	public abstract double getMassDiff();

	public abstract int getMiss();

	// Common getter methods
	public String getFileName() {
		return fileName;
	}

	public int getScan() {
		return scan;
	}

	public String getBaseSequence() {
		return baseSequence;
	}

	public String getFullSequence() {
		return fullSequence;
	}

	public String getModSeq() {
		return fullSequence;
	}

	public String getSequence() {
		return baseSequence;
	}

	public int getCharge() {
		return charge;
	}

	public String[] getProteins() {
		return proteins;
	}

	/**
	 * Sets the protein assignments for this PSM (e.g. after matching to a MAG database).
	 */
	public void setProteins(String[] proteins) {
		this.proteins = proteins;
	}

	public double getScore() {
		return score;
	}

	/** Legacy: see the class javadoc - a q-value for every current engine. */
	public double getPEP() {
		return pep;
	}

	/** Legacy: the same value as {@link #getPEP()}. */
	public double getqValue() {
		return pep;
	}

	// ------------------------------------------------------------------
	// Normalised accessors (QPX PSM export). Defaults mean "not available".
	// ------------------------------------------------------------------

	/** Search engine that produced the PSM, e.g. {@code MSFragger}, {@code Sage}, {@code pFind}. */
	public String getSearchEngine() {
		return null;
	}

	/** The engine's own decoy decision; see the class javadoc. */
	public boolean isDecoy() {
		return false;
	}

	/**
	 * True when every protein of the PSM is one of MetaPilot's shuffled
	 * ({@code shuf_}/{@code rev_shuf_}) genome-FDR entries - targets to the engine,
	 * decoys to the genome-level FDR.
	 */
	public boolean isGenomeFdrDecoy() {
		String[] pros = getProteins();
		if (pros == null || pros.length == 0) {
			return false;
		}
		for (String p : pros) {
			if (p == null) {
				return false;
			}
			String t = p.trim().toLowerCase(java.util.Locale.ROOT);
			if (!(t.startsWith("shuf_") || t.startsWith("rev_shuf_"))) {
				return false;
			}
		}
		return true;
	}

	/** Spectrum (run) file name without directory or extension. */
	public String getSpectrumFileName() {
		return stripPathAndExtension(fileName);
	}

	/** MS2 retention time in seconds; NaN if unknown. */
	public double getRtSeconds() {
		return Double.NaN;
	}

	/** Engine-predicted retention time in seconds; NaN if not predicted or not in seconds. */
	public double getPredictedRtSeconds() {
		return Double.NaN;
	}

	/** Precursor ion mobility; NaN if not measured. */
	public double getIonMobility() {
		return Double.NaN;
	}

	/** Experimental precursor neutral monoisotopic mass (Da); NaN if unknown. */
	public double getObservedNeutralMass() {
		return Double.NaN;
	}

	/**
	 * Theoretical neutral monoisotopic mass (Da) of the peptidoform returned by
	 * {@link #getProForma()}; NaN if unknown.
	 */
	public double getTheoreticalNeutralMass() {
		return Double.NaN;
	}

	/** Posterior error probability; NaN when the engine does not estimate one. */
	public double getPosteriorErrorProbability() {
		return Double.NaN;
	}

	/** PSM-level q-value; NaN when the engine does not report one. */
	public double getPsmQValue() {
		return Double.NaN;
	}

	/**
	 * Modifications in QPX form, or null when this PSM's notation cannot be parsed
	 * (the ProForma then falls back to the bare sequence and the export records that
	 * the modifications are unknown rather than absent).
	 */
	public List<QpxModification> getQpxModifications() {
		return null;
	}

	/** ProForma 2.0 peptidoform. */
	public String getProForma() {
		List<QpxModification> mods = getQpxModifications();
		return ProForma.build(getBaseSequence(), mods == null ? Collections.emptyList() : mods);
	}

	/** Engine scores beyond the core QPX columns. */
	public List<QpxScore> getAdditionalScores() {
		return Collections.emptyList();
	}

	/** Extra annotations for the QPX {@code cv_params} column (name -> value). */
	public Map<String, String> getCvParams() {
		return Collections.emptyMap();
	}

	protected static String stripPathAndExtension(String name) {
		if (name == null) {
			return null;
		}
		String n = name;
		int slash = Math.max(n.lastIndexOf('/'), n.lastIndexOf('\\'));
		if (slash >= 0) {
			n = n.substring(slash + 1);
		}
		String lower = n.toLowerCase(java.util.Locale.ROOT);
		for (String ext : new String[] { ".mzml", ".mzxml", ".mgf", ".raw", ".d", ".wiff", ".ms2", ".pf2" }) {
			if (lower.endsWith(ext)) {
				return n.substring(0, n.length() - ext.length());
			}
		}
		return n;
	}

}
