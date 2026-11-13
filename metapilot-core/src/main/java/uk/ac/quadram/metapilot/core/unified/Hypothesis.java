package uk.ac.quadram.metapilot.core.unified;

/**
 * One (spectrum, peptide) hypothesis of the unified identification
 * ({@code docs/developer/unified-rescoring-batchwise-fdr-design.md}), with every search that
 * proposed it merged in. Holds only what the identification needs, so the search rows can be
 * released once the candidate sets are built.
 */
public final class Hypothesis {

	/** How the hypothesis explains the precursor mass. Decided from mass alone, before any label is read. */
	public enum Explanation {
		/** Proposed by the closed search, or an open delta inside the zero family. */
		EXACT,
		/** An open delta on a retained delta-mass peak. */
		PEAK,
		/** An open delta on no retained peak. */
		UNEXPLAINED
	}

	final int spectrum;
	final String peptide;
	final boolean target;
	final int charge;
	final double precursorMass;
	final double retentionTime;
	final int numTolTerm;
	final int missedCleavages;
	final String[] proteins;

	boolean foundByClosed;
	boolean foundByOpen;
	/** Closed search: best (lowest) rank; Integer.MAX_VALUE when not proposed by the closed search. */
	int closedRank = Integer.MAX_VALUE;
	/** Closed search: mass difference of its best-ranked row (precursor - peptide incl. variable modifications). */
	double closedMassDiff = Double.NaN;
	/** Closed search: variable modifications, index 0 the N-terminus, i the residue i (1-based). */
	double[] closedMods;
	/** Open search: raw mass difference of its best-ranked row. */
	double openMassDiff = Double.NaN;
	/** Open search: best (lowest) rank of its rows; Integer.MAX_VALUE when not proposed by the open search. */
	int openRank = Integer.MAX_VALUE;
	/** Open search row index into the peak fit, or -1. */
	int openFitIndex = -1;

	Explanation explanation;
	/** Retained peak (batch) id for {@link Explanation#PEAK}, else -1. */
	int peakId = -1;
	/** Modification mass placed on the fragments: the isotope-normalised delta, 0 when none. */
	double fragmentDelta;
	int isotopeOffset;

	/** Annotation features ({@link SpectrumAnnotator.Annotation}), null until annotated. */
	SpectrumAnnotator.Annotation annotation;
	/** Dropped by neighbour filtering. */
	boolean neighbour;
	/** RESET coin, assigned once: a decoy that would be a training decoy if it wins its spectrum. */
	boolean resetTraining;

	/**
	 * Diagnostic only ({@code UnifiedIdentification.Config#engineFeatures}): the closed classifier's own
	 * inputs for this candidate's closed row, so the training machinery can be checked against the
	 * reference result. Not part of the method.
	 */
	float[] engine;
	String protein;

	/** DIA-NN backbone prediction against the observed spectrum (delta not placed); NaN when not predicted. */
	float predEntropy = Float.NaN;
	float predPearson = Float.NaN;
	float predObserved = Float.NaN;
	float predRtError = Float.NaN;
	double predIrt = Double.NaN;

	Hypothesis(int spectrum, String peptide, boolean target, int charge, double precursorMass, double retentionTime,
			int numTolTerm, int missedCleavages, String[] proteins) {
		this.spectrum = spectrum;
		this.peptide = peptide;
		this.target = target;
		this.charge = charge;
		this.precursorMass = precursorMass;
		this.retentionTime = retentionTime;
		this.numTolTerm = numTolTerm;
		this.missedCleavages = missedCleavages;
		this.proteins = proteins;
	}

	public String peptide() {
		return peptide;
	}

	public boolean isTarget() {
		return target;
	}

	public Explanation explanation() {
		return explanation;
	}

	/** Base sequence with I merged into L: the peptide unit. */
	static String unit(String peptide) {
		return peptide.toUpperCase(java.util.Locale.ROOT).replace('I', 'L');
	}
}
