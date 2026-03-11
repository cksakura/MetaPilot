package uk.ac.quadram.metapilot.task.io.pep;

/**
 * Abstract base class for PSM (Peptide Spectrum Match) data Provides common
 * interface for different PSM implementations
 */
public abstract class MetaPSM {

	protected String fileName;
	protected int scan;
	protected String baseSequence;
	protected String fullSequence;
	protected int charge;
	protected String[] proteins;
	protected double score;
	protected double pep; // PEP or q-value

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

	public double getPEP() {
		return pep;
	}

	public double getqValue() {
		return pep;
	}

}
