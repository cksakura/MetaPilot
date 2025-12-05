package uk.ac.quadram.metapilot.dbSearch;

/**
 * Abstract Peptide-Spectrum Match (PSM) Data Model
 * 
 * <h2>Purpose</h2>
 * This abstract class represents a peptide-spectrum match (PSM), which is the fundamental
 * unit of identification in mass spectrometry-based proteomics. It encapsulates all
 * essential information linking a theoretical peptide to an experimental spectrum.
 * 
 * <h2>PSM Components</h2>
 * <h3>Spectrum Information:</h3>
 * <ul>
 *   <li><strong>fileName:</strong> Source raw data file</li>
 *   <li><strong>scan:</strong> Scan number or spectrum identifier</li>
 *   <li><strong>rt:</strong> Retention time in minutes</li>
 *   <li><strong>charge:</strong> Precursor ion charge state</li>
 *   <li><strong>precorsorMr:</strong> Observed precursor mass</li>
 * </ul>
 * 
 * <h3>Peptide Information:</h3>
 * <ul>
 *   <li><strong>sequence:</strong> Peptide amino acid sequence</li>
 *   <li><strong>pepmass:</strong> Theoretical peptide mass</li>
 *   <li><strong>massdiff:</strong> Mass difference (observed - theoretical)</li>
 *   <li><strong>protein:</strong> Parent protein accession or identifier</li>
 * </ul>
 * 
 * <h3>Match Quality:</h3>
 * <ul>
 *   <li><strong>num_matched_ions:</strong> Number of matched fragment ions</li>
 *   <li><strong>hit_rank:</strong> Rank among candidate peptides</li>
 *   <li><strong>isTarget:</strong> Target vs decoy database hit</li>
 * </ul>
 * 
 * <h2>Target-Decoy Strategy</h2>
 * <p>The {@code isTarget} field is crucial for false discovery rate (FDR) estimation:</p>
 * <ul>
 *   <li><strong>Target Hits:</strong> Matches to real protein sequences</li>
 *   <li><strong>Decoy Hits:</strong> Matches to reversed/scrambled sequences</li>
 *   <li><strong>FDR Calculation:</strong> FDR = (2 × decoys) / (targets + decoys)</li>
 * </ul>
 * 
 * <h2>Mass Accuracy</h2>
 * <p>Mass difference (massdiff) represents measurement accuracy:</p>
 * <ul>
 *   <li><strong>Positive values:</strong> Observed mass > theoretical mass</li>
 *   <li><strong>Negative values:</strong> Observed mass < theoretical mass</li>
 *   <li><strong>Typical range:</strong> ±20 ppm for high-accuracy instruments</li>
 * </ul>
 * 
 * <h2>Search Engine Integration</h2>
 * <p>This class serves as the base for search engine-specific PSM implementations:</p>
 * <ul>
 *   <li>MaxQuant evidence and msms tables</li>
 *   <li>pFind spectrum identification results</li>
 *   <li>MSFragger pepXML output</li>
 *   <li>Sage TSV results</li>
 *   <li>DIA-NN report files</li>
 * </ul>
 * 
 * <h2>Quality Metrics</h2>
 * <p>Different search engines provide additional quality metrics that extend this base:</p>
 * <ul>
 *   <li><strong>Score-based:</strong> Xcorr, expect value, hyperscore</li>
 *   <li><strong>Probability-based:</strong> PeptideProphet probability, posterior error probability</li>
 *   <li><strong>Spectral-based:</strong> Spectral angle, dot product</li>
 * </ul>
 * 
 * <h2>Data Flow</h2>
 * <pre>{@code
 *   Raw Data → Search Engine → PSM → Filtering → Protein Inference
 * }</pre>
 * 
 * <h2>Usage Example</h2>
 * <pre>{@code
 *   // Create PSM instance (via concrete implementation)
 *   AbstractPSM psm = new MaxQuantPSM(
 *       "sample1.raw", 1234, 2, 800.4521, 15.67,
 *       "PEPTIDEK", 800.4512, 0.0009, 8, 1, 
 *       "sp|P12345|PROTEIN_HUMAN", true
 *   );
 *   
 *   // Access PSM properties
 *   double accuracy = psm.getMassdiff() / psm.getPepmass() * 1e6; // ppm
 *   boolean highConfidence = psm.getNum_matched_ions() >= 6;
 * }</pre>
 * 
 * <h2>Thread Safety</h2>
 * <p>This class is thread-safe for read operations after construction.
 * Modification operations should be synchronized externally if needed.</p>
 * 
 * @author Kai Cheng - Quadram Institute Bioscience
 * @version 2.0
 * @since 1.0
 * 
 * @see uk.ac.quadram.metapilot.dbSearch.MaxQuant.MaxQuantPSM
 * @see uk.ac.quadram.metapilot.dbSearch.pfind.PFindPSM
 * @see uk.ac.quadram.metapilot.dbSearch.sage.SagePSM
 */
public abstract class AbstractPSM {

	private String fileName;
	private int scan;
	private int charge;
	private double precorsorMr;
	private double rt;
	private String sequence;
	private double pepmass;
	private double massdiff;
	private int num_matched_ions;
	private int hit_rank;
	private String protein;
	private boolean isTarget;
	
	public AbstractPSM(String fileName, int scan, int charge, double precorsorMr, double rt, String sequence,
			double pepmass, double massdiff, int num_matched_ions, int hit_rank, String protein, boolean isTarget) {
		this.fileName = fileName;
		this.scan = scan;
		this.charge = charge;
		this.precorsorMr = precorsorMr;
		this.rt = rt;
		this.sequence = sequence;
		this.pepmass = pepmass;
		this.massdiff = massdiff;
		this.num_matched_ions = num_matched_ions;
		this.hit_rank = hit_rank;
		this.protein = protein;
		this.isTarget = isTarget;
	}

	public boolean isTarget() {
		return isTarget;
	}

	public void setTarget(boolean isTarget) {
		this.isTarget = isTarget;
	}

	public String getFileName() {
		return fileName;
	}

	public void setFileName(String fileName) {
		this.fileName = fileName;
	}

	public int getScan() {
		return scan;
	}

	public void setScan(int scan) {
		this.scan = scan;
	}

	public int getCharge() {
		return charge;
	}

	public void setCharge(int charge) {
		this.charge = charge;
	}

	public double getPrecorsorMr() {
		return precorsorMr;
	}

	public void setPrecorsorMr(double precorsorMr) {
		this.precorsorMr = precorsorMr;
	}

	public double getRt() {
		return rt;
	}

	public void setRt(double rt) {
		this.rt = rt;
	}

	public String getSequence() {
		return sequence;
	}

	public void setSequence(String sequence) {
		this.sequence = sequence;
	}

	public double getPepmass() {
		return pepmass;
	}

	public void setPepmass(double pepmass) {
		this.pepmass = pepmass;
	}

	public double getMassdiff() {
		return massdiff;
	}

	public void setMassdiff(double massdiff) {
		this.massdiff = massdiff;
	}

	public int getNum_matched_ions() {
		return num_matched_ions;
	}

	public void setNum_matched_ions(int num_matched_ions) {
		this.num_matched_ions = num_matched_ions;
	}

	public int getHit_rank() {
		return hit_rank;
	}

	public void setHit_rank(int hit_rank) {
		this.hit_rank = hit_rank;
	}

	public String getProtein() {
		return protein;
	}

	public void setProtein(String protein) {
		this.protein = protein;
	}

	

}
