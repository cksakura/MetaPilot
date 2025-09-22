package uk.ac.quadram.metapilot.dbSearch.fragpipe;

import uk.ac.quadram.metapilot.task.io.pep.MetaPSM;

/**
 * FragPipe PSM implementation extending MetaPSM
 */
public class FragPipePSM extends MetaPSM {
    
    private String spectrum;
    private double expMass;
    private double calMass;
    private double massDiff;
    private int miss;
    private double rt;
    private double probability;
    private double qValue;
    
    /**
     * Constructor for FragPipePSM
     * 
     * @param fileName      File name
     * @param spectrum      Spectrum identifier
     * @param sequence      Base peptide sequence
     * @param modSeq        Modified sequence
     * @param charge        Charge state
     * @param expMass       Experimental mass
     * @param calMass       Calculated mass
     * @param massDiff      Mass difference
     * @param miss          Number of missed cleavages
     * @param rt            Retention time
     * @param score         Hyperscore
     * @param probability   Probability value
     * @param qValue        Q-value
     * @param proteins      Protein identifiers
     */
    public FragPipePSM(String fileName, String spectrum, String sequence, String modSeq, int charge,
            double expMass, double calMass, double massDiff, int miss, double rt, double score,
            double probability, double qValue, String[] proteins) {
        super(fileName, sequence, modSeq, charge, proteins, score, qValue);
        this.spectrum = spectrum;
        this.expMass = expMass;
        this.calMass = calMass;
        this.massDiff = massDiff;
        this.miss = miss;
        this.rt = rt;
        this.probability = probability;
        this.qValue = qValue;
    }
    
    // Getters for FragPipe-specific fields
    public String getSpectrum() {
        return spectrum;
    }
    
    public double getProbability() {
        return probability;
    }
    
    public double getqValue() {
        return qValue;
    }
    
    // Implement abstract methods from MetaPSM
    @Override
    public double getExpMass() {
        return expMass;
    }
    
    @Override
    public double getRt() {
        return rt;
    }
    
    @Override
    public double getCalMass() {
        return calMass;
    }
    
    @Override
    public double getMassDiff() {
        return massDiff;
    }
    
    @Override
    public int getMiss() {
        return miss;
    }
    
    // Override getPEP to return qValue for FragPipe (similar to Sage)
    @Override
    public double getPEP() {
        return qValue;
    }
}
