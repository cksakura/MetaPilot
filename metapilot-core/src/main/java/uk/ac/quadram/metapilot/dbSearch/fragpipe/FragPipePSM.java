package uk.ac.quadram.metapilot.dbSearch.fragpipe;

import java.util.ArrayList;
import java.util.List;

import uk.ac.quadram.metapilot.task.io.pep.MetaPSM;
import uk.ac.quadram.metapilot.task.io.qpx.EngineModificationParser;
import uk.ac.quadram.metapilot.task.io.qpx.QpxModification;
import uk.ac.quadram.metapilot.task.io.qpx.QpxScore;

/**
 * FragPipe PSM implementation extending MetaPSM
 */
public class FragPipePSM extends MetaPSM {
    
    private String spectrum;
    private String mod;
    private double expMass;
    private double calMass;
    private double massDiff;
    private int miss;
    private double rt;
    private double probability;
    private double qValue;
    private boolean decoy;
    // Optional psm.tsv columns, read when present (NaN = column absent)
    private double expectation = Double.NaN;
    private double nextscore = Double.NaN;
    private double ionMobility = Double.NaN;
    
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
    public FragPipePSM(String fileName, int scan, String spectrum, String sequence, String modSeq, String mod, int charge,
            double expMass, double calMass, double massDiff, int miss, double rt, double score,
			double probability, double qValue, String[] proteins) {
		this(fileName, scan, spectrum, sequence, modSeq, mod, charge, expMass, calMass, massDiff, miss, rt, score,
				probability, qValue, proteins, false);
	}

	public FragPipePSM(String fileName, int scan, String spectrum, String sequence, String modSeq, String mod, int charge,
			double expMass, double calMass, double massDiff, int miss, double rt, double score,
			double probability, double qValue, String[] proteins, boolean decoy) {
		super(fileName, scan, sequence, modSeq, charge, proteins, score, qValue);
		this.spectrum = spectrum;
		this.mod = mod;
		this.expMass = expMass;
		this.calMass = calMass;
		this.massDiff = massDiff;
		this.miss = miss;
		this.rt = rt;
		this.probability = probability;
		this.qValue = qValue;
		this.decoy = decoy;
	}
    
    // Getters for FragPipe-specific fields
    public String getSpectrum() {
        return spectrum;
    }
    
    public String getMod() {
		return mod;
	}

	public double getProbability() {
        return probability;
    }
    
    public double getqValue() {
        return qValue;
    }

	public boolean isDecoy() {
		return decoy;
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

    public void setExpectation(double expectation) {
        this.expectation = expectation;
    }

    public void setNextscore(double nextscore) {
        this.nextscore = nextscore;
    }

    public void setIonMobility(double ionMobility) {
        this.ionMobility = ionMobility;
    }

    // ------------------------------------------------------------------
    // Normalised accessors (see MetaPSM). psm.tsv semantics: "Retention" is in
    // seconds; "Observed Mass" / "Calculated Peptide Mass" are neutral masses, the
    // calculated one including every assigned modification; "Probability" is the
    // posterior probability the PSM is correct (PeptideProphet or Percolator), so
    // PEP = 1 - Probability.
    // ------------------------------------------------------------------

    @Override
    public String getSearchEngine() {
        return "MSFragger";
    }

    @Override
    public String getSpectrumFileName() {
        // "Spectrum" is <run>.<scan>.<scan>.<charge>; the run name is authoritative
        // (the IGC reader's fileName is the experiment folder, the MAG reader's is the
        // pep.xml path).
        if (spectrum != null) {
            String[] parts = spectrum.split("\\.");
            if (parts.length >= 4) {
                StringBuilder sb = new StringBuilder(parts[0]);
                for (int i = 1; i < parts.length - 3; i++) {
                    sb.append('.').append(parts[i]);
                }
                return sb.toString();
            }
        }
        return super.getSpectrumFileName();
    }

    @Override
    public double getRtSeconds() {
        return rt;
    }

    @Override
    public double getIonMobility() {
        return ionMobility > 0 ? ionMobility : Double.NaN;
    }

    @Override
    public double getObservedNeutralMass() {
        return expMass;
    }

    @Override
    public double getTheoreticalNeutralMass() {
        return calMass;
    }

    @Override
    public double getPosteriorErrorProbability() {
        if (Double.isNaN(probability) || probability < 0 || probability > 1) {
            return Double.NaN;
        }
        return 1.0 - probability;
    }

    @Override
    public double getPsmQValue() {
        return qValue;
    }

    @Override
    public List<QpxModification> getQpxModifications() {
        return EngineModificationParser.fromFragPipeAssigned(getBaseSequence(), mod);
    }

    @Override
    public List<QpxScore> getAdditionalScores() {
        List<QpxScore> scores = new ArrayList<>(4);
        QpxScore.addIfFinite(scores, "msfragger_hyperscore", score, Boolean.TRUE);
        QpxScore.addIfFinite(scores, "msfragger_nextscore", nextscore, Boolean.TRUE);
        QpxScore.addIfFinite(scores, "msfragger_expect", expectation, Boolean.FALSE);
        QpxScore.addIfFinite(scores, "q_value", qValue, Boolean.FALSE);
        return scores;
    }
}
