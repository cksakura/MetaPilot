package uk.ac.quadram.metapilot.dbSearch.sage;

import java.util.ArrayList;
import java.util.List;

import uk.ac.quadram.metapilot.task.io.pep.MetaPSM;
import uk.ac.quadram.metapilot.task.io.qpx.EngineModificationParser;
import uk.ac.quadram.metapilot.task.io.qpx.QpxModification;
import uk.ac.quadram.metapilot.task.io.qpx.QpxScore;

public class SagePSM extends MetaPSM {

	private double expMass;
	private double calMass;
	private double massPPM;
	private double fragPPM;
	private int miss;
	private int specificity;
	private double hyperscore;
	private double disScore;
	private double rt;
	private boolean isTarget;
	private double peptideQ;
	private double proteinQ;
	private double ms2Intensity;
	// Optional results.sage.tsv columns (defaults when absent)
	private int rank = 1;
	/** Sage writes posterior_error as log10(PEP) - linear_discriminant.rs: kde.posterior_error(..).log10(). */
	private double posteriorErrorLog10 = Double.NaN;
	private double ionMobility = Double.NaN;

	public SagePSM(String fileName, int scan, String sequence, String modSeq, double expMass, double calMass,
			double massPPM, double fragPPM, int charge, int miss, int specificity, double rt, double qValue,
			double hyperscore, double disScore, String[] proteins, boolean isTarget, double peptideQ, double proteinQ,
			double ms2Intensity) {
		super(fileName, scan, sequence, modSeq, charge, proteins, hyperscore, qValue);
		this.expMass = expMass;
		this.calMass = calMass;
		this.massPPM = massPPM;
		this.fragPPM = fragPPM;
		this.miss = miss;
		this.specificity = specificity;
		this.rt = rt;
		this.hyperscore = hyperscore;
		this.disScore = disScore;
		this.isTarget = isTarget;
		this.peptideQ = peptideQ;
		this.proteinQ = proteinQ;
		this.ms2Intensity = ms2Intensity;
	}

	public String getProtein() {
		StringBuilder sb = new StringBuilder();
		for (String pro : proteins) {
			sb.append(pro).append(";");
		}
		if (sb.length() > 0) {
			return sb.substring(0, sb.length() - 1);
		} else {
			return "";
		}
	}

	public int getScan() {
		return scan;
	}

	public void setScan(int scan) {
		this.scan = scan;
	}

	@Override
	public double getExpMass() {
		return expMass;
	}

	public void setExpMass(double expMass) {
		this.expMass = expMass;
	}

	@Override
	public double getCalMass() {
		return calMass;
	}

	public void setCalMass(double calMass) {
		this.calMass = calMass;
	}

	public double getMassPPM() {
		return massPPM;
	}

	public void setMassPPM(double massPPM) {
		this.massPPM = massPPM;
	}

	public double getFragPPM() {
		return fragPPM;
	}

	public void setFragPPM(double fragPPM) {
		this.fragPPM = fragPPM;
	}

	@Override
	public int getMiss() {
		return miss;
	}

	public void setMiss(int miss) {
		this.miss = miss;
	}

	public int getSpecificity() {
		return specificity;
	}

	public void setSpecificity(int specificity) {
		this.specificity = specificity;
	}

	public double getHyperscore() {
		return hyperscore;
	}

	public void setHyperscore(double hyperscore) {
		this.hyperscore = hyperscore;
	}

	public double getDisScore() {
		return disScore;
	}

	public void setDisScore(double disScore) {
		this.disScore = disScore;
	}

	@Override
	public double getRt() {
		return rt;
	}

	public void setRt(double rt) {
		this.rt = rt;
	}

	/**
	 * Precursor mass difference in Da (expmass - calcmass), the unit every other
	 * engine's getMassDiff() uses. Previously returned precursor_ppm; no caller read
	 * it for Sage PSMs. The ppm value is {@link #getMassPPM()}.
	 */
	@Override
	public double getMassDiff() {
		return expMass - calMass;
	}

	public boolean isTarget() {
		return isTarget;
	}

	public void setTarget(boolean isTarget) {
		this.isTarget = isTarget;
	}

	public double getPeptideQ() {
		return peptideQ;
	}

	public void setPeptideQ(double peptideQ) {
		this.peptideQ = peptideQ;
	}

	public double getProteinQ() {
		return proteinQ;
	}

	public void setProteinQ(double proteinQ) {
		this.proteinQ = proteinQ;
	}

	public double getMs2Intensity() {
		return ms2Intensity;
	}

	public void setMs2Intensity(double ms2Intensity) {
		this.ms2Intensity = ms2Intensity;
	}

	public int getRank() {
		return rank;
	}

	public void setRank(int rank) {
		this.rank = rank;
	}

	public double getPosteriorErrorLog10() {
		return posteriorErrorLog10;
	}

	public void setPosteriorErrorLog10(double posteriorErrorLog10) {
		this.posteriorErrorLog10 = posteriorErrorLog10;
	}

	public void setIonMobility(double ionMobility) {
		this.ionMobility = ionMobility;
	}

	// ------------------------------------------------------------------
	// Normalised accessors (see MetaPSM). results.sage.tsv semantics: rt is in
	// minutes (measured: FragPipe "Retention" / Sage rt = 60.0 over 63,544 shared
	// scans of 12mix YS-194-1-1); expmass/calcmass are neutral masses with every
	// modification in the peptide string; spectrum_q is the PSM q-value.
	// ------------------------------------------------------------------

	@Override
	public String getSearchEngine() {
		return "Sage";
	}

	@Override
	public boolean isDecoy() {
		return !isTarget;
	}

	@Override
	public double getRtSeconds() {
		return rt * 60.0;
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
		if (Double.isNaN(posteriorErrorLog10)) {
			return Double.NaN;
		}
		return Math.min(1.0, Math.pow(10.0, posteriorErrorLog10));
	}

	@Override
	public double getPsmQValue() {
		return getqValue();
	}

	@Override
	public List<QpxModification> getQpxModifications() {
		return EngineModificationParser.fromMassBrackets(getFullSequence()).modifications;
	}

	@Override
	public List<QpxScore> getAdditionalScores() {
		List<QpxScore> scores = new ArrayList<>(6);
		QpxScore.addIfFinite(scores, "sage_hyperscore", hyperscore, Boolean.TRUE);
		QpxScore.addIfFinite(scores, "sage_discriminant_score", disScore, Boolean.TRUE);
		QpxScore.addIfFinite(scores, "q_value", getqValue(), Boolean.FALSE);
		QpxScore.addIfFinite(scores, "sage_peptide_q", peptideQ, Boolean.FALSE);
		QpxScore.addIfFinite(scores, "sage_protein_q", proteinQ, Boolean.FALSE);
		QpxScore.addIfFinite(scores, "rank", rank, Boolean.FALSE);
		return scores;
	}

}
