/**
 * 
 */
package uk.ac.quadram.metapilot.dbSearch.pepxml;

import uk.ac.quadram.metapilot.dbSearch.AbstractPSM;

/**
 * @author Kai Cheng
 *
 */
public class SequestHTPSM extends AbstractPSM {
	
	private double xcorr;
	private double probability;

	public SequestHTPSM(String fileName, int scan, int charge, double precorsorMr, double rt, String sequence,
			double pepmass, double massdiff, int num_matched_ions, int hit_rank, String protein, double xcorr,
			double probability, boolean isTarget) {
		super(fileName, scan, charge, precorsorMr, rt, sequence, pepmass, massdiff, num_matched_ions, hit_rank, protein, isTarget);
		this.xcorr = xcorr;
		this.probability = probability;
	}

	public double getXcorr() {
		return xcorr;
	}

	public void setXcorr(double xcorr) {
		this.xcorr = xcorr;
	}

	public double getProbability() {
		return probability;
	}

	public void setProbability(double probability) {
		this.probability = probability;
	}

}
