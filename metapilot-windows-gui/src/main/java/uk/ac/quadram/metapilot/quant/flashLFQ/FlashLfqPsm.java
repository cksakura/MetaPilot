package uk.ac.quadram.metapilot.quant.flashLFQ;

import uk.ac.quadram.metapilot.task.io.pep.MetaPSM;

public class FlashLfqPsm extends MetaPSM {
	
	private double expMass;
	private double rt;
	private double calMass;
	private double massDiff;
	private int miss;

	public FlashLfqPsm(String fileName, String baseSequence, String fullSequence, double expMass, double rt,
			int charge, String[] proteins, double score, double pEP, double calMass, double massDiff, int miss) {
		super(fileName, baseSequence, fullSequence, charge, proteins, score, pEP);
		this.expMass = expMass;
		this.rt = rt;
		this.calMass = calMass;
		this.massDiff = massDiff;
		this.miss = miss;
	}

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
}
