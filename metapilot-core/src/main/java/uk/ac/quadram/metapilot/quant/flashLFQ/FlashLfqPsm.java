package uk.ac.quadram.metapilot.quant.flashLFQ;

import java.util.ArrayList;
import java.util.List;

import uk.ac.quadram.metapilot.task.io.pep.MetaPSM;
import uk.ac.quadram.metapilot.task.io.qpx.EngineModificationParser;
import uk.ac.quadram.metapilot.task.io.qpx.QpxModification;
import uk.ac.quadram.metapilot.task.io.qpx.QpxScore;

/**
 * A PSM read back from the FlashLFQ generic identification file MetaPilot writes for
 * pFind ({@code OpenPsmTsvWriter.addPFindPSM}) - the only producer of the files
 * {@link FlashLfqPsmReader} reads (DDA IGC/MAG pFind). The normalised accessors
 * below therefore describe pFind output: "Scan Retention Time" in minutes (pFind's
 * seconds / 60), neutral masses, a pFind full sequence ({@code LG(Gly->Pro)H...}),
 * and a q-value in the "PEP" column. pFind estimates no PEP.
 */
public class FlashLfqPsm extends MetaPSM {

	private static final String PFIND_DECOY_PREFIX = "REV_";
	
	private double expMass;
	private double rt;
	private double calMass;
	private double massDiff;
	private int miss;

	public FlashLfqPsm(String fileName, int scan, String baseSequence, String fullSequence, double expMass, double rt,
			int charge, String[] proteins, double score, double pEP, double calMass, double massDiff, int miss) {
		super(fileName, scan, baseSequence, fullSequence, charge, proteins, score, pEP);
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

	// ------------------------------------------------------------------
	// Normalised accessors (see MetaPSM and the class javadoc)
	// ------------------------------------------------------------------

	@Override
	public String getSearchEngine() {
		return "pFind";
	}

	/** pFind decoys carry REV_ on every protein (a shared target/decoy peptide is a target). */
	@Override
	public boolean isDecoy() {
		String[] pros = getProteins();
		if (pros == null || pros.length == 0) {
			return false;
		}
		for (String p : pros) {
			if (p == null || !p.trim().startsWith(PFIND_DECOY_PREFIX)) {
				return false;
			}
		}
		return true;
	}

	@Override
	public double getRtSeconds() {
		return rt * 60.0;
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
	public double getPsmQValue() {
		return pep;
	}

	@Override
	public List<QpxModification> getQpxModifications() {
		return EngineModificationParser.fromPFind(getFullSequence()).modifications;
	}

	@Override
	public List<QpxScore> getAdditionalScores() {
		List<QpxScore> scores = new ArrayList<>(2);
		QpxScore.addIfFinite(scores, "pfind_raw_score", score, Boolean.TRUE);
		QpxScore.addIfFinite(scores, "q_value", pep, Boolean.FALSE);
		return scores;
	}
}
