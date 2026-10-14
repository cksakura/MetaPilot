package uk.ac.quadram.metapilot.dbSearch.MaxQuant;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import uk.ac.quadram.metapilot.task.io.pep.MetaPSM;
import uk.ac.quadram.metapilot.task.io.qpx.EngineModificationParser;
import uk.ac.quadram.metapilot.task.io.qpx.QpxModification;
import uk.ac.quadram.metapilot.task.io.qpx.QpxScore;

/**
 * A MaxQuant {@code msms.txt} row as a {@link MetaPSM}, so the IGC MaxQuant workflow
 * can reach the common PSM layer (and the QPX PSM export) like every other DDA
 * workflow. Until now nothing read {@code msms.txt} into the pipeline.
 * <p>
 * msms.txt semantics: "Mass" is the theoretical neutral monoisotopic peptide mass
 * with all modifications; "Mass error [ppm]" is (observed - theoretical) /
 * theoretical; "Retention time" is in minutes; "PEP" is MaxQuant's posterior error
 * probability; "Reverse" = "+" marks decoys. msms.txt carries no per-PSM q-value.
 */
public class MaxquantMetaPsm extends MetaPSM {

	private final MaxquantMSMS source;

	public MaxquantMetaPsm(MaxquantMSMS source) {
		super(source.getRaw(), source.getScan(), source.getSequence(), source.getModseq(), source.getCharge(),
				splitProteins(source.getProtein()), source.getScore(), source.getPEP());
		this.source = source;
	}

	/** Every PSM of {@code msms.txt}, decoys included (flagged by {@link #isDecoy()}). */
	public static List<MetaPSM> read(File msmsTxt) {
		List<MetaPSM> list = new ArrayList<>();
		MaxquantMSMS[] rows = new MaxquantMSMSReader(msmsTxt).getMsms();
		if (rows != null) {
			for (MaxquantMSMS row : rows) {
				if (row != null) {
					list.add(new MaxquantMetaPsm(row));
				}
			}
		}
		return list;
	}

	private static String[] splitProteins(String proteins) {
		if (proteins == null || proteins.trim().isEmpty()) {
			return new String[0];
		}
		return proteins.split(";");
	}

	public MaxquantMSMS getSource() {
		return source;
	}

	@Override
	public double getExpMass() {
		return getObservedNeutralMass();
	}

	@Override
	public double getRt() {
		return source.getRt();
	}

	@Override
	public double getCalMass() {
		return source.getMass();
	}

	@Override
	public double getMassDiff() {
		return getObservedNeutralMass() - source.getMass();
	}

	@Override
	public int getMiss() {
		return source.getMiss();
	}

	// ------------------------------------------------------------------
	// Normalised accessors (see MetaPSM)
	// ------------------------------------------------------------------

	@Override
	public String getSearchEngine() {
		return "MaxQuant";
	}

	@Override
	public boolean isDecoy() {
		return !source.isTarget();
	}

	@Override
	public double getRtSeconds() {
		return source.getRt() * 60.0;
	}

	@Override
	public double getObservedNeutralMass() {
		double mass = source.getMass();
		double ppm = source.getMassError();
		if (mass <= 0 || Double.isNaN(ppm)) {
			return Double.NaN;
		}
		return mass * (1.0 + ppm / 1e6);
	}

	@Override
	public double getTheoreticalNeutralMass() {
		return source.getMass() > 0 ? source.getMass() : Double.NaN;
	}

	@Override
	public double getPosteriorErrorProbability() {
		return source.getPEP();
	}

	@Override
	public List<QpxModification> getQpxModifications() {
		return EngineModificationParser.fromMaxQuant(source.getModseq()).modifications;
	}

	@Override
	public List<QpxScore> getAdditionalScores() {
		List<QpxScore> scores = new ArrayList<>(2);
		QpxScore.addIfFinite(scores, "andromeda_score", source.getScore(), Boolean.TRUE);
		QpxScore.addIfFinite(scores, "andromeda_delta_score", source.getDeltaScore(), Boolean.TRUE);
		return scores;
	}
}
