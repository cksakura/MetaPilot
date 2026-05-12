package uk.ac.quadram.metapilot.task.dia.par;

import uk.ac.quadram.metapilot.task.par.MetaBaseParameter;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;

/**
 * Parameter class for MAG-based DIA marker/seed workflows.
 *
 * Shared by:
 * - DIA Marker workflows (Focus / Comprehensive, and Deep Marker)
 * - DIA Seed iterative workflows (Focus / Comprehensive)
 * - Deep comprehensive DIA workflow
 *
 * Focus / Comprehensive / Deep are encoded via {@link #markerMode}.
 */
public class MetaParameterDiaMarker extends MetaParameterDia {

	/**
	 * Marker workflow mode persisted for UI/CLI round-trip.
	 * <p>
	 * Known values:
	 * - FOCUS
	 * - COMPREHENSIVE
	 * - DEEP (iterative rounds + combined search over all DIA files)
	 * - DEEP_PER_FILE (iterative workflow run independently per DIA file; results under marker/<rawname>/...)
	 */
	private String markerMode = "FOCUS";

	// Modification settings (from DiannParameter)
	private boolean mExcision = true;
	private boolean carbamiC = true;
	private boolean oxidationM = false;
	private boolean acetyProN = false;
	private boolean phosphoSTY = false;
	private boolean kGG = false;

	/**
	 * Deep iterative workflows: genomes included per round.
	 * <p>
	 * Default is 5 to match legacy behaviour. Set to 10 (etc.) to benchmark speed/sensitivity
	 * trade-offs for single-file or per-file modes.
	 */
	private int genomesPerRound = 5;

	public MetaParameterDiaMarker() {
		super();
		this.setWorkflowType(uk.ac.quadram.metapilot.workflow.WorkflowType.DIA_MARKER.getAnalysisType());
	}

	public MetaParameterDiaMarker(MetaBaseParameter metaPar, MagDbItem magDbItem, boolean librarySearch,
			boolean selfModelSearch, String[] library, boolean useModel, String[] models) {
		super(metaPar, magDbItem, librarySearch, selfModelSearch, library, useModel, models);
		this.setWorkflowType(uk.ac.quadram.metapilot.workflow.WorkflowType.DIA_MARKER.getAnalysisType());
		if (metaPar instanceof MetaParameterDia) {
			MetaParameterDia dia = (MetaParameterDia) metaPar;
			this.mExcision = dia.isMExcision();
			this.carbamiC = dia.isCarbamiC();
			this.oxidationM = dia.isOxidationM();
			this.acetyProN = dia.isAcetyProN();
			this.phosphoSTY = dia.isPhosphoSTY();
			this.kGG = dia.isKGG();
			this.setUsePeptideBasedRound(dia.isUsePeptideBasedRound());
			this.setPeptideDetectabilityThreshold(dia.getPeptideDetectabilityThreshold());
		}
	}

	// Getters and setters for modifications
	public boolean isMExcision() {
		return mExcision;
	}

	public void setMExcision(boolean mExcision) {
		this.mExcision = mExcision;
	}

	public boolean isCarbamiC() {
		return carbamiC;
	}

	public void setCarbamiC(boolean carbamiC) {
		this.carbamiC = carbamiC;
	}

	public boolean isOxidationM() {
		return oxidationM;
	}

	public void setOxidationM(boolean oxidationM) {
		this.oxidationM = oxidationM;
	}

	public boolean isAcetyProN() {
		return acetyProN;
	}

	public void setAcetyProN(boolean acetyProN) {
		this.acetyProN = acetyProN;
	}

	public boolean isPhosphoSTY() {
		return phosphoSTY;
	}

	public void setPhosphoSTY(boolean phosphoSTY) {
		this.phosphoSTY = phosphoSTY;
	}

	public boolean isKGG() {
		return kGG;
	}

	public void setKGG(boolean kGG) {
		this.kGG = kGG;
	}

	public String getMarkerMode() {
		return markerMode;
	}
	
	public void setMarkerMode(String markerMode) {
		this.markerMode = markerMode;
	}

	public int getGenomesPerRound() {
		return genomesPerRound;
	}

	public void setGenomesPerRound(int genomesPerRound) {
		this.genomesPerRound = genomesPerRound;
	}
}

