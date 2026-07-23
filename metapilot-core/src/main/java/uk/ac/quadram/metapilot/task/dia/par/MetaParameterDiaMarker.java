package uk.ac.quadram.metapilot.task.dia.par;

import uk.ac.quadram.metapilot.task.par.MetaBaseParameter;
import uk.ac.quadram.metapilot.task.dia.DiaRefinementMode;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;

/**
 * Parameter class for MAG-based DIA marker/seed workflows.
 *
 * Shared by:
 * - DIA Marker workflows (Core / Focus, Comprehensive, and Deep variants)
 * - DIA Seed iterative workflows (Core / Focus and Comprehensive)
 * - Seed Deep Together DIA workflow
 *
 * Refinement mode is encoded via {@link #markerMode}.
 */
public class MetaParameterDiaMarker extends MetaParameterDia {

	/**
	 * Marker workflow mode persisted for UI/CLI round-trip.
	 * <p>
	 * Known values:
	 * - FOCUS (legacy persisted value for Core)
	 * - COMPREHENSIVE
	 * - DEEP (legacy persisted value for Deep Together)
	 * - DEEP_PER_FILE (iterative workflow run independently per DIA file; results under marker/<rawname>/...)
	 */
	private String markerMode = "FOCUS";

	/**
	 * MAG registry identity fields used when {@link #getMagDbItem()} is null
	 * (e.g. HPC resolve against a site MAGDB JSON).
	 */
	private String magDbIdentifier = "MGYG";
	private String magDbCatalogueVersion = "";

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

	public static final String REPRESENTATIVE_CLUSTERING_OFF = "OFF";
	public static final String REPRESENTATIVE_CLUSTERING_CONSERVATIVE = "CONSERVATIVE";
	public static final String REPRESENTATIVE_CLUSTERING_BALANCED = "BALANCED";
	public static final String REPRESENTATIVE_CLUSTERING_AGGRESSIVE = "AGGRESSIVE";
	public static final String REPRESENTATIVE_CLUSTERING_CUSTOM = "CUSTOM";

	/**
	 * Optional DEEP_PER_FILE speed-up for large cohorts. OFF preserves the original
	 * behavior exactly.
	 */
	private String representativeClusteringPreset = REPRESENTATIVE_CLUSTERING_OFF;
	private double representativePeptideJaccardThreshold = 0.62;
	private double representativePeptideContainmentThreshold = 0.79;
	private int representativeGenomeTopN = 50;
	private int representativeMinMarkerPeptides = 100;
	private boolean representativeClusteringDryRun = false;

	public MetaParameterDiaMarker() {
		super();
		this.setWorkflowType(uk.ac.quadram.metapilot.workflow.WorkflowType.DIA_MARKER.getAnalysisType());
	}

	public MetaParameterDiaMarker(MetaBaseParameter metaPar, MagDbItem magDbItem, boolean librarySearch,
			boolean selfModelSearch, String[] library, boolean useModel, String[] models) {
		super(metaPar, magDbItem, librarySearch, selfModelSearch, library, useModel, models);
		this.setWorkflowType(uk.ac.quadram.metapilot.workflow.WorkflowType.DIA_MARKER.getAnalysisType());
		if (magDbItem != null && magDbItem.getCatalogueID() != null && !magDbItem.getCatalogueID().isBlank()) {
			this.setMicroDb(magDbItem.getCatalogueID());
		}
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

	@Override
	public void setMagDbItem(MagDbItem magDbItem) {
		super.setMagDbItem(magDbItem);
		if (magDbItem != null && magDbItem.getCatalogueID() != null && !magDbItem.getCatalogueID().isBlank()) {
			this.setMicroDb(magDbItem.getCatalogueID());
		}
		if (magDbItem != null) {
			if (magDbItem.getIdentifier() != null && !magDbItem.getIdentifier().isBlank()) {
				this.magDbIdentifier = magDbItem.getIdentifier();
			}
			if (magDbItem.getUsedVersion() != null
					&& magDbItem.getUsedVersion().getCatalogueVersion() != null
					&& !magDbItem.getUsedVersion().getCatalogueVersion().isBlank()) {
				this.magDbCatalogueVersion = magDbItem.getUsedVersion().getCatalogueVersion();
			}
		}
	}

	public String getMagDbIdentifier() {
		return magDbIdentifier == null || magDbIdentifier.isBlank() ? "MGYG" : magDbIdentifier;
	}

	public void setMagDbIdentifier(String magDbIdentifier) {
		this.magDbIdentifier = magDbIdentifier == null || magDbIdentifier.isBlank() ? "MGYG" : magDbIdentifier.trim();
	}

	public String getMagDbCatalogueVersion() {
		return magDbCatalogueVersion == null ? "" : magDbCatalogueVersion;
	}

	public void setMagDbCatalogueVersion(String magDbCatalogueVersion) {
		this.magDbCatalogueVersion = magDbCatalogueVersion == null ? "" : magDbCatalogueVersion.trim();
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
		this.markerMode = DiaRefinementMode.fromMarkerMode(markerMode).getPersistedValue();
	}

	public DiaRefinementMode getMarkerRefinementMode() {
		return DiaRefinementMode.fromMarkerMode(markerMode);
	}

	public void setMarkerRefinementMode(DiaRefinementMode mode) {
		this.markerMode = (mode != null ? mode : DiaRefinementMode.CORE).getPersistedValue();
	}

	public int getGenomesPerRound() {
		return genomesPerRound;
	}

	public void setGenomesPerRound(int genomesPerRound) {
		this.genomesPerRound = genomesPerRound;
	}

	public String getRepresentativeClusteringPreset() {
		return representativeClusteringPreset;
	}

	public void setRepresentativeClusteringPreset(String preset) {
		String value = preset == null ? REPRESENTATIVE_CLUSTERING_OFF : preset.trim().toUpperCase();
		switch (value) {
		case REPRESENTATIVE_CLUSTERING_CONSERVATIVE:
			this.representativeClusteringPreset = REPRESENTATIVE_CLUSTERING_CONSERVATIVE;
			this.representativePeptideJaccardThreshold = 0.62;
			this.representativePeptideContainmentThreshold = 0.79;
			break;
		case REPRESENTATIVE_CLUSTERING_BALANCED:
			this.representativeClusteringPreset = REPRESENTATIVE_CLUSTERING_BALANCED;
			this.representativePeptideJaccardThreshold = 0.60;
			this.representativePeptideContainmentThreshold = 0.78;
			break;
		case REPRESENTATIVE_CLUSTERING_AGGRESSIVE:
			this.representativeClusteringPreset = REPRESENTATIVE_CLUSTERING_AGGRESSIVE;
			this.representativePeptideJaccardThreshold = 0.58;
			this.representativePeptideContainmentThreshold = 0.75;
			break;
		case REPRESENTATIVE_CLUSTERING_CUSTOM:
			this.representativeClusteringPreset = REPRESENTATIVE_CLUSTERING_CUSTOM;
			break;
		default:
			this.representativeClusteringPreset = REPRESENTATIVE_CLUSTERING_OFF;
			break;
		}
	}

	public boolean isRepresentativeClusteringEnabled() {
		return !REPRESENTATIVE_CLUSTERING_OFF.equalsIgnoreCase(representativeClusteringPreset);
	}

	public double getRepresentativePeptideJaccardThreshold() {
		return representativePeptideJaccardThreshold;
	}

	public void setRepresentativePeptideJaccardThreshold(double threshold) {
		if (threshold > 0.0 && threshold <= 1.0) {
			this.representativePeptideJaccardThreshold = threshold;
		}
	}

	public double getRepresentativePeptideContainmentThreshold() {
		return representativePeptideContainmentThreshold;
	}

	public void setRepresentativePeptideContainmentThreshold(double threshold) {
		if (threshold > 0.0 && threshold <= 1.0) {
			this.representativePeptideContainmentThreshold = threshold;
		}
	}

	public int getRepresentativeGenomeTopN() {
		return representativeGenomeTopN;
	}

	public void setRepresentativeGenomeTopN(int genomeTopN) {
		if (genomeTopN > 0) {
			this.representativeGenomeTopN = genomeTopN;
		}
	}

	public int getRepresentativeMinMarkerPeptides() {
		return representativeMinMarkerPeptides;
	}

	public void setRepresentativeMinMarkerPeptides(int minMarkerPeptides) {
		if (minMarkerPeptides >= 0) {
			this.representativeMinMarkerPeptides = minMarkerPeptides;
		}
	}

	public boolean isRepresentativeClusteringDryRun() {
		return representativeClusteringDryRun;
	}

	public void setRepresentativeClusteringDryRun(boolean dryRun) {
		this.representativeClusteringDryRun = dryRun;
	}
}

