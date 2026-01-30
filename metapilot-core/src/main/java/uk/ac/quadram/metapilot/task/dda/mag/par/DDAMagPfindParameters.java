package uk.ac.quadram.metapilot.task.dda.mag.par;

import java.util.ArrayList;
import java.util.List;

import uk.ac.quadram.metapilot.dbSearch.MaxQuant.MaxquantModification;
import uk.ac.quadram.metapilot.task.io.par.MetaOpenPar;
import uk.ac.quadram.metapilot.task.par.MetaBaseParameter;
import uk.ac.quadram.metapilot.task.par.MetaConstants;
import uk.ac.quadram.metapilot.task.par.MetaData;
import uk.ac.quadram.metapilot.task.par.MetaPep2TaxaPar;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;

/**
 * Parameter container for DDA MAG pFind workflow
 * 
 * This class holds all configuration parameters for a DDA MAG pFind analysis,
 * enabling parameter persistence, reproducibility, and UI state management.
 * Extends MetaBaseParameter to provide unified parameter architecture.
 */
public class DDAMagPfindParameters extends MetaBaseParameter {

	// Input/Output settings - now uses String[] arrays instead of JavaFX row
	// objects
	// ExperimentData format: [filePath, experimentName, fraction]
	private List<String[]> experimentData = new ArrayList<>();

	// Additional pFind-specific parameters
	private MaxquantModification[] isobaric = new MaxquantModification[0];
	private boolean combineLabel = false;
	private double[][] isoCorFactor = new double[0][0];

	// Add direct MagDbItem storage
	private MagDbItem magDbItem;

	// Open search and taxonomic parameters
	private MetaOpenPar osp = MetaOpenPar.getDefault();
	private MetaPep2TaxaPar ptp = MetaPep2TaxaPar.getDefault();

	// Isobaric labeling - now uses String[] arrays instead of JavaFX row objects
	// IsobaricLabelData format: [labelName, minus2, minus1, plus1, plus2]
	private String isobaricType = "";
	private List<String[]> isobaricTableData = new ArrayList<>();

	// Metadata - now uses String[] arrays instead of JavaFX row objects
	// MetadataData format: [sampleName, labelType, age, gender, condition,
	// treatment, timePoint, customField1, customField2, customField3]
	private List<String[]> metadataTableData = new ArrayList<>();

	// MAG-specific fields
	private String magDbCatalogueId = "";
	private String magDbVersion = "";
	private String genomeSelectorMode = MagDbItem.MARKER2;

	// Isotopic labeling
	private List<String> lightLabels = new ArrayList<>();
	private List<String> mediumLabels = new ArrayList<>();
	private List<String> heavyLabels = new ArrayList<>();

	// Default constructor with sensible defaults for MAG pFind
	public DDAMagPfindParameters() {
		super("DDA_MAG_PFIND", new MetaData(), "", "", false, "", MetaConstants.HCD_FTMS, 1, 1, true);

		// Initialize pFind-specific defaults
		setEnzyme("Trypsin_P");
		setMissCleavages(2);
		setDigestMode("Specific");
		setMs2ScanMode(MetaConstants.HCD_FTMS);
		setThreadCount(Runtime.getRuntime().availableProcessors());
		setQuantitativeMode(MetaConstants.labelFree);
		// Initialize with default modifications for pFind
		setFixMods(new String[] { "Carbamidomethyl[C] 0" });
		setVariMods(new String[] { "Acetyl[ProteinN-term] 0", "Oxidation[M] 0" });
	}

	// Constructor for creating from GUI parameters
	public DDAMagPfindParameters(String analysisType, MetaData metadata, String result, String microDb,
			boolean appendHostDb, String hostDB, String ms2ScanMode, int threadCount) {
		super(analysisType, metadata, result, microDb, appendHostDb, hostDB, ms2ScanMode, threadCount, threadCount,
				true);

		// Set pFind-specific defaults
		setEnzyme("Trypsin_P");
		setMissCleavages(2);
		setDigestMode("Specific");
		setQuantitativeMode(MetaConstants.labelFree);
		setReduceDb(true);
		setSpectraCluster(false);
	}

	// Getters and setters for inherited functionality
	public List<String[]> getExperimentData() {
		return experimentData;
	}

	public void setExperimentData(List<String[]> experimentData) {
		this.experimentData = experimentData;
	}

	// Use inherited methods for modifications - convert between String[] and
	// List<String>
	public List<String> getFixedModifications() {
		String[] fixMods = getFixMods();
		return fixMods != null ? List.of(fixMods) : new ArrayList<>();
	}

	public void setFixedModifications(List<String> fixedModifications) {
		setFixMods(fixedModifications.toArray(new String[0]));
	}

	public List<String> getVariableModifications() {
		String[] variMods = getVariMods();
		return variMods != null ? List.of(variMods) : new ArrayList<>();
	}

	public void setVariableModifications(List<String> variableModifications) {
		setVariMods(variableModifications.toArray(new String[0]));
	}

	// Isobaric labeling getters and setters
	public String getIsobaricType() {
		return isobaricType;
	}

	public void setIsobaricType(String isobaricType) {
		this.isobaricType = isobaricType;
	}

	public List<String[]> getIsobaricTableData() {
		return isobaricTableData;
	}

	public void setIsobaricTableData(List<String[]> isobaricTableData) {
		this.isobaricTableData = isobaricTableData;
	}

	public MaxquantModification[] getIsobaric() {
		return isobaric;
	}

	public void setIsobaric(MaxquantModification[] isobaric) {
		this.isobaric = isobaric;
	}

	public boolean isCombineLabel() {
		return combineLabel;
	}

	public void setCombineLabel(boolean combineLabel) {
		this.combineLabel = combineLabel;
	}

	public double[][] getIsoCorFactor() {
		return isoCorFactor;
	}

	public void setIsoCorFactor(double[][] isoCorFactor) {
		this.isoCorFactor = isoCorFactor;
	}

	// Add getter and setter for MagDbItem
	public MagDbItem getMagDbItem() {
		return magDbItem;
	}

	public void setMagDbItem(MagDbItem magDbItem) {
		this.magDbItem = magDbItem;
		this.setMicroDb(magDbItem.getCurrentFile().getAbsolutePath());
	}

	// Metadata getters and setters
	public List<String[]> getMetadataTableData() {
		return metadataTableData;
	}

	public void setMetadataTableData(List<String[]> metadataTableData) {
		this.metadataTableData = metadataTableData;
	}

	// MAG-specific getters and setters
	public String getMagDbCatalogueId() {
		return magDbCatalogueId;
	}

	public void setMagDbCatalogueId(String magDbCatalogueId) {
		this.magDbCatalogueId = magDbCatalogueId;
	}

	public String getMagDbVersion() {
		return magDbVersion;
	}

	public void setMagDbVersion(String magDbVersion) {
		this.magDbVersion = magDbVersion;
	}

	public String getGenomeSelectorMode() {
		return genomeSelectorMode;
	}

	public void setGenomeSelectorMode(String genomeSelectorMode) {
		this.genomeSelectorMode = genomeSelectorMode;
	}

	// Isotopic labeling getters and setters
	public List<String> getLightLabels() {
		return lightLabels;
	}

	public void setLightLabels(List<String> lightLabels) {
		this.lightLabels = lightLabels;
	}

	public List<String> getMediumLabels() {
		return mediumLabels;
	}

	public void setMediumLabels(List<String> mediumLabels) {
		this.mediumLabels = mediumLabels;
	}

	public List<String> getHeavyLabels() {
		return heavyLabels;
	}

	public void setHeavyLabels(List<String> heavyLabels) {
		this.heavyLabels = heavyLabels;
	}

	// Open search and taxonomic parameters
	public MetaOpenPar getOsp() {
		return osp;
	}

	public void setOsp(MetaOpenPar osp) {
		this.osp = osp;
	}

	// Alias for getOsp() to match the expected method name
	public MetaOpenPar getMop() {
		return osp;
	}

	public MetaPep2TaxaPar getPtp() {
		return ptp;
	}

	public void setPtp(MetaPep2TaxaPar ptp) {
		this.ptp = ptp;
	}

	// Open search parameter
	public boolean isOpenSearch() {
		return osp != null && osp.getOpenPsmFDR() > 0;
	}

	// Isobaric tag getter
	public uk.ac.quadram.metapilot.core.mod.IsobaricTag getIsobaricTag() {
		// Return null for now - this needs to be implemented based on the isobaric data
		return null;
	}

	/**
	 * Convert digestion mode string to integer representation
	 */
	private int getDigestModeInt(String digestMode) {
		switch (digestMode) {
		case "Specific":
			return 1;
		case "Semi free N-term":
			return 2;
		case "Semi free C-term":
			return 3;
		case "Semi specific":
			return 4;
		case "Unspecific":
			return 5;
		case "None":
			return 6;
		default:
			return 1;
		}
	}
}
