package uk.ac.quadram.metapilot.task.dda.mag.par;

import java.util.ArrayList;
import java.util.List;

import uk.ac.quadram.metapilot.core.mod.IsobaricTag;
import uk.ac.quadram.metapilot.dbSearch.MaxQuant.MaxquantModification;
import uk.ac.quadram.metapilot.task.io.par.MetaOpenPar;
import uk.ac.quadram.metapilot.task.par.IsobaricLabelSettings;
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
	private String genomeSelectionStrategy = "Focus";

	// Isotopic labeling
	private List<String> lightLabels = new ArrayList<>();
	private List<String> mediumLabels = new ArrayList<>();
	private List<String> heavyLabels = new ArrayList<>();

	/** Pangenome variant check (UHGP-100); runs only for open search. */
	private final uk.ac.quadram.metapilot.core.variant.PangenomeVariantSettings pangenomeVariantSettings =
			new uk.ac.quadram.metapilot.core.variant.PangenomeVariantSettings();
	/**
	 * A pFind substitution within this many Da of a modification on the same
	 * residue is reported as MS2-undecidable. pFind's own isobaric aliases are
	 * undecidable regardless.
	 */
	private double undecidableToleranceDa = DEFAULT_UNDECIDABLE_TOLERANCE_DA;
	public static final double DEFAULT_UNDECIDABLE_TOLERANCE_DA = 0.005;

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
		// Initialize with default modifications and digest mode for pFind
		setFixMods(new String[] { "Carbamidomethyl[C] 0" });
		setVariMods(new String[] { "Acetyl[ProteinN-term] 0", "Oxidation[M] 0" });
		setDigestMode("Fully-Specific");
	}

	// Constructor for creating from GUI parameters
	public DDAMagPfindParameters(String analysisType, MetaData metadata, String result, String microDb,
			boolean appendHostDb, String hostDB, String ms2ScanMode, int threadCount) {
		super(analysisType, metadata, result, microDb, appendHostDb, hostDB, ms2ScanMode, threadCount, threadCount,
				true);

		// Set pFind-specific defaults
		setEnzyme("Trypsin_P");
		setMissCleavages(2);
		setDigestMode("Fully-Specific");
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

	/**
	 * Per-channel reporter labels for the resolved tag. An explicitly set array
	 * is used only when it matches the tag channel for channel; otherwise the
	 * labels are built from the tag.
	 */
	public MaxquantModification[] getIsobaric() {
		IsobaricTag tag = getIsobaricTag();
		if (tag == null) {
			return isobaric;
		}
		if (IsobaricLabelSettings.reporterLabelsFit(tag, isobaric)) {
			return isobaric;
		}
		return IsobaricLabelSettings.buildReporterLabels(tag);
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

	/**
	 * Isotope correction factors (nPlex x 4). An explicitly set array is used
	 * when its shape fits the tag; otherwise the values come from the isobaric
	 * table, with 0 (no correction) for channels that have no row.
	 */
	public double[][] getIsoCorFactor() {
		IsobaricTag tag = getIsobaricTag();
		if (tag == null) {
			return isoCorFactor;
		}
		if (IsobaricLabelSettings.correctionFactorsFit(tag, isoCorFactor)) {
			return isoCorFactor;
		}
		return IsobaricLabelSettings.buildCorrectionFactors(tag, isobaricTableData);
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

	public String getGenomeSelectionStrategy() {
		return genomeSelectionStrategy;
	}

	public void setGenomeSelectionStrategy(String genomeSelectionStrategy) {
		if (genomeSelectionStrategy == null || genomeSelectionStrategy.trim().isEmpty()) {
			this.genomeSelectionStrategy = "Focus";
		} else {
			this.genomeSelectionStrategy = genomeSelectionStrategy.trim();
		}
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

	/** Whether open search is selected (open=1 in pFind.cfg). When true, fixmod/selectmod are not written. */
	public boolean isOpenSearch() {
		return osp != null && osp.isOpenSearch();
	}

	public uk.ac.quadram.metapilot.core.variant.PangenomeVariantSettings getPangenomeVariantSettings() {
		return pangenomeVariantSettings;
	}

	public double getUndecidableToleranceDa() {
		return undecidableToleranceDa;
	}

	public void setUndecidableToleranceDa(double undecidableToleranceDa) {
		this.undecidableToleranceDa = undecidableToleranceDa;
	}

	/**
	 * The isobaric tag of this run, or null when the quantification mode is not
	 * isobaric. The GUI never sets the tag directly, so when none was set it is
	 * resolved from the selected type, the isobaric table, the MetaData label
	 * titles or the metadata table (see {@link IsobaricLabelSettings}).
	 */
	@Override
	public IsobaricTag getIsobaricTag() {
		return IsobaricLabelSettings.resolveTag(getQuantitativeMode(), super.getIsobaricTag(), isobaricType,
				isobaricTableData, getMetadata(), metadataTableData);
	}

	/**
	 * Maps the GUI radio text ("Isobaric", "Label-free", ...) onto the
	 * MetaConstants values DDAMagPfindTask compares against; with the raw text
	 * "Isobaric" neither the label-free nor the isobaric quantification branch
	 * ran.
	 */
	@Override
	public void setQuantitativeMode(String quantitativeMode) {
		super.setQuantitativeMode(IsobaricLabelSettings.normalizeQuantitativeMode(quantitativeMode));
	}

}
