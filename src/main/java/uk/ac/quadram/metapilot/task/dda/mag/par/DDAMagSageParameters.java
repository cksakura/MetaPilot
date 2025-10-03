package uk.ac.quadram.metapilot.task.dda.mag.par;

import java.util.ArrayList;
import java.util.List;

import javafx.scene.control.Label;
import uk.ac.quadram.metapilot.dbSearch.MaxQuant.MaxquantModification;
import uk.ac.quadram.metapilot.task.common.gui.MetaPilotIOPanelFX.ExperimentRow;
import uk.ac.quadram.metapilot.task.common.gui.IsobaricLabelingPanel.IsobaricLabelRow;
import uk.ac.quadram.metapilot.task.dda.igc.gui.MetaPilotDdaIgcMaxQuantPanel.MetadataRow;
import uk.ac.quadram.metapilot.task.io.par.MetaOpenPar;
import uk.ac.quadram.metapilot.task.par.MetaBaseParameter;
import uk.ac.quadram.metapilot.task.par.MetaConstants;
import uk.ac.quadram.metapilot.task.par.MetaData;
import uk.ac.quadram.metapilot.task.par.MetaPep2TaxaPar;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;

/**
 * Parameter container for DDA MAG Sage workflow
 * 
 * This class holds all configuration parameters for a DDA MAG Sage analysis,
 * enabling parameter persistence, reproducibility, and UI state management.
 * Extends MetaBaseParameter to provide unified parameter architecture.
 */
public class DDAMagSageParameters extends MetaBaseParameter {

	// Input/Output settings
	private List<ExperimentRow> experimentData = new ArrayList<>();

	// Additional Sage-specific parameters
	private MaxquantModification[] isobaric = new MaxquantModification[0];
	private boolean combineLabel = false;
	private double[][] isoCorFactor = new double[0][0];

	// Add direct MagDbItem storage
	private MagDbItem magDbItem;

	// Status label for UI feedback
	private Label statusLabel;

	// Open search and taxonomic parameters
	private MetaOpenPar osp = MetaOpenPar.getDefault();
	private MetaPep2TaxaPar ptp = MetaPep2TaxaPar.getDefault();

	// Isobaric labeling
	private String isobaricType = "";
	private List<IsobaricLabelRow> isobaricTableData = new ArrayList<>();

	// Metadata
	private List<MetadataRow> metadataTableData = new ArrayList<>();

	// MAG-specific fields
	private String magDbCatalogueId = "";
	private String magDbVersion = "";

	// Isotopic labeling
	private List<String> lightLabels = new ArrayList<>();
	private List<String> mediumLabels = new ArrayList<>();
	private List<String> heavyLabels = new ArrayList<>();

	// Sage-specific parameters
	private double precursorToleranceMin = -500.0; // da (for open search)
	private double precursorToleranceMax = 150.0; // da (for open search)
	private double fragmentTolerance = 20.0; // ppm
	private String precursorToleranceUnit = "ppm";
	private String fragmentToleranceUnit = "ppm";
	private int minPeptideLength = 7;
	private int maxPeptideLength = 50;
	private double minPeakIntensity = 0.0;
	private int maxVariableModifications = 2;
	private boolean decoysEnabled = true;
	private boolean generateDecoys = true;
	private String decoysTag = "rev_";
	private double maxQ = 0.01; // FDR threshold
	private boolean writeCheckpoint = false;
	private String outputDirectory = "";

	// Additional Sage parameters from real config file
	private double peptideMinMass = 500.0;
	private double peptideMaxMass = 5000.0;
	private String[] ionKinds = { "b", "y" };
	private int minIonIndex = 2;
	private int minPeaks = 15;
	private int maxPeaks = 150;
	private int maxFragmentCharge = 1;
	private int minMatchedPeaks = 4;
	private boolean predictRt = true;
	private int[] isotopeErrors = { 0, 2 };
	private boolean deisotope = true;

	// LFQ-specific parameters
	private String peakScoring = "Hybrid";
	private String integration = "Sum";
	private double spectralAngle = 0.6;
	private double lfqPpmTolerance = 5.0;

	// Default constructor with sensible defaults for MAG Sage
	public DDAMagSageParameters() {
		super(uk.ac.quadram.metapilot.MetaPilotApp.WorkflowType.DDA_MAG_SAGE.getAnalysisType(), new MetaData(), "", "",
				false, "", MetaConstants.HCD_FTMS, 1, 1, true);

		// Initialize Sage-specific defaults
		setEnzyme("Trypsin");
		setMissCleavages(2);
		setDigestMode("Specific");
		setThreadCount(Runtime.getRuntime().availableProcessors());
		setQuantitativeMode(MetaConstants.labelFree);
		setReduceDb(true);
		setSpectraCluster(false);

		// Initialize with default modifications for Sage
		setFixMods(new String[] { "Carbamidomethyl[C] 0" });
		setVariMods(new String[] { "Acetyl[ProteinN-term] 0", "Oxidation[M] 0" });

		// Initialize Sage-specific defaults
		precursorToleranceMin = -20.0;
		precursorToleranceMax = 20.0;
		fragmentTolerance = 20.0;
		precursorToleranceUnit = "ppm";
		fragmentToleranceUnit = "ppm";
		minPeptideLength = 7;
		maxPeptideLength = 50;
		minPeakIntensity = 0.0;
		maxVariableModifications = 2;
		decoysEnabled = true;
		generateDecoys = true;
		decoysTag = "rev_";
		maxQ = 0.01;
		writeCheckpoint = false;

		// Initialize additional Sage parameters with defaults from real config
		peptideMinMass = 500.0;
		peptideMaxMass = 5000.0;
		ionKinds = new String[] { "b", "y" };
		minIonIndex = 2;
		minPeaks = 15;
		maxPeaks = 150;
		maxFragmentCharge = 1;
		minMatchedPeaks = 4;
		predictRt = true;
		isotopeErrors = new int[] { 0, 2 };
		deisotope = true;

		// Initialize LFQ parameters
		peakScoring = "Hybrid";
		integration = "Sum";
		spectralAngle = 0.6;
		lfqPpmTolerance = 5.0;
	}

	// Constructor for creating from GUI parameters
	public DDAMagSageParameters(String analysisType, MetaData metadata, String result, String microDb,
			boolean appendHostDb, String hostDB, String ms2ScanMode, int threadCount) {
		super(analysisType, metadata, result, microDb, appendHostDb, hostDB, ms2ScanMode, threadCount, threadCount,
				true);

		// Set Sage-specific defaults
		setEnzyme("Trypsin");
		setMissCleavages(2);
		setDigestMode("Specific");
		setQuantitativeMode(MetaConstants.labelFree);
		setReduceDb(true);
		setSpectraCluster(false);

		// Initialize Sage-specific defaults
		precursorToleranceMin = -20.0;
		precursorToleranceMax = 20.0;
		fragmentTolerance = 20.0;
		precursorToleranceUnit = "ppm";
		fragmentToleranceUnit = "ppm";
		minPeptideLength = 7;
		maxPeptideLength = 50;
		minPeakIntensity = 0.0;
		maxVariableModifications = 2;
		decoysEnabled = true;
		generateDecoys = true;
		decoysTag = "rev_";
		maxQ = 0.01;
		writeCheckpoint = false;

		// Initialize additional Sage parameters with defaults from real config
		peptideMinMass = 500.0;
		peptideMaxMass = 5000.0;
		ionKinds = new String[] { "b", "y" };
		minIonIndex = 2;
		minPeaks = 15;
		maxPeaks = 150;
		maxFragmentCharge = 1;
		minMatchedPeaks = 4;
		predictRt = true;
		isotopeErrors = new int[] { 0, 2 };
		deisotope = true;

		// Initialize LFQ parameters
		peakScoring = "Hybrid";
		integration = "Sum";
		spectralAngle = 0.6;
		lfqPpmTolerance = 5.0;
	}

	// Getters and setters for inherited functionality
	public List<ExperimentRow> getExperimentData() {
		return experimentData;
	}

	public void setExperimentData(List<ExperimentRow> experimentData) {
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

	public List<IsobaricLabelRow> getIsobaricTableData() {
		return isobaricTableData;
	}

	public void setIsobaricTableData(List<IsobaricLabelRow> isobaricTableData) {
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

	// Status label getter and setter
	public Label getStatusLabel() {
		return statusLabel;
	}

	public void setStatusLabel(Label statusLabel) {
		this.statusLabel = statusLabel;
	}

	// Metadata getters and setters
	public List<MetadataRow> getMetadataTableData() {
		return metadataTableData;
	}

	public void setMetadataTableData(List<MetadataRow> metadataTableData) {
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

	// Sage-specific parameter getters and setters

	public double getPrecursorToleranceMin() {
		return precursorToleranceMin;
	}

	public void setPrecursorToleranceMin(double precursorToleranceMin) {
		this.precursorToleranceMin = precursorToleranceMin;
	}

	public double getPrecursorToleranceMax() {
		return precursorToleranceMax;
	}

	public void setPrecursorToleranceMax(double precursorToleranceMax) {
		this.precursorToleranceMax = precursorToleranceMax;
	}

	public double getFragmentTolerance() {
		return fragmentTolerance;
	}

	public void setFragmentTolerance(double fragmentTolerance) {
		this.fragmentTolerance = fragmentTolerance;
	}

	public String getPrecursorToleranceUnit() {
		return precursorToleranceUnit;
	}

	public void setPrecursorToleranceUnit(String precursorToleranceUnit) {
		this.precursorToleranceUnit = precursorToleranceUnit;
	}

	public String getFragmentToleranceUnit() {
		return fragmentToleranceUnit;
	}

	public void setFragmentToleranceUnit(String fragmentToleranceUnit) {
		this.fragmentToleranceUnit = fragmentToleranceUnit;
	}

	public int getMinPeptideLength() {
		return minPeptideLength;
	}

	public void setMinPeptideLength(int minPeptideLength) {
		this.minPeptideLength = minPeptideLength;
	}

	public int getMaxPeptideLength() {
		return maxPeptideLength;
	}

	public void setMaxPeptideLength(int maxPeptideLength) {
		this.maxPeptideLength = maxPeptideLength;
	}

	public double getMinPeakIntensity() {
		return minPeakIntensity;
	}

	public void setMinPeakIntensity(double minPeakIntensity) {
		this.minPeakIntensity = minPeakIntensity;
	}

	public int getMaxVariableModifications() {
		return maxVariableModifications;
	}

	public void setMaxVariableModifications(int maxVariableModifications) {
		this.maxVariableModifications = maxVariableModifications;
	}

	public boolean isDecoysEnabled() {
		return decoysEnabled;
	}

	public void setDecoysEnabled(boolean decoysEnabled) {
		this.decoysEnabled = decoysEnabled;
	}

	public boolean isGenerateDecoys() {
		return generateDecoys;
	}

	public void setGenerateDecoys(boolean generateDecoys) {
		this.generateDecoys = generateDecoys;
	}

	public String getDecoysTag() {
		return decoysTag;
	}

	public void setDecoysTag(String decoysTag) {
		this.decoysTag = decoysTag;
	}

	public double getMaxQ() {
		return maxQ;
	}

	public void setMaxQ(double maxQ) {
		this.maxQ = maxQ;
	}

	public boolean isWriteCheckpoint() {
		return writeCheckpoint;
	}

	public void setWriteCheckpoint(boolean writeCheckpoint) {
		this.writeCheckpoint = writeCheckpoint;
	}

	public String getOutputDirectory() {
		return outputDirectory;
	}

	public void setOutputDirectory(String outputDirectory) {
		this.outputDirectory = outputDirectory;
	}

	// Additional Sage parameter getters and setters
	public double getPeptideMinMass() {
		return peptideMinMass;
	}

	public void setPeptideMinMass(double peptideMinMass) {
		this.peptideMinMass = peptideMinMass;
	}

	public double getPeptideMaxMass() {
		return peptideMaxMass;
	}

	public void setPeptideMaxMass(double peptideMaxMass) {
		this.peptideMaxMass = peptideMaxMass;
	}

	public String[] getIonKinds() {
		return ionKinds;
	}

	public void setIonKinds(String[] ionKinds) {
		this.ionKinds = ionKinds;
	}

	public int getMinIonIndex() {
		return minIonIndex;
	}

	public void setMinIonIndex(int minIonIndex) {
		this.minIonIndex = minIonIndex;
	}

	public int getMinPeaks() {
		return minPeaks;
	}

	public void setMinPeaks(int minPeaks) {
		this.minPeaks = minPeaks;
	}

	public int getMaxPeaks() {
		return maxPeaks;
	}

	public void setMaxPeaks(int maxPeaks) {
		this.maxPeaks = maxPeaks;
	}

	public int getMaxFragmentCharge() {
		return maxFragmentCharge;
	}

	public void setMaxFragmentCharge(int maxFragmentCharge) {
		this.maxFragmentCharge = maxFragmentCharge;
	}

	public int getMinMatchedPeaks() {
		return minMatchedPeaks;
	}

	public void setMinMatchedPeaks(int minMatchedPeaks) {
		this.minMatchedPeaks = minMatchedPeaks;
	}

	public boolean isPredictRt() {
		return predictRt;
	}

	public void setPredictRt(boolean predictRt) {
		this.predictRt = predictRt;
	}

	public int[] getIsotopeErrors() {
		return isotopeErrors;
	}

	public void setIsotopeErrors(int[] isotopeErrors) {
		this.isotopeErrors = isotopeErrors;
	}

	public boolean isDeisotope() {
		return deisotope;
	}

	public void setDeisotope(boolean deisotope) {
		this.deisotope = deisotope;
	}

	// LFQ parameter getters and setters
	public String getPeakScoring() {
		return peakScoring;
	}

	public void setPeakScoring(String peakScoring) {
		this.peakScoring = peakScoring;
	}

	public String getIntegration() {
		return integration;
	}

	public void setIntegration(String integration) {
		this.integration = integration;
	}

	public double getSpectralAngle() {
		return spectralAngle;
	}

	public void setSpectralAngle(double spectralAngle) {
		this.spectralAngle = spectralAngle;
	}

	public double getLfqPpmTolerance() {
		return lfqPpmTolerance;
	}

	public void setLfqPpmTolerance(double lfqPpmTolerance) {
		this.lfqPpmTolerance = lfqPpmTolerance;
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
