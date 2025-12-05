package uk.ac.quadram.metapilot.task.dda.igc.par;

import java.util.ArrayList;
import java.util.List;

import uk.ac.quadram.metapilot.dbSearch.MaxQuant.MaxquantModification;
import uk.ac.quadram.metapilot.task.par.MetaBaseParameter;
import uk.ac.quadram.metapilot.task.par.MetaConstants;
import uk.ac.quadram.metapilot.task.par.MetaData;

/**
 * Parameter container for DDA IGC Sage workflow
 * 
 * This class holds all configuration parameters for a DDA IGC Sage analysis,
 * enabling parameter persistence, reproducibility, and UI state management.
 * Extends MetaBaseParameter to provide unified parameter architecture.
 */
public class DDAIgcSageParameters extends MetaBaseParameter {

	// Input/Output settings - ExperimentData format: [filePath, experimentName, fraction]
	private List<String[]> experimentData = new ArrayList<>();

	// Additional Sage-specific parameters
	private MaxquantModification[] isobaric = new MaxquantModification[0];
	private boolean combineLabel = false;
	private double[][] isoCorFactor = new double[0][0];

	// Isobaric labeling - IsobaricLabelData format: [labelName, minus2, minus1, plus1, plus2]
	private String isobaricType = "";
	private List<String[]> isobaricTableData = new ArrayList<>();

	// Metadata - MetadataData format: [sampleName, labelType, age, gender, condition, treatment, timePoint, customField1, customField2, customField3]
	private List<String[]> metadataTableData = new ArrayList<>();

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

	// LFQ-specific parameters
	private String peakScoring = "Hybrid";
	private String integration = "Sum";
	private double spectralAngle = 0.6;
	private double lfqPpmTolerance = 5.0;
	
	// Default constructor with sensible defaults for IGC Sage
	public DDAIgcSageParameters() {
		super("DDA_IGC_SAGE", new MetaData(), "", "",
				false, "", "HCD-FTMS", Runtime.getRuntime().availableProcessors(),
				Runtime.getRuntime().availableProcessors(), true);

		// Initialize Sage-specific defaults
		setEnzyme("Trypsin");
		setMissCleavages(2);
		setDigestMode("Specific");
		setMs2ScanMode(MetaConstants.HCD_FTMS);
		setThreadCount(Runtime.getRuntime().availableProcessors());
		setQuantitativeMode(MetaConstants.labelFree);
		setReduceDb(true);
		setSpectraCluster(false);

		// Initialize with default modifications for Sage
		setFixMods(new String[] { "Carbamidomethyl[C] 0" });
		setVariMods(new String[] { "Acetyl[ProteinN-term] 0", "Oxidation[M] 0" });
	}

	// Constructor for creating from GUI parameters
	public DDAIgcSageParameters(String analysisType, MetaData metadata, String result, String microDb,
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

	// Sage-specific methods
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

	public String getQuanMode() {
		return getQuantitativeMode();
	}

	public void setQuanMode(String quanMode) {
		setQuantitativeMode(quanMode);
	}

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

	public List<String[]> getMetadataTableData() {
		return metadataTableData;
	}

	public void setMetadataTableData(List<String[]> metadataTableData) {
		this.metadataTableData = metadataTableData;
	}

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

	public boolean isWriteCheckpoint() {
		return writeCheckpoint;
	}

	public void setWriteCheckpoint(boolean writeCheckpoint) {
		this.writeCheckpoint = writeCheckpoint;
	}

	/**
	 * Create a copy of this parameters object
	 */
	public DDAIgcSageParameters copy() {
		DDAIgcSageParameters copy = new DDAIgcSageParameters();

		// Copy inherited fields from MetaBaseParameter
		copy.setResult(this.getResult());
		copy.setMicroDb(this.getMicroDb());
		copy.setAppendHostDb(this.isAppendHostDb());
		copy.setHostDB(this.getHostDB());
		copy.setEnzyme(this.getEnzyme());
		copy.setMissCleavages(this.getMissCleavages());
		copy.setDigestMode(this.getDigestModeString());
		copy.setThreadCount(this.getThreadCount());
		copy.setMs2ScanMode(this.getMs2ScanMode());
		copy.setReduceDb(this.isReduceDb());
		copy.setSpectraCluster(this.isSpectraCluster());
		copy.setFixMods(this.getFixMods());
		copy.setVariMods(this.getVariMods());
		copy.setQuantitativeMode(this.getQuantitativeMode());
		copy.setIsobaricTag(this.getIsobaricTag());
		copy.setIsobaricRefTable(this.getIsobaricRefTable());
		copy.setIsotopicLabels(this.getIsotopicLabels());

		// Copy Sage-specific fields
		copy.experimentData = new ArrayList<>(this.experimentData);
		copy.isobaric = this.isobaric.clone();
		copy.combineLabel = this.combineLabel;
		copy.isoCorFactor = this.isoCorFactor.clone();

		copy.isobaricType = this.isobaricType;
		copy.isobaricTableData = new ArrayList<>(this.isobaricTableData);

		copy.metadataTableData = new ArrayList<>(this.metadataTableData);

		copy.lightLabels = new ArrayList<>(this.lightLabels);
		copy.mediumLabels = new ArrayList<>(this.mediumLabels);
		copy.heavyLabels = new ArrayList<>(this.heavyLabels);

		// Copy Sage-specific parameters
		copy.precursorToleranceMin = this.precursorToleranceMin;
		copy.precursorToleranceMax = this.precursorToleranceMax;
		copy.fragmentTolerance = this.fragmentTolerance;
		copy.precursorToleranceUnit = this.precursorToleranceUnit;
		copy.fragmentToleranceUnit = this.fragmentToleranceUnit;
		copy.minPeptideLength = this.minPeptideLength;
		copy.maxPeptideLength = this.maxPeptideLength;
		copy.minPeakIntensity = this.minPeakIntensity;
		copy.maxVariableModifications = this.maxVariableModifications;
		copy.decoysEnabled = this.decoysEnabled;
		copy.generateDecoys = this.generateDecoys;
		copy.decoysTag = this.decoysTag;
		copy.maxQ = this.maxQ;
		copy.writeCheckpoint = this.writeCheckpoint;

		return copy;
	}

	@Override
	public String toString() {
		return String.format(
				"DDAIgcSageParameters[enzyme=%s, missed_cleavages=%d, quant_mode=%s, "
						+ "output=%s, threads=%d, precursor_tolerance=%.1f-%.1f%s]",
				getEnzyme(), getMissCleavages(), getQuantitativeMode(), getResult(), getThreadCount(),
				precursorToleranceMin, precursorToleranceMax, precursorToleranceUnit);
	}
}
