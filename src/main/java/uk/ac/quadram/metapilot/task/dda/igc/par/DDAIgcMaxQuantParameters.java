package uk.ac.quadram.metapilot.task.dda.igc.par;

import java.util.ArrayList;
import java.util.List;

import uk.ac.quadram.metapilot.task.common.gui.MetaPilotIOPanelFX.ExperimentRow;
import uk.ac.quadram.metapilot.task.dda.igc.gui.MetaPilotDdaIgcMaxQuantPanel.MetadataRow;
import uk.ac.quadram.metapilot.task.par.MetaBaseParameter;
import uk.ac.quadram.metapilot.task.par.MetaConstants;
import uk.ac.quadram.metapilot.task.par.MetaData;
import uk.ac.quadram.metapilot.task.dda.igc.gui.MetaPilotDdaIgcMaxQuantPanel.IsobaricLabelRow;

/**
 * Parameter container for DDA IGC MaxQuant workflow
 * 
 * This class holds all configuration parameters for a DDA IGC MaxQuant
 * analysis, enabling parameter persistence, reproducibility, and UI state
 * management.
 */
public class DDAIgcMaxQuantParameters extends MetaBaseParameter {

	// Input/Output settings
	private String outputDirectory = "";
	private List<ExperimentRow> experimentData = new ArrayList<>();

	// Isotopic labeling
	private List<String> lightLabels = new ArrayList<>();
	private List<String> mediumLabels = new ArrayList<>();
	private List<String> heavyLabels = new ArrayList<>();

	// Isobaric labeling
	private String isobaricType = "";
	private List<IsobaricLabelRow> isobaricTableData = new ArrayList<>();

	// Taxonomic analysis
	private boolean eukaryota = false;
	private boolean bacteria = false;
	private boolean archaea = false;
	private boolean viruses = false;
	private boolean viroids = false;
	private boolean unclassified = false;
	private boolean other = false;
	private String ignoreRank = "Species";
	private boolean ignoreEnvironmentalSamples = false;

	// Metadata
	private List<MetadataRow> metadataTableData = new ArrayList<>();

	// Default constructor with sensible defaults
	public DDAIgcMaxQuantParameters() {
		super("DDA_IGC_MAXQUANT", new MetaData(), "", "", false, "", MetaConstants.FTMS, 1, 1, true);

		// Initialize with default values
		setMs2ScanMode(MetaConstants.HCD_FTMS);
		setThreadCount(Runtime.getRuntime().availableProcessors());
		setSpectraCluster(false);
	}

	public DDAIgcMaxQuantParameters(String analysisType, MetaData metadata, String result, String microDb,
			boolean appendHostDb, String hostDB, String ms2ScanMode, int threadCount) {
		super(analysisType, metadata, result, microDb, appendHostDb, hostDB, ms2ScanMode, threadCount, threadCount,
				true);
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

	public boolean isEukaryota() {
		return eukaryota;
	}

	public void setEukaryota(boolean eukaryota) {
		this.eukaryota = eukaryota;
	}

	public boolean isBacteria() {
		return bacteria;
	}

	public void setBacteria(boolean bacteria) {
		this.bacteria = bacteria;
	}

	public boolean isArchaea() {
		return archaea;
	}

	public void setArchaea(boolean archaea) {
		this.archaea = archaea;
	}

	public boolean isViruses() {
		return viruses;
	}

	public void setViruses(boolean viruses) {
		this.viruses = viruses;
	}

	public boolean isViroids() {
		return viroids;
	}

	public void setViroids(boolean viroids) {
		this.viroids = viroids;
	}

	public boolean isUnclassified() {
		return unclassified;
	}

	public void setUnclassified(boolean unclassified) {
		this.unclassified = unclassified;
	}

	public boolean isOther() {
		return other;
	}

	public void setOther(boolean other) {
		this.other = other;
	}

	public String getIgnoreRank() {
		return ignoreRank;
	}

	public void setIgnoreRank(String ignoreRank) {
		this.ignoreRank = ignoreRank;
	}

	public boolean isIgnoreEnvironmentalSamples() {
		return ignoreEnvironmentalSamples;
	}

	public void setIgnoreEnvironmentalSamples(boolean ignoreEnvironmentalSamples) {
		this.ignoreEnvironmentalSamples = ignoreEnvironmentalSamples;
	}

	public List<MetadataRow> getMetadataTableData() {
		return metadataTableData;
	}

	public void setMetadataTableData(List<MetadataRow> metadataTableData) {
		this.metadataTableData = metadataTableData;
	}

	/**
	 * Create a copy of this parameters object
	 */
	public DDAIgcMaxQuantParameters copy() {
		DDAIgcMaxQuantParameters copy = new DDAIgcMaxQuantParameters();

		// Copy all fields
		copy.outputDirectory = this.outputDirectory;
		copy.experimentData = new ArrayList<>(this.experimentData);
		copy.setReduceDb(this.isReduceDb());

		// Copy inherited fields from MetaBaseParameter
		copy.setMicroDb(this.getMicroDb());
		copy.setAppendHostDb(this.isAppendHostDb());
		copy.setHostDB(this.getHostDB());
		copy.setEnzyme(this.getEnzyme());
		copy.setMissCleavages(this.getMissCleavages());
		copy.setDigestMode(this.getDigestModeString());
		copy.setThreadCount(this.getThreadCount());
		copy.setMs2ScanMode(this.getMs2ScanMode());
		copy.setSpectraCluster(this.isSpectraCluster());
		copy.setFixMods(this.getFixMods());
		copy.setVariMods(this.getVariMods());

		copy.setQuantitativeMode(this.getQuantitativeMode());
		copy.lightLabels = new ArrayList<>(this.lightLabels);
		copy.mediumLabels = new ArrayList<>(this.mediumLabels);
		copy.heavyLabels = new ArrayList<>(this.heavyLabels);

		copy.isobaricType = this.isobaricType;
		copy.isobaricTableData = new ArrayList<>(this.isobaricTableData);

		copy.eukaryota = this.eukaryota;
		copy.bacteria = this.bacteria;
		copy.archaea = this.archaea;
		copy.viruses = this.viruses;
		copy.viroids = this.viroids;
		copy.unclassified = this.unclassified;
		copy.other = this.other;
		copy.ignoreRank = this.ignoreRank;
		copy.ignoreEnvironmentalSamples = this.ignoreEnvironmentalSamples;

		copy.metadataTableData = new ArrayList<>(this.metadataTableData);

		return copy;
	}

	public String toString() {
		return String.format(
				"DDAIgcMaxQuantParameters[enzyme=%s, missed_cleavages=%d, quant_mode=%s, " + "output=%s, threads=%d]",
				getEnzyme(), getMissCleavages(), getQuantitativeMode(), outputDirectory, getThreadCount());
	}
}
