package uk.ac.quadram.metapilot.task.dda.mag.par;

import java.util.ArrayList;
import java.util.List;

import javafx.scene.control.Label;
import uk.ac.quadram.metapilot.dbSearch.MaxQuant.MaxquantModification;
import uk.ac.quadram.metapilot.task.common.gui.MetaPilotIOPanelFX.ExperimentRow;
import uk.ac.quadram.metapilot.task.common.gui.IsobaricLabelingPanel.IsobaricLabelRow;
import uk.ac.quadram.metapilot.task.dda.igc.gui.MetaPilotDdaIgcMaxQuantPanel.MetadataRow;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;
import uk.ac.quadram.metapilot.task.par.MetaBaseParameter;
import uk.ac.quadram.metapilot.task.par.MetaConstants;
import uk.ac.quadram.metapilot.task.par.MetaData;
import uk.ac.quadram.metapilot.dbSearch.fragpipe.FragpipeWorkflow;

/**
 * Parameter container for DDA MAG FragPipe workflow
 * 
 * This class holds all configuration parameters for a DDA MAG FragPipe analysis,
 * enabling parameter persistence, reproducibility, and UI state management.
 * Extends MetaBaseParameter to provide unified parameter architecture.
 */
public class DDAMagFragPipeParameters extends MetaBaseParameter {
    
    // Input/Output settings
    private List<ExperimentRow> experimentData = new ArrayList<>();
    
    // Additional FragPipe-specific parameters
    private MaxquantModification[] isobaric = new MaxquantModification[0];
    private boolean combineLabel = false;
    private double[][] isoCorFactor = new double[0][0];
    
    // Add direct MagDbItem storage
    private MagDbItem magDbItem;
    
    // FragPipe workflow configuration
    private FragpipeWorkflow fragpipeWorkflow;
    
    // Status label for UI feedback
    private Label statusLabel;

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
    
    // MAG-specific fields
    private String magDbCatalogueId = "";
    private String magDbVersion = "";
    
    // Isotopic labeling
    private List<String> lightLabels = new ArrayList<>();
    private List<String> mediumLabels = new ArrayList<>();
    private List<String> heavyLabels = new ArrayList<>();

	// Default constructor with sensible defaults for MAG FragPipe
	public DDAMagFragPipeParameters() {
		super(uk.ac.quadram.metapilot.MetaPilotApp.WorkflowType.DDA_MAG_FRAGPIPE.getAnalysisType(), new MetaData(), "", "", false, "", "HCD-FTMS",
				Runtime.getRuntime().availableProcessors(), Runtime.getRuntime().availableProcessors(), true);

		// Initialize FragPipe-specific defaults
		setEnzyme("Trypsin");
		setMissCleavages(2);
		setDigestMode("Specific");
		setMs2ScanMode(MetaConstants.HCD_FTMS);
		setThreadCount(Runtime.getRuntime().availableProcessors());
		setQuantitativeMode(MetaConstants.labelFree);
		setReduceDb(true);
		setSpectraCluster(false);

		// Initialize with default modifications for FragPipe
		setFixMods(new String[] { "Carbamidomethyl[C] 0" });
		setVariMods(new String[] { "Acetyl[ProteinN-term] 0", "Oxidation[M] 0" });
	}

	// Constructor for creating from GUI parameters
	public DDAMagFragPipeParameters(String analysisType, MetaData metadata, String result, String microDb,
			boolean appendHostDb, String hostDB, String ms2ScanMode, int threadCount) {
		super(analysisType, metadata, result, microDb, appendHostDb, hostDB, ms2ScanMode, threadCount, threadCount,
				true);

		// Set FragPipe-specific defaults
		setEnzyme("Trypsin");
		setMissCleavages(2);
		setDigestMode("Specific");
		setQuantitativeMode(MetaConstants.labelFree);
		setReduceDb(true);
		setSpectraCluster(false);
	}
    
    // Getters and setters for inherited functionality
    public List<ExperimentRow> getExperimentData() { return experimentData; }
    public void setExperimentData(List<ExperimentRow> experimentData) { this.experimentData = experimentData; }
    
    // Use inherited methods for modifications - convert between String[] and List<String>
    public List<String> getFixedModifications() {
        String[] fixMods = getFixMods();
        return fixMods != null ? List.of(fixMods) : new ArrayList<>();
    }
    
    public void setFixedModifications(List<String> fixedModifications) {
        if (fixedModifications != null) {
            setFixMods(fixedModifications.toArray(new String[0]));
        }
    }
    
    public List<String> getVariableModifications() {
        String[] variMods = getVariMods();
        return variMods != null ? List.of(variMods) : new ArrayList<>();
    }
    
    public void setVariableModifications(List<String> variableModifications) {
        if (variableModifications != null) {
            setVariMods(variableModifications.toArray(new String[0]));
        }
    }
    
    // Isobaric labeling getters and setters
    public String getIsobaricType() { return isobaricType; }
    public void setIsobaricType(String isobaricType) { this.isobaricType = isobaricType; }
    
    public List<IsobaricLabelRow> getIsobaricTableData() { return isobaricTableData; }
    public void setIsobaricTableData(List<IsobaricLabelRow> isobaricTableData) { this.isobaricTableData = isobaricTableData; }
    
    public MaxquantModification[] getIsobaric() { return isobaric; }
    public void setIsobaric(MaxquantModification[] isobaric) { this.isobaric = isobaric; }
    
    public boolean isCombineLabel() { return combineLabel; }
    public void setCombineLabel(boolean combineLabel) { this.combineLabel = combineLabel; }
    
    public double[][] getIsoCorFactor() { return isoCorFactor; }
    public void setIsoCorFactor(double[][] isoCorFactor) { this.isoCorFactor = isoCorFactor; }
    
    // MagDbItem getters and setters
    public MagDbItem getMagDbItem() { return magDbItem; }
    public void setMagDbItem(MagDbItem magDbItem) { this.magDbItem = magDbItem; }
    
    // Status label getters and setters
    public Label getStatusLabel() { return statusLabel; }
    public void setStatusLabel(Label statusLabel) { this.statusLabel = statusLabel; }
    
    // Taxonomic analysis getters and setters
    public boolean isEukaryota() { return eukaryota; }
    public void setEukaryota(boolean eukaryota) { this.eukaryota = eukaryota; }
    
    public boolean isBacteria() { return bacteria; }
    public void setBacteria(boolean bacteria) { this.bacteria = bacteria; }
    
    public boolean isArchaea() { return archaea; }
    public void setArchaea(boolean archaea) { this.archaea = archaea; }
    
    public boolean isViruses() { return viruses; }
    public void setViruses(boolean viruses) { this.viruses = viruses; }
    
    public boolean isViroids() { return viroids; }
    public void setViroids(boolean viroids) { this.viroids = viroids; }
    
    public boolean isUnclassified() { return unclassified; }
    public void setUnclassified(boolean unclassified) { this.unclassified = unclassified; }
    
    public boolean isOther() { return other; }
    public void setOther(boolean other) { this.other = other; }
    
    public String getIgnoreRank() { return ignoreRank; }
    public void setIgnoreRank(String ignoreRank) { this.ignoreRank = ignoreRank; }
    
    public boolean isIgnoreEnvironmentalSamples() { return ignoreEnvironmentalSamples; }
    public void setIgnoreEnvironmentalSamples(boolean ignoreEnvironmentalSamples) { this.ignoreEnvironmentalSamples = ignoreEnvironmentalSamples; }
    
    // Metadata getters and setters
    public List<MetadataRow> getMetadataTableData() { return metadataTableData; }
    public void setMetadataTableData(List<MetadataRow> metadataTableData) { this.metadataTableData = metadataTableData; }
    
    // MAG-specific getters and setters
    public String getMagDbCatalogueId() { return magDbCatalogueId; }
    public void setMagDbCatalogueId(String magDbCatalogueId) { this.magDbCatalogueId = magDbCatalogueId; }
    
    public String getMagDbVersion() { return magDbVersion; }
    public void setMagDbVersion(String magDbVersion) { this.magDbVersion = magDbVersion; }
    
    // Isotopic labeling getters and setters
    public List<String> getLightLabels() { return lightLabels; }
    public void setLightLabels(List<String> lightLabels) { this.lightLabels = lightLabels; }
    
    public List<String> getMediumLabels() { return mediumLabels; }
    public void setMediumLabels(List<String> mediumLabels) { this.mediumLabels = mediumLabels; }
    
    public List<String> getHeavyLabels() { return heavyLabels; }
    public void setHeavyLabels(List<String> heavyLabels) { this.heavyLabels = heavyLabels; }
 
    // Isobaric tag getter
    public uk.ac.quadram.metapilot.core.mod.IsobaricTag getIsobaricTag() {
        // Return null for now - this needs to be implemented based on the isobaric data
        return null;
    }
    
    // FragPipe workflow getters and setters
    public FragpipeWorkflow getFragpipeWorkflow() { return fragpipeWorkflow; }
    public void setFragpipeWorkflow(FragpipeWorkflow fragpipeWorkflow) { this.fragpipeWorkflow = fragpipeWorkflow; }
    
    /**
     * Convert to MetaParameterMag for compatibility with existing MAG tasks
     */
    public uk.ac.quadram.metapilot.task.mag.par.MetaParameterMag toMetaParameterMag() {
        // Create MetaBaseParameter first
        uk.ac.quadram.metapilot.task.par.MetaBaseParameter metaBasePar = new uk.ac.quadram.metapilot.task.par.MetaBaseParameter(
            getWorkflowType().name(),
            getMetadata(),
            getResult(),
            getMicroDb(),
            isAppendHostDb(),
            getHostDB(),
            getMs2ScanMode(),
            getCoreCount(),
            getThreadCount(),
            isMetaWorkflow()
        );
        
        // Convert digestion mode string to integer
        int digestModeInt = getDigestModeInt(getDigestModeString());
        
        // Create MetaParameterMag
        return new uk.ac.quadram.metapilot.task.mag.par.MetaParameterMag(
            metaBasePar,
            getFixMods(),
            getVariMods(),
            getEnzyme(),
            getMissCleavages(),
            digestModeInt,
            getQuantitativeMode(),
            getIsobaric(),
            getIsobaricTag(),
            isCombineLabel(),
            getIsoCorFactor(),
            getMagDbCatalogueId(),
            getMagDbVersion()
        );
    }
    
    /**
     * Convert digestion mode string to integer representation
     */
    private int getDigestModeInt(String digestMode) {
        switch (digestMode) {
            case "Specific": return 1;
            case "Semi free N-term": return 2;
            case "Semi free C-term": return 3;
            case "Semi specific": return 4;
            case "Unspecific": return 5;
            case "None": return 6;
            default: return 1;
        }
    }
}
