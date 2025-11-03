package uk.ac.quadram.metapilot.task.dda.igc.par;

import java.util.ArrayList;
import java.util.List;

import javafx.scene.control.Label;
import uk.ac.quadram.metapilot.dbSearch.MaxQuant.MaxquantModification;
import uk.ac.quadram.metapilot.task.common.gui.MetaPilotIOPanelFX.ExperimentRow;
import uk.ac.quadram.metapilot.task.common.gui.IsobaricLabelingPanel.IsobaricLabelRow;
import uk.ac.quadram.metapilot.task.dda.igc.gui.MetaPilotDdaIgcMaxQuantPanel.MetadataRow;
import uk.ac.quadram.metapilot.task.par.MetaBaseParameter;
import uk.ac.quadram.metapilot.task.par.MetaConstants;
import uk.ac.quadram.metapilot.task.par.MetaData;
import uk.ac.quadram.metapilot.dbSearch.fragpipe.FragpipeWorkflow;

/**
 * Parameter container for DDA IGC FragPipe workflow
 * 
 * This class holds all configuration parameters for a DDA IGC FragPipe analysis,
 * enabling parameter persistence, reproducibility, and UI state management.
 * Extends MetaBaseParameter to provide unified parameter architecture.
 */
public class DDAIgcFragPipeParameters extends MetaBaseParameter {
    
    // Input/Output settings
    private List<ExperimentRow> experimentData = new ArrayList<>();
    
    // Additional FragPipe-specific parameters
    private MaxquantModification[] isobaric = new MaxquantModification[0];
    private boolean combineLabel = false;
    private double[][] isoCorFactor = new double[0][0];
    
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
    
    // Isotopic labeling
    private List<String> lightLabels = new ArrayList<>();
    private List<String> mediumLabels = new ArrayList<>();
    private List<String> heavyLabels = new ArrayList<>();

    // Default constructor with sensible defaults for IGC FragPipe
    public DDAIgcFragPipeParameters() {
        super(uk.ac.quadram.metapilot.MetaPilotApp.WorkflowType.DDA_IGC_FRAGPIPE.getAnalysisType(), new MetaData(), "", "", false, "", "HCD-FTMS",
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
    public DDAIgcFragPipeParameters(String analysisType, MetaData metadata, String result, String microDb,
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
        setFixMods(fixedModifications.toArray(new String[0]));
    }
    
    public List<String> getVariableModifications() {
        String[] variMods = getVariMods();
        return variMods != null ? List.of(variMods) : new ArrayList<>();
    }
    
    public void setVariableModifications(List<String> variableModifications) {
        setVariMods(variableModifications.toArray(new String[0]));
    }
    
    // FragPipe-specific methods
    public FragpipeWorkflow getFragpipeWorkflow() {
        return fragpipeWorkflow;
    }
    
    public void setFragpipeWorkflow(FragpipeWorkflow fragpipeWorkflow) {
        this.fragpipeWorkflow = fragpipeWorkflow;
    }
    
    public Label getStatusLabel() {
        return statusLabel;
    }
    
    public void setStatusLabel(Label statusLabel) {
        this.statusLabel = statusLabel;
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
    
    public String getQuanMode() {
        return getQuantitativeMode();
    }
    
    public void setQuanMode(String quanMode) {
        setQuantitativeMode(quanMode);
    }
    
    public String getIsobaricType() { return isobaricType; }
    public void setIsobaricType(String isobaricType) { this.isobaricType = isobaricType; }
    
    public List<IsobaricLabelRow> getIsobaricTableData() { return isobaricTableData; }
    public void setIsobaricTableData(List<IsobaricLabelRow> isobaricTableData) { 
        this.isobaricTableData = isobaricTableData; 
    }
    
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
    public void setIgnoreEnvironmentalSamples(boolean ignoreEnvironmentalSamples) { 
        this.ignoreEnvironmentalSamples = ignoreEnvironmentalSamples; 
    }
    
    public List<MetadataRow> getMetadataTableData() { return metadataTableData; }
    public void setMetadataTableData(List<MetadataRow> metadataTableData) { 
        this.metadataTableData = metadataTableData; 
    }
    
    public List<String> getLightLabels() { return lightLabels; }
    public void setLightLabels(List<String> lightLabels) { this.lightLabels = lightLabels; }
    
    public List<String> getMediumLabels() { return mediumLabels; }
    public void setMediumLabels(List<String> mediumLabels) { this.mediumLabels = mediumLabels; }
    
    public List<String> getHeavyLabels() { return heavyLabels; }
    public void setHeavyLabels(List<String> heavyLabels) { this.heavyLabels = heavyLabels; }
    
    /**
     * Create a copy of this parameters object
     */
    public DDAIgcFragPipeParameters copy() {
        DDAIgcFragPipeParameters copy = new DDAIgcFragPipeParameters();
        
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
        
        // Copy FragPipe-specific fields
        copy.experimentData = new ArrayList<>(this.experimentData);
        copy.isobaric = this.isobaric.clone();
        copy.combineLabel = this.combineLabel;
        copy.isoCorFactor = this.isoCorFactor.clone();
        copy.fragpipeWorkflow = this.fragpipeWorkflow;
        copy.statusLabel = this.statusLabel;
        
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
        
        copy.lightLabels = new ArrayList<>(this.lightLabels);
        copy.mediumLabels = new ArrayList<>(this.mediumLabels);
        copy.heavyLabels = new ArrayList<>(this.heavyLabels);
        
        return copy;
    }
    
    @Override
    public String toString() {
        return String.format("DDAIgcFragPipeParameters[enzyme=%s, missed_cleavages=%d, quant_mode=%s, " +
                           "output=%s, threads=%d]", 
                           getEnzyme(), getMissCleavages(), getQuantitativeMode(), getResult(), getThreadCount());
    }
}
