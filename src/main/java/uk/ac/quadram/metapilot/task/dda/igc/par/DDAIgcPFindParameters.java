package uk.ac.quadram.metapilot.task.dda.igc.par;

import java.util.ArrayList;
import java.util.List;

import uk.ac.quadram.metapilot.dbSearch.MaxQuant.MaxquantModification;
import uk.ac.quadram.metapilot.task.common.gui.MetaPilotIOPanelFX.ExperimentRow;
import uk.ac.quadram.metapilot.task.common.gui.IsobaricLabelingPanel.IsobaricLabelRow;
import uk.ac.quadram.metapilot.task.dda.igc.gui.MetaPilotDdaIgcMaxQuantPanel.MetadataRow;
import uk.ac.quadram.metapilot.task.par.MetaBaseParameter;
import uk.ac.quadram.metapilot.task.par.MetaConstants;
import uk.ac.quadram.metapilot.task.par.MetaData;
import uk.ac.quadram.metapilot.task.par.MetaPep2TaxaPar;
import uk.ac.quadram.metapilot.task.v2.par.MetaOpenPar;

/**
 * Parameter container for DDA IGC pFind workflow
 * 
 * This class holds all configuration parameters for a DDA IGC pFind analysis,
 * enabling parameter persistence, reproducibility, and UI state management.
 * Extends MetaBaseParameter to provide unified parameter architecture.
 */
public class DDAIgcPFindParameters extends MetaBaseParameter {
    
    // Input/Output settings
    private List<ExperimentRow> experimentData = new ArrayList<>();
    
    // Additional pFind-specific parameters
    private MaxquantModification[] isobaric = new MaxquantModification[0];
    private boolean combineLabel = false;
    private double[][] isoCorFactor = new double[0][0];
    
    // Open search and taxonomic parameters
    private MetaOpenPar osp = MetaOpenPar.getDefault();
    private MetaPep2TaxaPar ptp = MetaPep2TaxaPar.getDefault();
    
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
    
    // Default constructor with sensible defaults for pFind
    public DDAIgcPFindParameters() {
        super("DDA_IGC_PFIND", new MetaData(), "", "", false, "", MetaConstants.HCD_FTMS, 1, 1, true);
        
        // Initialize pFind-specific defaults
        setEnzyme("Trypsin_P");
        setMissCleavages(2);
        setDigestMode("Specific");
        setMs2ScanMode(MetaConstants.HCD_FTMS);
        setThreadCount(Runtime.getRuntime().availableProcessors());
        setQuantitativeMode(MetaConstants.labelFree);
        setReduceDb(true);
        setSpectraCluster(false);
        
        // Initialize with default modifications for pFind
        setFixMods(new String[]{"Carbamidomethyl[C] 0"});
        setVariMods(new String[]{"Acetyl[ProteinN-term] 0", "Oxidation[M] 0"});
    }
    
    // Constructor for creating from GUI parameters
    public DDAIgcPFindParameters(String analysisType, MetaData metadata, String result, String microDb,
            boolean appendHostDb, String hostDB, String ms2ScanMode, int threadCount) {
        super(analysisType, metadata, result, microDb, appendHostDb, hostDB, ms2ScanMode, threadCount, threadCount, true);
        
        // Set pFind-specific defaults
        setEnzyme("Trypsin_P");
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
    
    // pFind-specific methods (compatible with MetaParameterPFind)
    public boolean isOpenSearch() {
        return osp.isOpenSearch();
    }
    
    public double getOpenPsmFDR() {
        return osp.getOpenPsmFDR();
    }
    
    public double getOpenProFDR() {
        return osp.getOpenProFDR();
    }
    
    public MetaOpenPar getMop() {
        return osp;
    }
    
    public void setMop(MetaOpenPar mop) {
        this.osp = mop;
    }
    
    public MetaOpenPar getOsp() {
        return osp;
    }
    
    public void setOsp(MetaOpenPar osp) {
        this.osp = osp;
    }
    
    public MetaPep2TaxaPar getPtp() {
        return ptp;
    }
    
    public void setPtp(MetaPep2TaxaPar ptp) {
        this.ptp = ptp;
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
    
    /**
     * Create a copy of this parameters object
     */
    public DDAIgcPFindParameters copy() {
        DDAIgcPFindParameters copy = new DDAIgcPFindParameters();
        
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
        
        // Copy pFind-specific fields
        copy.experimentData = new ArrayList<>(this.experimentData);
        copy.isobaric = this.isobaric.clone();
        copy.combineLabel = this.combineLabel;
        copy.isoCorFactor = this.isoCorFactor.clone();
        copy.osp = this.osp;
        copy.ptp = this.ptp;
        
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
    
    @Override
    public String toString() {
        return String.format("DDAIgcPFindParameters[enzyme=%s, missed_cleavages=%d, quant_mode=%s, " +
                           "output=%s, threads=%d]", 
                           getEnzyme(), getMissCleavages(), getQuantitativeMode(), getResult(), getThreadCount());
    }
}
