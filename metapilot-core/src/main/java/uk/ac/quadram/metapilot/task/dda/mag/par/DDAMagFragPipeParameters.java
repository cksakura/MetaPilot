package uk.ac.quadram.metapilot.task.dda.mag.par;

import java.util.ArrayList;
import java.util.List;

import uk.ac.quadram.metapilot.dbSearch.MaxQuant.MaxquantModification;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;
import uk.ac.quadram.metapilot.task.par.MetaBaseParameter;
import uk.ac.quadram.metapilot.task.par.MetaConstants;
import uk.ac.quadram.metapilot.task.par.MetaData;

/**
 * Parameter container for DDA MAG FragPipe workflow
 * 
 * This class holds all configuration parameters for a DDA MAG FragPipe analysis,
 * enabling parameter persistence, reproducibility, and UI state management.
 * Extends MetaBaseParameter to provide unified parameter architecture.
 */
public class DDAMagFragPipeParameters extends MetaBaseParameter {
    
    // Input/Output settings - now uses String[] arrays instead of JavaFX row objects
    // ExperimentData format: [filePath, experimentName, fraction]
    private List<String[]> experimentData = new ArrayList<>();
    
    // Additional FragPipe-specific parameters
    private MaxquantModification[] isobaric = new MaxquantModification[0];
    private boolean combineLabel = false;
    private double[][] isoCorFactor = new double[0][0];
    
    // Add direct MagDbItem storage
    private MagDbItem magDbItem;

    // Isobaric labeling - now uses String[] arrays instead of JavaFX row objects
    // IsobaricLabelData format: [labelName, minus2, minus1, plus1, plus2]
    private String isobaricType = "";
    private List<String[]> isobaricTableData = new ArrayList<>();
    
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
    
    // Metadata - now uses String[] arrays instead of JavaFX row objects
    // MetadataData format: [sampleName, labelType, age, gender, condition, treatment, timePoint, customField1, customField2, customField3]
    private List<String[]> metadataTableData = new ArrayList<>();
    
    // MAG-specific fields
    private String magDbCatalogueId = "";
    /** Registry identifier (e.g. MGYG for MGnify, or a custom local-MAG identifier). */
    private String magDbIdentifier = "MGYG";
    private String magDbVersion = "";
    private String genomeSelectorMode = MagDbItem.MARKER2;
    private String genomeSelectionStrategy = "Focus";
    private int fragpipeParallelTasks = 0;

    // Search-mode configuration (Closed vs Open search). Closed is the default and
    // reproduces the historical MSFragger closed-search behavior exactly. Open search
    // is a "must" for the combined/quant stage when selected, but is independently
    // optional at the marker/genome-selection stage (see openSearchAtMarkerStage).
    private String searchMode = "Closed";
    /** Only meaningful when searchMode == "Open"; off by default. */
    private boolean openSearchAtMarkerStage = false;
    /** Power-user escape hatch: explicit .workflow file for the marker/sample-specific stage. */
    private String customBasicWorkflowPath = "";
    /** Power-user escape hatch: explicit .workflow file for the combined/quant stage. */
    private String customCombinedWorkflowPath = "";
    
    // Isotopic labeling
    private List<String> lightLabels = new ArrayList<>();
    private List<String> mediumLabels = new ArrayList<>();
    private List<String> heavyLabels = new ArrayList<>();

	// Default constructor with sensible defaults for MAG FragPipe
	public DDAMagFragPipeParameters() {
		super("DDA_MAG_FRAGPIPE", new MetaData(), "", "", false, "", "HCD-FTMS",
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
    public List<String[]> getExperimentData() { return experimentData; }
    public void setExperimentData(List<String[]> experimentData) { this.experimentData = experimentData; }
    
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
    
    public List<String[]> getIsobaricTableData() { return isobaricTableData; }
    public void setIsobaricTableData(List<String[]> isobaricTableData) { this.isobaricTableData = isobaricTableData; }
    
    public MaxquantModification[] getIsobaric() { return isobaric; }
    public void setIsobaric(MaxquantModification[] isobaric) { this.isobaric = isobaric; }
    
    public boolean isCombineLabel() { return combineLabel; }
    public void setCombineLabel(boolean combineLabel) { this.combineLabel = combineLabel; }
    
    public double[][] getIsoCorFactor() { return isoCorFactor; }
    public void setIsoCorFactor(double[][] isoCorFactor) { this.isoCorFactor = isoCorFactor; }
    
    // MagDbItem getters and setters
    public MagDbItem getMagDbItem() { return magDbItem; }
    public void setMagDbItem(MagDbItem magDbItem) {
        this.magDbItem = magDbItem;
        if (magDbItem != null) {
            if (magDbItem.getCatalogueID() != null && !magDbItem.getCatalogueID().isBlank()) {
                this.magDbCatalogueId = magDbItem.getCatalogueID();
            }
            if (magDbItem.getIdentifier() != null && !magDbItem.getIdentifier().isBlank()) {
                this.magDbIdentifier = magDbItem.getIdentifier();
            }
            if (magDbItem.getUsedVersion() != null
                    && magDbItem.getUsedVersion().getCatalogueVersion() != null
                    && !magDbItem.getUsedVersion().getCatalogueVersion().isBlank()) {
                this.magDbVersion = magDbItem.getUsedVersion().getCatalogueVersion();
            }
        }
    }
    
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
    public List<String[]> getMetadataTableData() { return metadataTableData; }
    public void setMetadataTableData(List<String[]> metadataTableData) { this.metadataTableData = metadataTableData; }
    
    // MAG-specific getters and setters
    public String getMagDbCatalogueId() { return magDbCatalogueId; }
    public void setMagDbCatalogueId(String magDbCatalogueId) { this.magDbCatalogueId = magDbCatalogueId; }

    public String getMagDbIdentifier() {
        return magDbIdentifier == null || magDbIdentifier.isBlank() ? "MGYG" : magDbIdentifier;
    }
    public void setMagDbIdentifier(String magDbIdentifier) {
        this.magDbIdentifier = magDbIdentifier == null || magDbIdentifier.isBlank() ? "MGYG" : magDbIdentifier.trim();
    }
    
    public String getMagDbVersion() { return magDbVersion; }
    public void setMagDbVersion(String magDbVersion) { this.magDbVersion = magDbVersion; }

    public String getSearchMode() { return searchMode; }
    public void setSearchMode(String searchMode) {
        this.searchMode = "Open".equalsIgnoreCase(searchMode) ? "Open" : "Closed";
    }
    public boolean isOpenSearch() { return "Open".equals(searchMode); }

    public boolean isOpenSearchAtMarkerStage() { return openSearchAtMarkerStage; }
    public void setOpenSearchAtMarkerStage(boolean openSearchAtMarkerStage) {
        this.openSearchAtMarkerStage = openSearchAtMarkerStage;
    }

    public String getCustomBasicWorkflowPath() { return customBasicWorkflowPath; }
    public void setCustomBasicWorkflowPath(String customBasicWorkflowPath) {
        this.customBasicWorkflowPath = customBasicWorkflowPath == null ? "" : customBasicWorkflowPath.trim();
    }

    public String getCustomCombinedWorkflowPath() { return customCombinedWorkflowPath; }
    public void setCustomCombinedWorkflowPath(String customCombinedWorkflowPath) {
        this.customCombinedWorkflowPath = customCombinedWorkflowPath == null ? "" : customCombinedWorkflowPath.trim();
    }

    public String getGenomeSelectorMode() { return genomeSelectorMode; }
    public void setGenomeSelectorMode(String genomeSelectorMode) {
        this.genomeSelectorMode = genomeSelectorMode != null && !genomeSelectorMode.trim().isEmpty()
                ? genomeSelectorMode.trim() : MagDbItem.MARKER2;
    }

    public String getGenomeSelectionStrategy() { return genomeSelectionStrategy; }
    public void setGenomeSelectionStrategy(String genomeSelectionStrategy) {
        if (genomeSelectionStrategy == null || genomeSelectionStrategy.trim().isEmpty()) {
            this.genomeSelectionStrategy = "Focus";
        } else {
            this.genomeSelectionStrategy = genomeSelectionStrategy.trim();
        }
    }

    public int getFragpipeParallelTasks() { return fragpipeParallelTasks; }
    public void setFragpipeParallelTasks(int fragpipeParallelTasks) {
        this.fragpipeParallelTasks = Math.max(0, fragpipeParallelTasks);
    }
    
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
