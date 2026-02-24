package uk.ac.quadram.metapilot.workflow;

/**
 * Enumeration of all MetaPilot workflow types.
 * 
 * This enum defines all available workflow types in MetaPilot, providing
 * a centralized, type-safe way to identify and manage different analysis
 * workflows across the application.
 * 
 * Workflow categories:
 * - DDA IGC: Data-Dependent Acquisition with Integrated Gene Catalog
 * - DDA MAG: Data-Dependent Acquisition with Metagenome-Assembled Genomes
 * - DIA: Data-Independent Acquisition
 * 
 * @author Kai Cheng - Quadram Institute Bioscience
 * @version 1.0.0-beta.6
 * @since 1.0.0-beta.3
 */
public enum WorkflowType {
    // DDA IGC Workflows
    DDA_IGC_MAXQUANT("DDA_IGC_MAXQUANT"),
    DDA_IGC_PFIND("DDA_IGC_PFIND"),
    DDA_IGC_FRAGPIPE("DDA_IGC_FRAGPIPE"),
    DDA_IGC_SAGE("DDA_IGC_SAGE"),
    
    // DDA MAG Workflows
    DDA_MAG_PFIND("DDA_MAG_PFIND"),
    DDA_MAG_FRAGPIPE("DDA_MAG_FRAGPIPE"),
    DDA_MAG_SAGE("DDA_MAG_SAGE"),
    
    // DIA Workflows
    /** Marker-based DIA workflow (Marker4 genome selection). */
    DIA_MARKER("DIA_FOCUS"),
    /** Seed/library-based DIA workflow (MAG library search). */
    DIA_SEED("DIA_NORMAL"),
    DIA_TARGET("DIA_TARGET"),
    DIA_REALTIME("DIA_REALTIME"),
    /** Deeper, more time-consuming comprehensive DIA analysis. */
    DIA_DEEP_COMPREHENSIVE("DIA_DEEP_COMPREHENSIVE");
    
    private final String analysisType;
    
    WorkflowType(String analysisType) {
        this.analysisType = analysisType;
    }
    
    /**
     * Get the analysis type string identifier for this workflow.
     * This is used for parameter persistence and workflow identification.
     * 
     * @return The analysis type string
     */
    public String getAnalysisType() {
        return analysisType;
    }
    
    /**
     * Get WorkflowType from its analysis type string.
     * 
     * @param analysisType The analysis type string
     * @return The corresponding WorkflowType, or null if not found
     */
    public static WorkflowType fromAnalysisType(String analysisType) {
        if (analysisType == null) {
            return null;
        }
        for (WorkflowType type : values()) {
            if (type.getAnalysisType().equals(analysisType)) {
                return type;
            }
        }
        return null;
    }
}

