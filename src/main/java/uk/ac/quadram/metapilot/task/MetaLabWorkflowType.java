package uk.ac.quadram.metapilot.task;

/**
 * MetaPilot Workflow Type Enumeration - Comprehensive Analysis Pipeline Definitions
 * 
 * <h2>Purpose</h2>
 * This enumeration defines all supported analysis workflows in MetaPilot, providing
 * a type-safe way to identify and configure different metaproteomics analysis
 * pipelines. Each workflow type represents a specific combination of search engine,
 * database type, and analysis methodology.
 * 
 * <h2>Workflow Categories</h2>
 * 
 * <h3>Classical Proteomics Workflows:</h3>
 * <ul>
 *   <li><strong>MaxQuantWorkflow:</strong> MaxQuant-based label-free quantification</li>
 *   <li><strong>MSFraggerWorkflow:</strong> MSFragger database search pipeline</li>
 *   <li><strong>pFindWorkflow:</strong> pFind database search and identification</li>
 * </ul>
 * 
 * <h3>Metaproteomics Analysis Workflows:</h3>
 * <ul>
 *   <li><strong>TaxonomyAnalysis:</strong> Taxonomic classification of identified proteins</li>
 *   <li><strong>FunctionalAnnotation:</strong> Functional annotation and pathway analysis</li>
 * </ul>
 * 
 * <h3>Human Gut Microbiome (HGM) Workflows:</h3>
 * <ul>
 *   <li><strong>pFindHGM:</strong> pFind workflow optimized for human gut microbiome analysis</li>
 * </ul>
 * 
 * <h3>Metagenome-Assembled Genomes (MAG) Workflows:</h3>
 * <h4>DDA (Data-Dependent Acquisition) MAG Workflows:</h4>
 * <ul>
 *   <li><strong>pFindMAG:</strong> pFind search against MAG protein databases</li>
 *   <li><strong>FragpipeIMMag:</strong> FragPipe analysis for MAG proteomes</li>
 *   <li><strong>SageMAG:</strong> Sage high-speed search for MAG databases</li>
 *   <li><strong>AlphapeptIMMag:</strong> AlphaPept workflow for MAG analysis</li>
 * </ul>
 * 
 * <h4>DIA (Data-Independent Acquisition) MAG Workflows:</h4>
 * <ul>
 *   <li><strong>DiaNNMAG:</strong> Standard DIA-NN analysis for MAG proteomes</li>
 *   <li><strong>DiaNNMAGRT:</strong> Real-time DIA-NN monitoring for MAG analysis</li>
 *   <li><strong>DiaNNMAGTG:</strong> Targeted DIA-NN analysis for specific MAG proteins</li>
 * </ul>
 * 
 * <h3>Combined Workflows:</h3>
 * <ul>
 *   <li><strong>DdaDiaMAG:</strong> Integrated DDA/DIA workflow for comprehensive MAG analysis</li>
 * </ul>
 * 
 * <h2>Workflow Characteristics</h2>
 * 
 * <h3>Database Types:</h3>
 * <ul>
 *   <li><strong>Classical:</strong> Standard protein databases (UniProt, RefSeq)</li>
 *   <li><strong>HGM:</strong> Human gut microbiome-specific databases</li>
 *   <li><strong>MAG:</strong> Metagenome-assembled genome protein collections</li>
 * </ul>
 * 
 * <h3>Acquisition Methods:</h3>
 * <ul>
 *   <li><strong>DDA:</strong> Data-dependent acquisition (shotgun proteomics)</li>
 *   <li><strong>DIA:</strong> Data-independent acquisition (SWATH-MS)</li>
 *   <li><strong>Mixed:</strong> Combined DDA/DIA approaches</li>
 * </ul>
 * 
 * <h3>Analysis Focus:</h3>
 * <ul>
 *   <li><strong>Identification:</strong> Protein and peptide identification</li>
 *   <li><strong>Quantification:</strong> Label-free protein quantification</li>
 *   <li><strong>Taxonomy:</strong> Species-level classification</li>
 *   <li><strong>Function:</strong> Biological pathway analysis</li>
 * </ul>
 * 
 * <h2>Search Engine Integration</h2>
 * 
 * <h3>Supported Search Engines:</h3>
 * <ul>
 *   <li><strong>MaxQuant:</strong> Comprehensive proteomics platform</li>
 *   <li><strong>pFind:</strong> Chinese-developed high-performance engine</li>
 *   <li><strong>MSFragger:</strong> Ultra-fast database search</li>
 *   <li><strong>Sage:</strong> Rust-based high-speed search</li>
 *   <li><strong>DIA-NN:</strong> Deep learning DIA analysis</li>
 *   <li><strong>AlphaPept:</strong> Python-based proteomics platform</li>
 * </ul>
 * 
 * <h2>Workflow Selection Guidelines</h2>
 * 
 * <h3>For Standard Proteomics:</h3>
 * <ul>
 *   <li>Use <strong>MaxQuantWorkflow</strong> for comprehensive analysis</li>
 *   <li>Use <strong>MSFraggerWorkflow</strong> for high-throughput studies</li>
 * </ul>
 * 
 * <h3>For Metaproteomics:</h3>
 * <ul>
 *   <li>Use <strong>pFindHGM</strong> for human gut microbiome samples</li>
 *   <li>Use <strong>MAG workflows</strong> for environmental/clinical samples</li>
 * </ul>
 * 
 * <h3>For DIA Analysis:</h3>
 * <ul>
 *   <li>Use <strong>DiaNNMAG</strong> for standard DIA metaproteomics</li>
 *   <li>Use <strong>DiaNNMAGRT</strong> for real-time monitoring</li>
 *   <li>Use <strong>DiaNNMAGTG</strong> for targeted analysis</li>
 * </ul>
 * 
 * <h2>Performance Characteristics</h2>
 * <ul>
 *   <li><strong>pFind/Sage:</strong> Fastest search engines</li>
 *   <li><strong>MaxQuant:</strong> Most comprehensive but slower</li>
 *   <li><strong>DIA-NN:</strong> Optimized for DIA data</li>
 *   <li><strong>MAG workflows:</strong> Designed for large metagenomic databases</li>
 * </ul>
 * 
 * <h2>Usage Examples</h2>
 * 
 * <h3>Workflow Selection:</h3>
 * <pre>{@code
 *   // Select workflow based on sample type
 *   MetaLabWorkflowType workflow;
 *   if (isHumanGutSample) {
 *       workflow = MetaLabWorkflowType.pFindHGM;
 *   } else if (isEnvironmentalSample) {
 *       workflow = MetaLabWorkflowType.pFindMAG;
 *   } else {
 *       workflow = MetaLabWorkflowType.MaxQuantWorkflow;
 *   }
 * }</pre>
 * 
 * <h3>Workflow Properties:</h3>
 * <pre>{@code
 *   MetaLabWorkflowType workflow = MetaLabWorkflowType.pFindMAG;
 *   int id = workflow.getId();
 *   String description = workflow.getDescription();
 *   
 *   // Parse workflow from ID
 *   MetaLabWorkflowType parsed = MetaLabWorkflowType.parse(6);
 * }</pre>
 * 
 * <h2>Future Extensions</h2>
 * <p>This enumeration is designed to be extensible for new workflows:</p>
 * <ul>
 *   <li>Additional search engines (Comet, X!Tandem)</li>
 *   <li>Specialized databases (plant, marine microbiomes)</li>
 *   <li>Novel acquisition methods (timsTOF, FAIMS)</li>
 *   <li>Integrated multi-omics workflows</li>
 * </ul>
 * 
 * @author Kai Cheng - Quadram Institute Bioscience
 * @version 2.0
 * @since 1.0
 * 
 * @see uk.ac.quadram.metapilot.task.MetaAbstractTask
 * @see uk.ac.quadram.metapilot.MetaPilotApp
 */
public enum MetaLabWorkflowType {

	MaxQuantWorkflow(0, "MaxQuant workflow"),

	MSFraggerWorkflow(1, "MSFragger workflow"),

	pFindWorkflow(2, "pFind workflow"),

	TaxonomyAnalysis(3, "Taxonomy analysis"),

	FunctionalAnnotation(4, "Functional annotation"),

	pFindHGM(5, "pFind human gut microbiome (HGM) workflow"),

	pFindMAG(6, "pFind metagenome-assembled genomes (MAG) workflow"),

	DiaNNMAG(7, "DiaNN metagenome-assembled genomes (MAG) workflow"),

	DdaDiaMAG(8, "DDA Dia workflow"),

	FragpipeIMMag(9, "Fragpipe metagenome-assembled genomes (MAG) workflow"),

	AlphapeptIMMag(10, "Alphapept metagenome-assembled genomes (MAG) workflow"),
	
	DiaNNMAGRT(11, "DiaNN metagenome-assembled genomes (MAG) realtime workflow"),
	
	DiaNNMAGTG(12, "DiaNN metagenome-assembled genomes (MAG) target workflow"),
	
	SageMAG(13, "Sage metagenome-assembled genomes (MAG) workflow");

	int id;
	String name;

	MetaLabWorkflowType(int id, String name) {
		this.id = id;
		this.name = name;
	}

	public int getId() {
		return id;
	}

	public String getDescription() {
		return name;
	}

	public static MetaLabWorkflowType parse(int id) {
		for (MetaLabWorkflowType tt : MetaLabWorkflowType.values()) {
			if (tt.getId() == id) {
				return tt;
			}
		}
		return null;
	}

}
