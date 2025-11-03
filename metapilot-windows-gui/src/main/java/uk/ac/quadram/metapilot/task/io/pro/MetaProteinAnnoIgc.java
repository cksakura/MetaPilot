package uk.ac.quadram.metapilot.task.io.pro;

/**
 * Represents a protein with functional annotations for IGC (Integrated Gene Catalog) workflows.
 * 
 * This class extends MetaProtein to include functional annotation information
 * specific to IGC database searches, including COG, NOG, KEGG, GO, and EC annotations.
 */
public class MetaProteinAnnoIgc extends MetaProtein {

    // Functional annotation fields
    private String[] cogAnnotations;
    private String[] nogAnnotations;
    private String keggAnnotation;
    private String goAnnotation;
    private String ecAnnotation;
    
    // Additional IGC-specific fields
    private String geneId;
    private String geneSymbol;
    private String geneDescription;
    private String organism;
    private String taxonomy;

    /**
     * Default constructor
     */
    public MetaProteinAnnoIgc() {
        super("", new double[0]);
    }

    /**
     * Constructor with basic protein information
     */
    public MetaProteinAnnoIgc(String name, double[] intensity) {
        super(name, intensity);
    }

    // Getters and setters for functional annotations
    public String[] getCogAnnotations() {
        return cogAnnotations;
    }

    public void setCogAnnotations(String[] cogAnnotations) {
        this.cogAnnotations = cogAnnotations;
    }

    public String[] getNogAnnotations() {
        return nogAnnotations;
    }

    public void setNogAnnotations(String[] nogAnnotations) {
        this.nogAnnotations = nogAnnotations;
    }

    public String getKeggAnnotation() {
        return keggAnnotation;
    }

    public void setKeggAnnotation(String keggAnnotation) {
        this.keggAnnotation = keggAnnotation;
    }

    public String getGoAnnotation() {
        return goAnnotation;
    }

    public void setGoAnnotation(String goAnnotation) {
        this.goAnnotation = goAnnotation;
    }

    public String getEcAnnotation() {
        return ecAnnotation;
    }

    public void setEcAnnotation(String ecAnnotation) {
        this.ecAnnotation = ecAnnotation;
    }

    // Getters and setters for IGC-specific fields
    public String getGeneId() {
        return geneId;
    }

    public void setGeneId(String geneId) {
        this.geneId = geneId;
    }

    public String getGeneSymbol() {
        return geneSymbol;
    }

    public void setGeneSymbol(String geneSymbol) {
        this.geneSymbol = geneSymbol;
    }

    public String getGeneDescription() {
        return geneDescription;
    }

    public void setGeneDescription(String geneDescription) {
        this.geneDescription = geneDescription;
    }

    public String getOrganism() {
        return organism;
    }

    public void setOrganism(String organism) {
        this.organism = organism;
    }

    public String getTaxonomy() {
        return taxonomy;
    }

    public void setTaxonomy(String taxonomy) {
        this.taxonomy = taxonomy;
    }
}
