package uk.ac.quadram.metapilot.core.function.v2;

import java.util.ArrayList;

/**
 * Enzyme Commission (EC) Number Classification System for Functional Annotation
 * 
 * <h2>Purpose</h2>
 * This class represents the Enzyme Commission (EC) number classification system,
 * which provides a standardized method for naming and classifying enzymes based
 * on their catalytic function. EC numbers are essential for functional annotation
 * in metaproteomics and systems biology applications.
 * 
 * <h2>EC Number System</h2>
 * <p>The EC number is a hierarchical classification scheme with four levels:</p>
 * <pre>{@code
 *   EC X.Y.Z.W
 *   │ │ │ └─ Serial number (individual enzyme)
 *   │ │ └─── Sub-subclass (specific bond/group)
 *   │ └───── Subclass (substrate type)
 *   └─────── Main class (reaction type)
 * }</pre>
 * 
 * <h3>Main Classes (First Digit):</h3>
 * <ul>
 *   <li><strong>EC 1.x.x.x:</strong> Oxidoreductases (oxidation-reduction reactions)</li>
 *   <li><strong>EC 2.x.x.x:</strong> Transferases (transfer of functional groups)</li>
 *   <li><strong>EC 3.x.x.x:</strong> Hydrolases (hydrolysis reactions)</li>
 *   <li><strong>EC 4.x.x.x:</strong> Lyases (addition/removal without hydrolysis)</li>
 *   <li><strong>EC 5.x.x.x:</strong> Isomerases (intramolecular rearrangements)</li>
 *   <li><strong>EC 6.x.x.x:</strong> Ligases (formation of bonds with ATP cleavage)</li>
 *   <li><strong>EC 7.x.x.x:</strong> Translocases (movement across membranes)</li>
 * </ul>
 * 
 * <h2>Data Model</h2>
 * <h3>Core Fields:</h3>
 * <ul>
 *   <li><strong>id:</strong> EC number identifier (e.g., "1.1.1.1")</li>
 *   <li><strong>de:</strong> Description/definition of enzyme function</li>
 *   <li><strong>ca:</strong> Catalytic activity description</li>
 *   <li><strong>an:</strong> Alternative names for the enzyme</li>
 * </ul>
 * 
 * <h2>Field Descriptions</h2>
 * 
 * <h3>Catalytic Activity (CA):</h3>
 * <p>Describes the chemical reaction catalyzed by the enzyme:</p>
 * <ul>
 *   <li>Substrate → Product format</li>
 *   <li>Cofactor requirements</li>
 *   <li>Reaction conditions</li>
 * </ul>
 * 
 * <h3>Alternative Names (AN):</h3>
 * <p>Common names and synonyms for the enzyme:</p>
 * <ul>
 *   <li>Systematic names</li>
 *   <li>Common names</li>
 *   <li>Obsolete names</li>
 *   <li>Abbreviations</li>
 * </ul>
 * 
 * <h2>Usage in Metaproteomics</h2>
 * 
 * <h3>Functional Classification:</h3>
 * <ul>
 *   <li><strong>Protein Annotation:</strong> Assign EC numbers to identified proteins</li>
 *   <li><strong>Pathway Analysis:</strong> Map enzymes to metabolic pathways</li>
 *   <li><strong>Comparative Analysis:</strong> Compare enzymatic activities across samples</li>
 *   <li><strong>Systems Biology:</strong> Understand metabolic capabilities</li>
 * </ul>
 * 
 * <h3>Database Integration:</h3>
 * <ul>
 *   <li><strong>KEGG:</strong> Kyoto Encyclopedia of Genes and Genomes</li>
 *   <li><strong>UniProt:</strong> Universal Protein Resource</li>
 *   <li><strong>BRENDA:</strong> Comprehensive enzyme database</li>
 *   <li><strong>MetaCyc:</strong> Metabolic pathway database</li>
 * </ul>
 * 
 * <h2>Examples</h2>
 * 
 * <h3>Common Enzymes:</h3>
 * <ul>
 *   <li><strong>EC 1.1.1.1:</strong> Alcohol dehydrogenase</li>
 *   <li><strong>EC 2.7.1.1:</strong> Hexokinase</li>
 *   <li><strong>EC 3.2.1.1:</strong> α-amylase</li>
 *   <li><strong>EC 4.1.1.1:</strong> Pyruvate decarboxylase</li>
 *   <li><strong>EC 5.3.1.1:</strong> Triose-phosphate isomerase</li>
 *   <li><strong>EC 6.2.1.1:</strong> Acetyl-CoA synthetase</li>
 * </ul>
 * 
 * <h3>Usage Example:</h3>
 * <pre>{@code
 *   // Create enzyme entry
 *   EnzymeCommission ec = new EnzymeCommission("1.1.1.1");
 *   ec.setDe("Alcohol dehydrogenase");
 *   ec.setCa("Primary alcohol + NAD+ = aldehyde + NADH + H+");
 *   ec.addAn("ADH");
 *   ec.addAn("alcohol:NAD+ oxidoreductase");
 *   
 *   // Access information
 *   String id = ec.getId();          // "1.1.1.1"
 *   String names = ec.getAnString(); // "ADH;alcohol:NAD+ oxidoreductase"
 * }</pre>
 * 
 * <h2>Data Sources</h2>
 * <ul>
 *   <li><strong>IUBMB:</strong> International Union of Biochemistry and Molecular Biology</li>
 *   <li><strong>ENZYME:</strong> ENZYME database at ExPASy</li>
 *   <li><strong>UniProt:</strong> Cross-referenced EC numbers</li>
 *   <li><strong>KEGG:</strong> KEGG ENZYME database</li>
 * </ul>
 * 
 * <h2>Quality Considerations</h2>
 * <ul>
 *   <li><strong>Specificity:</strong> More specific EC numbers provide better annotation</li>
 *   <li><strong>Evidence:</strong> Experimental vs computational predictions</li>
 *   <li><strong>Coverage:</strong> Not all proteins have assigned EC numbers</li>
 *   <li><strong>Updates:</strong> EC classification evolves with new discoveries</li>
 * </ul>
 * 
 * <h2>Integration with Other Classifications</h2>
 * <ul>
 *   <li><strong>COG:</strong> Clusters of Orthologous Groups</li>
 *   <li><strong>PFAM:</strong> Protein families</li>
 *   <li><strong>GO:</strong> Gene Ontology molecular function</li>
 *   <li><strong>InterPro:</strong> Integrated protein signature database</li>
 * </ul>
 * 
 * <h2>Thread Safety</h2>
 * <p><strong>Note:</strong> This class is NOT thread-safe. External synchronization
 * is required if accessed from multiple threads concurrently.</p>
 * 
 * @author Kai Cheng - Quadram Institute Bioscience
 * @version 2.0
 * @since 1.0
 * 
 * @see uk.ac.quadram.metapilot.core.function
 * @see uk.ac.quadram.metapilot.task.MetaAbstractTask
 */
public class EnzymeCommission {

	private String id;
	private String de;
	private String ca;
	private ArrayList<String> an;

	public EnzymeCommission(String id) {
		super();
		this.id = id;
		this.ca = "";
		this.an = new ArrayList<String>();
	}

	public EnzymeCommission(String id, String de, String ca, ArrayList<String> an) {
		super();
		this.id = id;
		this.de = de;
		this.ca = ca;
		this.an = an;
	}

	public String getId() {
		return id;
	}

	public void setId(String id) {
		this.id = id;
	}

	public String getDe() {
		return de;
	}

	public void setDe(String de) {
		this.de = de;
	}

	public String getCa() {
		return ca;
	}

	public void setCa(String ca) {
		this.ca = ca;
	}

	public ArrayList<String> getAn() {
		return an;
	}
	
	public String getAnString() {
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < an.size(); i++) {
			sb.append(an.get(i)).append(";");
		}
		if (sb.length() > 0) {
			sb.deleteCharAt(sb.length() - 1);
		}
		return sb.toString();
	}

	public void setAn(ArrayList<String> an) {
		this.an = an;
	}

	public void addCa(String caString) {
		this.ca += caString;
	}

	public void addAn(String anString) {
		this.an.add(anString);
	}
}
