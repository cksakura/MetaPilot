/**
 * 
 */
package uk.ac.quadram.metapilot.core.taxonomy;

import java.util.HashMap;

/**
 * Comprehensive Taxonomic Hierarchy Classification System for Metaproteomics
 * 
 * <h2>Purpose</h2>
 * This enumeration defines the complete taxonomic hierarchy used in biological
 * classification systems, specifically designed for metaproteomics applications.
 * It provides standardized taxonomic ranks essential for species identification,
 * phylogenetic analysis, and microbiome community structure analysis.
 * 
 * <h2>Taxonomic Hierarchy</h2>
 * <p>The taxonomic classification follows the standard biological hierarchy
 * established by Carl Linnaeus and refined by modern taxonomy:</p>
 * 
 * <h3>Complete Hierarchy (Top to Bottom):</h3>
 * <pre>{@code
 *   Root
 *   ├── Domain (Superkingdom)    [Main Rank 0]
 *   │   ├── Kingdom              [Main Rank 1]
 *   │   │   ├── Subkingdom
 *   │   │   ├── Superphylum
 *   │   │   ├── Phylum           [Main Rank 2]
 *   │   │   │   ├── Subphylum
 *   │   │   │   ├── Superclass
 *   │   │   │   ├── Class        [Main Rank 3]
 *   │   │   │   │   ├── Subclass
 *   │   │   │   │   ├── Infraclass
 *   │   │   │   │   ├── Superorder
 *   │   │   │   │   ├── Order    [Main Rank 4]
 *   │   │   │   │   │   ├── Suborder
 *   │   │   │   │   │   ├── Infraorder
 *   │   │   │   │   │   ├── Parvorder
 *   │   │   │   │   │   ├── Superfamily
 *   │   │   │   │   │   ├── Family [Main Rank 5]
 *   │   │   │   │   │   │   ├── Subfamily
 *   │   │   │   │   │   │   ├── Tribe
 *   │   │   │   │   │   │   ├── Subtribe
 *   │   │   │   │   │   │   ├── Genus   [Main Rank 6]
 *   │   │   │   │   │   │   │   ├── Subgenus
 *   │   │   │   │   │   │   │   ├── Species group
 *   │   │   │   │   │   │   │   ├── Species [Main Rank 7]
 *   │   │   │   │   │   │   │   │   ├── Species subgroup
 *   │   │   │   │   │   │   │   │   ├── Subspecies
 *   │   │   │   │   │   │   │   │   ├── Varietas
 *   │   │   │   │   │   │   │   │   └── Forma
 * }</pre>
 * 
 * <h2>Main Ranks vs Intermediate Ranks</h2>
 * 
 * <h3>Main Ranks (8 core levels):</h3>
 * <p>These are the primary taxonomic levels used in most classifications:</p>
 * <ul>
 *   <li><strong>Domain:</strong> Bacteria, Archaea, Eukaryota</li>
 *   <li><strong>Kingdom:</strong> Traditional kingdoms (Animalia, Plantae, etc.)</li>
 *   <li><strong>Phylum:</strong> Major body plan divisions</li>
 *   <li><strong>Class:</strong> Fundamental groupings within phyla</li>
 *   <li><strong>Order:</strong> Related families</li>
 *   <li><strong>Family:</strong> Related genera</li>
 *   <li><strong>Genus:</strong> Closely related species</li>
 *   <li><strong>Species:</strong> Individual species</li>
 * </ul>
 * 
 * <h3>Simplified 7-Rank System:</h3>
 * <p>For streamlined analysis, excludes Kingdom level:</p>
 * <ul>
 *   <li>Domain → Phylum → Class → Order → Family → Genus → Species</li>
 * </ul>
 * 
 * <h3>Intermediate Ranks:</h3>
 * <p>Additional subdivisions for complex taxonomies:</p>
 * <ul>
 *   <li><strong>Sub-ranks:</strong> Subkingdom, Subphylum, Subclass, etc.</li>
 *   <li><strong>Super-ranks:</strong> Superphylum, Superclass, Superorder, etc.</li>
 *   <li><strong>Infra-ranks:</strong> Infraclass, Infraorder</li>
 *   <li><strong>Special ranks:</strong> Tribe, Parvorder, Varietas, Forma</li>
 * </ul>
 * 
 * <h2>Metaproteomics Applications</h2>
 * 
 * <h3>Taxonomic Assignment:</h3>
 * <ul>
 *   <li><strong>Protein Classification:</strong> Assign proteins to taxonomic groups</li>
 *   <li><strong>LCA Analysis:</strong> Lowest Common Ancestor determination</li>
 *   <li><strong>Phylogenetic Profiling:</strong> Community composition analysis</li>
 *   <li><strong>Differential Analysis:</strong> Compare taxonomic abundances</li>
 * </ul>
 * 
 * <h3>Database Integration:</h3>
 * <ul>
 *   <li><strong>NCBI Taxonomy:</strong> Primary taxonomic database</li>
 *   <li><strong>GTDB:</strong> Genome Taxonomy Database</li>
 *   <li><strong>Silva:</strong> rRNA sequence database</li>
 *   <li><strong>UniProt:</strong> Protein taxonomic annotations</li>
 * </ul>
 * 
 * <h2>ID System</h2>
 * <h3>Rank IDs:</h3>
 * <ul>
 *   <li><strong>Root:</strong> -1 (top level)</li>
 *   <li><strong>Superkingdom:</strong> 0 to Species: 23</li>
 *   <li><strong>Subspecific:</strong> 24-27 (subspecies, varietas, forma)</li>
 * </ul>
 * 
 * <h3>Main Rank IDs:</h3>
 * <ul>
 *   <li>Main ranks have corresponding mainId values (0-7)</li>
 *   <li>Intermediate ranks have mainId = -1</li>
 * </ul>
 * 
 * <h2>Usage Examples</h2>
 * 
 * <h3>Basic Usage:</h3>
 * <pre>{@code
 *   // Get main taxonomic ranks
 *   TaxonomyRanks[] mainRanks = TaxonomyRanks.getMainRanks();
 *   
 *   // Check if rank is a main rank
 *   boolean isMain = TaxonomyRanks.Genus.isMainRank(); // true
 *   boolean isMain2 = TaxonomyRanks.Subgenus.isMainRank(); // false
 *   
 *   // Get rank from name
 *   TaxonomyRanks rank = TaxonomyRanks.getRankFromName("species");
 * }</pre>
 * 
 * <h3>Rank Mapping:</h3>
 * <pre>{@code
 *   // Create rank to ID mapping
 *   HashMap<String, Integer> rankMap = TaxonomyRanks.getRankIdMap();
 *   
 *   // Convert between ID and rank
 *   TaxonomyRanks genus = TaxonomyRanks.getRankFromId(20);
 *   String rankName = TaxonomyRanks.getRankNameFromId(23); // "species"
 * }</pre>
 * 
 * <h3>Hierarchical Analysis:</h3>
 * <pre>{@code
 *   // Analyze taxonomic resolution
 *   TaxonomyRanks[] ranks = TaxonomyRanks.getMainRanks7();
 *   for (TaxonomyRanks rank : ranks) {
 *       System.out.println(rank.getName() + " (ID: " + rank.getId() + ")");
 *   }
 * }</pre>
 * 
 * <h2>Biological Examples</h2>
 * 
 * <h3>Human Classification:</h3>
 * <ul>
 *   <li><strong>Domain:</strong> Eukaryota</li>
 *   <li><strong>Kingdom:</strong> Animalia</li>
 *   <li><strong>Phylum:</strong> Chordata</li>
 *   <li><strong>Class:</strong> Mammalia</li>
 *   <li><strong>Order:</strong> Primates</li>
 *   <li><strong>Family:</strong> Hominidae</li>
 *   <li><strong>Genus:</strong> Homo</li>
 *   <li><strong>Species:</strong> Homo sapiens</li>
 * </ul>
 * 
 * <h3>E. coli Classification:</h3>
 * <ul>
 *   <li><strong>Domain:</strong> Bacteria</li>
 *   <li><strong>Phylum:</strong> Proteobacteria</li>
 *   <li><strong>Class:</strong> Gammaproteobacteria</li>
 *   <li><strong>Order:</strong> Enterobacterales</li>
 *   <li><strong>Family:</strong> Enterobacteriaceae</li>
 *   <li><strong>Genus:</strong> Escherichia</li>
 *   <li><strong>Species:</strong> Escherichia coli</li>
 * </ul>
 * 
 * <h2>Quality Considerations</h2>
 * <ul>
 *   <li><strong>Resolution:</strong> Higher ranks more reliable than species-level</li>
 *   <li><strong>Database Coverage:</strong> Some organisms poorly represented</li>
 *   <li><strong>Phylogenetic Updates:</strong> Taxonomy evolves with new research</li>
 *   <li><strong>Assignment Confidence:</strong> Use appropriate confidence thresholds</li>
 * </ul>
 * 
 * <h2>Performance Considerations</h2>
 * <ul>
 *   <li><strong>Memory Efficiency:</strong> Enum provides memory-efficient storage</li>
 *   <li><strong>Lookup Speed:</strong> O(1) access for ID-based lookups</li>
 *   <li><strong>Mapping Operations:</strong> Pre-computed maps for frequent operations</li>
 * </ul>
 * 
 * <h2>Thread Safety</h2>
 * <p>This enum is thread-safe as all methods are stateless and operate on
 * immutable enum constants.</p>
 * 
 * @author Kai Cheng - Quadram Institute Bioscience
 * @version 2.0
 * @since 1.0
 * 
 * @see uk.ac.quadram.metapilot.core.taxonomy.TaxonomyDatabase
 * @see uk.ac.quadram.metapilot.core.taxonomy.TaxonomyParameter
 */
public enum TaxonomyRanks {

	Root(-1, "root", true, -1),

	Superkingdom(0, "domain", true, 0),

	Kingdom(1, "kingdom", true, 1),

	Subkingdom(2, "subkingdom", false, -1),

	Superphylum(3, "superphylum", false, -1),

	Phylum(4, "phylum", true, 2),

	Subphylum(5, "subphylum", false, -1),

	Superclass(6, "superclass", false, -1),

	Class(7, "class", true, 3),

	Subclass(8, "subclass", false, -1),

	Infraclass(9, "infraclass", false, -1),

	Superorder(10, "superorder", false, -1),

	Order(11, "order", true, 4),

	Suborder(12, "suborder", false, -1),

	Infraorder(13, "infraorder", false, -1),

	Parvorder(14, "parvorder", true, -1),

	Superfamily(15, "superfamily", false, -1),

	Family(16, "family", true, 5),

	Subfamily(17, "subfamily", false, -1),

	Tribe(18, "tribe", false, -1),

	Subtribe(19, "subtribe", false, -1),

	Genus(20, "genus", true, 6),

	Subgenus(21, "subgenus", false, -1),

	Species_group(22, "species group", false, -1),

	Species(23, "species", true, 7),

	Species_subgroup(24, "species subgroup", false, -1),

	Subspecies(25, "subspecies", false, -1),

	Varietas(26, "varietas", false, -1),

	Forma(27, "forma", false, -1),

	;

	protected static final String[] ranks = new String[] { "domain", "kingdom", "subkingdom", "superphylum",
			"phylum", "subphylum", "superclass", "class", "subclass", "infraclass", "superorder", "order", "suborder",
			"infraorder", "parvorder", "superfamily", "family", "subfamily", "tribe", "subtribe", "genus", "subgenus",
			"species group", "species", "species subgroup", "subspecies", "varietas", "forma" };

	private int id;
	private String name;
	private boolean mainRank;
	private int mainId;

	TaxonomyRanks(int id, String name, boolean mainRank, int mainId) {
		this.id = id;
		this.name = name;
		this.mainRank = mainRank;
		this.mainId = mainId;
	}

	public int getId() {
		return id;
	}

	public String getName() {
		return name;
	}

	public boolean isMainRank() {
		return mainRank;
	}

	public int getMainId() {
		return mainId;
	}

	public static HashMap<String, Integer> getRankIdMap() {
		HashMap<String, Integer> rankIdMap = new HashMap<String, Integer>();
		for (TaxonomyRanks tr : TaxonomyRanks.values()) {
			rankIdMap.put(tr.getName(), tr.getId());
		}
		return rankIdMap;
	}

	public static TaxonomyRanks getRankFromId(int id) {
		TaxonomyRanks[] values = TaxonomyRanks.values();
		if (id >= -1 && id <= 28) {
			return values[id + 1];
		} else {
			return null;
		}
	}

	public static TaxonomyRanks getRankFromName(String name) {
		TaxonomyRanks[] values = TaxonomyRanks.values();
		for (TaxonomyRanks tr : values) {
			if (tr.getName().equalsIgnoreCase(name)) {
				return tr;
			}
		}
		return null;
	}

	public static String getRankNameFromId(int id) {
		TaxonomyRanks[] values = TaxonomyRanks.values();
		if (id >= -1 && id <= 28) {
			return values[id + 1].getName();
		} else {
			return null;
		}
	}

	public static TaxonomyRanks[] getMainRanks() {
		TaxonomyRanks[] values = new TaxonomyRanks[] { Superkingdom, Kingdom, Phylum, Class, Order, Family, Genus,
				Species };
		return values;
	}
	
	public static TaxonomyRanks[] getMainRanks7() {
		TaxonomyRanks[] values = new TaxonomyRanks[] { Superkingdom, Phylum, Class, Order, Family, Genus, Species };
		return values;
	}

	public static int[] getMainRankIds() {
		TaxonomyRanks[] mainRanks = getMainRanks();
		int[] mainRankIds = new int[mainRanks.length];
		for (int i = 0; i < mainRankIds.length; i++) {
			mainRankIds[i] = mainRanks[i].getId();
		}
		return mainRankIds;
	}
}
