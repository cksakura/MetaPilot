package uk.ac.quadram.metapilot.core.taxonomy;

import java.io.IOException;
import java.sql.*;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * SQL-based implementation of TaxonomyDatabase that extends the original class
 * to provide the same interface but with SQL backend and smart caching
 */
public class TaxonomyDatabaseSQL extends TaxonomyDatabase {

	private static Logger LOGGER = LogManager.getLogger(TaxonomyDatabaseSQL.class);

	// SQL connection and prepared statements
	private Connection connection;
	private PreparedStatement getTaxonStmt;
	private PreparedStatement getTaxonByNameStmt;
	private PreparedStatement getTaxonByRankStmt;

	// Smart caching strategy
	private final Map<Integer, Taxon> hierarchyCache; // Always cached (main ranks)
	private final Map<Integer, Taxon> lruCache; // LRU cache for everything else
	private final Map<String, Integer> nameCache; // Name to ID mapping cache
	private final int MAX_LRU_SIZE = 100000;

	// Cache for taxonNameMap functionality
	private final Map<String, Integer> taxonNameMapCache;

	public TaxonomyDatabaseSQL(String connectionString) throws SQLException {
		// The parent constructor will be called automatically and will skip file reading
		// because of the instanceof check in the read() method
		
		// Initialize caches
		this.hierarchyCache = new ConcurrentHashMap<>();
		this.lruCache = new LinkedHashMap<Integer, Taxon>(16, 0.75f, true) {
			@Override
			protected boolean removeEldestEntry(Map.Entry<Integer, Taxon> eldest) {
				return size() > MAX_LRU_SIZE;
			}
		};
		this.nameCache = new ConcurrentHashMap<>();
		this.taxonNameMapCache = new ConcurrentHashMap<>();

		// Initialize SQL connection
		initializeDatabase(connectionString);

		// Build hierarchy cache (most important taxa)
		buildHierarchyCache();
	}

	private void initializeDatabase(String connectionString) throws SQLException {
		// Ensure the connection string has the proper jdbc:sqlite: prefix
		String jdbcUrl;
		if (connectionString.startsWith("jdbc:sqlite:")) {
			jdbcUrl = connectionString;
		} else {
			jdbcUrl = "jdbc:sqlite:" + connectionString;
		}
		
		LOGGER.info("Connecting to database: " + jdbcUrl);
		this.connection = DriverManager.getConnection(jdbcUrl);

		// Prepare statements
		this.getTaxonStmt = connection.prepareStatement(
				"SELECT taxon_id, scientific_name, rank, lineage, parent FROM taxonomy WHERE taxon_id = ?");

		this.getTaxonByNameStmt = connection.prepareStatement(
				"SELECT taxon_id, scientific_name, rank, lineage, parent FROM taxonomy WHERE scientific_name = ? AND rank = ?");

		this.getTaxonByRankStmt = connection.prepareStatement(
				"SELECT taxon_id, scientific_name, rank, lineage, parent FROM taxonomy WHERE rank = ?");
	}

	private void buildHierarchyCache() {
		try {
			LOGGER.info("Loading all taxonomy data into RAM for better performance and reliability...");
			
			// Load ALL taxonomy data into RAM at startup
			// This eliminates database query issues and invalid parent references
			try (Statement stmt = connection.createStatement()) {
				ResultSet rs = stmt.executeQuery("SELECT * FROM taxonomy ORDER BY taxon_id");
				
				int totalLoaded = 0;
				int invalidParentCount = 0;
				
				while (rs.next()) {
					Taxon taxon = createTaxonFromResultSet(rs);
					hierarchyCache.put(taxon.getId(), taxon);
					nameCache.put(taxon.getName(), taxon.getId());
					taxonNameMapCache.put(taxon.getName() + "_" + taxon.getRankId(), taxon.getId());
					totalLoaded++;
					
					// Track invalid parent references for debugging
					// Note: We have multiple root nodes: 10239 (viruses), 131567 (cellular organisms), unclassified, other entries
					if (taxon.getParentId() == -1) { // -1 means no parent (root node)
						LOGGER.debug("Root taxon: " + taxon.getId() + " (" + taxon.getName() + ")");
					} else if (taxon.getParentId() == 0) { // 0 indicates a data corruption issue
						invalidParentCount++;
						if (invalidParentCount <= 10) { // Log first 10 examples
							LOGGER.debug("Taxon " + taxon.getId() + " (" + taxon.getName() + ") has corrupted parent ID: 0 (should be -1 for no parent or valid parent ID)");
						}
					}
				}
				
				LOGGER.info("Loaded " + totalLoaded + " taxa into RAM");
				if (invalidParentCount > 0) {
					LOGGER.warn("Found " + invalidParentCount + " taxa with invalid parent IDs (parent <= 0)");
				}
			}

			LOGGER.info("Total taxonomy cache size: " + hierarchyCache.size() + " taxa");
			
			// Identify and log all root nodes
			identifyRootNodes();
			
			// Check what rank names are used in the database
			checkRankNamesInDatabase();
			
			// Validate taxonomy data integrity
			validateTaxonomyData();

		} catch (SQLException e) {
			LOGGER.error("Error building hierarchy cache", e);
		}
	}
	
	/**
	 * Identify and log all root nodes in the taxonomy
	 */
	private void identifyRootNodes() {
		LOGGER.info("Identifying root nodes in taxonomy...");
		
		int rootCount = 0;
		for (Taxon taxon : hierarchyCache.values()) {
			if (taxon.getParentId() == -1) {
				rootCount++;
				LOGGER.info("  Root " + rootCount + ": ID=" + taxon.getId() + ", Name=" + taxon.getName());
			}
		}
		
		LOGGER.info("Found " + rootCount + " root nodes in taxonomy");
	}
	
	/**
	 * Map unrecognized rank names to appropriate rank IDs
	 */
	private int mapUnrecognizedRank(String rankStr) {
		if (rankStr == null) return -1;
		
		String rank = rankStr.toLowerCase().trim();
		
		// Map common unrecognized ranks to known ranks
		switch (rank) {
		case "clade":
			// Clade is a general term - try to map based on context
			// For now, treat as intermediate rank (continue traversal)
			return -1;
		case "no rank":
		case "unclassified":
			// These should continue traversal but not be stored
			return -1;
		default:
			// For other unrecognized ranks, return -1 (intermediate rank)
			return -1;
		}
	}
	
	/**
	 * Check what rank names are actually used in the database
	 */
	private void checkRankNamesInDatabase() {
		LOGGER.info("Checking rank names used in database...");
		
		// Check what rank names exist in the database
		HashMap<String, Integer> rankNameCounts = new HashMap<String, Integer>();
		for (Taxon taxon : hierarchyCache.values()) {
			String rankName = taxon.getRank();
			if (rankName != null) {
				rankNameCounts.put(rankName, rankNameCounts.getOrDefault(rankName, 0) + 1);
			}
		}
		
		// Log the most common rank names
		LOGGER.info("Most common rank names in database:");
		rankNameCounts.entrySet().stream()
			.sorted((e1, e2) -> e2.getValue().compareTo(e1.getValue()))
			.limit(15)
			.forEach(entry -> LOGGER.info("  " + entry.getKey() + ": " + entry.getValue() + " taxa"));
		
		// Check specifically for domain vs superkingdom
		if (rankNameCounts.containsKey("domain")) {
			LOGGER.info("Database uses 'domain' instead of 'superkingdom' - " + rankNameCounts.get("domain") + " taxa");
		}
		if (rankNameCounts.containsKey("superkingdom")) {
			LOGGER.info("Database uses 'superkingdom' - " + rankNameCounts.get("superkingdom") + " taxa");
		}
		
		// Check for unrecognized ranks
		HashMap<String, Integer> unrecognizedRanks = new HashMap<String, Integer>();
		for (String rankName : rankNameCounts.keySet()) {
			if (!TaxonomyRanks.getRankIdMap().containsKey(rankName)) {
				unrecognizedRanks.put(rankName, rankNameCounts.get(rankName));
			}
		}
		
		if (!unrecognizedRanks.isEmpty()) {
			LOGGER.warn("Found " + unrecognizedRanks.size() + " unrecognized rank names:");
			unrecognizedRanks.entrySet().stream()
				.sorted((e1, e2) -> e2.getValue().compareTo(e1.getValue()))
				.forEach(entry -> LOGGER.warn("  " + entry.getKey() + ": " + entry.getValue() + " taxa"));
		}
	}
	
	/**
	 * Validate taxonomy data integrity and report issues
	 */
	private void validateTaxonomyData() {
		LOGGER.info("Validating taxonomy data integrity...");
		
		int totalTaxa = hierarchyCache.size();
		int taxaWithInvalidParents = 0;
		int taxaWithValidLineage = 0;
		int taxaWithBrokenLineage = 0;
		
		// Sample some common taxa for validation
		int[] testTaxa = {9606, 10090, 10116, 3702, 559292}; // Human, Mouse, Rat, Arabidopsis, S. cerevisiae
		
		for (int testTaxonId : testTaxa) {
			Taxon taxon = hierarchyCache.get(testTaxonId);
			if (taxon != null) {
				LOGGER.info("Testing lineage for taxon " + testTaxonId + " (" + taxon.getName() + ")");
				LOGGER.info("  Taxon rank ID: " + taxon.getRankId() + ", Rank name: " + taxon.getRank() + ", Parent ID: " + taxon.getParentId());
				int[] lineage = getMainParentTaxonIds(taxon);
				LOGGER.info("  Lineage array: " + java.util.Arrays.toString(lineage));
				
				int validIds = 0;
				int missingIds = 0;
				for (int id : lineage) {
					if (id > 0) {
						validIds++;
					} else if (id == -1) {
						missingIds++;
					}
				}
				if (validIds > 0) {
					taxaWithValidLineage++;
					LOGGER.info("  Result: " + validIds + " valid lineage entries, " + missingIds + " missing ranks");
				} else {
					taxaWithBrokenLineage++;
					LOGGER.warn("  Result: BROKEN lineage (no valid entries found) - all entries are -1");
				}
			} else {
				LOGGER.warn("Test taxon " + testTaxonId + " not found in cache");
			}
		}
		
		// Count taxa with invalid parents
		for (Taxon taxon : hierarchyCache.values()) {
			if (taxon.getParentId() == 0) { // 0 indicates data corruption
				taxaWithInvalidParents++;
			}
			// Note: -1 is valid for multiple root nodes (viruses, cellular organisms, unclassified, other)
		}
		
		LOGGER.info("Taxonomy validation results:");
		LOGGER.info("  Total taxa: " + totalTaxa);
		LOGGER.info("  Taxa with invalid parents: " + taxaWithInvalidParents);
		LOGGER.info("  Test taxa with valid lineage: " + taxaWithValidLineage);
		LOGGER.info("  Test taxa with broken lineage: " + taxaWithBrokenLineage);
		
		if (taxaWithBrokenLineage > 0) {
			LOGGER.warn("Some taxa have broken lineages - this may cause issues in taxonomic analysis");
		}
	}

	private Taxon createTaxonFromResultSet(ResultSet rs) throws SQLException {
		int taxonId = rs.getInt("taxon_id");
		String name = rs.getString("scientific_name");
		String rankStr = rs.getString("rank");
		String lineage = rs.getString("lineage");
		
		// Handle NULL parent values properly - use getObject() to detect NULL
		Object parentObj = rs.getObject("parent");
		int parent = -1; // Default for no parent
		if (parentObj != null) {
			parent = ((Number) parentObj).intValue();
		}

		// Parse rank
		int rank = -1;
		if (rankStr != null && !rankStr.trim().isEmpty()) {
			if (TaxonomyRanks.getRankIdMap().containsKey(rankStr)) {
				rank = TaxonomyRanks.getRankIdMap().get(rankStr);
			} else {
				// Handle unrecognized rank names by mapping them to appropriate IDs
				rank = mapUnrecognizedRank(rankStr);
				if (rank == -1) {
					// Log unrecognized ranks for debugging
					if (taxonId == 9606 || taxonId == 10090 || taxonId == 10116 || taxonId == 3702 || taxonId == 559292) {
						LOGGER.debug("Unrecognized rank '" + rankStr + "' for taxon " + taxonId + " (" + name + ")");
					}
				}
			}
		}

		// Parse root type from lineage
		RootType rootType = null;
		if (lineage != null && lineage.trim().length() > 0) {
			String[] lineageParts = lineage.split("[; ]");
			if (lineageParts.length > 0) {
				String rootTypeStr = lineageParts[0].trim().replace(" ", "_");
				try {
					rootType = RootType.valueOf(rootTypeStr);
				} catch (IllegalArgumentException e) {
					rootType = RootType.getRootType(taxonId);
				}
			}
		} else {
			rootType = RootType.getRootType(taxonId);
		}

		return new Taxon(taxonId, parent, rank, name, rootType);
	}

	// Override the core data access methods to use SQL + caching

	@Override
	public Taxon getTaxonFromId(int taxonId) {
		// Since all data is now loaded into RAM, check cache directly
		if (hierarchyCache.containsKey(taxonId)) {
			LOGGER.debug("Found taxon ID " + taxonId + " in RAM cache");
			return hierarchyCache.get(taxonId);
		}

		// Taxon not found in cache (shouldn't happen if database is complete)
		LOGGER.warn("No taxon found in RAM cache for ID: " + taxonId + " (database may be incomplete)");
		return null;
	}

	@Override
	public Taxon getTaxonFromName(String name_rankId) {
		// Since all data is now loaded into RAM, check cache directly
		if (taxonNameMapCache.containsKey(name_rankId)) {
			return getTaxonFromId(taxonNameMapCache.get(name_rankId));
		}

		// Parse name and rank from the key
		String[] parts = name_rankId.split("_");
		if (parts.length >= 2) {
			String name = parts[0];
			int rankId = Integer.parseInt(parts[1]);
			return getTaxonFromName(name, rankId);
		}

		return null;
	}

	@Override
	public Taxon getTaxonFromName(String name, int rankId) {
		// Since all data is now loaded into RAM, check cache directly
		String key = name + "_" + rankId;
		if (taxonNameMapCache.containsKey(key)) {
			return getTaxonFromId(taxonNameMapCache.get(key));
		}

		// Taxon not found in cache (shouldn't happen if database is complete)
		LOGGER.debug("No taxon found in RAM cache for name: " + name + " with rank ID: " + rankId);
		return null;
	}

	// Override HashMap access methods to prevent NullPointerException

	@Override
	public HashMap<Integer, Taxon> getTaxonMap() {
		// Return a view that combines all caches
		HashMap<Integer, Taxon> combinedMap = new HashMap<>();
		combinedMap.putAll(hierarchyCache);
		synchronized (lruCache) {
			combinedMap.putAll(lruCache);
		}
		return combinedMap;
	}

	@Override
	public HashMap<String, Integer> getTaxonNameMap() {
		// Return a copy of the name cache
		return new HashMap<>(taxonNameMapCache);
	}

	// Override methods that directly access taxonMap to use our getTaxonFromId

	@Override
	public int getMainParentId(int taxonId) {
		Taxon taxon = getTaxonFromId(taxonId);
		if (taxon == null)
			return -1;

		int parentId = taxon.getParentId();
		if (parentId == -1)
			return -1;

		Taxon parent = getTaxonFromId(parentId);
		if (parent == null)
			return -1;

		if (TaxonomyRanks.getRankFromId(parent.getRankId()).isMainRank()) {
			return parent.getId();
		} else {
			return getMainParentId(parent.getId());
		}
	}

	@Override
	public int getMainRankId(int taxonId) {
		Taxon taxon = getTaxonFromId(taxonId);
		if (taxon == null)
			return -1;

		int rankId = taxon.getRankId();
		if (taxon.isMainRank()) {
			return rankId;
		} else {
			return getMainRankId(taxon.getParentId());
		}
	}

	@Override
	public int getTaxonId4Rank(int taxonId, int rankId) {
		Taxon taxon = getTaxonFromId(taxonId);
		if (taxon == null)
			return -1;

		if (taxon.getRankId() == rankId) {
			return taxonId;
		} else {
			return getTaxonId4Rank(taxon.getParentId(), rankId);
		}
	}

	@Override
	public Taxon getMainParent(int taxonId) {
		return getMainParent(getTaxonFromId(taxonId));
	}

	@Override
	public Taxon getMainParent(Taxon taxon) {
		if (taxon == null)
			return null;
		if (taxon.getId() == 131567)
			return taxon;
		if (taxon.isMainRank())
			return taxon;
		else
			return getMainParent(getTaxonFromId(taxon.getParentId()));
	}

	@Override
	public int[] getMainParentTaxonIds(int taxonId) {
		Taxon taxon = getTaxonFromId(taxonId);
		if (taxon == null)
			return null;
		else
			return getMainParentTaxonIds(taxon);
	}

	@Override
	public int[] getMainParentTaxonIds(Taxon taxon) {
		int[] taxonomy = new int[8];
		Arrays.fill(taxonomy, -1);

		int traversalDepth = 0;
		while (taxon != null && traversalDepth < 50) { // Prevent infinite loops
			int rankid = taxon.getRankId();
			LOGGER.debug("Traversal depth " + traversalDepth + ": Processing taxon " + taxon.getId() + 
						" (" + taxon.getName() + ") with rank ID " + rankid);
			
			// Map rank IDs to array positions based on TaxonomyRanks enum
			// Main ranks: Domain(0), Kingdom(1), Phylum(4), Class(7), Order(11), Family(16), Genus(20), Species(23)
			switch (rankid) {
			case -1: // Root or unrecognized rank - continue traversal
				LOGGER.debug("Root or unrecognized rank " + rankid + " for taxon " + taxon.getId() + ", continuing traversal");
				break;
			case 0: // Domain
				taxonomy[0] = taxon.getId();
				break;
			case 1: // Kingdom
				taxonomy[1] = taxon.getId();
				break;
			case 4: // Phylum
				taxonomy[2] = taxon.getId();
				break;
			case 7: // Class
				taxonomy[3] = taxon.getId();
				break;
			case 11: // Order
				taxonomy[4] = taxon.getId();
				break;
			case 16: // Family
				taxonomy[5] = taxon.getId();
				break;
			case 20: // Genus
				taxonomy[6] = taxon.getId();
				break;
			case 23: // Species
				taxonomy[7] = taxon.getId();
				break;
			default:
				// Intermediate ranks - continue traversal but don't store
				LOGGER.debug("Intermediate rank " + rankid + " for taxon " + taxon.getId() + ", continuing traversal");
				break;
			}
			int parentId = taxon.getParentId();
			LOGGER.debug("Looking up parent taxon ID: " + parentId + " for taxon: " + taxon.getId());
			
			// Handle root nodes (parent ID = -1) and invalid parent IDs (parent ID = 0)
			if (parentId == -1) {
				LOGGER.debug("Reached root taxon " + taxon.getId() + " (" + taxon.getName() + "), stopping traversal");
				break;
			} else if (parentId == 0) {
				LOGGER.debug("Taxon " + taxon.getId() + " has corrupted parent ID 0, stopping traversal");
				break;
			}
			
			// Since all data is now in RAM, get taxon directly from cache
			taxon = hierarchyCache.get(parentId);
			if (taxon == null) {
				LOGGER.warn("Parent taxon ID " + parentId + " not found in cache, stopping hierarchy traversal at depth " + traversalDepth);
				break;
			}
			traversalDepth++;
		}

		if (traversalDepth >= 50) {
			LOGGER.warn("Hierarchy traversal stopped due to depth limit (50 levels)");
		}

		return taxonomy;
	}

	@Override
	public String[] getTaxonArrays(int taxonId) {
		Taxon taxon = getTaxonFromId(taxonId);
		if (taxon == null)
			return null;
		else
			return getTaxonArrays(taxon);
	}

	@Override
	public String[] getTaxonArrays(Taxon taxon) {
		String[] taxonomy = new String[8];
		
		// Get the lineage IDs first
		int[] lineageIds = getMainParentTaxonIds(taxon);
		
		// Convert IDs to names, using empty string for missing ranks (-1)
		for (int i = 0; i < lineageIds.length && i < taxonomy.length; i++) {
			if (lineageIds[i] > 0) {
				Taxon lineageTaxon = getTaxonFromId(lineageIds[i]);
				if (lineageTaxon != null) {
					taxonomy[i] = lineageTaxon.getName();
				} else {
					taxonomy[i] = ""; // Missing taxon
				}
			} else {
				taxonomy[i] = ""; // Missing rank (-1)
			}
		}
		
		return taxonomy;
	}

	// Override getRankNameMap to work with our caching system
	@Override
	public HashMap<String, Integer>[] getRankNameMap() {
		@SuppressWarnings("unchecked")
		HashMap<String, Integer>[] maps = new HashMap[8];
		for (int i = 0; i < maps.length; i++) {
			maps[i] = new HashMap<String, Integer>();
		}

		// Use our combined cache
		HashMap<Integer, Taxon> combinedMap = getTaxonMap();
		for (Integer id : combinedMap.keySet()) {
			Taxon taxon = combinedMap.get(id);
			int rankid = taxon.getRankId();

			switch (rankid) {
			case 0:
				maps[0].put(taxon.getName(), id);
				break;
			case 1:
				maps[1].put(taxon.getName(), id);
				break;
			case 4:
				maps[2].put(taxon.getName(), id);
				break;
			case 7:
				maps[3].put(taxon.getName(), id);
				break;
			case 11:
				maps[4].put(taxon.getName(), id);
				break;
			case 16:
				maps[5].put(taxon.getName(), id);
				break;
			case 20:
				maps[6].put(taxon.getName(), id);
				break;
			case 23:
				maps[7].put(taxon.getName(), id);
				break;
			}
		}

		return maps;
	}

	// Debug method to check database contents
	public void debugDatabaseContents() {
		try {
			LOGGER.info("=== TaxonomyDatabaseSQL Debug Information ===");
			
			// Check if table exists
			try (Statement stmt = connection.createStatement()) {
				ResultSet rs = stmt.executeQuery("SELECT name FROM sqlite_master WHERE type='table' AND name='taxonomy'");
				if (rs.next()) {
					LOGGER.info("Table 'taxonomy' exists in database");
				} else {
					LOGGER.error("Table 'taxonomy' does NOT exist in database");
					return;
				}
			}
			
			// Check table schema
			try (Statement stmt = connection.createStatement()) {
				ResultSet rs = stmt.executeQuery("PRAGMA table_info(taxonomy)");
				LOGGER.info("Table schema:");
				while (rs.next()) {
					LOGGER.info("  Column: " + rs.getString("name") + " Type: " + rs.getString("type"));
				}
			}
			
			// Check total row count
			try (Statement stmt = connection.createStatement()) {
				ResultSet rs = stmt.executeQuery("SELECT COUNT(*) as count FROM taxonomy");
				if (rs.next()) {
					LOGGER.info("Total rows in taxonomy table: " + rs.getInt("count"));
				}
			}
			
			// Check first few rows
			try (Statement stmt = connection.createStatement()) {
				ResultSet rs = stmt.executeQuery("SELECT * FROM taxonomy LIMIT 5");
				LOGGER.info("First 5 rows in taxonomy table:");
				while (rs.next()) {
					LOGGER.info("  ID: " + rs.getInt("taxon_id") + 
							   ", Name: " + rs.getString("scientific_name") + 
							   ", Rank: " + rs.getString("rank") + 
							   ", Parent: " + rs.getInt("parent"));
				}
			}
			
			// Check if taxon ID 1 exists
			try (Statement stmt = connection.createStatement()) {
				ResultSet rs = stmt.executeQuery("SELECT * FROM taxonomy WHERE taxon_id = 1");
				if (rs.next()) {
					LOGGER.info("Taxon ID 1 found: " + rs.getString("scientific_name"));
				} else {
					LOGGER.warn("Taxon ID 1 NOT found in database");
				}
			}
			
			// Check taxon ID 131567 (cellular organisms)
			try (Statement stmt = connection.createStatement()) {
				ResultSet rs = stmt.executeQuery("SELECT * FROM taxonomy WHERE taxon_id = 131567");
				if (rs.next()) {
					LOGGER.info("Taxon ID 131567 (cellular organisms) found: " + rs.getString("scientific_name") + 
							   ", Parent: " + rs.getInt("parent"));
				} else {
					LOGGER.warn("Taxon ID 131567 (cellular organisms) NOT found in database");
				}
			}
			
			// Check for taxons with missing parents (should only be 131567)
			try (Statement stmt = connection.createStatement()) {
				ResultSet rs = stmt.executeQuery("SELECT taxon_id, scientific_name, parent FROM taxonomy WHERE parent IS NULL OR parent = ''");
				LOGGER.info("Taxons with missing parents:");
				int count = 0;
				while (rs.next()) {
					count++;
					LOGGER.info("  ID: " + rs.getInt("taxon_id") + ", Name: " + rs.getString("scientific_name") + 
							   ", Parent: " + rs.getString("parent"));
				}
				LOGGER.info("Total taxons with missing parents: " + count);
			}
			
			// Check for taxons with invalid parent IDs (parent doesn't exist)
			try (Statement stmt = connection.createStatement()) {
				ResultSet rs = stmt.executeQuery(
					"SELECT t1.taxon_id, t1.scientific_name, t1.parent " +
					"FROM taxonomy t1 " +
					"WHERE t1.parent IS NOT NULL AND t1.parent != '' " +
					"AND NOT EXISTS (SELECT 1 FROM taxonomy t2 WHERE t2.taxon_id = t1.parent) " +
					"LIMIT 10"
				);
				LOGGER.info("Sample taxons with invalid parent IDs (parent doesn't exist):");
				int count = 0;
				while (rs.next()) {
					count++;
					LOGGER.info("  ID: " + rs.getInt("taxon_id") + ", Name: " + rs.getString("scientific_name") + 
							   ", Invalid Parent: " + rs.getInt("parent"));
				}
				LOGGER.info("Found " + count + " taxons with invalid parent IDs (showing first 10)");
			}
			
			// Check for taxons with parent ID = 0 (invalid)
			try (Statement stmt = connection.createStatement()) {
				ResultSet rs = stmt.executeQuery(
					"SELECT taxon_id, scientific_name, parent, rank " +
					"FROM taxonomy " +
					"WHERE parent = 0 " +
					"LIMIT 10"
				);
				LOGGER.info("Sample taxons with parent ID = 0 (invalid):");
				int count = 0;
				while (rs.next()) {
					count++;
					LOGGER.info("  ID: " + rs.getInt("taxon_id") + ", Name: " + rs.getString("scientific_name") + 
							   ", Rank: " + rs.getString("rank") + ", Parent: " + rs.getInt("parent"));
				}
				LOGGER.info("Found " + count + " taxons with parent ID = 0 (showing first 10)");
			}
			
			// Check for taxons with parent ID = -1 (invalid)
			try (Statement stmt = connection.createStatement()) {
				ResultSet rs = stmt.executeQuery(
					"SELECT taxon_id, scientific_name, parent, rank " +
					"FROM taxonomy " +
					"WHERE parent = -1 " +
					"LIMIT 10"
				);
				LOGGER.info("Sample taxons with parent ID = -1 (invalid):");
				int count = 0;
				while (rs.next()) {
					count++;
					LOGGER.info("  ID: " + rs.getInt("taxon_id") + ", Name: " + rs.getString("scientific_name") + 
							   ", Rank: " + rs.getString("rank") + ", Parent: " + rs.getInt("parent"));
				}
				LOGGER.info("Found " + count + " taxons with parent ID = -1 (showing first 10)");
			}
			
			// Check parent ID distribution
			try (Statement stmt = connection.createStatement()) {
				ResultSet rs = stmt.executeQuery(
					"SELECT parent, COUNT(*) as count " +
					"FROM taxonomy " +
					"WHERE parent IS NOT NULL " +
					"GROUP BY parent " +
					"ORDER BY count DESC " +
					"LIMIT 20"
				);
				LOGGER.info("Top 20 parent IDs by frequency:");
				while (rs.next()) {
					int parentId = rs.getInt("parent");
					int count = rs.getInt("count");
					LOGGER.info("  Parent ID " + parentId + ": " + count + " children");
				}
			}
			
			// Check minimum and maximum taxon IDs
			try (Statement stmt = connection.createStatement()) {
				ResultSet rs = stmt.executeQuery("SELECT MIN(taxon_id) as min_id, MAX(taxon_id) as max_id FROM taxonomy");
				if (rs.next()) {
					LOGGER.info("Taxon ID range: " + rs.getInt("min_id") + " to " + rs.getInt("max_id"));
				}
			}
			
			LOGGER.info("=== End Debug Information ===");
			
		} catch (SQLException e) {
			LOGGER.error("Error during database debug", e);
		}
	}

	// Clean up resources
	public void close() throws SQLException {
		if (getTaxonStmt != null)
			getTaxonStmt.close();
		if (getTaxonByNameStmt != null)
			getTaxonByNameStmt.close();
		if (getTaxonByRankStmt != null)
			getTaxonByRankStmt.close();
		if (connection != null)
			connection.close();
	}

	/**
	 * Debug method to investigate specific taxon lineage issues
	 */
	public void debugTaxonLineage(int taxonId) {
		LOGGER.info("=== Debugging lineage for taxon " + taxonId + " ===");
		
		Taxon taxon = getTaxonFromId(taxonId);
		if (taxon == null) {
			LOGGER.warn("Taxon " + taxonId + " not found in database");
			return;
		}
		
		LOGGER.info("Taxon " + taxonId + ": " + taxon.getName() + " (rank: " + taxon.getRankId() + ")");
		LOGGER.info("Parent ID: " + taxon.getParentId());
		
		// Trace the lineage step by step
		int depth = 0;
		Taxon currentTaxon = taxon;
		
		while (currentTaxon != null && depth < 20) {
			LOGGER.info("  Depth " + depth + ": ID=" + currentTaxon.getId() + 
					   ", Name=" + currentTaxon.getName() + 
					   ", Rank=" + currentTaxon.getRankId() + 
					   ", Parent=" + currentTaxon.getParentId());
			
			int parentId = currentTaxon.getParentId();
			if (parentId <= 0) {
				LOGGER.warn("    Invalid parent ID " + parentId + " - stopping traversal");
				break;
			}
			
			currentTaxon = getTaxonFromId(parentId);
			if (currentTaxon == null) {
				LOGGER.warn("    Parent taxon " + parentId + " not found in database - stopping traversal");
				break;
			}
			
			depth++;
		}
		
		// Now try the getMainParentTaxonIds method
		LOGGER.info("Testing getMainParentTaxonIds method:");
		int[] lineage = getMainParentTaxonIds(taxon);
		LOGGER.info("Lineage array: " + java.util.Arrays.toString(lineage));
		
		int validCount = 0;
		for (int i = 0; i < lineage.length; i++) {
			if (lineage[i] > 0) {
				validCount++;
				Taxon lineageTaxon = getTaxonFromId(lineage[i]);
				if (lineageTaxon != null) {
					LOGGER.info("  Position " + i + ": ID=" + lineage[i] + ", Name=" + lineageTaxon.getName());
				}
			}
		}
		
		LOGGER.info("Valid lineage entries: " + validCount + " out of " + lineage.length);
		LOGGER.info("=== End lineage debug for taxon " + taxonId + " ===");
	}

	/**
	 * Static method to import taxonomy data from TSV file into SQLite database
	 * 
	 * @param tsvFilePath Path to the taxonomy TSV file
	 * @param dbFilePath  Path to the SQLite database file
	 * @throws SQLException if database operations fail
	 * @throws IOException  if file operations fail
	 */
	public static void importTaxonomyFromTSV(String tsvFilePath, String dbFilePath) throws SQLException, IOException {
		LOGGER.info("Starting import of taxonomy data from " + tsvFilePath + " to " + dbFilePath);

		Connection conn = null;
		PreparedStatement insertStmt = null;

		try {
			// Create database connection
			conn = DriverManager.getConnection("jdbc:sqlite:" + dbFilePath);
			conn.setAutoCommit(false); // Use transactions for better performance

			// Create taxonomy table if it doesn't exist
			String createTableSQL = """
					CREATE TABLE IF NOT EXISTS taxonomy (
					    taxon_id INTEGER PRIMARY KEY,
					    scientific_name TEXT NOT NULL,
					    lineage TEXT,
					    rank TEXT,
					    parent INTEGER,
					    common_name TEXT,
					    links TEXT
					)
					""";

			try (Statement stmt = conn.createStatement()) {
				stmt.execute(createTableSQL);
			}

			// Create indexes for better query performance
			String[] indexSQLs = { "CREATE INDEX IF NOT EXISTS idx_taxonomy_rank ON taxonomy(rank)",
					"CREATE INDEX IF NOT EXISTS idx_taxonomy_parent ON taxonomy(parent)",
					"CREATE INDEX IF NOT EXISTS idx_taxonomy_name ON taxonomy(scientific_name)",
					"CREATE INDEX IF NOT EXISTS idx_taxonomy_name_rank ON taxonomy(scientific_name, rank)" };

			try (Statement stmt = conn.createStatement()) {
				for (String indexSQL : indexSQLs) {
					stmt.execute(indexSQL);
				}
			}

			// Clear existing data
			try (Statement stmt = conn.createStatement()) {
				stmt.execute("DELETE FROM taxonomy");
			}

			// Prepare insert statement
			String insertSQL = """
					INSERT INTO taxonomy (taxon_id, scientific_name, lineage, rank, parent, common_name, links)
					VALUES (?, ?, ?, ?, ?, ?, ?)
					""";
			insertStmt = conn.prepareStatement(insertSQL);

			// Read and process TSV file
			try (java.io.BufferedReader reader = new java.io.BufferedReader(new java.io.FileReader(tsvFilePath))) {
				String line = reader.readLine(); // Skip header
				if (line == null) {
					throw new IOException("Empty file or missing header");
				}

				// Parse header to get column indices
				String[] headers = line.split("\t");
				int taxonIdCol = -1, scientificNameCol = -1, lineageCol = -1, rankCol = -1, parentCol = -1,
						commonNameCol = -1, linksCol = -1;

				for (int i = 0; i < headers.length; i++) {
					String header = headers[i].trim();
					switch (header) {
					case "Taxon Id":
						taxonIdCol = i;
						break;
					case "Scientific name":
						scientificNameCol = i;
						break;
					case "Lineage":
						lineageCol = i;
						break;
					case "Rank":
						rankCol = i;
						break;
					case "Parent":
						parentCol = i;
						break;
					case "Common name":
						commonNameCol = i;
						break;
					case "Links":
						linksCol = i;
						break;
					}
				}

				if (taxonIdCol == -1 || scientificNameCol == -1) {
					throw new IOException("Required columns 'Taxon Id' and 'Scientific name' not found in header");
				}

				LOGGER.info("Column mapping: Taxon Id=" + taxonIdCol + ", Scientific name=" + scientificNameCol
						+ ", Lineage=" + lineageCol + ", Rank=" + rankCol + ", Parent=" + parentCol);

				int batchSize = 10000;
				int recordCount = 0;
				int batchCount = 0;

				// Process data rows
				while ((line = reader.readLine()) != null) {
					String[] fields = line.split("\t");

					if (fields.length <= Math.max(taxonIdCol, scientificNameCol)) {
						LOGGER.warn("Skipping malformed line: " + line);
						continue;
					}

					try {
						// Parse taxon ID
						int taxonId = Integer.parseInt(fields[taxonIdCol].trim());

						// Parse scientific name
						String scientificName = fields[scientificNameCol].trim();
						if (scientificName.isEmpty()) {
							LOGGER.warn("Skipping record with empty scientific name: " + taxonId);
							continue;
						}

						// Parse other fields (handle missing columns gracefully)
						String lineage = (lineageCol >= 0 && lineageCol < fields.length) ? fields[lineageCol].trim()
								: null;
						String rank = (rankCol >= 0 && rankCol < fields.length) ? fields[rankCol].trim() : null;
						String parentStr = (parentCol >= 0 && parentCol < fields.length) ? fields[parentCol].trim()
								: null;
						String commonName = (commonNameCol >= 0 && commonNameCol < fields.length)
								? fields[commonNameCol].trim()
								: null;
						String links = (linksCol >= 0 && linksCol < fields.length) ? fields[linksCol].trim() : null;

						// Parse parent ID
						Integer parent = null;
						if (parentStr != null && !parentStr.isEmpty()) {
							try {
								parent = Integer.parseInt(parentStr);
								// Log some sample parent relationships for debugging
								if (recordCount < 10 || taxonId == 131567) {
									LOGGER.debug("Taxon " + taxonId + " (" + scientificName + ") has parent " + parent);
								}
							} catch (NumberFormatException e) {
								LOGGER.warn("Invalid parent ID for taxon " + taxonId + ": " + parentStr);
							}
						} else {
							// Log taxons without parents (should only be root)
							if (recordCount < 10 || taxonId == 131567) {
								LOGGER.debug("Taxon " + taxonId + " (" + scientificName + ") has no parent (root node)");
							}
						}

						// Set parameters
						insertStmt.setInt(1, taxonId);
						insertStmt.setString(2, scientificName);
						insertStmt.setString(3, lineage);
						insertStmt.setString(4, rank);
						insertStmt.setObject(5, parent);
						insertStmt.setString(6, commonName);
						insertStmt.setString(7, links);

						insertStmt.addBatch();
						recordCount++;

						// Commit batch periodically
						if (recordCount % batchSize == 0) {
							insertStmt.executeBatch();
							conn.commit();
							batchCount++;
							LOGGER.info("Processed " + recordCount + " records (" + batchCount + " batches)");
						}

					} catch (NumberFormatException e) {
						LOGGER.warn("Skipping record with invalid taxon ID: " + fields[taxonIdCol]);
						continue;
					} catch (Exception e) {
						LOGGER.error("Error processing record: " + line, e);
						continue;
					}
				}

				// Execute remaining batch
				if (recordCount % batchSize != 0) {
					insertStmt.executeBatch();
					conn.commit();
				}

				LOGGER.info("Import completed successfully. Total records imported: " + recordCount);

			}

		} finally {
			// Clean up resources
			if (insertStmt != null) {
				try {
					insertStmt.close();
				} catch (SQLException e) {
					LOGGER.error("Error closing insert statement", e);
				}
			}
			if (conn != null) {
				try {
					conn.setAutoCommit(true);
					conn.close();
				} catch (SQLException e) {
					LOGGER.error("Error closing database connection", e);
				}
			}
		}
	}

	/**
	 * Convenience method to import the default taxonomy file to the default
	 * database
	 */
	public static void importDefaultTaxonomy() throws SQLException, IOException {
		importTaxonomyFromTSV("tests/taxonomy_all_2025_09_17.tsv", "tests/func_def.db");
	}
	
	public static void main(String[] args) {
		try {
			importTaxonomyFromTSV("C:\\Users\\daz25hak\\MetaPilot\\MetaPilot\\tests\\taxonomy_all_2025_09_17.tsv", 
					"C:\\Users\\daz25hak\\MetaPilot\\MetaPilot\\tests\\func_def.db");
		} catch (SQLException | IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}
