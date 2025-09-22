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
			// Cache all main taxonomic ranks (most commonly accessed)
			String[] mainRanks = { "superkingdom", "phylum", "class", "order", "family", "genus", "species" };

			for (String rank : mainRanks) {
				getTaxonByRankStmt.setString(1, rank);
				ResultSet rs = getTaxonByRankStmt.executeQuery();

				while (rs.next()) {
					Taxon taxon = createTaxonFromResultSet(rs);
					hierarchyCache.put(taxon.getId(), taxon);
					nameCache.put(taxon.getName(), taxon.getId());
					taxonNameMapCache.put(taxon.getName() + "_" + taxon.getRankId(), taxon.getId());
				}
				rs.close();
			}

			LOGGER.info("Built hierarchy cache with " + hierarchyCache.size() + " main rank taxa");

			// Also cache taxons that are commonly referenced as parents
			// This helps with hierarchy traversal
			try (Statement stmt = connection.createStatement()) {
				// Cache taxons that are referenced as parents by many other taxons
				ResultSet rs = stmt.executeQuery(
					"SELECT t.* FROM taxonomy t " +
					"WHERE t.taxon_id IN (" +
					"  SELECT parent FROM taxonomy " +
					"  WHERE parent IS NOT NULL " +
					"  GROUP BY parent " +
					"  HAVING COUNT(*) > 100" +
					")"
				);
				
				int parentCacheCount = 0;
				while (rs.next()) {
					Taxon taxon = createTaxonFromResultSet(rs);
					if (!hierarchyCache.containsKey(taxon.getId())) {
						hierarchyCache.put(taxon.getId(), taxon);
						parentCacheCount++;
					}
				}
				LOGGER.info("Added " + parentCacheCount + " frequently referenced parent taxons to hierarchy cache");
			}

			LOGGER.info("Total hierarchy cache size: " + hierarchyCache.size() + " taxa");

		} catch (SQLException e) {
			LOGGER.error("Error building hierarchy cache", e);
		}
	}

	private Taxon createTaxonFromResultSet(ResultSet rs) throws SQLException {
		int taxonId = rs.getInt("taxon_id");
		String name = rs.getString("scientific_name");
		String rankStr = rs.getString("rank");
		String lineage = rs.getString("lineage");
		int parent = rs.getInt("parent");

		// Parse rank
		int rank = -1;
		if (rankStr != null && TaxonomyRanks.getRankIdMap().containsKey(rankStr)) {
			rank = TaxonomyRanks.getRankIdMap().get(rankStr);
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
		// Check hierarchy cache first (fastest)
		if (hierarchyCache.containsKey(taxonId)) {
			LOGGER.debug("Found taxon ID " + taxonId + " in hierarchy cache");
			return hierarchyCache.get(taxonId);
		}

		// Check LRU cache
		synchronized (lruCache) {
			if (lruCache.containsKey(taxonId)) {
				LOGGER.debug("Found taxon ID " + taxonId + " in LRU cache");
				return lruCache.get(taxonId);
			}
		}

		// Query database
		try {
			LOGGER.debug("Querying database for taxon ID: " + taxonId);
			getTaxonStmt.setInt(1, taxonId);
			ResultSet rs = getTaxonStmt.executeQuery();

			if (rs.next()) {
				Taxon taxon = createTaxonFromResultSet(rs);
				LOGGER.debug("Found taxon ID " + taxonId + " in database: " + taxon.getName());

				// Add to LRU cache
				synchronized (lruCache) {
					lruCache.put(taxonId, taxon);
				}

				// Update name cache
				nameCache.put(taxon.getName(), taxonId);
				taxonNameMapCache.put(taxon.getName() + "_" + taxon.getRankId(), taxonId);

				rs.close();
				return taxon;
			} else {
				LOGGER.warn("No taxon found in database for ID: " + taxonId);
			}
			rs.close();

		} catch (SQLException e) {
			LOGGER.error("Error querying taxon by ID: " + taxonId, e);
		}

		return null;
	}

	@Override
	public Taxon getTaxonFromName(String name_rankId) {
		// Check cache first
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
		// Check cache first
		String key = name + "_" + rankId;
		if (taxonNameMapCache.containsKey(key)) {
			return getTaxonFromId(taxonNameMapCache.get(key));
		}

		// Query database
		try {
			String rankStr = TaxonomyRanks.getRankFromId(rankId).toString().toLowerCase();
			getTaxonByNameStmt.setString(1, name);
			getTaxonByNameStmt.setString(2, rankStr);
			ResultSet rs = getTaxonByNameStmt.executeQuery();

			if (rs.next()) {
				Taxon taxon = createTaxonFromResultSet(rs);

				// Add to caches
				synchronized (lruCache) {
					lruCache.put(taxon.getId(), taxon);
				}
				nameCache.put(name, taxon.getId());
				taxonNameMapCache.put(key, taxon.getId());

				rs.close();
				return taxon;
			}
			rs.close();

		} catch (SQLException e) {
			LOGGER.error("Error querying taxon by name: " + name, e);
		}

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
			
			switch (rankid) {
			case -1:
				break;
			case 0:
				taxonomy[0] = taxon.getId();
				break;
			case 1:
				taxonomy[1] = taxon.getId();
				break;
			case 4:
				taxonomy[2] = taxon.getId();
				break;
			case 7:
				taxonomy[3] = taxon.getId();
				break;
			case 11:
				taxonomy[4] = taxon.getId();
				break;
			case 16:
				taxonomy[5] = taxon.getId();
				break;
			case 20:
				taxonomy[6] = taxon.getId();
				break;
			case 23:
				taxonomy[7] = taxon.getId();
				break;
			}
			int parentId = taxon.getParentId();
			LOGGER.debug("Looking up parent taxon ID: " + parentId + " for taxon: " + taxon.getId());
			
			// Check if parent is in cache first
			boolean inHierarchyCache = hierarchyCache.containsKey(parentId);
			boolean inLruCache = lruCache.containsKey(parentId);
			LOGGER.debug("Parent " + parentId + " in hierarchy cache: " + inHierarchyCache + 
						", in LRU cache: " + inLruCache);
			
			taxon = getTaxonFromId(parentId);
			if (taxon == null) {
				LOGGER.warn("Parent taxon ID " + parentId + " not found in database, stopping hierarchy traversal at depth " + traversalDepth);
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
		while (taxon != null) {
			int rankid = taxon.getRankId();
			switch (rankid) {
			case 0:
				taxonomy[0] = taxon.getName();
				break;
			case 1:
				taxonomy[1] = taxon.getName();
				break;
			case 4:
				taxonomy[2] = taxon.getName();
				break;
			case 7:
				taxonomy[3] = taxon.getName();
				break;
			case 11:
				taxonomy[4] = taxon.getName();
				break;
			case 16:
				taxonomy[5] = taxon.getName();
				break;
			case 20:
				taxonomy[6] = taxon.getName();
				break;
			case 23:
				taxonomy[7] = taxon.getName();
				break;
			}
			int parentId = taxon.getParentId();
			taxon = getTaxonFromId(parentId);
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
