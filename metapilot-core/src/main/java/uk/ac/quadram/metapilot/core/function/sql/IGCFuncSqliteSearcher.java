package uk.ac.quadram.metapilot.core.function.sql;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.sql.Connection;
import java.sql.DatabaseMetaData;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.Set;

import org.dom4j.DocumentException;

import uk.ac.quadram.metapilot.core.MGYG.db.ProDbConfig;
import uk.ac.quadram.metapilot.core.function.v2.EnzymeCommission;
import uk.ac.quadram.metapilot.core.function.v2.GoObo;
import uk.ac.quadram.metapilot.task.io.pro.MetaProtein;
import uk.ac.quadram.metapilot.task.io.pro.MetaProteinAnnoEggNog;
import uk.ac.quadram.metapilot.task.io.pro.MetaProteinXMLReader2;
import uk.ac.quadram.metapilot.task.io.pro.MetaProteinXMLWriter2;
import uk.ac.quadram.metapilot.task.par.MetaConstants;

public class IGCFuncSqliteSearcher {

	protected String dbPath;
	private String[] dbNames = new String[] { "mouse", "human", "dog", "rat", "pig", "chicken" };

	protected String defDbPath;
	protected String cogDbName = "cog_20_def_tab";
	protected String nogDbName = "nog_e5_og_annotations";
	protected String ecDbName = "enzyme";
	protected String goDbName = "go_obo";
	protected String keggDbName = "kegg_pathway";

	protected HashMap<String, String[]> usedCogMap;
	protected HashMap<String, String[]> usedNogMap;
	protected HashMap<String, String> usedKeggMap;
	protected HashMap<String, GoObo> usedGoMap;
	protected HashMap<String, EnzymeCommission> usedEcMap;

	public IGCFuncSqliteSearcher(String dbpath, String defDbPath) throws SQLException, NumberFormatException, IOException {

		this.dbPath = dbpath;
		this.defDbPath = defDbPath;

		this.usedCogMap = new HashMap<String, String[]>();
		this.usedNogMap = new HashMap<String, String[]>();
		this.usedKeggMap = new HashMap<String, String>();
		this.usedGoMap = new HashMap<String, GoObo>();
		this.usedEcMap = new HashMap<String, EnzymeCommission>();
	}

	public MetaProteinAnnoEggNog[] match(MetaProtein[] metapros) throws SQLException, IOException {

		Connection dbConn = DriverManager.getConnection("jdbc:sqlite:" + this.dbPath);
		Statement dbStmt = dbConn.createStatement();

		String dbname = null;
		L: for (int i = 0; i < metapros.length; i++) {
			for (int j = 0; j < dbNames.length; j++) {
				String sql = "select * from " + dbNames[j] + " where query_name=" + "\"" + metapros[i].getName() + "\"";
				ResultSet rs = dbStmt.executeQuery(sql);
				if (rs.next()) {
					dbname = dbNames[j];
					break L;
				}
			}
		}

		ArrayList<MetaProteinAnnoEggNog> list = new ArrayList<MetaProteinAnnoEggNog>();

		if (dbname != null) {

			HashSet<String> proNameSet = new HashSet<String>();
			for (int i = 0; i < metapros.length; i++) {
				proNameSet.add(metapros[i].getName());
			}

			HashSet<String> cogSet = new HashSet<String>();
			HashSet<String> nogSet = new HashSet<String>();
			HashSet<String> keggSet = new HashSet<String>();
			HashSet<String> goSet = new HashSet<String>();
			HashSet<String> ecSet = new HashSet<String>();

			HashMap<String, String[]> proMap = new HashMap<String, String[]>();
			StringBuilder prosb = new StringBuilder();
			prosb.append("(");
			for (String pro : proNameSet) {
				prosb.append("\"").append(pro).append("\"").append(",");
			}
			prosb.deleteCharAt(prosb.length() - 1);
			prosb.append(")");

			StringBuilder sb = new StringBuilder();
			sb.append("select * from ").append(dbname).append(" where query_name IN ").append(prosb);

			ResultSet rs = dbStmt.executeQuery(sb.toString());
			while (rs.next()) {
				String[] content = new String[ProDbConfig.PRO2FUNTITLE_STRINGS.length];
				for (int i = 0; i < content.length; i++) {
					content[i] = rs.getString(ProDbConfig.PRO2FUNTITLE_STRINGS[i]);
				}
				proMap.put(content[0], content);

				String[] gos = content[5].split(",");
				for (int i = 0; i < gos.length; i++) {
					goSet.add(gos[i]);
				}

				String[] ecs = content[6].split(",");
				for (int i = 0; i < ecs.length; i++) {
					ecSet.add(ecs[i]);
				}

				String[] kos = content[7].split(",");
				for (int i = 0; i < kos.length; i++) {
					keggSet.add(kos[i]);
				}

				String[] pathways = content[8].split(",");
				for (int i = 0; i < pathways.length; i++) {
					keggSet.add(pathways[i]);
				}

				String[] cogs = content[16].split(",");
				for (int i = 0; i < cogs.length; i++) {
					cogSet.add(cogs[i]);
				}

				String[] nogs = content[17].split(",");
				for (int i = 0; i < nogs.length; i++) {
					nogSet.add(nogs[i] + "@" + content[2]);
				}
			}

			for (int i = 0; i < metapros.length; i++) {
				String proName = metapros[i].getName();
				if (proMap.containsKey(proName)) {
					String[] content = proMap.get(proName);
					MetaProteinAnnoEggNog mpa = new MetaProteinAnnoEggNog(metapros[i], content[2], content[3], content[4],
							content[5], content[6], content[7], content[8], content[9], content[10], content[11],
							content[12], content[13], content[14], content[15], content[16], content[17]);
					list.add(mpa);
				}
			}

			getReference(cogSet, nogSet, goSet, ecSet, keggSet);
		}

		dbStmt.close();
		dbConn.close();

		MetaProteinAnnoEggNog[] mpas = list.toArray(new MetaProteinAnnoEggNog[list.size()]);

		return mpas;
	}
	
	protected void getReference(HashSet<String> cogSet, HashSet<String> nogSet, HashSet<String> goSet,
			HashSet<String> ecSet, HashSet<String> keggSet) throws SQLException {

		try (Connection defDbConn = DriverManager.getConnection("jdbc:sqlite:" + this.defDbPath)) {
			if (hasTable(defDbConn, "annotation_term")) {
				getReferenceFromNormalizedDb(defDbConn, cogSet, nogSet, goSet, ecSet, keggSet);
			} else {
				getReferenceFromLegacyDb(defDbConn, cogSet, nogSet, goSet, ecSet, keggSet);
			}
		}
	}

	private void getReferenceFromNormalizedDb(Connection conn, HashSet<String> cogSet, HashSet<String> nogSet,
			HashSet<String> goSet, HashSet<String> ecSet, HashSet<String> keggSet) throws SQLException {

		for (ReferenceRow row : queryTerms(conn, "COG", cogSet)) {
			usedCogMap.put(row.id, new String[] { row.category, row.description });
		}

		for (ReferenceRow row : queryNogTerms(conn, nogSet)) {
			usedNogMap.put(row.id, new String[] { row.category, row.description });
		}

		for (ReferenceRow row : queryEcTerms(conn, ecSet)) {
			ArrayList<String> anList = new ArrayList<String>();
			for (String an : row.synonyms.split(";")) {
				if (an != null && an.length() > 0) {
					anList.add(an);
				}
			}
			usedEcMap.put(row.id, new EnzymeCommission(row.id, row.description, row.reaction, anList));
		}

		for (ReferenceRow row : queryTerms(conn, "GO", goSet)) {
			usedGoMap.put(row.id, new GoObo(row.id, firstNonEmpty(row.name, row.description), row.namespace));
		}

		for (String kegg : cleanIds(keggSet)) {
			ReferenceRow row = queryKeggTerm(conn, kegg);
			if (row != null) {
				usedKeggMap.put(kegg, row.description);
			}
		}
	}

	private void getReferenceFromLegacyDb(Connection conn, HashSet<String> cogSet, HashSet<String> nogSet,
			HashSet<String> goSet, HashSet<String> ecSet, HashSet<String> keggSet) throws SQLException {

		try (Statement stmt = conn.createStatement()) {
			for (String id : cleanIds(cogSet)) {
				try (ResultSet rs = stmt.executeQuery(
						"select id, category, description from " + cogDbName + " where id=" + sqlQuote(id))) {
					while (rs.next()) {
						usedCogMap.put(rs.getString("id"),
								new String[] { rs.getString("category"), rs.getString("description") });
					}
				}
			}

			for (String id : cleanIds(nogSet)) {
				try (ResultSet rs = stmt.executeQuery(
						"select id, category, description from " + nogDbName + " where id=" + sqlQuote(id))) {
					while (rs.next()) {
						usedNogMap.put(rs.getString("id"),
								new String[] { rs.getString("category"), rs.getString("description") });
					}
				}
			}

			for (String id : cleanIds(ecSet)) {
				try (ResultSet rs = stmt.executeQuery(
						"select id, de, ca, an from " + ecDbName + " where id=" + sqlQuote(id))) {
					while (rs.next()) {
						ArrayList<String> anList = new ArrayList<String>();
						String anString = rs.getString("an");
						if (anString != null) {
							for (String an : anString.split(";")) {
								if (an.length() > 0) anList.add(an);
							}
						}
						String ecId = rs.getString("id");
						usedEcMap.put(ecId, new EnzymeCommission(ecId, rs.getString("de"), rs.getString("ca"), anList));
					}
				}
			}

			for (String id : cleanIds(goSet)) {
				try (ResultSet rs = stmt.executeQuery(
						"select id, name, namespace from " + goDbName + " where id=" + sqlQuote(id))) {
					while (rs.next()) {
						usedGoMap.put(rs.getString("id"),
								new GoObo(rs.getString("id"), rs.getString("name"), rs.getString("namespace")));
					}
				}
			}

			for (String id : cleanIds(keggSet)) {
				try (ResultSet rs = stmt.executeQuery(
						"select name, description from " + keggDbName + " where name=" + sqlQuote(id))) {
					while (rs.next()) {
						usedKeggMap.put(id, rs.getString("description"));
					}
				}
			}
		}
	}

	private static ArrayList<ReferenceRow> queryTerms(Connection conn, String source, Set<String> ids)
			throws SQLException {
		ArrayList<ReferenceRow> rows = new ArrayList<ReferenceRow>();
		String sql = "select source, term_id, name, description, namespace, category from annotation_term "
				+ "where source=? and term_id=?";
		try (PreparedStatement ps = conn.prepareStatement(sql)) {
			for (String id : cleanIds(ids)) {
				ps.setString(1, source);
				ps.setString(2, id);
				try (ResultSet rs = ps.executeQuery()) {
					while (rs.next()) {
						rows.add(ReferenceRow.fromTerm(rs, id));
					}
				}
			}
		}
		return rows;
	}

	private static ArrayList<ReferenceRow> queryNogTerms(Connection conn, Set<String> ids) throws SQLException {
		ArrayList<ReferenceRow> rows = new ArrayList<ReferenceRow>();
		if (hasTable(conn, "eggnog_og_annotation")) {
			String sql = "select id, category, description from eggnog_og_annotation where id=?";
			try (PreparedStatement ps = conn.prepareStatement(sql)) {
				for (String id : cleanIds(ids)) {
					ps.setString(1, id);
					try (ResultSet rs = ps.executeQuery()) {
						while (rs.next()) {
							ReferenceRow row = new ReferenceRow();
							row.id = rs.getString("id");
							row.category = rs.getString("category");
							row.description = rs.getString("description");
							rows.add(row);
						}
					}
				}
			}
		} else {
			rows.addAll(queryLegacyLikeRows(conn, "nog_e5_og_annotations", ids));
		}
		return rows;
	}

	private static ArrayList<ReferenceRow> queryEcTerms(Connection conn, Set<String> ids) throws SQLException {
		ArrayList<ReferenceRow> rows = new ArrayList<ReferenceRow>();
		String sql = "select t.term_id, t.description, "
				+ "coalesce(ca.property_value, '') as reaction, coalesce(an.property_value, '') as synonyms "
				+ "from annotation_term t "
				+ "left join annotation_property ca on ca.source=t.source and ca.term_id=t.term_id and ca.property_key='reaction' "
				+ "left join annotation_property an on an.source=t.source and an.term_id=t.term_id and an.property_key='synonyms' "
				+ "where t.source='EC' and t.term_id=?";
		try (PreparedStatement ps = conn.prepareStatement(sql)) {
			for (String id : cleanIds(ids)) {
				ps.setString(1, id);
				try (ResultSet rs = ps.executeQuery()) {
					while (rs.next()) {
						ReferenceRow row = new ReferenceRow();
						row.id = rs.getString("term_id");
						row.description = rs.getString("description");
						row.reaction = rs.getString("reaction");
						row.synonyms = rs.getString("synonyms");
						rows.add(row);
					}
				}
			}
		}
		return rows;
	}

	private static ReferenceRow queryTermByIdOrAlias(Connection conn, String source, String id) throws SQLException {
		String sql = "select t.source, t.term_id, t.name, t.description, t.namespace, t.category "
				+ "from annotation_term t where t.source=? and t.term_id=? "
				+ "union all "
				+ "select t.source, t.term_id, t.name, t.description, t.namespace, t.category "
				+ "from annotation_alias a join annotation_term t on t.source=a.source and t.term_id=a.canonical_id "
				+ "where a.source=? and a.alias=? limit 1";
		try (PreparedStatement ps = conn.prepareStatement(sql)) {
			ps.setString(1, source);
			ps.setString(2, id);
			ps.setString(3, source);
			ps.setString(4, id);
			try (ResultSet rs = ps.executeQuery()) {
				if (rs.next()) {
					return ReferenceRow.fromTerm(rs, id);
				}
			}
		}
		return null;
	}

	private static ReferenceRow queryKeggTerm(Connection conn, String id) throws SQLException {
		String normalized = id == null ? "" : id.trim();
		if (normalized.length() == 0) {
			return null;
		}

		ArrayList<String> sources = new ArrayList<String>();
		if (normalized.startsWith("ko:") || normalized.matches("^K\\d{5}$")) {
			sources.add("KEGG_ko");
		} else if (normalized.startsWith("path:") || normalized.startsWith("map") || normalized.startsWith("ko")) {
			sources.add("KEGG_Pathway");
		} else if (normalized.startsWith("md:") || normalized.matches("^M\\d{5}$")) {
			sources.add("KEGG_Module");
		} else if (normalized.startsWith("rn:") || normalized.matches("^R\\d{5}$")) {
			sources.add("KEGG_Reaction");
		} else if (normalized.startsWith("rc:") || normalized.matches("^RC\\d{5}$")) {
			sources.add("KEGG_rclass");
		}

		sources.add("KEGG_ko");
		sources.add("KEGG_Pathway");
		sources.add("KEGG_Module");
		sources.add("KEGG_Reaction");
		sources.add("KEGG_rclass");

		LinkedHashSet<String> uniqueSources = new LinkedHashSet<String>(sources);
		for (String source : uniqueSources) {
			ReferenceRow row = queryTermByIdOrAlias(conn, source, normalized);
			if (row != null) {
				return row;
			}
		}
		return null;
	}

	private static ArrayList<ReferenceRow> queryLegacyLikeRows(Connection conn, String table, Set<String> ids)
			throws SQLException {
		ArrayList<ReferenceRow> rows = new ArrayList<ReferenceRow>();
		String sql = "select id, category, description from " + table + " where id=?";
		try (PreparedStatement ps = conn.prepareStatement(sql)) {
			for (String id : cleanIds(ids)) {
				ps.setString(1, id);
				try (ResultSet rs = ps.executeQuery()) {
					while (rs.next()) {
						ReferenceRow row = new ReferenceRow();
						row.id = rs.getString("id");
						row.category = rs.getString("category");
						row.description = rs.getString("description");
						rows.add(row);
					}
				}
			}
		}
		return rows;
	}

	private static Set<String> cleanIds(Set<String> ids) {
		if (ids == null || ids.isEmpty()) {
			return Collections.emptySet();
		}
		LinkedHashSet<String> clean = new LinkedHashSet<String>();
		for (String id : ids) {
			if (id == null) continue;
			String value = id.trim();
			if (value.length() == 0 || value.equals("-")) continue;
			clean.add(value);
		}
		return clean;
	}

	private static boolean hasTable(Connection conn, String name) throws SQLException {
		DatabaseMetaData meta = conn.getMetaData();
		try (ResultSet rs = meta.getTables(null, null, name, new String[] { "TABLE" })) {
			if (rs.next()) return true;
		}
		try (ResultSet rs = meta.getTables(null, null, name, new String[] { "VIEW" })) {
			return rs.next();
		}
	}

	private static String firstNonEmpty(String first, String second) {
		if (first != null && first.length() > 0) return first;
		return second == null ? "" : second;
	}

	private static String sqlQuote(String value) {
		return "'" + value.replace("'", "''") + "'";
	}

	private static final class ReferenceRow {
		String id = "";
		String name = "";
		String description = "";
		String namespace = "";
		String category = "";
		String reaction = "";
		String synonyms = "";

		static ReferenceRow fromTerm(ResultSet rs, String requestedId) throws SQLException {
			ReferenceRow row = new ReferenceRow();
			row.id = firstNonEmpty(requestedId, rs.getString("term_id"));
			row.name = rs.getString("name");
			row.description = rs.getString("description");
			row.namespace = rs.getString("namespace");
			row.category = rs.getString("category");
			return row;
		}
	}
	
	public HashMap<String, String[]> getUsedCogMap() {
		return usedCogMap;
	}

	public HashMap<String, String[]> getUsedNogMap() {
		return usedNogMap;
	}

	public HashMap<String, String> getUsedKeggMap() {
		return usedKeggMap;
	}

	public HashMap<String, GoObo> getUsedGoMap() {
		return usedGoMap;
	}

	public HashMap<String, EnzymeCommission> getUsedEcMap() {
		return usedEcMap;
	}

	@SuppressWarnings("unused")
	private void test(String in)
			throws IOException, NumberFormatException, SQLException, DocumentException {
		ArrayList<MetaProtein> list = new ArrayList<MetaProtein>();
		BufferedReader reader = new BufferedReader(new FileReader(in));
		String line = reader.readLine();
		String[] title = line.split("\t");
		String[] fileNames = new String[title.length - 10];
		System.arraycopy(title, 6, fileNames, 0, fileNames.length);

		int count = 1;
		while ((line = reader.readLine()) != null) {
			String[] cs = line.split("\t");
			if (!cs[1].contains("ECOLI")) {
				double[] intensity = new double[cs.length - 10];
				for (int i = 0; i < intensity.length; i++) {
					intensity[i] = Double.parseDouble(cs[i + 6]);
				}
				MetaProtein mp = new MetaProtein(count++, 1, cs[1], Integer.parseInt(cs[2]), Integer.parseInt(cs[3]),
						Double.parseDouble(cs[4]), intensity);
				list.add(mp);
			}
		}
		reader.close();

		MetaProtein[] mps = list.toArray(new MetaProtein[list.size()]);
		MetaProteinAnnoEggNog[] metapros = this.match(mps);

		File folder = (new File(in)).getParentFile();

		MetaProteinXMLWriter2 writer = new MetaProteinXMLWriter2(new File(folder, "function.xml"), MetaConstants.pFind,
				MetaConstants.labelFree, fileNames);

		writer.addProteins(metapros, this.getUsedCogMap(), this.getUsedNogMap(), this.getUsedKeggMap(),
				this.getUsedGoMap(), this.getUsedEcMap());

		writer.close();

		MetaProteinXMLReader2 funReader = new MetaProteinXMLReader2(new File(folder, "function.xml"));
		funReader.exportTsv(new File(folder, "function.tsv"));
	}
}
