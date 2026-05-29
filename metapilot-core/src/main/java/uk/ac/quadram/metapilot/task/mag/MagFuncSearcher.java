package uk.ac.quadram.metapilot.task.mag;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import uk.ac.quadram.metapilot.core.MGYG.db.ProDbConfig;
import uk.ac.quadram.metapilot.core.function.sql.IGCFuncSqliteSearcher;
import uk.ac.quadram.metapilot.core.taxonomy.TaxonomyRanks;
import uk.ac.quadram.metapilot.core.tools.ResourceLoader;
import uk.ac.quadram.metapilot.task.io.pro.MetaProtein;

public class MagFuncSearcher extends IGCFuncSqliteSearcher {

	private String genome;
	private MagDbItem usedMagDbItem;

	public MagFuncSearcher(String dbpath, String defDbPath, MagDbItem usedMagDbItem)
			throws SQLException, NumberFormatException, IOException {
		super(dbpath, defDbPath);
		this.usedMagDbItem = usedMagDbItem;
	}

	public MagFuncSearcher(String dbpath, String defDbPath, MagDbItem usedMagDbItem, String genome)
			throws SQLException, NumberFormatException, IOException {
		super(dbpath, defDbPath);
		this.genome = genome;
		this.usedMagDbItem = usedMagDbItem;
	}

	public HashMap<String, String[]> matchGenome(String[] genomes) throws IOException {

		HashMap<String, String> gnMap = new HashMap<String, String>();
		BufferedReader gtdbNcbiReader = new BufferedReader(
				new InputStreamReader(ResourceLoader.load("gtdb_ncbi"), "UTF-8"));
		String gnLine = gtdbNcbiReader.readLine();
		while ((gnLine = gtdbNcbiReader.readLine()) != null) {
			String[] cs = gnLine.split("\t");
			gnMap.put(cs[0], cs[1]);
		}
		gtdbNcbiReader.close();

		File genomeTaxaFile = usedMagDbItem.getOriginalTaxaFile();
		HashMap<String, String[]> genomeTaxaMap = new HashMap<String, String[]>();
		for (int i = 0; i < genomes.length; i++) {
			genomeTaxaMap.put(genomes[i], new String[ProDbConfig.genome2TaxaGTDBNCBITitle.length]);
		}

		if (genomeTaxaFile.exists()) {
			BufferedReader reader = new BufferedReader(new FileReader(genomeTaxaFile));
			String line = reader.readLine();
			int lineageId = -1;
			String[] title = line.split("\t");
			for (int i = 0; i < title.length; i++) {
				if (title[i].equals("Lineage") || title[i].equals("Species_designation")) {
					lineageId = i;
				}
			}
			while ((line = reader.readLine()) != null) {
				String[] cs = line.split("\t");
				if (genomeTaxaMap.containsKey(cs[0])) {

					String gtdb = cs[lineageId];
					ArrayList<String> list = new ArrayList<String>();
					list.add(cs[0]);

					String rank = "";
					String[] lineage = gtdb.split(";");
					lineage = normalizeLineage(lineage);
					for (int i = lineage.length - 1; i >= 0; i--) {
						if (lineage[i].length() > 3) {
							list.add(lineage[i].substring(3));

							char rr = lineage[i].charAt(0);

							switch (rr) {
							case 'd': {
								rank = TaxonomyRanks.Superkingdom.getName();
								break;
							}
							case 'p': {
								rank = TaxonomyRanks.Phylum.getName();
								break;
							}
							case 'c': {
								rank = TaxonomyRanks.Class.getName();
								break;
							}
							case 'o': {
								rank = TaxonomyRanks.Order.getName();
								break;
							}
							case 'f': {
								rank = TaxonomyRanks.Family.getName();
								break;
							}
							case 'g': {
								rank = TaxonomyRanks.Genus.getName();
								break;
							}
							case 's': {
								rank = TaxonomyRanks.Species.getName();
								break;
							}
							default:
								break;
							}

							break;
						}
					}
					list.add(rank);

					for (int i = 0; i < lineage.length; i++) {
						String name = "";
						if (lineage[i].length() > 3) {
							name = lineage[i].substring(3);
						}
						if (i == 0 && name.length() == 0) {
							name = inferDomain(lineage);
						}
						if (i == 0) {
							list.add(name);
							list.add("");
						} else {
							list.add(name);
						}
					}

					if (gnMap.containsKey(gtdb)) {
						String ncbi = gnMap.get(gtdb);
						rank = "";
						lineage = ncbi.split(";");
						lineage = normalizeLineage(lineage);
						for (int i = lineage.length - 1; i >= 0; i--) {
							if (lineage[i].length() > 3) {
								list.add(lineage[i].substring(3));

								char rr = lineage[i].charAt(0);

								switch (rr) {
								case 'd': {
									rank = TaxonomyRanks.Superkingdom.getName();
									break;
								}
								case 'p': {
									rank = TaxonomyRanks.Phylum.getName();
									break;
								}
								case 'c': {
									rank = TaxonomyRanks.Class.getName();
									break;
								}
								case 'o': {
									rank = TaxonomyRanks.Order.getName();
									break;
								}
								case 'f': {
									rank = TaxonomyRanks.Family.getName();
									break;
								}
								case 'g': {
									rank = TaxonomyRanks.Genus.getName();
									break;
								}
								case 's': {
									rank = TaxonomyRanks.Species.getName();
									break;
								}
								default:
									break;
								}

								break;
							}
						}
						list.add(rank);

						for (int i = 0; i < lineage.length; i++) {
							String name = "";
							if (lineage[i].length() > 3) {
								name = lineage[i].substring(3);
							}
							if (i == 0 && name.length() == 0) {
								name = inferDomain(lineage);
							}
							if (i == 0) {
								list.add(name);
								list.add("");
							} else {
								list.add(name);
							}
						}
					} else {
						for (int i = 0; i < 10; i++) {
							list.add("");
						}
					}

					genomeTaxaMap.put(cs[0], list.toArray(new String[list.size()]));
				}
			}
			reader.close();
		} else {
			throw new IOException("Taxonomic annotation file " + genomeTaxaFile + " was not found.");
		}

		return genomeTaxaMap;
	}

	private static String[] normalizeLineage(String[] lineage) {
		if (lineage == null || lineage.length == 0) {
			return new String[0];
		}
		String domain = inferDomain(lineage);
		String first = lineage[0] == null ? "" : lineage[0].trim();
		if (first.length() > 0 && Character.toLowerCase(first.charAt(0)) == 'd') {
			if (first.length() <= 3 || first.substring(3).trim().isEmpty()) {
				String[] normalized = lineage.clone();
				normalized[0] = "d__" + domain;
				return normalized;
			}
			return lineage;
		}
		String[] normalized = new String[lineage.length + 1];
		normalized[0] = "d__" + domain;
		System.arraycopy(lineage, 0, normalized, 1, lineage.length);
		return normalized;
	}

	private static String inferDomain(String[] lineage) {
		for (String taxon : lineage) {
			if (taxon == null || taxon.length() < 3) {
				continue;
			}
			char prefix = Character.toLowerCase(taxon.charAt(0));
			if (prefix == 'd' && taxon.length() > 3 && !taxon.substring(3).trim().isEmpty()) {
				return taxon.substring(3).trim();
			}
			if (prefix == 'p' || prefix == 'c' || prefix == 'o' || prefix == 'f' || prefix == 'g' || prefix == 's') {
				return "Bacteria";
			}
		}
		return "";
	}

	public MetaProteinAnnoMag[] match(MetaProtein[] metapros) throws IOException, SQLException {

		if (genome == null) {
			return null;
		}

		ArrayList<MetaProteinAnnoMag> list = new ArrayList<MetaProteinAnnoMag>();
		HashSet<String> cogSet = new HashSet<String>();
		HashSet<String> nogSet = new HashSet<String>();
		HashSet<String> keggSet = new HashSet<String>();
		HashSet<String> goSet = new HashSet<String>();
		HashSet<String> ecSet = new HashSet<String>();

		File originalDbFolder = usedMagDbItem.getOriginalDbFolder();
		File eggnogFolder = usedMagDbItem.getOriginalFuncFolder();

		HashMap<String, MetaProtein> proMap = new HashMap<String, MetaProtein>();
		for (MetaProtein mp : metapros) {
			proMap.put(mp.getName(), mp);
		}

		HashMap<String, String> nameMap = new HashMap<String, String>();

		File fastaFile = new File(originalDbFolder, genome + ".faa");
		if (!fastaFile.exists()) {
			fastaFile = new File(originalDbFolder, genome + ".1.faa");
		}
		if (fastaFile.exists()) {
			BufferedReader reader = new BufferedReader(new FileReader(fastaFile));
			String line = null;
			while ((line = reader.readLine()) != null) {
				if (line.startsWith(">")) {
					int id = line.indexOf(" ");
					String proIdString = line.substring(1, id);
					String proNameString = line.substring(id + 1);

					if (proMap.containsKey(proIdString)) {
						nameMap.put(proIdString, proNameString);
					}
				}
			}
			reader.close();
		}

		File eggnogFile = null;
		File[] eggnogFiles = eggnogFolder.listFiles();
		if (eggnogFiles != null) {
			for (int i = 0; i < eggnogFiles.length; i++) {
				String name = eggnogFiles[i].getName();
				if (name.startsWith(genome)) {
					eggnogFile = eggnogFiles[i];
				}
			}
		}

		if (eggnogFile != null && eggnogFile.exists()) {
			BufferedReader reader = new BufferedReader(new FileReader(eggnogFile));
			reader.readLine();
			String line;
			while ((line = reader.readLine()) != null) {
				String[] content = line.split("\t");
				if (proMap.containsKey(content[0])) {
					MetaProtein mp = proMap.get(content[0]);

					String[] cognogs = content[MetaProteinAnnoMag.nogID].split(",");
					String cogRoot = "";
					String cogBacteria = "";

					StringBuilder nogSb = new StringBuilder();

					for (int i = 0; i < cognogs.length; i++) {
						if (cognogs[i].startsWith("COG")) {

							if (cognogs[i].endsWith("root")) {
								cogRoot = cognogs[i].split("@")[0];
							} else if (cognogs[i].endsWith("Bacteria")) {
								cogBacteria = cognogs[i].split("@")[0];
							}
						} else {
							String nog = cognogs[i].split("\\|")[0];
							nogSb.append(cognogs[i]).append(",");
							nogSet.add(nog);
						}
					}

					String cog = "";
					if (cogRoot.length() > 0) {
						if (cogBacteria.length() > 0) {
							if (cogRoot.equals(cogBacteria)) {
								cog = cogRoot;
								cogSet.add(cogRoot);
							} else {
								cog = cogRoot + "," + cogBacteria;
								cogSet.add(cogRoot);
								cogSet.add(cogBacteria);
							}
						} else {
							cog = cogRoot;
							cogSet.add(cogRoot);
						}
					} else {
						if (cogBacteria.length() > 0) {
							cog = cogBacteria;
							cogSet.add(cogBacteria);
						}
					}

					String nog = "";
					if (nogSb.length() > 0) {
						nog = nogSb.substring(0, nogSb.length() - 1);
					}

					String[] gos = content[MetaProteinAnnoMag.goID].split(",");
					for (int i = 0; i < gos.length; i++) {
						goSet.add(gos[i]);
					}

					String[] ecs = content[MetaProteinAnnoMag.ecID].split(",");
					for (int i = 0; i < ecs.length; i++) {
						ecSet.add(ecs[i]);
					}

					String[] kos = content[MetaProteinAnnoMag.koID].split(",");
					for (int i = 0; i < kos.length; i++) {
						keggSet.add(kos[i]);
					}

					String[] pathways = content[MetaProteinAnnoMag.pathwayID].split(",");
					for (int i = 0; i < pathways.length; i++) {
						keggSet.add(pathways[i]);
					}

					MetaProteinAnnoMag mpa = new MetaProteinAnnoMag(mp, content, nameMap.get(content[0]), cog, nog);
					list.add(mpa);
				}

			}
			reader.close();
		} else {
			throw new IOException("Functional annotation file " + eggnogFile + " was not found.");
		}

		getReference(cogSet, nogSet, goSet, ecSet, keggSet);

		return list.toArray(new MetaProteinAnnoMag[list.size()]);
	}

}
