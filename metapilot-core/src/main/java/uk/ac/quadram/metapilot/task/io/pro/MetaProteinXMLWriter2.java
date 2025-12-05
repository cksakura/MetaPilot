package uk.ac.quadram.metapilot.task.io.pro;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.dom4j.Document;
import org.dom4j.DocumentFactory;
import org.dom4j.DocumentHelper;
import org.dom4j.Element;
import org.dom4j.io.OutputFormat;
import org.dom4j.io.XMLWriter;

import uk.ac.quadram.metapilot.core.function.v2.EnzymeCommission;
import uk.ac.quadram.metapilot.core.function.v2.GoObo;
import uk.ac.quadram.metapilot.core.taxonomy.Taxon;
import uk.ac.quadram.metapilot.task.hgm.MetaProteinAnnoHgm;
import uk.ac.quadram.metapilot.task.mag.MetaProteinAnnoMag;
import uk.ac.quadram.metapilot.task.par.MetaConstants;

public class MetaProteinXMLWriter2 {

	protected String searchType;
	protected String quanType;
	protected File output;
	protected Document document;
	protected XMLWriter writer;
	protected Element root;
	protected Element proRoot;
	protected Element taxonRoot;

	private HashSet<Integer> taxonIdSet;

	protected static final Logger LOGGER = LogManager.getLogger(MetaProteinXMLWriter2.class);

	public MetaProteinXMLWriter2(String output, String searchType, String quanType, String[] fileNames) {
		this(new File(output), searchType, quanType, fileNames);
		// TODO Auto-generated constructor stub
	}

	public MetaProteinXMLWriter2(File output, String searchType, String quanType, String[] fileNames) {
		// TODO Auto-generated constructor stub

		this.searchType = searchType;
		this.quanType = quanType;
		this.output = output;
		this.document = DocumentHelper.createDocument();
		OutputStreamWriter bufferedWriter = null;
		try {
			bufferedWriter = new OutputStreamWriter(new FileOutputStream(output), "UTF8");
		} catch (UnsupportedEncodingException | FileNotFoundException e) {
			// TODO Auto-generated catch block
			LOGGER.error("Error in writing proteins to " + output, e);
		}
		this.writer = new XMLWriter(bufferedWriter, OutputFormat.createPrettyPrint());
		this.root = DocumentFactory.getInstance().createElement("Protein_function");
		this.document.setRootElement(root);

		StringBuilder sb = new StringBuilder();
		for (String name : fileNames) {
			sb.append(name).append(";");
		}
		if (sb.length() > 0) {
			sb.deleteCharAt(sb.length() - 1);
		}

		this.taxonIdSet = new HashSet<Integer>();

		this.root.addAttribute("fileNames", sb.toString());
		this.root.addAttribute("searchType", searchType);
		this.root.addAttribute("quanMode", quanType);

		this.proRoot = DocumentFactory.getInstance().createElement("Proteins");
		this.root.add(proRoot);

		this.taxonRoot = DocumentFactory.getInstance().createElement("Taxa");
		this.root.add(taxonRoot);
	}

	public void addProteins(MetaProteinAnnoEggNog[] proteins, HashMap<String, String[]> usedCogMap,
			HashMap<String, String[]> usedNogMap, HashMap<String, String> usedKeggMap, HashMap<String, GoObo> usedGoMap,
			HashMap<String, EnzymeCommission> usedEcMap) {
		this.root.addAttribute("Protein_count", String.valueOf(proteins.length));

		Element cogRoot = DocumentFactory.getInstance().createElement("COGs");
		for (String key : usedCogMap.keySet()) {
			Element e = DocumentFactory.getInstance().createElement("COG");
			e.addAttribute("id", key);
			e.addAttribute("category", usedCogMap.get(key)[0]);
			e.addAttribute("function", usedCogMap.get(key)[1]);
			cogRoot.add(e);
		}
		this.root.add(cogRoot);
		this.root.addAttribute("COG_count", String.valueOf(usedCogMap.size()));

		Element nogRoot = DocumentFactory.getInstance().createElement("NOGs");
		for (String key : usedNogMap.keySet()) {
			Element e = DocumentFactory.getInstance().createElement("NOG");
			e.addAttribute("id", key);
			e.addAttribute("category", usedNogMap.get(key)[0]);
			e.addAttribute("function", usedNogMap.get(key)[1]);
			nogRoot.add(e);
		}
		this.root.add(nogRoot);
		this.root.addAttribute("NOG_count", String.valueOf(usedNogMap.size()));

		Element keggRoot = DocumentFactory.getInstance().createElement("KEGGs");
		for (String key : usedKeggMap.keySet()) {
			Element e = DocumentFactory.getInstance().createElement("KEGG");
			e.addAttribute("id", key);
			e.addAttribute("function", usedKeggMap.get(key));
			keggRoot.add(e);
		}
		this.root.add(keggRoot);
		this.root.addAttribute("KEGG_count", String.valueOf(usedKeggMap.size()));

		Element goRoot = DocumentFactory.getInstance().createElement("GOs");
		for (String key : usedGoMap.keySet()) {
			Element e = DocumentFactory.getInstance().createElement("GO");
			GoObo goObo = usedGoMap.get(key);

			e.addAttribute("id", goObo.getId());
			e.addAttribute("name", goObo.getName());
			e.addAttribute("namespace", goObo.getNamespace());
			goRoot.add(e);
		}
		this.root.add(goRoot);
		this.root.addAttribute("GO_count", String.valueOf(usedGoMap.size()));

		Element ecRoot = DocumentFactory.getInstance().createElement("ECs");
		for (String key : usedEcMap.keySet()) {
			Element e = DocumentFactory.getInstance().createElement("EC");
			EnzymeCommission ec = usedEcMap.get(key);

			e.addAttribute("id", ec.getId());
			e.addAttribute("de", ec.getDe());
			e.addAttribute("ca", ec.getCa());

			ArrayList<String> anlist = ec.getAn();
			StringBuilder ansb = new StringBuilder();
			for (String an : anlist) {
				ansb.append(an).append(";");
			}
			if (ansb.length() > 0) {
				ansb.deleteCharAt(ansb.length() - 1);
			}
			e.addAttribute("an", ansb.toString());

			ecRoot.add(e);
		}
		this.root.add(ecRoot);
		this.root.addAttribute("EC_count", String.valueOf(usedEcMap.size()));

		for (MetaProteinAnnoEggNog mp : proteins) {
			this.addProtein(mp);
		}
	}

	private void addProtein(MetaProteinAnnoEggNog mpa) {

		Element ePro = DocumentFactory.getInstance().createElement("Protein");
		MetaProtein protein = mpa.getPro();
		ePro.addAttribute("Name", protein.getName());
		if (this.searchType.equals(MetaConstants.maxQuant)) {

			ePro.addAttribute("Group_ID", String.valueOf(protein.getGroupId()));
			ePro.addAttribute("Protein_ID", String.valueOf(protein.getProteinId()));
			ePro.addAttribute("Peptide_count", String.valueOf(protein.getPepCount()));
			ePro.addAttribute("Razor_unique_peptide_count", String.valueOf(protein.getRazorUniPepCount()));
			ePro.addAttribute("Unique_peptide_count", String.valueOf(protein.getUniPepCount()));
			ePro.addAttribute("E_value", String.valueOf(protein.getEvalue()));

			int[] ms2Counts = protein.getMs2Counts();
			StringBuilder ms2CountSb = new StringBuilder();
			for (int ms2c : ms2Counts) {
				ms2CountSb.append(ms2c).append("_");
			}
			if (ms2CountSb.length() > 0) {
				ms2CountSb.deleteCharAt(ms2CountSb.length() - 1);
			}
			ePro.addAttribute("ms2Counts", ms2CountSb.toString());

		} else if (this.searchType.equals(MetaConstants.openSearch)) {

			ePro.addAttribute("Group_ID", String.valueOf(protein.getGroupId()));
			ePro.addAttribute("Protein_ID", String.valueOf(protein.getProteinId()));
			ePro.addAttribute("Peptide_count", String.valueOf(protein.getPepCount()));
			ePro.addAttribute("PSM_count", String.valueOf(protein.getPsmCount()));
			ePro.addAttribute("Score", String.valueOf(protein.getScore()));

		} else if (this.searchType.equals(MetaConstants.pFind) || this.searchType.equals(MetaConstants.fragpipeIGC)) {

			ePro.addAttribute("Group_ID", String.valueOf(protein.getGroupId()));
			ePro.addAttribute("Protein_ID", String.valueOf(protein.getProteinId()));
			ePro.addAttribute("Peptide_count", String.valueOf(protein.getPepCount()));
			ePro.addAttribute("PSM_count", String.valueOf(protein.getPsmCount()));
			ePro.addAttribute("Score", String.valueOf(protein.getScore()));

		} else if (this.searchType.equals(MetaConstants.pFindHap) || this.searchType.equals(MetaConstants.pFindMag)
				|| this.searchType.equals(MetaConstants.fragpipeMAG)
				|| this.searchType.equals(MetaConstants.alphapeptMAG)) {

			ePro.addAttribute("Group_ID", String.valueOf(protein.getGroupId()));
			ePro.addAttribute("Protein_ID", String.valueOf(protein.getProteinId()));
			ePro.addAttribute("Peptide_count", String.valueOf(protein.getPepCount()));
			ePro.addAttribute("PSM_count", String.valueOf(protein.getPsmCount()));
			ePro.addAttribute("Score", String.valueOf(protein.getScore()));
		}

		double[] intensities = protein.getIntensities();
		StringBuilder intensitySb = new StringBuilder();
		for (double inten : intensities) {
			intensitySb.append(inten).append("_");
		}
		if (intensitySb.length() > 0) {
			intensitySb.deleteCharAt(intensitySb.length() - 1);
		}

		ePro.addAttribute("intensities", intensitySb.toString());

		if (searchType.equals(MetaConstants.pFindHap)) {
			ePro.addAttribute("Pro_name", ((MetaProteinAnnoHgm) mpa).getPro_name());
			ePro.addAttribute("seed", ((MetaProteinAnnoHgm) mpa).getSeed());
			ePro.addAttribute("KofamKOALA_KO", ((MetaProteinAnnoHgm) mpa).getKofamKOALA_KO());
			ePro.addAttribute("description", ((MetaProteinAnnoHgm) mpa).getDescription());
			ePro.addAttribute("pfam", ((MetaProteinAnnoHgm) mpa).getPfam());
		} else if (searchType.equals(MetaConstants.pFindMag) || searchType.equals(MetaConstants.fragpipeMAG)
				|| this.searchType.equals(MetaConstants.alphapeptMAG)) {
			ePro.addAttribute("Pro_name", ((MetaProteinAnnoMag) mpa).getPro_name());
			ePro.addAttribute("seed", ((MetaProteinAnnoMag) mpa).getSeed());
			ePro.addAttribute("description", ((MetaProteinAnnoMag) mpa).getDescription());
			ePro.addAttribute("pfam", ((MetaProteinAnnoMag) mpa).getPfam());
		} else if (searchType.equals(MetaConstants.fragpipeIGC)) {

		}

		ePro.addAttribute("Tax_id", mpa.tax_id);
		ePro.addAttribute("Best_tax_name", mpa.best_tax_level);
		ePro.addAttribute("Preferred_name", mpa.prefered_name);
		ePro.addAttribute("Gene_Ontology", mpa.getGene_Ontology());
		ePro.addAttribute("EC", mpa.getEC());
		ePro.addAttribute("KEGG_ko", mpa.getKEGG_ko());
		ePro.addAttribute("KEGG_Pathway", mpa.getKEGG_Pathway());
		ePro.addAttribute("KEGG_Module", mpa.getKEGG_Module());
		ePro.addAttribute("KEGG_Reaction", mpa.getKEGG_Reaction());
		ePro.addAttribute("KEGG_rclass", mpa.getKEGG_rclass());
		ePro.addAttribute("BRITE", mpa.getBRITE());
		ePro.addAttribute("KEGG_TC", mpa.getKEGG_TC());
		ePro.addAttribute("CAZy", mpa.getCAZy());
		ePro.addAttribute("BiGG_Reaction", mpa.getBiGG_Reaction());
		ePro.addAttribute("COG", mpa.getCog());
		ePro.addAttribute("NOG", mpa.getNog());

		this.proRoot.add(ePro);
	}
	
	public boolean addTaxon(Taxon taxon) {
		int id = taxon.getId();
		if (taxonIdSet.contains(id)) {
			return false;
		} else {
			taxonIdSet.add(id);

			String name = taxon.getName();
			int parentId = taxon.getParentId();
			int rankId = taxon.getRankId();
			int rootTypeId = taxon.getRootType().getId();

			Element eTaxon = DocumentFactory.getInstance().createElement("Taxon");
			eTaxon.addAttribute("id", String.valueOf(id));
			eTaxon.addAttribute("name", name);
			eTaxon.addAttribute("parentId", String.valueOf(parentId));
			eTaxon.addAttribute("rankId", String.valueOf(rankId));
			eTaxon.addAttribute("rootTypeId", String.valueOf(rootTypeId));

			int[] mainParentIds = taxon.getMainParentIds();

			if (mainParentIds != null) {
				StringBuilder mpsb = new StringBuilder();
				for (int mpid : mainParentIds) {
					mpsb.append(mpid).append(",");
				}
				if (mpsb.length() > 0) {
					mpsb.deleteCharAt(mpsb.length() - 1);
				}
				eTaxon.addAttribute("mainParentIds", mpsb.toString());
			}

			taxonRoot.add(eTaxon);
			return true;
		}
	}

	public void close() {
		try {
			this.writer.write(document);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			LOGGER.error("Error in writing proteins to " + output, e);
		}
		try {
			this.writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			LOGGER.error("Error in writing proteins to " + output, e);
		}
	}
}
