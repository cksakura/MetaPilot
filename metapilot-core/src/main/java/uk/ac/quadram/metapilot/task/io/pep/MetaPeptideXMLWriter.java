/**
 * 
 */
package uk.ac.quadram.metapilot.task.io.pep;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.HashSet;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.dom4j.Document;
import org.dom4j.DocumentFactory;
import org.dom4j.DocumentHelper;
import org.dom4j.Element;
import org.dom4j.io.OutputFormat;
import org.dom4j.io.XMLWriter;

import uk.ac.quadram.metapilot.core.taxonomy.Taxon;

/**
 * @author Kai Cheng
 *
 */
public class MetaPeptideXMLWriter {

	private File output;
	private Document document;
	private XMLWriter writer;
	private Element root;
	private Element peptideRoot;
	private Element taxonRoot;

	private int peptideCount = 0;
	private int taxonCount = 0;

	private HashSet<Integer> taxonIdSet;

	private static final Logger LOGGER = LogManager.getLogger(MetaPeptideXMLWriter.class);

	public MetaPeptideXMLWriter(String output, String searchEngine, String quanType, String... fileNames) {
		this(new File(output), searchEngine, quanType, fileNames);
	}

	public MetaPeptideXMLWriter(File output, String searchEngine, String quanType, String... fileNames) {

		this.output = output;
		this.document = DocumentHelper.createDocument();

		OutputStreamWriter bufferedWriter;
		try {
			bufferedWriter = new OutputStreamWriter(new FileOutputStream(output), "UTF8");
			this.writer = new XMLWriter(bufferedWriter, OutputFormat.createPrettyPrint());
		} catch (IOException e) {
			// TODO Auto-generated catch block
			LOGGER.error("Error in writing taxonomy information to " + output, e);
		}

		this.root = DocumentFactory.getInstance().createElement("Peptide_taxonomy");
		this.document.setRootElement(root);

		root.addAttribute("searchEngine", searchEngine);
		root.addAttribute("quantitativeMethodId", quanType);

		StringBuilder sb = new StringBuilder();
		for (String name : fileNames) {
			sb.append(name).append(";");
		}
		sb.deleteCharAt(sb.length() - 1);

		root.addAttribute("fileNames", sb.toString());

		peptideRoot = DocumentFactory.getInstance().createElement("Peptides");
		root.add(peptideRoot);

		taxonRoot = DocumentFactory.getInstance().createElement("Taxons");
		root.add(taxonRoot);

		this.taxonIdSet = new HashSet<Integer>();
	}

	public void addPeptide(MetaPeptide pep) {

		Element ePep = pep.getXmlPepElement();
		peptideRoot.add(ePep);
		peptideCount++;
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
			taxonCount++;
			return true;
		}
	}

	public void close() {
		root.addAttribute("peptideCount", String.valueOf(peptideCount));
		root.addAttribute("taxonCount", String.valueOf(taxonCount));
		try {
			writer.write(document);
			writer.close();

		} catch (IOException e) {
			// TODO Auto-generated catch block
			LOGGER.error("Error in writing taxonomy information to " + output, e);
		}
	}

}