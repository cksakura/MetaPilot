/**
 * 
 */
package uk.ac.quadram.metapilot.mdb.unipept;

import java.io.File;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import uk.ac.quadram.metapilot.core.taxonomy.Taxon;
import uk.ac.quadram.metapilot.core.taxonomy.TaxonomyDatabase;
import uk.ac.quadram.metapilot.dbSearch.MaxQuant.MaxquantPepReader;
import uk.ac.quadram.metapilot.task.io.pep.MetaPeptide;
import uk.ac.quadram.metapilot.task.io.pep.AbstractMetaPeptideReader;
import uk.ac.quadram.metapilot.task.io.pep.MetaPeptideXMLWriter;
import uk.ac.quadram.metapilot.task.par.MetaConstants;

/**
 * @author Kai Cheng
 *
 */
public class UnipeptLCAGetter {

	private AbstractMetaPeptideReader reader;
	private MetaPeptideXMLWriter writer;
	private Appendable log;
	private File outFile;

	private TaxonomyDatabase taxDb;
	private MetaPeptide[] peps;

	private SimpleDateFormat format = new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss");

	private String taskName = "Unipept LCA calculation task";
	private static Logger LOGGER = LogManager.getLogger(UnipeptLCAGetter.class);

	public void getLCA() {

		ArrayList<String> seqlist = new ArrayList<String>();
		for (int i = 0; i < this.peps.length; i++) {
			if (peps[i].isUse()) {
				seqlist.add(peps[i].getSequence());
			}
		}
		String[] sequences = seqlist.toArray(new String[seqlist.size()]);

		HashMap<String, UnipLCAResult> lcaResultMap = UnipeptRequester.pept2lca(sequences);

		for (MetaPeptide metapep : this.peps) {
			if (!metapep.isUse()) {
				continue;
			}
			String sequence = metapep.getSequence();
			if (lcaResultMap.containsKey(sequence)) {
				UnipLCAResult lcaResult = lcaResultMap.get(sequence);
				int taxId = lcaResult.getTaxon().getId();
				Taxon lca = taxDb.getTaxonFromId(taxId);

				if (lca == null) {
					continue;
				}
				if (lca.getMainParentIds() == null) {
					lca.setMainParentIds(taxDb.getMainParentTaxonIds(taxId));
				}

				int[] mainTaxonIds = lca.getMainParentIds();

				for (int parentId : mainTaxonIds) {
					Taxon parentTaxon = this.taxDb.getTaxonFromId(parentId);
					if (parentTaxon != null) {
						if (parentTaxon.getMainParentIds() == null) {
							parentTaxon.setMainParentIds(taxDb.getMainParentTaxonIds(parentId));
						}
						writer.addTaxon(parentTaxon);
					}
				}

				metapep.setLca(lca);
				writer.addPeptide(metapep);
			}
		}

		try {
			log.append(format.format(new Date()) + "\t" + "LCA analysis by Unipept started" + "\n");
		} catch (IOException e) {
			// TODO Auto-generated catch block
			LOGGER.error("Error in writing log file", e);
		}

		writer.close();
	}

}
