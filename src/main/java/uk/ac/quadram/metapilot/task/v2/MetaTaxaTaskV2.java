/**
 * 
 */
package uk.ac.quadram.metapilot.task.v2;

import java.io.File;
import java.io.FileNotFoundException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.ExecutionException;

import javax.swing.JProgressBar;
import javax.swing.SwingWorker;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import uk.ac.quadram.metapilot.core.taxonomy.RootType;
import uk.ac.quadram.metapilot.core.taxonomy.Taxon;
import uk.ac.quadram.metapilot.core.taxonomy.TaxonomyDatabase;
import uk.ac.quadram.metapilot.core.taxonomy.TaxonomyRanks;
import uk.ac.quadram.metapilot.dbSearch.MaxQuant.MaxquantPepReader;
import javafx.scene.control.ProgressBar;
import uk.ac.quadram.metapilot.mdb.pep.Pep2TaxaDatabase;
import uk.ac.quadram.metapilot.mdb.unipept.UnipLCAResult;
import uk.ac.quadram.metapilot.mdb.unipept.UnipeptRequester;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.task.MetaAbstractTask;
import uk.ac.quadram.metapilot.task.MetaReportTask;
import uk.ac.quadram.metapilot.task.io.MetaAlgorithm;
import uk.ac.quadram.metapilot.task.io.MetaBiomJsonHandler;
import uk.ac.quadram.metapilot.task.io.MetaTreeHandler;
import uk.ac.quadram.metapilot.task.io.pep.MetaPeptide;
import uk.ac.quadram.metapilot.task.io.pep.MetaPeptideXMLReader;
import uk.ac.quadram.metapilot.task.io.pep.MetaPeptideXMLWriter;
import uk.ac.quadram.metapilot.task.par.MetaBaseParameter;
import uk.ac.quadram.metapilot.task.par.MetaConstants;
import uk.ac.quadram.metapilot.task.par.MetaPep2TaxaPar;
import uk.ac.quadram.metapilot.task.v2.par.MetaSourcesV2;
import uk.ac.quadram.metapilot.task.v2.par.MetaParameterMQ;

/**
 * @author Kai Cheng
 *
 */
public class MetaTaxaTaskV2 extends MetaAbstractTask {

	private MetaPep2TaxaPar mptp;
	private String taxDbString;
	private static Logger LOGGER = LogManager.getLogger(MetaTaxaTaskV2.class);
	private static final String taskName = "Taxonomy analysis";
	private static final SimpleDateFormat format = new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss");

	public MetaTaxaTaskV2(MetaBaseParameter metaPar, MetaResourceParameters mrp, ProgressBar progressBar) {
		super(metaPar, mrp, progressBar);
	}

	protected void initial() {
		this.mptp = ((MetaParameterMQ) metaPar).getPtp();
	}

	@Override
	protected Boolean call() throws Exception {
		// TODO Auto-generated method stub

		boolean finished1 = false;
		boolean finished2 = false;

		updateMessage(taskName);

		LOGGER.info(taskName + ": start");
		System.out.println(format.format(new Date()) + "\t" + taskName + ": start");

		File summaryMetaFile = metaPar.getMetafile();
		if (summaryMetaFile == null || !summaryMetaFile.exists()) {
			summaryMetaFile = new File(metaPar.getDbSearchResultFile(), "metainfo.tsv");
			try {
				metaPar.getMetadata().exportMetadata(summaryMetaFile);
			} catch (FileNotFoundException e) {
				// TODO Auto-generated catch block

				LOGGER.error(taskName + ": error in exporting metadata information", e);
				System.out.println(
						format.format(new Date()) + "\t" + taskName + ": error in exporting metadata information");

				return false;
			}

			metaPar.setMetafile(summaryMetaFile);
		}

		if (mptp.isUnipept()) {
			finished2 = unipeptCalculate();
		}

		updateProgress(85, 100);

		if (finished1 || finished2) {

			LOGGER.info(taskName + ": finished");
			System.out.println(format.format(new Date()) + "\t" + taskName + ": finished");

			return true;
		} else {

			LOGGER.info(taskName + ": failed");
			System.out.println(format.format(new Date()) + "\t" + taskName + ": failed");

			return false;
		}
	}

	private boolean unipeptCalculate() {

		boolean finished = false;

		File taxFile;
		if (metaPar.getDbSearchResultFile() != null) {
			taxFile = new File(metaPar.getDbSearchResultFile(), "taxonomy_analysis");
		} else {
			taxFile = new File(metaPar.getResult());
		}

		if (!taxFile.exists()) {
			taxFile.mkdir();
		}

		File taxResultFile = new File(taxFile, MetaAlgorithm.Unipept.getName() + ".taxonomy.xml");

		if (taxResultFile.exists() && taxResultFile.length() > 0) {

			LOGGER.info(taskName + ": taxonomy result file " + taxResultFile + " was found");

			MetaPeptideXMLReader reader = new MetaPeptideXMLReader(taxResultFile);

			File refinedTaxon = new File(taxFile, MetaAlgorithm.Unipept.getName() + ".taxa.refine.csv");
			reader.exportCsv(taxFile, MetaAlgorithm.Unipept, mptp.getLeastPepCount());

			MetaPeptide[] peps = reader.getPeptides();
			Taxon[] taxons = reader.getTaxons();
			String[] expNames = reader.getFileNames();

			File megan = new File(taxFile, MetaAlgorithm.Unipept.getName() + ".biom");
			if (!megan.exists() || megan.length() == 0) {
				MetaBiomJsonHandler.export(peps, taxons, expNames, megan.getAbsolutePath());
			}

			File iMetaLab = new File(taxFile, MetaAlgorithm.Unipept.getName() + ".iMetaLab.tree.csv");
			if (!iMetaLab.exists() || iMetaLab.length() == 0) {
				MetaTreeHandler.export(peps, taxons, expNames, iMetaLab);
			}

			File jsOutput = new File(metaPar.getReportDir() + "\\reports\\metamap\\data\\tree.js");
			if (!jsOutput.getParentFile().exists()) {
				jsOutput.getParentFile().mkdirs();
			}

			if (!jsOutput.exists() || jsOutput.length() == 0) {
				MetaTreeHandler.exportJS(peps, taxons, expNames, jsOutput);
			}

			if (megan.exists() && iMetaLab.exists() && jsOutput.exists()) {
				finished = true;
			}

			File reportDir = new File(metaPar.getReportDir() + "\\reports\\rmd_html");

			if (!reportDir.exists()) {
				reportDir.mkdirs();
			}

			MetaReportTask task = new MetaReportTask(metaPar.getThreadCount());

			File report_taxonomy_summary = new File(reportDir, "report_taxonomy_summary.html");
			if (!report_taxonomy_summary.exists() || report_taxonomy_summary.length() == 0) {
				if (metaPar.getMetafile() != null && metaPar.getMetafile().exists()) {
					task.addTask(MetaReportTask.taxon, refinedTaxon, report_taxonomy_summary, metaPar.getMetafile());
				} else {
					task.addTask(MetaReportTask.taxon, refinedTaxon, report_taxonomy_summary);
				}
			}

			task.run();

			try {
				finished = task.get();
			} catch (InterruptedException | ExecutionException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}

			return finished;

		} else {

			TaxonomyDatabase taxDb = new TaxonomyDatabase(taxDbString);
			String[] fileNames = null;

			String searchType = "";
			MetaPeptide[] peps = null;

			switch (metaPar.getWorkflowType()) {

			case MaxQuantWorkflow: {

				File maxquantTxt = new File(metaPar.getDbSearchResultFile() + "\\combined\\txt\\peptides.txt");

				if (maxquantTxt.exists()) {

					searchType = MetaConstants.maxQuant;
					MaxquantPepReader reader = new MaxquantPepReader(maxquantTxt, ((MetaParameterMQ) metaPar));
					fileNames = ((MaxquantPepReader) reader).getIntensityTitle();
					peps = reader.getMetaPeptides();

					MetaPeptideXMLWriter writer = new MetaPeptideXMLWriter(taxResultFile, searchType,
							((MetaParameterMQ) metaPar).getQuanMode(), fileNames);

					ArrayList<String> seqlist = new ArrayList<String>();
					for (int i = 0; i < peps.length; i++) {
						if (peps[i].isUse()) {
							seqlist.add(peps[i].getSequence());
						}
					}
					String[] sequences = seqlist.toArray(new String[seqlist.size()]);

					HashMap<String, UnipLCAResult> lcaResultMap = UnipeptRequester.pept2lca(sequences);

					for (MetaPeptide metapep : peps) {
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
								Taxon parentTaxon = taxDb.getTaxonFromId(parentId);
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

					writer.close();
				}

				break;
			}

			case TaxonomyAnalysis: {

				searchType = ((MetaParameterMQ) metaPar).getQuanResultType();
				fileNames = ((MetaParameterMQ) metaPar).getIntensityTitles();
				peps = ((MetaParameterMQ) metaPar).getPeps();

				MetaPeptideXMLWriter writer = new MetaPeptideXMLWriter(taxResultFile, searchType,
						((MetaParameterMQ) metaPar).getQuanMode(), fileNames);

				ArrayList<String> seqlist = new ArrayList<String>();
				for (int i = 0; i < peps.length; i++) {
					if (searchType.equals(MetaConstants.maxQuant) && !peps[i].isUse()) {
						continue;
					}
					seqlist.add(peps[i].getSequence());
				}
				String[] sequences = seqlist.toArray(new String[seqlist.size()]);

				HashMap<String, UnipLCAResult> lcaResultMap = UnipeptRequester.pept2lca(sequences);

				for (MetaPeptide metapep : peps) {
					if (searchType.equals(MetaConstants.maxQuant) && !metapep.isUse()) {
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
							Taxon parentTaxon = taxDb.getTaxonFromId(parentId);
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

				writer.close();

				break;
			}

			default:
				break;
			}
		}

		if (!taxResultFile.exists() || taxResultFile.length() == 0) {
			return false;
		}

		MetaPeptideXMLReader reader = new MetaPeptideXMLReader(taxResultFile);
		File refinedTaxon = new File(taxFile, MetaAlgorithm.Unipept.getName() + ".taxa.refine.csv");
		reader.exportCsv(taxFile, MetaAlgorithm.Unipept, mptp.getLeastPepCount());

		MetaPeptide[] peps = reader.getPeptides();
		Taxon[] taxons = reader.getTaxons();
		String[] expNames = reader.getFileNames();

		File megan = new File(taxFile, MetaAlgorithm.Unipept.getName() + ".biom");
		if (!megan.exists() || megan.length() == 0) {
			MetaBiomJsonHandler.export(peps, taxons, expNames, megan.getAbsolutePath());
		}

		File iMetaLab = new File(taxFile, MetaAlgorithm.Unipept.getName() + ".iMetaLab.tree.csv");
		if (!iMetaLab.exists() || iMetaLab.length() == 0) {
			MetaTreeHandler.export(peps, taxons, expNames, iMetaLab);
		}

		File jsOutput = new File(metaPar.getReportDir() + "\\reports\\metamap\\data\\tree.js");
		if (!jsOutput.exists() || jsOutput.length() == 0) {
			MetaTreeHandler.exportJS(peps, taxons, expNames, jsOutput);
		} else {
			jsOutput = new File(metaPar.getReportDir() + "\\reports\\metamap\\data\\tree2.js");
			MetaTreeHandler.exportJS(peps, taxons, expNames, jsOutput);
		}

		if (megan.exists() && iMetaLab.exists() && jsOutput.exists()) {
			finished = true;
		}

		File reportDir = new File(metaPar.getReportDir() + "\\reports\\rmd_html");

		if (!reportDir.exists()) {
			reportDir.mkdirs();
		}

		MetaReportTask task = new MetaReportTask(metaPar.getThreadCount());

		File report_taxonomy_summary = new File(reportDir, "report_taxonomy_summary.html");
		if (!report_taxonomy_summary.exists() || report_taxonomy_summary.length() == 0) {
			if (this.metaPar.getMetafile() != null && this.metaPar.getMetafile().exists()) {
				task.addTask(MetaReportTask.taxon, refinedTaxon, report_taxonomy_summary, metaPar.getMetafile());
			} else {
				task.addTask(MetaReportTask.taxon, refinedTaxon, report_taxonomy_summary);
			}
		}

		task.run();

		try {
			finished = task.get();
		} catch (InterruptedException | ExecutionException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		return finished;
	}

	@Override
	protected String getTaskName() {
		// TODO Auto-generated method stub
		return taskName;
	}

	@Override
	protected Logger getLogger() {
		// TODO Auto-generated method stub
		return LOGGER;
	}
}
