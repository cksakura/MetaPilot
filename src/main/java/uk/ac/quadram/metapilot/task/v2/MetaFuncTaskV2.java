/**
 * 
 */
package uk.ac.quadram.metapilot.task.v2;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;

import javax.swing.JProgressBar;
import javax.swing.SwingWorker;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.dom4j.Document;
import org.dom4j.DocumentException;
import org.dom4j.Element;
import org.dom4j.io.SAXReader;

import uk.ac.quadram.metapilot.core.function.COGFinder;
import uk.ac.quadram.metapilot.core.function.CategoryFinder;
import uk.ac.quadram.metapilot.core.function.GOBPFinder;
import uk.ac.quadram.metapilot.core.function.GOCCFinder;
import uk.ac.quadram.metapilot.core.function.GOMFFinder;
import uk.ac.quadram.metapilot.core.function.KEGGFinder;
import uk.ac.quadram.metapilot.core.function.KOFinder;
import uk.ac.quadram.metapilot.core.function.NOGFinder;
import uk.ac.quadram.metapilot.core.function.sql.IGCFuncSqliteSearcher;
import uk.ac.quadram.metapilot.dbSearch.MaxQuant.MaxquantProReader;
import javafx.scene.control.ProgressBar;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.task.MetaAbstractTask;
import uk.ac.quadram.metapilot.task.MetaReportTask;
import uk.ac.quadram.metapilot.task.io.pro.MetaProtein;
import uk.ac.quadram.metapilot.task.io.pro.MetaProteinAnno1;
import uk.ac.quadram.metapilot.task.io.pro.MetaProteinXMLWriter;
import uk.ac.quadram.metapilot.task.io.pro.pfind.MetaProteinAnnoEggNog;
import uk.ac.quadram.metapilot.task.io.pro.pfind.MetaProteinXMLReader2;
import uk.ac.quadram.metapilot.task.io.pro.pfind.MetaProteinXMLWriter2;
import uk.ac.quadram.metapilot.task.par.MetaBaseParameter;
import uk.ac.quadram.metapilot.task.par.MetaConstants;
import uk.ac.quadram.metapilot.task.v2.par.MetaSourcesV2;
import uk.ac.quadram.metapilot.task.v2.par.MetaParameterMQ;

/**
 * @author Kai Cheng
 *
 */
public class MetaFuncTaskV2 extends MetaAbstractTask {

	private static Logger LOGGER = LogManager.getLogger(MetaFuncTaskV2.class);
	private static final String taskName = "Functional annotation";
	private String searchType;
	private String quanType;

	public MetaFuncTaskV2(MetaBaseParameter metaPar, MetaResourceParameters mrp, ProgressBar progressBar) {
		super(metaPar, mrp, progressBar);
	}


	
	@Override
	protected Boolean call() throws Exception {
		// TODO Auto-generated method stub

		boolean finished = false;

		updateMessage(taskName);

		LOGGER.info(taskName + ": start");
		System.out.println(format.format(new Date()) + "\t" + taskName + ": start");

		File funcResultFile = null;

		if (((MetaParameterMQ) metaPar).getQuickOutput() != null) {
			funcResultFile = new File(((MetaParameterMQ) metaPar).getQuickOutput());
			if (!funcResultFile.exists()) {
				funcResultFile.mkdir();
			}

		} else {
			File summaryMetaFile;
			if (metaPar.getDbSearchResultFile() != null) {
				funcResultFile = new File(metaPar.getDbSearchResultFile(), "functional_annotation");
				summaryMetaFile = new File(metaPar.getDbSearchResultFile(), "metainfo.tsv");
			} else {
				funcResultFile = new File(metaPar.getResult(), "functional_annotation");
				summaryMetaFile = new File(metaPar.getResult(), "metainfo.tsv");
			}

			if (!funcResultFile.exists()) {
				funcResultFile.mkdir();
			}

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

		File funcResultXml = new File(funcResultFile, "functions.xml");

		if (funcResultXml.exists() && funcResultXml.length() > 0) {

			LOGGER.info(taskName + ": functional annotation result file " + funcResultXml + " was found");
			System.out.println(format.format(new Date()) + "\t" + taskName + ": functional annotation result file "
					+ funcResultXml + " was found");

			updateProgress(95, 100);

			File funTsv = new File(funcResultFile, "functions.tsv");
			File html = new File(funcResultFile, "functions.html");

			SAXReader reader = new SAXReader();
			Document document = reader.read(funcResultXml);
			Element root = document.getRootElement();

			if (root.element("COGs") != null && root.element("NOGs") != null && root.element("KEGGs") != null
					&& root.element("ECs") != null && root.element("GOs") != null) {

				MetaProteinXMLReader2 funReader = new MetaProteinXMLReader2(funcResultXml);
				if (!funTsv.exists() || funTsv.length() == 0) {
					LOGGER.info(taskName + ": writing results to " + funTsv);
					System.out.println(format.format(new Date()) + "\t" + taskName + ": writing results to " + funTsv);

					try {
						funReader.exportTsv(funTsv);
					} catch (IOException e) {
						LOGGER.error(taskName + ": writing results to " + funTsv + " failed", e);
						System.err.println(format.format(new Date()) + "\t" + taskName + ": writing results to "
								+ funTsv + " failed");

						return false;
					}

					LOGGER.info(taskName + ": writing results to " + funTsv + " finished");
					System.out.println(format.format(new Date()) + "\t" + taskName + ": writing results to " + funTsv
							+ " finished");
				} else {
					LOGGER.info(taskName + ": result file " + funTsv + " already exist");
					System.out.println(
							format.format(new Date()) + "\t" + taskName + ": result file " + funTsv + " already exist");
				}

				updateProgress(96, 100);

				if (!html.exists() || html.length() == 0) {
					LOGGER.info(taskName + ": writing results to " + html);
					System.out.println(format.format(new Date()) + "\t" + taskName + ": writing results to " + html);

					try {
						funReader.exportHtml(html);
					} catch (IOException e) {
						LOGGER.error(taskName + ": writing results to " + html + " failed", e);
						System.err.println(format.format(new Date()) + "\t" + taskName + ": writing results to " + html
								+ " failed");

						return false;
					}

					LOGGER.info(taskName + ": writing results to " + html + " finished");
					System.out.println(
							format.format(new Date()) + "\t" + taskName + ": writing results to " + html + " finished");

				} else {
					LOGGER.info(taskName + ": result file " + funTsv + " already exist");
					System.out.println(
							format.format(new Date()) + "\t" + taskName + ": result file " + funTsv + " already exist");
				}
			}

			updateProgress(97, 100);

			if (funTsv.exists() && html.exists()) {
				finished = true;
			}

			File reportDir = new File(metaPar.getReportDir() + "\\reports\\rmd_html");

			if (!reportDir.exists()) {
				reportDir.mkdirs();
			}

			MetaReportTask task = new MetaReportTask(metaPar.getThreadCount());

			File report_function_summary = new File(reportDir, "report_function_summary.html");
			if (!report_function_summary.exists() || report_function_summary.length() == 0) {

				if (this.metaPar.getMetafile() != null && this.metaPar.getMetafile().exists()) {

					try {
						task.addTask(MetaReportTask.function, funTsv, report_function_summary, metaPar.getMetafile());
					} catch (Exception e) {
						LOGGER.error(taskName + ": writing report to " + report_function_summary + " failed", e);
						System.err.println(format.format(new Date()) + "\t" + taskName + ": writing report to "
								+ report_function_summary + " failed");

						return false;
					}

				} else {
					try {
						task.addTask(MetaReportTask.function, funTsv, report_function_summary);
					} catch (Exception e) {
						LOGGER.error(taskName + ": writing report to " + report_function_summary + " failed", e);
						System.err.println(format.format(new Date()) + "\t" + taskName + ": writing report to "
								+ report_function_summary + " failed");

						return false;
					}
				}
			} else {
				LOGGER.info(taskName + ": report file " + report_function_summary + " already exist");
				System.out.println(format.format(new Date()) + "\t" + taskName + ": report file "
						+ report_function_summary + " already exist");
			}

			task.run();

			finished = task.get();

			if (finished) {
				updateProgress(100, 100);
			}

			return finished;
		}

		MetaProtein[] metapros = null;
		String[] expTitles = null;

		switch (metaPar.getWorkflowType()) {
		case FunctionalAnnotation: {
			metapros = ((MetaParameterMQ) metaPar).getPros();
			expTitles = ((MetaParameterMQ) metaPar).getIntensityTitles();
			searchType = ((MetaParameterMQ) metaPar).getQuanResultType();
			quanType = ((MetaParameterMQ) metaPar).getQuanMode();
			break;
		}
		case MaxQuantWorkflow: {

			File maxquantTxt = new File(metaPar.getDbSearchResultFile() + "\\combined\\txt\\proteinGroups.txt");
			if (maxquantTxt.exists()) {
				MaxquantProReader reader = new MaxquantProReader(maxquantTxt, (MetaParameterMQ) metaPar);
				metapros = reader.getMetaProteins();
				expTitles = reader.getIntensityTitle();
				searchType = MetaConstants.maxQuant;
				quanType = reader.getQuanMode();

			} else {
				LOGGER.error(taskName + ": protein identification result was not found, task failed");
				System.err.println(format.format(new Date()) + "\t" + taskName
						+ ": protein identification result was not found, task failed");

				return false;
			}

			break;
		}

		case MSFraggerWorkflow: {

			File funTsv = new File(funcResultFile, "functions.tsv");
			File html = new File(funcResultFile, "functions.html");

			try {

				MetaProteinXMLReader2 funReader = new MetaProteinXMLReader2(funcResultXml);

				LOGGER.info(taskName + ": writing results to " + funTsv);
				System.out.println(format.format(new Date()) + "\t" + taskName + ": writing results to " + funTsv);

				try {
					funReader.exportTsv(funTsv);
				} catch (IOException e) {
					LOGGER.error(taskName + ": writing results to " + funTsv + " failed", e);
					System.err.println(
							format.format(new Date()) + "\t" + taskName + ": writing results to " + funTsv + " failed");

					return false;
				}

				LOGGER.info(taskName + ": writing results to " + funTsv + " finished");
				System.out.println(
						format.format(new Date()) + "\t" + taskName + ": writing results to " + funTsv + " finished");

				LOGGER.info(taskName + ": writing results to " + html);
				System.out.println(format.format(new Date()) + "\t" + taskName + ": writing results to " + html);

				try {
					funReader.exportHtml(html);
				} catch (IOException e) {
					LOGGER.error(taskName + ": writing results to " + html + " failed", e);
					System.err.println(
							format.format(new Date()) + "\t" + taskName + ": writing results to " + html + " failed");

					return false;
				}

				LOGGER.info(taskName + ": writing results to " + html + " finished");
				System.out.println(
						format.format(new Date()) + "\t" + taskName + ": writing results to " + html + " finished");

			} catch (DocumentException e) {
				LOGGER.error(taskName + ": error in reading function result file " + funcResultXml, e);
				System.err.println(format.format(new Date()) + "\t" + taskName
						+ ": error in reading function result file " + funcResultXml);

				return false;
					}

		updateProgress(100, 100);

		return finished;
		}

		default:
			break;
		}

		if (metapros == null) {
			LOGGER.error(taskName + ": no proteins were found");
			LOGGER.error(taskName + ": task failed");

			System.err.println(format.format(new Date()) + "\t" + taskName + ": no proteins were found");
			System.err.println(format.format(new Date()) + "\t" + taskName + ": task failed");

			return false;
		}

		updateProgress(86, 100);

		File funTsv = new File(funcResultFile, "functions.tsv");


		boolean match = this.funcAnnoSqlite(metapros, expTitles, funcResultXml);
		if (!match) {
			LOGGER.error(taskName + ": error in searching functional annotation database");
			System.err.println(format.format(new Date()) + "\t" + taskName
					+ ": error in searching functional annotation database");
			return false;
		}

		MetaProteinXMLReader2 funReader = new MetaProteinXMLReader2(funcResultXml);

		System.out.println(format.format(new Date()) + "\t" + taskName + ": writing results to " + funTsv);

		funReader.exportTsv(funTsv);

		File html = new File(funcResultFile, "functions.html");

		System.out.println(format.format(new Date()) + "\t" + taskName + ": writing results to " + html);

		funReader.exportHtml(html);

		if (funTsv.exists() && html.exists()) {
			finished = true;
		}
	

		File reportDir = new File(metaPar.getReportDir() + "\\reports\\rmd_html");

		if (!reportDir.exists()) {
			reportDir.mkdirs();
		}

		MetaReportTask task = new MetaReportTask(metaPar.getThreadCount());

		File report_function_summary = new File(reportDir, "report_function_summary.html");
		if (!report_function_summary.exists() || report_function_summary.length() == 0) {

			if (this.metaPar.getMetadata() == null) {
				task.addTask(MetaReportTask.function, funTsv, report_function_summary);
			} else {
				task.addTask(MetaReportTask.function, funTsv, report_function_summary, metaPar.getMetafile());
			}
		}

		task.run();

		finished = task.get();

		if (finished) {
			updateProgress(100, 100);
		}

		return finished;
	}

	private boolean funcAnnoSqlite(MetaProtein[] proteins, String[] fileNames, File final_pro_xml_file) {

		return true;
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
