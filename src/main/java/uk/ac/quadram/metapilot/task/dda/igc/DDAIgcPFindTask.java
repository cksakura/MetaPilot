package uk.ac.quadram.metapilot.task.dda.igc;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileFilter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.nio.file.StandardWatchEventKinds;
import java.nio.file.WatchEvent;
import java.nio.file.WatchKey;
import java.nio.file.WatchService;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.ExecutionException;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import uk.ac.quadram.metapilot.core.mod.PostTransModification;
import uk.ac.quadram.metapilot.core.taxonomy.Taxon;
import uk.ac.quadram.metapilot.core.taxonomy.TaxonomyDatabase;
import uk.ac.quadram.metapilot.core.tools.FormatTool;
import uk.ac.quadram.metapilot.dbSearch.dbReducer.DBReducerTask;
import uk.ac.quadram.metapilot.dbSearch.pfind.PFindIsobaricTask;
import uk.ac.quadram.metapilot.dbSearch.pfind.PFindLabelFreeTask;
import uk.ac.quadram.metapilot.dbSearch.pfind.PFindPeptide;
import uk.ac.quadram.metapilot.dbSearch.pfind.PFindProtein;
import uk.ac.quadram.metapilot.dbSearch.pfind.PFindResultReader;
import uk.ac.quadram.metapilot.dbSearch.pfind.PFindTask;
import uk.ac.quadram.metapilot.dbSearch.pfind.PFindResultReader.PeptideResult;
import javafx.scene.control.ProgressBar;
import uk.ac.quadram.metapilot.mdb.unipept.UnipLCAResult;
import uk.ac.quadram.metapilot.mdb.unipept.UnipeptRequester;
import uk.ac.quadram.metapilot.quant.flashLFQ.FlashLfqTask;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters.ToolType;
import uk.ac.quadram.metapilot.task.MetaPilotTask;
import uk.ac.quadram.metapilot.task.MetaPythonReportCopyTask;
import uk.ac.quadram.metapilot.task.MetaPythonReportTask;
import uk.ac.quadram.metapilot.task.dda.PFindSpectraConvertTask;
import uk.ac.quadram.metapilot.task.dda.igc.par.DDAIgcPFindParameters;
import uk.ac.quadram.metapilot.task.io.MetaAlgorithm;
import uk.ac.quadram.metapilot.task.io.MetaBiomJsonHandler;
import uk.ac.quadram.metapilot.task.io.MetaTreeHandler;
import uk.ac.quadram.metapilot.task.io.pep.MetaPeptide;
import uk.ac.quadram.metapilot.task.io.pep.MetaPeptideXMLReader;
import uk.ac.quadram.metapilot.task.io.pep.MetaPeptideXMLWriter;
import uk.ac.quadram.metapilot.task.io.pro.MetaProtein;
import uk.ac.quadram.metapilot.task.par.MetaConstants;
import uk.ac.quadram.metapilot.task.par.MetaPep2TaxaPar;
import uk.ac.quadram.metapilot.task.pfind.MetaDbReducerTask;

import uk.ac.quadram.metapilot.task.v2.par.MetaOpenPar;
import uk.ac.quadram.metapilot.task.par.MetaData;

public class DDAIgcPFindTask extends MetaPilotTask {

	private static final Logger LOGGER = LogManager.getLogger(DDAIgcPFindTask.class);
	private static final String taskName = "DDA IGC database search by pFind workflow";

	protected File rawDir;
	protected File quan_pep_file;
	protected File quan_pro_file;
	protected HashMap<String, String> expNameMap;
	protected PFindResultReader pFindResultReader;
	protected HashMap<String, PostTransModification> ptmMap;

	public DDAIgcPFindTask(DDAIgcPFindParameters metaPar, MetaResourceParameters mrp, ProgressBar progressBar,
			ProgressBar secondBar) {
		super(metaPar, mrp, progressBar, secondBar);
		
		// Debug logging for metadata
		LOGGER.info("{}: Constructor called with parameters", taskName);
		if (metaPar != null) {
			LOGGER.info("{}: Parameters object is not null", taskName);
			if (metaPar.getMetadata() != null) {
				MetaData metadata = metaPar.getMetadata();
				String[] rawFiles = metadata.getRawFiles();
				String[] expNames = metadata.getExpNames();
				LOGGER.info("{}: Metadata contains {} raw files", taskName, rawFiles != null ? rawFiles.length : 0);
				if (rawFiles != null) {
					for (int i = 0; i < rawFiles.length; i++) {
						LOGGER.info("{}: Raw file {}: {} -> {}", taskName, i, rawFiles[i], expNames != null && i < expNames.length ? expNames[i] : "N/A");
					}
				}
			} else {
				LOGGER.warn("{}: Metadata is null in parameters", taskName);
			}
		} else {
			LOGGER.error("{}: Parameters object is null!", taskName);
		}
	}

	@Override
	public void forceStop() {
		// TODO Auto-generated method stub

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

	@Override
	protected Boolean call() throws Exception {
		// Update progress at start
		updateProgress(0, 100);
		updateMessage("Starting DDA IGC pFind workflow...");

		updateMessage(getTaskName());

		LOGGER.info(getTaskName() + ": start");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": start");

		// Update progress before spectra conversion
		updateProgress(5, 100);
		updateMessage("Starting spectra conversion...");

		PFindSpectraConvertTask spTask = new PFindSpectraConvertTask(metaPar, mrp, progressBar);
		LOGGER.info("{}: PFindSpectraConvertTask created, executing...", taskName);
		spTask.execute();
		try {
			boolean spectraConverted = spTask.get(); // Wait for completion and get result
			if (!spectraConverted) {
				LOGGER.error(getTaskName() + ": error in pFind spectrum convert workflow");
				System.err.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": error in pFind spectrum convert workflow");
				return false;
			}
		} catch (InterruptedException | ExecutionException e) {
			LOGGER.error(getTaskName() + ": error in pFind spectrum convert workflow", e);
			System.err.println(format.format(new Date()) + "\t" + getTaskName() + ": error in pFind spectrum convert workflow: " + e.getMessage());
			return false;
		}

		// Update progress after spectra conversion
		updateProgress(20, 100);
		updateMessage("Spectra conversion completed, starting database reduction...");

		LOGGER.info(getTaskName() + ": pFind spectrum convert workflow: finished");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": pFind spectrum convert workflow: finished");

		// Update progress before database reduction
		updateProgress(25, 100);
		updateMessage("Creating sample-specific database...");

		MetaDbReducerTask dbrTask = new MetaDbReducerTask(metaPar, mrp, secondBar);
		LOGGER.info("{}: MetaDbReducerTask created, executing...", taskName);
		dbrTask.execute();
		try {
			boolean dbReduced = dbrTask.get(); // Wait for completion and get result
			if (!dbReduced) {
				LOGGER.error(getTaskName() + ": error in sample-specific database generation");
				System.err.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": error in sample-specific database generation");
				return false;
			}
		} catch (InterruptedException | ExecutionException e) {
			LOGGER.error(getTaskName() + ": error in sample-specific database generation", e);
			System.err.println(format.format(new Date()) + "\t" + getTaskName() + ": error in sample-specific database generation: " + e.getMessage());
			return false;
		}

		// Update progress after database reduction
		updateProgress(40, 100);
		updateMessage("Database reduction completed, starting identification...");

		LOGGER.info(getTaskName() + ": sample specific database generation: finished");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": sample specific database generation: finished");

		File pfindSearchFile = new File(metaPar.getResult(), "pfind_search");
		if (!pfindSearchFile.exists()) {
			pfindSearchFile.mkdir();
			LOGGER.info("  Created pFind search directory: " + pfindSearchFile.getAbsolutePath());
		} else {
			LOGGER.info("  pFind search directory already exists: " + pfindSearchFile.getAbsolutePath());
		}

		metaPar.setDbSearchResultFile(pfindSearchFile);
		
		boolean iden = iden(spTask.getMgfFiles(), spTask.getPf2Files());
		if (!iden) {
			LOGGER.error(getTaskName() + ": error in peptide and protein identification");
			System.err.println(
					format.format(new Date()) + "\t" + getTaskName() + ": error in peptide and protein identification");
			return false;
		}

		// Update progress after identification
		updateProgress(60, 100);
		updateMessage("Identification completed, starting quantification...");

		updateProgress(70, 100);

		String quanMode = metaPar.getQuantitativeMode();
		boolean quant = false;
		if (quanMode.equals(MetaConstants.labelFree)) {
			quant = lfQuant();
			if (!quant) {
				LOGGER.error(getTaskName() + ": error in peptide and protein quantification");
				System.err.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": error in peptide and protein quantification");
				return false;
			}
		} else if (quanMode.equals(MetaConstants.isobaricLabel)) {
			quant = lfQuant();
			if (!quant) {
				LOGGER.error(getTaskName() + ": error in peptide and protein quantification");
				System.err.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": error in peptide and protein quantification");
				return false;
			}

			updateProgress(75, 100);

			quant = isobaricQuant();
			if (!quant) {
				LOGGER.error(getTaskName() + ": error in peptide and protein quantification");
				System.err.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": error in peptide and protein quantification");
				return false;
			}
		}

		// Update progress before report generation
		updateProgress(80, 100);
		updateMessage("Quantification completed, generating reports...");

		updateProgress(80, 100);

		exportReport();

		// Update progress to completion
		updateProgress(100, 100);
		updateMessage("Workflow completed successfully!");

		updateProgress(100, 100);

		return true;
	}

	protected boolean iden(String[] mgfFiles, String[] pf2Files) {

		updateMessage(getTaskName() + ": peptide and protein identification");

		LOGGER.info(getTaskName() + ": peptide and protein identification started");
		System.out.println(
				format.format(new Date()) + "\t" + getTaskName() + ": peptide and protein identification started");

		// Update progress at start of identification
		updateProgress(45, 100);
		updateMessage("Starting peptide and protein identification...");

		File pFindBin = (new File(mrp.getToolPath(ToolType.PFIND))).getParentFile();
		PFindTask pFindTask = new PFindTask(pFindBin, metaPar.getThreadCount());
		HashMap<String, PostTransModification> ptmMap = pFindTask.getPtmMap();
		MetaOpenPar mosp = ((DDAIgcPFindParameters) metaPar).getMop();

		String ms2Type = ((DDAIgcPFindParameters) metaPar).getMs2ScanMode();
		boolean useDBReducer = false;
		File dbReducerExe = new File(mrp.getToolPath(ToolType.DBREDUCER));
		DBReducerTask dbReducerTask = null;
		if (dbReducerExe.exists() && ms2Type.equals(MetaConstants.HCD_FTMS)
				&& ((DDAIgcPFindParameters) metaPar).getQuanMode().equals(MetaConstants.labelFree)) {
			useDBReducer = true;
			dbReducerTask = new DBReducerTask(dbReducerExe);
		}

		PFindResultReader pFindResultReader = new PFindResultReader(metaPar.getDbSearchResultFile(),
				metaPar.getSpectraFile(), ptmMap);

		if (pFindResultReader.hasResult()) {

			LOGGER.info(getTaskName() + ": peptide and protein identification results have been found in "
					+ metaPar.getDbSearchResultFile() + ", go to next step");
			System.out.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": peptide and protein identification results have been found in "
					+ metaPar.getDbSearchResultFile() + ", go to next step");
			
			// Update progress if results already exist
			updateProgress(55, 100);
			updateMessage("Identification results found, proceeding...");
		} else {
			// Update progress before database search
			updateProgress(50, 100);
			updateMessage("Performing database search...");
			String db = metaPar.getCurrentDb();
			if (useDBReducer && dbReducerTask != null) {
				dbReducerTask.setTotalThread(metaPar.getThreadCount());
				File dbReducerFolder = new File(metaPar.getDbSearchResultFile().getAbsolutePath(), "DBReducer");
				dbReducerFolder.mkdir();

				dbReducerTask.searchMgf(mgfFiles, dbReducerFolder.getAbsolutePath(), db, ms2Type, mosp.getOpenPsmFDR(),
						((DDAIgcPFindParameters) this.metaPar).getFixMods(),
						((DDAIgcPFindParameters) this.metaPar).getVariMods());

				dbReducerTask.run(1, mgfFiles.length * 48);

				File dbReducerFile = new File(dbReducerFolder, "DBReducer.fasta");
				if (dbReducerFile.exists() && dbReducerFile.length() > 0) {

					db = dbReducerFile.getAbsolutePath();

					LOGGER.info(getTaskName() + ": reduced database was generated in " + dbReducerFile);
					System.out.println(format.format(new Date()) + "\t" + getTaskName()
							+ ": reduced database was generated in " + dbReducerFile);
				} else {
					LOGGER.info(getTaskName() + ": reduced database was not found");
					System.out.println(
							format.format(new Date()) + "\t" + getTaskName() + ": reduced database was not found");
				}
			}

			pFindTask.setTotalThread(metaPar.getThreadCount());

			pFindTask.searchPF2(pf2Files, metaPar.getDbSearchResultFile().getAbsolutePath(), db, ms2Type,
					((DDAIgcPFindParameters) metaPar).isOpenSearch(), true, mosp.getOpenPsmFDR(),
					((DDAIgcPFindParameters) this.metaPar).getFixMods(), ((DDAIgcPFindParameters) this.metaPar).getVariMods(),
					((DDAIgcPFindParameters) this.metaPar).getEnzyme(),
					((DDAIgcPFindParameters) this.metaPar).getMissCleavages(),
					((DDAIgcPFindParameters) this.metaPar).getDigestMode(),
					((DDAIgcPFindParameters) this.metaPar).getIsobaricTag());

			pFindTask.run(1, pf2Files.length * 48);

		}

		if (!pFindResultReader.hasResult()) {
			LOGGER.error(getTaskName() + ": peptide and protein identification failed");
			System.err.println(
					format.format(new Date()) + "\t" + getTaskName() + ": peptide and protein identification failed");
			return false;
		}

		// Update progress before reading results
		updateProgress(55, 100);
		updateMessage("Reading identification results...");

		PFindTask.cleanPFindResultFolder(metaPar.getDbSearchResultFile(), true);

		try {
			pFindResultReader.read();
		} catch (NumberFormatException | IOException e) {
			// TODO Auto-generated catch block
			LOGGER.error(getTaskName() + ": error in reading the database search result files");
			System.err.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": error in reading the database search result files");
			return false;
		}

		// Update progress after identification completion
		updateProgress(60, 100);
		updateMessage("Identification completed successfully!");

		LOGGER.info(getTaskName() + ": peptide and protein identification finished");
		System.out.println(
				format.format(new Date()) + "\t" + getTaskName() + ": peptide and protein identification finished");

		return true;
	}

	protected boolean lfQuant() {
		updateMessage(getTaskName() + ": peptide and protein quantification");

		LOGGER.info(getTaskName() + ": peptide and protein label-free quantification started");
		System.out.println(format.format(new Date()) + "\t" + getTaskName()
				+ ": peptide and protein label-free quantification started");

		// Update progress at start of label-free quantification
		updateProgress(65, 100);
		updateMessage("Starting label-free quantification...");

		String flashLFQ = mrp.getToolPath(ToolType.FLASHLFQ);
		FlashLfqTask flashLfqTask = new FlashLfqTask(flashLFQ, metaPar.getThreadCount());

		PFindLabelFreeTask pfindLFTask = new PFindLabelFreeTask(flashLfqTask, metaPar.getMetadata(),
				metaPar.getDbSearchResultFile(), rawDir, metaPar.getSpectraFile());
		boolean quan = pfindLFTask.run();

		if (!quan) {

			LOGGER.info(getTaskName() + ": error in peptide and protein label-free quantification");
			System.out.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": error in peptide and protein label-free quantification");

			return false;
		}

		this.quan_pep_file = pfindLFTask.getQuan_pep_file();
		this.quan_pro_file = pfindLFTask.getQuan_pro_file_razor();

		// Update progress after label-free quantification completion
		updateProgress(70, 100);
		updateMessage("Label-free quantification completed!");

		LOGGER.info(getTaskName() + ": peptide and protein quantification label-free finished");
		System.out.println(format.format(new Date()) + "\t" + getTaskName()
				+ ": peptide and protein label-free quantification finished");

		return true;
	}

	protected boolean isobaricQuant() {
		updateMessage(getTaskName() + ": peptide and protein quantification");

		LOGGER.info(getTaskName() + ": peptide and protein isobaric quantification started");
		System.out.println(format.format(new Date()) + "\t" + getTaskName()
				+ ": peptide and protein isobaric quantification started");

		// Update progress at start of isobaric quantification
		updateProgress(75, 100);
		updateMessage("Starting isobaric quantification...");

		File pFindBin = (new File(mrp.getToolPath(ToolType.PFIND))).getParentFile();
		PFindTask pFindTask = new PFindTask(pFindBin, metaPar.getThreadCount());

		PFindResultReader reader = new PFindResultReader(metaPar.getDbSearchResultFile(), metaPar.getSpectraFile(),
				pFindTask.getPtmMap());
		PFindIsobaricTask pfindIsobaricTask = new PFindIsobaricTask(reader, this.metaPar.getSpectraFile(),
				metaPar.getDbSearchResultFile(), ((DDAIgcPFindParameters) this.metaPar).getIsobaricTag(),
				((DDAIgcPFindParameters) this.metaPar).getIsobaric(),
				((DDAIgcPFindParameters) this.metaPar).getIsoCorFactor(),
				((DDAIgcPFindParameters) this.metaPar).isCombineLabel(), metaPar.getMetadata());

		try {
			pfindIsobaricTask.run();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			LOGGER.error(getTaskName() + ": error in quantifying PSMs", e);
			System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": error in quantifying PSMs");
		}

		this.quan_pep_file = pfindIsobaricTask.getQuanPepFile();
		this.quan_pro_file = pfindIsobaricTask.getQuanRazorProFile();

		if (quan_pep_file.exists() && quan_pro_file.exists()) {
			LOGGER.info(getTaskName() + ": quantification results have been fount in " + this.metaPar.getSpectraFile());
			System.out.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": quantification results have been fount in " + this.metaPar.getSpectraFile());
			
			// Update progress after isobaric quantification completion
			updateProgress(80, 100);
			updateMessage("Isobaric quantification completed!");
			
			return true;
		}

		if (!this.quan_pep_file.exists()) {
			return false;
		}

		if (!this.quan_pro_file.exists()) {
			return false;
		}

		LOGGER.info(getTaskName() + ": peptide and protein isobaric quantification finished");
		System.out.println(format.format(new Date()) + "\t" + getTaskName()
				+ ": peptide and protein isobaric quantification finished");

		return true;
	}

	protected boolean exportReport() {

		updateMessage(getTaskName() + ": exporting report");

		LOGGER.info(getTaskName() + ": exporting report started");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": exporting report started");

		File reportDir = new File(this.metaPar.getResult(), "report");
		if (!reportDir.exists()) {
			reportDir.mkdirs();
		}

		MetaPythonReportCopyTask reportCopyTask = new MetaPythonReportCopyTask(new File(metaPar.getResult()));
		Thread reportCopyThread = new Thread(reportCopyTask);
		reportCopyThread.start();
		try {
			reportCopyThread.join();
			if (!reportCopyTask.get()) {
				LOGGER.error(getTaskName() + ": error in copying report template");
				System.err.println(
						format.format(new Date()) + "\t" + getTaskName() + ": error in copying report template");
				return false;
			}
		} catch (InterruptedException | ExecutionException e) {
			LOGGER.error(getTaskName() + ": error in copying report template", e);
			System.err.println(format.format(new Date()) + "\t" + getTaskName() + ": error in copying report template: "
					+ e.getMessage());
			return false;
		}

		metaPar.setReportDir(reportDir);

		File reportHtmlDir = new File(reportDir + "\\reports\\rmd_html");

		if (!reportDir.exists()) {
			reportDir.mkdirs();
		}

		File summaryMetaFile = new File(this.metaPar.getResult(), "metainfo.tsv");
		try {
			this.metaPar.getMetadata().exportMetadata(summaryMetaFile);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		MetaPythonReportTask task = new MetaPythonReportTask(metaPar.getThreadCount());

		boolean summayrFinished = this.exportSummary(task);

		updateProgress(82, 100);

		String[] rawFiles = this.metaPar.getMetadata().getRawFiles();
		String[] quanExpNames = new String[rawFiles.length];
		for (int i = 0; i < rawFiles.length; i++) {
			quanExpNames[i] = rawFiles[i].substring(rawFiles[i].lastIndexOf("\\") + 1, rawFiles[i].lastIndexOf("."));
		}
		PFindResultReader pFindResultReader = new PFindResultReader(metaPar.getDbSearchResultFile(),
				metaPar.getSpectraFile(), ptmMap);
		try {
			pFindResultReader.read();
		} catch (NumberFormatException | IOException e) {
			// TODO Auto-generated catch block
			LOGGER.error(getTaskName() + ": error in reading the database search result files");
			System.err.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": error in reading the database search result files");
			return false;
		}

		this.final_pep_txt = new File(this.metaPar.getDbSearchResultFile(), "final_peptides.tsv");

		if (!final_pep_txt.exists() || final_pep_txt.length() == 0) {

			LOGGER.info(getTaskName() + ": exporting peptide report started");
			System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": exporting peptide report started");

			try (BufferedReader quanPepReader = new BufferedReader(new FileReader(this.quan_pep_file));
					PrintWriter writer = new PrintWriter(this.final_pep_txt)) {

				String line = quanPepReader.readLine();
				String[] title = line.split("\t");
				int seqId = -1;

				int[] intenId = new int[quanExpNames.length];
				int[] idenId = new int[quanExpNames.length];

				for (int i = 0; i < title.length; i++) {
					if (title[i].equals("Sequence")) {
						seqId = i;
					} else if (title[i].startsWith("Intensity_")) {
						for (int j = 0; j < quanExpNames.length; j++) {
							if (quanExpNames[j].equals(title[i].substring("Intensity_".length()))) {
								intenId[j] = i;
							}
						}
					} else if (title[i].startsWith("Detection Type_")) {
						for (int j = 0; j < quanExpNames.length; j++) {
							if (quanExpNames[j].equals(title[i].substring("Detection Type_".length()))) {
								idenId[j] = i;
							}
						}
					}
				}

				HashMap<String, double[]> pepIntensityMap = new HashMap<String, double[]>();
				HashMap<String, boolean[]> idenTypeMap = new HashMap<String, boolean[]>();

				while ((line = quanPepReader.readLine()) != null) {
					String[] cs = line.split("\t");
					double[] intensity = new double[quanExpNames.length];
					for (int i = 0; i < intensity.length; i++) {
						intensity[i] = Double.parseDouble(cs[intenId[i]]);
					}

					if (pepIntensityMap.containsKey(cs[seqId])) {
						double[] total = pepIntensityMap.get(cs[seqId]);
						for (int i = 0; i < intensity.length; i++) {
							total[i] += intensity[i];
						}
						pepIntensityMap.put(cs[seqId], total);
					} else {
						pepIntensityMap.put(cs[seqId], intensity);
					}

					boolean[] idenType = new boolean[quanExpNames.length];
					for (int i = 0; i < idenType.length; i++) {
						if (cs[idenId[i]].equals("MSMS")) {
							idenType[i] = true;
						} else {
							idenType[i] = false;
						}
					}
					idenTypeMap.put(cs[seqId], idenType);
				}
				quanPepReader.close();

				StringBuilder titlesb = new StringBuilder();
				titlesb.append("Sequence").append("\t");
				titlesb.append("Base sequence").append("\t");
				titlesb.append("Length").append("\t");
				titlesb.append("Missed cleavages").append("\t");
				titlesb.append("Mass").append("\t");
				titlesb.append("Proteins").append("\t");
				titlesb.append("Charges").append("\t");
				titlesb.append("Score").append("\t");
				titlesb.append("Final score").append("\t");
				titlesb.append("Q-value").append("\t");

				for (int i = 0; i < quanExpNames.length; i++) {
					if (this.expNameMap.containsKey(quanExpNames[i])) {
						titlesb.append("Identification type " + expNameMap.get(quanExpNames[i])).append("\t");
					} else {
						titlesb.append("Identification type " + quanExpNames[i]).append("\t");
					}
				}

				titlesb.append("Intensity").append("\t");

				for (int i = 0; i < quanExpNames.length; i++) {
					if (this.expNameMap.containsKey(quanExpNames[i])) {
						titlesb.append("Intensity " + expNameMap.get(quanExpNames[i])).append("\t");
					} else {
						titlesb.append("Intensity " + quanExpNames[i]).append("\t");
					}
				}

				titlesb.append("Reverse").append("\t");
				titlesb.append("Potential contaminant").append("\t");
				titlesb.append("id").append("\t");
				titlesb.append("Protein group IDs").append("\t");
				titlesb.append("MS/MS Count");

				writer.println(titlesb);

				PeptideResult[] pepresult = pFindResultReader.getPeptideResults();
				HashMap<String, PFindProtein> proteinMap = pFindResultReader.getProteinMap();

				for (int i = 0; i < pepresult.length; i++) {

					String sequence = pepresult[i].getSequence();
					StringBuilder sb = new StringBuilder();
					sb.append(sequence).append("\t");
					sb.append(pepresult[i].getBaseSeq()).append("\t");
					sb.append(pepresult[i].getLength()).append("\t");
					sb.append(pepresult[i].getMiss()).append("\t");
					sb.append(pepresult[i].getMass()).append("\t");

					String[] proteins = pepresult[i].getProteins();
					if (proteins.length > 0) {
						for (String pro : proteins) {
							sb.append(pro).append(";");
						}
						sb.deleteCharAt(sb.length() - 1);
					}
					sb.append("\t");

					int[] charge = pepresult[i].getCharge();
					StringBuilder chargesb = new StringBuilder();
					for (int j = 0; j < charge.length; j++) {
						if (charge[j] == 1) {
							chargesb.append(j + 1).append(";");
						}
					}
					if (chargesb.length() > 1) {
						sb.append(chargesb.subSequence(0, chargesb.length() - 1)).append("\t");
					} else {
						sb.append("2").append("\t");
					}

					sb.append(pepresult[i].getRawScore()).append("\t");
					sb.append(pepresult[i].getFinalScore()).append("\t");
					sb.append(pepresult[i].getQValue()).append("\t");

					if (idenTypeMap.containsKey(sequence)) {
						boolean[] idenType = idenTypeMap.get(sequence);
						for (int j = 0; j < idenType.length; j++) {
							if (idenType[j]) {
								sb.append("By MS/MS").append("\t");
							} else {
								sb.append("\t");
							}
						}
					} else {
						for (int j = 0; j < quanExpNames.length; j++) {
							sb.append("\t");
						}
					}

					double totalIntensity = 0;

					if (pepIntensityMap.containsKey(sequence)) {
						double[] intensity = pepIntensityMap.get(sequence);
						for (int j = 0; j < intensity.length; j++) {
							totalIntensity += intensity[j];
						}

						sb.append((int) totalIntensity).append("\t");

						for (int j = 0; j < intensity.length; j++) {
							sb.append((int) intensity[j]).append("\t");
						}
					} else {
						sb.append("0\t");

						for (int j = 0; j < quanExpNames.length; j++) {
							sb.append("0\t");
						}
					}

					if (pepresult[i].isTarget()) {
						sb.append("").append("\t");
					} else {
						sb.append("+").append("\t");
					}

					sb.append("").append("\t");

					sb.append(i).append("\t");

					if (proteins.length > 0) {
						for (String pro : proteins) {
							if (proteinMap.containsKey(pro)) {
								sb.append(proteinMap.get(pro).getId()).append(";");
							}
						}
						sb.deleteCharAt(sb.length() - 1);
					}
					sb.append("\t");

					sb.append(pepresult[i].getMs2Count());

					writer.println(sb);
				}

				writer.close();

			} catch (IOException e) {
				// TODO Auto-generated catch block
				LOGGER.error(getTaskName() + ": error in reading peptides from " + this.quan_pep_file.getName(), e);
				System.err.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": error in reading peptides from " + this.quan_pep_file.getName());
			}

			LOGGER.info(getTaskName() + ": peptide report has been exported to " + final_pep_txt);
			System.out.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": peptide report has been exported to " + final_pep_txt);

			if (!final_pep_txt.exists() || final_pep_txt.length() == 0) {
				LOGGER.error(getTaskName() + ": error in writing the peptides to " + this.final_pep_txt);
				System.err.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": error in writing the peptides to " + this.final_pep_txt);
				return false;
			}
		}

		boolean pepReportOk = false;
		if (!final_pep_txt.exists() || final_pep_txt.length() == 0) {
			return false;
		} else {
			File report_peptides_summary = new File(reportHtmlDir, "report_peptides_summary.html");
			if (!report_peptides_summary.exists() || report_peptides_summary.length() == 0) {

				if (summaryMetaFile.exists()) {
					pepReportOk = task.addTask(MetaPythonReportTask.PEPTIDE, final_pep_txt, report_peptides_summary,
							summaryMetaFile);
				} else {
					pepReportOk = task.addTask(MetaPythonReportTask.PEPTIDE, final_pep_txt, report_peptides_summary);
				}
			}
		}

		updateProgress(84, 100);

		this.final_pro_txt = new File(this.metaPar.getResult(), "final_proteins.tsv");

		if (!final_pro_txt.exists() || final_pro_txt.length() == 0) {

			LOGGER.info(getTaskName() + ": exporting protein report started");
			System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": exporting protein report started");

			try (BufferedReader quanProReader = new BufferedReader(new FileReader(this.quan_pro_file));
					PrintWriter writer = new PrintWriter(this.final_pro_txt)) {

				String line = quanProReader.readLine();
				String[] title = line.split("\t");
				int proId = -1;
				ArrayList<String> tlist = new ArrayList<String>();
				for (int i = 0; i < title.length; i++) {
					if (title[i].equals("Protein Groups")) {
						proId = i;
					} else if (title[i].startsWith("Intensity")) {
						tlist.add(title[i]);
					}
				}

				String[] expNames = tlist.toArray(new String[tlist.size()]);

				StringBuilder titlesb = new StringBuilder();
				titlesb.append("Protein IDs").append("\t");
				titlesb.append("Majority protein IDs").append("\t");
				titlesb.append("Peptide counts (all)").append("\t");
				titlesb.append("Peptide counts (razor)").append("\t");
				titlesb.append("Number of proteins").append("\t");
				titlesb.append("Peptides").append("\t");
				titlesb.append("Score").append("\t");
				titlesb.append("Intensity").append("\t");

				int[] fileIds = new int[expNames.length];

				if (((DDAIgcPFindParameters) this.metaPar).getQuanMode().equals(MetaConstants.labelFree)) {

					for (int i = 0; i < expNames.length; i++) {
						for (int j = 0; j < title.length; j++) {
							if (expNames[i].equals(title[j])) {
								fileIds[i] = j;
							}
						}
					}

					for (int i = 0; i < expNames.length; i++) {
						String name = expNames[i].substring(expNames[i].indexOf("_") + 1);
						if (this.expNameMap.containsKey(name)) {
							titlesb.append("Intensity " + expNameMap.get(name)).append("\t");
						} else {
							titlesb.append("Intensity " + name).append("\t");
						}
					}

				} else if (((DDAIgcPFindParameters) this.metaPar).getQuanMode().equals(MetaConstants.isobaricLabel)) {

					for (int i = 0; i < expNames.length; i++) {
						for (int j = 0; j < title.length; j++) {
							if (expNames[i].equals(title[j])) {
								fileIds[i] = j;
							}
						}
					}

					for (int i = 0; i < expNames.length; i++) {
						String name = expNames[i].substring(expNames[i].indexOf("_") + 1);
						titlesb.append("Intensity " + name).append("\t");
					}
				}

				titlesb.append("Reverse").append("\t");
				titlesb.append("Potential contaminant").append("\t");
				titlesb.append("id").append("\t");
				titlesb.append("Peptide IDs").append("\t");
				titlesb.append("Peptide is razor");

				writer.println(titlesb);

				PFindProtein[] proteins = this.pFindResultReader.getProteins();
				HashMap<String, HashSet<String>> razorProPepMap = this.pFindResultReader.getRazorProPepMap();
				HashMap<String, HashSet<String>> allProPepMap = this.pFindResultReader.getProPepMap();
				PeptideResult[] peptideResults = this.pFindResultReader.getPeptideResults();

				HashMap<String, double[]> proIntensityMap = new HashMap<String, double[]>();
				for (int i = 0; i < proteins.length; i++) {
					proIntensityMap.put(proteins[i].getName(), null);
				}

				while ((line = quanProReader.readLine()) != null) {

					String[] cs = line.split("\t");

					if (proIntensityMap.containsKey(cs[proId])) {
						double[] intensity = new double[expNames.length];
						for (int i = 0; i < fileIds.length; i++) {
							intensity[i] = Double.parseDouble(cs[fileIds[i]]);
						}
						proIntensityMap.put(cs[proId], intensity);
					}
				}

				quanProReader.close();

				DecimalFormat df = FormatTool.getDFE4();
				for (int i = 0; i < proteins.length; i++) {

					StringBuilder sb = new StringBuilder();
					String pro = proteins[i].getName();
					double[] intensity = proIntensityMap.get(pro);

					if (intensity == null) {
						intensity = new double[expNames.length];
						Arrays.fill(intensity, 0.0);
					}

					HashSet<String> sameSet = proteins[i].getSameSet();
					HashSet<String> subSet = proteins[i].getSubSet();
					HashSet<Integer> pepIdSet = proteins[i].getPepIdSet();

					Integer[] pepIds = pepIdSet.toArray(new Integer[pepIdSet.size()]);
					Arrays.sort(pepIds);
					boolean[] isRazorPep = new boolean[pepIds.length];

					int proCount = 1 + sameSet.size() + subSet.size();
					int[] pepCountAll = new int[proCount];
					int[] pepCountRazor = new int[proCount];

					for (int j = 0; j < pepIds.length; j++) {
						String sequence = peptideResults[pepIds[j]].getBaseSeq();
						pepCountAll[0]++;
						if (razorProPepMap.get(pro).contains(sequence)) {
							pepCountRazor[0]++;
							isRazorPep[j] = true;
						}
					}

					int proNameId = 1;
					StringBuilder sameSb = new StringBuilder();
					for (String p : sameSet) {
						sameSb.append(";").append(p);
						for (int j = 0; j < pepIds.length; j++) {
							String sequence = peptideResults[pepIds[j]].getBaseSeq();
							if (allProPepMap.get(p).contains(sequence)) {
								pepCountAll[proNameId]++;
								if (isRazorPep[j]) {
									pepCountRazor[proNameId]++;
								}
							}
						}
						proNameId++;
					}

					StringBuilder subSb = new StringBuilder();
					for (String p : subSet) {
						subSb.append(";").append(p);
						for (int j = 0; j < pepIds.length; j++) {
							String sequence = peptideResults[pepIds[j]].getBaseSeq();
							if (allProPepMap.get(p).contains(sequence)) {
								pepCountAll[proNameId]++;
								if (isRazorPep[j]) {
									pepCountRazor[proNameId]++;
								}
							}
						}
						proNameId++;
					}

					sb.append(pro).append(sameSb).append(subSb).append("\t");
					sb.append(pro).append(sameSb).append("\t");
					for (int j = 0; j < pepCountAll.length; j++) {
						sb.append(pepCountAll[j]).append(";");
					}
					sb.append("\t");

					for (int j = 0; j < pepCountRazor.length; j++) {
						sb.append(pepCountRazor[j]).append(";");
					}
					sb.append("\t");

					sb.append(proCount).append("\t");
					sb.append(pepIds.length).append("\t");

					sb.append(FormatTool.getDF2().format(proteins[i].getScore())).append("\t");

					double totalIntensity = 0;

					for (int j = 0; j < intensity.length; j++) {
						totalIntensity += intensity[j];
					}

					sb.append(df.format(totalIntensity)).append("\t");
					for (int j = 0; j < intensity.length; j++) {
						sb.append(df.format(intensity[j])).append("\t");
					}

					if (pro.startsWith(PFindResultReader.REV)) {
						sb.append("+").append("\t");
					} else {
						sb.append("\t");
					}
					sb.append("\t");
					sb.append(proteins[i].getId()).append("\t");

					for (int j = 0; j < pepIds.length; j++) {
						sb.append(pepIds[j]).append(";");
					}
					sb.append("\t");

					for (int j = 0; j < isRazorPep.length; j++) {
						sb.append(isRazorPep[j]).append(";");
					}

					// currently no reverse proteins are exported
					if (!pro.startsWith(PFindResultReader.REV)) {
						writer.println(sb);
					}
				}

				writer.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				LOGGER.error(
						getTaskName() + ": error in reading quantified proteins from " + this.quan_pro_file.getName(),
						e);
				System.err.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": error in reading quantified proteins from " + this.quan_pro_file.getName());
			}

			if (final_pro_txt == null || final_pro_txt.length() == 0) {
				LOGGER.error(getTaskName() + ": error in writing final protein result to " + this.final_pro_txt);
				System.err.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": error in writing final protein result to " + this.final_pro_txt);
			}

			LOGGER.info(getTaskName() + ": protein report has been exported to " + final_pro_txt.getName());
			System.out.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": protein report has been exported to " + final_pro_txt.getName());
		}

		boolean proReportOk = false;
		if (!final_pro_txt.exists() || final_pro_txt.length() == 0) {
			return false;
		} else {
			File report_proteinGroups_summary = new File(reportHtmlDir, "report_proteinGroups_summary.html");
			if (!report_proteinGroups_summary.exists() || report_proteinGroups_summary.length() == 0) {
				if (summaryMetaFile.exists()) {
					proReportOk = task.addTask(MetaPythonReportTask.PROTEIN, final_pro_txt, report_proteinGroups_summary,
							summaryMetaFile);
				} else {
					proReportOk = task.addTask(MetaPythonReportTask.PROTEIN, final_pro_txt, report_proteinGroups_summary);
				}
			}
		}

		updateProgress(86, 100);

		LOGGER.info(getTaskName() + ": taxonomy analysis started");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": taxonomy analysis started");

		boolean finish = false;

		File taxFile = new File(metaPar.getDbSearchResultFile(), "taxonomy_analysis");
		if (!taxFile.exists()) {
			taxFile.mkdir();
		}

		HashMap<String, double[]> pepIntensityMap = new HashMap<String, double[]>();
		String[] expNames = null;
		try (BufferedReader quanPepReader = new BufferedReader(new FileReader(this.quan_pep_file))) {

			String line = quanPepReader.readLine();
			String[] title = line.split("\t");
			int seqId = -1;
			ArrayList<String> list = new ArrayList<String>();
			for (int i = 0; i < title.length; i++) {
				if (title[i].equals("Sequence")) {
					seqId = i;
				} else if (title[i].startsWith("Intensity_")) {
					list.add(title[i]);
				}
			}

			expNames = list.toArray(new String[list.size()]);

			if (((DDAIgcPFindParameters) this.metaPar).getQuanMode().equals(MetaConstants.labelFree)) {
				Arrays.sort(expNames);
			}

			int[] intenId = new int[list.size()];
			for (int i = 0; i < expNames.length; i++) {
				for (int j = 0; j < title.length; j++) {
					if (expNames[i].equals(title[j])) {
						intenId[i] = j;
					}
				}
			}

			while ((line = quanPepReader.readLine()) != null) {
				String[] cs = line.split("\t");
				double[] intensity = new double[expNames.length];
				for (int i = 0; i < intensity.length; i++) {
					intensity[i] = Double.parseDouble(cs[intenId[i]]);
				}
				pepIntensityMap.put(cs[seqId], intensity);
			}
			quanPepReader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			LOGGER.error("Error in reading peptides from " + this.quan_pep_file, e);
			System.err.println(format.format(new Date()) + "\t" + getTaskName() + ": error in reading peptides from "
					+ this.quan_pep_file);
			return false;
		}

		boolean executeUnipept = false;
		MetaPep2TaxaPar mptp = ((DDAIgcPFindParameters) metaPar).getPtp();
		File taxResultFile = new File(taxFile, MetaAlgorithm.Unipept.getName() + ".taxonomy.xml");
		if (!taxResultFile.exists() || taxResultFile.length() == 0) {
			executeUnipept = true;
		} else {

			MetaPeptideXMLReader reader = new MetaPeptideXMLReader(taxResultFile);

			File refinedTaxon = new File(taxFile, MetaAlgorithm.Builtin.getName() + ".taxa.refine.csv");
			reader.exportCsv(taxFile, MetaAlgorithm.Unipept, mptp.getLeastPepCount());

			MetaPeptide[] peps = reader.getPeptides();
			Taxon[] taxons = reader.getTaxons();

			File megan = new File(taxFile, MetaAlgorithm.Unipept.getName() + ".biom");
			if (!megan.exists() || megan.length() == 0) {
				MetaBiomJsonHandler.export(peps, taxons, expNames, megan.getAbsolutePath());
			}

			File allPepTaxa = new File(taxFile, MetaAlgorithm.Unipept.getName() + ".allPepTaxa.csv");
			if (!allPepTaxa.exists() || allPepTaxa.length() == 0) {
				reader.exportPeptideTaxaAll(allPepTaxa);
			}

			File iMetaLab = new File(taxFile, MetaAlgorithm.Unipept.getName() + ".iMetaLab.tree.csv");
			if (!iMetaLab.exists() || iMetaLab.length() == 0) {
				MetaTreeHandler.export(peps, taxons, expNames, iMetaLab);
			}

			File jsOutput = new File(metaPar.getReportDir() + "\\reports\\metamap\\data\\tree.js");
			if (!jsOutput.exists() || jsOutput.length() == 0) {
				MetaTreeHandler.exportJS(peps, taxons, expNames, iMetaLab);
			} else {
				jsOutput = new File(metaPar.getReportDir() + "\\reports\\metamap\\data\\tree2.js");
				MetaTreeHandler.exportJS(peps, taxons, expNames, jsOutput);
			}

			if (megan.exists() && allPepTaxa.exists() && iMetaLab.exists() && jsOutput.exists()) {
				finish = true;
			}

			File report_taxonomy_summary = new File(reportHtmlDir, "report_taxonomy_summary.html");
			if (!report_taxonomy_summary.exists() || report_taxonomy_summary.length() == 0) {
				if (report_taxonomy_summary.exists()) {
					finish = task.addTask(MetaPythonReportTask.TAXONOMY, refinedTaxon, report_taxonomy_summary,
							report_taxonomy_summary);

				} else {
					finish = task.addTask(MetaPythonReportTask.TAXONOMY, refinedTaxon, report_taxonomy_summary);
				}
			}
		}

		updateProgress(88, 100);

		HashMap<String, ArrayList<PeptideResult>> baseSeqMap = this.pFindResultReader.getBasePeptideMap();
		TaxonomyDatabase td = new TaxonomyDatabase(mrp.getToolPath(ToolType.FUNCDB));
		HashSet<String> sequenceSet = new HashSet<String>();
		sequenceSet.addAll(baseSeqMap.keySet());

		updateProgress(90, 100);

		if (executeUnipept) {

			String[] sequences = baseSeqMap.keySet().toArray(new String[baseSeqMap.size()]);
			Arrays.sort(sequences);

			HashMap<String, UnipLCAResult> lcaResultMap = UnipeptRequester.pept2lca(sequences);

			MetaPeptideXMLWriter writer = new MetaPeptideXMLWriter(taxResultFile, MetaConstants.pFind,
					MetaConstants.labelFree, expNames);

			for (int i = 0; i < sequences.length; i++) {

				ArrayList<PFindResultReader.PeptideResult> list = baseSeqMap.get(sequences[i]);

				String[] modseqs = new String[list.size()];
				double[] scores = new double[list.size()];
				int[] psmCounts = new int[list.size()];
				double[][] intensities = new double[list.size()][];

				int miss = 0;
				double mass = Double.MAX_VALUE;
				double score = 0;

				HashSet<String> proset = new HashSet<String>();
				for (int j = 0; j < list.size(); j++) {
					PeptideResult pr = list.get(j);
					modseqs[j] = pr.getSequence();
					scores[j] = pr.getRawScore();
					psmCounts[j] = pr.getMs2Count();

					intensities[j] = new double[expNames.length];
					if (pepIntensityMap.containsKey(pr.getSequence())) {
						intensities[j] = pepIntensityMap.get(pr.getSequence());
					}

					if (pr.getMiss() > miss) {
						miss = pr.getMiss();
					}
					if (pr.getMass() < mass) {
						mass = pr.getMass();
					}
					if (pr.getRawScore() > score) {
						score = pr.getRawScore();
					}

					String[] pros = pr.getProteins();
					for (int k = 0; k < pros.length; k++) {
						proset.add(pros[k]);
					}
				}
				String pro = "";

				StringBuilder prosb = new StringBuilder();
				for (String p : proset) {
					prosb.append(p).append(";");
				}
				if (prosb.length() > 0) {
					prosb = prosb.deleteCharAt(prosb.length() - 1);
				}
				pro = prosb.toString();

				if (lcaResultMap.containsKey(sequences[i])) {
					UnipLCAResult lcaResult = lcaResultMap.get(sequences[i]);
					int taxId = lcaResult.getTaxon().getId();
					Taxon lca = td.getTaxonFromId(taxId);

					if (lca != null) {

						PFindPeptide op = new PFindPeptide(sequences[i], miss, score, pro, modseqs, scores, psmCounts,
								intensities);

						op.setLca(lca);

						if (lca.getMainParentIds() == null) {
							lca.setMainParentIds(td.getMainParentTaxonIds(taxId));
						}

						int[] mainTaxonIds = lca.getMainParentIds();

						for (int parentId : mainTaxonIds) {
							Taxon parentTaxon = td.getTaxonFromId(parentId);
							if (parentTaxon != null) {
								if (parentTaxon.getMainParentIds() == null) {
									parentTaxon.setMainParentIds(td.getMainParentTaxonIds(parentId));
								}
								writer.addTaxon(parentTaxon);
							}
						}

						op.setLca(lca);
						writer.addPeptide(op);
					}
				}
			}
			writer.close();

			if (taxResultFile.exists() && taxResultFile.length() > 0) {

				MetaPeptideXMLReader reader = new MetaPeptideXMLReader(taxResultFile);

				File refinedTaxon = new File(taxFile, MetaAlgorithm.Builtin.getName() + ".taxa.refine.csv");
				reader.exportCsv(taxFile, MetaAlgorithm.Unipept, mptp.getLeastPepCount());

				MetaPeptide[] peps = reader.getPeptides();
				Taxon[] taxons = reader.getTaxons();

				File megan = new File(taxFile, MetaAlgorithm.Unipept.getName() + ".biom");
				if (!megan.exists() || megan.length() == 0) {
					MetaBiomJsonHandler.export(peps, taxons, expNames, megan.getAbsolutePath());
				}

				File allPepTaxa = new File(taxFile, MetaAlgorithm.Unipept.getName() + ".allPepTaxa.csv");
				if (!allPepTaxa.exists() || allPepTaxa.length() == 0) {
					reader.exportPeptideTaxaAll(allPepTaxa);
				}

				File iMetaLab = new File(taxFile, MetaAlgorithm.Unipept.getName() + ".iMetaLab.tree.csv");
				if (!iMetaLab.exists() || iMetaLab.length() == 0) {
					MetaTreeHandler.export(peps, taxons, expNames, iMetaLab);
				}

				File jsOutput = new File(metaPar.getReportDir() + "\\reports\\metamap\\data\\tree.js");
				if (!jsOutput.exists() || jsOutput.length() == 0) {
					MetaTreeHandler.exportJS(peps, taxons, expNames, iMetaLab);
				} else {
					jsOutput = new File(metaPar.getReportDir() + "\\reports\\metamap\\data\\tree2.js");
					MetaTreeHandler.exportJS(peps, taxons, expNames, jsOutput);
				}

				if (megan.exists() && allPepTaxa.exists() && iMetaLab.exists() && jsOutput.exists()) {
					finish = true;
				}

				File report_taxonomy_summary = new File(reportHtmlDir, "report_taxonomy_summary.html");
				if (!report_taxonomy_summary.exists() || report_taxonomy_summary.length() == 0) {
					if (summaryMetaFile.exists()) {
						finish = task.addTask(MetaPythonReportTask.TAXONOMY, refinedTaxon, report_taxonomy_summary,
								summaryMetaFile);

					} else {
						finish = task.addTask(MetaPythonReportTask.TAXONOMY, refinedTaxon, report_taxonomy_summary);
					}
				}
			}
		}

		if (finish) {
			LOGGER.info(getTaskName() + ": taxonomy analysis finished");
			System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": taxonomy analysis finished");
		} else {
			LOGGER.error(getTaskName() + ": taxonomy analysis failed");
			System.err.println(format.format(new Date()) + "\t" + getTaskName() + ": taxonomy analysis failed");
			return false;
		}

		LOGGER.info(getTaskName() + ": functional annotation started");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": functional annotation started");

		File funcFile = new File(metaPar.getDbSearchResultFile(), "functional_annotation");
		if (!funcFile.exists()) {
			funcFile.mkdir();
		}

		File final_pro_xml_file = new File(funcFile, "functions.xml");

		if (final_pro_xml_file.exists() && final_pro_xml_file.length() > 0) {

			LOGGER.info(getTaskName() + ": proteins with functional annotations have been exported to "
					+ final_pro_xml_file.getName());
			System.out.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": proteins with functional annotations have been exported to " + final_pro_xml_file.getName());
		} else {

			PFindProtein[] proteins = this.pFindResultReader.getProteins();
			HashMap<String, PFindProtein> promap = new HashMap<String, PFindProtein>();
			for (int i = 0; i < proteins.length; i++) {
				String pro = proteins[i].getName();
				promap.put(pro, proteins[i]);
			}

			ArrayList<MetaProtein> list = new ArrayList<MetaProtein>();
			String[] fileNames = null;
			try (BufferedReader quanProReader = new BufferedReader(new FileReader(this.final_pro_txt))) {

				String line = quanProReader.readLine();
				String[] title = line.split("\t");
				int id = -1;
				int proId = -1;

				ArrayList<String> tlist = new ArrayList<String>();
				for (int i = 0; i < title.length; i++) {
					if (title[i].startsWith("Intensity ")) {
						tlist.add(title[i]);
					} else if (title[i].equals("id")) {
						id = i;
					} else if (title[i].equals("Majority protein IDs")) {
						proId = i;
					}
				}

				fileNames = tlist.toArray(new String[tlist.size()]);

				int[] fileIds = new int[fileNames.length];
				for (int i = 0; i < fileNames.length; i++) {
					for (int j = 0; j < title.length; j++) {
						if (fileNames[i].equals(title[j])) {
							fileIds[i] = j;
						}
					}
				}

				while ((line = quanProReader.readLine()) != null) {

					String[] cs = line.split("\t");
					String[] pros = cs[proId].split(";");

					double[] intensity = new double[fileIds.length];
					for (int j = 0; j < fileIds.length; j++) {
						intensity[j] = Double.parseDouble(cs[fileIds[j]]);
					}

					if (promap.containsKey(pros[0])) {
						PFindProtein protein = promap.get(pros[0]);
						for (int i = 0; i < pros.length; i++) {
							MetaProtein mp = new MetaProtein(Integer.parseInt(cs[id]), i + 1, pros[i],
									protein.getPepCount(), protein.getSpCount(), protein.getScore(), intensity);
							list.add(mp);
						}
					}
				}
				quanProReader.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				LOGGER.error(getTaskName() + ": error in reading quantified proteins from " + this.quan_pro_file, e);
				System.err.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": error in reading quantified proteins from " + this.quan_pro_file);
			}

			updateProgress(94, 100);
		}

		LOGGER.info(getTaskName() + ": proteins with functional annotations have been exported to "
				+ final_pro_xml_file.getName());
		System.out.println(format.format(new Date()) + "\t" + getTaskName()
				+ ": proteins with functional annotations have been exported to " + final_pro_xml_file.getName());

		updateProgress(96, 100);

		File funcResultFile = final_pro_xml_file.getParentFile();
		File funTsv = new File(funcResultFile, "functions.tsv");
		File funReportTsv = new File(funcResultFile, "functions_report.tsv");
		File report_function_summary = new File(reportHtmlDir, "report_function_summary.html");
		if (!report_function_summary.exists() || report_function_summary.length() == 0) {
			if (summaryMetaFile.exists()) {
				finish = task.addTask(MetaPythonReportTask.FUNCTION, funReportTsv, report_function_summary, summaryMetaFile);
			} else {
				finish = task.addTask(MetaPythonReportTask.FUNCTION, funReportTsv, report_function_summary);
			}
		}

		LOGGER.info(getTaskName() + ": functional annotation finished");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": functional annotation finished");

		if (!summayrFinished || !pepReportOk || !proReportOk) {
			updateProgress(99, 100);

			LOGGER.info(getTaskName() + ": error in generating the reports");
			System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": error in generating the reports");

			return false;
		}

		task.run();

		try {
			if (task.get()) {

				LOGGER.info(getTaskName() + ": report generation successes");
				System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": report generation successes");

				updateProgress(100, 100);
				return true;
			}

		} catch (InterruptedException | ExecutionException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		return false;
	}

	protected boolean exportSummary(MetaPythonReportTask task) {

		File final_summary_txt = new File(this.metaPar.getResult(), "final_summary.tsv");
		try {

			PrintWriter countWriter = new PrintWriter(final_summary_txt);
			StringBuilder titlesb = new StringBuilder();
			titlesb.append("Raw file").append("\t");
			titlesb.append("Experiment").append("\t");
			titlesb.append("MS/MS").append("\t");
			titlesb.append("MS/MS Identified").append("\t");
			titlesb.append("MS/MS Identified [%]").append("\t");
			titlesb.append("Peptide Sequences Identified");

			countWriter.println(titlesb);

			int totalms2count = 0;
			int totalpsmCount = 0;
			int totalUnipepCount = 0;

			BufferedReader idReader = new BufferedReader(
					new FileReader(new File(this.metaPar.getDbSearchResultFile(), "id.tsv")));
			HashMap<String, int[]> idMap = new HashMap<String, int[]>();
			String line = idReader.readLine();
			while ((line = idReader.readLine()) != null) {
				String[] cs = line.split("\t");
				idMap.put(cs[0], new int[] { Integer.parseInt(cs[1]), Integer.parseInt(cs[2]) });
			}
			idReader.close();

			File spectraDirFile = metaPar.getSpectraFile();
			File rtFile = new File(spectraDirFile, "rt.txt");
			File ms2CountFile = new File(spectraDirFile, "ms2Count.txt");

			BufferedReader ms2CountReader = new BufferedReader(new FileReader(ms2CountFile));
			line = ms2CountReader.readLine();

			while ((line = ms2CountReader.readLine()) != null) {
				String[] cs = line.split("\t");
				String rawName = cs[0];

				if (idMap.containsKey(rawName)) {

					StringBuilder sb = new StringBuilder();
					sb.append(cs[0]).append("\t");

					if (this.expNameMap.containsKey(cs[0])) {
						sb.append(expNameMap.get(cs[0])).append("\t");
					} else {
						sb.append(cs[0]).append("\t");
					}

					int ms2count = Integer.parseInt(cs[1]);

					int[] idCounts = idMap.get(rawName);

					sb.append(ms2count).append("\t");
					sb.append(idCounts[1]).append("\t");

					double ratio = 0;
					if (ms2count > 0) {
						ratio = (double) idCounts[1] / (double) ms2count * 100.0;
					}

					sb.append(FormatTool.getDF2().format(ratio)).append("\t");

					sb.append(idCounts[0]);

					countWriter.println(sb);

					totalms2count += ms2count;
					totalpsmCount += idCounts[1];
					totalUnipepCount += idCounts[0];
				}
			}

			ms2CountReader.close();

			StringBuilder sb = new StringBuilder();
			sb.append("Total").append("\t");
			sb.append("\t");
			sb.append(totalms2count).append("\t");
			sb.append(totalpsmCount).append("\t");

			double ratio = 0;
			if (totalms2count > 0) {
				ratio = (double) totalpsmCount / (double) totalms2count * 100.0;
			}

			sb.append(FormatTool.getDF2().format(ratio)).append("\t");

			sb.append(totalUnipepCount);

			countWriter.print(sb);

			countWriter.close();

		} catch (IOException e) {
			// TODO Auto-generated catch block
			LOGGER.error(getTaskName() + ": error in writing PSM information to " + final_summary_txt.getName(), e);
			System.err.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": error in writing PSM information to " + final_summary_txt.getName());
		}
		File summaryMetaFile = new File(this.metaPar.getDbSearchResultFile(), "metainfo.tsv");
		try {
			this.metaPar.getMetadata().exportMetadata(summaryMetaFile);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		if (!final_summary_txt.exists() || final_summary_txt.length() == 0) {
			return false;

		} else {
			boolean finish = false;
			File report_ID_summary = new File(metaPar.getReportDir(), "report_ID_summary.html");
			if (summaryMetaFile != null && summaryMetaFile.exists()) {
				finish = task.addTask(MetaPythonReportTask.ID_SUMMARY, final_summary_txt, report_ID_summary, summaryMetaFile);
			} else {
				finish = task.addTask(MetaPythonReportTask.ID_SUMMARY, final_summary_txt, report_ID_summary);
			}

			return finish;
		}
	}

	// The inner class that implements the watch service loop
	public class WatchServiceLoop implements Runnable {

		private WatchService watchService;
		private FileFilter filter;
		private int totalCount = 0;
		private int currentCount = 0;
		private int barStart;
		private int barEnd;

		public WatchServiceLoop(WatchService watchService, FileFilter filter, int totalCount, int currentCount,
				int barStart, int barEnd) {
			this.watchService = watchService;
			this.filter = filter;
			this.totalCount = totalCount;
			this.currentCount = currentCount;
			this.barStart = barStart;
			this.barEnd = barEnd;
		}

		@SuppressWarnings("unchecked")
		@Override
		public void run() {

			while (currentCount < totalCount) {
				try {
					// Wait for a watch key to be available
					WatchKey key = watchService.take();
					// Iterate over the watch events
					for (WatchEvent<?> event : key.pollEvents()) {
						// Get the event kind and context
						WatchEvent.Kind<?> kind = event.kind();
						WatchEvent<Path> pathEvent = (WatchEvent<Path>) event;
						Path fileName = pathEvent.context();
						// If the event is a creation of a psm.tsv file, increment the count and update
						// the progress bar
						if (kind == StandardWatchEventKinds.ENTRY_CREATE && filter.accept(fileName.toFile())) {
							currentCount++;
							int percentage = (int) ((double) currentCount / (double) totalCount
									* (double) (barEnd - barStart));

							DDAIgcPFindTask.this.updateProgress(barStart + percentage, 100);

							System.out.println(format.format(new Date()) + "\t" + getTaskName()
									+ ": WatchService find one finished task at " + fileName.toFile());

							LOGGER.info(
									getTaskName() + ": WatchService find one finished task at " + fileName.toFile());
						}
					}
					// Reset the key and continue the loop
					key.reset();
				} catch (InterruptedException e) {
					System.out.println(format.format(new Date()) + "\t" + getTaskName()
							+ ": error in create a WatchService for separate search");

					LOGGER.error(getTaskName() + ": error in create a WatchService for separate search", e);
				}
			}
			// Close the watch service when done
			try {
				watchService.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}

	@Override
	protected String getVersion() {
		// TODO Auto-generated method stub
		return null;
	}
}
