package uk.ac.quadram.metapilot.task.dda.igc;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.concurrent.ExecutionException;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import uk.ac.quadram.metapilot.core.mod.PostTransModification;
import uk.ac.quadram.metapilot.core.mod.IsobaricTag;
import uk.ac.quadram.metapilot.dbSearch.dbReducer.DBReducerTask;
import uk.ac.quadram.metapilot.dbSearch.pfind.PFindIsobaricTask;
import uk.ac.quadram.metapilot.dbSearch.pfind.PFindLabelFreeTask;
import uk.ac.quadram.metapilot.dbSearch.pfind.PFindResultReader;
import uk.ac.quadram.metapilot.dbSearch.pfind.PFindTask;
import javafx.scene.control.ProgressBar;
import uk.ac.quadram.metapilot.quant.flashLFQ.FlashLfqTask;
import uk.ac.quadram.metapilot.quant.flashLFQ.FlashLfqQuanPepReader;
import uk.ac.quadram.metapilot.quant.flashLFQ.FlashLfqQuanPeptide;
import uk.ac.quadram.metapilot.quant.flashLFQ.FlashLfqPsmReader;
import uk.ac.quadram.metapilot.quant.flashLFQ.FlashLfqPsm;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters.ToolType;
import uk.ac.quadram.metapilot.task.igc.MetaPilotIgcTask;
import uk.ac.quadram.metapilot.task.dda.MetaDbReducerTask;
import uk.ac.quadram.metapilot.task.dda.PFindSpectraConvertTask;
import uk.ac.quadram.metapilot.task.dda.igc.par.DDAIgcPFindParameters;
import uk.ac.quadram.metapilot.task.io.par.MetaOpenPar;
import uk.ac.quadram.metapilot.task.par.MetaConstants;
import uk.ac.quadram.metapilot.task.par.MetaData;

public class DDAIgcPFindTask extends MetaPilotIgcTask {

	private static final Logger LOGGER = LogManager.getLogger(DDAIgcPFindTask.class);
	private static final String taskName = "DDA IGC database search by pFind workflow";

	protected HashMap<String, PostTransModification> ptmMap;

	// Fields to track running subtasks
	private PFindSpectraConvertTask spectraConvertTask;
	private MetaDbReducerTask metaDbReducerTask;
	private PFindTask pFindTask;
	private DBReducerTask dbReducerTask;
	private Thread watchServiceThread;

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
				if (rawFiles != null && expNames != null) {
					this.expNameMap = new HashMap<String, String>();
					this.quanExpNames = new String[rawFiles.length];
					for (int i = 0; i < rawFiles.length; i++) {
						LOGGER.info("{}: Raw file {}: {} -> {}", taskName, i, rawFiles[i],
								expNames != null && i < expNames.length ? expNames[i] : "N/A");

						this.quanExpNames[i] = rawFiles[i].substring(rawFiles[i].lastIndexOf("\\") + 1,
								rawFiles[i].lastIndexOf("."));
						this.expNameMap.put(quanExpNames[i], expNames[i]);
					}
					this.rawDir = (new File(rawFiles[0])).getParentFile();
				} else {
					LOGGER.warn("{}: Raw files are null in parameters", taskName);
				}
			} else {
				LOGGER.warn("{}: Metadata is null in parameters", taskName);
			}
		} else {
			LOGGER.error("{}: Parameters object is null!", taskName);
		}

		File pFindBin = (new File(mrp.getToolPath(ToolType.PFIND))).getParentFile();
		pFindTask = new PFindTask(pFindBin, metaPar.getThreadCount());
		ptmMap = pFindTask.getPtmMap();

		File dbReducerExe = new File(mrp.getToolPath(ToolType.DBREDUCER));
		dbReducerTask = new DBReducerTask(dbReducerExe);
	}

	@Override
	public void forceStop() {
		LOGGER.info(getTaskName() + ": Force stop requested");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": Force stop requested");

		// Cancel spectra conversion task if running
		if (spectraConvertTask != null && !spectraConvertTask.isDone()) {
			LOGGER.info(getTaskName() + ": Cancelling spectra conversion task");
			spectraConvertTask.cancel(true);
		}

		// Cancel database reduction task if running
		if (metaDbReducerTask != null && !metaDbReducerTask.isDone()) {
			LOGGER.info(getTaskName() + ": Cancelling database reduction task");
			metaDbReducerTask.cancel(true);
		}

		// Cancel pFind task if running
		if (pFindTask != null) {
			LOGGER.info(getTaskName() + ": Stopping pFind identification task");
			// PFindTask doesn't have forceStop method, but we can interrupt any running
			// threads
			// The task will be stopped when this main task is cancelled
		}

		// Cancel DBReducer task if running
		if (dbReducerTask != null) {
			LOGGER.info(getTaskName() + ": Stopping DBReducer task");
			// DBReducerTask doesn't have forceStop method, but we can interrupt any running
			// threads
			// The task will be stopped when this main task is cancelled
		}

		// Cancel watch service thread if running
		if (watchServiceThread != null && watchServiceThread.isAlive()) {
			LOGGER.info(getTaskName() + ": Interrupting watch service thread");
			watchServiceThread.interrupt();
		}

		// Cancel this task itself
		this.cancel(true);

		LOGGER.info(getTaskName() + ": Force stop completed");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": Force stop completed");
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
	protected String getVersion() {
		// TODO Auto-generated method stub
		return "1.0.0";
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

		spectraConvertTask = new PFindSpectraConvertTask(metaPar, mrp, secondBar);
		LOGGER.info("{}: PFindSpectraConvertTask created, executing...", taskName);
		spectraConvertTask.execute();
		try {
			boolean spectraConverted = spectraConvertTask.get(); // Wait for completion and get result
			if (!spectraConverted) {
				LOGGER.error(getTaskName() + ": error in pFind spectrum convert workflow");
				System.err.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": error in pFind spectrum convert workflow");
				return false;
			}
		} catch (InterruptedException | ExecutionException e) {
			LOGGER.error(getTaskName() + ": error in pFind spectrum convert workflow", e);
			System.err.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": error in pFind spectrum convert workflow: " + e.getMessage());
			return false;
		}

		// Update progress after spectra conversion
		updateProgress(20, 100);
		updateMessage("Spectra conversion completed, starting database reduction...");

		LOGGER.info(getTaskName() + ": pFind spectrum convert workflow: finished");
		System.out.println(
				format.format(new Date()) + "\t" + getTaskName() + ": pFind spectrum convert workflow: finished");

		// Update progress before database reduction
		updateProgress(25, 100);
		updateMessage("Creating sample-specific database...");

		metaDbReducerTask = new MetaDbReducerTask(metaPar, mrp, secondBar);
		LOGGER.info("{}: MetaDbReducerTask created, executing...", taskName);
		metaDbReducerTask.execute();
		try {
			boolean dbReduced = metaDbReducerTask.get(); // Wait for completion and get result
			if (!dbReduced) {
				LOGGER.error(getTaskName() + ": error in sample-specific database generation");
				System.err.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": error in sample-specific database generation");
				return false;
			}
		} catch (InterruptedException | ExecutionException e) {
			LOGGER.error(getTaskName() + ": error in sample-specific database generation", e);
			System.err.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": error in sample-specific database generation: " + e.getMessage());
			return false;
		}

		// Update progress after database reduction
		updateProgress(40, 100);
		updateMessage("Database reduction completed, starting identification...");

		LOGGER.info(getTaskName() + ": sample specific database generation: finished");
		System.out.println(
				format.format(new Date()) + "\t" + getTaskName() + ": sample specific database generation: finished");

		File pfindSearchFile = new File(metaPar.getResult(), "pfind_IGC_search");
		if (!pfindSearchFile.exists()) {
			pfindSearchFile.mkdir();
			LOGGER.info("  Created pFind search directory: " + pfindSearchFile.getAbsolutePath());
		} else {
			LOGGER.info("  pFind search directory already exists: " + pfindSearchFile.getAbsolutePath());
		}

		metaPar.setDbSearchResultFile(pfindSearchFile);
		
		this.combinedFile = new File(pfindSearchFile, "combined");
		if (!combinedFile.exists()) {
			combinedFile.mkdir();
			LOGGER.info("  Created combined search directory: " + combinedFile.getAbsolutePath());
		} else {
			LOGGER.info("  Combined search directory already exists: " + combinedFile.getAbsolutePath());
		}

		boolean iden = iden(spectraConvertTask.getMgfFiles(), spectraConvertTask.getPf2Files());
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

		HashMap<String, PostTransModification> ptmMap = pFindTask.getPtmMap();
		MetaOpenPar mosp = ((DDAIgcPFindParameters) metaPar).getMop();

		String ms2Type = ((DDAIgcPFindParameters) metaPar).getMs2ScanMode();
		boolean useDBReducer = false;

		// dbReducerTask2 is already declared as a field
		if (ms2Type.equals(MetaConstants.HCD_FTMS)
				&& ((DDAIgcPFindParameters) metaPar).getQuanMode().equals(MetaConstants.labelFree)) {
			useDBReducer = true;
		}

		PFindResultReader pFindResultReader = new PFindResultReader(combinedFile, metaPar.getSpectraFile(), ptmMap);

		if (pFindResultReader.hasResult()) {

			LOGGER.info(getTaskName() + ": peptide and protein identification results have been found in "
					+ combinedFile + ", go to next step");
			System.out.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": peptide and protein identification results have been found in " + combinedFile
					+ ", go to next step");

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

			pFindTask.searchPF2(pf2Files, combinedFile.getAbsolutePath(), db, ms2Type,
					((DDAIgcPFindParameters) metaPar).isOpenSearch(), true, mosp.getOpenPsmFDR(),
					((DDAIgcPFindParameters) this.metaPar).getFixMods(),
					((DDAIgcPFindParameters) this.metaPar).getVariMods(),
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

		PFindTask.cleanPFindResultFolder(combinedFile, true);

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

		PFindLabelFreeTask pfindLFTask = new PFindLabelFreeTask(flashLfqTask, metaPar.getMetadata(), combinedFile,
				rawDir, metaPar.getSpectraFile());
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

		PFindResultReader reader = new PFindResultReader(combinedFile, metaPar.getSpectraFile(), pFindTask.getPtmMap());
		PFindIsobaricTask pfindIsobaricTask = new PFindIsobaricTask(reader, this.metaPar.getSpectraFile(), combinedFile,
				((DDAIgcPFindParameters) this.metaPar).getIsobaricTag(),
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

			LOGGER.info(getTaskName() + ": exporting peptide report started");
			System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": exporting peptide report started");

		// Build file list and mapping
		ArrayList<String> fileList = new ArrayList<String>();
		HashMap<String, Integer> fileIdMap = new HashMap<String, Integer>();
		File expDesignFile = new File(rawDir, "ExperimentalDesign.tsv");

		if (expDesignFile.exists()) {
			try (BufferedReader reader = new BufferedReader(new FileReader(expDesignFile))) {
				String line = reader.readLine();
				while ((line = reader.readLine()) != null) {
					String[] cs = line.split("\t");
					fileIdMap.put(cs[0], fileIdMap.size());
					fileList.add(cs[0]);
				}
				reader.close();
			} catch (IOException e) {
				LOGGER.error(getTaskName() + ": error in reading the experimental design from " + expDesignFile, e);
				System.err.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": error in reading the experimental design from " + expDesignFile);
			}
		} else {
			File[] files = metaPar.getDbSearchResultFile().listFiles();
			for (int i = 0; i < files.length; i++) {
				if (files[i].isDirectory() && (new File(files[i], "pFind-Filtered.spectra")).exists()) {
					fileIdMap.put(files[i].getName(), fileIdMap.size());
					fileList.add(files[i].getName());
				}
			}
		}

		// Read PSMs from combined results (following DDAMagPfindTask pattern)
		FlashLfqPsm[] psms;
		try {
			FlashLfqPsmReader flashLfqPsmReader = new FlashLfqPsmReader(new File(combinedFile, "psms.tsv"));
			psms = flashLfqPsmReader.getPSMs();
			} catch (IOException e) {
			LOGGER.error(getTaskName() + ": error in reading the PSMs from " + (new File(combinedFile, "psms.tsv")), e);
			System.err.println(format.format(new Date()) + "\t" + getTaskName() + ": error in reading the PSMs from "
					+ (new File(combinedFile, "psms.tsv")));
			return false;
		}

		// Read quantified peptides (following DDAMagPfindTask pattern)
		FlashLfqQuanPepReader quanPepReader;
		String quanMode = ((DDAIgcPFindParameters) this.metaPar).getQuanMode();
		if (quanMode.equals(MetaConstants.labelFree)) {
			quanPepReader = new FlashLfqQuanPepReader(this.quan_pep_file, fileList);
		} else if (quanMode.equals(MetaConstants.isobaricLabel)) {
			IsobaricTag isoTag = ((DDAIgcPFindParameters) this.metaPar).getIsobaricTag();
			quanPepReader = new FlashLfqQuanPepReader(this.quan_pep_file, fileList, isoTag.getReportTitle());
		} else {
			LOGGER.error(getTaskName() + ": unknown quantitative type");
			System.err.println(format.format(new Date()) + "\t" + getTaskName() + ": unknown quantitative type");
			return false;
		}

		FlashLfqQuanPeptide[] quanPeps = quanPepReader.getMetaPeptides();

		// Build intensity title array
		String[] intensityTitle;
		if (quanMode.equals(MetaConstants.labelFree)) {
			intensityTitle = fileList.toArray(String[]::new);
		} else if (quanMode.equals(MetaConstants.isobaricLabel)) {
			IsobaricTag isoTag = ((DDAIgcPFindParameters) this.metaPar).getIsobaricTag();
			String[] repoTitle = isoTag.getReportTitle();
			intensityTitle = new String[fileList.size() * repoTitle.length];
			for (int i = 0; i < fileList.size(); i++) {
				for (int j = 0; j < repoTitle.length; j++) {
					intensityTitle[i * repoTitle.length + j] = fileList.get(i) + "_" + repoTitle[j];
				}
			}
						} else {
			LOGGER.error(getTaskName() + ": unknown quantitative type");
			return false;
		}

		// Create and execute the post-analysis task directly
		DDAIgcPostAnalysisTask postAnalysisTask = new DDAIgcPostAnalysisTask(
			metaPar, mrp, secondBar, quanPeps, psms, fileIdMap, fileList, quanMode, intensityTitle, getTaskName(), this);
		
		// Execute the task and wait for completion
		Thread taskThread = new Thread(postAnalysisTask);
		taskThread.setDaemon(true);
		taskThread.start();
		
		try {
			taskThread.join();
			boolean success = postAnalysisTask.get();
			
            if (success) {
                LOGGER.info(getTaskName() + ": Post-analysis completed successfully");
                return true;
					} else {
				LOGGER.error(getTaskName() + ": Post-analysis task failed");
			return false;
			}
				} catch (InterruptedException e) {
			LOGGER.error(getTaskName() + ": Post-analysis task was interrupted", e);
			Thread.currentThread().interrupt();
			return false;
		} catch (ExecutionException e) {
			LOGGER.error(getTaskName() + ": Post-analysis task failed with exception", e);
			return false;
		}
	}
}
