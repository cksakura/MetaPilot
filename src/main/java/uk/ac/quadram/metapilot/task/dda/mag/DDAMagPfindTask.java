package uk.ac.quadram.metapilot.task.dda.mag;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileFilter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.math.BigInteger;
import java.nio.file.FileSystems;
import java.nio.file.Path;
import java.nio.file.StandardWatchEventKinds;
import java.nio.file.WatchEvent;
import java.nio.file.WatchKey;
import java.nio.file.WatchService;
import java.sql.SQLException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.ExecutionException;

import org.apache.commons.lang3.math.NumberUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.dom4j.DocumentException;

import uk.ac.quadram.metapilot.core.mod.IsobaricTag;
import uk.ac.quadram.metapilot.core.taxonomy.TaxonomyRanks;
import uk.ac.quadram.metapilot.core.tools.FormatTool;
import uk.ac.quadram.metapilot.dbSearch.dbReducer.DBReducerTask;
import uk.ac.quadram.metapilot.dbSearch.open.io.OpenPsmTsvWriter;
import uk.ac.quadram.metapilot.dbSearch.pfind.PFindIsobaricTask;
import uk.ac.quadram.metapilot.dbSearch.pfind.PFindPSM;
import uk.ac.quadram.metapilot.dbSearch.pfind.PFindPsmReader;
import uk.ac.quadram.metapilot.dbSearch.pfind.PFindResultReader;
import uk.ac.quadram.metapilot.dbSearch.pfind.PFindTask;
import uk.ac.quadram.metapilot.gui.JavaFXProgressBar;
import uk.ac.quadram.metapilot.quant.flashLFQ.FlashLfqPsmReader;
import uk.ac.quadram.metapilot.quant.flashLFQ.FlashLfqQuanPepReader;
import uk.ac.quadram.metapilot.quant.flashLFQ.FlashLfqQuanPeptide;
import uk.ac.quadram.metapilot.quant.flashLFQ.FlashLfqTask;
import uk.ac.quadram.metapilot.quant.flashLFQ.FlashLfqPsmReader.FlashLfqPsm;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters.ToolType;
import uk.ac.quadram.metapilot.task.MetaPythonReportTask;
import uk.ac.quadram.metapilot.task.MetaReportCopyTask;
import uk.ac.quadram.metapilot.task.dda.PFindSpectraConvertTask;
import uk.ac.quadram.metapilot.task.io.MetaTreeHandler;
import uk.ac.quadram.metapilot.task.io.pro.MetaProtein;
import uk.ac.quadram.metapilot.task.io.pro.pfind.MetaProteinXMLReader2;
import uk.ac.quadram.metapilot.task.mag.HMGenomeHtmlWriter;
import uk.ac.quadram.metapilot.task.mag.MagFuncSearcher;
import uk.ac.quadram.metapilot.task.mag.MetaPilotMagTask;
import uk.ac.quadram.metapilot.task.mag.MetaProteinAnnoMag;
import uk.ac.quadram.metapilot.task.mag.par.MetaParameterMag;
import uk.ac.quadram.metapilot.task.par.MetaBaseParameter;
import uk.ac.quadram.metapilot.task.par.MetaConstants;
import uk.ac.quadram.metapilot.task.pfind.par.MetaParameterPFind;
import uk.ac.quadram.metapilot.task.v2.par.MetaOpenPar;

public class DDAMagPfindTask extends MetaPilotMagTask {

	private String[] mgfFiles;
	private String[] pf2Files;
	private File dbReducerExe;
	private File pFindExe;
	private boolean useDBReducer = false;
	
	
	
	private static final Logger LOGGER = LogManager.getLogger(DDAMagPfindTask.class);
	private static final String taskName = "DDA MAG database search by pFind workflow";

	public DDAMagPfindTask(MetaBaseParameter metaPar, MetaResourceParameters mrp, JavaFXProgressBar progressBar,
			JavaFXProgressBar secondBar) {
		super(metaPar, mrp, progressBar, secondBar);
		// TODO Auto-generated constructor stub

		
		// TODO Auto-generated method stub
		String[] rawFiles = metaPar.getMetadata().getRawFiles();
		String ms2Type = this.metaPar.getMs2ScanMode();
		String suffix = "";
		if (ms2Type.equals(MetaConstants.HCD_FTMS)) {
			suffix = "_HCDFT";
		} else if (ms2Type.equals(MetaConstants.HCD_ITMS)) {
			suffix = "_HCDIT";
		} else if (ms2Type.equals(MetaConstants.CID_FTMS)) {
			suffix = "_CIDFT";
		} else if (ms2Type.equals(MetaConstants.CID_ITMS)) {
			suffix = "_CIDIT";
		}
		this.mgfFiles = new String[rawFiles.length];
		this.pf2Files = new String[rawFiles.length];
		this.expNameMap = new HashMap<String, String>();
		File spectraDirFile = metaPar.getSpectraFile();
		for (int i = 0; i < rawFiles.length; i++) {
			File rawi = new File(rawFiles[i]);
			this.rawDir = rawi.getParentFile();
			String rawName = rawi.getName();
			rawName = rawName.substring(0, rawName.length() - ".raw".length());
			File mgfFile = new File(spectraDirFile, rawName + suffix + ".mgf");
			File pf2 = new File(spectraDirFile, rawName + suffix + ".pf2");
			mgfFiles[i] = mgfFile.getAbsolutePath();
			pf2Files[i] = pf2.getAbsolutePath();

			this.separateResultFolders = new File[rawFiles.length];
			this.separateResultFolders[i] = new File(metaPar.getDbSearchResultFile(), rawName);
			if (!this.separateResultFolders[i].exists()) {
				this.separateResultFolders[i].mkdir();
			}
		}

		this.dbReducerExe = new File(mrp.getToolPath(ToolType.DBREDUCER));
		this.pFindExe = new File(mrp.getToolPath(ToolType.PFIND));
		if (dbReducerExe.getName().equals("DBReducer.exe") && dbReducerExe.exists()
				&& metaPar.getMs2ScanMode().equals(MetaConstants.HCD_FTMS)) {
			useDBReducer = true;
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
		// TODO Auto-generated method stub
		progressBar.setString(getTaskName());
		setProgress(0);

		LOGGER.info(getTaskName() + ": start");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": start");

		PFindSpectraConvertTask spTask = new PFindSpectraConvertTask(metaPar, mrp, secondBar);
		spTask.execute();
		boolean convert = spTask.getValue();
		if (!convert) {
			LOGGER.error(getTaskName() + ": error in spectra format conversion");
			System.err
					.println(format.format(new Date()) + "\t" + getTaskName() + ": error in spectra format conversion");
			return false;
		}

		LOGGER.info(getTaskName() + ": searching high abundant database started");
		System.out.println(
				format.format(new Date()) + "\t" + getTaskName() + ": searching high abundant database started");

		boolean searchHap = this.separateSearchHap();

		if (searchHap) {
			LOGGER.info(getTaskName() + ": searching high abundant database finished");
			System.out.println(
					format.format(new Date()) + "\t" + getTaskName() + ": searching high abundant database finished");
		} else {
			LOGGER.error(getTaskName() + ": error in searching high abundant database");
			System.err.println(
					format.format(new Date()) + "\t" + getTaskName() + ": error in searching high abundant database");

			return false;
		}

		LOGGER.info(getTaskName() + ": extracting proteins from the HAP result started");
		System.out.println(
				format.format(new Date()) + "\t" + getTaskName() + ": extracting proteins from the HAP result started");

		boolean extract = this.fastaExtract4Iterator();

		if (extract) {
			LOGGER.info(getTaskName() + ": extracting proteins from the HAP result finished");
			System.out.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": extracting proteins from the HAP result finished");
		} else {
			LOGGER.error(getTaskName() + ": error in extracting proteins from the HAP result");
			System.err.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": error in extracting proteins from the HAP result");

			return false;
		}

		setProgress(40);
		
		progressBar.setString(getTaskName() + ": peptide and protein identification");

		LOGGER.info(getTaskName() + ": peptide and protein identification started");
		System.out.println(
				format.format(new Date()) + "\t" + getTaskName() + ": peptide and protein identification started");

		boolean search = separateSearch();

		if (search) {
			LOGGER.info(getTaskName() + ": database search through the sample-specific database finished");
			System.out.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": database search through the sample-specific database finished");
		} else {
			LOGGER.error(getTaskName() + ": error in database search");
			System.err.println(format.format(new Date()) + "\t" + getTaskName() + ": error in database search");

			return false;
		}

		boolean combine = combineResult();

		if (combine) {
			LOGGER.info(getTaskName() + ": combining of the separate search results finished");
			System.out.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": combining of the separate search results finished");
		} else {
			LOGGER.error(getTaskName() + ": error in the combination of the separate search results");
			System.err.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": error in the combination of the separate search results");
			return false;
		}

		LOGGER.info(getTaskName() + ": peptide and protein identification finished");
		System.out.println(
				format.format(new Date()) + "\t" + getTaskName() + ": peptide and protein identification finished");

		setProgress(70);
		
		boolean quant = false;
		if (((MetaParameterPFind) metaPar).getQuanMode().equals(MetaConstants.labelFree)) {
			quant = lfQuant();
			if (!quant) {
				LOGGER.error(getTaskName() + ": error in peptide and protein quantification");
				System.err.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": error in peptide and protein quantification");
				return false;
			}
		} else if (((MetaParameterPFind) metaPar).getQuanMode().equals(MetaConstants.isobaricLabel)) {
			quant = lfQuant();
			if (!quant) {
				LOGGER.error(getTaskName() + ": error in peptide and protein quantification");
				System.err.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": error in peptide and protein quantification");
				return false;
			}

			setProgress(75);

			quant = isobaricQuant();
			if (!quant) {
				LOGGER.error(getTaskName() + ": error in peptide and protein quantification");
				System.err.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": error in peptide and protein quantification");
				return false;
			}
		}

		setProgress(80);

		exportReport();

		setProgress(100);
		
		return true;
	}

	protected boolean separateSearchHap() {

		int taskCount = 0;
		this.separateResultFolders = new File[rawFiles.length];
		for (int i = 0; i < rawFiles.length; i++) {
			String name = rawFiles[i].substring(rawFiles[i].lastIndexOf("\\") + 1, rawFiles[i].lastIndexOf("."));

			this.separateResultFolders[i] = new File(metaPar.getDbSearchResultFile(), name);
			if (!this.separateResultFolders[i].exists()) {
				this.separateResultFolders[i].mkdir();
			}

			File hapFile = new File(separateResultFolders[i], "hap");
			if (!hapFile.exists()) {
				hapFile.mkdirs();
			}

			PFindResultReader reader = new PFindResultReader(hapFile, null);
			if (reader.hasResult()) {
				System.out.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": pFind search result already exists in " + separateResultFolders[i].getAbsolutePath());

				LOGGER.info(getTaskName() + ": pFind search result already exists in "
						+ separateResultFolders[i].getAbsolutePath());
			} else {
				taskCount++;
			}
		}
		MetaOpenPar mosp = ((MetaParameterMag) metaPar).getMop();
		String ms2Type = ((MetaParameterMag) metaPar).getMs2ScanMode();
		if (taskCount > 0) {
			int totalThreadCount = metaPar.getThreadCount();
			int threadCount = totalThreadCount > taskCount ? totalThreadCount / taskCount : 1;
			int threadPoolCount = totalThreadCount > taskCount ? taskCount : totalThreadCount;

			if (this.useDBReducer) {
				DBReducerTask task = new DBReducerTask(dbReducerExe);
				task.setTotalThread(threadCount);
				for (int i = 0; i < this.mgfFiles.length; i++) {
					File mgfFile = new File(mgfFiles[i]);
					File hapFile = new File(separateResultFolders[i], "hap");
					if (mgfFile.exists()) {
						task.searchMgf(new String[] { mgfFiles[i] }, hapFile.getAbsolutePath(),
								this.magDb.getHapFastaFile().getAbsolutePath(), ms2Type, mosp.getOpenPsmFDR(),
								((MetaParameterMag) this.metaPar).getFixMods(),
								((MetaParameterMag) this.metaPar).getVariMods());
					} else {

						System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": spectra file "
								+ mgfFiles[i] + " was not found");

						LOGGER.error(getTaskName() + ": spectra file " + mgfFile + " was not found");

						return false;
					}
				}
				task.run(threadPoolCount, separateResultFolders.length * 3);

			} else {
				PFindTask task = new PFindTask(pFindExe.getParentFile());
				task.setTotalThread(threadCount);

				for (int i = 0; i < this.pf2Files.length; i++) {
					File hapFile = new File(separateResultFolders[i], "hap");
					File pf2 = new File(pf2Files[i]);
					if (pf2.exists()) {
						task.searchPF2(new String[] { pf2.getAbsolutePath() }, hapFile.getAbsolutePath(),
								this.magDb.getHapFastaFile().getAbsolutePath(), ms2Type,
								((MetaParameterMag) metaPar).isOpenSearch(), true, mosp.getOpenPsmFDR(),
								((MetaParameterMag) this.metaPar).getFixMods(),
								((MetaParameterMag) this.metaPar).getVariMods(),
								((MetaParameterMag) this.metaPar).getEnzyme(),
								((MetaParameterMag) this.metaPar).getMissCleavages(),
								((MetaParameterMag) this.metaPar).getDigestMode(),
								((MetaParameterMag) this.metaPar).getIsobaricTag());
					} else {

					}
				}
				task.run(threadPoolCount, separateResultFolders.length * 3);
			}
		}

		FileFilter filter = new FileFilter() {
			@Override
			public boolean accept(File file) {
				return file.getName().equals("pFind-Filtered.spectra");
			}
		};

		// The WatchService for monitoring the folders
		WatchService watchService = null;
		// Create the WatchService that monitors the folders
		try {
			watchService = FileSystems.getDefault().newWatchService();
			// Register each folder with the watch service
			for (File folder : separateResultFolders) {
				File hapFile = new File(folder, "hap");
				Path path = hapFile.toPath();
				path.register(watchService, StandardWatchEventKinds.ENTRY_CREATE);
			}
		} catch (IOException e) {
			System.out.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": error in create a WatchService for separate search");

			LOGGER.error(getTaskName() + ": error in create a WatchService for separate search", e);
		}

		new Thread(new WatchServiceLoop(watchService, filter, taskCount, separateResultFolders.length, 10, 38)).start();

		for (int i = 0; i < this.separateResultFolders.length; i++) {
			File hapFile = new File(separateResultFolders[i], "hap");
			PFindTask.cleanPFindResultFolder(hapFile, true);

			PFindResultReader reader = new PFindResultReader(hapFile, null);
			if (!reader.hasResult()) {
				System.out.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": pFind search result was not found from " + hapFile.getAbsolutePath());

				LOGGER.error(getTaskName() + ": pFind search result was not found from " + hapFile.getAbsolutePath());

				return false;
			}
		}

		return true;
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

							DDAMagPfindTask.this.setProgress(barStart + percentage);

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

	protected boolean fastaExtract4Iterator() {

		for (int i = 0; i < this.separateResultFolders.length; i++) {

			File hapFile = new File(separateResultFolders[i], "hap");
			File dbFile = new File(hapFile, separateResultFolders[i].getName() + ".fasta");

			if (dbFile.exists() && dbFile.length() > 0) {

				System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": sample-specific database for "
						+ separateResultFolders[i].getName() + " already exists in " + dbFile.getAbsolutePath());

				LOGGER.info(getTaskName() + ": sample-specific database for " + separateResultFolders[i].getName()
						+ " already exists in " + dbFile.getAbsolutePath());

				continue;
			}

			File spectraFile = new File(hapFile, "pFind-Filtered.spectra");

			HashMap<String, HashSet<String>> genomeMap = new HashMap<String, HashSet<String>>();
			HashSet<String> pepSet = new HashSet<String>();

			try (BufferedReader hapReader = new BufferedReader(new FileReader(spectraFile))) {
				String line = hapReader.readLine();
				while ((line = hapReader.readLine()) != null) {
					String[] cs = line.split("\t");
					String[] pros = cs[12].split("/");
					for (int j = 0; j < pros.length; j++) {
						if (pros[j].startsWith(magDb.getIdentifier())) {
							String genome = pros[j].split("_")[0];
							if (genomeMap.containsKey(genome)) {
								genomeMap.get(genome).add(cs[0]);
							} else {
								HashSet<String> set = new HashSet<String>();
								set.add(cs[0]);
								genomeMap.put(genome, set);
							}
						}
					}
					pepSet.add(cs[0]);
				}
				hapReader.close();
			} catch (IOException e) {
				System.out.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": error in reading the pFind search result file from " + spectraFile.getAbsolutePath());

				LOGGER.error(getTaskName() + ": error in reading the pFind search result file from "
						+ spectraFile.getAbsolutePath());

				return false;
			}

			System.out.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": start generating sample-specific database for " + separateResultFolders[i].getName());

			LOGGER.info(getTaskName() + ": start generating sample-specific database for "
					+ separateResultFolders[i].getName());

			File hapGenomeFile = new File(hapFile, "hap_genomes.tsv");
			HashSet<String> genomeSet = refineGenomes(genomeMap, pepSet.size(), psmCoverThreshold, hapGenomeFile);

			try (PrintWriter writer = new PrintWriter(dbFile)) {
				File dbFolder = this.magDb.getOriginalDbFolder();
				for (String genome : genomeSet) {
					File genomeDbFile = new File(dbFolder, genome + ".faa");
					if (!genomeDbFile.exists()) {
						genomeDbFile = new File(dbFolder, genome + ".fasta");
					}
					if (genomeDbFile.exists()) {
						BufferedReader readeri = new BufferedReader(new FileReader(genomeDbFile));
						String linei = null;
						while ((linei = readeri.readLine()) != null) {
							writer.println(linei);
						}
						readeri.close();
					}
				}

				if (metaPar.isAppendHostDb()) {
					BufferedReader hostReader = new BufferedReader(new FileReader(this.metaPar.getHostDB()));
					String linei = null;
					while ((linei = hostReader.readLine()) != null) {
						writer.println(linei);
					}
					hostReader.close();
				}

				writer.close();
			} catch (IOException e) {
				System.out.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": error in writing the sample-specific database to " + dbFile.getAbsolutePath());

				LOGGER.error(getTaskName() + ": error in writing the sample-specific database to "
						+ dbFile.getAbsolutePath());

				return false;
			}

			System.out.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": generating sample-specific database for " + separateResultFolders[i].getName() + " finished");

			LOGGER.info(getTaskName() + ": generating sample-specific database for "
					+ separateResultFolders[i].getName() + " finished");
		}

		return true;
	}

	protected boolean separateSearch() {

		int taskCount = 0;
		for (int i = 0; i < this.separateResultFolders.length; i++) {

			PFindResultReader reader = new PFindResultReader(separateResultFolders[i], null);
			if (reader.hasResult()) {
				System.out.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": pFind search result already exists in " + separateResultFolders[i].getAbsolutePath());

				LOGGER.info(getTaskName() + ": pFind search result already exists in "
						+ separateResultFolders[i].getAbsolutePath());
			} else {
				File hapFile = new File(separateResultFolders[i], "hap");
				File db = new File(hapFile, separateResultFolders[i].getName() + ".fasta");
				if (!db.exists()) {
					System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": sample-specific database"
							+ db + " was not found");

					LOGGER.error(getTaskName() + ": sample-specific database" + db + " was not found");

					return false;
				}
				taskCount++;
			}
		}
		MetaOpenPar mosp = ((MetaParameterMag) metaPar).getMop();
		String ms2Type = ((MetaParameterMag) metaPar).getMs2ScanMode();
		if (taskCount > 0) {
			int totalThreadCount = metaPar.getThreadCount();
			int threadCount = totalThreadCount > taskCount ? totalThreadCount / taskCount : 1;
			int threadPoolCount = totalThreadCount > taskCount ? taskCount : totalThreadCount;

			if (this.useDBReducer) {
				DBReducerTask task = new DBReducerTask(dbReducerExe);
				task.setTotalThread(threadCount);

				for (int i = 0; i < this.mgfFiles.length; i++) {
					File dbReducerFile = new File(separateResultFolders[i], "DBReducer");
					if (!dbReducerFile.exists()) {
						dbReducerFile.mkdir();
					}
					File reducedDb = new File(dbReducerFile, "DBReducer.fasta");
					if (!reducedDb.exists()) {
						File mgf = new File(mgfFiles[i]);
						if (mgf.exists()) {
							File hapFile = new File(separateResultFolders[i], "hap");
							File db = new File(hapFile, separateResultFolders[i].getName() + ".fasta");
							if (db.length() > 0) {
								task.searchMgf(new String[] { mgf.getAbsolutePath() }, dbReducerFile.getAbsolutePath(),
										db.getAbsolutePath(), ms2Type, mosp.getOpenPsmFDR(),
										((MetaParameterMag) this.metaPar).getFixMods(),
										((MetaParameterMag) this.metaPar).getVariMods());
							} else {
								System.out.println(format.format(new Date()) + "\t" + getTaskName()
										+ ": no protein sequences were in " + db
										+ ", which means no peptides/proteins can be identifed from "
										+ separateResultFolders[i].getName());

								LOGGER.error(getTaskName() + ": no protein sequences were in " + db
										+ ", which means no peptides/proteins can be identifed from "
										+ separateResultFolders[i].getName());
							}
						} else {

							System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": spectra file "
									+ mgf + " was not found");

							LOGGER.error(getTaskName() + ": spectra file " + mgf + " was not found");

							return false;
						}
					}
				}
				task.run(threadPoolCount, separateResultFolders.length * 3);

				PFindTask pFindTask = new PFindTask(pFindExe.getParentFile());
				pFindTask.setTotalThread(threadCount);

				for (int i = 0; i < this.pf2Files.length; i++) {
					File pf2 = new File(pf2Files[i]);
					File dbReducerFile = new File(separateResultFolders[i], "DBReducer");
					File db = new File(dbReducerFile, "DBReducer.fasta");

					if (pf2.exists()) {

						if (db.exists() && db.length() > 0) {
							pFindTask.searchPF2(new String[] { pf2.getAbsolutePath() },
									separateResultFolders[i].getAbsolutePath(), db.getAbsolutePath(), ms2Type,
									((MetaParameterMag) this.metaPar).isOpenSearch(), true, mosp.getOpenPsmFDR(),
									((MetaParameterMag) this.metaPar).getFixMods(),
									((MetaParameterMag) this.metaPar).getVariMods(),
									((MetaParameterMag) this.metaPar).getEnzyme(),
									((MetaParameterMag) this.metaPar).getMissCleavages(),
									((MetaParameterMag) this.metaPar).getDigestMode(),
									((MetaParameterMag) this.metaPar).getIsobaricTag());

						} else {

							System.out.println(
									format.format(new Date()) + "\t" + getTaskName() + ": no protein sequences were in "
											+ db + ", which means no peptides/proteins can be identifed from "
											+ separateResultFolders[i].getName());

							LOGGER.error(getTaskName() + ": no protein sequences were in " + db
									+ ", which means no peptides/proteins can be identifed from "
									+ separateResultFolders[i].getName());
						}
					}
				}

				pFindTask.run(threadPoolCount, this.separateResultFolders.length * 10);

			} else {
				PFindTask pFindTask = new PFindTask(pFindExe.getParentFile());
				pFindTask.setTotalThread(threadCount);

				for (int i = 0; i < this.pf2Files.length; i++) {
					File pf2 = new File(pf2Files[i]);
					File hapFile = new File(separateResultFolders[i], "hap");
					File db = new File(hapFile, separateResultFolders[i].getName() + ".fasta");
					if (pf2.exists()) {
						if (db.exists() && db.length() > 0) {
							pFindTask.searchPF2(new String[] { pf2.getAbsolutePath() },
									separateResultFolders[i].getAbsolutePath(), db.getAbsolutePath(), ms2Type,
									((MetaParameterMag) this.metaPar).isOpenSearch(), true, mosp.getOpenPsmFDR(),
									((MetaParameterMag) this.metaPar).getFixMods(),
									((MetaParameterMag) this.metaPar).getVariMods(),
									((MetaParameterMag) this.metaPar).getEnzyme(),
									((MetaParameterMag) this.metaPar).getMissCleavages(),
									((MetaParameterMag) this.metaPar).getDigestMode(),
									((MetaParameterMag) this.metaPar).getIsobaricTag());
						} else {

							System.out.println(
									format.format(new Date()) + "\t" + getTaskName() + ": no protein sequences were in "
											+ db + ", which means no peptides/proteins can be identifed from "
											+ separateResultFolders[i].getName());

							LOGGER.error(getTaskName() + ": no protein sequences were in " + db
									+ ", which means no peptides/proteins can be identifed from "
									+ separateResultFolders[i].getName());
						}
					}
				}

				pFindTask.run(threadPoolCount, this.separateResultFolders.length * 10);
			}
		}

		FileFilter filter = new FileFilter() {
			@Override
			public boolean accept(File file) {
				return file.getName().equals("pFind-Filtered.spectra");
			}
		};

		// The WatchService for monitoring the folders
		WatchService watchService = null;
		// Create the WatchService that monitors the folders
		try {
			watchService = FileSystems.getDefault().newWatchService();
			// Register each folder with the watch service
			for (File folder : separateResultFolders) {
				Path path = folder.toPath();
				path.register(watchService, StandardWatchEventKinds.ENTRY_CREATE);
			}
		} catch (IOException e) {
			e.printStackTrace();
		}

		new Thread(new WatchServiceLoop(watchService, filter, taskCount, separateResultFolders.length, 40, 68)).start();

		for (int i = 0; i < this.separateResultFolders.length; i++) {
			PFindTask.cleanPFindResultFolder(separateResultFolders[i], true);

			File hapFile = new File(separateResultFolders[i], "hap");
			if (hapFile.exists()) {
				PFindTask.cleanPFindResultFolder(hapFile, true);
			}

			File dbReduceFile = new File(separateResultFolders[i], "DBReducer");
			if (dbReduceFile.exists()) {
				PFindTask.cleanPFindResultFolder(dbReduceFile, true);
			}
		}

		return true;
	}

	protected boolean lfQuant() {
		
		progressBar.setString(getTaskName() + ": peptide and protein quantification");

		LOGGER.info(getTaskName() + ": peptide and protein label-free quantification started");
		System.out.println(format.format(new Date()) + "\t" + getTaskName()
				+ ": peptide and protein label-free quantification started");

		String flashLFQ = mrp.getToolPath(ToolType.FLASHLFQ);
		FlashLfqTask flashLfqTask = new FlashLfqTask(flashLFQ, metaPar.getThreadCount());
		File rawDir = (new File(rawFiles[0])).getParentFile();
		this.quan_pep_file = new File(metaPar.getDbSearchResultFile(), "QuantifiedPeptides.tsv");
		if (!quan_pep_file.exists() || quan_pep_file.length() == 0) {

			File msms = new File(metaPar.getDbSearchResultFile(), "psms.tsv");
			boolean lfqQuan = flashLfqTask.run(metaPar.getMetadata(), msms.getAbsolutePath(), rawDir.getAbsolutePath(), true, true);
			if (!lfqQuan) {

				LOGGER.info(getTaskName() + ": error in peptide and protein label-free quantification");
				System.out.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": error in peptide and protein label-free quantification");

				return false;

			} else {
				if (!quan_pep_file.exists()) {

					LOGGER.info(getTaskName() + ": error in peptide and protein label-free quantification");
					System.out.println(format.format(new Date()) + "\t" + getTaskName()
							+ ": error in peptide and protein label-free quantification");

					return false;

				}
			}
		} else {
			LOGGER.info(getTaskName() + ": peptide quantification result has been found in " + quan_pep_file
					+ ", go to next step");
			System.out.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": peptide quantification result has been found in " + quan_pep_file + ", go to next step");
		}

		LOGGER.info(getTaskName() + ": peptide and protein quantification label-free finished");
		System.out.println(format.format(new Date()) + "\t" + getTaskName()
				+ ": peptide and protein label-free quantification finished");

		return true;
	}

	protected boolean isobaricQuant() {

		progressBar.setString(getTaskName() + ": peptide and protein quantification");

		LOGGER.info(getTaskName() + ": peptide and protein isobaric quantification started");
		System.out.println(format.format(new Date()) + "\t" + getTaskName()
				+ ": peptide and protein isobaric quantification started");

		ArrayList<PFindPSM> psmList = new ArrayList<PFindPSM>();
		for (int i = 0; i < this.separateResultFolders.length; i++) {
			PFindPsmReader psmReader = new PFindPsmReader(new File(this.separateResultFolders[i], "pFind.spectra"));
			PFindPSM[] psms = psmReader.getPsms();
			for (int j = 0; j < psms.length; j++) {
				psmList.add(psms[j]);
			}
		}
		PFindPSM[] psms = psmList.toArray(new PFindPSM[psmList.size()]);
		PFindIsobaricTask pfindIsobaricTask = new PFindIsobaricTask(psms, this.metaPar.getSpectraFile(),
				metaPar.getDbSearchResultFile(), ((MetaParameterPFind) this.metaPar).getIsobaricTag(),
				((MetaParameterPFind) this.metaPar).getIsobaric(),
				((MetaParameterPFind) this.metaPar).getIsoCorFactor(),
				((MetaParameterPFind) this.metaPar).isCombineLabel(), metaPar.getMetadata());

		try {
			pfindIsobaricTask.runOnlyPeptide();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			LOGGER.error(getTaskName() + ": error in quantifying PSMs", e);
			System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": error in quantifying PSMs");
		}

		 File quan_pep_file = pfindIsobaricTask.getQuanPepFile();

		if (quan_pep_file.exists()) {
			LOGGER.info(getTaskName() + ": quantification results have been fount in " + metaPar.getDbSearchResultFile());
			System.out.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": quantification results have been fount in " + metaPar.getDbSearchResultFile());
			return true;
		}

		if (!quan_pep_file.exists()) {
			return false;
		}

		LOGGER.info(getTaskName() + ": peptide and protein isobaric quantification finished");
		System.out.println(format.format(new Date()) + "\t" + getTaskName()
				+ ": peptide and protein isobaric quantification finished");

		return true;
	}

	protected boolean combineResult() {
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": result combination started");
		LOGGER.info(getTaskName() + ": result combination started");

		File msms = new File(metaPar.getDbSearchResultFile(), "psms.tsv");

		if (msms.exists() && msms.length() > 200) {

			LOGGER.info(getTaskName() + ": peptide and protein identification results have been found in "
					+ metaPar.getDbSearchResultFile() + ", go to next step");
			System.out.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": peptide and protein identification results have been found in "
					+ metaPar.getDbSearchResultFile() + ", go to next step");

			return true;
		}

		File rtFile = new File(metaPar.getSpectraFile(), "rt.txt");
		HashMap<String, Double> rtMap = new HashMap<String, Double>();
		try (BufferedReader reader = new BufferedReader(new FileReader(rtFile))) {
			String line = reader.readLine();
			while ((line = reader.readLine()) != null) {
				String[] cs = line.split("\t");
				rtMap.put(cs[0] + "_" + cs[1], Double.parseDouble(cs[2]));
			}
			reader.close();
		} catch (IOException e) {
			System.out.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": error in reading the retention time information from " + rtFile);
			LOGGER.error(getTaskName() + ": error in reading the retention time information from " + rtFile);
		}

		ArrayList<PFindPSM[]> psmList = new ArrayList<PFindPSM[]>();
		HashMap<String, HashSet<String>> pepProMap = new HashMap<String, HashSet<String>>();

		for (int i = 0; i < this.separateResultFolders.length; i++) {
			PFindPsmReader psmReader = new PFindPsmReader(new File(this.separateResultFolders[i], "pFind.spectra"));

			PFindPSM[] psms = psmReader.getPsms();

			for (int j = 0; j < psms.length; j++) {
				String sequence = psms[j].getSequence();

				if (sequence.contains("X")) {
					continue;
				}

				String[] proteins = psms[j].getProteins();
				HashSet<String> proSet;
				if (pepProMap.containsKey(sequence)) {
					proSet = pepProMap.get(sequence);
				} else {
					proSet = new HashSet<String>();
					pepProMap.put(sequence, proSet);
				}
				for (int k = 0; k < proteins.length; k++) {
					proSet.add(proteins[k]);
				}
			}

			psmList.add(psms);
		}

		OpenPsmTsvWriter msmsWriter = new OpenPsmTsvWriter(msms, OpenPsmTsvWriter.pFindFormat);
		for (int i = 0; i < psmList.size(); i++) {
			PFindPSM[] psms = psmList.get(i);
			for (int j = 0; j < psms.length; j++) {
				String sequence = psms[j].getSequence();

				if (pepProMap.containsKey(sequence)) {
					HashSet<String> proSet = pepProMap.get(sequence);
					String[] proteins = proSet.toArray(new String[proSet.size()]);
					psms[j].setProteins(proteins);

					String scanKey = psms[j].getFileName() + "_" + psms[j].getScan();

					if (rtMap.containsKey(scanKey)) {
						psms[j].setRt(rtMap.get(scanKey));
						msmsWriter.addPFindPSM(psms[j]);
					}
				}
			}
		}

		msmsWriter.close();

		if (msms.exists() && msms.length() > 200) {

			System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": result combination finished");
			LOGGER.info(getTaskName() + ": result combination finished");

			return true;
		} else {

			System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": result combination failed");
			LOGGER.info(getTaskName() + ": result combination failed");

			return false;
		}
	}
	
	protected boolean exportReport() {

		LOGGER.info(getTaskName() + ": exporting peptide report started");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": exporting peptide report started");

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

		HashMap<String, FlashLfqPsm> flashLfqPsmMap = new HashMap<String, FlashLfqPsm>();
		HashMap<String, int[]> pepPsmCountMap = new HashMap<String, int[]>();
		HashMap<String, int[]> pepChargeMap = new HashMap<String, int[]>();
		HashMap<String, double[]> pepExpectMap = new HashMap<String, double[]>();
		HashMap<String, Integer> psmCountMap = new HashMap<String, Integer>();

		int[] filePsmCount = new int[fileIdMap.size()];
		HashSet<String>[] fileSequenceSets = new HashSet[fileIdMap.size()];
		for (int i = 0; i < fileSequenceSets.length; i++) {
			fileSequenceSets[i] = new HashSet<String>();
		}

		try {
			FlashLfqPsmReader flashLfqPsmReader = new FlashLfqPsmReader(new File(metaPar.getDbSearchResultFile(), "psms.tsv"));
			FlashLfqPsm[] psms = flashLfqPsmReader.getPSMs();
			for (FlashLfqPsm psm : psms) {
				String sequence = psm.getBaseSequence();
				String fileName = psm.getFileName();
				String modseq = psm.getFullSequence();
				int charge = psm.getCharge();
				double PEP = psm.getPEP();
				int fileId = fileIdMap.get(fileName);

				int[] psmCount;
				if (pepPsmCountMap.containsKey(modseq)) {
					psmCount = pepPsmCountMap.get(modseq);
				} else {
					psmCount = new int[fileIdMap.size()];
					pepPsmCountMap.put(modseq, psmCount);
				}
				if (fileIdMap.containsKey(fileName)) {
					psmCount[fileId]++;
				}

				filePsmCount[fileId]++;
				fileSequenceSets[fileId].add(sequence);

				int[] charges;
				if (pepChargeMap.containsKey(modseq)) {
					charges = pepChargeMap.get(modseq);
				} else {
					charges = new int[7];
					pepChargeMap.put(modseq, charges);
				}
				if (charge < charges.length) {
					charges[charge] = 1;
				}
				double[] expect;
				if (pepExpectMap.containsKey(modseq)) {
					expect = pepExpectMap.get(modseq);
				} else {
					expect = new double[fileIdMap.size()];
					Arrays.fill(expect, 1.0);
					pepExpectMap.put(modseq, expect);
				}
				if (PEP < expect[fileId]) {
					expect[fileId] = PEP;
				}

				if (flashLfqPsmMap.containsKey(modseq)) {
					if (PEP < flashLfqPsmMap.get(modseq).getPEP()) {
						flashLfqPsmMap.put(modseq, psm);
					}
				} else {
					flashLfqPsmMap.put(modseq, psm);
				}

				String[] pros = psm.getProteins();
				for (String pro : pros) {
					if (psmCountMap.containsKey(pro)) {
						psmCountMap.put(pro, psmCountMap.get(pro) + 1);
					} else {
						psmCountMap.put(pro, 1);
					}

					if (pro.startsWith(this.magDb.getIdentifier())) {
						String genome = pro.split("_")[0];
						if (psmCountMap.containsKey(genome)) {
							psmCountMap.put(genome, psmCountMap.get(genome) + 1);
						} else {
							psmCountMap.put(genome, 1);
						}
					}
				}
			}

		} catch (IOException e) {
			// TODO Auto-generated catch block
			LOGGER.error(
					getTaskName() + ": error in reading the PSMs from  " + (new File(metaPar.getDbSearchResultFile(), "psms.tsv")), e);
			System.err.println(format.format(new Date()) + "\t" + getTaskName() + ": error in reading the PSMs from  "
					+ (new File(metaPar.getDbSearchResultFile(), "psms.tsv")));
		}

		FlashLfqQuanPepReader quanPepReader;

		String quanMode = ((MetaParameterMag) this.metaPar).getQuanMode();
		if (quanMode.equals(MetaConstants.labelFree)) {
			quanPepReader = new FlashLfqQuanPepReader(this.quan_pep_file, fileList);
		} else if (quanMode.equals(MetaConstants.isobaricLabel)) {
			IsobaricTag isoTag = ((MetaParameterMag) this.metaPar).getIsobaricTag();
			quanPepReader = new FlashLfqQuanPepReader(this.quan_pep_file, fileList, isoTag.getReportTitle());
		} else {
			LOGGER.error(getTaskName() + ": unknown quantitative type");
			System.err.println(format.format(new Date()) + "\t" + getTaskName() + ": unknown quantitative type");
			return false;
		}

		FlashLfqQuanPeptide[] quanPeps = quanPepReader.getMetaPeptides();
		for (int i = 0; i < quanPeps.length; i++) {
			String seq = quanPeps[i].getModSeq();
			if (pepExpectMap.containsKey(seq)) {
				quanPeps[i].setIndividualPEP(pepExpectMap.get(seq));
			} else {
				double[] expect = new double[fileIdMap.size()];
				Arrays.fill(expect, 1.0);
				quanPeps[i].setIndividualPEP(expect);
			}
		}

		HashMap<String, Double> genomeScoreMap = this.refineGenomes(quanPeps);

		String[] razorProteins = this.refineProteins(quanPeps, genomeScoreMap);

		setProgress(82);

		HashMap<String, double[]> genomeIntensityMap = new HashMap<String, double[]>();
		HashMap<String, double[]> proteinIntensityMap = new HashMap<String, double[]>();
		HashMap<String, Double> proteinTotalIntenMap = new HashMap<String, Double>();

		HashSet<String> pepSequenceSet = new HashSet<String>();
		HashMap<String, HashSet<String>> totalPepCountMap = new HashMap<String, HashSet<String>>();
		HashMap<String, HashSet<String>> razorPepCountMap = new HashMap<String, HashSet<String>>();

		ArrayList<FlashLfqQuanPeptide> peptideList = new ArrayList<FlashLfqQuanPeptide>();

		for (int i = 0; i < quanPeps.length; i++) {
			if (razorProteins[i] != null) {
				String stripeSequence = quanPeps[i].getSequence();

				pepSequenceSet.add(stripeSequence);
				peptideList.add(quanPeps[i]);

				if (razorPepCountMap.containsKey(razorProteins[i])) {
					razorPepCountMap.get(razorProteins[i]).add(stripeSequence);
				} else {
					HashSet<String> set = new HashSet<String>();
					set.add(stripeSequence);
					razorPepCountMap.put(razorProteins[i], set);
				}

				if (razorProteins[i].startsWith(magDb.getIdentifier())) {
					String genome = razorProteins[i].split("_")[0];
					if (razorPepCountMap.containsKey(genome)) {
						razorPepCountMap.get(genome).add(stripeSequence);
					} else {
						HashSet<String> set = new HashSet<String>();
						set.add(stripeSequence);
						razorPepCountMap.put(genome, set);
					}
				}
			}
		}

		for (int i = 0; i < quanPeps.length; i++) {
			String[] proteins = quanPeps[i].getProteins();
			int count = 0;
			for (int j = 0; j < proteins.length; j++) {
				if (razorPepCountMap.containsKey(proteins[j])) {
					count++;
				}
			}

			if (count > 0) {
				String stripeSequence = quanPeps[i].getSequence();
				double[] pepIntensity = new double[quanPeps[i].getIntensity().length];
				for (int j = 0; j < pepIntensity.length; j++) {
					pepIntensity[j] = quanPeps[i].getIntensity()[j] / (double) count;
				}

				for (int j = 0; j < proteins.length; j++) {
					if (razorPepCountMap.containsKey(proteins[j])) {
						if (proteins[j].startsWith(magDb.getIdentifier())) {
							String genome = proteins[j].split("_")[0];
							if (genomeScoreMap.containsKey(genome)) {
								double[] genomeIntensity;
								if (genomeIntensityMap.containsKey(genome)) {
									genomeIntensity = genomeIntensityMap.get(genome);
								} else {
									genomeIntensity = new double[pepIntensity.length];
									genomeIntensityMap.put(genome, genomeIntensity);
								}
								for (int k = 0; k < pepIntensity.length; k++) {
									genomeIntensity[k] += pepIntensity[k];
								}

								double proteinTotalInten;
								double[] proteinIntensity;
								if (proteinIntensityMap.containsKey(proteins[j])) {
									proteinIntensity = proteinIntensityMap.get(proteins[j]);
									proteinTotalInten = proteinTotalIntenMap.get(proteins[j]);
								} else {
									proteinIntensity = new double[pepIntensity.length];
									proteinIntensityMap.put(proteins[j], proteinIntensity);
									proteinTotalInten = 0.0;
								}
								for (int k = 0; k < pepIntensity.length; k++) {
									proteinIntensity[k] += pepIntensity[k];
									proteinTotalInten += pepIntensity[k];
								}
								proteinTotalIntenMap.put(proteins[j], proteinTotalInten);

								if (totalPepCountMap.containsKey(genome)) {
									totalPepCountMap.get(genome).add(stripeSequence);
								} else {
									HashSet<String> set = new HashSet<String>();
									set.add(stripeSequence);
									totalPepCountMap.put(genome, set);
								}

								if (totalPepCountMap.containsKey(proteins[j])) {
									totalPepCountMap.get(proteins[j]).add(stripeSequence);
								} else {
									HashSet<String> set = new HashSet<String>();
									set.add(stripeSequence);
									totalPepCountMap.put(proteins[j], set);
								}
							}
						} else {
							double proteinTotalInten;
							double[] proteinIntensity;
							if (proteinIntensityMap.containsKey(proteins[j])) {
								proteinIntensity = proteinIntensityMap.get(proteins[j]);
								proteinTotalInten = proteinTotalIntenMap.get(proteins[j]);
							} else {
								proteinIntensity = new double[pepIntensity.length];
								proteinIntensityMap.put(proteins[j], proteinIntensity);
								proteinTotalInten = 0.0;
							}
							for (int k = 0; k < pepIntensity.length; k++) {
								proteinIntensity[k] += pepIntensity[k];
								proteinTotalInten += pepIntensity[k];
							}
							proteinTotalIntenMap.put(proteins[j], proteinTotalInten);

							if (totalPepCountMap.containsKey(proteins[j])) {
								totalPepCountMap.get(proteins[j]).add(stripeSequence);
							} else {
								HashSet<String> set = new HashSet<String>();
								set.add(stripeSequence);
								totalPepCountMap.put(proteins[j], set);
							}
						}
					}
				}
			}
		}

		String[] proteins = proteinTotalIntenMap.keySet().toArray(String[]::new);
		Arrays.sort(proteins, new Comparator<String>() {

			@Override
			public int compare(String o1, String o2) {
				// TODO Auto-generated method stub
				if (proteinTotalIntenMap.get(o1) > proteinTotalIntenMap.get(o2))
					return -1;
				if (proteinTotalIntenMap.get(o1) < proteinTotalIntenMap.get(o2))
					return 1;
				return 0;
			}
		});

		HashMap<String, Integer> proteinIdMap = new HashMap<String, Integer>();
		for (int i = 0; i < proteins.length; i++) {
			proteinIdMap.put(proteins[i], i);
		}

		try (PrintWriter idWriter = new PrintWriter(new File(this.metaPar.getResult(), "id.tsv"))) {
			idWriter.println("Raw file\tUnique peptide count\tPSM count");

			String[] rawNames = fileIdMap.keySet().toArray(new String[fileIdMap.size()]);
			Arrays.sort(rawNames, new Comparator<String>() {

				@Override
				public int compare(String o1, String o2) {
					// TODO Auto-generated method stub
					return fileIdMap.get(o1) - fileIdMap.get(o2);
				}
			});

			for (int i = 0; i < rawNames.length; i++) {
				idWriter.println(rawNames[i] + "\t" + fileSequenceSets[i].size() + "\t" + filePsmCount[i]);
			}
			idWriter.close();
		} catch (IOException e) {
			System.out.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": error in writing identification information to "
					+ new File(this.metaPar.getResult(), "id.tsv"));
			LOGGER.error(getTaskName() + ": error in writing identification information to "
					+ new File(this.metaPar.getResult(), "id.tsv"), e);
		}

		setProgress(85);

		HashMap<String, ArrayList<Integer>> proPepIdMap = new HashMap<String, ArrayList<Integer>>();
		HashMap<String, ArrayList<Integer>> proPepRazorIdMap = new HashMap<String, ArrayList<Integer>>();
		HashMap<String, Double> proScoreMap = new HashMap<String, Double>();
		try {
			final_pep_txt = new File(this.metaPar.getResult(), "final_peptides.tsv");

			PrintWriter writer = new PrintWriter(final_pep_txt);

			StringBuilder titlesb = new StringBuilder();
			titlesb.append("Sequence").append("\t");
			titlesb.append("Base sequence").append("\t");
			titlesb.append("Length").append("\t");
			titlesb.append("Missed cleavages").append("\t");
			titlesb.append("Mass").append("\t");
			titlesb.append("Proteins").append("\t");
			titlesb.append("Charges").append("\t");
			titlesb.append("Score").append("\t");
			titlesb.append("PEP").append("\t");

			for (int i = 0; i < fileList.size(); i++) {
				titlesb.append("Identification type " + fileList.get(i)).append("\t");
			}

			for (int i = 0; i < fileList.size(); i++) {
				titlesb.append("MS/MS count " + fileList.get(i)).append("\t");
			}

			titlesb.append("Intensity").append("\t");

			if (quanMode.equals(MetaConstants.labelFree)) {
				for (int i = 0; i < fileList.size(); i++) {
					titlesb.append("Intensity " + fileList.get(i)).append("\t");
				}
			} else if (quanMode.equals(MetaConstants.isobaricLabel)) {
				IsobaricTag isoTag = ((MetaParameterMag) this.metaPar).getIsobaricTag();
				String[] repoTitle = isoTag.getReportTitle();
				for (int i = 0; i < fileList.size(); i++) {
					for (int j = 0; j < repoTitle.length; j++) {
						titlesb.append("Intensity " + fileList.get(i) + "_" + repoTitle[j]).append("\t");
					}
				}
			}

			titlesb.append("Reverse").append("\t");
			titlesb.append("Potential contaminant").append("\t");
			titlesb.append("id").append("\t");
			titlesb.append("Protein group IDs").append("\t");
			titlesb.append("MS/MS Count");

			writer.println(titlesb);

			int id = 0;
			for (FlashLfqQuanPeptide peptide : peptideList) {

				String fullSeq = peptide.getModSeq();
				String sequence = peptide.getSequence();

				FlashLfqPsm psm = flashLfqPsmMap.get(fullSeq);
				if (psm == null) {
					continue;
				}

				StringBuilder sb = new StringBuilder();
				sb.append(fullSeq).append("\t");
				sb.append(sequence).append("\t");
				sb.append(sequence.length()).append("\t");
				sb.append(psm.getMiss()).append("\t");
				sb.append(psm.getCalMass()).append("\t");

				String[] pros = peptide.getProteins();

				for (String pro : pros) {
					if (proteinIdMap.containsKey(pro)) {
						sb.append(pro).append(";");
					}
				}
				sb.deleteCharAt(sb.length() - 1);
				sb.append("\t");

				int[] charge = pepChargeMap.get(fullSeq);
				StringBuilder chargesb = new StringBuilder();
				for (int j = 0; j < charge.length; j++) {
					if (charge[j] == 1) {
						chargesb.append(j).append(";");
					}
				}
				if (chargesb.length() > 1) {
					sb.append(chargesb.subSequence(0, chargesb.length() - 1)).append("\t");
				} else {
					sb.append("2").append("\t");
				}

				sb.append(psm.getScore()).append("\t");
				sb.append(psm.getPEP()).append("\t");

				int totalMsCount = 0;
				int[] psmCount = pepPsmCountMap.get(fullSeq);
				double[] intensity = peptide.getIntensity();
				for (int j = 0; j < psmCount.length; j++) {
					if (psmCount[j] > 0) {
						totalMsCount += psmCount[j];
						sb.append("By MS/MS").append("\t");
					} else {
						if (intensity[j] > 0) {
							sb.append("MBR\t");
						} else {
							sb.append("NotDetected\t");
						}
					}
				}

				for (int j = 0; j < psmCount.length; j++) {
					sb.append(psmCount[j]).append("\t");
				}

				double totalIntensity = 0;
				for (int j = 0; j < intensity.length; j++) {
					totalIntensity += intensity[j];
				}

				sb.append(BigInteger.valueOf((long) Math.floor(totalIntensity))).append("\t");

				for (int j = 0; j < intensity.length; j++) {
					sb.append(BigInteger.valueOf((long) Math.floor(intensity[j]))).append("\t");
				}

				sb.append("").append("\t");
				sb.append("").append("\t");

				sb.append(id).append("\t");

				for (String pro : pros) {
					if (proteinIdMap.containsKey(pro)) {
						sb.append(proteinIdMap.get(pro)).append(";");

						if (proPepIdMap.containsKey(pro)) {
							proPepIdMap.get(pro).add(id);
						} else {
							ArrayList<Integer> list = new ArrayList<Integer>();
							list.add(id);
							proPepIdMap.put(pro, list);
						}

						if (proScoreMap.containsKey(pro)) {
							proScoreMap.put(pro, proScoreMap.get(pro) * psm.getPEP());
						} else {
							proScoreMap.put(pro, psm.getPEP());
						}

						HashSet<String> rSet = razorPepCountMap.get(pro);
						if (proPepRazorIdMap.containsKey(pro)) {
							if (rSet.contains(sequence)) {
								proPepRazorIdMap.get(pro).add(1);
							} else {
								proPepRazorIdMap.get(pro).add(0);
							}
						} else {
							ArrayList<Integer> list = new ArrayList<Integer>();
							if (rSet.contains(sequence)) {
								list.add(1);
							} else {
								list.add(0);
							}
							proPepRazorIdMap.put(pro, list);
						}
					}
				}
				sb.deleteCharAt(sb.length() - 1);
				sb.append("\t");

				sb.append(totalMsCount);

				writer.println(sb);
				id++;
			}

			writer.close();

		} catch (Exception e) {
			// TODO Auto-generated catch block
			LOGGER.error(getTaskName() + ": error in writing the peptides to  " + final_pep_txt, e);
			System.err.println(format.format(new Date()) + "\t" + getTaskName() + ": error in writing the peptides to  "
					+ final_pep_txt);
		}

		LOGGER.info(getTaskName() + ": peptide report has been exported to " + final_pep_txt);
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": peptide report has been exported to "
				+ final_pep_txt);

		if (!final_pep_txt.exists() || final_pep_txt.length() == 0) {
			return false;
		}

		setProgress(87);

		LOGGER.info(getTaskName() + ": exporting protein report started");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": exporting protein report started");

		ArrayList<MetaProtein> proteinList = new ArrayList<MetaProtein>();
		try {

			this.final_pro_txt = new File(this.metaPar.getResult(), "final_proteins.tsv");

			PrintWriter writer = new PrintWriter(this.final_pro_txt);

			StringBuilder titlesb = new StringBuilder();
			titlesb.append("Protein IDs").append("\t");
//			titlesb.append("Majority protein IDs").append("\t");
			titlesb.append("Peptide counts (all)").append("\t");
			titlesb.append("Peptide counts (razor)").append("\t");
			titlesb.append("Number of PSMs").append("\t");
			titlesb.append("PEP").append("\t");
			titlesb.append("Intensity").append("\t");

			if (quanMode.equals(MetaConstants.labelFree)) {
				for (int i = 0; i < fileList.size(); i++) {
					titlesb.append("Intensity " + fileList.get(i)).append("\t");
				}
			} else if (quanMode.equals(MetaConstants.isobaricLabel)) {
				IsobaricTag isoTag = ((MetaParameterMag) this.metaPar).getIsobaricTag();
				String[] repoTitle = isoTag.getReportTitle();
				for (int i = 0; i < fileList.size(); i++) {
					for (int j = 0; j < repoTitle.length; j++) {
						titlesb.append("Intensity " + fileList.get(i) + "_" + repoTitle[j]).append("\t");
					}
				}
			}

			titlesb.append("Reverse").append("\t");
			titlesb.append("Potential contaminant").append("\t");
			titlesb.append("id").append("\t");
			titlesb.append("Peptide IDs").append("\t");
			titlesb.append("Peptide is razor");

			writer.println(titlesb);

			DecimalFormat df = FormatTool.getDFE4();
			for (int i = 0; i < proteins.length; i++) {

				StringBuilder sb = new StringBuilder();
				sb.append(proteins[i]).append("\t");
				sb.append(totalPepCountMap.get(proteins[i]).size()).append("\t");
				sb.append(razorPepCountMap.get(proteins[i]).size()).append("\t");
				sb.append(psmCountMap.get(proteins[i])).append("\t");
				sb.append(proScoreMap.get(proteins[i])).append("\t");

				double totalIntensity = proteinTotalIntenMap.get(proteins[i]);
				sb.append(df.format(totalIntensity)).append("\t");

				double[] intensity = proteinIntensityMap.get(proteins[i]);
				for (int j = 0; j < intensity.length; j++) {
					sb.append(df.format(intensity[j])).append("\t");
				}
				sb.append("\t");
				sb.append("\t");
				sb.append(i).append("\t");

				ArrayList<Integer> pepIdList = proPepIdMap.get(proteins[i]);
				ArrayList<Integer> pepRazorList = proPepRazorIdMap.get(proteins[i]);

				for (int j = 0; j < pepIdList.size(); j++) {
					sb.append(pepIdList.get(j)).append(";");
				}
				sb.deleteCharAt(sb.length() - 1);
				sb.append("\t");

				for (int j = 0; j < pepRazorList.size(); j++) {
					if (pepRazorList.get(j) == 1) {
						sb.append("true;");
					} else {
						sb.append("false;");
					}
				}
				sb.deleteCharAt(sb.length() - 1);
				writer.println(sb);

				MetaProtein mp = new MetaProtein(i, 0, proteins[i], totalPepCountMap.get(proteins[i]).size(),
						psmCountMap.get(proteins[i]), proScoreMap.get(proteins[i]), intensity);
				proteinList.add(mp);
			}

			writer.close();

		} catch (IOException e) {
			// TODO Auto-generated catch block
			LOGGER.error(getTaskName() + ": error in writing final protein result to " + this.final_pro_txt.getName(),
					e);
			System.err.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": error in writing final protein result to " + this.final_pro_txt.getName());
		}

		LOGGER.info(getTaskName() + ": protein report has been exported to " + final_pro_txt.getName());
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": protein report has been exported to "
				+ final_pro_txt.getName());

		setProgress(89);

		LOGGER.info(getTaskName() + ": functional annotation started");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": functional annotation started");

		File funcFile = new File(metaPar.getResult(), "functional_annotation");
		if (!funcFile.exists()) {
			funcFile.mkdir();
		}

		File funTsv = new File(funcFile, "functions.tsv");
		File funReportTsv = new File(funcFile, "functions_report.tsv");
		File final_pro_xml_file = new File(funcFile, "functions.xml");

		HashMap<String, Integer> proTaxIdMap = new HashMap<String, Integer>();
		MetaProtein[] metapros = proteinList.toArray(new MetaProtein[proteinList.size()]);

		String[] intensityTitle;
		if (quanMode.equals(MetaConstants.labelFree)) {
			intensityTitle = fileList.toArray(String[]::new);
		} else if (quanMode.equals(MetaConstants.isobaricLabel)) {
			IsobaricTag isoTag = ((MetaParameterMag) this.metaPar).getIsobaricTag();
			String[] repoTitle = isoTag.getReportTitle();
			intensityTitle = new String[fileList.size() * repoTitle.length];

			for (int i = 0; i < fileList.size(); i++) {
				for (int j = 0; j < repoTitle.length; j++) {
					intensityTitle[i * repoTitle.length + j] = fileList.get(i) + "_" + repoTitle[j];
				}
			}
		} else {
			return false;
		}

		MetaProteinAnnoMag[] mpas = funcAnnoSqlite(metapros, intensityTitle);
		for (int i = 0; i < mpas.length; i++) {
			String name = mpas[i].getPro().getName();
			if (NumberUtils.isCreatable(mpas[i].getTax_id())) {
				proTaxIdMap.put(name, Integer.parseInt(mpas[i].getTax_id()));
			} else {
				proTaxIdMap.put(name, 131567);
			}
		}

		MetaProteinXMLReader2 funReader = null;
		try {
			funReader = new MetaProteinXMLReader2(final_pro_xml_file);
		} catch (DocumentException e) {
			// TODO Auto-generated catch block
			LOGGER.error(getTaskName() + ": error in reading protein function information from "
					+ final_pro_xml_file.getName(), e);
			System.err.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": error in reading protein function information from " + this.quan_pro_file.getName());
			return false;
		}

		if (!funTsv.exists() || !funReportTsv.exists()) {
			System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": writing results to "
					+ funTsv.getName() + " started");
			try {
				funReader.exportTsv(funTsv);
				funReader.exportTsvReport(funReportTsv);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				LOGGER.error(getTaskName() + ": error in writing protein function information to " + funTsv.getName(),
						e);
				System.err.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": error in writing protein function information to " + funTsv.getName());
			}
			System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": writing results to "
					+ funTsv.getName() + " finished");
		}

		File html = new File(funcFile, "functions.html");
		if (!html.exists()) {
			System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": writing results to "
					+ html.getName() + " started");
			try {
				funReader.exportHtml(html);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				LOGGER.error(getTaskName() + ": error in writing protein function information to " + html.getName(), e);
				System.err.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": error in writing protein function information to " + html.getName());
			}
			System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": writing results to "
					+ html.getName() + " finished");
		}

		setProgress(91);

		// taxonomy analysis

		LOGGER.info(getTaskName() + ": taxonomy analysis started");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": taxonomy analysis started");

		File taxFile = new File(metaPar.getResult(), "taxonomy_analysis");
		if (!taxFile.exists()) {
			taxFile.mkdir();
		}

		MagFuncSearcher searcher = null;
		try {
			searcher = new MagFuncSearcher(metaPar.getMicroDb(), mrp.getToolPath(ToolType.FUNCDB));

		} catch (NumberFormatException | SQLException | IOException e) {
			// TODO Auto-generated catch block
			LOGGER.error(getTaskName() + ": error in genome analysis, task failed", e);
			System.err.println(
					format.format(new Date()) + "\t" + getTaskName() + ": error in genome analysis, task failed");

			return false;
		}

		String[] genomes = genomeIntensityMap.keySet().toArray(new String[genomeIntensityMap.size()]);
		HashMap<String, String[]> genomeTaxaMap = null;
		try {
			genomeTaxaMap = searcher.matchGenome(genomes);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			LOGGER.error(getTaskName() + ": Error in searching the genome2taxa database", e);
			System.err.println(
					format.format(new Date()) + "\t" + getTaskName() + ": error in searching the genome2taxa database");
			return false;
		}

		File genomeHtmlFile = new File(taxFile, "Genome.html");
		try {
			HMGenomeHtmlWriter htmlWriter = new HMGenomeHtmlWriter(genomeHtmlFile);
			htmlWriter.write(genomeTaxaMap);
			htmlWriter.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			LOGGER.error("Error in writing genome information to " + genomeHtmlFile, e);
			System.err.println(
					format.format(new Date()) + "\t" + "Error in writing genome information to " + genomeHtmlFile);
			return false;
		}

		TaxonomyRanks[] mainRanks = TaxonomyRanks.getMainRanks();
		double[] cellularIntensity = new double[intensityTitle.length];
		HashMap<String, String[]>[] taxRankMaps = new HashMap[mainRanks.length];
		HashMap<String, double[]>[] taxIntensityMaps = new HashMap[mainRanks.length];
		HashMap<String, Double>[] taxScoreMaps = new HashMap[mainRanks.length];
		for (int i = 0; i < taxRankMaps.length; i++) {
			taxRankMaps[i] = new HashMap<String, String[]>();
			taxIntensityMaps[i] = new HashMap<String, double[]>();
			taxScoreMaps[i] = new HashMap<String, Double>();
		}

		File genomeFile = new File(taxFile, "Genome.tsv");
		PrintWriter genomeWriter = null;
		try {
			genomeWriter = new PrintWriter(genomeFile);
			StringBuilder genomeTitlesb = new StringBuilder();
			genomeTitlesb.append("Genome\t");
			genomeTitlesb.append("Minus_log(p-value)\t");
			genomeTitlesb.append("PSM count\t");
			genomeTitlesb.append("Peptide count\t");
			genomeTitlesb.append("Razor peptide count\t");
			genomeTitlesb.append("GTDB_Name\t");
			genomeTitlesb.append("GTDB_Rank\t");

			for (int i = 0; i < mainRanks.length; i++) {
				genomeTitlesb.append("GTDB_").append(mainRanks[i].getName()).append("\t");
			}

			genomeTitlesb.append("NCBI_Name\t");
			genomeTitlesb.append("NCBI_Rank");

			for (int i = 0; i < mainRanks.length; i++) {
				genomeTitlesb.append("\t").append("NCBI_").append(mainRanks[i].getName());
			}

			for (String exp : intensityTitle) {
				if (exp.startsWith("Intensity ")) {
					genomeTitlesb.append("\t").append(exp);
				} else {
					genomeTitlesb.append("\tIntensity ").append(exp);
				}
			}

			genomeWriter.println(genomeTitlesb);

			setProgress(94);

			for (int i = 0; i < genomes.length; i++) {
				String[] taxaContent = genomeTaxaMap.get(genomes[i]);
				if (genomeIntensityMap.containsKey(genomes[i])) {
					double[] genomeIntensity = genomeIntensityMap.get(genomes[i]);
					for (int j = 0; j < genomeIntensity.length; j++) {
						cellularIntensity[j] += genomeIntensity[j];
					}

					double score = genomeScoreMap.get(genomes[i]);
					StringBuilder genomeTaxaSb = new StringBuilder();
					genomeTaxaSb.append(genomes[i]).append("\t");

					genomeTaxaSb.append(FormatTool.getDF2().format(-Math.log10(score))).append("\t");

					if (psmCountMap.containsKey(genomes[i])) {
						genomeTaxaSb.append(psmCountMap.get(genomes[i])).append("\t");
					} else {
						genomeTaxaSb.append("0").append("\t");
					}
					if (totalPepCountMap.containsKey(genomes[i])) {
						genomeTaxaSb.append(totalPepCountMap.get(genomes[i]).size()).append("\t");
					} else {
						genomeTaxaSb.append("0").append("\t");
					}
					if (razorPepCountMap.containsKey(genomes[i])) {
						genomeTaxaSb.append(razorPepCountMap.get(genomes[i]).size());
					} else {
						genomeTaxaSb.append("0");
					}

					for (int j = 1; j < taxaContent.length; j++) {
						genomeTaxaSb.append("\t").append(taxaContent[j]);
					}

					for (int j = 0; j < genomeIntensity.length; j++) {
						genomeTaxaSb.append("\t").append(genomeIntensity[j]);
					}
					genomeWriter.println(genomeTaxaSb);

					for (int j = 0; j < taxRankMaps.length; j++) {
						String taxon = taxaContent[j + 3];
						if (taxon == null || taxon.length() == 0) {
							continue;
						}
						if (!taxRankMaps[j].containsKey(taxon)) {
							String[] lineage = new String[j + 1];
							System.arraycopy(taxaContent, 3, lineage, 0, lineage.length);
							taxRankMaps[j].put(taxon, lineage);
						}

						double[] taxIntensity;
						if (taxIntensityMaps[j].containsKey(taxon)) {
							taxIntensity = taxIntensityMaps[j].get(taxon);
						} else {
							taxIntensity = new double[genomeIntensity.length];
						}

						for (int k = 0; k < taxIntensity.length; k++) {
							taxIntensity[k] += genomeIntensity[k];
						}
						taxIntensityMaps[j].put(taxon, taxIntensity);

						if (taxScoreMaps[j].containsKey(taxon)) {
							taxScoreMaps[j].put(taxon, taxScoreMaps[j].get(taxon) - Math.log10(score));
						} else {
							taxScoreMaps[j].put(taxon, -Math.log10(score));
						}
					}
				}
			}
			genomeWriter.close();

		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			LOGGER.error("Error in writing genome information to " + genomeFile, e);
			System.err
					.println(format.format(new Date()) + "\t" + "Error in writing genome information to " + genomeFile);
			return false;
		}

		File taxaFile = new File(taxFile, "Taxa.tsv");
		File taxaReportFile = new File(taxFile, "Taxa_report.tsv");
		PrintWriter taxaWriter = null;
		PrintWriter taxaReportWriter = null;
		try {
			taxaWriter = new PrintWriter(taxaFile);
			taxaReportWriter = new PrintWriter(taxaReportFile);

			StringBuilder taxaTitlesb1 = new StringBuilder();
			taxaTitlesb1.append("Name\t");
			taxaTitlesb1.append("Minus_log(p-value)\t");
			taxaTitlesb1.append("Rank");

			for (int i = 0; i < mainRanks.length; i++) {
				taxaTitlesb1.append("\t").append(mainRanks[i].getName());
			}
			for (String exp : intensityTitle) {
				if (exp.startsWith("Intensity ")) {
					taxaTitlesb1.append("\t").append(exp);
				} else {
					taxaTitlesb1.append("\tIntensity ").append(exp);
				}
			}
			StringBuilder taxaTitlesb2 = new StringBuilder();
			taxaTitlesb2.append("Name\t");
			taxaTitlesb2.append("Rank");

			for (int i = 0; i < mainRanks.length; i++) {
				taxaTitlesb2.append("\t").append(mainRanks[i].getName());
			}
			for (String exp : intensityTitle) {
				if (exp.startsWith("Intensity ")) {
					taxaTitlesb2.append("\t").append(exp);
				} else {
					taxaTitlesb2.append("\tIntensity ").append(exp);
				}
			}
			taxaWriter.println(taxaTitlesb1);
			taxaReportWriter.println(taxaTitlesb2);

			for (int i = 0; i < taxRankMaps.length; i++) {
				for (String taxon : taxRankMaps[i].keySet()) {
					String[] taxaContent = taxRankMaps[i].get(taxon);
					double[] intensity = taxIntensityMaps[i].get(taxon);

					StringBuilder taxaSb = new StringBuilder();
					taxaSb.append(taxon).append("\t");

					double taxonScore = taxScoreMaps[i].get(taxon);
					taxaSb.append(FormatTool.getDF2().format(taxonScore)).append("\t");

					taxaSb.append(mainRanks[i].getName());

					for (int j = 0; j < mainRanks.length; j++) {
						if (j < taxaContent.length) {
							taxaSb.append("\t").append(taxaContent[j]);
						} else {
							taxaSb.append("\t");
						}
					}

					for (int j = 0; j < intensity.length; j++) {
						taxaSb.append("\t").append(intensity[j]);
					}

					taxaWriter.println(taxaSb);
				}
			}

			taxaWriter.close();

			for (int i = 0; i < taxRankMaps.length; i++) {
				for (String taxon : taxRankMaps[i].keySet()) {
					String[] taxaContent = taxRankMaps[i].get(taxon);
					double[] intensity = taxIntensityMaps[i].get(taxon);

					StringBuilder taxaSb = new StringBuilder();
					taxaSb.append(taxon).append("\t");
					taxaSb.append(mainRanks[i].getName());

					for (int j = 0; j < mainRanks.length; j++) {
						if (j < taxaContent.length) {
							taxaSb.append("\t").append(taxaContent[j]);
						} else {
							taxaSb.append("\t");
						}
					}
					for (int j = 0; j < intensity.length; j++) {
						taxaSb.append("\t").append(intensity[j]);
					}

					taxaReportWriter.println(taxaSb);
				}
			}

			taxaReportWriter.close();

		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			LOGGER.error("Error in writing taxa information to " + taxaFile, e);
			System.err.println(format.format(new Date()) + "\t" + "Error in writing taxa information to " + taxaFile);
		}

		LOGGER.info(getTaskName() + ": taxonomy analysis finished");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": taxonomy analysis finished");

		setProgress(95);

		progressBar.setString(getTaskName() + ": exporting report");

		LOGGER.info(getTaskName() + ": exporting report started");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": exporting report started");

		File reportDir = new File(this.metaPar.getResult(), "report");
		if (!reportDir.exists()) {
			reportDir.mkdirs();
		}

		MetaReportCopyTask.copy(reportDir, getVersion());

		metaPar.setReportDir(reportDir);

		File reportHtmlDir = new File(reportDir + "\\reports\\rmd_html");

		if (!reportHtmlDir.exists()) {
			reportHtmlDir.mkdirs();
		}

		File jsOutput = new File(metaPar.getReportDir() + "\\reports\\metamap\\data\\tree.js");
		if (!jsOutput.getParentFile().exists()) {
			jsOutput.getParentFile().mkdirs();
		}

		if (!jsOutput.exists() || jsOutput.length() == 0) {
			MetaTreeHandler.exportJS(taxRankMaps, taxIntensityMaps, cellularIntensity, quanExpNames, jsOutput);
		}

		MetaPythonReportTask task = new MetaPythonReportTask(metaPar.getThreadCount());

		File summaryMetaFile = new File(this.metaPar.getResult(), "metainfo.tsv");
		try {
			metaPar.getMetadata().exportMetadata(summaryMetaFile);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		this.exportSummary(task);

		if (!final_pep_txt.exists() || final_pep_txt.length() == 0) {
			return false;
		} else {
			File report_peptides_summary = new File(reportHtmlDir, "report_peptides_summary.html");
			if (!report_peptides_summary.exists() || report_peptides_summary.length() == 0) {
				if (summaryMetaFile.exists()) {
					task.addTask(MetaPythonReportTask.PEPTIDE, final_pep_txt, report_peptides_summary, summaryMetaFile);
				} else {
					task.addTask(MetaPythonReportTask.PEPTIDE, final_pep_txt, report_peptides_summary);
				}
			}
		}

		if (!final_pro_txt.exists() || final_pro_txt.length() == 0) {
			return false;
		} else {
			File report_proteinGroups_summary = new File(reportHtmlDir, "report_proteinGroups_summary.html");
			if (!report_proteinGroups_summary.exists() || report_proteinGroups_summary.length() == 0) {
				if (summaryMetaFile.exists()) {
					task.addTask(MetaPythonReportTask.PROTEIN, final_pro_txt, report_proteinGroups_summary, summaryMetaFile);
				} else {
					task.addTask(MetaPythonReportTask.PROTEIN, final_pro_txt, report_proteinGroups_summary);
				}
			}
		}

		if (!taxaReportFile.exists() || taxaReportFile.length() == 0) {
			return false;
		} else {
			File report_taxonomy_summary = new File(reportHtmlDir, "report_taxonomy_summary.html");
			if (!report_taxonomy_summary.exists() || report_taxonomy_summary.length() == 0) {

				if (summaryMetaFile != null && summaryMetaFile.exists()) {
					task.addTask(MetaPythonReportTask.TAXONOMY, taxaReportFile, report_taxonomy_summary, summaryMetaFile);
				} else {
					task.addTask(MetaPythonReportTask.TAXONOMY, taxaReportFile, report_taxonomy_summary);
				}
			}
		}

		if (!funReportTsv.exists() || funReportTsv.length() == 0) {
			return false;
		} else {
			File report_function_summary = new File(reportHtmlDir, "report_function_summary.html");
			if (!report_function_summary.exists() || report_function_summary.length() == 0) {

				if (summaryMetaFile != null && summaryMetaFile.exists()) {
					task.addTask(MetaPythonReportTask.FUNCTION, funReportTsv, report_function_summary, summaryMetaFile);
				} else {
					task.addTask(MetaPythonReportTask.FUNCTION, funReportTsv, report_function_summary);
				}
			}
		}

		task.run();

		try {
			boolean finish = task.get();
			if (finish) {

				LOGGER.info(getTaskName() + ": report generation successes");
				System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": report generation successes");

				setProgress(100);
				return true;
			}

		} catch (InterruptedException | ExecutionException e) {
			// TODO Auto-generated catch block
			LOGGER.error(getTaskName() + ": error in the report generation", e);
			System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": error in the report generation");
		}

		return true;
	}

	@Override
	protected String getVersion() {
		// TODO Auto-generated method stub
		return null;
	}
}
