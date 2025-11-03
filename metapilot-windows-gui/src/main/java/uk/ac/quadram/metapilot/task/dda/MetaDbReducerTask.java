package uk.ac.quadram.metapilot.task.dda;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileFilter;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

import org.apache.commons.io.FileUtils;
import org.apache.logging.log4j.Logger;

import uk.ac.quadram.metapilot.core.prodb.FastaReader;
import uk.ac.quadram.metapilot.core.prodb.ProteinItem;
import uk.ac.quadram.metapilot.dbSearch.PepProGroup;
import uk.ac.quadram.metapilot.dbSearch.dbReducer.DBReducerTask;
import uk.ac.quadram.metapilot.dbSearch.pfind.PFindPSM;
import uk.ac.quadram.metapilot.dbSearch.pfind.PFindPsmReader;
import uk.ac.quadram.metapilot.dbSearch.pfind.PFindTask;
import javafx.scene.control.ProgressBar;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters.ToolType;
import uk.ac.quadram.metapilot.spectra.cluster.ClusterTask;
import uk.ac.quadram.metapilot.task.MetaDbCreateTask;
import uk.ac.quadram.metapilot.task.par.MetaBaseParameter;
import uk.ac.quadram.metapilot.task.par.MetaConstants;

/**
 * Generate sample-specific database in pFind workflow
 * 
 * @author Kai Cheng
 * @version 2025.03.20
 */
public class MetaDbReducerTask extends MetaDbCreateTask {

	private String[] mgfFiles;
	private String[] pf2Files;
	private File dbReducerExe;
	private File pFindExe;
	private boolean useDBReducer = false;
	private double psmFdr = 0.01;

	public MetaDbReducerTask(MetaBaseParameter metaPar, MetaResourceParameters mrp, ProgressBar progressBar,
			Logger logger) {
		super(metaPar, mrp, progressBar, logger);

		logger.info("{}: MetaDbReducerTask constructor called", taskName);

		// Log metadata information
		if (metaPar.getMetadata() != null) {
			String[] rawFiles = metaPar.getMetadata().getRawFiles();
			logger.info("{}: Constructor - Found {} raw files in metadata", taskName,
					rawFiles != null ? rawFiles.length : 0);

			if (rawFiles != null) {
				for (int i = 0; i < rawFiles.length; i++) {
					File rawFile = new File(rawFiles[i]);
					logger.debug("{}: Constructor - Raw file {}: {} (exists: {}, size: {} bytes)", taskName, i,
							rawFiles[i], rawFile.exists(), rawFile.length());
				}
			}
		} else {
			logger.error("{}: Constructor - Metadata is null in MetaBaseParameter", taskName);
		}

		// Log tool paths
		this.dbReducerExe = new File(mrp.getToolPath(ToolType.DBREDUCER));
		this.pFindExe = new File(mrp.getToolPath(ToolType.PFIND));

		logger.info("{}: DBReducer executable: {} (exists: {})", taskName, dbReducerExe.getAbsolutePath(),
				dbReducerExe.exists());
		logger.info("{}: pFind executable: {} (exists: {})", taskName, pFindExe.getAbsolutePath(), pFindExe.exists());

		if (dbReducerExe.getName().equals("DBReducer.exe") && dbReducerExe.exists()
				&& metaPar.getMs2ScanMode().equals(MetaConstants.HCD_FTMS)) {
			useDBReducer = true;
			logger.info("{}: Using DBReducer for database reduction", taskName);
		} else {
			logger.info("{}: Using pFind for database reduction (DBReducer conditions not met)", taskName);
			logger.debug("{}: DBReducer conditions - name: {}, exists: {}, ms2Mode: {}", taskName,
					dbReducerExe.getName(), dbReducerExe.exists(), metaPar.getMs2ScanMode());
		}

		// Initialize file arrays directly in constructor
		logger.info("{}: Initializing file arrays in constructor", taskName);

		// Log MS2 scan mode
		String ms2Type = this.metaPar.getMs2ScanMode();
		logger.info("{}: MS2 scan mode: {}", taskName, ms2Type);

		// Log spectra file directory
		File spectraDirFile = metaPar.getSpectraFile();
		logger.info("{}: Spectra directory: {} (exists: {}, isDirectory: {})", taskName,
				spectraDirFile.getAbsolutePath(), spectraDirFile.exists(), spectraDirFile.isDirectory());

		// Initialize file arrays with safety checks
		if (metaPar.getMetadata() == null) {
			logger.error("{}: Cannot initialize file arrays - metadata is null", taskName);
			this.mgfFiles = new String[0];
			this.pf2Files = new String[0];
			return;
		}

		String[] rawFiles = metaPar.getMetadata().getRawFiles();
		if (rawFiles == null || rawFiles.length == 0) {
			logger.error("{}: Cannot initialize file arrays - raw files array is null or empty", taskName);
			this.mgfFiles = new String[0];
			this.pf2Files = new String[0];
			return;
		}

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

		logger.info("{}: Using suffix: '{}' for MS2 type: {}", taskName, suffix, ms2Type);

		this.mgfFiles = new String[rawFiles.length];
		this.pf2Files = new String[rawFiles.length];

		logger.info("{}: Initializing {} MGF and PF2 files", taskName, rawFiles.length);

		for (int i = 0; i < rawFiles.length; i++) {
			File rawi = new File(rawFiles[i]);
			String rawName = rawi.getName();
			rawName = rawName.substring(0, rawName.length() - ".raw".length());

			File mgfFile = new File(spectraDirFile, rawName + ".mgf");
			File pf2 = new File(spectraDirFile, rawName + suffix + ".pf2");

			mgfFiles[i] = mgfFile.getAbsolutePath();
			pf2Files[i] = pf2.getAbsolutePath();

			logger.info(
					"{}: File {} - Raw: {}, MGF: {} (exists: {}, size: {} bytes), PF2: {} (exists: {}, size: {} bytes)",
					taskName, i, rawName, mgfFile.getName(), mgfFile.exists(), mgfFile.length(), pf2.getName(),
					pf2.exists(), pf2.length());
		}

		logger.info("{}: Constructor completed. MGF files: {}, PF2 files: {}", taskName, mgfFiles.length,
				pf2Files.length);
	}

	/**
	 * Clustering strategy
	 * 
	 * @throws Exception
	 */
	protected void cluster() {

		updateMessage(taskName + ": spectra clustering...");

		logger.info(taskName + ": spectra cluster started");
		System.out.println(format.format(new Date()) + "\t" + taskName + ": spectra cluster started");

		// Update progress at start of clustering
		updateProgress(20, 100);
		updateMessage("Starting spectra clustering...");

		// Test CDF resource loading to ensure the resource issue is resolved
		try {
			logger.info("{}: Testing CDF resource loading...", taskName);
			uk.ac.ebi.pride.spectracluster.cdf.CumulativeDistributionFunction cdf = uk.ac.ebi.pride.spectracluster.cdf.CumulativeDistributionFunctionFactory
					.getDefaultCumlativeDistributionFunctionForSimilarityMetric(
							uk.ac.ebi.pride.spectracluster.similarity.CombinedFisherIntensityTest.class);
			logger.info("{}: CDF loaded successfully: {}", taskName, cdf);
		} catch (Exception e) {
			logger.warn("{}: CDF loading test failed, but continuing with clustering: {}", taskName, e.getMessage());
		}

		// Log the MGF files that will be used for clustering
		logger.info("{}: Clustering {} MGF files", taskName, mgfFiles != null ? mgfFiles.length : 0);
		if (mgfFiles == null || mgfFiles.length == 0) {
			logger.error("{}: No MGF files available for clustering! Check if raw files were properly loaded.",
					taskName);
			logger.error("{}: This usually means the MetaData object was not properly set on the parameters.",
					taskName);
			return;
		}

		for (int i = 0; i < mgfFiles.length; i++) {
			File mgfFile = new File(mgfFiles[i]);
			logger.info("{}: MGF file {}: {} (exists: {}, size: {} bytes)", taskName, i, mgfFiles[i], mgfFile.exists(),
					mgfFile.length());
		}

		// Update progress before creating cluster directory
		updateProgress(22, 100);
		updateMessage("Preparing clustering directory...");

		File clusterDir = new File(dbcFile, "clustered_spectra");
		File clusterMgf = new File(clusterDir, "cluster.mgf");

		logger.info("{}: Cluster directory: {} (exists: {})", taskName, clusterDir.getAbsolutePath(),
				clusterDir.exists());
		logger.info("{}: Target cluster MGF: {} (exists: {}, size: {} bytes)", taskName, clusterMgf.getAbsolutePath(),
				clusterMgf.exists(), clusterMgf.length());

		boolean created = clusterDir.mkdir();
		logger.info("{}: Created cluster directory: {}", taskName, created);

		if (clusterMgf.exists()) {
			logger.info(taskName + ": clustered MS/MS spectra file already existed in " + clusterMgf);
			System.out.println(format.format(new Date()) + "\t" + taskName
					+ ": clustered MS/MS spectra file already existed in " + clusterMgf);
		} else {

			// Update progress before starting clustering
			updateProgress(25, 100);
			updateMessage("Clustering MS/MS spectra...");

			logger.info(taskName + ": clustering MS/MS spectra to " + clusterMgf + " started");
			System.out.println(format.format(new Date()) + "\t" + taskName + ": clustering MS/MS spectra to "
					+ clusterMgf + " started");

			logger.info("{}: Creating ClusterTask with threadCount={}, fastMode=true, isobaricTag={}", taskName,
					metaPar.getThreadCount(), metaPar.getIsobaricTag());

			ClusterTask clusterTask = new ClusterTask(metaPar.getThreadCount(), true, metaPar.getIsobaricTag(), logger);

			logger.info("{}: Starting clustering with ClusterTask.runCli()", taskName);
			clusterTask.runCli(mgfFiles, clusterDir, clusterMgf);

			logger.info(taskName + ": clustering MS/MS spectra to " + clusterMgf + " finished");
			System.out.println(format.format(new Date()) + "\t" + taskName + ": clustering MS/MS spectra to "
					+ clusterMgf + " finished");
		}

		// Update progress after clustering completion
		updateProgress(30, 100);
		updateMessage("Spectra clustering completed!");

		logger.info(taskName + ": spectra cluster finished");
		System.out.println(format.format(new Date()) + "\t" + taskName + ": spectra cluster finished");

		updateProgress(30, 100);
	}

	/**
	 * Clustering strategy
	 * 
	 * 
	 */
	protected File[] clusterSearch() {

		updateMessage(taskName + ": searching clustered spectra against the original database");

		logger.info(taskName + ": searching clustered spectra against the original database started");
		System.out.println(format.format(new Date()) + "\t" + taskName
				+ ": searching clustered spectra in original database started");

		// Update progress at start of search
		updateProgress(32, 100);
		updateMessage("Searching clustered spectra against database...");

		File clusterDir = new File(dbcFile, "clustered_spectra");

		File[] mgfs = clusterDir.listFiles(new FileFilter() {

			@Override
			public boolean accept(File arg0) {
				// TODO Auto-generated method stub
				if (arg0.getName().endsWith("mgf") || arg0.getName().endsWith("MGF"))
					return true;

				return false;
			}
		});

		// Safety check: if no MGF files were found, return empty array
		if (mgfs == null || mgfs.length == 0) {
			logger.warn("{}: No MGF files found in cluster directory: {}", taskName, clusterDir.getAbsolutePath());
			logger.warn("{}: This may indicate that clustering failed or no spectra were generated.", taskName);
			return new File[0];
		}

		// Update progress before starting search
		updateProgress(34, 100);
		updateMessage("Found " + mgfs.length + " MGF files, starting database search...");

		if (useDBReducer) {
			File reducedDbFile = new File(dbcFile, "DBReducer.fasta");
			if (reducedDbFile.exists() && reducedDbFile.length() > 0) {

				logger.info(taskName + ": refined database already existed in " + reducedDbFile);
				System.out.println(format.format(new Date()) + "\t" + taskName
						+ ": refined database already existed in " + reducedDbFile);

			} else {
				if (reducedDbFile.exists() && reducedDbFile.length() == 0) {
					FileUtils.deleteQuietly(reducedDbFile);
				}

				logger.info(taskName + ": generating refined database from the original database started");
				System.out.println(format.format(new Date()) + "\t" + taskName
						+ ": generating refined database from the original database started");

				String ms2Type = metaPar.getMs2ScanMode();

				String[] mgfFilePaths = new String[mgfs.length];
				for (int i = 0; i < mgfFilePaths.length; i++) {
					mgfFilePaths[i] = mgfs[i].getAbsolutePath();
				}

				logger.info("{}: Starting DBReducer search with {} MGF files", taskName, mgfFilePaths.length);

				DBReducerTask task = new DBReducerTask(dbReducerExe);
				task.setTotalThread(metaPar.getThreadCount());

				task.searchMgf(mgfFilePaths, dbcFile.getAbsolutePath(), originalFasta, ms2Type, psmFdr,
						metaPar.getFixMods(), metaPar.getVariMods());
				task.run(1, mgfFilePaths.length * 48);
			}

			if (reducedDbFile.exists() && reducedDbFile.length() > 0) {
				logger.info(taskName + ": generating refined database from the original database finished");
				System.out.println(format.format(new Date()) + "\t" + taskName
						+ ": generating refined database from the original database finished");
			}

			PFindTask.cleanPFindResultFolder(dbcFile, true);

			updateProgress(35, 100);

			return new File[] { reducedDbFile };

		} else {
			File pSpecResult = new File(dbcFile, "pFind.spectra");

			if (pSpecResult.exists() && pSpecResult.length() > 0) {

				logger.info(taskName + ": search results already existed in " + pSpecResult);
				System.out.println(format.format(new Date()) + "\t" + taskName + ": search results already existed in "
						+ pSpecResult);

			} else {

				if (pSpecResult.exists() && pSpecResult.length() == 0) {
					FileUtils.deleteQuietly(pSpecResult);
				}

				logger.info(taskName + ": searching clustered spectra in original database started");
				System.out.println(format.format(new Date()) + "\t" + taskName
						+ ": searching clustered spectra in original database started");

				String ms2Type = metaPar.getMs2ScanMode();

				String[] mgfFilePaths = new String[mgfs.length];
				for (int i = 0; i < mgfFilePaths.length; i++) {
					mgfFilePaths[i] = mgfs[i].getAbsolutePath();
				}

				logger.info("{}: Starting pFind search with {} MGF files", taskName, mgfFilePaths.length);

				File pFindBin = pFindExe.getParentFile();

				PFindTask task = new PFindTask(pFindBin);
				task.setTotalThread(metaPar.getThreadCount());

				task.searchMgf(mgfFilePaths, dbcFile.getAbsolutePath(), originalFasta, ms2Type, false, false, psmFdr,
						metaPar.getFixMods(), metaPar.getVariMods(), metaPar.getEnzyme(), metaPar.getMissCleavages(),
						metaPar.getDigestMode(), metaPar.getIsobaricTag(), true);
				task.run(1, mgfFilePaths.length * 48);
			}

			if (pSpecResult.exists() && pSpecResult.length() > 0) {
				logger.info(taskName + ": searching clustered spectra against the original database finished");
				System.out.println(format.format(new Date()) + "\t" + taskName
						+ ": searching clustered spectra against the original database finished");
			}

			PFindTask.cleanPFindResultFolder(pSpecResult.getParentFile(), true);

			updateProgress(35, 100);

			return new File[] { pSpecResult };
		}
	}

	@Override
	protected File[] metaProIQ() {
		// TODO Auto-generated method stub

		updateMessage(taskName + ": searching clustered spectra in original database");

		logger.info(taskName + ": MetaProIQ started");
		System.out.println(format.format(new Date()) + "\t" + taskName + ": MetaProIQ started");

		// Safety check: ensure we have MGF files to work with
		if (this.mgfFiles == null || this.mgfFiles.length == 0) {
			logger.error("{}: No MGF files available for MetaProIQ workflow! Check if raw files were properly loaded.",
					taskName);
			logger.error("{}: This usually means the MetaData object was not properly set on the parameters.",
					taskName);
			return new File[0];
		}

		String ms2Type = metaPar.getMs2ScanMode();

		File firstSearch = new File(dbcFile, "first_search");
		if (!firstSearch.exists()) {
			firstSearch.mkdir();
		}

		File secondSearch = new File(dbcFile, "second_search");
		if (!secondSearch.exists()) {
			secondSearch.mkdir();
		}

		File[] secondResults = new File[this.mgfFiles.length];

		int firstTaskCount = 0;
		for (int i = 0; i < this.mgfFiles.length; i++) {
			File mgfFilei = new File(this.mgfFiles[i]);
			String name = mgfFilei.getName();
			name = name.substring(0, name.lastIndexOf("."));

			File firstFasta = new File(firstSearch, name + ".fasta");
			File firstDir = new File(firstSearch, name);
			if (!firstDir.exists()) {
				firstDir.mkdir();
			}

			if (!firstFasta.exists() || firstFasta.length() == 0) {
				File pSpecResult = new File(firstDir, "pFind.spectra");
				if (!pSpecResult.exists() || pSpecResult.length() == 0) {
					firstTaskCount++;
				}
			}
		}

		if (firstTaskCount > 0) {

			int totalThreadCount = metaPar.getThreadCount();
			int threadCount = totalThreadCount > firstTaskCount ? totalThreadCount / firstTaskCount : 1;
			int threadPoolCount = totalThreadCount > firstTaskCount ? firstTaskCount : totalThreadCount;

			if (this.useDBReducer) {
				DBReducerTask task = new DBReducerTask(dbReducerExe);
				task.setTotalThread(threadCount);

				for (int i = 0; i < this.mgfFiles.length; i++) {

					File specFilei = new File(this.mgfFiles[i]);
					String name = specFilei.getName();
					name = name.substring(0, name.lastIndexOf("."));

					File firstFasta = new File(firstSearch, name + ".fasta");
					File firstDir = new File(firstSearch, name);
					if (!firstDir.exists()) {
						firstDir.mkdir();
					}

					if (!firstFasta.exists() || firstFasta.length() == 0) {

						File pSpecResult = new File(firstDir, "pFind.spectra");

						if (!pSpecResult.exists() || pSpecResult.length() == 0) {

							logger.info(taskName + ": first search of " + name + " in original database started");
							System.out.println(format.format(new Date()) + "\t" + taskName + ": first search of " + name
									+ " in original database started...");

							task.searchMgf(new String[] { this.mgfFiles[i] }, firstDir.getAbsolutePath(), originalFasta,
									ms2Type, 0.01, metaPar.getFixMods(), metaPar.getVariMods());
						}
					}
				}

				task.run(threadPoolCount, firstTaskCount * 12);

			} else {

				File pFindBin = pFindExe.getParentFile();
				PFindTask task = new PFindTask(pFindBin);
				task.setTotalThread(threadCount);

				for (int i = 0; i < this.mgfFiles.length; i++) {

					File mgfFilei = new File(this.mgfFiles[i]);
					String name = mgfFilei.getName();
					name = name.substring(0, name.lastIndexOf("."));

					File firstFasta = new File(firstSearch, name + ".fasta");
					File firstDir = new File(firstSearch, name);
					if (!firstDir.exists()) {
						firstDir.mkdir();
					}

					if (!firstFasta.exists() || firstFasta.length() == 0) {

						File pSpecResult = new File(firstDir, "pFind.spectra");

						if (!pSpecResult.exists() || pSpecResult.length() == 0) {

							logger.info(taskName + ": first search of " + name + " in original database started");
							System.out.println(format.format(new Date()) + "\t" + taskName + ": first search of " + name
									+ " in original database started...");

							if ((new File(this.pf2Files[i])).exists()) {
								task.searchPF2(new String[] { this.pf2Files[i] }, firstDir.getAbsolutePath(),
										originalFasta, ms2Type, false, false, 0.01, metaPar.getFixMods(),
										metaPar.getVariMods(), metaPar.getEnzyme(), metaPar.getMissCleavages(),
										metaPar.getDigestMode(), metaPar.getIsobaricTag(), true);
							} else {

								task.searchMgf(new String[] { this.mgfFiles[i] }, firstDir.getAbsolutePath(),
										originalFasta, ms2Type, false, false, 0.01, metaPar.getFixMods(),
										metaPar.getVariMods(), metaPar.getEnzyme(), metaPar.getMissCleavages(),
										metaPar.getDigestMode(), metaPar.getIsobaricTag(), true);
							}
						}
					}
				}

				task.run(threadPoolCount, firstTaskCount * 12);
			}
		}

		int secondTaskCount = 0;
		for (int i = 0; i < this.mgfFiles.length; i++) {

			File specFilei = new File(this.mgfFiles[i]);
			String name = specFilei.getName();
			name = name.substring(0, name.lastIndexOf("."));

			File firstFasta = new File(firstSearch, name + ".fasta");
			File firstDir = new File(firstSearch, name);
			if (!firstDir.exists()) {
				firstDir.mkdir();
			}

			File pSpecResult = new File(firstDir, "pFind.spectra");
			if (pSpecResult.exists() && pSpecResult.length() > 0) {
				logger.info(taskName + ": first search of " + name + " in original database finished");
				System.out.println(format.format(new Date()) + "\t" + taskName
						+ ": searching clustered spectra in original database finished...");
			} else {
				logger.error(taskName + ": error in searching " + this.mgfFiles[i] + " in first search, task failed");
				System.err.println(format.format(new Date()) + "\t" + taskName + ": error in searching "
						+ this.mgfFiles[i] + " in first search, task failed");
			}

			int firstProCount = 0;

			File dbReducerFasta = new File(firstDir, "DBReducer.fasta");
			if (dbReducerFasta.exists()) {
				try (BufferedReader reduceDbReader = new BufferedReader(new FileReader(dbReducerFasta))) {
					PrintWriter dbWriter = new PrintWriter(firstFasta);
					String line = null;
					while ((line = reduceDbReader.readLine()) != null) {
						dbWriter.println(line);
						if (line.startsWith(">")) {
							firstProCount++;
						}
					}
					reduceDbReader.close();
					dbWriter.close();
				} catch (IOException e) {
					logger.error(taskName + ": error in reading the database file " + dbReducerFasta, e);
					System.err.println(format.format(new Date()) + "\t" + taskName
							+ ": error in reading the database file " + dbReducerFasta);
				}
			} else {
				firstProCount = this.writeDatabase(originalFasta, firstFasta.getAbsolutePath(), false, pSpecResult);
			}

			if (!firstFasta.exists() || firstFasta.length() == 0 || firstProCount == 0) {
				logger.error(taskName + ": error in creating sample-specific database based on " + this.mgfFiles[i]
						+ " in first search");
				System.err.println(format.format(new Date()) + "\t" + taskName
						+ ": error in creating sample-specific database based on " + this.mgfFiles[i]
						+ " in first search");
				continue;

			} else {
				logger.info(taskName + ": first search of " + name + " in original database finished, " + firstProCount
						+ " proteins added to " + firstFasta);
				System.out.println(format.format(new Date()) + "\t" + taskName + ": first search of " + name
						+ " in original database finished, " + firstProCount + " proteins added to " + firstFasta);

				File secondDir = new File(secondSearch, name);
				if (!secondDir.exists()) {
					secondDir.mkdir();
				}

				secondResults[i] = new File(secondDir, "pFind-Filtered.spectra");

				if (secondResults[i].exists() && secondResults[i].length() > 0) {
					logger.info(taskName + ": second search results " + secondResults[i] + " already existed");
					System.out.println(format.format(new Date()) + "\t" + taskName + ": second search results "
							+ secondResults[i] + " already existed");
					continue;
				} else {
					secondTaskCount++;
				}
			}
		}

		if (secondTaskCount > 0) {

			File pFindBin = pFindExe.getParentFile();
			int totalThreadCount = metaPar.getThreadCount();
			int threadCount = totalThreadCount > secondTaskCount ? totalThreadCount / secondTaskCount : 1;
			int threadPoolCount = totalThreadCount > secondTaskCount ? secondTaskCount : totalThreadCount;

			if (this.useDBReducer) {
				DBReducerTask task = new DBReducerTask(dbReducerExe);
				task.setTotalThread(threadCount);

				for (int i = 0; i < this.mgfFiles.length; i++) {

					File specFilei = new File(this.mgfFiles[i]);
					String name = specFilei.getName();
					name = name.substring(0, name.lastIndexOf("."));

					File firstFasta = new File(firstSearch, name + ".fasta");

					File secondDir = new File(secondSearch, name);
					if (!secondDir.exists()) {
						secondDir.mkdir();
					}

					secondResults[i] = new File(secondDir, "DBReducer.fasta");

					if (secondResults[i].exists() && secondResults[i].length() > 0) {
						logger.info(taskName + ": refined database " + secondResults[i] + " already existed");
						System.out.println(format.format(new Date()) + "\t" + taskName + ": refined database "
								+ secondResults[i] + " already existed");

						continue;
					}

					logger.info(taskName + ": second search of " + name + " against the original database started");
					System.out.println(format.format(new Date()) + "\t" + taskName + ": second search of " + name
							+ " against the original database started...");

					task.searchMgf(new String[] { this.mgfFiles[i] }, secondDir.getAbsolutePath(),
							firstFasta.getAbsolutePath(), ms2Type, 0.01, metaPar.getFixMods(), metaPar.getVariMods());
				}

				task.run(threadPoolCount, secondTaskCount * 6);

			} else {
				PFindTask task = new PFindTask(pFindBin);
				task.setTotalThread(1);

				for (int i = 0; i < this.mgfFiles.length; i++) {

					File mgfFilei = new File(this.mgfFiles[i]);
					String name = mgfFilei.getName();
					name = name.substring(0, name.lastIndexOf("."));

					File firstFasta = new File(firstSearch, name + ".fasta");

					File secondDir = new File(secondSearch, name);
					if (!secondDir.exists()) {
						secondDir.mkdir();
					}

					secondResults[i] = new File(secondDir, "pFind-Filtered.spectra");

					if (secondResults[i].exists() && secondResults[i].length() > 0) {
						logger.info(taskName + ": second search results " + secondResults[i] + " already existed");
						System.out.println(format.format(new Date()) + "\t" + taskName + ": second search results "
								+ secondResults[i] + " already existed");

						continue;
					}

					logger.info(taskName + ": second search of " + name + " against the original database started");
					System.out.println(format.format(new Date()) + "\t" + taskName + ": second search of " + name
							+ " against the original database started...");

					if ((new File(this.pf2Files[i])).exists()) {
						task.searchPF2(new String[] { this.pf2Files[i] }, secondDir.getAbsolutePath(),
								firstFasta.getAbsolutePath(), ms2Type, false, true, 0.01, metaPar.getFixMods(),
								metaPar.getVariMods(), metaPar.getEnzyme(), metaPar.getMissCleavages(),
								metaPar.getDigestMode(), metaPar.getIsobaricTag(), true);
					} else {
						task.searchMgf(new String[] { this.mgfFiles[i] }, secondDir.getAbsolutePath(),
								firstFasta.getAbsolutePath(), ms2Type, false, true, 0.01, metaPar.getFixMods(),
								metaPar.getVariMods(), metaPar.getEnzyme(), metaPar.getMissCleavages(),
								metaPar.getDigestMode(), metaPar.getIsobaricTag(), true);
					}
				}

				task.run(threadPoolCount, secondTaskCount * 6);
			}
		}

		for (int i = 0; i < this.mgfFiles.length; i++) {

			File mgfFilei = new File(this.mgfFiles[i]);
			String name = mgfFilei.getName();
			name = name.substring(0, name.lastIndexOf("."));

			File firstFasta = new File(firstSearch, name + ".fasta");

			if (secondResults[i].exists() && secondResults[i].length() > 0) {
				logger.info(taskName + ": second search of " + name + " in " + firstFasta.getName() + " finished");
				System.out.println(format.format(new Date()) + "\t" + taskName + ": second search of " + name + " in "
						+ firstFasta.getName() + " finished");
			} else {
				logger.error(taskName + ": second search of " + name + " in " + firstFasta.getName() + " failed");
				System.err.println(format.format(new Date()) + "\t" + taskName + ": second search of " + name + " in "
						+ firstFasta.getName() + " failed");
			}
		}

		File[] firstFiles = firstSearch.listFiles();
		for (int i = 0; i < firstFiles.length; i++) {
			if (firstFiles[i].isDirectory()) {
				PFindTask.cleanPFindResultFolder(firstFiles[i], true);
			}
		}

		File[] secondFiles = secondSearch.listFiles();
		for (int i = 0; i < secondFiles.length; i++) {
			if (secondFiles[i].isDirectory()) {
				PFindTask.cleanPFindResultFolder(secondFiles[i], true);
			}
		}

		logger.info(taskName + ": MetaProIQ finished");

		System.out.println(format.format(new Date()) + "\t" + taskName + ":");

		updateProgress(30, 100);

		return secondResults;
	}

	@Override
	protected int writeDatabase(boolean cluster) {
		// TODO Auto-generated method stub

		if (cluster) {
			if (this.useDBReducer) {
				int count = 0;
				try (BufferedReader reader = new BufferedReader(new FileReader(clusterResults[0]));
						PrintWriter writer = new PrintWriter(ssdbFile)) {
					String line = null;
					while ((line = reader.readLine()) != null) {
						writer.println(line);
						if (line.startsWith(">")) {
							count++;
						}
					}
					reader.close();
					writer.close();
				} catch (IOException e) {
					logger.error(taskName + ": error in reading the DBReducer result file " + clusterResults[0], e);
					System.err.println(format.format(new Date()) + "\t" + taskName
							+ ": error in reading the DBReducer result file " + clusterResults[0]);
				}
				return count;
			} else {
				return writeDatabase(originalFasta, ssdbFile.getAbsolutePath(), false, clusterResults);
			}
		} else {
			if (this.useDBReducer) {
				HashMap<String, ProteinItem> proMap = new HashMap<String, ProteinItem>();
				for (File reduceDb : clusterResults) {
					try {
						FastaReader fr = new FastaReader(reduceDb);
						ProteinItem[] pis = fr.getItems();
						for (int i = 0; i < pis.length; i++) {
							String ref = pis[i].getRef();
							proMap.put(ref, pis[i]);
						}
					} catch (IOException e) {
						logger.error(taskName + ": error in reading the DBReducer result file " + reduceDb, e);
						System.err.println(format.format(new Date()) + "\t" + taskName
								+ ": error in reading the DBReducer result file " + reduceDb);
					}
				}

				try (PrintWriter writer = new PrintWriter(ssdbFile)) {
					for (String ref : proMap.keySet()) {
						writer.println(">" + ref);
						writer.println(proMap.get(ref).getSequence());
					}
					writer.close();
				} catch (IOException e) {
					// TODO Auto-generated catch block
					logger.error(taskName + ": error in writing the protein sequences to " + ssdbFile, e);
					System.err.println(format.format(new Date()) + "\t" + taskName
							+ ": error in writing the protein sequences to " + ssdbFile);
				}

				return proMap.size();
			} else {
				return writeDatabase(originalFasta, ssdbFile.getAbsolutePath(), true, clusterResults);
			}
		}
	}

	private HashSet<String> filterRedundantPros(HashMap<String, HashSet<String>> pepProMap,
			HashMap<String, HashSet<String>> proPepMap) throws IOException {

		String[] proteins = proPepMap.keySet().toArray(new String[proPepMap.size()]);
		Arrays.sort(proteins, (g1, g2) -> {
			int s1 = proPepMap.get(g1).size();
			int s2 = proPepMap.get(g2).size();
			return s2 - s1;
		});

		int proCount = 0;
		for (int iteratorCount = 0; iteratorCount < 10; iteratorCount++) {
			HashSet<String> totalSet = new HashSet<String>();
			HashMap<String, Integer> orderMap = new HashMap<String, Integer>();
			int count = 0;
			for (int i = 0; i < proteins.length; i++) {
				int currentSize = totalSet.size();
				totalSet.addAll(proPepMap.get(proteins[i]));
				int increaseSize = totalSet.size() - currentSize;
				orderMap.put(proteins[i], increaseSize);

				if (totalSet.size() == pepProMap.size()) {
					count = i;
				}
			}

			if (proCount == count) {
				break;
			} else {
				Arrays.sort(proteins, (g1, g2) -> {
					if (orderMap.get(g2) == orderMap.get(g1)) {
						return proPepMap.get(g2).size() - proPepMap.get(g1).size();
					} else {
						return orderMap.get(g2) - orderMap.get(g1);
					}
				});

				proCount = count;

				ArrayList<String> list = new ArrayList<String>();
				for (int i = 0; i < proteins.length; i++) {
					if (orderMap.get(proteins[i]) == 0) {
						break;
					} else {
						list.add(proteins[i]);
					}
				}
				proteins = list.toArray(String[]::new);
			}
		}

		HashSet<String> proSet = new HashSet<String>();
		HashSet<String> totalPepSet = new HashSet<String>();
		for (int i = 0; i < proteins.length; i++) {
			proSet.add(proteins[i]);
			totalPepSet.addAll(proPepMap.get(proteins[i]));
			if (totalPepSet.size() == pepProMap.size()) {
				break;
			}
		}

		return proSet;
	}

	/**
	 * Time consuming
	 * 
	 * @param pepProMap
	 * @param proPepMap
	 * @return
	 * @throws IOException
	 * @deprecated
	 */
	@SuppressWarnings("unused")
	private HashSet<String> filterRedundantProsNotused(HashMap<String, HashSet<String>> pepProMap,
			HashMap<String, HashSet<String>> proPepMap) throws IOException {

		HashSet<String> uniquePepSet = new HashSet<String>();
		HashSet<String> uniqueProSet = new HashSet<String>();

		Iterator<String> pepit = pepProMap.keySet().iterator();
		while (pepit.hasNext()) {
			String pep = pepit.next();
			HashSet<String> proSet = pepProMap.get(pep);

			if (proSet.size() == 1) {
				uniquePepSet.add(pep);

				for (String pro : proSet) {
					uniqueProSet.add(pro);

					HashSet<String> pepSet = proPepMap.get(pro);
					for (String subpep : pepSet) {
						uniquePepSet.add(subpep);
					}
				}
			}
		}

		for (String uniquePep : uniquePepSet) {
			HashSet<String> proSet = pepProMap.get(uniquePep);

			pepProMap.remove(uniquePep);

			for (String pro : proSet) {
				HashSet<String> pepSet = proPepMap.get(pro);
				pepSet.remove(uniquePep);
			}
		}

		HashSet<String> usedPepSet = new HashSet<String>();

		try {
			for (String pep : pepProMap.keySet()) {
				if (!usedPepSet.contains(pep)) {

					usedPepSet.add(pep);

					PepProGroup ppgroup = new PepProGroup();

					addPep(pep, usedPepSet, pepProMap, proPepMap, ppgroup);

					ppgroup.refine(proPepMap);

					uniqueProSet.addAll(ppgroup.getProCountMap().keySet());
				}
			}
		} catch (Exception e) {

		} finally {

		}

		return uniqueProSet;
	}

	private void addPep(String pep, HashSet<String> usedPepSet, HashMap<String, HashSet<String>> pepProMap,
			HashMap<String, HashSet<String>> proPepMap, PepProGroup ppgroup) {

		HashSet<String> pros = pepProMap.get(pep);
		ppgroup.addPep(pep, pros);

		for (String pro : pros) {
			HashSet<String> peps = proPepMap.get(pro);
			for (String subpep : peps) {
				if (!ppgroup.containsPep(subpep)) {
					usedPepSet.add(subpep);
					addPep(subpep, usedPepSet, pepProMap, proPepMap, ppgroup);
				}
			}
		}
	}

	/**
	 * 
	 * @param originalDb
	 * @param ssDb
	 * @param targetDecoy
	 * @param searchResult
	 * @return
	 */
	private int writeDatabase(String originalDb, String ssDb, boolean targetDecoy, File... searchResult) {

		HashSet<String> set = new HashSet<String>();

		for (File file : searchResult) {

			PFindPsmReader reader = new PFindPsmReader(file);
			PFindPSM[] psms = reader.getPsms();

			if (targetDecoy) {
				for (int i = 0; i < psms.length; i++) {
					if (psms[i].isTarget()) {
						String[] proteins = psms[i].getProteins();
						for (int j = 0; j < proteins.length; j++) {
							set.add(proteins[j]);
						}
					}
				}
			} else {

				HashMap<String, HashSet<String>> pepProMap = new HashMap<String, HashSet<String>>();
				HashMap<String, HashSet<String>> proPepMap = new HashMap<String, HashSet<String>>();

				for (int i = 0; i < psms.length; i++) {
					if (psms[i].isTarget()) {
						String sequence = psms[i].getSequence();
						String[] pros = psms[i].getProteins();
						for (String pro : pros) {
							if (pepProMap.containsKey(sequence)) {
								pepProMap.get(sequence).add(pro);
							} else {
								HashSet<String> proSet = new HashSet<String>();
								proSet.add(pro);
								pepProMap.put(sequence, proSet);
							}

							if (proPepMap.containsKey(pro)) {
								proPepMap.get(pro).add(sequence);
							} else {
								HashSet<String> pepSet = new HashSet<String>();
								pepSet.add(sequence);
								proPepMap.put(pro, pepSet);
							}
						}
					}
				}

				if (proPepMap.size() > 1000000) {
					try {
						set = this.filterRedundantPros(pepProMap, proPepMap);
					} catch (IOException e) {
						// TODO Auto-generated catch block
						logger.error(taskName + ": error in filtering the redundant proteins", e);
						System.err.println(format.format(new Date()) + "\t" + taskName
								+ ": error in filtering the redundant proteins");
					}
				} else {
					set.addAll(proPepMap.keySet());
				}
			}
		}

		int count = 0;

		try (BufferedReader fastaReader = new BufferedReader(new FileReader(originalDb));
				PrintWriter writer = new PrintWriter(ssDb)) {
			boolean write = false;
			String line = null;
			while ((line = fastaReader.readLine()) != null) {
				if (line.startsWith(">")) {
					String ref;
					int id = line.indexOf(" ");
					if (id > 1) {
						ref = line.substring(1, id);
					} else {
						ref = line.substring(1);
					}

					if (set.contains(ref)) {
						write = true;
						writer.println(line);
						count++;
					} else {
						write = false;
					}
				} else {
					if (write) {
						writer.println(line);
					}
				}
			}
			fastaReader.close();
			writer.close();

		} catch (IOException e) {
			// TODO Auto-generated catch block
			logger.error(taskName + ": error in writing protein database to " + ssDb, e);
			System.err.println(
					format.format(new Date()) + "\t" + taskName + ": error in writing protein database to " + ssDb);
		}

		return count;
	}

}
