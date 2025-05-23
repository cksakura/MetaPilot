package uk.ac.quadram.metapilot.task.pfind;

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

import uk.ac.quadram.metapilot.core.prodb.FastaReader;
import uk.ac.quadram.metapilot.core.prodb.ProteinItem;
import uk.ac.quadram.metapilot.dbSearch.PepProGroup;
import uk.ac.quadram.metapilot.dbSearch.dbReducer.DBReducerTask;
import uk.ac.quadram.metapilot.dbSearch.pfind.PFindPSM;
import uk.ac.quadram.metapilot.dbSearch.pfind.PFindPsmReader;
import uk.ac.quadram.metapilot.dbSearch.pfind.PFindTask;
import uk.ac.quadram.metapilot.gui.JavaFXProgressBar;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters.ToolType;
import uk.ac.quadram.metapilot.spectra.cluster.ClusterTask;
import uk.ac.quadram.metapilot.task.MetaDbCreateTask;
import uk.ac.quadram.metapilot.task.MetaLabWorkflowType;
import uk.ac.quadram.metapilot.task.par.MetaConstants;
import uk.ac.quadram.metapilot.task.pfind.par.MetaParameterPFind;

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

	public MetaDbReducerTask(MetaParameterPFind metaPar, MetaResourceParameters mrp, JavaFXProgressBar progressBar) {
		super(metaPar, mrp, progressBar);
		this.dbReducerExe = new File(mrp.getToolPath(ToolType.DBREDUCER));
		this.pFindExe = new File(mrp.getToolPath(ToolType.PFIND));
		if (dbReducerExe.getName().equals("DBReducer.exe") && dbReducerExe.exists()
				&& metaPar.getMs2ScanMode().equals(MetaConstants.HCD_FTMS)) {
			useDBReducer = true;
		}
	}

	@Override
	protected void initial() {
		// TODO Auto-generated method stub
		this.mscp = metaPar.getMscp();
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
		File spectraDirFile = metaPar.getSpectraFile();
		for (int i = 0; i < rawFiles.length; i++) {
			File rawi = new File(rawFiles[i]);
			String rawName = rawi.getName();
			rawName = rawName.substring(0, rawName.length() - ".raw".length());
			File mgfFile = new File(spectraDirFile, rawName + suffix + ".mgf");
			File pf2 = new File(spectraDirFile, rawName + suffix + ".pf2");
			mgfFiles[i] = mgfFile.getAbsolutePath();
			pf2Files[i] = pf2.getAbsolutePath();
		}
	}

	/**
	 * Clustering strategy
	 * 
	 * @throws Exception
	 */
	protected void cluster() {

		progressBar.setString(taskName + ": spectra clustering...");

		LOGGER.info(taskName + ": spectra cluster started");
		System.out.println(format.format(new Date()) + "\t" + taskName + ": spectra cluster started");

		File clusterDir = new File(dbcFile, "clustered_spectra");
		File clusterMgf = new File(clusterDir, "cluster.mgf");
		clusterDir.mkdir();

		if (clusterMgf.exists()) {
			LOGGER.info(taskName + ": clustered MS/MS spectra file already existed in " + clusterMgf);
			System.out.println(format.format(new Date()) + "\t" + taskName
					+ ": clustered MS/MS spectra file already existed in " + clusterMgf);
		} else {

			LOGGER.info(taskName + ": clustering MS/MS spectra to " + clusterMgf + " started");
			System.out.println(format.format(new Date()) + "\t" + taskName + ": clustering MS/MS spectra to "
					+ clusterMgf + " started");

			ClusterTask clusterTask = new ClusterTask(metaPar.getThreadCount(), true,
					((MetaParameterPFind) metaPar).getIsobaricTag(), MetaLabWorkflowType.pFindWorkflow.name());

			clusterTask.runCli(mgfFiles, clusterDir, clusterMgf);

			LOGGER.info(taskName + ": clustering MS/MS spectra to " + clusterMgf + " finished");
			System.out.println(format.format(new Date()) + "\t" + taskName + ": clustering MS/MS spectra to "
					+ clusterMgf + " finished");
		}

		LOGGER.info(taskName + ": spectra cluster finished");
		System.out.println(format.format(new Date()) + "\t" + taskName + ": spectra cluster finished");

		setProgress(30);
	}

	/**
	 * Clustering strategy
	 * 
	 * 
	 */
	protected File[] clusterSearch() {

		progressBar.setString(taskName + ": searching clustered spectra against the original database");

		LOGGER.info(taskName + ": searching clustered spectra against the original database started");
		System.out.println(format.format(new Date()) + "\t" + taskName
				+ ": searching clustered spectra in original database started");

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

		if (useDBReducer) {
			File reducedDbFile = new File(dbcFile, "DBReducer.fasta");
			if (reducedDbFile.exists() && reducedDbFile.length() > 0) {

				LOGGER.info(taskName + ": refined database already existed in " + reducedDbFile);
				System.out.println(format.format(new Date()) + "\t" + taskName
						+ ": refined database already existed in " + reducedDbFile);

			} else {
				if (reducedDbFile.exists() && reducedDbFile.length() == 0) {
					FileUtils.deleteQuietly(reducedDbFile);
				}

				LOGGER.info(taskName + ": generating refined database from the original database started");
				System.out.println(format.format(new Date()) + "\t" + taskName
						+ ": generating refined database from the original database started");

				String ms2Type = metaPar.getMs2ScanMode();

				String[] mgfFilePaths = new String[mgfs.length];
				for (int i = 0; i < mgfFilePaths.length; i++) {
					mgfFilePaths[i] = mgfs[i].getAbsolutePath();
				}

				DBReducerTask task = new DBReducerTask(dbReducerExe);
				task.setTotalThread(metaPar.getThreadCount());

				task.searchMgf(mgfFilePaths, dbcFile.getAbsolutePath(), originalFasta, ms2Type, psmFdr,
						((MetaParameterPFind) this.metaPar).getFixMods(),
						((MetaParameterPFind) this.metaPar).getVariMods());
				task.run(1, mgfFilePaths.length * 48);
			}

			if (reducedDbFile.exists() && reducedDbFile.length() > 0) {
				LOGGER.info(taskName + ": generating refined database from the original database finished");
				System.out.println(format.format(new Date()) + "\t" + taskName
						+ ": generating refined database from the original database finished");
			}

			PFindTask.cleanPFindResultFolder(dbcFile, true);

			setProgress(35);

			return new File[] { reducedDbFile };

		} else {
			File pSpecResult = new File(dbcFile, "pFind.spectra");

			if (pSpecResult.exists() && pSpecResult.length() > 0) {

				LOGGER.info(taskName + ": search results already existed in " + pSpecResult);
				System.out.println(format.format(new Date()) + "\t" + taskName + ": search results already existed in "
						+ pSpecResult);

			} else {

				if (pSpecResult.exists() && pSpecResult.length() == 0) {
					FileUtils.deleteQuietly(pSpecResult);
				}

				LOGGER.info(taskName + ": searching clustered spectra in original database started");
				System.out.println(format.format(new Date()) + "\t" + taskName
						+ ": searching clustered spectra in original database started");

				String ms2Type = metaPar.getMs2ScanMode();

				String[] mgfFilePaths = new String[mgfs.length];
				for (int i = 0; i < mgfFilePaths.length; i++) {
					mgfFilePaths[i] = mgfs[i].getAbsolutePath();
				}

				File pFindBin = pFindExe.getParentFile();

				PFindTask task = new PFindTask(pFindBin);
				task.setTotalThread(metaPar.getThreadCount());

				task.searchMgf(mgfFilePaths, dbcFile.getAbsolutePath(), originalFasta, ms2Type, mscp.isOpen(), false,
						psmFdr, ((MetaParameterPFind) this.metaPar).getFixMods(),
						((MetaParameterPFind) this.metaPar).getVariMods(),
						((MetaParameterPFind) this.metaPar).getEnzyme(),
						((MetaParameterPFind) this.metaPar).getMissCleavages(),
						((MetaParameterPFind) this.metaPar).getDigestMode(),
						((MetaParameterPFind) this.metaPar).getIsobaricTag(), true);
				task.run(1, mgfFilePaths.length * 48);
			}

			if (pSpecResult.exists() && pSpecResult.length() > 0) {
				LOGGER.info(taskName + ": searching clustered spectra against the original database finished");
				System.out.println(format.format(new Date()) + "\t" + taskName
						+ ": searching clustered spectra against the original database finished");
			}

			PFindTask.cleanPFindResultFolder(pSpecResult.getParentFile(), true);

			setProgress(35);

			return new File[] { pSpecResult };
		}
	}

	@Override
	protected File[] metaProIQ() {
		// TODO Auto-generated method stub

		progressBar.setString(taskName + ": searching clustered spectra in original database");

		LOGGER.info(taskName + ": MetaProIQ started");
		System.out.println(format.format(new Date()) + "\t" + taskName + ": MetaProIQ started");

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

							LOGGER.info(taskName + ": first search of " + name + " in original database started");
							System.out.println(format.format(new Date()) + "\t" + taskName + ": first search of " + name
									+ " in original database started...");

							task.searchMgf(new String[] { this.mgfFiles[i] }, firstDir.getAbsolutePath(), originalFasta,
									ms2Type, 0.01, ((MetaParameterPFind) this.metaPar).getFixMods(),
									((MetaParameterPFind) this.metaPar).getVariMods());
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

							LOGGER.info(taskName + ": first search of " + name + " in original database started");
							System.out.println(format.format(new Date()) + "\t" + taskName + ": first search of " + name
									+ " in original database started...");

							if ((new File(this.pf2Files[i])).exists()) {
								task.searchPF2(new String[] { this.pf2Files[i] }, firstDir.getAbsolutePath(),
										originalFasta, ms2Type, mscp.isOpen(), false, 0.01,
										((MetaParameterPFind) this.metaPar).getFixMods(),
										((MetaParameterPFind) this.metaPar).getVariMods(),
										((MetaParameterPFind) this.metaPar).getEnzyme(),
										((MetaParameterPFind) this.metaPar).getMissCleavages(),
										((MetaParameterPFind) this.metaPar).getDigestMode(),
										((MetaParameterPFind) this.metaPar).getIsobaricTag(), true);
							} else {

								task.searchMgf(new String[] { this.mgfFiles[i] }, firstDir.getAbsolutePath(),
										originalFasta, ms2Type, mscp.isOpen(), false, 0.01,
										((MetaParameterPFind) this.metaPar).getFixMods(),
										((MetaParameterPFind) this.metaPar).getVariMods(),
										((MetaParameterPFind) this.metaPar).getEnzyme(),
										((MetaParameterPFind) this.metaPar).getMissCleavages(),
										((MetaParameterPFind) this.metaPar).getDigestMode(),
										((MetaParameterPFind) this.metaPar).getIsobaricTag(), true);
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
				LOGGER.info(taskName + ": first search of " + name + " in original database finished");
				System.out.println(format.format(new Date()) + "\t" + taskName
						+ ": searching clustered spectra in original database finished...");
			} else {
				LOGGER.error(taskName + ": error in searching " + this.mgfFiles[i] + " in first search, task failed");
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
					LOGGER.error(taskName + ": error in reading the database file " + dbReducerFasta, e);
					System.err.println(format.format(new Date()) + "\t" + taskName
							+ ": error in reading the database file " + dbReducerFasta);
				}
			} else {
				firstProCount = this.writeDatabase(originalFasta, firstFasta.getAbsolutePath(), false, pSpecResult);
			}

			if (!firstFasta.exists() || firstFasta.length() == 0 || firstProCount == 0) {
				LOGGER.error(taskName + ": error in creating sample-specific database based on " + this.mgfFiles[i]
						+ " in first search");
				System.err.println(format.format(new Date()) + "\t" + taskName
						+ ": error in creating sample-specific database based on " + this.mgfFiles[i]
						+ " in first search");
				continue;

			} else {
				LOGGER.info(taskName + ": first search of " + name + " in original database finished, " + firstProCount
						+ " proteins added to " + firstFasta);
				System.out.println(format.format(new Date()) + "\t" + taskName + ": first search of " + name
						+ " in original database finished, " + firstProCount + " proteins added to " + firstFasta);

				File secondDir = new File(secondSearch, name);
				if (!secondDir.exists()) {
					secondDir.mkdir();
				}

				secondResults[i] = new File(secondDir, "pFind-Filtered.spectra");

				if (secondResults[i].exists() && secondResults[i].length() > 0) {
					LOGGER.info(taskName + ": second search results " + secondResults[i] + " already existed");
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
						LOGGER.info(taskName + ": refined database " + secondResults[i] + " already existed");
						System.out.println(format.format(new Date()) + "\t" + taskName + ": refined database "
								+ secondResults[i] + " already existed");

						continue;
					}

					LOGGER.info(taskName + ": second search of " + name + " against the original database started");
					System.out.println(format.format(new Date()) + "\t" + taskName + ": second search of " + name
							+ " against the original database started...");

					task.searchMgf(new String[] { this.mgfFiles[i] }, secondDir.getAbsolutePath(),
							firstFasta.getAbsolutePath(), ms2Type, 0.01,
							((MetaParameterPFind) this.metaPar).getFixMods(),
							((MetaParameterPFind) this.metaPar).getVariMods());
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
						LOGGER.info(taskName + ": second search results " + secondResults[i] + " already existed");
						System.out.println(format.format(new Date()) + "\t" + taskName + ": second search results "
								+ secondResults[i] + " already existed");

						continue;
					}

					LOGGER.info(taskName + ": second search of " + name + " against the original database started");
					System.out.println(format.format(new Date()) + "\t" + taskName + ": second search of " + name
							+ " against the original database started...");

					if ((new File(this.pf2Files[i])).exists()) {
						task.searchPF2(new String[] { this.pf2Files[i] }, secondDir.getAbsolutePath(),
								firstFasta.getAbsolutePath(), ms2Type, mscp.isOpen(), true, 0.01,
								((MetaParameterPFind) this.metaPar).getFixMods(),
								((MetaParameterPFind) this.metaPar).getVariMods(),
								((MetaParameterPFind) this.metaPar).getEnzyme(),
								((MetaParameterPFind) this.metaPar).getMissCleavages(),
								((MetaParameterPFind) this.metaPar).getDigestMode(),
								((MetaParameterPFind) this.metaPar).getIsobaricTag(), true);
					} else {
						task.searchMgf(new String[] { this.mgfFiles[i] }, secondDir.getAbsolutePath(),
								firstFasta.getAbsolutePath(), ms2Type, mscp.isOpen(), true, 0.01,
								((MetaParameterPFind) this.metaPar).getFixMods(),
								((MetaParameterPFind) this.metaPar).getVariMods(),
								((MetaParameterPFind) this.metaPar).getEnzyme(),
								((MetaParameterPFind) this.metaPar).getMissCleavages(),
								((MetaParameterPFind) this.metaPar).getDigestMode(),
								((MetaParameterPFind) this.metaPar).getIsobaricTag(), true);
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
				LOGGER.info(taskName + ": second search of " + name + " in " + firstFasta.getName() + " finished");
				System.out.println(format.format(new Date()) + "\t" + taskName + ": second search of " + name + " in "
						+ firstFasta.getName() + " finished");
			} else {
				LOGGER.error(taskName + ": second search of " + name + " in " + firstFasta.getName() + " failed");
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

		LOGGER.info(taskName + ": MetaProIQ finished");

		System.out.println(format.format(new Date()) + "\t" + taskName + ":");

		setProgress(30);

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
					LOGGER.error(taskName + ": error in reading the DBReducer result file " + clusterResults[0], e);
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
						LOGGER.error(taskName + ": error in reading the DBReducer result file " + reduceDb, e);
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
					LOGGER.error(taskName + ": error in writing the protein sequences to " + ssdbFile, e);
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
						LOGGER.error(taskName + ": error in filtering the redundant proteins", e);
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
			LOGGER.error(taskName + ": error in writing protein database to " + ssDb, e);
			System.err.println(
					format.format(new Date()) + "\t" + taskName + ": error in writing protein database to " + ssDb);
		}

		return count;
	}
}
