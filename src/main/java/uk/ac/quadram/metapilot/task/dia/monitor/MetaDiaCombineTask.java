package uk.ac.quadram.metapilot.task.dia.monitor;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Random;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import javafx.scene.control.ProgressBar;
import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNTask;
import uk.ac.quadram.metapilot.dbSearch.diann.DiannParameter;
import uk.ac.quadram.metapilot.gui.JavaFXProgressBar;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters.ToolType;
import uk.ac.quadram.metapilot.task.dia.DiaLibSearchPar;
import uk.ac.quadram.metapilot.task.dia.MetaDiaTask;
import uk.ac.quadram.metapilot.task.dia.monitor.par.MetaParameterDiaRT;
import uk.ac.quadram.metapilot.task.dia.par.MetaParameterDia;
import uk.ac.quadram.metapilot.task.mag.par.MetaParameterMag;
import uk.ac.quadram.metapilot.task.par.MetaData;

public class MetaDiaCombineTask extends MetaDiaTask {

	private HashMap<String, String[]> rawResultFileMap;
	private static final Logger LOGGER = LogManager.getLogger(MetaDiaMonitorTask.class);
	private static final String taskName = "Realtime DIA data analysis";

	public MetaDiaCombineTask(MetaParameterDiaRT metaPar, MetaResourceParameters mrp, ProgressBar bar1,
			ProgressBar bar2, HashMap<String, String[]> rawResultFileMap) {
		super(metaPar, mrp, bar1, bar2);
		this.rawResultFileMap = rawResultFileMap;
	}

	protected void initial() {

		this.magDb = ((MetaParameterMag) metaPar).getUsedMagDbItem();
		this.threadCount = metaPar.getThreadCount();

		this.diaNNTask = new DiaNNTask(mrp.getToolPath(ToolType.DIANN));
		this.diaLibSearchPar = new DiaLibSearchPar();
		this.diaLibSearchPar.setThreadCount(metaPar.getThreadCount());

		this.diannPar = new DiannParameter();
		this.diannPar.setThreads(metaPar.getThreadCount());

		setTaskCount();
	}

	@Override
	protected Boolean call() throws IOException {

		boolean iden = iden();
		if (!iden) {
			LOGGER.error(getTaskName() + ": error in peptide and protein identification");
			System.err.println(
					format.format(new Date()) + "\t" + getTaskName() + ": error in peptide and protein identification");
			return false;
		}

		updateProgress(90, 100);

		exportReport();

		cleanTempFiles();

		updateProgress(100, 100);

		return true;
	}

	protected boolean iden() {

		updateMessage(getTaskName() + ": peptide and protein identification");
		LOGGER.info(getTaskName() + ": peptide and protein identification started");
		System.out.println(
				format.format(new Date()) + "\t" + getTaskName() + ": peptide and protein identification started");

		this.resultFolderFile = new File(metaPar.getResult(), "combined");
		if (!this.resultFolderFile.exists()) {
			resultFolderFile.mkdirs();
		}
		this.expNameMap = new HashMap<String, String>();

		this.rawFiles = new String[rawResultFileMap.size()];
		this.expNames = new String[rawResultFileMap.size()];
		this.quanExpNames = new String[rawResultFileMap.size()];
		this.separateResultFolders = new File[rawResultFileMap.size()];

		int fileid = 0;
		for (String key : rawResultFileMap.keySet()) {
			String[] files = rawResultFileMap.get(key);
			expNames[fileid] = key;
			quanExpNames[fileid] = key;
			rawFiles[fileid] = files[0];
			separateResultFolders[fileid] = new File(files[1]);
			fileid++;
		}
		this.metadata = new MetaData(this.rawFiles, this.expNames);
		metaPar.setMetadata(metadata);
		
		this.diaFiles = rawFiles;

		File magDbFile = this.magDb.getCurrentFile();
		File libFolderFile = new File(magDbFile, "libraries");
		String[] usedLibs = ((MetaParameterDia) metaPar).getLibrary();
		File hostTsvFile = null;
		for (int i = 0; i < usedLibs.length; i++) {
			if (usedLibs[i].equals("host")) {

				File hostFile = new File(metaPar.getResult(), "host");
				if (!hostFile.exists()) {
					hostFile.mkdir();
				}

				hostTsvFile = new File(hostFile, "host.tsv");
				if (!hostTsvFile.exists()) {
					File libFile = new File(new File(libFolderFile, usedLibs[i]), usedLibs[i] + ".speclib");

					diaLibSearchPar.setMbr(false);
					diaLibSearchPar.setDiaFiles(diaFiles);
					diaLibSearchPar.setThreadCount(threadCount);

					diannPar.setThreads(threadCount);
					diannPar.setMatrices(false);
					diannPar.setRelaxed_prot_inf(true);
					diannPar.setReanalyse(false);

					diaNNTask.addTask(diannPar, diaLibSearchPar, hostTsvFile.getAbsolutePath(),
							libFile.getAbsolutePath());
					diaNNTask.run(diaFiles.length * 12);
				}

				break;
			}
		}

		diaLibSearchPar.setMbr(true);
		diaLibSearchPar.setDiaFiles(diaFiles);
		diaLibSearchPar.setThreadCount(threadCount);

		diannPar.setThreads(threadCount);
		diannPar.setMatrices(true);
		diannPar.setRelaxed_prot_inf(true);
		diannPar.setReanalyse(true);

		this.quan_pep_file = new File(resultFolderFile, "combined.pr_matrix.tsv");
		File combineFile = new File(resultFolderFile, "combined.tsv");

		if (!quan_pep_file.exists() || !combineFile.exists()) {

			LOGGER.info(getTaskName() + ": library search started");
			System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": library search started");

			HashSet<String> pepSet = new HashSet<String>();
			for (int i = 0; i < separateResultFolders.length; i++) {

				File firstFile = new File(new File(separateResultFolders[i], "firstSearch"), "firstSearch.tsv");
				if (firstFile.exists()) {
					try (BufferedReader reader = new BufferedReader(new FileReader(firstFile))) {
						String pline = reader.readLine();
						String[] title = pline.split("\t");
						int seqid = -1;
						for (int j = 0; j < title.length; j++) {
							if (title[j].equals("Stripped.Sequence")) {
								seqid = j;
							}
						}
						while ((pline = reader.readLine()) != null) {
							String[] cs = pline.split("\t");
							pepSet.add(cs[seqid]);
						}
						reader.close();
					} catch (IOException e) {
						LOGGER.error(getTaskName() + ": error in reading the library search result file " + firstFile,
								e);
						System.err.println(format.format(new Date()) + "\t" + getTaskName()
								+ ": error in reading the library search result file " + firstFile);

						return false;
					}
				}

				File secondFile = new File(new File(separateResultFolders[i], "secondSearch"), "secondSearch.tsv");
				if (secondFile.exists()) {
					try (BufferedReader reader = new BufferedReader(new FileReader(secondFile))) {
						String pline = reader.readLine();
						String[] title = pline.split("\t");
						int seqid = -1;
						for (int j = 0; j < title.length; j++) {
							if (title[j].equals("Stripped.Sequence")) {
								seqid = j;
							}
						}
						while ((pline = reader.readLine()) != null) {
							String[] cs = pline.split("\t");
							pepSet.add(cs[seqid]);
						}
						reader.close();
					} catch (IOException e) {
						LOGGER.error(getTaskName() + ": error in reading the library search result file " + secondFile,
								e);
						System.err.println(format.format(new Date()) + "\t" + getTaskName()
								+ ": error in reading the library search result file " + secondFile);

						return false;
					}
				}

				File thirdFile = new File(new File(separateResultFolders[i], "thirdSearch"), "thirdSearch.tsv");
				if (thirdFile.exists()) {
					try (BufferedReader reader = new BufferedReader(new FileReader(thirdFile))) {
						String pline = reader.readLine();
						String[] title = pline.split("\t");
						int seqid = -1;
						for (int j = 0; j < title.length; j++) {
							if (title[j].equals("Stripped.Sequence")) {
								seqid = j;
							}
						}
						while ((pline = reader.readLine()) != null) {
							String[] cs = pline.split("\t");
							pepSet.add(cs[seqid]);
						}
						reader.close();
					} catch (IOException e) {
						LOGGER.error(getTaskName() + ": error in reading the library search result file " + thirdFile,
								e);
						System.err.println(format.format(new Date()) + "\t" + getTaskName()
								+ ": error in reading the library search result file " + thirdFile);

						return false;
					}
				}

				File lastFile = new File(separateResultFolders[i], separateResultFolders[i].getName() + ".tsv");
				if (!lastFile.exists()) {
					lastFile = new File(separateResultFolders[i], "lastSearch.tsv");
				}
				if (lastFile.exists()) {
					try (BufferedReader reader = new BufferedReader(new FileReader(lastFile))) {
						String pline = reader.readLine();
						String[] title = pline.split("\t");
						int seqid = -1;
						for (int j = 0; j < title.length; j++) {
							if (title[j].equals("Stripped.Sequence")) {
								seqid = j;
							}
						}
						while ((pline = reader.readLine()) != null) {
							String[] cs = pline.split("\t");
							pepSet.add(cs[seqid]);
						}
						reader.close();
					} catch (IOException e) {
						LOGGER.error(getTaskName() + ": error in reading the library search result file " + lastFile,
								e);
						System.err.println(format.format(new Date()) + "\t" + getTaskName()
								+ ": error in reading the library search result file " + lastFile);

						return false;
					}
				}
			}

			updateProgress(10, 100);
			
			if (hostTsvFile != null && hostTsvFile.exists()) {
				try (BufferedReader reader = new BufferedReader(new FileReader(hostTsvFile))) {
					String pline = reader.readLine();
					String[] title = pline.split("\t");
					int seqid = -1;
					for (int j = 0; j < title.length; j++) {
						if (title[j].equals("Stripped.Sequence")) {
							seqid = j;
						}
					}
					while ((pline = reader.readLine()) != null) {
						String[] cs = pline.split("\t");
						pepSet.add(cs[seqid]);
					}
					reader.close();
				} catch (IOException e) {
					LOGGER.error(getTaskName() + ": error in reading the library search result file " + hostTsvFile, e);
					System.err.println(format.format(new Date()) + "\t" + getTaskName()
							+ ": error in reading the library search result file " + hostTsvFile);

					return false;
				}
			}

			File fastaFile = new File(resultFolderFile, "combined.fasta");
			if (!fastaFile.exists() || fastaFile.length() == 0) {
				try (PrintWriter fastaWriter = new PrintWriter(fastaFile)) {

					int id = 0;
					for (String pep : pepSet) {

						fastaWriter.println(">pro" + id + " pro" + id);
						fastaWriter.println(pep);
						id++;
					}
					fastaWriter.close();
				} catch (IOException e) {
					LOGGER.error(getTaskName() + ": error in writing the protein database file " + fastaFile, e);
					System.err.println(format.format(new Date()) + "\t" + getTaskName()
							+ ": error in writing the protein database file " + fastaFile);

					return false;
				}
			}

			File combineLibFile = new File(resultFolderFile, "combined.speclib");

			if (!combineLibFile.exists() || combineLibFile.length() == 0) {

				File tempLibFile = new File(resultFolderFile, "combined.predicted.speclib");
				if (tempLibFile.exists()) {
					tempLibFile.renameTo(combineLibFile);
				} else {
					diaNNTask.addTask(diannPar, fastaFile.getAbsolutePath(), combineLibFile.getAbsolutePath(), false);
					diaNNTask.run(2);

					if (tempLibFile.exists()) {
						tempLibFile.renameTo(combineLibFile);
					} else {

						LOGGER.error(getTaskName() + ": error in writing the combined library file " + combineLibFile);
						System.err.println(format.format(new Date()) + "\t" + getTaskName()
								+ ": error in writing the combined library file " + combineLibFile);

						return false;
					}
				}

			} else {

				LOGGER.info(getTaskName() + ": find the combined library file in " + combineLibFile
						+ ", go to the next step");
				System.out.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": find the combined library file in " + combineLibFile + ", go to the next step");
			}

			diaNNTask.addTask(diannPar, diaFiles, combineFile.getAbsolutePath(), combineLibFile.getAbsolutePath());
			diaNNTask.run(diaFiles.length * 4);
			
			updateProgress(50, 100);

			File shufFastaFile = new File(resultFolderFile, "shuffled.fasta");
			File shufPepFile = new File(resultFolderFile, "shuffledPeptide.tsv");

			if (!shufFastaFile.exists() || shufFastaFile.length() == 0) {
				Random random = new Random();

				try (PrintWriter shufFastaWriter = new PrintWriter(shufFastaFile)) {

					PrintWriter shufPepWriter = new PrintWriter(shufPepFile);
					shufPepWriter.println("Shuffled peptide\tOriginal peptide");

					int id = 0;
					for (String pep : pepSet) {
						List<Character> sequenceChars = new ArrayList<Character>();
						for (int i = 0; i < pep.length() - 1; i++) {
							sequenceChars.add(pep.charAt(i));
						}
						Collections.shuffle(sequenceChars, random);

						StringBuilder shuffledSequence = new StringBuilder(sequenceChars.size() + 1);
						for (Character c : sequenceChars) {
							shuffledSequence.append(c);
						}
						shuffledSequence.append(pep.charAt(pep.length() - 1));
						String shufPep = shuffledSequence.toString();

						shufFastaWriter.println(">pro" + id + " pro" + id);
						shufFastaWriter.println(pep);

						shufFastaWriter.println(">shuf_pro" + id + " shuf_pro" + id);
						shufFastaWriter.println(shufPep);

						shufPepWriter.println(shufPep + "\t" + pep);
						id++;
					}
					shufFastaWriter.close();
					shufPepWriter.close();

				} catch (IOException e) {
					LOGGER.error(getTaskName() + ": error in writing the shuffled database file " + shufFastaFile, e);
					System.err.println(format.format(new Date()) + "\t" + getTaskName()
							+ ": error in writing the shuffled database file " + shufFastaFile);

					return false;
				}
			} else {
				LOGGER.info(getTaskName() + ": find the shuffled database file in " + shufFastaFile
						+ ", go to the next step");
				System.out.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": find the shuffled database file in " + shufFastaFile + ", go to the next step");
			}

			File shufLibFile = new File(resultFolderFile, "shuffled.speclib");

			if (!shufLibFile.exists() || shufLibFile.length() == 0) {

				File tempLibFile = new File(resultFolderFile, "shuffled.predicted.speclib");
				if (tempLibFile.exists()) {
					tempLibFile.renameTo(shufLibFile);
				} else {
					diaNNTask.addTask(diannPar, shufFastaFile.getAbsolutePath(), shufLibFile.getAbsolutePath(), false);
					diaNNTask.run(2);

					if (tempLibFile.exists()) {
						tempLibFile.renameTo(shufLibFile);
					} else {

						LOGGER.error(getTaskName() + ": error in writing the shuffled library file " + shufLibFile);
						System.err.println(format.format(new Date()) + "\t" + getTaskName()
								+ ": error in writing the shuffled library file " + shufLibFile);

						return false;
					}
				}

			} else {

				LOGGER.info(
						getTaskName() + ": find the shuffled library file in " + shufLibFile + ", go to the next step");
				System.out.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": find the shuffled library file in " + shufLibFile + ", go to the next step");
			}

			File shufPsmFile = new File(resultFolderFile, "shuffled.tsv");
			if (!shufPsmFile.exists() || shufPsmFile.length() == 0) {
				diaNNTask.addTask(diannPar, diaFiles, shufPsmFile.getAbsolutePath(), shufLibFile.getAbsolutePath());
				diaNNTask.run(diaFiles.length * 4);
			}

			LOGGER.info(getTaskName() + ": library search finished");
			System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": library search finished");

		} else {

			updateProgress(50, 100);
			
			HashSet<String> pepSet = new HashSet<String>();
			try (BufferedReader reader = new BufferedReader(new FileReader(combineFile))) {
				String line = reader.readLine();
				String[] title = line.split("\t");
				int seqid = -1;
				for (int i = 0; i < title.length; i++) {
					if (title[i].equals("Stripped.Sequence")) {
						seqid = i;
					}
				}

				while ((line = reader.readLine()) != null) {
					String[] cs = line.split("\t");
					pepSet.add(cs[seqid]);
				}
				reader.close();

			} catch (IOException e) {
				LOGGER.error(getTaskName() + ": error in reading the combined result file " + combineFile, e);
				System.err.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": error in reading the combined result file " + combineFile);

				return false;
			}

			File shufFastaFile = new File(resultFolderFile, "shuffled.fasta");
			File shufPepFile = new File(resultFolderFile, "shuffledPeptide.tsv");

			if (!shufFastaFile.exists() || shufFastaFile.length() == 0) {
				Random random = new Random();

				try (PrintWriter shufFastaWriter = new PrintWriter(shufFastaFile)) {

					PrintWriter shufPepWriter = new PrintWriter(shufPepFile);
					shufPepWriter.println("Shuffled peptide\tOriginal peptide");

					int id = 0;
					for (String pep : pepSet) {
						List<Character> sequenceChars = new ArrayList<Character>();
						for (int i = 0; i < pep.length() - 1; i++) {
							sequenceChars.add(pep.charAt(i));
						}
						Collections.shuffle(sequenceChars, random);

						StringBuilder shuffledSequence = new StringBuilder(sequenceChars.size() + 1);
						for (Character c : sequenceChars) {
							shuffledSequence.append(c);
						}
						shuffledSequence.append(pep.charAt(pep.length() - 1));
						String shufPep = shuffledSequence.toString();

						shufFastaWriter.println(">pro" + id + " pro" + id);
						shufFastaWriter.println(pep);

						shufFastaWriter.println(">shuf_pro" + id + " shuf_pro" + id);
						shufFastaWriter.println(shufPep);

						shufPepWriter.println(shufPep + "\t" + pep);
						id++;
					}
					shufFastaWriter.close();
					shufPepWriter.close();

				} catch (IOException e) {
					LOGGER.error(getTaskName() + ": error in writing the shuffled database file " + shufFastaFile, e);
					System.err.println(format.format(new Date()) + "\t" + getTaskName()
							+ ": error in writing the shuffled database file " + shufFastaFile);

					return false;
				}
			} else {
				LOGGER.info(getTaskName() + ": find the shuffled database file in " + shufFastaFile
						+ ", go to the next step");
				System.out.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": find the shuffled database file in " + shufFastaFile + ", go to the next step");
			}
			
			File shufPsmFile = new File(resultFolderFile, "shuffled.tsv");
			if (!shufPsmFile.exists() || shufPsmFile.length() == 0) {
				File shufLibFile = new File(resultFolderFile, "shuffled.speclib");
				if (!shufLibFile.exists() || shufLibFile.length() == 0) {

					File tempLibFile = new File(resultFolderFile, "shuffled.predicted.speclib");
					if (tempLibFile.exists()) {
						tempLibFile.renameTo(shufLibFile);
					} else {
						diaNNTask.addTask(diannPar, shufFastaFile.getAbsolutePath(), shufLibFile.getAbsolutePath(),
								false);
						diaNNTask.run(2);

						if (tempLibFile.exists()) {
							tempLibFile.renameTo(shufLibFile);
						} else {

							LOGGER.error(getTaskName() + ": error in writing the shuffled library file " + shufLibFile);
							System.err.println(format.format(new Date()) + "\t" + getTaskName()
									+ ": error in writing the shuffled library file " + shufLibFile);

							return false;
						}
					}

				} else {

					LOGGER.info(getTaskName() + ": find the shuffled library file in " + shufLibFile
							+ ", go to the next step");
					System.out.println(format.format(new Date()) + "\t" + getTaskName()
							+ ": find the shuffled library file in " + shufLibFile + ", go to the next step");
				}

				diaNNTask.addTask(diannPar, diaFiles, shufPsmFile.getAbsolutePath(), shufLibFile.getAbsolutePath());
				diaNNTask.run(diaFiles.length * 4);
			}
		}

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
