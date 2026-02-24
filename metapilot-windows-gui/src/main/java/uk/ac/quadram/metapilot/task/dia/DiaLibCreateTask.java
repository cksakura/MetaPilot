package uk.ac.quadram.metapilot.task.dia;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.apache.commons.io.FileUtils;
import org.apache.logging.log4j.Logger;

import uk.ac.quadram.metapilot.dbSearch.deepDetect.DeepDetectParameter;
import uk.ac.quadram.metapilot.dbSearch.deepDetect.DeepDetectTask;
import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNParquetReader;
import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNTask;
import uk.ac.quadram.metapilot.dbSearch.diann.DiannParameter;
import javafx.scene.control.ProgressBar;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters.ToolType;
import uk.ac.quadram.metapilot.task.MetaAbstractTask;
import uk.ac.quadram.metapilot.task.dia.par.MetaParameterDia;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;

/**
 * DIA library build task (Migrated to JavaFX Task)
 * 
 * @author Kai Cheng
 *
 */
public class DiaLibCreateTask extends MetaAbstractTask {

	private static String taskName = "MAG spectra libarary create task";
	private SimpleDateFormat format = new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss");

	private int taskType;
	private MagDbItem magDbItem;
	private String filePath;
	private String enzyme;
	private String libName;
	private int miss;

	public static final int fromLibrary = 0;
	public static final int fromProtein = 1;
	public static final int fromPeptide = 2;

	private boolean matchToMag = true;

	public DiaLibCreateTask(MetaParameterDia metaPar, MetaResourceParameters mrp, ProgressBar progressBar,
			String filePath, ArrayList<String> pathList, int buildType, Logger logger) {
		super(metaPar, mrp, progressBar, logger);
		this.magDbItem = metaPar.getMagDbItem();
		this.filePath = filePath;
		this.enzyme = pathList.get(0);
		this.libName = pathList.get(1);
		this.miss = Integer.parseInt(pathList.get(2));
		if (pathList.size() > 3) {
			this.matchToMag = Boolean.parseBoolean(pathList.get(3));
		}

		if (filePath.endsWith(".speclib")) {
			this.taskType = fromLibrary;
		} else if (filePath.endsWith(".faa") || filePath.endsWith(".fasta")) {
			this.taskType = fromProtein;
		} else if (filePath.endsWith(".tsv") || filePath.endsWith(".csv")) {
			this.taskType = fromPeptide;
		}
	}

	@Override
	protected Boolean call() throws Exception {

		File libFolder = new File(magDbItem.getCurrentFile(), "libraries");
		if (!libFolder.exists()) {
			libFolder.mkdir();
		}

		File libFile = new File(libFolder, libName);
		if (!libFile.exists()) {
			libFile.mkdir();
		}

		DiaNNTask diaNNTask = new DiaNNTask(mrp.getToolPath(ToolType.DIANN));
		DiannParameter diannPar = new DiannParameter();
		diannPar.setThreads(metaPar.getThreadCount());
		diannPar.setEnzyme(enzyme);
		diannPar.setMissed_cleavages(miss);
		if (metaPar instanceof MetaParameterDia) {
			MetaParameterDia diaPar = (MetaParameterDia) metaPar;
			diannPar.setM_excision(diaPar.isMExcision());
			diannPar.setCarbami_C(diaPar.isCarbamiC());
			diannPar.setOxidation_M(diaPar.isOxidationM());
			diannPar.setAcety_proN(diaPar.isAcetyProN());
			diannPar.setPhospho_STY(diaPar.isPhosphoSTY());
			diannPar.setK_GG(diaPar.isKGG());
		}

		File destLibFile = new File(libFile, libName + ".speclib");
		File parquetFile = new File(libFile, libName + ".parquet");

		if (taskType == fromLibrary) {
			if (!destLibFile.exists()) {
				FileUtils.copyFile(new File(filePath), destLibFile);
			}
		} else if (taskType == fromProtein) {

			File pepFile = new File(libFile, libName + ".predicted.tsv");
			if (!pepFile.exists()) {
				DeepDetectTask deepDetectTask = new DeepDetectTask(metaPar, mrp, progressBar, null);
				DeepDetectParameter deepDetectPar = new DeepDetectParameter();
				deepDetectPar.setMissCleavage(miss);
				if (enzyme.equals("Lys-C")) {
					deepDetectPar.setProtease("LysC");
				} else if (enzyme.equals("Chymotrypsin")) {
					deepDetectPar.setProtease("Chymotrypsin");
				} else if (enzyme.equals("AspN")) {
					deepDetectPar.setProtease("AspN");
				} else if (enzyme.equals("GluC")) {
					deepDetectPar.setProtease("GluC");
				} else {
					deepDetectPar.setProtease("Trypsin");
				}

				deepDetectPar.setInput(new File(filePath));
				deepDetectPar.setOutput(pepFile);

				deepDetectTask.addTask(deepDetectPar);
				deepDetectTask.run(2);
			}

			HashMap<String, Double> pepScoreMap = new HashMap<String, Double>();
			HashMap<String, HashSet<String>> pepProMap = new HashMap<String, HashSet<String>>();
			try (BufferedReader reader = new BufferedReader(new FileReader(pepFile))) {
				String pline = reader.readLine();
				while ((pline = reader.readLine()) != null) {
					String[] cs = pline.split("\t");
					double score = Double.parseDouble(cs[2]);
					if (pepScoreMap.containsKey(cs[1])) {
						if (score > pepScoreMap.get(cs[1])) {
							pepScoreMap.put(cs[1], score);
						}
					} else {
						pepScoreMap.put(cs[1], score);
					}
					if (pepProMap.containsKey(cs[1])) {
						pepProMap.get(cs[1]).add(cs[0]);
					} else {
						HashSet<String> proSet = new HashSet<String>();
						proSet.add(cs[0]);
						pepProMap.put(cs[1], proSet);
					}
				}
				reader.close();
			} catch (IOException e) {
				logger.error(taskName + ": error in reading " + pepFile, e);
				System.err.println(format.format(new Date()) + "\t" + taskName + ": error in reading " + pepFile);
			}

			File pepDbFile = new File(libFile, libName + ".fasta");
			if (!pepDbFile.exists()) {
				List<Map.Entry<String, Double>> pepList = new ArrayList<>(pepScoreMap.entrySet());
				pepList.sort(Map.Entry.comparingByValue(Collections.reverseOrder()));
				String[] sortedPeps = pepList.stream().map(Map.Entry::getKey).toArray(String[]::new);

				try (PrintWriter writer = new PrintWriter(pepDbFile)) {
					for (int i = 0; i < sortedPeps.length; i++) {
						double score = pepScoreMap.get(sortedPeps[i]);
						if (score < 0.85) {
							if (i > 2000000) {
								break;
							}
						} else {
							if (i > 3000000) {
								break;
							}
						}

						HashSet<String> proSet = pepProMap.get(sortedPeps[i]);
						if (proSet != null && proSet.size() > 0) {
							for (String pro : proSet) {
								writer.println(">" + pro);
								writer.println(sortedPeps[i]);
							}
						}
					}

					writer.close();
				} catch (IOException e) {
					logger.error(taskName + ": error in writing " + pepDbFile, e);
					System.err.println(format.format(new Date()) + "\t" + taskName + ": error in reading " + pepDbFile);
				}
			}

			if (!destLibFile.exists()) {
				File tempLibFile = new File(libFile, libName + ".predicted.speclib");
				if (tempLibFile.exists()) {
					tempLibFile.renameTo(destLibFile);
				} else {
					diaNNTask.addTask(diannPar, pepDbFile.getAbsolutePath(), destLibFile.getAbsolutePath(), false);
					diaNNTask.run(2);

					if (tempLibFile.exists()) {
						tempLibFile.renameTo(destLibFile);
					} else {
						return false;
					}
				}
			}
		} else if (taskType == fromPeptide) {

			File pepDbFile = new File(libFile, libName + ".fasta");
			if (!pepDbFile.exists()) {
				HashSet<String> pepSet = new HashSet<String>();
				try (BufferedReader reader = new BufferedReader(new FileReader(filePath))) {
					String pline = reader.readLine();
					String[] title = pline.split("\t");
					int seqid = -1;
					for (int i = 0; i < title.length; i++) {
						if (title[i].equals("Sequence")) {
							seqid = i;
						} else if (title[i].equals("Base sequence")) {
							seqid = i;
						}
					}

					if (seqid == -1) {

						logger.error(taskName + ": there is no \"Sequence\" column in the file " + filePath);
						System.out.println(format.format(new Date()) + "\t" + taskName
								+ ": there is no \"Sequence\" column in the file " + filePath);

						return false;
					}

					while ((pline = reader.readLine()) != null) {
						String[] cs = pline.split("\t");
						int c1 = 0;
						int c2 = 0;
						boolean start = true;
						StringBuilder sb = new StringBuilder();
						for (int i = 0; i < cs[seqid].length(); i++) {
							char aa = cs[seqid].charAt(i);
							if (aa == '(' || aa == '[') {
								start = false;
								c1++;
							} else if (aa == ')' || aa == ']') {
								start = true;
								c2++;
							} else {
								if (start && aa >= 'A' && aa <= 'Z') {
									sb.append(aa);
								}
							}
						}

						if (c1 == c2) {
							String seq = sb.toString();
							pepSet.add(seq);
						}
					}
					reader.close();
				} catch (IOException e) {

					logger.error(taskName + ": error in reading " + filePath, e);
					System.err.println(format.format(new Date()) + "\t" + taskName + ": error in reading " + filePath);

					return false;
				}

				String[] sortedPeps = pepSet.toArray(String[]::new);

				try (PrintWriter writer = new PrintWriter(pepDbFile)) {
					for (int i = 0; i < sortedPeps.length; i++) {
						writer.println(">pro" + i + " pro" + i);
						writer.println(sortedPeps[i]);
					}
					writer.close();
				} catch (IOException e) {
					logger.error(taskName + ": error in writing " + pepDbFile, e);
					System.err.println(format.format(new Date()) + "\t" + taskName + ": error in reading " + pepDbFile);
				}
			}

			if (!destLibFile.exists()) {
				File tempLibFile = new File(libFile, libName + ".predicted.speclib");
				if (tempLibFile.exists()) {
					tempLibFile.renameTo(destLibFile);
				} else {
					diaNNTask.addTask(diannPar, pepDbFile.getAbsolutePath(), destLibFile.getAbsolutePath(), false);
					diaNNTask.run(2);

					if (tempLibFile.exists()) {
						tempLibFile.renameTo(destLibFile);
					} else {
						return false;
					}
				}
			}
		} else {
			return false;
		}

		if (!parquetFile.exists()) {
			diaNNTask.addTask(diannPar, destLibFile.getAbsolutePath(), parquetFile.getAbsolutePath());
			diaNNTask.run(1);
		}

		if (!parquetFile.exists()) {
			return false;
		}

		HashSet<String> pepSet = DiaNNParquetReader.getPeptideSet(parquetFile.getAbsolutePath());
		File infoFile = new File(libFile, libName + "_info.txt");
		try (PrintWriter writer = new PrintWriter(infoFile)) {
			writer.println("Data set size:\t" + pepSet.size());
			writer.println("Enzyme:\t" + enzyme);
			writer.println("Miss cleavages:\t" + miss);
			writer.println("MatchToMag:\t" + matchToMag);
			writer.close();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			logger.error(taskName + ": error in writing the information to " + infoFile, e);
			System.err.println(
					format.format(new Date()) + "\t" + taskName + ": error in writing the information to " + infoFile);
		}

		magDbItem.addPepLib(libName, pepSet.size(), enzyme, miss, matchToMag);
		return true;
	}

	public void forceStop() {
		// Cancel this task
		this.cancel();
	}

	@Override
	protected String getTaskName() {
		// TODO Auto-generated method stub
		return taskName;
	}
}
