package uk.ac.quadram.metapilot.task.dia;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashSet;

import org.apache.commons.io.FileUtils;
import org.apache.logging.log4j.Logger;

import uk.ac.quadram.metapilot.core.dia.PeptideListLibrarySupport;
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
	private String lastFailureReason;
	/** fromPeptide only: what went into the library, written to _info.txt. */
	private java.util.List<String> peptideLibraryInfo;

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

		// Use buildType from UI so "Protein Database" is respected (e.g. .fa files); do not infer from extension only.
		this.taskType = buildType;
	}

	@Override
	protected Boolean call() throws Exception {
		lastFailureReason = null;

		if (magDbItem == null) {
			return fail("MAG database item is null.", null);
		}
		if (magDbItem.getCurrentFile() == null) {
			return fail("MAG database current folder is null for database: " + magDbItem.getCatalogueID(), null);
		}

		File libFolder = new File(magDbItem.getCurrentFile(), "libraries");
		if (!libFolder.exists()) {
			if (!libFolder.mkdir()) {
				return fail("Failed to create libraries folder: " + libFolder.getAbsolutePath(), null);
			}
		}

		File libFile = new File(libFolder, libName);
		if (!libFile.exists()) {
			if (!libFile.mkdir()) {
				return fail("Failed to create library folder: " + libFile.getAbsolutePath(), null);
			}
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
			if (!destLibFile.exists()) {
				File tempLibFile = new File(libFile, libName + ".predicted.speclib");
				if (tempLibFile.exists()) {
					tempLibFile.renameTo(destLibFile);
				} else {
					diaNNTask.addTask(diannPar, filePath, destLibFile.getAbsolutePath(), false);
					diaNNTask.run(2);

					if (tempLibFile.exists()) {
						tempLibFile.renameTo(destLibFile);
					} else if (!destLibFile.exists()) {
						return fail(
								"DIA-NN library generation did not produce expected file: "
										+ destLibFile.getAbsolutePath(),
								null);
					}
				}
			}
		} else if (taskType == fromPeptide) {

			// Peptide entries: no protein N-terminal rules, missed cleavages that keep the peptides intact, sorted FASTA
			// (see PeptideListLibrarySupport).
			File pepDbFile = new File(libFile, libName + ".fasta");
			PeptideListLibrarySupport.PeptideList peptideList;
			try {
				if (!pepDbFile.exists()) {
					peptideList = PeptideListLibrarySupport.readPeptideTable(new File(filePath));
					PeptideListLibrarySupport.writePeptideFasta(peptideList, pepDbFile);
				} else {
					peptideList = PeptideListLibrarySupport.readPeptideFasta(pepDbFile);
				}
			} catch (IOException e) {
				return fail("Error reading the peptide list or writing the peptide FASTA: " + filePath, e);
			}
			miss = PeptideListLibrarySupport.configureForPeptideLibrary(diannPar, miss, peptideList);
			peptideLibraryInfo = PeptideListLibrarySupport.describe(peptideList, diannPar);
			for (String line : peptideLibraryInfo) {
				logger.info(taskName + ": " + line.replace('\t', ' '));
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
						return fail(
								"DIA-NN library generation did not produce expected file: "
										+ tempLibFile.getAbsolutePath(),
								null);
					}
				}
			}
		} else {
			return fail("Unsupported library build type: " + taskType, null);
		}

		if (!parquetFile.exists()) {
			diaNNTask.addTask(diannPar, destLibFile.getAbsolutePath(), parquetFile.getAbsolutePath());
			diaNNTask.run(1);
		}

		if (!parquetFile.exists()) {
			return fail("DIA-NN parquet conversion failed. Missing output: " + parquetFile.getAbsolutePath(), null);
		}

		HashSet<String> pepSet = DiaNNParquetReader.getPeptideSet(parquetFile.getAbsolutePath());
		File infoFile = new File(libFile, libName + "_info.txt");
		try (PrintWriter writer = new PrintWriter(infoFile)) {
			writer.println("Data set size:\t" + pepSet.size());
			writer.println("Enzyme:\t" + enzyme);
			writer.println("Miss cleavages:\t" + miss);
			writer.println("MatchToMag:\t" + matchToMag);
			if (peptideLibraryInfo != null) {
				for (String line : peptideLibraryInfo) {
					writer.println(line);
				}
			}
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

	public String getLastFailureReason() {
		return lastFailureReason;
	}

	private Boolean fail(String message, Exception e) {
		lastFailureReason = message;
		if (e != null) {
			logger.error(taskName + ": " + message, e);
		} else {
			logger.error(taskName + ": " + message);
		}
		System.err.println(format.format(new Date()) + "\t" + taskName + ": " + message);
		return false;
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
