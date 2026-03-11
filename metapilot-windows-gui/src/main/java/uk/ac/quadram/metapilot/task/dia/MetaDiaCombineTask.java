package uk.ac.quadram.metapilot.task.dia;

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
import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNParquetReader;
import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNTask;
import uk.ac.quadram.metapilot.dbSearch.diann.DiannParameter;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters.ToolType;
import uk.ac.quadram.metapilot.task.dia.par.MetaParameterDia;
import uk.ac.quadram.metapilot.task.dia.par.MetaParameterDiaRT;
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

		this.magDb = ((MetaParameterDiaRT) metaPar).getMagDbItem();
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

		HashSet<String> peptideSet = new HashSet<String>();
		for (String baseName : rawResultFileMap.keySet()) {
			File perFileDir = new File(rawResultFileMap.get(baseName)[1]);
			File finalParquet = new File(perFileDir, baseName + ".parquet");
			if (!finalParquet.exists()) {
				LOGGER.warn(getTaskName() + ": missing per-file final parquet: {}", finalParquet.getAbsolutePath());
				continue;
			}
			try {
				peptideSet.addAll(DiaNNParquetReader.getPeptideSet(finalParquet.getAbsolutePath()));
			} catch (IOException e) {
				LOGGER.error(getTaskName() + ": failed reading {}", finalParquet.getAbsolutePath(), e);
				return false;
			}
		}
		if (peptideSet.isEmpty()) {
			LOGGER.error(getTaskName() + ": no peptides found from per-file final parquets");
			return false;
		}

		updateProgress(10, 100);

		File fastaFile = new File(resultFolderFile, "combined.fasta");
		if (!fastaFile.exists() || fastaFile.length() == 0) {
			try (PrintWriter fastaWriter = new PrintWriter(fastaFile)) {
				int id = 0;
				for (String pep : peptideSet) {
					fastaWriter.println(">pro" + id + " pro" + id);
					fastaWriter.println(pep);
					id++;
				}
			} catch (IOException e) {
				LOGGER.error(getTaskName() + ": error writing {}", fastaFile.getAbsolutePath(), e);
				return false;
			}
		}

		File combineLibFile = new File(resultFolderFile, "combined.speclib");
		if (!combineLibFile.exists() || combineLibFile.length() == 0) {
			File tempLibFile = new File(resultFolderFile, "combined.predicted.speclib");
			if (tempLibFile.exists()) {
				tempLibFile.renameTo(combineLibFile);
			} else {
				diannPar.setThreads(threadCount);
				diannPar.setExportDecoy(true);
				diaNNTask.addTask(diannPar, fastaFile.getAbsolutePath(), combineLibFile.getAbsolutePath(), false);
				diaNNTask.run(2);
				if (!combineLibFile.exists() && tempLibFile.exists()) {
					tempLibFile.renameTo(combineLibFile);
				}
			}
		}
		if (!combineLibFile.exists()) {
			LOGGER.error(getTaskName() + ": failed generating combined.speclib");
			return false;
		}

		File combineParquet = new File(resultFolderFile, "combined.parquet");
		this.quan_pep_file = new File(resultFolderFile, "combined.pr_matrix.tsv");
		this.quan_pro_file = new File(resultFolderFile, "combined.pg_matrix.tsv");

		diaLibSearchPar.setMbr(true);
		diaLibSearchPar.setDiaFiles(diaFiles);
		diaLibSearchPar.setThreadCount(threadCount);
		diaLibSearchPar.setqValue(0.01);

		diannPar.setThreads(threadCount);
		diannPar.setMatrices(true);
		diannPar.setRelaxed_prot_inf(true);
		diannPar.setReanalyse(true);
		diannPar.setExportDecoy(true);

		if (!combineParquet.exists() || !quan_pep_file.exists()) {
			diaNNTask.addTask(diannPar, diaLibSearchPar, combineParquet.getAbsolutePath(), combineLibFile.getAbsolutePath());
			diaNNTask.run(Math.max(4, diaFiles.length * 12));
		}

		if (!combineParquet.exists()) {
			LOGGER.error(getTaskName() + ": final combined search failed - combined.parquet not found");
			return false;
		}

		updateProgress(50, 100);
		LOGGER.info(getTaskName() + ": library search finished");
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
