package uk.ac.quadram.metapilot.task.dia;

import java.io.File;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;

import javafx.scene.control.ProgressBar;

import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNPrecursor;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.task.dia.par.MetaParameterDia;
import uk.ac.quadram.metapilot.task.dia.par.MetaParameterDiaTg;
import uk.ac.quadram.metapilot.task.mag.MagDbItem.MagPepLib;
import uk.ac.quadram.metapilot.task.mag.MarkerGenomeSelector;

/**
 * DIA Seed workflow task (Direct identification, no iterative search).
 * Behaviour depends on selected library types (Microbiome vs Standalone, from MagPepLib.isMatchToMag()):
 * <ul>
 *   <li>All Standalone: search with all libraries, then stop after identification (no export – no MAG info).</li>
 *   <li>All Microbiome: search with all libraries and run full workflow (iden → export report, etc.).</li>
 *   <li>Mixed: search microbiome libraries first with --out-lib; full workflow continues (next step TBD by user).</li>
 * </ul>
 */
public class MetaDiaSeedDirectTask extends MetaDiaTask {

	private static final String TASK_NAME = "DIA Seed analysis";
	private static final SimpleDateFormat format = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");

	/** When true, call() stops after iden() without running exportReport(). */
	private boolean stopAfterIden;

	public MetaDiaSeedDirectTask(MetaParameterDia metaPar, MetaResourceParameters mrp, ProgressBar bar1, ProgressBar bar2) {
		super(metaPar, mrp, bar1, bar2);
		initResultFolderIfNeeded();
	}

	@Override
	protected String getTaskName() {
		return TASK_NAME;
	}

	@Override
	protected void initResultFolderIfNeeded() {
		if (metaPar == null) {
			return;
		}
		if (metaPar.getDbSearchResultFile() != null) {
			this.resultFolderFile = metaPar.getDbSearchResultFile();
			return;
		}
		// Direct seed library search: seed_direct under result
		File base = new File(metaPar.getResult() != null ? metaPar.getResult() : ".");
		File dir = new File(base, "seed_direct");
		if (!dir.exists()) {
			dir.mkdirs();
			logger.info("Created result directory for DIA Seed Direct workflow: {}", dir.getAbsolutePath());
		}
		metaPar.setDbSearchResultFile(dir);
		this.resultFolderFile = dir;
	}

	/**
	 * Runs the same steps as MetaDiaTask.call() but skips export (and later steps) when stopAfterIden is true.
	 */
	@Override
	protected Boolean call() throws Exception {
		updateMessage(getTaskName());
		updateProgress(0, 100);

		logger.info(getTaskName() + ": started");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": started");

		updateProgress(5, 100);
		updateSecondaryProgress(ProgressBar.INDETERMINATE_PROGRESS);

		if (!copyReportTemplatesInBackground()) {
			logger.error(getTaskName() + ": Failed to copy report templates");
			return false;
		}

		logger.info("Step 1: Starting spectra format conversion...");
		updateMessage(getTaskName() + ": Converting spectra format...");
		updateProgress(10, 100);
		boolean convert = spectraConvert();
		if (!convert) {
			logger.error(getTaskName() + ": error in spectra format conversion");
			updateMessage(getTaskName() + ": Failed - Spectra conversion error");
			updateSecondaryProgress(0.0);
			return false;
		}

		logger.info("Step 2: Starting sample-specific library generation...");
		updateProgress(15, 100);
		updateMessage(getTaskName() + ": Creating library...");
		boolean createSSDB = createSSLibrary();
		if (!createSSDB) {
			logger.error(getTaskName() + ": error in creating sample-specific library");
			updateMessage(getTaskName() + ": Failed - sample-specific library creation error");
			updateSecondaryProgress(0.0);
			return false;
		}

		logger.info("Step 3: Starting identification...");
		updateProgress(20, 100);
		updateMessage(getTaskName() + ": Running identification...");
		updateSecondaryProgress(ProgressBar.INDETERMINATE_PROGRESS);
		boolean iden = iden();
		if (!iden) {
			logger.error(getTaskName() + ": error in identification");
			updateMessage(getTaskName() + ": Failed - Identification error");
			updateSecondaryProgress(0.0);
			return false;
		}

		if (stopAfterIden) {
			updateProgress(100, 100);
			updateMessage(getTaskName() + ": Identification completed (no report – standalone libraries only).");
			updateSecondaryProgress(0.0);
			logger.info(getTaskName() + ": finished after identification (standalone libraries, no export).");
			System.out.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": finished after identification (standalone libraries, no export).");
			return true;
		}

		// Step 4: Export Report
		logger.info("Step 4: Starting report generation...");
		updateProgress(85, 100);
		updateMessage(getTaskName() + ": Generating reports...");
		updateSecondaryProgress(ProgressBar.INDETERMINATE_PROGRESS);
		boolean exportReport = exportReport();
		if (!exportReport) {
			logger.error(getTaskName() + ": error in report generation");
			updateMessage(getTaskName() + ": Failed - Report generation error");
			updateSecondaryProgress(0.0);
			return false;
		}

		updateProgress(100, 100);
		updateMessage(getTaskName() + ": Completed successfully!");
		updateSecondaryProgress(0.0);
		logger.info(getTaskName() + ": finished - all workflow steps completed successfully");
		System.out.println(format.format(new Date()) + "\t" + getTaskName()
				+ ": finished - all workflow steps completed successfully");
		return true;
	}

	/** Use new genome/protein scoring scheme (S, p_emp, q_value) like MetaDiaMarkerTask. */
	@Override
	protected boolean exportReport() {
		return super.exportReport();
	}

	/**
	 * Map peptides to real genomes using MAG predicted TSVs so genome/protein scoring uses MAG protein IDs.
	 */
	@Override
	protected DiaNNPrecursor[] mapPrecursorsToGenomes(DiaNNPrecursor[] precursors) {
		if (precursors == null || precursors.length == 0 || magDb == null) {
			return precursors;
		}
		File predictedFolder = new File(magDb.getCurrentFile(), "predicted");
		if (!predictedFolder.isDirectory()) {
			logger.warn(getTaskName() + ": predicted folder not found at " + predictedFolder.getAbsolutePath() + ", skipping peptide-to-genome mapping");
			return precursors;
		}
		List<DiaNNPrecursor> mapped = MarkerGenomeSelector.setProteinInfoFromPredictedTsVs(precursors, predictedFolder);
		if (mapped == null || mapped.isEmpty()) {
			logger.warn(getTaskName() + ": no precursors mapped to MAG proteins, using empty set");
			return new DiaNNPrecursor[0];
		}
		return mapped.toArray(new DiaNNPrecursor[0]);
	}

	/** Delegate downstream report generation to the same implementation as MetaDiaMarkerTask. */
	@Override
	protected boolean exportReportDownstreamFromNewScoring(java.util.HashMap<String, double[]> genomeScoreMap, java.util.HashMap<String, double[]> proteinScoreMap) {
		return MetaDiaMarkerTask.runExportReportDownstream(this, genomeScoreMap, proteinScoreMap);
	}

	@Override
	protected boolean iden() {
		logger.info("Starting identification process...");
		updateMessage(getTaskName() + ": peptide and protein identification");

		diaLibSearchPar.setMbr(true);
		diaLibSearchPar.setDiaFiles(diaFiles);
		diaLibSearchPar.setThreadCount(threadCount);
		diaLibSearchPar.setqValue(0.01);

		diannPar.setThreads(threadCount);
		diannPar.setMatrices(true);
		diannPar.setRelaxed_prot_inf(true);
		diannPar.setReanalyse(true);

		File combinedFile = new File(metaPar.getDbSearchResultFile(), "combined");
		if (!combinedFile.exists()) {
			combinedFile.mkdir();
		}
		// Required for exportReport() / getExportQuanPepFile()
		this.quan_pep_file = new File(combinedFile, "combined.pr_matrix.tsv");
		this.quan_pro_file = new File(combinedFile, "combined.pg_matrix.tsv");
		File combineParquetFile = new File(combinedFile, "combined.parquet");

		if (!quan_pep_file.exists() || !combineParquetFile.exists()) {

			logger.info(getTaskName() + ": library search started");

			updateProgress(10, 100);
			diannPar.setExportDecoy(true);

			if (metaPar instanceof MetaParameterDiaTg) {
				// Single external library: no MAG, stop after iden.
				String libraryPath = ((MetaParameterDiaTg) metaPar).getLibraryPath();
				diaNNTask.addTask(diannPar, diaLibSearchPar, combineParquetFile.getAbsolutePath(),
						new String[] { libraryPath });
				diaNNTask.run(diaFiles.length * 12);
				stopAfterIden = true;
				return true;
			}

			// MetaParameterDia: resolve paths and classify by matchToMag
			File magDbFile = magDb.getCurrentFile();
			File libFolderFile = new File(magDbFile, "libraries");
			String[] usedLibs = ((MetaParameterDia) metaPar).getLibrary();
			MagPepLib[] allLibs = magDb.getMagPepLibs();

			String[] libPath = new String[usedLibs.length];
			boolean[] matchToMag = new boolean[usedLibs.length];
			for (int i = 0; i < usedLibs.length; i++) {
				File directLibFile = new File(usedLibs[i]);
				if (directLibFile.exists()) {
					libPath[i] = directLibFile.getAbsolutePath();
					matchToMag[i] = false; // external path → treat as standalone
				} else {
					File libFolder = new File(libFolderFile, usedLibs[i]);
					File speclibFile = new File(libFolder, usedLibs[i] + ".speclib");
					File tsvFile = new File(libFolder, usedLibs[i] + ".tsv");
					File parquetFile = new File(libFolder, usedLibs[i] + ".parquet");

					File chosen = null;
					if (speclibFile.exists()) {
						chosen = speclibFile;
					} else if (tsvFile.exists()) {
						chosen = tsvFile;
					} else if (parquetFile.exists()) {
						chosen = parquetFile;
					}

					if (chosen == null) {
						logger.error(getTaskName() + ": no library file (.speclib/.tsv/.parquet) found for library '{}'", usedLibs[i]);
						return false;
					}

					libPath[i] = chosen.getAbsolutePath();

					boolean micro = false;
					for (MagPepLib lib : allLibs) {
						if (lib.getLibName().equals(usedLibs[i])) {
							micro = lib.isMatchToMag();
							break;
						}
					}
					matchToMag[i] = micro;
				}
			}

			int microCount = 0;
			for (boolean b : matchToMag) {
				if (b) microCount++;
			}
			boolean allStandalone = (microCount == 0);
			boolean allMicrobiome = (microCount == usedLibs.length);

			if (allStandalone) {
				diaNNTask.addTask(diannPar, diaLibSearchPar, combineParquetFile.getAbsolutePath(), libPath);
				diaNNTask.run(diaFiles.length * 12);
				stopAfterIden = true;
				logger.info(getTaskName() + ": all libraries standalone – will stop after identification.");
			} else if (allMicrobiome) {
				diaNNTask.addTask(diannPar, diaLibSearchPar, combineParquetFile.getAbsolutePath(), libPath);
				diaNNTask.run(diaFiles.length * 12);
				stopAfterIden = false;
			} else {
				// Mixed: search microbiome libraries only with --out-lib
				List<String> microPaths = new ArrayList<>();
				for (int i = 0; i < usedLibs.length; i++) {
					if (matchToMag[i]) {
						microPaths.add(libPath[i]);
					}
				}
				String[] microPathArr = microPaths.toArray(new String[0]);
				diaNNTask.addTask(diannPar, diaLibSearchPar, combineParquetFile.getAbsolutePath(),
						false, true, true, microPathArr);
				diaNNTask.run(diaFiles.length * 12);
				stopAfterIden = false;
				logger.info(getTaskName() + ": mixed libraries – searched microbiome libraries with --out-lib.");
			}
		}

		return true;
	}
}
