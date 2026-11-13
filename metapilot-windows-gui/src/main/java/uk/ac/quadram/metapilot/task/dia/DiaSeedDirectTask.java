package uk.ac.quadram.metapilot.task.dia;

import java.io.File;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;

import javafx.scene.control.ProgressBar;

import uk.ac.quadram.metapilot.core.dia.DiaHostLibrarySupport;
import uk.ac.quadram.metapilot.core.dia.MagPredictedPeptideProteinMapper;
import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNPrecursor;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters.ToolType;
import uk.ac.quadram.metapilot.task.dia.par.MetaParameterDia;
import uk.ac.quadram.metapilot.task.dia.par.MetaParameterDiaTg;
import uk.ac.quadram.metapilot.task.mag.MagDbItem.MagPepLib;

/**
 * DIA Seed workflow task (Direct identification, no iterative search).
 * Behaviour depends on selected library types (Microbiome vs Standalone, from MagPepLib.isMatchToMag()):
 * <ul>
 *   <li>All Standalone: search with all libraries, then stop after identification (no export - no MAG info).</li>
 *   <li>All Microbiome: search with all libraries and run full workflow (iden -> export report, etc.).</li>
 *   <li>Mixed: one search with the microbiome and the standalone (host) libraries together, then the full workflow;
 *   the export removes host peptides (host takes shared peptides), as in every DIA workflow (DiaHostLibrarySupport).</li>
 * </ul>
 */
public class DiaSeedDirectTask extends MetaDiaTask {

	private static final String TASK_NAME = "DIA Seed analysis";
	private static final SimpleDateFormat format = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");

	/** When true, call() stops after iden() without running exportReport(). */
	private boolean stopAfterIden;

	public DiaSeedDirectTask(MetaParameterDia metaPar, MetaResourceParameters mrp, ProgressBar bar1, ProgressBar bar2) {
		super(withMode(metaPar), mrp, bar1, bar2);
		initResultFolderIfNeeded();
	}

	private static MetaParameterDia withMode(MetaParameterDia metaPar) {
		if (metaPar != null) {
			metaPar.setSeedRefinementMode(DiaRefinementMode.DIRECT);
		}
		return metaPar;
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
			updateMessage(getTaskName() + ": Identification completed (no report - standalone libraries only).");
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
		updateSecondaryProgress(1.0);
		logger.info(getTaskName() + ": finished - all workflow steps completed successfully");
		System.out.println(format.format(new Date()) + "\t" + getTaskName()
				+ ": finished - all workflow steps completed successfully");
		return true;
	}

	/** Use new genome/protein scoring scheme (S, p_emp, q_value) like DiaMarkerCoreTask. */
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
			throw new IllegalStateException("Seed Direct requires predicted mapping resources: " + predictedFolder);
		}
        return MagPredictedPeptideProteinMapper.mapForFinalisation(precursors, predictedFolder, null);
	}

	/** Delegate downstream report generation to the same implementation as DiaMarkerCoreTask. */
	@Override
	protected boolean exportReportDownstreamFromNewScoring(java.util.HashMap<String, double[]> genomeScoreMap, java.util.HashMap<String, double[]> proteinScoreMap) {
		return DiaMarkerCoreTask.runExportReportDownstream(this, genomeScoreMap, proteinScoreMap);
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
			}

			if (!stopAfterIden) {
				// MetaParameterDia: resolve paths and classify by matchToMag
				String[] usedLibs = ((MetaParameterDia) metaPar).getLibrary();
				MagPepLib[] allLibs = magDb.getMagPepLibs();

				String[] libPath = new String[usedLibs.length];
				boolean[] matchToMag = new boolean[usedLibs.length];
				for (int i = 0; i < usedLibs.length; i++) {
					File directLibFile = new File(usedLibs[i]);
					if (directLibFile.exists()) {
						libPath[i] = directLibFile.getAbsolutePath();
						matchToMag[i] = false; // external path -> treat as standalone
					} else {
						File chosen = magDb.resolveSpectralLibraryFile(usedLibs[i]);
						if (chosen == null || !chosen.exists()) {
							logger.error(getTaskName()
									+ ": no library file (.speclib/.tsv/.parquet) found for library '{}'", usedLibs[i]);
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
					if (b)
						microCount++;
				}
				boolean allStandalone = (microCount == 0);
				boolean allMicrobiome = (microCount == usedLibs.length);

				if (allStandalone) {
					diaNNTask.addTask(diannPar, diaLibSearchPar, combineParquetFile.getAbsolutePath(), libPath);
					diaNNTask.run(diaFiles.length * 12);
					stopAfterIden = true;
					logger.info(getTaskName() + ": all libraries standalone - will stop after identification.");
				} else if (allMicrobiome) {
					diaNNTask.addTask(diannPar, diaLibSearchPar, combineParquetFile.getAbsolutePath(), libPath);
					diaNNTask.run(diaFiles.length * 12);
					stopAfterIden = false;
				} else {
					// Mixed: Direct has one search, which is the final combined search, so the standalone (host)
					// libraries are searched together with the microbiome libraries (all as .parquet).
					List<File> microLibs = new ArrayList<>();
					List<File> hostLibs = new ArrayList<>();
					for (int i = 0; i < usedLibs.length; i++) {
						File lib = DiaHostLibrarySupport.parquetInLibraryFolder(new File(libPath[i]));
						if (matchToMag[i]) {
							microLibs.add(lib);
						} else if (!lib.getName().toLowerCase().endsWith(".parquet")) {
							logger.error(getTaskName() + ": standalone (host) libraries must be .parquet libraries: {}",
									libPath[i]);
							return false;
						} else {
							hostLibs.add(lib);
						}
					}
					String[] mixedLibs;
					try {
						File diannExe = MetaDiaMarkerCoreRunnerBridge.resolveGuiToolFile(mrp, ToolType.DIANN);
						mixedLibs = DiaHostLibrarySupport.finalSearchLibraries(microLibs, hostLibs, combinedFile,
								diannPar, diannExe, null, MetaDiaMarkerCoreRunnerBridge.inferCondaPrefixFromPython(
										MetaDiaMarkerCoreRunnerBridge.resolveGuiToolFile(mrp, ToolType.PYTHON)));
					} catch (Exception e) {
						logger.error(getTaskName() + ": preparing the combined search libraries failed", e);
						return false;
					}
					diaNNTask.addTask(diannPar, diaLibSearchPar, combineParquetFile.getAbsolutePath(), mixedLibs);
					diaNNTask.run(diaFiles.length * 12);
					stopAfterIden = false;
					logger.info(getTaskName() + ": mixed libraries - microbiome and {} host libraries searched together.",
							hostLibs.size());
				}
			}
		}

		if (!quan_pep_file.exists() || !combineParquetFile.exists()) {
			logCombinedSearchFailed(combineParquetFile, "library search");
			return false;
		}
		return true;
	}
}
