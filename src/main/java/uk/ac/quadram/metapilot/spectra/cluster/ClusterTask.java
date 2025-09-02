/**
 * 
 */
package uk.ac.quadram.metapilot.spectra.cluster;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import java.util.List;

import org.apache.commons.io.FileUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import uk.ac.quadram.metapilot.core.mod.IsobaricTag;
import uk.ac.ebi.pride.spectracluster.binning.BinningSpectrumConverter;
import uk.ac.ebi.pride.spectracluster.cli.SpectraClusterCliMain;
import uk.ac.ebi.pride.spectracluster.clustering.BinaryClusterFileReference;
import uk.ac.ebi.pride.spectracluster.clustering.BinaryFileClusterer;
import uk.ac.ebi.pride.spectracluster.conversion.MergingCGFConverter;
import uk.ac.ebi.pride.spectracluster.implementation.ClusteringSettings;
import uk.ac.ebi.pride.spectracluster.implementation.SpectraClusterStandalone;
import uk.ac.ebi.pride.spectracluster.merging.BinaryFileMergingClusterer;
import uk.ac.ebi.pride.spectracluster.util.Defaults;
import uk.ac.ebi.pride.spectracluster.util.function.peak.BinnedHighestNPeakFunction;
import uk.ac.ebi.pride.spectracluster.util.function.peak.HighestNPeakFunction;
import uk.ac.ebi.pride.spectracluster.util.function.spectrum.RemoveReporterIonPeaksFunction;

/**
 * @author Kai Cheng
 * @version 2025.03.20
 */
public class ClusterTask {

	private int threadCount;
	private boolean fastMode;
	private IsobaricTag isobaricTag;

	private float startThreshold = 0.999F;
	private float endThreshold = 0.99f;
	private int rounds = 4;
	private int nHighestPeaks = 70;

	private static Logger LOGGER = LogManager.getLogger(ClusterTask.class);
	private static final String taskName = "Spectra clustering";
	private static final SimpleDateFormat format = new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss");

	public ClusterTask(int threadCount, boolean fastMode) {
		this.threadCount = threadCount;
		this.fastMode = fastMode;
		this.isobaricTag = null;
		LOGGER.info("{}: ClusterTask initialized with threadCount={}, fastMode={}, isobaricTag=null", 
			taskName, threadCount, fastMode);
	}

	public ClusterTask(int threadCount, boolean fastMode, IsobaricTag isobaricTag) {
		this.threadCount = threadCount;
		this.fastMode = fastMode;
		this.isobaricTag = isobaricTag;
		
		this.startThreshold = 0.999F;
		this.endThreshold = 0.99f;
		LOGGER.info("{}: ClusterTask initialized with threadCount={}, fastMode={}, isobaricTag={}, startThreshold={}, endThreshold={}", 
			taskName, threadCount, fastMode, isobaricTag != null ? isobaricTag.getName() : "null", startThreshold, endThreshold);
	}

	public void runCli(String[] fileNames, File clusterDir, File resultMgf) {
		runCli(fileNames, clusterDir, resultMgf, false);
	}

	public void runCli(String[] fileNames, File clusterDir, File resultMgf, boolean singleFile) {
		LOGGER.info("{}: Starting CLI-based clustering with {} input files", taskName, fileNames.length);
		LOGGER.info("{}: Cluster directory: {}", taskName, clusterDir.getAbsolutePath());
		LOGGER.info("{}: Result MGF file: {}", taskName, resultMgf.getAbsolutePath());
		LOGGER.info("{}: Single file mode: {}", taskName, singleFile);
		
		// Test CDF loading to ensure the resource issue is resolved
		try {
			LOGGER.info("{}: Testing CDF resource loading...", taskName);
			uk.ac.ebi.pride.spectracluster.cdf.CumulativeDistributionFunction cdf = 
				uk.ac.ebi.pride.spectracluster.cdf.CumulativeDistributionFunctionFactory
					.getDefaultCumlativeDistributionFunctionForSimilarityMetric(
						uk.ac.ebi.pride.spectracluster.similarity.CombinedFisherIntensityTest.class);
			LOGGER.info("{}: CDF loaded successfully: {}", taskName, cdf);
		} catch (Exception e) {
			LOGGER.warn("{}: CDF loading test failed, but continuing with clustering: {}", taskName, e.getMessage());
		}

		// Validate input files
		if (!validateInputFiles(fileNames)) {
			LOGGER.error("{}: Input file validation failed, aborting clustering", taskName);
			return;
		}

		if (!clusterDir.exists()) {
			LOGGER.info("{}: Creating cluster directory: {}", taskName, clusterDir.getAbsolutePath());
			boolean created = clusterDir.mkdirs();
			if (!created) {
				LOGGER.error("{}: Failed to create cluster directory: {}", taskName, clusterDir.getAbsolutePath());
				return;
			}
		} else {
			LOGGER.info("{}: Cluster directory already exists: {}", taskName, clusterDir.getAbsolutePath());
		}

		ArrayList<String> list = new ArrayList<String>();

		for (int i = 0; i < fileNames.length; i++) {
			list.add(fileNames[i]);
			LOGGER.debug("{}: Added input file {}: {}", taskName, i, fileNames[i]);
		}

		if (fastMode) {
			list.add("-fast_mode");
			LOGGER.info("{}: Added fast mode flag", taskName);
		}

		list.add("-filter");
		list.add("mz_150");
		LOGGER.debug("{}: Added mz_150 filter", taskName);

		list.add("-rounds");
		list.add("4");
		LOGGER.debug("{}: Added 4 rounds parameter", taskName);

		list.add("-binary_directory");
		list.add(clusterDir.getAbsolutePath());
		LOGGER.debug("{}: Added binary directory: {}", taskName, clusterDir.getAbsolutePath());

		// FILTER PEAKS PER MZ new BinnedHighestNPeakFunction(10, 100, 30)
		list.add("-x_filter_peaks_mz");
		LOGGER.debug("{}: Added peak filtering parameter", taskName);

		list.add("-threshold_start");
		list.add(String.valueOf(startThreshold));
		LOGGER.debug("{}: Added start threshold: {}", taskName, startThreshold);

		list.add("-threshold_end");
		list.add(String.valueOf(endThreshold));
		LOGGER.debug("{}: Added end threshold: {}", taskName, endThreshold);

		File cfgFile = new File(resultMgf.getParent(), "cluster.cfg");
		list.add("-output_path");
		list.add(cfgFile.getAbsolutePath());
		LOGGER.debug("{}: Added output path: {}", taskName, cfgFile.getAbsolutePath());

		list.add("-major_peak_jobs");
		list.add(String.valueOf(threadCount));
		LOGGER.debug("{}: Added thread count: {}", taskName, threadCount);

		if (isobaricTag != null) {
			String name = isobaricTag.getName();
			LOGGER.info("{}: Processing isobaric tag: {}", taskName, name);
			if (name.startsWith("iTRAQ")) {
				list.add("-remove_reporters");
				list.add("itraq");
				LOGGER.debug("{}: Added iTRAQ reporter removal", taskName);
			} else {
				list.add("-remove_reporters");
				list.add("all");
				LOGGER.debug("{}: Added all reporter removal", taskName);
			}
		} else {
			LOGGER.debug("{}: No isobaric tag specified, skipping reporter removal", taskName);
		}

		String[] args = list.toArray(new String[list.size()]);
		LOGGER.info("{}: Executing SpectraClusterCliMain with {} arguments", taskName, args.length);
		LOGGER.debug("{}: CLI arguments: {}", taskName, String.join(" ", args));
		
		try {
			SpectraClusterCliMain.main(args);
			LOGGER.info("{}: SpectraClusterCliMain completed successfully", taskName);
		} catch (Exception e) {
			LOGGER.error("{}: SpectraClusterCliMain failed with exception", taskName, e);
			return;
		}

		LOGGER.info("{}: Running garbage collection", taskName);
		System.gc();

		LOGGER.info("{}: Processing generated files in cluster directory", taskName);
		File[] generatedFiles = clusterDir.listFiles();
		if (generatedFiles == null) {
			LOGGER.error("{}: Failed to list files in cluster directory", taskName);
			return;
		}
		
		LOGGER.info("{}: Found {} files in cluster directory", taskName, generatedFiles.length);
		for (File file : generatedFiles) {
			LOGGER.debug("{}: Generated file: {} (size: {} bytes, isDirectory: {})", 
				taskName, file.getName(), file.length(), file.isDirectory());
		}
		
		File cgfFile = null;
		for (int i = 0; i < generatedFiles.length; i++) {
			if (generatedFiles[i].isDirectory()) {
				LOGGER.debug("{}: Processing directory: {}", taskName, generatedFiles[i].getName());
				File[] subFiles = generatedFiles[i].listFiles();
				if (subFiles != null) {
					for (int j = 0; j < subFiles.length; j++) {
						LOGGER.debug("{}: Deleting sub-file: {}", taskName, subFiles[j].getName());
						FileUtils.deleteQuietly(subFiles[j]);
					}
				}
				LOGGER.debug("{}: Deleting directory: {}", taskName, generatedFiles[i].getName());
				FileUtils.deleteQuietly(generatedFiles[i]);
			} else if (generatedFiles[i].isFile()) {
				if (generatedFiles[i].getName().endsWith("cgf")) {
					cgfFile = generatedFiles[i];
					LOGGER.info("{}: Found CGF file: {} (size: {} bytes)", taskName, cgfFile.getName(), cgfFile.length());
				} else {
					LOGGER.debug("{}: Deleting non-CGF file: {}", taskName, generatedFiles[i].getName());
					FileUtils.deleteQuietly(generatedFiles[i]);
				}
			}
		}

		if (cgfFile == null) {
			LOGGER.error("{}: No CGF file found after clustering, cannot proceed with conversion", taskName);
			return;
		}

		int clusterCount = 0;
		LOGGER.info("{}: Converting CGF file to MGF: {} -> {}", taskName, cgfFile.getAbsolutePath(), resultMgf.getAbsolutePath());

		try {
			clusterCount = CGF2MGF.convert2Mgf(fileNames, cgfFile, resultMgf, false);
			LOGGER.info("{}: CGF to MGF conversion completed successfully, {} spectra clustered", taskName, clusterCount);
		} catch (IOException e) {
			LOGGER.error("{}: Error in converting the CGF file to MGF file", taskName, e);
			System.err.println(
					format.format(new Date()) + "\t" + taskName + ": error in converting the CGF file to MGF file.");
		}

		// Validate output file
		if (resultMgf.exists()) {
			LOGGER.info("{}: Result MGF file created successfully: {} (size: {} bytes)", 
				taskName, resultMgf.getAbsolutePath(), resultMgf.length());
		} else {
			LOGGER.error("{}: Result MGF file was not created", taskName);
		}

		System.out.println(
				format.format(new Date()) + "\t" + taskName + ": clustered MS/MS spectra count: " + clusterCount);
		LOGGER.info(taskName + ": clustered MS/MS spectra count: " + clusterCount);
	}

	/**
	 * Validates input files before processing
	 */
	private boolean validateInputFiles(String[] fileNames) {
		LOGGER.info("{}: Validating {} input files", taskName, fileNames.length);
		
		if (fileNames == null || fileNames.length == 0) {
			LOGGER.error("{}: No input files provided", taskName);
			return false;
		}

		for (int i = 0; i < fileNames.length; i++) {
			String fileName = fileNames[i];
			LOGGER.debug("{}: Validating input file {}: {}", taskName, i, fileName);
			
			if (fileName == null || fileName.trim().isEmpty()) {
				LOGGER.error("{}: Input file {} is null or empty", taskName, i);
				return false;
			}

			File file = new File(fileName);
			LOGGER.debug("{}: Input file {}: absolute path = {}, exists = {}, isFile = {}, canRead = {}, size = {} bytes", 
				taskName, i, file.getAbsolutePath(), file.exists(), file.isFile(), file.canRead(), file.length());

			if (!file.exists()) {
				LOGGER.error("{}: Input file {} does not exist: {}", taskName, i, file.getAbsolutePath());
				return false;
			}

			if (!file.isFile()) {
				LOGGER.error("{}: Input file {} is not a regular file: {}", taskName, i, file.getAbsolutePath());
				return false;
			}

			if (!file.canRead()) {
				LOGGER.error("{}: Input file {} is not readable: {}", taskName, i, file.getAbsolutePath());
				return false;
			}

			if (file.length() == 0) {
				LOGGER.error("{}: Input file {} is empty: {}", taskName, i, file.getAbsolutePath());
				return false;
			}

			// Check if it's a valid MGF file by looking for BEGIN IONS
			if (!validateMGFContent(file)) {
				LOGGER.error("{}: Input file {} does not appear to be a valid MGF file: {}", taskName, i, file.getAbsolutePath());
				return false;
			}
		}

		LOGGER.info("{}: All {} input files validated successfully", taskName, fileNames.length);
		return true;
	}

	/**
	 * Validates that a file contains valid MGF content
	 */
	private boolean validateMGFContent(File file) {
		LOGGER.debug("{}: Validating MGF content for file: {}", taskName, file.getName());
		
		try (BufferedReader reader = new BufferedReader(new FileReader(file))) {
			String line;
			int lineCount = 0;
			int ionCount = 0;
			boolean foundBeginIons = false;
			
			while ((line = reader.readLine()) != null && lineCount < 1000) { // Check first 1000 lines
				lineCount++;
				if (line.trim().startsWith("BEGIN IONS")) {
					foundBeginIons = true;
					ionCount++;
				}
			}
			
			LOGGER.debug("{}: File {} validation - lines checked: {}, BEGIN IONS found: {}, ion count: {}", 
				taskName, file.getName(), lineCount, foundBeginIons, ionCount);
			
			return foundBeginIons && ionCount > 0;
		} catch (IOException e) {
			LOGGER.error("{}: Error reading file {} for MGF validation", taskName, file.getName(), e);
			return false;
		}
	}

	/**
	 * @deprecated
	 * @param fileNames
	 * @param clusterDir
	 * @param log
	 * @param resultName
	 */
	public void clustering(String[] fileNames, File clusterDir, Appendable log, String resultName) {
		LOGGER.info("{}: Starting deprecated clustering method with {} input files", taskName, fileNames.length);
		LOGGER.info("{}: Cluster directory: {}, result name: {}", taskName, clusterDir.getAbsolutePath(), resultName);

		try {
			log.append(format.format(new Date()) + "\tclustering started\n");
		} catch (IOException e) {
			LOGGER.error("{}: Error in writing log file", taskName, e);
		}
		
		if (!clusterDir.exists()) {
			LOGGER.info("{}: Creating cluster directory: {}", taskName, clusterDir.getAbsolutePath());
			boolean created = clusterDir.mkdirs();
			if (!created) {
				LOGGER.error("{}: Failed to create cluster directory: {}", taskName, clusterDir.getAbsolutePath());
				return;
			}
		} else {
			LOGGER.info("{}: Cluster directory already exists: {}", taskName, clusterDir.getAbsolutePath());
		}

		int totalCount = 0;
		LOGGER.info("{}: Counting MS/MS spectra in input files", taskName);
		
		for (int i = 0; i < fileNames.length; i++) {
			File filei = new File(fileNames[i]);
			LOGGER.debug("{}: Processing input file {}: {} (size: {} bytes)", taskName, i, filei.getName(), filei.length());
			
			int count = 0;
			String line = null;
			
			try (BufferedReader reader = new BufferedReader(new FileReader(filei))) {
				while ((line = reader.readLine()) != null) {
					if (line.startsWith("BEGIN IONS")) {
						count++;
					}
				}
				reader.close();
				LOGGER.debug("{}: File {} contains {} MS/MS spectra", taskName, filei.getName(), count);
			} catch (IOException e) {
				LOGGER.error("{}: Error in reading mgf file {}", taskName, filei.getName(), e);
			}

			totalCount += count;

			try {
				log.append(format.format(new Date()) + "\t" + filei.getName() + " MS/MS count: " + count + "\n");
			} catch (IOException e) {
				LOGGER.error("{}: Error in writing log file", taskName, e);
			}
		}
		
		LOGGER.info("{}: Total MS/MS count across all files: {}", taskName, totalCount);
		
		try {
			log.append(format.format(new Date()) + "\tTotal MS/MS count: " + totalCount + "\n");
			log.append(format.format(new Date()) + "\tGenerating binning spectrum started\n");
		} catch (IOException e) {
			LOGGER.error("{}: Error in writing log file", taskName, e);
		}

		LOGGER.info("{}: Creating BinningSpectrumConverter with threadCount={}, fastMode={}", taskName, threadCount, fastMode);
		BinningSpectrumConverter binningSpectrumConverter = new BinningSpectrumConverter(clusterDir, threadCount, fastMode);

		try {
			LOGGER.info("{}: Processing peaklist files for binning", taskName);
			binningSpectrumConverter.processPeaklistFiles(fileNames);
			LOGGER.info("{}: Peaklist processing completed successfully", taskName);
		} catch (Exception e) {
			LOGGER.error("{}: Error in clustering spectra during binning", taskName, e);
		}
		
		List<BinaryClusterFileReference> binaryFiles = binningSpectrumConverter.getWrittenFiles();
		LOGGER.info("{}: Generated {} binary files from binning", taskName, binaryFiles.size());
		for (BinaryClusterFileReference ref : binaryFiles) {
			LOGGER.debug("{}: Binary file: {}", taskName, ref.getResultFile().getName());
		}

		File clusteringResultDirectory = createTemporaryDirectory("clustering_results", clusterDir);
		File tmpClusteringResultDirectory = createTemporaryDirectory("clustering_results", clusterDir);
		List<Float> thresholds = generateThreshold(0.999f, 0.99f, 4);
		
		LOGGER.info("{}: Created clustering directories: {} and {}", taskName, 
			clusteringResultDirectory.getAbsolutePath(), tmpClusteringResultDirectory.getAbsolutePath());
		LOGGER.info("{}: Generated {} clustering thresholds: {}", taskName, thresholds.size(), thresholds);

		try {
			log.append(format.format(new Date()) + "\tThresholds: ");
			for (Float f : thresholds) {
				log.append(String.valueOf(f)).append("\t");
			}
			log.append("\n");
			log.append(format.format(new Date()) + "\tClustering MS/MS spectra started\n");
		} catch (IOException e) {
			LOGGER.error("{}: Error in writing log file", taskName, e);
		}

		LOGGER.info("{}: Creating BinaryFileClusterer with threadCount={}, fastMode={}", taskName, threadCount, fastMode);
		BinaryFileClusterer binaryFileClusterer = new BinaryFileClusterer(threadCount, clusteringResultDirectory,
				thresholds, fastMode, tmpClusteringResultDirectory, null);

		try {
			LOGGER.info("{}: Starting binary file clustering with {} files", taskName, binaryFiles.size());
			binaryFileClusterer.clusterFiles(binaryFiles);
			LOGGER.info("{}: Binary file clustering completed successfully", taskName);
		} catch (Exception e) {
			LOGGER.error("{}: Error in clustering spectra during binary clustering", taskName, e);
		}

		List<BinaryClusterFileReference> clusteredFiles = binaryFileClusterer.getResultFiles();
		clusteredFiles = new ArrayList<BinaryClusterFileReference>(clusteredFiles);
		Collections.sort(clusteredFiles);
		LOGGER.info("{}: Binary clustering generated {} result files", taskName, clusteredFiles.size());

		File mergedResultsDirectory = createTemporaryDirectory("merged_results", clusterDir);
		LOGGER.info("{}: Created merged results directory: {}", taskName, mergedResultsDirectory.getAbsolutePath());
		
		BinaryFileMergingClusterer mergingClusterer = new BinaryFileMergingClusterer(threadCount,
				mergedResultsDirectory, thresholds, fastMode, Defaults.getDefaultPrecursorIonTolerance(), true,
				mergedResultsDirectory);

		try {
			log.append(format.format(new Date()) + "\tMerging clusters started\n");
		} catch (IOException e) {
			LOGGER.error("{}: Error in writing log file", taskName, e);
		}

		File combinedResultFile = null;
		try {
			combinedResultFile = File.createTempFile("combined_clustering_results", ".cgf", clusterDir);
			LOGGER.info("{}: Created combined result file: {}", taskName, combinedResultFile.getAbsolutePath());
		} catch (IOException e) {
			LOGGER.error("{}: Error in creating temp file", taskName, e);
			return;
		}
		
		MergingCGFConverter mergingCGFConverter = new MergingCGFConverter(combinedResultFile, true);
		mergingClusterer.addListener(mergingCGFConverter);

		try {
			LOGGER.info("{}: Starting merging clustering with {} files", taskName, clusteredFiles.size());
			mergingClusterer.clusterFiles(clusteredFiles);
			LOGGER.info("{}: Merging clustering completed successfully", taskName);
		} catch (Exception e) {
			LOGGER.error("{}: Error in clustering spectra during merging", taskName, e);
		}

		LOGGER.info("{}: Running garbage collection", taskName);
		System.gc();

		LOGGER.info("{}: Processing generated files in cluster directory", taskName);
		File[] generatedFiles = clusterDir.listFiles();
		if (generatedFiles == null) {
			LOGGER.error("{}: Failed to list files in cluster directory", taskName);
			return;
		}
		
		LOGGER.info("{}: Found {} files in cluster directory", taskName, generatedFiles.length);
		for (File file : generatedFiles) {
			LOGGER.debug("{}: Generated file: {} (size: {} bytes, isDirectory: {})", 
				taskName, file.getName(), file.length(), file.isDirectory());
		}
		
		File cgfFile = null;
		for (int i = 0; i < generatedFiles.length; i++) {
			if (generatedFiles[i].isDirectory()) {
				LOGGER.debug("{}: Processing directory: {}", taskName, generatedFiles[i].getName());
				File[] subFiles = generatedFiles[i].listFiles();
				if (subFiles != null) {
					for (int j = 0; j < subFiles.length; j++) {
						LOGGER.debug("{}: Deleting sub-file: {}", taskName, subFiles[j].getName());
						FileUtils.deleteQuietly(subFiles[j]);
					}
				}
				LOGGER.debug("{}: Deleting directory: {}", taskName, generatedFiles[i].getName());
				FileUtils.deleteQuietly(generatedFiles[i]);
			} else if (generatedFiles[i].isFile()) {
				if (generatedFiles[i].getName().endsWith("cgf")) {
					cgfFile = new File(clusterDir, "cluster.cgf");
					LOGGER.info("{}: Renaming CGF file {} to {}", taskName, generatedFiles[i].getName(), cgfFile.getName());
					boolean renamed = generatedFiles[i].renameTo(cgfFile);
					if (!renamed) {
						LOGGER.error("{}: Failed to rename CGF file", taskName);
					}
				} else {
					LOGGER.debug("{}: Deleting non-CGF file: {}", taskName, generatedFiles[i].getName());
					FileUtils.deleteQuietly(generatedFiles[i]);
				}
			}
		}

		if (cgfFile == null) {
			LOGGER.error("{}: No CGF file found after clustering, cannot proceed with conversion", taskName);
			return;
		}

		File mgfFile = new File(clusterDir, resultName);
		LOGGER.info("{}: Converting CGF to MGF: {} -> {}", taskName, cgfFile.getAbsolutePath(), mgfFile.getAbsolutePath());
		
		int clusterCount = CGF2MGF.convertMultiFile(cgfFile, mgfFile);
		LOGGER.info("{}: CGF to MGF conversion completed, {} spectra clustered", taskName, clusterCount);

		try {
			log.append(format.format(new Date()) + "\tClustered MS/MS count: " + clusterCount + "\n");
			log.append(format.format(new Date()) + "\tClustering MS/MS spectra finished\n");
		} catch (IOException e) {
			LOGGER.error("{}: Error in writing log file", taskName, e);
		}
		
		// Validate output file
		if (mgfFile.exists()) {
			LOGGER.info("{}: Result MGF file created successfully: {} (size: {} bytes)", 
				taskName, mgfFile.getAbsolutePath(), mgfFile.length());
		} else {
			LOGGER.error("{}: Result MGF file was not created", taskName);
		}
	}

	/**
	 * @deprecated
	 * @param fileNames
	 * @param clusterDir
	 * @param resultMgf
	 */
	public void clusterStandAlone(String[] fileNames, File clusterDir, File resultMgf) {
		LOGGER.info("{}: Starting standalone clustering with {} input files", taskName, fileNames.length);
		LOGGER.info("{}: Cluster directory: {}, result MGF: {}", taskName, clusterDir.getAbsolutePath(), resultMgf.getAbsolutePath());

		if (!clusterDir.exists()) {
			LOGGER.info("{}: Creating cluster directory: {}", taskName, clusterDir.getAbsolutePath());
			boolean created = clusterDir.mkdirs();
			if (!created) {
				LOGGER.error("{}: Failed to create cluster directory: {}", taskName, clusterDir.getAbsolutePath());
				return;
			}
		} else {
			LOGGER.info("{}: Cluster directory already exists: {}", taskName, clusterDir.getAbsolutePath());
		}

		LOGGER.info("{}: Creating SpectraClusterStandalone instance", taskName);
		SpectraClusterStandalone spectraClusterStandalone = null;
		try {
			spectraClusterStandalone = new SpectraClusterStandalone(clusterDir);
			LOGGER.info("{}: SpectraClusterStandalone created successfully", taskName);
		} catch (Exception e) {
			LOGGER.error("{}: Error in creating SpectraClusterStandalone", taskName, e);
			System.err.println(format.format(new Date()) + "\t" + taskName + ": error in clustering spectra");
			return;
		}

		LOGGER.info("{}: Configuring clustering settings", taskName);
		ClusteringSettings
				.addIntitalSpectrumFilter(ClusteringSettings.SPECTRUM_FILTER.getFilterForName("immonium_ions").filter);
		LOGGER.debug("{}: Added immonium ions filter", taskName);

		ClusteringSettings
				.addIntitalSpectrumFilter(ClusteringSettings.SPECTRUM_FILTER.getFilterForName("mz_200").filter);
		LOGGER.debug("{}: Added mz_200 filter", taskName);

		ClusteringSettings.setLoadingSpectrumFilter(new HighestNPeakFunction(nHighestPeaks));
		LOGGER.debug("{}: Set loading spectrum filter to HighestNPeakFunction with {} peaks", taskName, nHighestPeaks);

		ClusteringSettings.setLoadingSpectrumFilter(new BinnedHighestNPeakFunction(10, 100, 30));
		LOGGER.debug("{}: Set loading spectrum filter to BinnedHighestNPeakFunction(10, 100, 30)", taskName);

		spectraClusterStandalone.setParallelJobs(threadCount);
		LOGGER.info("{}: Set parallel jobs to {}", taskName, threadCount);

		List<Float> thresholds = spectraClusterStandalone.generateClusteringThresholds(startThreshold, endThreshold, rounds);
		LOGGER.info("{}: Generated {} clustering thresholds: {}", taskName, thresholds.size(), thresholds);

		spectraClusterStandalone.setUseFastMode(fastMode);
		LOGGER.info("{}: Set fast mode to {}", taskName, fastMode);

		if (isobaricTag == null) {
			LOGGER.debug("{}: No isobaric tag specified, skipping reporter removal", taskName);
		} else {
			String name = isobaricTag.getName();
			LOGGER.info("{}: Processing isobaric tag: {}", taskName, name);
			if (name.startsWith("iTRAQ")) {
				spectraClusterStandalone.setReporterType(RemoveReporterIonPeaksFunction.REPORTER_TYPE.ITRAQ);
				LOGGER.debug("{}: Set reporter type to ITRAQ", taskName);
			} else {
				spectraClusterStandalone.setReporterType(RemoveReporterIonPeaksFunction.REPORTER_TYPE.ALL);
				LOGGER.debug("{}: Set reporter type to ALL", taskName);
			}
		}

		ArrayList<File> peaklistFiles = new ArrayList<File>();
		for (int i = 0; i < fileNames.length; i++) {
			File file = new File(fileNames[i]);
			peaklistFiles.add(file);
			LOGGER.debug("{}: Added peaklist file {}: {} (size: {} bytes)", taskName, i, file.getName(), file.length());
		}
		
		LOGGER.info("{}: Added {} peaklist files for clustering", taskName, peaklistFiles.size());
		LOGGER.debug("{}: Loading spectrum filter class: {}", taskName, ClusteringSettings.getLoadingSpectrumFilter().getClass());
		
		File finalResultFile = new File(clusterDir, "cluster.cgf");
		LOGGER.info("{}: Final result file will be: {}", taskName, finalResultFile.getAbsolutePath());
		
		try {
			LOGGER.info("{}: Starting clustering of {} peaklist files", taskName, peaklistFiles.size());
			spectraClusterStandalone.clusterPeaklistFiles(peaklistFiles, thresholds, finalResultFile);
			LOGGER.info("{}: Clustering completed successfully", taskName);
		} catch (Exception e) {
			LOGGER.error("{}: Error in clustering spectra during standalone clustering", taskName, e);
			System.err.println(format.format(new Date()) + "\t" + taskName + ": error in clustering spectra");
		}

		LOGGER.info("{}: Running garbage collection", taskName);
		System.gc();

		LOGGER.info("{}: Processing generated files in cluster directory", taskName);
		File[] generatedFiles = clusterDir.listFiles();
		if (generatedFiles == null) {
			LOGGER.error("{}: Failed to list files in cluster directory", taskName);
			return;
		}
		
		LOGGER.info("{}: Found {} files in cluster directory", taskName, generatedFiles.length);
		for (File file : generatedFiles) {
			LOGGER.debug("{}: Generated file: {} (size: {} bytes, isDirectory: {})", 
				taskName, file.getName(), file.length(), file.isDirectory());
		}
		
		File cgfFile = null;
		for (int i = 0; i < generatedFiles.length; i++) {
			if (generatedFiles[i].isDirectory()) {
				LOGGER.debug("{}: Processing directory: {}", taskName, generatedFiles[i].getName());
				File[] subFiles = generatedFiles[i].listFiles();
				if (subFiles != null) {
					for (int j = 0; j < subFiles.length; j++) {
						LOGGER.debug("{}: Deleting sub-file: {}", taskName, subFiles[j].getName());
						FileUtils.deleteQuietly(subFiles[j]);
					}
				}
				LOGGER.debug("{}: Deleting directory: {}", taskName, generatedFiles[i].getName());
				FileUtils.deleteQuietly(generatedFiles[i]);
			} else if (generatedFiles[i].isFile()) {
				if (generatedFiles[i].getName().endsWith("cgf")) {
					cgfFile = generatedFiles[i];
					LOGGER.info("{}: Found CGF file: {} (size: {} bytes)", taskName, cgfFile.getName(), cgfFile.length());
				} else {
					LOGGER.debug("{}: Deleting non-CGF file: {}", taskName, generatedFiles[i].getName());
					FileUtils.deleteQuietly(generatedFiles[i]);
				}
			}
		}

		if (cgfFile == null) {
			LOGGER.error("{}: No CGF file found after clustering, cannot proceed with conversion", taskName);
			return;
		}

		LOGGER.info("{}: Converting CGF to MGF: {} -> {}", taskName, cgfFile.getAbsolutePath(), resultMgf.getAbsolutePath());
		int clusterCount = CGF2MGF.convertMultiFile(cgfFile, resultMgf);
		LOGGER.info("{}: CGF to MGF conversion completed, {} spectra clustered", taskName, clusterCount);

		System.out.println(
				format.format(new Date()) + "\t" + taskName + ": clustered MS/MS spectra count: " + clusterCount);
		LOGGER.info(taskName + ": clustered MS/MS spectra count: " + clusterCount);
		
		// Validate output file
		if (resultMgf.exists()) {
			LOGGER.info("{}: Result MGF file created successfully: {} (size: {} bytes)", 
				taskName, resultMgf.getAbsolutePath(), resultMgf.length());
		} else {
			LOGGER.error("{}: Result MGF file was not created", taskName);
		}
	}

	/**
	 * Can't remove report ions from isobaric labeling quantification, not used
	 * 
	 * @deprecated
	 * @param fileNames
	 * @param clusterDir
	 * @param resultMgf
	 */
	public void clustering(String[] fileNames, File clusterDir, File resultMgf) {
		LOGGER.info("{}: Starting deprecated clustering method (no reporter removal) with {} input files", taskName, fileNames.length);
		LOGGER.info("{}: Cluster directory: {}, result MGF: {}", taskName, clusterDir.getAbsolutePath(), resultMgf.getAbsolutePath());

		if (!clusterDir.exists()) {
			LOGGER.info("{}: Creating cluster directory: {}", taskName, clusterDir.getAbsolutePath());
			boolean created = clusterDir.mkdirs();
			if (!created) {
				LOGGER.error("{}: Failed to create cluster directory: {}", taskName, clusterDir.getAbsolutePath());
				return;
			}
		} else {
			LOGGER.info("{}: Cluster directory already exists: {}", taskName, clusterDir.getAbsolutePath());
		}

		LOGGER.info("{}: Creating BinningSpectrumConverter with threadCount={}, fastMode={}", taskName, threadCount, fastMode);
		BinningSpectrumConverter binningSpectrumConverter = new BinningSpectrumConverter(clusterDir, threadCount, fastMode);

		try {
			LOGGER.info("{}: Processing peaklist files for binning", taskName);
			binningSpectrumConverter.processPeaklistFiles(fileNames);
			LOGGER.info("{}: Peaklist processing completed successfully", taskName);
		} catch (Exception e) {
			LOGGER.error("{}: Error in clustering spectra during binning", taskName, e);
		}
		
		List<BinaryClusterFileReference> binaryFiles = binningSpectrumConverter.getWrittenFiles();
		LOGGER.info("{}: Generated {} binary files from binning", taskName, binaryFiles.size());
		for (BinaryClusterFileReference ref : binaryFiles) {
			LOGGER.debug("{}: Binary file: {}", taskName, ref.getResultFile().getName());
		}

		File clusteringResultDirectory = createTemporaryDirectory("clustering_results", clusterDir);
		File tmpClusteringResultDirectory = createTemporaryDirectory("clustering_results", clusterDir);
		List<Float> thresholds = generateThreshold(startThreshold, endThreshold, rounds);

		BinaryFileClusterer binaryFileClusterer = new BinaryFileClusterer(threadCount, clusteringResultDirectory,
				thresholds, fastMode, tmpClusteringResultDirectory, null);

		try {
			binaryFileClusterer.clusterFiles(binaryFiles);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			LOGGER.error("Error in clutering spectra", e);
		}

		List<BinaryClusterFileReference> clusteredFiles = binaryFileClusterer.getResultFiles();
		clusteredFiles = new ArrayList<BinaryClusterFileReference>(clusteredFiles);
		Collections.sort(clusteredFiles);

		File mergedResultsDirectory = createTemporaryDirectory("merged_results", clusterDir);
		BinaryFileMergingClusterer mergingClusterer = new BinaryFileMergingClusterer(threadCount,
				mergedResultsDirectory, thresholds, fastMode, Defaults.getDefaultPrecursorIonTolerance(), true,
				mergedResultsDirectory);

		File combinedResultFile = null;
		try {
			combinedResultFile = File.createTempFile("combined_clustering_results", ".cgf", clusterDir);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			LOGGER.error("Error in creating temp file", e);
		}
		MergingCGFConverter mergingCGFConverter = new MergingCGFConverter(combinedResultFile, true);
		mergingClusterer.addListener(mergingCGFConverter);

		try {
			mergingClusterer.clusterFiles(clusteredFiles);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			LOGGER.error("Error in clutering spectra", e);
		}

		System.gc();

		File[] generatedFiles = clusterDir.listFiles();
		File cgfFile = null;
		for (int i = 0; i < generatedFiles.length; i++) {
			if (generatedFiles[i].isDirectory()) {
				File[] subFiles = generatedFiles[i].listFiles();
				for (int j = 0; j < subFiles.length; j++) {
					FileUtils.deleteQuietly(subFiles[j]);
				}
				FileUtils.deleteQuietly(generatedFiles[i]);

			} else if (generatedFiles[i].isFile()) {
				if (generatedFiles[i].getName().endsWith("cgf")) {
					cgfFile = new File(clusterDir, "cluster.cgf");
					generatedFiles[i].renameTo(cgfFile);
				} else {
					FileUtils.deleteQuietly(generatedFiles[i]);
				}
			}
		}

		int clusterCount = CGF2MGF.convertMultiFile(cgfFile, resultMgf);

		LOGGER.info(taskName + ": Clustered MS/MS spectra count: " + clusterCount);
	}

	private File createTemporaryDirectory(String prefix, File tmpDirectory) {
		LOGGER.debug("{}: Creating temporary directory with prefix '{}' in {}", taskName, prefix, 
			tmpDirectory != null ? tmpDirectory.getAbsolutePath() : "default temp location");
		
		File tmpFile = null;

		if (tmpDirectory != null && tmpDirectory.isDirectory()) {
			try {
				tmpFile = File.createTempFile(prefix, "", tmpDirectory);
				LOGGER.debug("{}: Created temp file in specified directory: {}", taskName, tmpFile.getAbsolutePath());
			} catch (IOException e) {
				LOGGER.error("{}: Error in creating temp file in specified directory", taskName, e);
				System.err.println(format.format(new Date()) + "\t" + taskName + ": error in creating temp file.");
			}
		} else {
			try {
				tmpFile = File.createTempFile(prefix, "");
				LOGGER.debug("{}: Created temp file in default location: {}", taskName, tmpFile.getAbsolutePath());
			} catch (IOException e) {
				LOGGER.error("{}: Error in creating temp file in default location", taskName, e);
				System.err.println(format.format(new Date()) + "\t" + taskName + ": error in creating temp file.");
			}
		}

		if (tmpFile != null && !tmpFile.delete()) {
			LOGGER.error("{}: Error in deleting temp file: {}", taskName, tmpFile.getAbsolutePath());
			System.err.println(format.format(new Date()) + "\t" + taskName + ": error in deleting temp file.");
			return tmpFile;
		}

		if (tmpFile != null && !tmpFile.mkdir()) {
			LOGGER.error("{}: Error in creating temp directory: {}", taskName, tmpFile.getAbsolutePath());
			System.err.println(format.format(new Date()) + "\t" + taskName + ": error in creating temp file.");
		} else if (tmpFile != null) {
			LOGGER.debug("{}: Successfully created temporary directory: {}", taskName, tmpFile.getAbsolutePath());
		}

		return tmpFile;
	}

	private List<Float> generateThreshold(float startThreshold, float endThreshold, int rounds) {
		LOGGER.debug("{}: Generating {} thresholds from {} to {}", taskName, rounds, startThreshold, endThreshold);
		
		List<Float> thresholds = new ArrayList<Float>(rounds);
		float stepSize = (startThreshold - endThreshold) / (rounds - 1);

		for (int i = 0; i < rounds; i++) {
			float threshold = startThreshold - (stepSize * i);
			thresholds.add(threshold);
			LOGGER.debug("{}: Threshold {}: {}", taskName, i, threshold);
		}

		LOGGER.debug("{}: Generated {} thresholds: {}", taskName, thresholds.size(), thresholds);
		return thresholds;
	}

}
