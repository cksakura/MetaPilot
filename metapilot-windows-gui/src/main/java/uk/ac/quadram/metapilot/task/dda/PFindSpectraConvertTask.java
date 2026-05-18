package uk.ac.quadram.metapilot.task.dda;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintWriter;
import java.nio.file.*;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Date;
import java.util.Locale;
import java.util.concurrent.atomic.AtomicInteger;

import org.apache.logging.log4j.Logger;

import uk.ac.quadram.metapilot.core.runner.MzmlSummaryFileWriter;
import uk.ac.quadram.metapilot.dbSearch.pfind.PFindTask;
import javafx.scene.control.ProgressBar;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters.ToolType;
import uk.ac.quadram.metapilot.task.MetaAbstractTask;
import uk.ac.quadram.metapilot.task.par.MetaBaseParameter;
import uk.ac.quadram.metapilot.task.par.MetaConstants;

public class PFindSpectraConvertTask extends MetaAbstractTask {

	private String[] mgfFiles;
	private String[] pf2Files;

	// Progress monitoring
	private WatchService watchService;
	private Path spectraDirPath;
	private File spectraDirFile;
	private AtomicInteger completedFiles;
	private int totalFiles;
	private volatile boolean watchServiceRunning = false;

	private static final String taskName = "pFind spectrum convert workflow";

	public PFindSpectraConvertTask(MetaBaseParameter metaPar, MetaResourceParameters mrp, ProgressBar progressBar,
			Logger logger) {
		super(metaPar, mrp, progressBar, logger);
		// TODO Auto-generated constructor stub
	}

	@Override
	protected String getTaskName() {
		// TODO Auto-generated constructor stub
		return taskName;
	}

	@Override
	protected Logger getLogger() {
		// TODO Auto-generated constructor stub
		return logger;
	}

	@Override
	protected Boolean call() throws Exception {
		// TODO Auto-generated constructor stub

		updateMessage(getTaskName() + ": extracting spectra");

		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": extracting spectra started");
		logger.info(getTaskName() + ": extracting spectra started");

		// Log detailed information about the task
		logger.info("PFindSpectraConvertTask parameters:");
		logger.info("  Spectra directory: " + metaPar.getSpectraFile());
		logger.info("  MS2 scan mode: " + metaPar.getMs2ScanMode());
		logger.info("  Thread count: " + metaPar.getThreadCount());
		if (metaPar.getMetadata() != null) {
			String[] rawFiles = metaPar.getMetadata().getRawFiles();
			logger.info("  Raw files count: " + (rawFiles != null ? rawFiles.length : 0));
		}
		logger.info("  pFind path: " + mrp.getToolPath(MetaResourceParameters.ToolType.PFIND));

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

		this.spectraDirFile = metaPar.getSpectraFile();
		File rtFile = new File(spectraDirFile, "rt.txt");
		File ms2CountFile = new File(spectraDirFile, "ms2Count.txt");
		String[] rawFiles = metaPar.getMetadata().getRawFiles();
		if (MzmlSummaryFileWriter.allMzml(rawFiles)) {
			MzmlSummaryFileWriter.ensure(spectraDirFile, rawFiles);
		}

		String[] rawFileNames = new String[rawFiles.length];
		String[] ms2FileNames = new String[rawFiles.length];
		this.mgfFiles = new String[rawFiles.length];
		this.pf2Files = new String[rawFiles.length];

		// Initialize progress tracking
		totalFiles = rawFiles.length;
		completedFiles = new AtomicInteger(0);

		// Count existing MGF files and set initial progress
		int existingMgfCount = countExistingMgfFiles();
		completedFiles.set(existingMgfCount);

		// Set initial progress based on existing files
		if (totalFiles > 0) {
			double initialProgress = (double) existingMgfCount / totalFiles;
			updateProgress(initialProgress, 1.0);
			updateMessage(getTaskName() + ": Found " + existingMgfCount + "/" + totalFiles + " existing MGF files");
		}

		ArrayList<String> extractRawList = new ArrayList<String>();
		logger.info("Checking existing spectra files...");
		for (int i = 0; i < rawFiles.length; i++) {
			File rawi = new File(rawFiles[i]);
			String rawName = rawi.getName();

			rawFileNames[i] = rawName.substring(0, rawName.lastIndexOf("."));

			File pf2 = new File(spectraDirFile, rawFileNames[i] + suffix + ".pf2");
			File ms2 = new File(spectraDirFile, rawFileNames[i] + ".ms2");

			ms2FileNames[i] = ms2.getAbsolutePath();

			logger.info("  Raw file " + (i + 1) + ": " + rawName);
			logger.info("    Expected pf2 file: " + pf2.getAbsolutePath());
			logger.info("    Expected ms2 file: " + ms2.getAbsolutePath());
			logger.info("    pf2 exists: " + pf2.exists() + " (size: " + (pf2.exists() ? pf2.length() : 0) + ")");
			logger.info("    ms2 exists: " + ms2.exists() + " (size: " + (ms2.exists() ? ms2.length() : 0) + ")");

			if (!pf2.exists() || pf2.length() == 0) {
				extractRawList.add(rawFiles[i]);
				logger.info("    -> Will extract spectra for this file");
			} else {
				logger.info("    -> Spectra already exist, skipping extraction");
			}

			pf2Files[i] = pf2.getAbsolutePath();
		}

		logger.info("Files to extract: " + extractRawList.size());
		if (extractRawList.size() > 0) {
			logger.info("Starting spectra extraction for " + extractRawList.size() + " files...");
			for (int i = 0; i < extractRawList.size(); i++) {
				logger.info("  Extract file " + (i + 1) + ": " + extractRawList.get(i));
			}

			// Start WatchService to monitor MGF file generation
			startWatchServiceLoop();

			File pFindBin = (new File(mrp.getToolPath(ToolType.PFIND))).getParentFile();
			logger.info("pFind binary directory: " + pFindBin.getAbsolutePath());
			logger.info("pFind binary exists: " + pFindBin.exists());

			PFindTask pFindTask = new PFindTask(pFindBin, metaPar.getThreadCount());
			String[] extractRaws = extractRawList.toArray(new String[extractRawList.size()]);
			logger.info("Calling pFindTask.extract() with spectra directory: " + spectraDirFile.getAbsolutePath());
			pFindTask.extract(spectraDirFile.getAbsolutePath(), extractRaws);
			logger.info("Calling pFindTask.run() with " + extractRaws.length + " files");
			pFindTask.run(extractRaws.length);
			logger.info("pFindTask.run() completed");

			// Wait for WatchService to detect all files or timeout
			updateMessage(getTaskName() + ": Waiting for PF2 file generation to complete...");
			int timeoutSeconds = 300; // 5 minutes timeout
			int waitedSeconds = 0;
			while (completedFiles.get() < totalFiles && waitedSeconds < timeoutSeconds) {
				try {
					Thread.sleep(1000);
					waitedSeconds++;
					if (waitedSeconds % 10 == 0) {
						logger.info("Waiting for MGF files... (" + completedFiles.get() + "/" + totalFiles
								+ " completed, " + waitedSeconds + "s elapsed)");
					}
				} catch (InterruptedException e) {
					logger.warn("Wait interrupted");
					Thread.currentThread().interrupt();
					break;
				}
			}

			// Stop WatchService
			stopWatchService();

			if (completedFiles.get() >= totalFiles) {
				logger.info("All MGF files have been generated successfully");
				updateMessage(getTaskName() + ": All " + totalFiles + " MGF files generated successfully");
			} else {
				logger.warn("Timeout waiting for MGF files. Generated " + completedFiles.get() + "/" + totalFiles
						+ " files");
				updateMessage(getTaskName() + ": Timeout - Generated " + completedFiles.get() + "/" + totalFiles
						+ " MGF files");
			}

			// Check which extracted raws still have no valid pf2; retry those with
			// extractWithoutPParse2Plus
			ArrayList<String> unconvertedRawList = new ArrayList<String>();
			for (String rawPath : extractRawList) {
				File rawFile = new File(rawPath);
				String rawName = rawFile.getName();
				rawName = rawName.substring(0, rawName.lastIndexOf("."));
				File pf2 = new File(spectraDirFile, rawName + suffix + ".pf2");
				if (!pf2.exists() || pf2.length() == 0) {
					unconvertedRawList.add(rawPath);
					logger.info("No valid pf2 for " + rawPath + ", will retry with extractWithoutPParse2Plus");
				}
			}
			if (!unconvertedRawList.isEmpty()) {
				logger.info("Retrying " + unconvertedRawList.size()
						+ " unconverted file(s) with extractWithoutPParse2Plus");
				updateMessage(getTaskName() + ": Retrying " + unconvertedRawList.size()
						+ " file(s) with legacy extract (no pParse2Plus)...");
				startWatchServiceLoop();
				String[] unconvertedRaws = unconvertedRawList.toArray(new String[unconvertedRawList.size()]);
				pFindTask.extractWithoutPParse2Plus(spectraDirFile.getAbsolutePath(), unconvertedRaws);
				pFindTask.run(unconvertedRaws.length);
				waitedSeconds = 0;
				while (completedFiles.get() < totalFiles && waitedSeconds < timeoutSeconds) {
					try {
						Thread.sleep(1000);
						waitedSeconds++;
						if (waitedSeconds % 10 == 0) {
							logger.info("Waiting for PF2 after retry... (" + completedFiles.get() + "/" + totalFiles
									+ " completed, " + waitedSeconds + "s elapsed)");
						}
					} catch (InterruptedException e) {
						logger.warn("Wait interrupted");
						Thread.currentThread().interrupt();
						break;
					}
				}
				stopWatchService();
				if (completedFiles.get() >= totalFiles) {
					logger.info("All spectra files generated after retry with extractWithoutPParse2Plus");
				} else {
					logger.warn("Still missing some files after retry: " + completedFiles.get() + "/" + totalFiles);
				}
			}

			// Final check: if any spectra are still missing, log everything and return false
			ArrayList<String> stillMissingRawList = new ArrayList<String>();
			for (int i = 0; i < rawFiles.length; i++) {
				File pf2 = new File(spectraDirFile, rawFileNames[i] + suffix + ".pf2");
				if (!pf2.exists() || pf2.length() == 0) {
					stillMissingRawList.add(rawFiles[i]);
				}
			}
			if (!stillMissingRawList.isEmpty()) {
				logger.error(getTaskName() + ": " + stillMissingRawList.size() + " raw file(s) still have no valid PF2 spectra after extract and retry");
				logger.error("Spectra directory: " + spectraDirFile.getAbsolutePath());
				logger.error("Total raw files: " + rawFiles.length + ", missing: " + stillMissingRawList.size());
				logger.error("Raw file(s) that may have a problem:");
				for (String rawPath : stillMissingRawList) {
					File rawFile = new File(rawPath);
					String rn = rawFile.getName();
					rn = rn.substring(0, rn.lastIndexOf("."));
					File expectedPf2 = new File(spectraDirFile, rn + suffix + ".pf2");
					File expectedMs2 = new File(spectraDirFile, rn + ".ms2");
					logger.error("  - Raw: " + rawPath);
					logger.error("    Expected PF2: " + expectedPf2.getAbsolutePath() + " (exists: " + expectedPf2.exists() + ", size: " + (expectedPf2.exists() ? expectedPf2.length() : 0) + ")");
					logger.error("    Expected MS2:  " + expectedMs2.getAbsolutePath() + " (exists: " + expectedMs2.exists() + ", size: " + (expectedMs2.exists() ? expectedMs2.length() : 0) + ")");
				}
				System.err.println(format.format(new Date()) + "\t" + getTaskName() + ": SPECTRA CONVERSION FAILED - " + stillMissingRawList.size() + " raw file(s) have no valid PF2. See log for details.");
				System.err.println(getTaskName() + ": Raw file(s) that may have a problem:");
				for (String rawPath : stillMissingRawList) {
					System.err.println("  - " + rawPath);
				}
				updateMessage(getTaskName() + ": Failed - " + stillMissingRawList.size() + " raw file(s) have no valid spectra (see log)");
				return false;
			}
		} else {
			System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": spectra files have been found in "
					+ spectraDirFile);
			logger.info(getTaskName() + ": spectra files have been found in " + spectraDirFile);
			logger.info("All required spectra files already exist");
		}

		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": extracting spectra finished");
		logger.info(getTaskName() + ": extracting spectra finished");

		if (!rtFile.exists() || rtFile.length() == 0 || !ms2CountFile.exists() || ms2CountFile.length() == 0) {

			System.out.println(
					format.format(new Date()) + "\t" + getTaskName() + ": writing retention time information started");
			logger.info(getTaskName() + ": writing retention time information started");

			PrintWriter rtWriter = new PrintWriter(rtFile);
			PrintWriter ms2CountWriter = new PrintWriter(ms2CountFile);

			rtWriter.println("FileName\tScannum\tRt");
			ms2CountWriter.println("FileName\tMS2 count");

			DecimalFormat df = null;

			for (int i = 0; i < rawFiles.length; i++) {
				File rawi = new File(rawFiles[i]);
				String rawName = rawi.getName();
				rawName = rawName.substring(0, rawName.length() - ".raw".length());

				File ms2File = new File(spectraDirFile, rawName + ".ms2");
				File mgfFile = new File(spectraDirFile, rawName + ".mgf");
				mgfFiles[i] = mgfFile.getAbsolutePath();

				PrintWriter writer = new PrintWriter(mgfFile);
				int scannum = 0;
				double premz = 0.0;
				double rtsecond = 0.0;
				boolean begin = false;
				int count = 0;

				BufferedReader reader = new BufferedReader(new FileReader(ms2File));
				String line = null;
				while ((line = reader.readLine()) != null) {
					String[] cs = line.split("\t");
					if (cs.length == 1) {
						writer.println(line.replaceAll(",", "."));
					} else {

						if (cs[0].equals("S")) {
							if (begin) {
								writer.println("END IONS");
								begin = false;
							}
							writer.println("BEGIN IONS");

							if (df == null) {
								df = new DecimalFormat();
								df.applyPattern("#,##0.00");
								if (cs[3].contains(",")) {
									df.setDecimalFormatSymbols(new DecimalFormatSymbols(Locale.GERMANY));
								} else {
									df.setDecimalFormatSymbols(new DecimalFormatSymbols(Locale.US));
								}
							}

							scannum = Integer.parseInt(cs[1]);
							try {
								premz = df.parse(cs[3]).doubleValue();
							} catch (ParseException e) {
								// TODO Auto-generated catch block
								System.err.println(format.format(new Date()) + "\t" + getTaskName()
										+ ": error in parsing the number from " + line);
								logger.error(getTaskName() + ": error in parsing the number from " + line, e);
								reader.close();
								writer.close();
								rtWriter.close();
								return false;
							}
							count++;

						} else if (cs[0].equals("Z")) {
							writer.println("TITLE=" + rawName + "." + scannum + "." + scannum + "." + cs[1]);
							writer.println("RTINSECONDS=" + rtsecond);
							writer.println("PEPMASS=" + premz);
							writer.println("CHARGE=" + cs[1] + "+");
							begin = true;
						} else if (cs[0].equals("I")) {
							if (cs[1].equals("RetTime")) {
								try {
									if (df != null) {
										rtsecond = df.parse(cs[2]).doubleValue();
									}
								} catch (ParseException e) {
									// TODO Auto-generated catch block
									System.err.println(format.format(new Date()) + "\t" + getTaskName()
											+ ": error in parsing the number from " + line);
									logger.error(getTaskName() + ": error in parsing the number from " + line, e);
									reader.close();
									writer.close();
									rtWriter.close();
									return false;
								}

								rtWriter.println(rawName + "\t" + scannum + "\t" + rtsecond / 60.0);
							}
						}
					}
				}
				reader.close();

				if (begin) {
					writer.println("END IONS");
				}
				writer.close();

				ms2CountWriter.println(rawName + "\t" + count);
			}

			rtWriter.close();
			ms2CountWriter.close();

			System.out.println(
					format.format(new Date()) + "\t" + getTaskName() + ": writing retention time information finished");
			logger.info(getTaskName() + ": writing retention time information finished");
		} else {
			DecimalFormat format = null;
			for (int i = 0; i < rawFiles.length; i++) {
				File rawi = new File(rawFiles[i]);
				String rawName = rawi.getName();
				rawName = rawName.substring(0, rawName.length() - ".raw".length());

				File ms2File = new File(spectraDirFile, rawName + ".ms2");
				File mgfFile = new File(spectraDirFile, rawName + ".mgf");
				mgfFiles[i] = mgfFile.getAbsolutePath();

				if (!mgfFile.exists()) {
					PrintWriter writer = new PrintWriter(mgfFile);
					int scannum = 0;
					double premz = 0.0;
					double rtsecond = 0.0;
					boolean begin = false;

					BufferedReader reader = new BufferedReader(new FileReader(ms2File));
					String line = null;
					while ((line = reader.readLine()) != null) {
						String[] cs = line.split("\t");
						if (cs.length == 1) {
							writer.println(line.replaceAll(",", "."));
						} else {

							if (cs[0].equals("S")) {
								if (begin) {
									writer.println("END IONS");
									begin = false;
								}
								writer.println("BEGIN IONS");

								if (format == null) {
									format = new DecimalFormat();
									format.applyPattern("#,##0.00");
									if (cs[3].contains(",")) {
										format.setDecimalFormatSymbols(new DecimalFormatSymbols(Locale.GERMANY));
									} else {
										format.setDecimalFormatSymbols(new DecimalFormatSymbols(Locale.US));
									}
								}

								scannum = Integer.parseInt(cs[1]);
								try {
									premz = format.parse(cs[3]).doubleValue();
								} catch (ParseException e) {
									// TODO Auto-generated catch block
									System.err.println(format.format(new Date()) + "\t" + getTaskName()
											+ ": error in parsing the number from " + line);
									logger.error(getTaskName() + ": error in parsing the number from " + line, e);
									reader.close();
									writer.close();
									return false;
								}
							} else if (cs[0].equals("Z")) {
								writer.println("TITLE=" + rawName + "." + scannum + "." + scannum + "." + cs[1]);
								writer.println("RTINSECONDS=" + rtsecond);
								writer.println("PEPMASS=" + premz);
								writer.println("CHARGE=" + cs[1] + "+");
								begin = true;
							} else if (cs[0].equals("I")) {
								if (cs[1].equals("RetTime")) {
									try {
										if (format != null) {
											rtsecond = format.parse(cs[2]).doubleValue();
										}
									} catch (ParseException e) {
										// TODO Auto-generated catch block
										System.err.println(format.format(new Date()) + "\t" + getTaskName()
												+ ": error in parsing the number from " + line);
										logger.error(getTaskName() + ": error in parsing the number from " + line, e);
										reader.close();
										writer.close();
										return false;
									}
								}
							}
						}
					}
					reader.close();

					if (begin) {
						writer.println("END IONS");
					}
					writer.close();
				}
			}
		}

		// Ensure WatchService is stopped
		stopWatchService();

		// Stop the progress bar on successful completion
		updateProgress(100, 100);
		updateMessage(getTaskName() + ": Spectra conversion completed successfully!");

		return true;
	}

	public String[] getMgfFiles() {
		return mgfFiles;
	}

	public String[] getPf2Files() {
		return pf2Files;
	}

	/**
	 * Starts a WatchService to monitor the spectra directory for .mgf file
	 * generation and updates progress accordingly.
	 */
	private void startWatchServiceLoop() {
		try {
			watchService = FileSystems.getDefault().newWatchService();
			spectraDirPath = spectraDirFile.toPath();
			spectraDirPath.register(watchService, StandardWatchEventKinds.ENTRY_CREATE,
					StandardWatchEventKinds.ENTRY_MODIFY);

			watchServiceRunning = true;
			logger.info("WatchService started for directory: " + spectraDirPath);

			// Start the watch service in a separate thread
			Thread watchThread = new Thread(() -> {
				try {
					while (watchServiceRunning && !Thread.currentThread().isInterrupted()) {
						WatchKey key = watchService.take();

						for (WatchEvent<?> event : key.pollEvents()) {
							WatchEvent.Kind<?> kind = event.kind();

							if (kind == StandardWatchEventKinds.OVERFLOW) {
								continue;
							}

							@SuppressWarnings("unchecked")
							WatchEvent<Path> ev = (WatchEvent<Path>) event;
							Path filename = ev.context();

							// Check if the created/modified file is a .mgf file
							if (filename.toString().toLowerCase().endsWith(".mgf")
									|| filename.toString().toLowerCase().endsWith(".ms2")) {

								Path fullPath = spectraDirPath.resolve(filename);
								File mgfFile = fullPath.toFile();

								// Wait a bit for the file to be fully written
								Thread.sleep(1000);

								// Check if file exists and has content
								if (mgfFile.exists() && mgfFile.length() > 0) {
									int currentCompleted = completedFiles.incrementAndGet();
									updateProgress((double) currentCompleted, (double) totalFiles * 2.0);
									updateMessage(getTaskName() + ": Generated " + currentCompleted + "/"
											+ totalFiles * 2 + " MGF/MS2 files");

									logger.info("MGF/MS2 file detected: " + filename + " (Progress: " + currentCompleted
											+ "/" + totalFiles * 2 + ")");

									// Check if all files are completed
									if (currentCompleted >= totalFiles * 2.0) {
										logger.info("All MGF files have been generated. Stopping WatchService.");
										stopWatchService();
										break;
									}
								}
							}
						}

						boolean valid = key.reset();
						if (!valid) {
							logger.warn("WatchKey is no longer valid");
							break;
						}
					}
				} catch (InterruptedException e) {
					logger.info("WatchService interrupted");
					Thread.currentThread().interrupt();
				} catch (Exception e) {
					logger.error("Error in WatchService loop", e);
				}
			});

			watchThread.setDaemon(true);
			watchThread.setName("PFindSpectraConvert-WatchService");
			watchThread.start();

		} catch (Exception e) {
			logger.error("Failed to start WatchService", e);
		}
	}

	/**
	 * Stops the WatchService loop
	 */
	private void stopWatchService() {
		watchServiceRunning = false;
		if (watchService != null) {
			try {
				watchService.close();
				logger.info("WatchService stopped");
			} catch (Exception e) {
				logger.error("Error stopping WatchService", e);
			}
		}
	}

	/**
	 * Counts existing .mgf files in the spectra directory and sets initial progress
	 */
	private int countExistingMgfFiles() {
		int existingCount = 0;
		String[] rawFiles = metaPar.getMetadata().getRawFiles();

		for (int i = 0; i < rawFiles.length; i++) {
			File rawi = new File(rawFiles[i]);
			String rawName = rawi.getName();
			rawName = rawName.substring(0, rawName.lastIndexOf("."));

			File mgfFile = new File(spectraDirFile, rawName + ".mgf");
			if (mgfFile.exists() && mgfFile.length() > 0) {
				existingCount++;
				logger.info("Existing MGF file found: " + mgfFile.getName());
			}

			File ms2File = new File(spectraDirFile, rawName + ".ms2");
			if (ms2File.exists() && ms2File.length() > 0) {
				existingCount++;
				logger.info("Existing MS2 file found: " + ms2File.getName());
			}
		}

		logger.info("Found " + existingCount + " existing MGF/MS2 files");
		return existingCount;
	}

}
