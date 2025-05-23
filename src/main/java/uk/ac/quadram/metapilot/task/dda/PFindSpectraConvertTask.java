package uk.ac.quadram.metapilot.task.dda;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Date;
import java.util.Locale;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import uk.ac.quadram.metapilot.dbSearch.pfind.PFindTask;
import uk.ac.quadram.metapilot.gui.JavaFXProgressBar;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters.ToolType;
import uk.ac.quadram.metapilot.task.MetaAbstractTask;
import uk.ac.quadram.metapilot.task.par.MetaBaseParameter;
import uk.ac.quadram.metapilot.task.par.MetaConstants;

public class PFindSpectraConvertTask extends MetaAbstractTask {
	
	private String[] mgfFiles;
	private String[] pf2Files;

	private static final Logger LOGGER = LogManager.getLogger(PFindSpectraConvertTask.class);
	private static final String taskName = "pFind spectrum convert workflow";

	public PFindSpectraConvertTask(MetaBaseParameter metaPar, MetaResourceParameters mrp,
			JavaFXProgressBar progressBar) {
		super(metaPar, mrp, progressBar);
		// TODO Auto-generated constructor stub
	}

	@Override
	protected void initial() {
		// TODO Auto-generated method stub

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

	@Override
	protected Boolean call() throws Exception {
		// TODO Auto-generated method stub

		progressBar.setString(getTaskName() + ": extracting spectra");

		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": extracting spectra started");
		LOGGER.info(getTaskName() + ": extracting spectra started");

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

		File spectraDirFile = metaPar.getSpectraFile();
		File rtFile = new File(spectraDirFile, "rt.txt");
		File ms2CountFile = new File(spectraDirFile, "ms2Count.txt");
		String[] rawFiles = metaPar.getMetadata().getRawFiles();

		String[] rawFileNames = new String[rawFiles.length];
		String[] ms2FileNames = new String[rawFiles.length];
		this.mgfFiles = new String[rawFiles.length];
		this.pf2Files = new String[rawFiles.length];

		ArrayList<String> extractRawList = new ArrayList<String>();
		for (int i = 0; i < rawFiles.length; i++) {
			File rawi = new File(rawFiles[i]);
			String rawName = rawi.getName();

			rawFileNames[i] = rawName.substring(0, rawName.lastIndexOf("."));

			File pf2 = new File(spectraDirFile, rawFileNames[i] + suffix + ".pf2");
			File ms2 = new File(spectraDirFile, rawFileNames[i] + ".ms2");

			ms2FileNames[i] = ms2.getAbsolutePath();

			if (!pf2.exists() || pf2.length() == 0) {
				extractRawList.add(rawFiles[i]);
			}
			
			pf2Files[i] = pf2.getAbsolutePath();
		}

		File pFindBin = (new File(mrp.getToolPath(ToolType.PFIND))).getParentFile();
		PFindTask pFindTask = new PFindTask(pFindBin, metaPar.getThreadCount());
		if (extractRawList.size() > 0) {
			String[] extractRaws = extractRawList.toArray(new String[extractRawList.size()]);
			pFindTask.extract(spectraDirFile.getAbsolutePath(), extractRaws);
			pFindTask.run(extractRaws.length);
		} else {
			System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": spectra files have been found in "
					+ spectraDirFile);
			LOGGER.info(getTaskName() + ": spectra files have been found in " + spectraDirFile);
			return false;
		}

		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": extracting spectra finished");
		LOGGER.info(getTaskName() + ": extracting spectra finished");

		if (!rtFile.exists() || rtFile.length() == 0 || !ms2CountFile.exists() || ms2CountFile.length() == 0) {

			System.out.println(
					format.format(new Date()) + "\t" + getTaskName() + ": writing retention time information started");
			LOGGER.info(getTaskName() + ": writing retention time information started");

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
				File mgfFile = new File(spectraDirFile, rawName + suffix + ".mgf");
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
								LOGGER.error(getTaskName() + ": error in parsing the number from " + line, e);
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
									LOGGER.error(getTaskName() + ": error in parsing the number from " + line, e);
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
			LOGGER.info(getTaskName() + ": writing retention time information finished");
		} else {
			DecimalFormat format = null;
			for (int i = 0; i < rawFiles.length; i++) {
				File rawi = new File(rawFiles[i]);
				String rawName = rawi.getName();
				rawName = rawName.substring(0, rawName.length() - ".raw".length());

				File ms2File = new File(spectraDirFile, rawName + ".ms2");
				File mgfFile = new File(spectraDirFile, rawName + suffix + ".mgf");
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
									LOGGER.error(getTaskName() + ": error in parsing the number from " + line, e);
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
										LOGGER.error(getTaskName() + ": error in parsing the number from " + line, e);
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
		return true;
	}

	public String[] getMgfFiles() {
		return mgfFiles;
	}

	public String[] getPf2Files() {
		return pf2Files;
	}

	
}
