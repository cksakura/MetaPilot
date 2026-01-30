package uk.ac.quadram.metapilot.task.dda;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Date;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import org.apache.logging.log4j.Logger;

import javafx.scene.control.ProgressBar;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters.ToolType;
import uk.ac.quadram.metapilot.spectra.Spectrum;
import uk.ac.quadram.metapilot.spectra.io.MSConvertor;
import uk.ac.quadram.metapilot.spectra.io.MgfReader;
import uk.ac.quadram.metapilot.task.MetaAbstractTask;
import uk.ac.quadram.metapilot.task.par.MetaBaseParameter;

public class MSConvertTask extends MetaAbstractTask {

	private static final String taskName = "pFind spectrum convert workflow";

	public MSConvertTask(MetaBaseParameter metaPar, MetaResourceParameters mrp, ProgressBar progressBar,
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
	protected Boolean call() throws Exception {
		// TODO Auto-generated constructor stub

		updateMessage(getTaskName() + ": extracting spectra");

		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": extracting spectra started");
		logger.info(getTaskName() + ": extracting spectra started");

		File spectraDirFile = metaPar.getSpectraFile();
		File rtFile = new File(spectraDirFile, "rt.txt");
		File ms2CountFile = new File(spectraDirFile, "ms2Count.txt");
		String[] rawFiles = metaPar.getMetadata().getRawFiles();

		String[] rawFileNames = new String[rawFiles.length];
		String[] mgfFileNames = new String[rawFiles.length];

		ArrayList<String> extractRawList = new ArrayList<String>();
		for (int i = 0; i < rawFiles.length; i++) {

			File rawi = new File(rawFiles[i]);
			String rawName = rawi.getName();

			rawFileNames[i] = rawName.substring(0, rawName.lastIndexOf("."));

			File mgf = new File(spectraDirFile, rawFileNames[i] + ".mgf");

			mgfFileNames[i] = mgf.getAbsolutePath();

			if (!mgf.exists() || mgf.length() == 0) {
				extractRawList.add(rawFiles[i]);
			}
		}

		if (extractRawList.size() > 0) {

			File msconvertFile = new File(mrp.getToolPath(ToolType.MSCONVERT));
			if (msconvertFile.exists()) {

				ExecutorService es = Executors.newFixedThreadPool(metaPar.getThreadCount());

				for (int i = 0; i < extractRawList.size(); i++) {

					String rawi = extractRawList.get(i);

					Runnable runable = new Runnable() {
						@Override
						public void run() {
							// TODO Auto-generated method stub

							if (rawi.endsWith(".mzml") || rawi.endsWith(".mzML")) {
								MSConvertor.rawMzml2mgf(msconvertFile.getAbsolutePath(), rawi,
										spectraDirFile.getAbsolutePath());
							} else if (rawi.endsWith(".d") || rawi.endsWith(".D")) {
								MSConvertor.d2mgf(msconvertFile.getAbsolutePath(), rawi,
										spectraDirFile.getAbsolutePath());
							} else {
								System.out.println(
										format.format(new Date()) + "\t" + getTaskName() + ": unknown spectra format");
								logger.info(getTaskName() + ": unknown spectra format");
							}
						}
					};
					es.submit(runable);
				}

				es.shutdown();

				try {
					boolean finish = es.awaitTermination(extractRawList.size() * 4, TimeUnit.HOURS);
					if (finish) {
						System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": msconvert finished");
						logger.info(getTaskName() + ": msconvert finished");
					} else {
						System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": msconvert failed");
						logger.info(getTaskName() + ": msconvert failed");
					}

				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": msconvert failed");
					logger.error(getTaskName() + ": msconvert failed", e);
				}

			} else {
				System.out.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": msconvert was not found, please " + "put \"ProteoWizard\" into the resource folder");
				logger.info(getTaskName() + ": msconvert was not found, please "
						+ "put \"ProteoWizard\" into the resource folder");
			}

		} else {
			System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": spectra files have been found in "
					+ spectraDirFile);
			logger.info(getTaskName() + ": spectra files have been found in " + spectraDirFile);
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

			for (int i = 0; i < mgfFileNames.length; i++) {

				int count = 0;
				MgfReader reader = new MgfReader(mgfFileNames[i]);
				Spectrum sp = null;
				while ((sp = reader.getNextSpectrum()) != null) {
					count++;
					int scannum = sp.getScannum();
					double rtsecond = sp.getRt();
					rtWriter.println(rawFileNames[i] + "\t" + scannum + "\t" + rtsecond);
				}
				reader.close();

				ms2CountWriter.println(rawFileNames[i] + "\t" + count);
			}

			rtWriter.close();
			ms2CountWriter.close();

			System.out.println(
					format.format(new Date()) + "\t" + getTaskName() + ": writing retention time information finished");
			logger.info(getTaskName() + ": writing retention time information finished");
		}

		return true;
	}

}
