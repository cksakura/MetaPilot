package uk.ac.quadram.metapilot.core.runner;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Writes rt.txt and ms2Count.txt from mzML inputs for summary report generation.
 */
public final class MzmlSummaryFileWriter {

	private static final Pattern SCAN_PATTERN = Pattern.compile("scan=(\\d+)");
	private static final Pattern RT_PATTERN = Pattern.compile("name=\"scan start time\"\\s+value=\"([^\"]+)\"");
	private static final Pattern MS_LEVEL_PATTERN = Pattern.compile("name=\"ms level\"\\s+value=\"(\\d+)\"");

	private MzmlSummaryFileWriter() {
	}

	public static boolean allMzml(String[] files) {
		if (files == null || files.length == 0) {
			return false;
		}
		for (String file : files) {
			if (file == null || !file.toLowerCase().endsWith(".mzml")) {
				return false;
			}
		}
		return true;
	}

	public static boolean ensure(File spectraDir, String[] mzmlFiles) throws IOException {
		if (spectraDir == null || mzmlFiles == null || mzmlFiles.length == 0) {
			return false;
		}
		File rtFile = new File(spectraDir, "rt.txt");
		File ms2CountFile = new File(spectraDir, "ms2Count.txt");
		if (rtFile.exists() && rtFile.length() > 0 && ms2CountFile.exists() && ms2CountFile.length() > 0) {
			return true;
		}
		if (!spectraDir.exists() && !spectraDir.mkdirs()) {
			throw new IOException("Cannot create spectra directory: " + spectraDir.getAbsolutePath());
		}

		try (PrintWriter rtWriter = new PrintWriter(rtFile);
				PrintWriter ms2CountWriter = new PrintWriter(ms2CountFile)) {
			rtWriter.println("FileName\tScannum\tRt");
			ms2CountWriter.println("FileName\tMS2 count");
			for (String mzml : mzmlFiles) {
				File mzmlFile = new File(mzml);
				if (!mzmlFile.exists() || mzmlFile.length() == 0) {
					return false;
				}
				int count = writeOneMzml(rtWriter, mzmlFile);
				ms2CountWriter.println(baseNameNoExt(mzmlFile) + "\t" + count);
			}
		}
		return true;
	}

	private static int writeOneMzml(PrintWriter rtWriter, File mzmlFile) throws IOException {
		String fileName = baseNameNoExt(mzmlFile);
		int count = 0;
		int scannum = 0;
		boolean isMS2 = false;

		try (BufferedReader reader = new BufferedReader(new FileReader(mzmlFile))) {
			String line;
			while ((line = reader.readLine()) != null) {
				line = line.trim();
				if (line.startsWith("<spectrum index=")) {
					Matcher matcher = SCAN_PATTERN.matcher(line);
					if (matcher.find()) {
						scannum = Integer.parseInt(matcher.group(1));
					}
					isMS2 = false;
				} else if (!isMS2 && (isMsnSpectrum(line) || isMsLevel2(line))) {
					count++;
					isMS2 = true;
				} else if (line.startsWith("<cvParam") && line.contains("name=\"scan start time\"")) {
					Matcher matcher = RT_PATTERN.matcher(line);
					if (matcher.find() && isMS2) {
						rtWriter.println(fileName + "\t" + scannum + "\t" + matcher.group(1));
					}
				}
			}
		}
		return count;
	}

	private static boolean isMsnSpectrum(String line) {
		return line.contains("name=\"MSn spectrum\"");
	}

	private static boolean isMsLevel2(String line) {
		Matcher matcher = MS_LEVEL_PATTERN.matcher(line);
		return matcher.find() && "2".equals(matcher.group(1));
	}

	private static String baseNameNoExt(File file) {
		String name = file.getName();
		String lower = name.toLowerCase();
		int dot = lower.lastIndexOf(".mzml");
		if (dot > 0) {
			return name.substring(0, dot);
		}
		dot = name.lastIndexOf('.');
		return dot > 0 ? name.substring(0, dot) : name;
	}
}
