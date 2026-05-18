package uk.ac.quadram.metapilot.dbSearch.diann;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Locale;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import uk.ac.quadram.metapilot.quant.flashLFQ.FlashLfqQuanProReader;
import uk.ac.quadram.metapilot.quant.flashLFQ.FlashLfqQuanProtein;

public class DiannPgMatrixReader extends FlashLfqQuanProReader {

	private boolean isCsv;
	private static Logger LOGGER = LogManager.getLogger(DiannPgMatrixReader.class);

	public DiannPgMatrixReader(String in) {
		this(new File(in));
		// TODO Auto-generated constructor stub
	}

	public DiannPgMatrixReader(File in) {
		super(in);
		// TODO Auto-generated constructor stub
	}

	protected void parseTitle() {
		// TODO Auto-generated method stub
		intensityTitleIdMap = new HashMap<String, Integer>();
		ArrayList<String> intensityTitleList = new ArrayList<String>();

		String line = null;
		try {
			line = reader.readLine();
			String fileName = this.getFile().getName();
			if (fileName.endsWith("tsv") || fileName.endsWith("TSV")) {
				this.title = line.split("\t");
				this.isCsv = false;
			} else if (fileName.endsWith("csv") || fileName.endsWith("CSV")) {
				this.title = line.split(",");
				this.isCsv = true;
			} else {

			}

			for (int i = 0; i < title.length; i++) {
				if (title[i].equals("Protein.Group") || title[i].equals("ProteinGroup")) {
					proteinId = i;
				} else if (isIntensityColumn(title[i])) {
					String runName = normalizeRunName(title[i]);
					intensityTitleIdMap.put(runName, i);
					intensityTitleList.add(runName);
				}
			}

			this.intensityTitles = intensityTitleList.toArray(new String[intensityTitleList.size()]);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			LOGGER.error("Error in reading FlashLFQ quantification result file " + super.getFile(), e);
		}
	}

	protected FlashLfqQuanProtein parse() {
		if (line.trim().length() == 0) {
			return null;
		}

		String[] cs;
		if (this.isCsv) {
			cs = line.split(",");
		} else {
			cs = line.split("\t");
		}

		double[] intensity = new double[intensityTitles.length];

		for (int i = 0; i < intensity.length; i++) {
			int id = this.intensityTitleIdMap.get(intensityTitles[i]);
			if (id > -1 && id < cs.length) {
				if (cs[id].trim().length() > 0) {
					intensity[i] = Double.parseDouble(cs[id]);
				} else {
					intensity[i] = 0;
				}

			} else {
				intensity[i] = 0;
			}
		}

		String proname = cs[proteinId].split(";")[0];
		FlashLfqQuanProtein protein = new FlashLfqQuanProtein(proname, intensity);

		return protein;
	}

	private static boolean isIntensityColumn(String title) {
		if (title == null) {
			return false;
		}
		String lower = title.toLowerCase(Locale.ROOT);
		return lower.endsWith(".d") || lower.endsWith(".dia") || lower.endsWith(".raw")
				|| lower.endsWith(".mzml");
	}

	private static String normalizeRunName(String runName) {
		if (runName == null) {
			return "";
		}
		String name = runName.trim();
		int slash = Math.max(name.lastIndexOf('\\'), name.lastIndexOf('/'));
		if (slash >= 0 && slash < name.length() - 1) {
			name = name.substring(slash + 1);
		}
		String lower = name.toLowerCase(Locale.ROOT);
		if (lower.endsWith(".raw.dia")) {
			return name.substring(0, name.length() - 8);
		}
		if (lower.endsWith(".mzml.dia")) {
			return name.substring(0, name.length() - 9);
		}
		if (lower.endsWith(".d.dia")) {
			return name.substring(0, name.length() - 6);
		}
		if (lower.endsWith(".mzml")) {
			return name.substring(0, name.length() - 5);
		}
		if (lower.endsWith(".raw")) {
			return name.substring(0, name.length() - 4);
		}
		if (lower.endsWith(".dia")) {
			return name.substring(0, name.length() - 4);
		}
		if (lower.endsWith(".d")) {
			return name.substring(0, name.length() - 2);
		}
		return name;
	}

}
