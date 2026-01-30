package uk.ac.quadram.metapilot.quant.flashLFQ;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import uk.ac.quadram.metapilot.dbSearch.pfind.PFindResultReader;
import uk.ac.quadram.metapilot.task.io.pep.AbstractMetaPeptideReader;
import uk.ac.quadram.metapilot.task.par.MetaConstants;

public class FlashLfqQuanPepReader extends AbstractMetaPeptideReader {

	protected BufferedReader reader;
	protected static int sequenceId = -1;
	protected static int modSeqId = -1;
	protected static int proteinId = -1;

	protected String[] intensityTitles;
	protected HashMap<String, Integer> intensityTitleIdMap;
	protected String[] idenTypeTitles;
	protected HashMap<String, Integer> idenTypeTitleIdMap;
	protected String[] title;
	protected String line;

	private static Logger LOGGER = LogManager.getLogger(FlashLfqQuanPepReader.class);

	public FlashLfqQuanPepReader(String quanPep) {
		this(new File(quanPep));
	}

	public FlashLfqQuanPepReader(File quanPep) {
		super(quanPep);
		try {
			this.reader = new BufferedReader(new FileReader(quanPep));
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			LOGGER.error("Error in reading FlashLFQ search result file " + quanPep, e);
		}
		this.parseTitle();
	}
	
	public FlashLfqQuanPepReader(String quanPep, ArrayList<String> fileList) {
		this(new File(quanPep), fileList);
	}

	public FlashLfqQuanPepReader(File quanPep, ArrayList<String> fileList) {
		super(quanPep);
		try {
			this.reader = new BufferedReader(new FileReader(quanPep));
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			LOGGER.error("Error in reading FlashLFQ search result file " + quanPep, e);
		}
		this.parseTitle(fileList);
	}
	
	public FlashLfqQuanPepReader(String quanPep, ArrayList<String> fileList, String[] reportTitle) {
		this(new File(quanPep), fileList, reportTitle);
	}

	public FlashLfqQuanPepReader(File quanPep, ArrayList<String> fileList, String[] reportTitle) {
		super(quanPep);
		try {
			this.reader = new BufferedReader(new FileReader(quanPep));
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			LOGGER.error("Error in reading FlashLFQ search result file " + quanPep, e);
		}
		this.parseTitle(fileList, reportTitle);
	}

	protected void parseTitle() {
		// TODO Auto-generated method stub
		intensityTitleIdMap = new HashMap<String, Integer>();
		ArrayList<String> intensityTitleList = new ArrayList<String>();
		idenTypeTitleIdMap = new HashMap<String, Integer>();
		ArrayList<String> idenTypeList = new ArrayList<String>();

		String line = null;
		try {
			line = reader.readLine();
			this.title = line.split("\t");
			LOGGER.info("FlashLfqQuanPepReader: header has {} columns; Intensity_ columns: {}; Detection Type_ columns: {}",
					title.length,
					java.util.Arrays.stream(title).filter(t -> t != null && t.startsWith("Intensity_")).toList(),
					java.util.Arrays.stream(title).filter(t -> t != null && t.startsWith("Detection Type_")).toList());
			for (int i = 0; i < title.length; i++) {
				if (title[i].equals("Sequence")) {
					modSeqId = i;
				} else if (title[i].equals("Base Sequence")) {
					sequenceId = i;
				} else if (title[i].equals("Protein Groups")) {
					proteinId = i;
				} else if (title[i].startsWith("Intensity_")) {
					intensityTitleIdMap.put(title[i], i);
					intensityTitleList.add(title[i]);
				} else if (title[i].startsWith("Detection Type_")) {
					idenTypeTitleIdMap.put(title[i], i);
					idenTypeList.add(title[i]);
				}
			}

			if (sequenceId == -1) {
				sequenceId = modSeqId;
			}

			this.intensityTitles = intensityTitleList.toArray(String[]::new);
			this.idenTypeTitles = idenTypeList.toArray(String[]::new);
			LOGGER.info("FlashLfqQuanPepReader: parsed {} intensity columns, {} idenType columns (parseTitle no fileList)",
					intensityTitles.length, idenTypeTitles.length);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			LOGGER.error("Error in reading FlashLFQ quantification result file " + super.getFile(), e);
		}
	}
	
	protected void parseTitle(ArrayList<String> fileList) {
		// TODO Auto-generated method stub
		intensityTitleIdMap = new HashMap<String, Integer>();
		idenTypeTitleIdMap = new HashMap<String, Integer>();
		this.intensityTitles = new String[fileList.size()];
		this.idenTypeTitles = new String[fileList.size()];

		String line = null;
		try {
			line = reader.readLine();
			this.title = line.split("\t");
			LOGGER.info("FlashLfqQuanPepReader: fileList size={}, names={}", fileList.size(), fileList);
			java.util.List<String> intensityCols = new java.util.ArrayList<>();
			java.util.List<String> detectionTypeCols = new java.util.ArrayList<>();
			for (int i = 0; i < title.length; i++) {
				if (title[i] != null && title[i].startsWith("Intensity_")) {
					intensityCols.add(title[i]);
				}
				if (title[i] != null && title[i].startsWith("Detection Type_")) {
					detectionTypeCols.add(title[i]);
				}
			}
			LOGGER.info("FlashLfqQuanPepReader: header Intensity_ columns ({}): {}", intensityCols.size(), intensityCols);
			LOGGER.info("FlashLfqQuanPepReader: header Detection Type_ columns ({}): {}", detectionTypeCols.size(), detectionTypeCols);

			for (int i = 0; i < title.length; i++) {
				if (title[i].equals("Sequence")) {
					modSeqId = i;
				} else if (title[i].equals("Base Sequence")) {
					sequenceId = i;
				} else if (title[i].equals("Protein Groups")) {
					proteinId = i;
				} else if (title[i].startsWith("Intensity_")) {
					intensityTitleIdMap.put(title[i], i);
					for (int j = 0; j < fileList.size(); j++) {
						if (title[i].equals("Intensity_" + fileList.get(j))) {
							this.intensityTitles[j] = title[i];
						}
					}
					// Fallback: match by file name only (header may have path: Intensity_/path/to/file.raw)
					String intensitySuffix = title[i].substring("Intensity_".length());
					for (int j = 0; j < fileList.size(); j++) {
						if (this.intensityTitles[j] == null
								&& (intensitySuffix.equals(fileList.get(j))
										|| intensitySuffix.replace('\\', '/').endsWith("/" + fileList.get(j)))) {
							this.intensityTitles[j] = title[i];
							LOGGER.info("FlashLfqQuanPepReader: matched intensity column '{}' to fileList[{}]='{}' by path fallback",
									title[i], j, fileList.get(j));
							break;
						}
					}
				} else if (title[i].startsWith("Detection Type_")) {
					idenTypeTitleIdMap.put(title[i], i);
					for (int j = 0; j < fileList.size(); j++) {
						if (title[i].equals("Detection Type_" + fileList.get(j))) {
							this.idenTypeTitles[j] = title[i];
						}
					}
					// Fallback: match by file name only
					String detectionSuffix = title[i].substring("Detection Type_".length());
					for (int j = 0; j < fileList.size(); j++) {
						if (this.idenTypeTitles[j] == null
								&& (detectionSuffix.equals(fileList.get(j))
										|| detectionSuffix.replace('\\', '/').endsWith("/" + fileList.get(j)))) {
							this.idenTypeTitles[j] = title[i];
							LOGGER.info("FlashLfqQuanPepReader: matched Detection Type column '{}' to fileList[{}]='{}' by path fallback",
									title[i], j, fileList.get(j));
							break;
						}
					}
				}
			}

			if (sequenceId == -1) {
				sequenceId = modSeqId;
			}

			int matchedIntensity = 0;
			int matchedIdenType = 0;
			for (int j = 0; j < fileList.size(); j++) {
				if (this.intensityTitles[j] != null) matchedIntensity++;
				if (this.idenTypeTitles[j] != null) matchedIdenType++;
			}
			LOGGER.info("FlashLfqQuanPepReader: after parseTitle(fileList): matched {}/{} intensity columns, {}/{} Detection Type columns",
					matchedIntensity, fileList.size(), matchedIdenType, fileList.size());
			if (matchedIntensity < fileList.size() || matchedIdenType < fileList.size()) {
				for (int j = 0; j < fileList.size(); j++) {
					if (this.intensityTitles[j] == null) {
						LOGGER.warn("FlashLfqQuanPepReader: no header match for Intensity_ fileList[{}]='{}'; expected column 'Intensity_{}'",
								j, fileList.get(j), fileList.get(j));
					}
					if (this.idenTypeTitles[j] == null) {
						LOGGER.warn("FlashLfqQuanPepReader: no header match for Detection Type_ fileList[{}]='{}'; expected 'Detection Type_{}'",
								j, fileList.get(j), fileList.get(j));
					}
				}
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			LOGGER.error("Error in reading FlashLFQ quantification result file " + super.getFile(), e);
		}
	}

	protected void parseTitle(ArrayList<String> fileList, String[] reportTitle) {
		// TODO Auto-generated method stub
		intensityTitleIdMap = new HashMap<String, Integer>();
		idenTypeTitleIdMap = new HashMap<String, Integer>();
		this.intensityTitles = new String[fileList.size() * reportTitle.length];
		this.idenTypeTitles = new String[fileList.size()];

		String line = null;
		try {
			line = reader.readLine();
			this.title = line.split("\t");
			java.util.List<String> intensityCols = new java.util.ArrayList<>();
			for (int i = 0; i < title.length; i++) {
				if (title[i] != null && title[i].startsWith("Intensity_")) intensityCols.add(title[i]);
			}
			LOGGER.info("FlashLfqQuanPepReader: parseTitle(fileList, reportTitle): fileList size={}, reportTitle length={}; header Intensity_ columns ({}): {}",
					fileList.size(), reportTitle.length, intensityCols.size(), intensityCols);
			for (int i = 0; i < title.length; i++) {
				if (title[i].equals("Sequence")) {
					modSeqId = i;
				} else if (title[i].equals("Base Sequence")) {
					sequenceId = i;
				} else if (title[i].equals("Protein Groups")) {
					proteinId = i;
				} else if (title[i].startsWith("Intensity_")) {
					intensityTitleIdMap.put(title[i], i);
					for (int j = 0; j < fileList.size(); j++) {
						String fileJ = fileList.get(j);
						for (int k = 0; k < reportTitle.length; k++) {
							if (title[i].equals("Intensity_" + fileJ + "_" + reportTitle[k])) {
								this.intensityTitles[j * reportTitle.length + k] = title[i];
							}
						}
					}
				} else if (title[i].startsWith("Detection Type_")) {
					idenTypeTitleIdMap.put(title[i], i);
					for (int j = 0; j < fileList.size(); j++) {
						if (title[i].equals("Detection Type_" + fileList.get(j))) {
							this.idenTypeTitles[j] = title[i];
						}
					}
				}
			}

			if (sequenceId == -1) {
				sequenceId = modSeqId;
			}
			int matchedIntensity = 0;
			for (int i = 0; i < this.intensityTitles.length; i++) {
				if (this.intensityTitles[i] != null) matchedIntensity++;
			}
			int matchedIdenType = 0;
			for (int j = 0; j < this.idenTypeTitles.length; j++) {
				if (this.idenTypeTitles[j] != null) matchedIdenType++;
			}
			LOGGER.info("FlashLfqQuanPepReader: after parseTitle(fileList, reportTitle): matched {}/{} intensity columns, {}/{} Detection Type columns",
					matchedIntensity, this.intensityTitles.length, matchedIdenType, fileList.size());
		} catch (IOException e) {
			// TODO Auto-generated catch block
			LOGGER.error("Error in reading FlashLFQ quantification result file " + super.getFile(), e);
		}
	}
	
	protected FlashLfqQuanPeptide parse() {
		String[] cs = line.split("\t");
		String[] ps = proteinId == -1 ? new String[] {} : cs[proteinId].split(";");
		String[] proteins = new String[ps.length];
		boolean decoy = false;
		for (int i = 0; i < proteins.length; i++) {
			proteins[i] = ps[i];
			if (proteins[i].startsWith(PFindResultReader.REV)) {
				decoy = true;
			}
		}

		double[] intensity = new double[intensityTitles.length];
		for (int i = 0; i < intensity.length; i++) {
			if (intensityTitles[i] == null) {
				continue;
			}
			Integer id = this.intensityTitleIdMap.get(intensityTitles[i]);
			if (id != null && id > -1 && id < cs.length) {
				try {
					intensity[i] = Double.parseDouble(cs[id]);
				} catch (NumberFormatException e) {
					LOGGER.warn("FlashLfqQuanPepReader: could not parse intensity at column {} ('{}') for title '{}'",
							id, id < cs.length ? cs[id] : "?", intensityTitles[i]);
				}
			}
		}

		int[] idenType = new int[idenTypeTitles.length];
		for (int i = 0; i < idenType.length; i++) {
			if (idenTypeTitles[i] == null) {
				continue;
			}
			Integer id = this.idenTypeTitleIdMap.get(idenTypeTitles[i]);
			if (id != null && id > -1 && id < cs.length) {
				if (cs[id].equals("MSMS")) {
					idenType[i] = 1;
				} else if (cs[id].equals("MBR")) {
					idenType[i] = 2;
				} else if (cs[id].equals("NotDetected")) {
					idenType[i] = 0;
				}
			}
		}

		FlashLfqQuanPeptide peptide = new FlashLfqQuanPeptide(cs[sequenceId], cs[modSeqId], proteins, intensity,
				idenType, decoy);

		return peptide;
	}

	@Override
	public String getQuanMode() {
		// TODO Auto-generated method stub
		return MetaConstants.labelFree;
	}

	@Override
	public FlashLfqQuanPeptide[] getMetaPeptides() {
		// TODO Auto-generated method stub
		ArrayList<FlashLfqQuanPeptide> list = new ArrayList<FlashLfqQuanPeptide>();
		boolean loggedFirstRow = false;
		try {
			while ((line = reader.readLine()) != null) {
				FlashLfqQuanPeptide pep = parse();
				if (pep != null) {
					list.add(pep);
					if (!loggedFirstRow) {
						double[] firstIntensity = pep.getIntensity();
						int[] firstIdenType = pep.getIdenType();
						java.util.List<Double> intList = new java.util.ArrayList<>();
						for (double v : firstIntensity) intList.add(v);
						java.util.List<Integer> typeList = new java.util.ArrayList<>();
						for (int v : firstIdenType) typeList.add(v);
						LOGGER.info("FlashLfqQuanPepReader: first peptide intensity ({} values): {}; idenType ({} values): {}",
								firstIntensity.length, intList, firstIdenType.length, typeList);
						loggedFirstRow = true;
					}
				}
			}
			reader.close();

		} catch (IOException e) {
			// TODO Auto-generated catch block
			LOGGER.error("Error in reading peptide from " + this.getFile(), e);
		}

		FlashLfqQuanPeptide[] peps = list.toArray(new FlashLfqQuanPeptide[list.size()]);
		LOGGER.info("FlashLfqQuanPepReader: read {} peptides from {}", peps.length, this.getFile());
		return peps;
	}

	@Override
	public Object[] getTitleObjs() {
		// TODO Auto-generated method stub
		Object[] titleObjs;
		if (proteinId > -1) {
			titleObjs = new Object[3 + intensityTitles.length];
			titleObjs[0] = this.title[sequenceId];
			titleObjs[1] = this.title[proteinId];
			titleObjs[2] = "Intensity";
			for (int i = 0; i < intensityTitles.length; i++) {
				titleObjs[i + 3] = intensityTitles[i];
			}
		} else {
			titleObjs = new Object[2 + intensityTitles.length];
			titleObjs[0] = this.title[sequenceId];
			titleObjs[1] = "Intensity";
			for (int i = 0; i < intensityTitles.length; i++) {
				titleObjs[i + 2] = intensityTitles[i];
			}
		}

		return titleObjs;
	}

	@Override
	public String[] getIntensityTitle() {
		// TODO Auto-generated method stub
		return this.intensityTitles;
	}
	
	public static void main(String[] args) {
		
	}
}
