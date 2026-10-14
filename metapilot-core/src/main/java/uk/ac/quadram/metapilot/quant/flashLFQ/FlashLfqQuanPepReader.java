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
	// Instance fields: these were static, so a reader inherited column positions
	// from whichever file was parsed before it (e.g. "Base Sequence" missing in the
	// new file kept the previous file's index).
	protected int sequenceId = -1;
	protected int modSeqId = -1;
	protected int proteinId = -1;

	protected String[] intensityTitles;
	protected HashMap<String, Integer> intensityTitleIdMap;
	protected String[] idenTypeTitles;
	protected HashMap<String, Integer> idenTypeTitleIdMap;
	protected String[] title;
	protected String line;

	private static Logger LOGGER = LogManager.getLogger(FlashLfqQuanPepReader.class);

	/*
	 * Column-match and content diagnostics. Callers (e.g. the DDA MAG task) report
	 * these through their own task logger: this class's LOGGER does not reach the
	 * per-task log, so a 0-column match used to be invisible there. A reader whose
	 * expected intensity columns are absent from the header still returns peptides,
	 * with all-zero intensity arrays, so callers must check getMatchedIntensityCount().
	 */
	/** Number of header columns starting with "Intensity_". */
	protected int headerIntensityColumnCount;
	/** First few header "Intensity_" column names, for diagnostics. */
	protected ArrayList<String> headerIntensityExamples = new ArrayList<String>();
	/** Intensity columns the caller expects (fileList.size(), or fileList.size() * reportTitle.length). */
	protected int expectedIntensityCount;
	/** Expected intensity columns that were found in the header. */
	protected int matchedIntensityCount;
	/** Expected intensity column names with no header match (all of them, in expected order). */
	protected ArrayList<String> unmatchedIntensityTitles = new ArrayList<String>();
	protected int expectedIdenTypeCount;
	protected int matchedIdenTypeCount;
	/** Set by getMetaPeptides(); -1 until then. */
	protected int peptideCount = -1;
	/** Peptides with at least one intensity value > 0; set by getMetaPeptides(), -1 until then. */
	protected int peptidesWithIntensity = -1;

	public FlashLfqQuanPepReader(String quanPep) {
		this(new File(quanPep));
	}

	public FlashLfqQuanPepReader(File quanPep) {
		super(quanPep);
		try {
			this.reader = new BufferedReader(new FileReader(quanPep));
		} catch (FileNotFoundException e) {
			LOGGER.error("Error in reading FlashLFQ search result file " + quanPep, e);
			throw new RuntimeException("Cannot open peptide matrix file: " + (quanPep != null ? quanPep.getAbsolutePath() : "null"), e);
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

	/**
	 * Reader for a peptide table whose intensity columns are named explicitly:
	 * column k is {@code "Intensity_" + intensityColumnSuffixes[k]} and the
	 * returned intensity arrays follow that order. Detection Type columns are
	 * still matched per raw file in {@code fileList}. Use it for tables whose
	 * columns are not a raw-file (x report title) grid, e.g. isobaric tables
	 * written with combined label columns; the writer should supply the suffixes
	 * (see PFindIsobaricTask.intensityColumnSuffixes).
	 */
	public static FlashLfqQuanPepReader withIntensityColumns(File quanPep, ArrayList<String> fileList,
			String[] intensityColumnSuffixes) {
		return new FlashLfqQuanPepReader(quanPep, fileList, intensityColumnSuffixes, true);
	}

	private FlashLfqQuanPepReader(File quanPep, ArrayList<String> fileList, String[] intensityColumnSuffixes,
			boolean explicitColumns) {
		super(quanPep);
		try {
			this.reader = new BufferedReader(new FileReader(quanPep));
		} catch (FileNotFoundException e) {
			LOGGER.error("Error in reading peptide quantification file " + quanPep, e);
			throw new RuntimeException("Cannot open peptide matrix file: "
					+ (quanPep != null ? quanPep.getAbsolutePath() : "null"), e);
		}
		this.parseTitleExplicit(fileList, intensityColumnSuffixes);
	}

	protected void parseTitleExplicit(ArrayList<String> fileList, String[] intensityColumnSuffixes) {
		intensityTitleIdMap = new HashMap<String, Integer>();
		idenTypeTitleIdMap = new HashMap<String, Integer>();
		this.intensityTitles = new String[intensityColumnSuffixes.length];
		this.idenTypeTitles = new String[fileList.size()];
		try {
			String line = reader.readLine();
			this.title = line != null ? line.split("\t") : new String[0];
			HashMap<String, Integer> expectedIndex = new HashMap<String, Integer>();
			for (int k = 0; k < intensityColumnSuffixes.length; k++) {
				expectedIndex.put("Intensity_" + intensityColumnSuffixes[k], k);
			}
			for (int i = 0; i < title.length; i++) {
				if (title[i].equals("Sequence")) {
					modSeqId = i;
				} else if (title[i].equals("Base Sequence")) {
					sequenceId = i;
				} else if (title[i].equals("Protein Groups")) {
					proteinId = i;
				} else if (title[i].startsWith("Intensity_")) {
					intensityTitleIdMap.put(title[i], i);
					Integer k = expectedIndex.get(title[i]);
					if (k != null) {
						this.intensityTitles[k] = title[i];
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
			for (int k = 0; k < intensityTitles.length; k++) {
				if (intensityTitles[k] != null) {
					matchedIntensity++;
				} else {
					unmatchedIntensityTitles.add("Intensity_" + intensityColumnSuffixes[k]);
				}
			}
			int matchedIdenType = 0;
			for (String t : idenTypeTitles) {
				if (t != null) {
					matchedIdenType++;
				}
			}
			recordHeaderIntensityColumns();
			this.expectedIntensityCount = intensityTitles.length;
			this.matchedIntensityCount = matchedIntensity;
			this.expectedIdenTypeCount = fileList.size();
			this.matchedIdenTypeCount = matchedIdenType;
			LOGGER.info("FlashLfqQuanPepReader: explicit intensity columns: matched {}/{}, Detection Type {}/{}",
					matchedIntensity, intensityTitles.length, matchedIdenType, fileList.size());
			if (matchedIntensity == 0) {
				LOGGER.warn("FlashLfqQuanPepReader: NO expected intensity column found - every peptide will have zero intensity. {}",
						getIntensityMatchSummary());
			}
		} catch (IOException e) {
			LOGGER.error("Error in reading peptide quantification file " + super.getFile(), e);
		}
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
			recordHeaderIntensityColumns();
			this.expectedIntensityCount = intensityTitles.length;
			this.matchedIntensityCount = intensityTitles.length;
			this.expectedIdenTypeCount = idenTypeTitles.length;
			this.matchedIdenTypeCount = idenTypeTitles.length;
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
			recordHeaderIntensityColumns();
			this.expectedIntensityCount = fileList.size();
			this.matchedIntensityCount = matchedIntensity;
			this.expectedIdenTypeCount = fileList.size();
			this.matchedIdenTypeCount = matchedIdenType;
			for (int j = 0; j < fileList.size(); j++) {
				if (this.intensityTitles[j] == null) {
					this.unmatchedIntensityTitles.add("Intensity_" + fileList.get(j));
				}
			}
			LOGGER.info("FlashLfqQuanPepReader: after parseTitle(fileList): matched {}/{} intensity columns, {}/{} Detection Type columns",
					matchedIntensity, fileList.size(), matchedIdenType, fileList.size());
			if (matchedIntensity == 0) {
				LOGGER.warn("FlashLfqQuanPepReader: NO expected intensity column found - every peptide will have zero intensity. {}",
						getIntensityMatchSummary());
			}
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
			recordHeaderIntensityColumns();
			this.expectedIntensityCount = this.intensityTitles.length;
			this.matchedIntensityCount = matchedIntensity;
			this.expectedIdenTypeCount = fileList.size();
			this.matchedIdenTypeCount = matchedIdenType;
			for (int j = 0; j < fileList.size(); j++) {
				for (int k = 0; k < reportTitle.length; k++) {
					if (this.intensityTitles[j * reportTitle.length + k] == null) {
						this.unmatchedIntensityTitles.add("Intensity_" + fileList.get(j) + "_" + reportTitle[k]);
					}
				}
			}
			LOGGER.info("FlashLfqQuanPepReader: after parseTitle(fileList, reportTitle): matched {}/{} intensity columns, {}/{} Detection Type columns",
					matchedIntensity, this.intensityTitles.length, matchedIdenType, fileList.size());
			if (matchedIntensity == 0) {
				LOGGER.warn("FlashLfqQuanPepReader: NO expected isobaric intensity column found - every peptide will have zero intensity. {}",
						getIntensityMatchSummary());
			}
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
		int withIntensity = 0;
		for (FlashLfqQuanPeptide pep : peps) {
			double[] values = pep.getIntensity();
			if (values == null) {
				continue;
			}
			for (double v : values) {
				if (v > 0) {
					withIntensity++;
					break;
				}
			}
		}
		this.peptideCount = peps.length;
		this.peptidesWithIntensity = withIntensity;
		LOGGER.info("FlashLfqQuanPepReader: read {} peptides from {} ({} with intensity > 0)", peps.length,
				this.getFile(), withIntensity);
		if (peps.length > 0 && withIntensity == 0) {
			LOGGER.warn("FlashLfqQuanPepReader: all {} peptides have zero intensity. {}", peps.length,
					getIntensityMatchSummary());
		}
		return peps;
	}

	private void recordHeaderIntensityColumns() {
		this.headerIntensityColumnCount = 0;
		this.headerIntensityExamples.clear();
		if (this.title == null) {
			return;
		}
		for (String t : this.title) {
			if (t != null && t.startsWith("Intensity_")) {
				this.headerIntensityColumnCount++;
				if (this.headerIntensityExamples.size() < 3) {
					this.headerIntensityExamples.add(t);
				}
			}
		}
	}

	public int getExpectedIntensityCount() {
		return expectedIntensityCount;
	}

	public int getMatchedIntensityCount() {
		return matchedIntensityCount;
	}

	public int getHeaderIntensityColumnCount() {
		return headerIntensityColumnCount;
	}

	public ArrayList<String> getUnmatchedIntensityTitles() {
		return unmatchedIntensityTitles;
	}

	public int getMatchedIdenTypeCount() {
		return matchedIdenTypeCount;
	}

	public int getExpectedIdenTypeCount() {
		return expectedIdenTypeCount;
	}

	/** Number of peptides read by getMetaPeptides(); -1 before it is called. */
	public int getPeptideCount() {
		return peptideCount;
	}

	/** Peptides with at least one intensity > 0, from getMetaPeptides(); -1 before it is called. */
	public int getPeptidesWithIntensity() {
		return peptidesWithIntensity;
	}

	/**
	 * One-line description of how the expected intensity columns matched this
	 * file's header, for task logs. Lists the first few unmatched expected names
	 * next to the first few header names, so a naming mismatch (or the wrong file)
	 * is visible at a glance.
	 */
	public String getIntensityMatchSummary() {
		StringBuilder sb = new StringBuilder();
		sb.append("file=").append(this.getFile() != null ? this.getFile().getAbsolutePath() : "null");
		sb.append("; header Intensity_ columns=").append(headerIntensityColumnCount);
		sb.append("; expected intensity columns matched=").append(matchedIntensityCount).append('/')
				.append(expectedIntensityCount);
		sb.append("; Detection Type columns matched=").append(matchedIdenTypeCount).append('/')
				.append(expectedIdenTypeCount);
		if (peptideCount >= 0) {
			sb.append("; peptides=").append(peptideCount).append(", with intensity > 0=")
					.append(peptidesWithIntensity);
		}
		if (!unmatchedIntensityTitles.isEmpty()) {
			sb.append("; first expected-but-missing: ")
					.append(unmatchedIntensityTitles.subList(0, Math.min(3, unmatchedIntensityTitles.size())));
			sb.append("; first header Intensity_ columns: ").append(headerIntensityExamples);
		}
		return sb.toString();
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
