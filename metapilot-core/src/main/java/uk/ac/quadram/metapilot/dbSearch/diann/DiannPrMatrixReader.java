package uk.ac.quadram.metapilot.dbSearch.diann;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Locale;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import uk.ac.quadram.metapilot.quant.flashLFQ.FlashLfqQuanPepReader;
import uk.ac.quadram.metapilot.quant.flashLFQ.FlashLfqQuanPeptide;

public class DiannPrMatrixReader extends FlashLfqQuanPepReader {

	private static int chargeId = -1;
	private boolean isCsv;
	private static Logger LOGGER = LogManager.getLogger(DiannPrMatrixReader.class);
	
	private HashMap<String, Integer> runIdMap;

	private static final class PrecursorSummary {
		String[] proteins;
		String razorProtein;
		int miss;
		double mass;
		double pep = 1;
		double qValue = 1;
		int[] spCount;
		double[] individualPep;
	}

	public DiannPrMatrixReader(String quanPep) {
		this(new File(quanPep));
		// TODO Auto-generated constructor stub
	}

	public DiannPrMatrixReader(File quanPep) {
		super(quanPep);
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

			this.runIdMap = new HashMap<String, Integer>();

			for (int i = 0; i < title.length; i++) {
				if (title[i].equals("Stripped.Sequence") || title[i].equals("StrippedSequence")) {
					sequenceId = i;
				} else if (title[i].equals("Modified.Sequence") || title[i].equals("ModifiedSequence")) {
					modSeqId = i;
				} else if (title[i].equals("Precursor.Charge") || title[i].equals("PrecursorCharge")) {
					chargeId = i;
				} else if (title[i].equals("Protein.Group")) {
					proteinId = i;
				} else if (isIntensityColumn(title[i])) {
					String runName = normalizeRunName(title[i]);
					int runId = intensityTitleList.size();
					intensityTitleIdMap.put(runName, i);
					intensityTitleList.add(runName);
					putRunAlias(runIdMap, title[i], runId);
					putRunAlias(runIdMap, runName, runId);
				}
			}

			this.intensityTitles = intensityTitleList.toArray(new String[intensityTitleList.size()]);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			LOGGER.error("Error in reading FlashLFQ quantification result file " + super.getFile(), e);
		}
	}

	protected FlashLfqQuanPeptide parse() {
		if (line.trim().length() == 0) {
			return null;
		}

		String[] cs;
		if (this.isCsv) {
			cs = line.split(",");
		} else {
			cs = line.split("\t");
		}

		String[] ps = proteinId == -1 ? new String[] {} : cs[proteinId].split(";");
		String[] proteins = new String[ps.length];
		for (int i = 0; i < proteins.length; i++) {
			proteins[i] = ps[i];
		}

		int[] idenType = new int[intensityTitles.length];
		double[] intensity = new double[intensityTitles.length];

		for (int i = 0; i < intensity.length; i++) {
			int id = this.intensityTitleIdMap.get(intensityTitles[i]);
			if (id > -1 && id < cs.length) {
				if (cs[id].trim().length() > 0) {
					intensity[i] = Double.parseDouble(cs[id]);
				} else {
					intensity[i] = 0;
				}
				idenType[i] = (intensity[i] > 0) ? 1 : 0; // MS/MS when quantified, else unmatched
			} else {
				intensity[i] = 0;
				idenType[i] = 0;
			}
		}

		int charge = chargeId > -1 ? Integer.parseInt(cs[chargeId]) : 1;
		FlashLfqQuanPeptide peptide = new FlashLfqQuanPeptide(cs[sequenceId], cs[modSeqId], proteins, intensity,
				idenType, false);
		peptide.setCharges(new int[] { charge });

		return peptide;
	}

	@Override
	public FlashLfqQuanPeptide[] getMetaPeptides() {
		// TODO Auto-generated method stub
		HashMap<String, FlashLfqQuanPeptide> map = new HashMap<String, FlashLfqQuanPeptide>();
		ArrayList<FlashLfqQuanPeptide> list = new ArrayList<FlashLfqQuanPeptide>();
		try {
			while ((line = reader.readLine()) != null) {
				FlashLfqQuanPeptide pep = parse();
				if (pep != null) {

					String modSeq = pep.getModSeq();
					if (map.containsKey(modSeq)) {
						FlashLfqQuanPeptide pep0 = map.get(modSeq);
						int[] charge0 = pep0.getCharges();
						int[] charge = pep.getCharges();

						int[] charge1 = new int[charge0.length + 1];
						System.arraycopy(charge0, 0, charge1, 0, charge0.length);
						charge1[charge0.length] = charge[0];
						pep0.setCharges(charge1);

						double[] intensity0 = pep0.getIntensity();
						double[] intensity = pep.getIntensity();
						int[] idenType = new int[intensity0.length];
						for (int i = 0; i < intensity0.length; i++) {
							intensity0[i] += intensity[i];
							idenType[i] = (intensity0[i] > 0) ? 1 : 0; // MS/MS when quantified, else unmatched
						}
						pep0.setIdenType(idenType);

					} else {
						list.add(pep);
						map.put(modSeq, pep);
						double[] intensity = pep.getIntensity();
						int[] idenType = new int[intensity.length];
						for (int i = 0; i < intensity.length; i++) {
							idenType[i] = (intensity[i] > 0) ? 1 : 0; // MS/MS when quantified, else unmatched
						}
						pep.setIdenType(idenType);
					}
				}
			}
			reader.close();

		} catch (IOException e) {
			// TODO Auto-generated catch block
			LOGGER.error("Error in reading peptide from " + this.getFile(), e);
		}

		FlashLfqQuanPeptide[] peps = list.toArray(new FlashLfqQuanPeptide[list.size()]);
		return peps;
	}

	public FlashLfqQuanPeptide[] getMetaPeptides(HashSet<String> usedPepSet) {
		// TODO Auto-generated method stub
		HashMap<String, FlashLfqQuanPeptide> map = new HashMap<String, FlashLfqQuanPeptide>();
		ArrayList<FlashLfqQuanPeptide> list = new ArrayList<FlashLfqQuanPeptide>();
		try {
			while ((line = reader.readLine()) != null) {
				FlashLfqQuanPeptide pep = parse();
				if (pep != null) {
					String stripeSeq = pep.getSequence();
					String modSeq = pep.getModSeq();
					if (usedPepSet.contains(stripeSeq)) {
						if (map.containsKey(modSeq)) {
							FlashLfqQuanPeptide pep0 = map.get(modSeq);
							int[] charge0 = pep0.getCharges();
							int[] charge = pep.getCharges();

							int[] charge1 = new int[charge0.length + 1];
							System.arraycopy(charge0, 0, charge1, 0, charge0.length);
							charge1[charge0.length] = charge[0];
							pep0.setCharges(charge1);

							double[] intensity0 = pep0.getIntensity();
							double[] intensity = pep.getIntensity();
							int[] idenType = new int[intensity0.length];
							for (int i = 0; i < intensity0.length; i++) {
								intensity0[i] += intensity[i];
								idenType[i] = (intensity0[i] > 0) ? 1 : 0; // MS/MS when quantified, else unmatched
							}
							pep0.setIdenType(idenType);
						} else {
							list.add(pep);
							map.put(modSeq, pep);
							double[] intensity = pep.getIntensity();
							int[] idenType = new int[intensity.length];
							for (int i = 0; i < intensity.length; i++) {
								idenType[i] = (intensity[i] > 0) ? 1 : 0; // MS/MS when quantified, else unmatched
							}
							pep.setIdenType(idenType);
						}
					}
				}
			}
			reader.close();

		} catch (IOException e) {
			// TODO Auto-generated catch block
			LOGGER.error("Error in reading peptide from " + this.getFile(), e);
		}

		FlashLfqQuanPeptide[] peps = list.toArray(new FlashLfqQuanPeptide[list.size()]);		
		return peps;
	}
	
	public FlashLfqQuanPeptide[] getMetaPeptides(DiaNNPrecursor[] diaNNPrecursors) {
		// TODO Auto-generated method stub

		HashMap<String, PrecursorSummary> ppMap = buildPrecursorSummaryMap(diaNNPrecursors);

		HashMap<String, FlashLfqQuanPeptide> map = new HashMap<String, FlashLfqQuanPeptide>();
		ArrayList<FlashLfqQuanPeptide> list = new ArrayList<FlashLfqQuanPeptide>();
		try {
			while ((line = reader.readLine()) != null) {
				FlashLfqQuanPeptide peptide = parse();
				if (peptide != null) {
					String modSeq = peptide.getModSeq();
					if (ppMap.containsKey(modSeq)) {
						if (map.containsKey(modSeq)) {
							FlashLfqQuanPeptide pep0 = map.get(modSeq);
							int[] charge0 = pep0.getCharges();
							int[] charge = peptide.getCharges();

							int[] charge1 = new int[charge0.length + 1];
							System.arraycopy(charge0, 0, charge1, 0, charge0.length);
							charge1[charge0.length] = charge[0];
							pep0.setCharges(charge1);

							double[] intensity0 = pep0.getIntensity();
							double[] intensity = peptide.getIntensity();
							for (int i = 0; i < intensity0.length; i++) {
								intensity0[i] += intensity[i];
							}

							int[] spCount = pep0.getMs2Counts();
							int[] idenType = pep0.getIdenType();
							for (int i = 0; i < spCount.length; i++) {
								if (spCount[i] > 0) {
									idenType[i] = 1; // MS/MS
								} else if (intensity0[i] > 0) {
									idenType[i] = 2; // MBR
								} else {
									idenType[i] = 0; // unmatched
								}
							}
							peptide.setIdenType(idenType);

						} else {
							list.add(peptide);
							map.put(modSeq, peptide);

							PrecursorSummary summary = ppMap.get(modSeq);
							int[] spCount = summary.spCount.clone();
							int[] idenType = new int[intensityTitleIdMap.size()];
							double[] indiPEP = summary.individualPep.clone();
							double[] intensity = peptide.getIntensity();

							peptide.setProteins(summary.proteins);
							peptide.setRazorProtein(summary.razorProtein);
							peptide.setMissCleave(summary.miss);
							peptide.setMass(summary.mass);
							peptide.setPEP(summary.pep);
							peptide.setScore(summary.qValue);

							peptide.setMs2Counts(spCount);
							peptide.setIndividualPEP(indiPEP);
							
							for (int i = 0; i < spCount.length; i++) {
								if (spCount[i] > 0) {
									idenType[i] = 1; // MS/MS
								} else if (intensity[i] > 0) {
									idenType[i] = 2; // MBR
								} else {
									idenType[i] = 0; // unmatched
								}
							}
							peptide.setIdenType(idenType);
						}
					}
				}
			}
			reader.close();

		} catch (IOException e) {
			// TODO Auto-generated catch block
			LOGGER.error("Error in reading peptide from " + this.getFile(), e);
		}

		FlashLfqQuanPeptide[] peps = list.toArray(new FlashLfqQuanPeptide[list.size()]);
		return peps;
	}

	private HashMap<String, PrecursorSummary> buildPrecursorSummaryMap(DiaNNPrecursor[] diaNNPrecursors) {
		HashMap<String, PrecursorSummary> ppMap = new HashMap<String, PrecursorSummary>();
		if (diaNNPrecursors == null) {
			return ppMap;
		}
		int sampleCount = intensityTitleIdMap != null ? intensityTitleIdMap.size() : 0;
		for (DiaNNPrecursor pp : diaNNPrecursors) {
			if (pp == null) {
				continue;
			}
			String modseq = pp.getModSeqString();
			if (modseq == null) {
				continue;
			}
			PrecursorSummary summary = ppMap.get(modseq);
			if (summary == null) {
				summary = new PrecursorSummary();
				summary.spCount = new int[sampleCount];
				summary.individualPep = new double[sampleCount];
				Arrays.fill(summary.individualPep, 1.0);
				ppMap.put(modseq, summary);
			}
			Integer runIdObj = findRunId(pp.getRunString());
			if (runIdObj != null) {
				int runId = runIdObj;
				if (runId >= 0 && runId < summary.spCount.length) {
					summary.spCount[runId]++;
					if (pp.getPEP() < summary.individualPep[runId]) {
						summary.individualPep[runId] = pp.getPEP();
					}
				}
			}
			if (pp.getPEP() < summary.pep) {
				summary.pep = pp.getPEP();
			}
			if (pp.getGlobalQvalue() < summary.qValue) {
				summary.qValue = pp.getGlobalQvalue();
			}
			summary.miss = pp.getMissCleavage();
			summary.mass = pp.getMass();
			summary.proteins = pp.getProGroups();
			summary.razorProtein = pp.getRazorProtein();
		}
		return ppMap;
	}

	private Integer findRunId(String runName) {
		if (runIdMap == null || runName == null) {
			return null;
		}
		Integer runId = runIdMap.get(runName);
		if (runId != null) {
			return runId;
		}
		return runIdMap.get(normalizeRunName(runName));
	}

	private static void putRunAlias(HashMap<String, Integer> map, String runName, int runId) {
		if (map == null || runName == null) {
			return;
		}
		String trimmed = runName.trim();
		if (!trimmed.isEmpty()) {
			map.putIfAbsent(trimmed, runId);
		}
		String normalized = normalizeRunName(trimmed);
		if (!normalized.isEmpty()) {
			map.putIfAbsent(normalized, runId);
		}
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
