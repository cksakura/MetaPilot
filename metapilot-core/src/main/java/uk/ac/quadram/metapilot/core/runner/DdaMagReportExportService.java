package uk.ac.quadram.metapilot.core.runner;

import uk.ac.quadram.metapilot.core.tools.FormatTool;
import uk.ac.quadram.metapilot.task.io.pep.MetaPSM;
import uk.ac.quadram.metapilot.task.io.pep.MetaPeptide;
import uk.ac.quadram.metapilot.task.par.MetaConstants;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.math.BigInteger;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 * Shared DDA MAG TSV report writer.
 */
public final class DdaMagReportExportService {

	private DdaMagReportExportService() {
	}

	public static final class Result {
		public File idTsv;
		public File finalSummaryTsv;
		public File finalPeptidesTsv;
		public File razorProteinsTsv;
	}

	public static void writeIdentificationSummary(File combinedDir, File resultDir, File spectraDir,
			HashMap<String, Integer> fileIdMap, int[] filePsmCount, HashSet<String>[] fileSequenceSets,
			Map<String, String> expNameMap) throws IOException {
		if (combinedDir == null)
			throw new IOException("combinedDir is null");
		if (resultDir == null)
			throw new IOException("resultDir is null");
		combinedDir.mkdirs();
		resultDir.mkdirs();

		File idTsv = new File(combinedDir, "id.tsv");
		try (PrintWriter idWriter = new PrintWriter(idTsv)) {
			idWriter.println("Raw file\tUnique peptide count\tPSM count");
			String[] rawNames = fileIdMap.keySet().toArray(new String[0]);
			Arrays.sort(rawNames, Comparator.comparingInt(fileIdMap::get));
			for (String rawName : rawNames) {
				int idx = fileIdMap.get(rawName);
				int unique = idx >= 0 && idx < fileSequenceSets.length && fileSequenceSets[idx] != null
						? fileSequenceSets[idx].size()
						: 0;
				int psm = idx >= 0 && idx < filePsmCount.length ? filePsmCount[idx] : 0;
				idWriter.println(rawName + "\t" + unique + "\t" + psm);
			}
		}

		File finalSummary = new File(resultDir, "final_summary.tsv");
		try (PrintWriter writer = new PrintWriter(finalSummary)) {
			writer.println(
					"Raw file\tExperiment\tMS/MS\tMS/MS Identified\tMS/MS Identified [%]\tPeptide Sequences Identified");
			HashMap<String, int[]> idMap = readIdMap(idTsv);
			File ms2CountFile = spectraDir != null ? new File(spectraDir, "ms2Count.txt") : null;
			if (ms2CountFile != null && ms2CountFile.exists()) {
				writeSummaryFromMs2Count(writer, ms2CountFile, idMap, expNameMap);
			} else {
				writeSummaryFromIdentifications(writer, idMap, expNameMap);
			}
		}
	}

	public static void writePeptideReport(File out, ArrayList<MetaPeptide> peptideList, HashMap<String, MetaPSM> psmMap,
			HashMap<String, int[]> pepPsmCountMap, HashMap<String, int[]> pepChargeMap,
			HashMap<String, Integer> proteinIdMap, ArrayList<String> fileList, String quanMode, String[] intensityTitle)
			throws IOException {
		writePeptideReport(out, peptideList, psmMap, pepPsmCountMap, pepChargeMap, proteinIdMap, fileList, quanMode,
				intensityTitle, null);
	}

	public static void writePeptideReport(File out, ArrayList<MetaPeptide> peptideList, HashMap<String, MetaPSM> psmMap,
			HashMap<String, int[]> pepPsmCountMap, HashMap<String, int[]> pepChargeMap,
			HashMap<String, Integer> proteinIdMap, ArrayList<String> fileList, String quanMode, String[] intensityTitle,
			HashMap<String, String> razorProteinMap)
			throws IOException {
		if (peptideList == null)
			peptideList = new ArrayList<>();
		if (psmMap == null)
			psmMap = new HashMap<>();
		if (pepPsmCountMap == null)
			pepPsmCountMap = new HashMap<>();
		if (pepChargeMap == null)
			pepChargeMap = new HashMap<>();
		if (proteinIdMap == null)
			proteinIdMap = new HashMap<>();
		if (fileList == null)
			fileList = new ArrayList<>();
		if (intensityTitle == null)
			intensityTitle = new String[0];
		if (out.getParentFile() != null)
			out.getParentFile().mkdirs();
		try (PrintWriter writer = new PrintWriter(out)) {
			StringBuilder titlesb = new StringBuilder();
			titlesb.append("Sequence\tBase sequence\tLength\tMissed cleavages\tMass\tRazor protein\tMapped proteins\tCharges\tScore\tPEP\t");
			for (String file : fileList)
				titlesb.append("Identification type ").append(file).append("\t");
			for (String file : fileList)
				titlesb.append("MS/MS count ").append(file).append("\t");
			titlesb.append("Intensity\t");
			if (MetaConstants.labelFree.equals(quanMode)) {
				for (String title : intensityTitle)
					titlesb.append("Intensity ").append(title).append("\t");
			}
			titlesb.append("Reverse\tPotential contaminant\tid\tProtein group IDs\tMS/MS Count");
			writer.println(titlesb);

			int id = 0;
			for (MetaPeptide peptide : peptideList) {
				if (peptide == null)
					continue;
				String fullSeq = peptide.getModSeq();
				String sequence = peptide.getSequence();
				MetaPSM psm = psmMap.get(fullSeq);
				if (psm == null)
					continue;

				StringBuilder sb = new StringBuilder();
				sb.append(fullSeq).append("\t");
				sb.append(sequence).append("\t");
				sb.append(sequence != null ? sequence.length() : 0).append("\t");
				sb.append(psm.getMiss()).append("\t");
				sb.append(psm.getCalMass()).append("\t");
				sb.append(resolveRazorProtein(razorProteinMap, fullSeq, sequence)).append("\t");

				String[] proteins = peptide.getProteins();
				if (proteins != null) {
					for (String protein : proteins) {
						if (protein != null)
							sb.append(protein).append(";");
					}
					removeTrailingSemicolon(sb);
				}
				sb.append("\t");

				appendCharges(sb, pepChargeMap.get(fullSeq));
				sb.append(psm.getScore()).append("\t");
				sb.append(psm.getPEP()).append("\t");

				int totalMsCount = 0;
				int[] psmCount = pepPsmCountMap.get(fullSeq);
				double[] intensity = peptide.getIntensity();
				if (intensity == null)
					intensity = new double[0];
				if (psmCount == null)
					psmCount = new int[fileList.size()];
				for (int j = 0; j < fileList.size(); j++) {
					int count = j < psmCount.length ? psmCount[j] : 0;
					if (count > 0) {
						totalMsCount += count;
						sb.append("By MS/MS\t");
					} else if (j < intensity.length && intensity[j] > 0) {
						sb.append("MBR\t");
					} else {
						sb.append("NotDetected\t");
					}
				}
				for (int j = 0; j < fileList.size(); j++) {
					sb.append(j < psmCount.length ? psmCount[j] : 0).append("\t");
				}

				double totalIntensity = 0;
				for (double value : intensity)
					totalIntensity += value;
				sb.append(BigInteger.valueOf((long) Math.floor(totalIntensity))).append("\t");
				if (MetaConstants.labelFree.equals(quanMode)) {
					for (int j = 0; j < intensityTitle.length; j++) {
						double value = j < intensity.length ? intensity[j] : 0;
						sb.append(BigInteger.valueOf((long) Math.floor(value))).append("\t");
					}
				}

				sb.append("\t\t");
				sb.append(id).append("\t");
				if (proteins != null) {
					for (String protein : proteins) {
						if (proteinIdMap.containsKey(protein))
							sb.append(proteinIdMap.get(protein)).append(";");
					}
					removeTrailingSemicolon(sb);
				}
				sb.append("\t").append(totalMsCount);
				writer.println(sb);
				id++;
			}
		}
	}

	private static String resolveRazorProtein(HashMap<String, String> razorProteinMap, String modSeq, String sequence) {
		if (razorProteinMap == null) {
			return "";
		}
		String razor = modSeq != null ? razorProteinMap.get(modSeq) : null;
		if ((razor == null || razor.isEmpty()) && sequence != null) {
			razor = razorProteinMap.get(sequence);
		}
		return razor != null ? razor : "";
	}

	public static void writeProteinReport(File out, String[] proteins, HashMap<String, double[]> proteinIntensityMap,
			HashMap<String, Double> proteinTotalIntenMap, HashMap<String, HashSet<String>> totalPepCountMap,
			HashMap<String, HashSet<String>> razorPepCountMap, HashMap<String, Integer> psmCountMap, String quanMode,
			String[] intensityTitle, ArrayList<MetaPeptide> peptideList) throws IOException {
		writeProteinReport(out, proteins, proteinIntensityMap, proteinTotalIntenMap, totalPepCountMap,
				razorPepCountMap, psmCountMap, quanMode, intensityTitle, peptideList, null);
	}

	public static void writeProteinReport(File out, String[] proteins, HashMap<String, double[]> proteinIntensityMap,
			HashMap<String, Double> proteinTotalIntenMap, HashMap<String, HashSet<String>> totalPepCountMap,
			HashMap<String, HashSet<String>> razorPepCountMap, HashMap<String, Integer> psmCountMap, String quanMode,
			String[] intensityTitle, ArrayList<MetaPeptide> peptideList, HashMap<String, double[]> proteinScoreMap)
			throws IOException {
		if (proteins == null)
			proteins = new String[0];
		if (proteinIntensityMap == null)
			proteinIntensityMap = new HashMap<>();
		if (proteinTotalIntenMap == null)
			proteinTotalIntenMap = new HashMap<>();
		if (totalPepCountMap == null)
			totalPepCountMap = new HashMap<>();
		if (razorPepCountMap == null)
			razorPepCountMap = new HashMap<>();
		if (psmCountMap == null)
			psmCountMap = new HashMap<>();
		if (intensityTitle == null)
			intensityTitle = new String[0];
		if (proteinScoreMap == null)
			proteinScoreMap = new HashMap<>();
		if (out.getParentFile() != null)
			out.getParentFile().mkdirs();
		try (PrintWriter writer = new PrintWriter(out)) {
			StringBuilder titlesb = new StringBuilder();
			titlesb.append(
					"Protein IDs\tPeptide counts (all)\tPeptide counts (razor)\tNumber of PSMs\tprotein_score_S\tprotein_p_emp\tprotein_q_value\tIntensity\t");
			if (MetaConstants.labelFree.equals(quanMode) || MetaConstants.isobaricLabel.equals(quanMode)) {
				for (String title : intensityTitle)
					titlesb.append("Intensity ").append(title).append("\t");
			}
			titlesb.append("Reverse\tPotential contaminant\tid\tPeptide IDs\tPeptide is razor");
			writer.println(titlesb);

			DecimalFormat df = FormatTool.getDFE4();
			for (int i = 0; i < proteins.length; i++) {
				String protein = proteins[i];
				HashSet<String> totalPepSet = totalPepCountMap.getOrDefault(protein, new HashSet<>());
				HashSet<String> razorPepSet = razorPepCountMap.getOrDefault(protein, new HashSet<>());
				StringBuilder sb = new StringBuilder();
				sb.append(protein).append("\t");
				sb.append(totalPepSet.size()).append("\t");
				sb.append(razorPepSet.size()).append("\t");
				sb.append(psmCountMap.getOrDefault(protein, 0)).append("\t");
				double[] scores = proteinScoreMap.get(protein);
				if (scores != null && scores.length >= 3) {
					sb.append(df.format(scores[0])).append("\t");
					sb.append(df.format(scores[1])).append("\t");
					sb.append(df.format(scores[2])).append("\t");
				} else {
					sb.append("\t\t\t");
				}
				double totalIntensity = proteinTotalIntenMap.getOrDefault(protein, 0.0);
				sb.append(df.format(totalIntensity)).append("\t");
				double[] intensity = proteinIntensityMap.getOrDefault(protein, new double[0]);
				if (MetaConstants.labelFree.equals(quanMode) || MetaConstants.isobaricLabel.equals(quanMode)) {
					for (int j = 0; j < intensityTitle.length; j++) {
						double value = j < intensity.length ? intensity[j] : 0;
						sb.append(df.format(value)).append("\t");
					}
				}
				sb.append("\t\t").append(i).append("\t");
				appendPeptideIdsAndRazorFlags(sb, totalPepSet, razorPepSet, peptideList);
				writer.println(sb);
			}
		}
	}

	public static void writeAllProteinReport(File out, String[] proteins, HashMap<String, double[]> proteinIntensityMap,
			HashMap<String, Double> proteinTotalIntenMap, HashMap<String, HashSet<String>> totalPepCountMap,
			HashMap<String, HashSet<String>> razorPepCountMap, HashMap<String, Integer> psmCountMap, String quanMode,
			String[] intensityTitle, ArrayList<MetaPeptide> peptideList) throws IOException {
		writeAllProteinReport(out, proteins, proteinIntensityMap, proteinTotalIntenMap, totalPepCountMap,
				razorPepCountMap, psmCountMap, quanMode, intensityTitle, peptideList, null);
	}

	public static void writeAllProteinReport(File out, String[] proteins, HashMap<String, double[]> proteinIntensityMap,
			HashMap<String, Double> proteinTotalIntenMap, HashMap<String, HashSet<String>> totalPepCountMap,
			HashMap<String, HashSet<String>> razorPepCountMap, HashMap<String, Integer> psmCountMap, String quanMode,
			String[] intensityTitle, ArrayList<MetaPeptide> peptideList, HashMap<String, double[]> proteinScoreMap)
			throws IOException {
		if (proteins == null)
			proteins = new String[0];
		if (proteinIntensityMap == null)
			proteinIntensityMap = new HashMap<>();
		if (proteinTotalIntenMap == null)
			proteinTotalIntenMap = new HashMap<>();
		if (totalPepCountMap == null)
			totalPepCountMap = new HashMap<>();
		if (razorPepCountMap == null)
			razorPepCountMap = new HashMap<>();
		if (psmCountMap == null)
			psmCountMap = new HashMap<>();
		if (intensityTitle == null)
			intensityTitle = new String[0];
		if (proteinScoreMap == null)
			proteinScoreMap = new HashMap<>();
		if (out.getParentFile() != null)
			out.getParentFile().mkdirs();

		try (PrintWriter writer = new PrintWriter(out)) {
			StringBuilder header = new StringBuilder();
			header.append("Protein IDs\tPeptide counts (all)\tPeptide counts (razor)\tNumber of PSMs\t");
			header.append("protein_score_S\tprotein_p_emp\tprotein_q_value\tIntensity\t");
			if (MetaConstants.labelFree.equals(quanMode) || MetaConstants.isobaricLabel.equals(quanMode)) {
				for (String title : intensityTitle) {
					header.append("Intensity ").append(title).append("\t");
				}
			}
			header.append("Reverse\tPotential contaminant\tid\tPeptide IDs\tPeptide is razor");
			writer.println(header);

			DecimalFormat df = FormatTool.getDFE4();
			for (int i = 0; i < proteins.length; i++) {
				String protein = proteins[i];
				HashSet<String> totalPepSet = totalPepCountMap.getOrDefault(protein, new HashSet<>());
				HashSet<String> razorPepSet = razorPepCountMap.getOrDefault(protein, new HashSet<>());
				StringBuilder sb = new StringBuilder();
				sb.append(protein).append("\t");
				sb.append(totalPepSet.size()).append("\t");
				sb.append(razorPepSet.size()).append("\t");
				sb.append(psmCountMap.getOrDefault(protein, 0)).append("\t");
				double[] scores = proteinScoreMap.get(protein);
				if (scores != null && scores.length >= 3) {
					sb.append(df.format(scores[0])).append("\t");
					sb.append(df.format(scores[1])).append("\t");
					sb.append(df.format(scores[2])).append("\t");
				} else {
					sb.append("\t\t\t");
				}
				double totalIntensity = proteinTotalIntenMap.getOrDefault(protein, 0.0);
				sb.append(df.format(totalIntensity)).append("\t");
				double[] intensity = proteinIntensityMap.getOrDefault(protein, new double[0]);
				if (MetaConstants.labelFree.equals(quanMode) || MetaConstants.isobaricLabel.equals(quanMode)) {
					for (int j = 0; j < intensityTitle.length; j++) {
						double value = j < intensity.length ? intensity[j] : 0;
						sb.append(df.format(value)).append("\t");
					}
				}
				sb.append("\t\t").append(i).append("\t");
				appendPeptideIdsAndRazorFlags(sb, totalPepSet, razorPepSet, peptideList);
				writer.println(sb);
			}
		}
	}

	private static LinkedHashMap<String, int[]> readIdMap(File idTsv) throws IOException {
		LinkedHashMap<String, int[]> idMap = new LinkedHashMap<>();
		try (BufferedReader reader = new BufferedReader(new FileReader(idTsv))) {
			String line = reader.readLine();
			while ((line = reader.readLine()) != null) {
				String[] cs = line.split("\t");
				if (cs.length >= 3) {
					idMap.put(cs[0], new int[] { Integer.parseInt(cs[1]), Integer.parseInt(cs[2]) });
				}
			}
		}
		return idMap;
	}

	private static void writeSummaryFromMs2Count(PrintWriter writer, File ms2CountFile, Map<String, int[]> idMap,
			Map<String, String> expNameMap) throws IOException {
		HashMap<String, Integer> ms2CountMap = readMs2CountMap(ms2CountFile);
		int totalMs2Count = 0;
		int totalPsmCount = 0;
		int totalUniquePepCount = 0;
		for (Map.Entry<String, int[]> entry : idMap.entrySet()) {
			String rawName = entry.getKey();
			int[] idCounts = entry.getValue();
			Integer ms2Count = ms2CountMap.get(normalizeSampleName(rawName));
			if (ms2Count == null) {
				ms2Count = idCounts[1];
			}
			writer.println(summaryLine(rawName, expNameMap, ms2Count, idCounts[1], idCounts[0]));
			totalMs2Count += ms2Count;
			totalPsmCount += idCounts[1];
			totalUniquePepCount += idCounts[0];
		}
		writer.print(summaryLine("Total", null, totalMs2Count, totalPsmCount, totalUniquePepCount));
	}

	private static HashMap<String, Integer> readMs2CountMap(File ms2CountFile) throws IOException {
		HashMap<String, Integer> ms2CountMap = new HashMap<>();
		try (BufferedReader reader = new BufferedReader(new FileReader(ms2CountFile))) {
			String line = reader.readLine();
			while ((line = reader.readLine()) != null) {
				String[] cs = line.split("\t");
				if (cs.length < 2) {
					continue;
				}
				try {
					ms2CountMap.put(normalizeSampleName(cs[0]), Integer.parseInt(cs[1]));
				} catch (NumberFormatException e) {
					// Ignore malformed rows and let the caller fall back to identified counts.
				}
			}
		}
		return ms2CountMap;
	}

	private static void writeSummaryFromIdentifications(PrintWriter writer, Map<String, int[]> idMap,
			Map<String, String> expNameMap) {
		int totalPsmCount = 0;
		int totalUniquePepCount = 0;
		for (Map.Entry<String, int[]> entry : idMap.entrySet()) {
			String rawName = entry.getKey();
			int[] idCounts = entry.getValue();
			writer.println(summaryLine(rawName, expNameMap, idCounts[1], idCounts[1], idCounts[0]));
			totalPsmCount += idCounts[1];
			totalUniquePepCount += idCounts[0];
		}
		writer.print(summaryLine("Total", null, totalPsmCount, totalPsmCount, totalUniquePepCount));
	}

	private static String summaryLine(String rawName, Map<String, String> expNameMap, int ms2Count, int psmCount,
			int uniquePepCount) {
		String expName = resolveExperimentName(rawName, expNameMap);
		double ratio = ms2Count > 0 ? (double) psmCount / (double) ms2Count * 100.0 : 0;
		return rawName + "\t" + expName + "\t" + ms2Count + "\t" + psmCount + "\t" + FormatTool.getDF2().format(ratio)
				+ "\t" + uniquePepCount;
	}

	private static String resolveExperimentName(String rawName, Map<String, String> expNameMap) {
		if ("Total".equals(rawName)) {
			return "";
		}
		if (expNameMap == null) {
			return rawName;
		}
		String expName = expNameMap.get(rawName);
		if (expName != null) {
			return expName;
		}
		String normalized = normalizeSampleName(rawName);
		for (Map.Entry<String, String> entry : expNameMap.entrySet()) {
			if (normalizeSampleName(entry.getKey()).equals(normalized)) {
				return entry.getValue();
			}
		}
		return rawName;
	}

	private static String normalizeSampleName(String name) {
		if (name == null) {
			return "";
		}
		String fileName = new File(name.trim()).getName();
		int dot = fileName.lastIndexOf('.');
		return dot > 0 ? fileName.substring(0, dot) : fileName;
	}

	private static void appendCharges(StringBuilder sb, int[] charge) {
		if (charge == null) {
			sb.append("2\t");
			return;
		}
		StringBuilder chargesb = new StringBuilder();
		for (int j = 0; j < charge.length; j++) {
			if (charge[j] == 1)
				chargesb.append(j).append(";");
		}
		if (chargesb.length() > 1)
			sb.append(chargesb.substring(0, chargesb.length() - 1)).append("\t");
		else
			sb.append("2\t");
	}

	private static void appendPeptideIdsAndRazorFlags(StringBuilder sb, HashSet<String> totalPepSet,
			HashSet<String> razorPepSet, ArrayList<MetaPeptide> peptideList) {
		StringBuilder idSb = new StringBuilder();
		StringBuilder razorSb = new StringBuilder();
		if (peptideList != null) {
			for (int id = 0; id < peptideList.size(); id++) {
				String seq = peptideList.get(id).getSequence();
				if (totalPepSet.contains(seq)) {
					if (idSb.length() > 0) {
						idSb.append(";");
						razorSb.append(";");
					}
					idSb.append(id);
					razorSb.append(razorPepSet.contains(seq) ? "true" : "false");
				}
			}
		}
		sb.append(idSb).append("\t").append(razorSb);
	}

	private static void removeTrailingSemicolon(StringBuilder sb) {
		if (sb.length() > 0 && sb.charAt(sb.length() - 1) == ';') {
			sb.deleteCharAt(sb.length() - 1);
		}
	}
}
