package uk.ac.quadram.metapilot.task.igc;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.math.BigInteger;
import java.util.concurrent.ExecutionException;
import java.util.stream.Collectors;

import uk.ac.quadram.metapilot.core.math.MathTool;
import uk.ac.quadram.metapilot.core.taxonomy.Taxon;
import uk.ac.quadram.metapilot.core.taxonomy.TaxonomyDatabaseSQL;
import uk.ac.quadram.metapilot.core.taxonomy.TaxonomyRanks;
import uk.ac.quadram.metapilot.core.tools.FormatTool;
import uk.ac.quadram.metapilot.dbSearch.pfind.PFindResultReader;
import javafx.scene.control.ProgressBar;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters.ToolType;
import uk.ac.quadram.metapilot.task.MetaAbstractTask;
import uk.ac.quadram.metapilot.task.MetaPilotTask;
import uk.ac.quadram.metapilot.task.MetaPythonReportTask;
import uk.ac.quadram.metapilot.task.io.ReportManifestGenerator;
import uk.ac.quadram.metapilot.task.io.pep.MetaPeptide;
import uk.ac.quadram.metapilot.task.io.pep.MetaPSM;
import uk.ac.quadram.metapilot.task.par.MetaBaseParameter;
import uk.ac.quadram.metapilot.mdb.unipept.UnipeptRequester;
import uk.ac.quadram.metapilot.mdb.unipept.UnipProResult;
import uk.ac.quadram.metapilot.mdb.unipept.UnipProInfo;

import org.apache.commons.math3.stat.inference.TTest;

/**
 * Abstract base class for DDA IGC (Integrated Gene Catalog) workflow tasks.
 * 
 * This class provides common functionality and fields shared across all DDA IGC
 * workflows, including pFind, MaxQuant, Sage, and FragPipe implementations. It
 * handles common operations such as result directory creation, report
 * generation, and taxonomic analysis.
 */
public abstract class MetaPilotIgcTask extends MetaPilotTask {

	// Logger is inherited from MetaPilotTask as an instance field

	// Common fields for all IGC workflows
	protected String[] rawFiles;
	protected String[] rawFileNames;
	protected String[] quanExpNames;
	protected File[] separateResultFolders;
	protected double psmCoverThreshold = 0.95;
	protected File spectraDirFile;
	protected File rawDir;
	protected File quan_pep_file;
	protected File quan_pro_file;
	protected File final_pep_txt;
	protected File final_pro_txt;
	protected HashMap<String, String> expNameMap;
	protected File summaryMetaFile;
	protected File combinedFile;

	protected MetaPeptide[] metaPeptides;

	public MetaPilotIgcTask(MetaBaseParameter metaPar, MetaResourceParameters mrp, ProgressBar progressBar,
			ProgressBar secondBar) {
		super(metaPar, mrp, progressBar, secondBar);
		// Report templates will be copied in background during call() method
		// This prevents UI freeze during task construction
		logger.info(getTaskName() + ": Constructor completed, report templates will be copied during execution");
	}

	/**
	 * Helper method to refine proteins using Occam's razor principle (similar to
	 * DDA MAG workflow).
	 * 
	 * @param proteinToPeptides Map from protein names to their associated peptides
	 * @param totalPeptides     Total number of unique peptides
	 * @param threshold         Threshold for protein selection (e.g., 0.95 for 95%
	 *                          peptide coverage)
	 * @return Set of selected protein names
	 */
	protected boolean generateIdentificationSummary(HashMap<String, Integer> fileIdMap, int[] filePsmCount,
			HashSet<String>[] fileSequenceSets, String taskName) {

		try (PrintWriter idWriter = new PrintWriter(new File(combinedFile, "id.tsv"))) {
			idWriter.println("Raw file\tUnique peptide count\tPSM count");

			String[] rawNames = fileIdMap.keySet().toArray(new String[fileIdMap.size()]);
			Arrays.sort(rawNames, new Comparator<String>() {
				@Override
				public int compare(String o1, String o2) {
					return fileIdMap.get(o1) - fileIdMap.get(o2);
				}
			});

			for (int i = 0; i < rawNames.length; i++) {
				idWriter.println(rawNames[i] + "\t" + fileSequenceSets[i].size() + "\t" + filePsmCount[i]);
			}
		} catch (IOException e) {
			logger.error(taskName + ": error in writing identification information", e);
			System.err.println(
					format.format(new Date()) + "\t" + taskName + ": error in writing identification information");
			return false;
		}

		File final_summary_txt = new File(this.metaPar.getDbSearchResultFile(), "final_summary.tsv");
		try {

			PrintWriter countWriter = new PrintWriter(final_summary_txt);
			StringBuilder titlesb = new StringBuilder();
			titlesb.append("Raw file").append("\t");
			titlesb.append("Experiment").append("\t");
			titlesb.append("MS/MS").append("\t");
			titlesb.append("MS/MS Identified").append("\t");
			titlesb.append("MS/MS Identified [%]").append("\t");
			titlesb.append("Peptide Sequences Identified");

			countWriter.println(titlesb);

			int totalms2count = 0;
			int totalpsmCount = 0;
			int totalUnipepCount = 0;

			BufferedReader idReader = new BufferedReader(new FileReader(new File(combinedFile, "id.tsv")));
			HashMap<String, int[]> idMap = new HashMap<String, int[]>();
			String line = idReader.readLine();
			while ((line = idReader.readLine()) != null) {
				String[] cs = line.split("\t");
				idMap.put(cs[0], new int[] { Integer.parseInt(cs[1]), Integer.parseInt(cs[2]) });
			}
			idReader.close();

			File ms2CountFile = new File(metaPar.getSpectraFile(), "ms2Count.txt");
			BufferedReader ms2CountReader = new BufferedReader(new FileReader(ms2CountFile));
			line = ms2CountReader.readLine();

			while ((line = ms2CountReader.readLine()) != null) {
				String[] cs = line.split("\t");
				String rawName = cs[0];

				if (idMap.containsKey(rawName)) {

					StringBuilder sb = new StringBuilder();
					sb.append(cs[0]).append("\t");

					if (this.expNameMap.containsKey(cs[0])) {
						sb.append(expNameMap.get(cs[0])).append("\t");
					} else {
						sb.append(cs[0]).append("\t");
					}

					int ms2count = Integer.parseInt(cs[1]);

					int[] idCounts = idMap.get(rawName);

					sb.append(ms2count).append("\t");
					sb.append(idCounts[1]).append("\t");

					double ratio = 0;
					if (ms2count > 0) {
						ratio = (double) idCounts[1] / (double) ms2count * 100.0;
					}

					sb.append(FormatTool.getDF2().format(ratio)).append("\t");

					sb.append(idCounts[0]);

					countWriter.println(sb);

					totalms2count += ms2count;
					totalpsmCount += idCounts[1];
					totalUnipepCount += idCounts[0];
				}
			}

			ms2CountReader.close();

			StringBuilder sb = new StringBuilder();
			sb.append("Total").append("\t");
			sb.append("\t");
			sb.append(totalms2count).append("\t");
			sb.append(totalpsmCount).append("\t");

			double ratio = 0;
			if (totalms2count > 0) {
				ratio = (double) totalpsmCount / (double) totalms2count * 100.0;
			}

			sb.append(FormatTool.getDF2().format(ratio)).append("\t");

			sb.append(totalUnipepCount);

			countWriter.print(sb);

			countWriter.close();

		} catch (IOException e) {
			// TODO Auto-generated catch block
			logger.error(getTaskName() + ": error in writing PSM information to " + final_summary_txt.getName(), e);
			System.err.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": error in writing PSM information to " + final_summary_txt.getName());
		}

		if (!final_summary_txt.exists() || final_summary_txt.length() == 0) {
			return false;

		} else {
			return true;
		}
	}

	protected boolean generatePeptideReport(ArrayList<MetaPeptide> peptideList, HashMap<String, MetaPSM> psmMap,
			HashMap<String, int[]> pepPsmCountMap, HashMap<String, int[]> pepChargeMap,
			HashMap<String, Integer> proteinIdMap, ArrayList<String> fileList, String quanMode, String[] intensityTitle,
			String taskName) {
		return generatePeptideReport(peptideList, psmMap, pepPsmCountMap, pepChargeMap, proteinIdMap, fileList,
				quanMode, intensityTitle, taskName, null);
	}

	protected boolean generatePeptideReport(ArrayList<MetaPeptide> peptideList, HashMap<String, MetaPSM> psmMap,
			HashMap<String, int[]> pepPsmCountMap, HashMap<String, int[]> pepChargeMap,
			HashMap<String, Integer> proteinIdMap, ArrayList<String> fileList, String quanMode, String[] intensityTitle,
			String taskName, HashMap<String, String> razorProteinMap) {

		logger.info(taskName + ": exporting peptide report started");
		System.out.println(format.format(new Date()) + "\t" + taskName + ": exporting peptide report started");

		try {
			final_pep_txt = new File(this.metaPar.getDbSearchResultFile(), "final_peptides.tsv");
			PrintWriter writer = new PrintWriter(final_pep_txt);

			StringBuilder titlesb = new StringBuilder();
			titlesb.append("Sequence").append("\t");
			titlesb.append("Base sequence").append("\t");
			titlesb.append("Length").append("\t");
			titlesb.append("Missed cleavages").append("\t");
			titlesb.append("Mass").append("\t");
			titlesb.append("Razor protein").append("\t");
			titlesb.append("Mapped proteins").append("\t");
			titlesb.append("Charges").append("\t");
			titlesb.append("Score").append("\t");
			titlesb.append("PEP").append("\t");

			for (int i = 0; i < fileList.size(); i++) {
				titlesb.append("Identification type " + fileList.get(i)).append("\t");
			}

			for (int i = 0; i < fileList.size(); i++) {
				titlesb.append("MS/MS count " + fileList.get(i)).append("\t");
			}

			titlesb.append("Intensity").append("\t");

			for (int i = 0; i < intensityTitle.length; i++) {
				titlesb.append("Intensity " + intensityTitle[i]).append("\t");
			}

			titlesb.append("Reverse").append("\t");
			titlesb.append("Potential contaminant").append("\t");
			titlesb.append("id").append("\t");
			titlesb.append("Protein group IDs").append("\t");
			titlesb.append("Total MS/MS count");

			writer.println(titlesb);

			int id = 0;
			for (MetaPeptide peptide : peptideList) {
				String fullSeq = peptide.getModSeq();
				MetaPSM psm = psmMap.get(fullSeq);

				StringBuilder sb = new StringBuilder();
				sb.append(psm.getFullSequence()).append("\t");
				sb.append(psm.getBaseSequence()).append("\t");
				sb.append(psm.getBaseSequence().length()).append("\t");
				sb.append(psm.getMiss()).append("\t");

				double mass = peptide.getMass();
				if (mass == 0 && psm != null) {
					mass = psm.getCalMass();
				}
				sb.append(mass).append("\t");
				sb.append(resolveRazorProtein(razorProteinMap, fullSeq, psm.getBaseSequence())).append("\t");

				String[] pros = peptide.getProteins();
				if (pros != null) {
					StringBuilder proSb = new StringBuilder();
					for (String pro : pros) {
						proSb.append(pro).append(";");
					}
					if (proSb.length() > 0) {
						proSb.deleteCharAt(proSb.length() - 1);
					}
					sb.append(proSb);
				}
				sb.append("\t");

				int[] charge = pepChargeMap.get(fullSeq);
				if (charge != null) {
					StringBuilder chargesb = new StringBuilder();
					for (int j = 0; j < charge.length; j++) {
						if (charge[j] == 1) {
							chargesb.append(j).append(";");
						}
					}
					if (chargesb.length() > 0) {
						if (chargesb.charAt(chargesb.length() - 1) == ';') {
							chargesb.deleteCharAt(chargesb.length() - 1);
						}
						sb.append(chargesb).append("\t");
					} else {
						sb.append("2\t");
					}
				} else {
					sb.append("2\t");
				}

				sb.append(psm.getScore()).append("\t");
				sb.append(psm.getPEP()).append("\t");

				int totalMsCount = 0;
				int[] psmCount = pepPsmCountMap.get(fullSeq);
				if (psmCount != null) {
					double[] intensity = peptide.getIntensity();
					for (int j = 0; j < psmCount.length; j++) {
						if (psmCount[j] > 0) {
							totalMsCount += psmCount[j];
							sb.append("By MS/MS").append("\t");
						} else {
							if (intensity[j] > 0) {
								sb.append("MBR\t");
							} else {
								sb.append("NotDetected\t");
							}
						}
					}

					for (int j = 0; j < psmCount.length; j++) {
						sb.append(psmCount[j]).append("\t");
					}
				} else {
					for (int j = 0; j < fileList.size(); j++) {
						sb.append("NotDetected\t");
					}
					for (int j = 0; j < fileList.size(); j++) {
						sb.append(0).append("\t");
					}
				}

				double totalIntensity = 0;
				double[] intensity = peptide.getIntensity();
				for (int j = 0; j < intensity.length; j++) {
					totalIntensity += intensity[j];
				}

				sb.append(BigInteger.valueOf((long) Math.floor(totalIntensity))).append("\t");

				for (int j = 0; j < intensity.length; j++) {
					sb.append(BigInteger.valueOf((long) Math.floor(intensity[j]))).append("\t");
				}

				sb.append("").append("\t");
				sb.append("").append("\t");
				sb.append(id).append("\t");

				StringBuilder groupSb = new StringBuilder();
				for (String pro : pros) {
					if (proteinIdMap.containsKey(pro)) {
						groupSb.append(proteinIdMap.get(pro)).append(";");
					}
				}
				if (groupSb.length() > 0) {
					groupSb.deleteCharAt(groupSb.length() - 1);
				}
				sb.append(groupSb).append("\t");

				sb.append(totalMsCount);

				writer.println(sb);
				id++;
			}

			writer.close();

		} catch (Exception e) {
			logger.error(taskName + ": error in writing the peptides to " + final_pep_txt, e);
			System.err.println(format.format(new Date()) + "\t" + taskName + ": error in writing the peptides to "
					+ final_pep_txt);
			return false;
		}

		logger.info(taskName + ": peptide report has been exported to " + final_pep_txt);
		System.out.println(
				format.format(new Date()) + "\t" + taskName + ": peptide report has been exported to " + final_pep_txt);

		if (!final_pep_txt.exists() || final_pep_txt.length() == 0) {
			return false;
		}

		return true;
	}

	private String resolveRazorProtein(HashMap<String, String> razorProteinMap, String modSeq, String sequence) {
		if (razorProteinMap == null) {
			return "";
		}
		String razor = modSeq != null ? razorProteinMap.get(modSeq) : null;
		if ((razor == null || razor.isEmpty()) && sequence != null) {
			razor = razorProteinMap.get(sequence);
		}
		return razor != null ? razor : "";
	}

	protected boolean generateProteinReport(String[] proteins, HashMap<String, double[]> proteinIntensityMap,
			HashMap<String, Double> proteinTotalIntensityMap, HashMap<String, HashSet<String>> totalPepCountMap,
			HashMap<String, HashSet<String>> razorPepCountMap, HashMap<String, Integer> proteinPsmCountMap,
			HashMap<String, Integer> proteinIdMap, String quanMode, String[] intensityTitle, String taskName) {

		logger.info(taskName + ": exporting protein report started");
		System.out.println(format.format(new Date()) + "\t" + taskName + ": exporting protein report started");

		try {
			final_pro_txt = new File(this.metaPar.getDbSearchResultFile(), "razor_proteins.tsv");
			PrintWriter writer = new PrintWriter(final_pro_txt);

			StringBuilder titlesb = new StringBuilder();
			titlesb.append("Protein IDs").append("\t");
			titlesb.append("Peptide counts (all)").append("\t");
			titlesb.append("Peptide counts (razor)").append("\t");
			titlesb.append("Number of PSMs").append("\t");
			titlesb.append("PEP").append("\t");
			titlesb.append("Intensity").append("\t");

			for (int i = 0; i < intensityTitle.length; i++) {
				titlesb.append(intensityTitle[i]).append("\t");
			}

			titlesb.append("Reverse").append("\t");
			titlesb.append("Potential contaminant").append("\t");
			titlesb.append("id").append("\t");

			writer.println(titlesb);

			for (String protein : proteins) {
				StringBuilder sb = new StringBuilder();
				sb.append(protein).append("\t");

				int totalPepCount = 0;
				if (totalPepCountMap.containsKey(protein)) {
					totalPepCount = totalPepCountMap.get(protein).size();
				}
				sb.append(totalPepCount).append("\t");

				int razorPepCount = 0;
				if (razorPepCountMap.containsKey(protein)) {
					razorPepCount = razorPepCountMap.get(protein).size();
				}
				sb.append(razorPepCount).append("\t");

				int psmCount = 0;
				if (proteinPsmCountMap.containsKey(protein)) {
					psmCount = proteinPsmCountMap.get(protein);
				}
				sb.append(psmCount).append("\t");

				double totalIntensity = proteinTotalIntensityMap.containsKey(protein)
						? proteinTotalIntensityMap.get(protein)
						: 0.0;
				sb.append(0.0).append("\t");
				sb.append(BigInteger.valueOf((long) Math.floor(totalIntensity))).append("\t");

				double[] intensity = proteinIntensityMap.get(protein);
				if (intensity != null) {
					for (int i = 0; i < intensity.length; i++) {
						sb.append(BigInteger.valueOf((long) Math.floor(intensity[i]))).append("\t");
					}
				} else {
					for (int i = 0; i < intensityTitle.length; i++) {
						sb.append(0).append("\t");
					}
				}

				sb.append("").append("\t");
				sb.append("").append("\t");
				sb.append(proteinIdMap.get(protein));

				writer.println(sb);
			}

			writer.close();

		} catch (Exception e) {
			logger.error(taskName + ": error in writing the proteins to " + final_pro_txt, e);
			System.err.println(format.format(new Date()) + "\t" + taskName + ": error in writing the proteins to "
					+ final_pro_txt);
			return false;
		}

		logger.info(taskName + ": protein report has been exported to " + final_pro_txt);
		System.out.println(
				format.format(new Date()) + "\t" + taskName + ": protein report has been exported to " + final_pro_txt);

		if (!final_pro_txt.exists() || final_pro_txt.length() == 0) {
			return false;
		}

		return true;
	}

	/**
	 * Export an additional all-protein report (no razor-only restriction) using the
	 * same table schema as razor_proteins.tsv.
	 */
	protected boolean generateAllProteinReport(String[] proteins, HashMap<String, double[]> proteinIntensityMap,
			HashMap<String, Double> proteinTotalIntensityMap, HashMap<String, HashSet<String>> totalPepCountMap,
			HashMap<String, HashSet<String>> razorPepCountMap, HashMap<String, Integer> proteinPsmCountMap,
			String quanMode, String[] intensityTitle, String taskName) {

		File outFile = new File(this.metaPar.getDbSearchResultFile(), "all_proteins.tsv");
		logger.info(taskName + ": exporting all-protein report started");
		System.out.println(format.format(new Date()) + "\t" + taskName + ": exporting all-protein report started");

		try (PrintWriter writer = new PrintWriter(outFile)) {
			StringBuilder titlesb = new StringBuilder();
			titlesb.append("Protein IDs").append("\t");
			titlesb.append("Peptide counts (all)").append("\t");
			titlesb.append("Peptide counts (razor)").append("\t");
			titlesb.append("Number of PSMs").append("\t");
			titlesb.append("PEP").append("\t");
			titlesb.append("Intensity").append("\t");

			for (int i = 0; i < intensityTitle.length; i++) {
				titlesb.append(intensityTitle[i]).append("\t");
			}

			titlesb.append("Reverse").append("\t");
			titlesb.append("Potential contaminant").append("\t");
			titlesb.append("id");

			writer.println(titlesb);

			for (int idx = 0; idx < proteins.length; idx++) {
				String protein = proteins[idx];
				StringBuilder sb = new StringBuilder();
				sb.append(protein).append("\t");

				int totalPepCount = totalPepCountMap.containsKey(protein) ? totalPepCountMap.get(protein).size() : 0;
				sb.append(totalPepCount).append("\t");

				int razorPepCount = razorPepCountMap.containsKey(protein) ? razorPepCountMap.get(protein).size() : 0;
				sb.append(razorPepCount).append("\t");

				int psmCount = proteinPsmCountMap.containsKey(protein) ? proteinPsmCountMap.get(protein) : 0;
				sb.append(psmCount).append("\t");

				double totalIntensity = proteinTotalIntensityMap.containsKey(protein) ? proteinTotalIntensityMap.get(protein)
						: 0.0;
				sb.append(0.0).append("\t");
				sb.append(BigInteger.valueOf((long) Math.floor(totalIntensity))).append("\t");

				double[] intensity = proteinIntensityMap.get(protein);
				if (intensity != null) {
					for (int i = 0; i < intensity.length; i++) {
						sb.append(BigInteger.valueOf((long) Math.floor(intensity[i]))).append("\t");
					}
				} else {
					for (int i = 0; i < intensityTitle.length; i++) {
						sb.append(0).append("\t");
					}
				}

				sb.append("").append("\t");
				sb.append("").append("\t");
				sb.append(idx);

				writer.println(sb);
			}
		} catch (Exception e) {
			logger.error(taskName + ": error in writing all proteins to " + outFile, e);
			System.err.println(
					format.format(new Date()) + "\t" + taskName + ": error in writing all proteins to " + outFile);
			return false;
		}

		logger.info(taskName + ": all-protein report has been exported to " + outFile);
		System.out.println(format.format(new Date()) + "\t" + taskName + ": all-protein report has been exported to "
				+ outFile);

		return outFile.exists() && outFile.length() > 0;
	}

	protected boolean performFunctionalAnnotation(String[] proteins, HashMap<String, double[]> proteinIntensityMap,
			HashMap<String, Double> proteinTotalIntensityMap, HashMap<String, HashSet<String>> totalPepCountMap,
			HashMap<String, UnipProInfo> proteinInfoMap, String[] intensityTitle, String taskName) {

		logger.info(taskName + ": functional annotation started");
		System.out.println(format.format(new Date()) + "\t" + taskName + ": functional annotation started");

		File funcFile = new File(metaPar.getDbSearchResultFile(), "functional_annotation");
		if (!funcFile.exists()) {
			funcFile.mkdir();
		}

		// Write functions.tsv directly from UnipProInfo data (no XML generation needed)
		File funTsv = new File(funcFile, "functions.tsv");
		try (PrintWriter writer = new PrintWriter(funTsv)) {
			// Build header
			StringBuilder header = new StringBuilder();
			header.append("Protein ID\t");
			header.append("Protein Name\t");
			header.append("EC Number\t");
			header.append("GO Terms\t");
			header.append("Taxonomic ID\t");
			header.append("Peptide Count\t");
			header.append("Total Intensity");

			// Add intensity columns
			for (String intensity : intensityTitle) {
				if (intensity.startsWith("Intensity ")) {
					header.append("\t").append(intensity);
				} else {
					header.append("\tIntensity ").append(intensity);
				}
			}

			writer.println(header);

			// Write protein data
			for (String protein : proteins) {
				UnipProInfo info = proteinInfoMap.get(protein);
				if (info == null) {
					continue; // Skip proteins without UnipProInfo
				}

				StringBuilder sb = new StringBuilder();
				sb.append(protein).append("\t");
				sb.append(info.getName() != null ? info.getName() : "").append("\t");

				// Extract EC numbers (if available)
				String ecNumbers = "";
				if (info.getEc() != null && info.getEc().length > 0) {
					ecNumbers = String.join(";", info.getEc());
				}
				sb.append(ecNumbers).append("\t");

				// Extract GO terms (if available)
				String goTerms = "";
				if (info.getGo() != null && info.getGo().length > 0) {
					goTerms = String.join(";", info.getGo());
				}
				sb.append(goTerms).append("\t");

				sb.append(info.getTaxId()).append("\t");

				// Add peptide count and intensity data
				if (totalPepCountMap.containsKey(protein)) {
					sb.append(totalPepCountMap.get(protein).size()).append("\t");
				} else {
					sb.append("0").append("\t");
				}

				if (proteinTotalIntensityMap.containsKey(protein)) {
					sb.append(proteinTotalIntensityMap.get(protein));
				} else {
					sb.append("0");
				}

				// Add individual intensity values
				if (proteinIntensityMap.containsKey(protein)) {
					double[] intensities = proteinIntensityMap.get(protein);
					for (double intensity : intensities) {
						sb.append("\t").append(intensity);
					}
				}

				writer.println(sb);
			}

		} catch (IOException e) {
			logger.error(taskName + ": error in writing functional annotation to " + funTsv.getName(), e);
			System.err.println(format.format(new Date()) + "\t" + taskName
					+ ": error in writing functional annotation to " + funTsv.getName());
			return false;
		}

		logger.info(taskName + ": functional annotation completed - wrote " + proteins.length + " proteins to "
				+ funTsv.getName());
		System.out.println(format.format(new Date()) + "\t" + taskName + ": functional annotation completed");

		return true;
	}

	protected boolean performTaxonomyAnalysis(String[] proteins, HashMap<String, double[]> proteinIntensityMap,
			HashMap<String, Integer> proteinPsmCountMap, HashMap<String, HashSet<String>> totalPepCountMap,
			HashMap<String, HashSet<String>> razorPepCountMap, HashMap<String, UnipProInfo> proteinInfoMap,
			String[] intensityTitle, String taskName) {

		logger.info(taskName + ": taxonomy analysis started");
		System.out.println(format.format(new Date()) + "\t" + taskName + ": taxonomy analysis started");

		File taxFile = new File(metaPar.getDbSearchResultFile(), "taxonomy_analysis");
		if (!taxFile.exists()) {
			taxFile.mkdir();
		}

		// Build protein-to-taxon mapping from UnipProInfo
		HashMap<String, Integer> proteinToTaxonMap = new HashMap<String, Integer>();
		HashMap<Integer, HashSet<String>> taxonToProteins = new HashMap<Integer, HashSet<String>>();

		for (String protein : proteins) {
			UnipProInfo info = proteinInfoMap.get(protein);
			if (info != null && info.getTaxId() > 0) {
				int taxonId = info.getTaxId();
				proteinToTaxonMap.put(protein, taxonId);

				if (!taxonToProteins.containsKey(taxonId)) {
					taxonToProteins.put(taxonId, new HashSet<String>());
				}
				taxonToProteins.get(taxonId).add(protein);
			}
		}

		logger.info(taskName + ": Mapped " + proteinToTaxonMap.size() + " proteins to " + taxonToProteins.size()
				+ " taxonomic IDs");

		// Get taxonomic lineage information for each taxon ID (similar to MAG workflow)
		HashMap<Integer, Taxon> taxonInfoMap = new HashMap<Integer, Taxon>();
		TaxonomyDatabaseSQL taxDb = null;
		try {
			taxDb = new TaxonomyDatabaseSQL(mrp.getToolPath(ToolType.FUNCDB));
			for (Integer taxonId : taxonToProteins.keySet()) {
				Taxon taxon = taxDb.getTaxonFromId(taxonId);
				if (taxon != null) {
					taxonInfoMap.put(taxonId, taxon);
				}
			}
		} catch (Exception e) {
			logger.error(taskName + ": Error accessing taxonomy database", e);
			System.err.println(format.format(new Date()) + "\t" + taskName + ": Error accessing taxonomy database");
			return false;
		}

		// Build taxonomic rank maps (similar to MAG workflow but protein-based)
		TaxonomyRanks[] mainRanks = TaxonomyRanks.getMainRanks();
		HashMap<String, String[]>[] taxRankMaps = new HashMap[mainRanks.length];
		HashMap<String, double[]>[] taxIntensityMaps = new HashMap[mainRanks.length];
		HashMap<String, Integer>[] taxPsmCountMaps = new HashMap[mainRanks.length];
		HashMap<String, HashSet<String>>[] taxPepCountMaps = new HashMap[mainRanks.length];
		HashMap<String, Integer>[] taxProteinCountMaps = new HashMap[mainRanks.length];

		for (int i = 0; i < taxRankMaps.length; i++) {
			taxRankMaps[i] = new HashMap<String, String[]>();
			taxIntensityMaps[i] = new HashMap<String, double[]>();
			taxPsmCountMaps[i] = new HashMap<String, Integer>();
			taxPepCountMaps[i] = new HashMap<String, HashSet<String>>();
			taxProteinCountMaps[i] = new HashMap<String, Integer>();
		}

		// Process each protein and aggregate to taxonomic ranks (like MAG workflow)
		for (String protein : proteins) {
			if (!proteinToTaxonMap.containsKey(protein)) {
				continue; // Skip proteins without taxonomic information
			}

			Integer taxonId = proteinToTaxonMap.get(protein);
			Taxon taxon = taxonInfoMap.get(taxonId);
			if (taxon == null) {
				continue;
			}

			// Get protein data
			double[] proteinIntensity = proteinIntensityMap.get(protein);
			Integer psmCount = proteinPsmCountMap.get(protein);
			HashSet<String> pepCount = totalPepCountMap.get(protein);

			if (proteinIntensity == null || psmCount == null || pepCount == null) {
				continue;
			}

			// Get taxonomic lineage using the same approach as existing IGC workflow
			if (taxon.getMainParentIds() == null) {
				try {
					taxon.setMainParentIds(taxDb.getMainParentTaxonIds(taxonId));
				} catch (Exception e) {
					logger.warn(taskName + ": Could not get lineage for taxon " + taxonId);
					continue;
				}
			}

			int[] mainTaxonIds = taxon.getMainParentIds();
			if (mainTaxonIds == null || mainTaxonIds.length == 0) {
				continue;
			}

			// Process lineage - handle missing ranks gracefully
			// Note: -1 values indicate missing intermediate ranks, not errors
			ArrayList<Integer> validTaxonIds = new ArrayList<Integer>();
			for (int lineageTaxonId : mainTaxonIds) {
				if (lineageTaxonId > 0) {
					validTaxonIds.add(lineageTaxonId);
				} else if (lineageTaxonId == -1) {
					logger.debug(taskName + ": Missing rank in lineage for taxon " + taxonId + " (sparse lineage)");
				} else if (lineageTaxonId == 0) {
					logger.debug(
							taskName + ": Invalid taxon ID 0 in lineage for taxon " + taxonId + " (data corruption)");
				}
			}

			if (validTaxonIds.isEmpty()) {
				logger.warn(taskName + ": No valid taxon IDs found in lineage for taxon " + taxonId
						+ " (lineage contained only missing ranks or invalid IDs), skipping");

				// Debug this specific taxon to understand the issue
				if (taxonId == 9606) { // Homo sapiens
					logger.warn(taskName + ": Debugging Homo sapiens lineage issue - this should not happen!");
					// Enable detailed debugging to see what's going wrong
					((TaxonomyDatabaseSQL) taxDb).debugTaxonLineage(taxonId);
				}
				continue; // Skip if no valid taxon IDs in lineage
			}

			// Convert back to array for processing
			mainTaxonIds = validTaxonIds.stream().mapToInt(Integer::intValue).toArray();

			// Process each taxonomic rank in the lineage
			for (int rankIndex = 0; rankIndex < mainTaxonIds.length && rankIndex < mainRanks.length; rankIndex++) {
				int parentTaxonId = mainTaxonIds[rankIndex];
				Taxon parentTaxon = taxonInfoMap.get(parentTaxonId);
				if (parentTaxon == null) {
					// Get parent taxon if not already cached (reuse existing taxDb)
					try {
						parentTaxon = taxDb.getTaxonFromId(parentTaxonId);
						if (parentTaxon != null) {
							taxonInfoMap.put(parentTaxonId, parentTaxon);
						}
					} catch (Exception e) {
						continue;
					}
				}

				if (parentTaxon == null) {
					continue;
				}

				String taxonName = parentTaxon.getName();
				if (taxonName == null || taxonName.isEmpty() || taxonName.trim().isEmpty()) {
					// Skip taxa with empty names - don't aggregate intensities under ""
					continue;
				}

				// Build lineage for this rank
				String[] lineage = new String[rankIndex + 1];
				for (int i = 0; i <= rankIndex && i < mainTaxonIds.length; i++) {
					int lineageTaxonId = mainTaxonIds[i];
					if (lineageTaxonId > 0) { // Only process valid taxon IDs
						Taxon lineageTaxon = taxonInfoMap.get(lineageTaxonId);
						if (lineageTaxon != null) {
							String lineageTaxonName = lineageTaxon.getName();
							if (lineageTaxonName != null && !lineageTaxonName.isEmpty()) {
								lineage[i] = lineageTaxonName;
							} else {
								lineage[i] = ""; // Missing rank name
							}
						} else {
							lineage[i] = ""; // Missing taxon
						}
					} else {
						lineage[i] = ""; // Missing rank (-1)
					}
				}

				if (!taxRankMaps[rankIndex].containsKey(taxonName)) {
					taxRankMaps[rankIndex].put(taxonName, lineage);
					taxIntensityMaps[rankIndex].put(taxonName, new double[proteinIntensity.length]);
					taxPsmCountMaps[rankIndex].put(taxonName, 0);
					taxPepCountMaps[rankIndex].put(taxonName, new HashSet<String>());
					taxProteinCountMaps[rankIndex].put(taxonName, 0);
				}

				// Aggregate intensity (like MAG workflow line 2212-2214)
				double[] taxIntensity = taxIntensityMaps[rankIndex].get(taxonName);
				for (int i = 0; i < proteinIntensity.length; i++) {
					taxIntensity[i] += proteinIntensity[i];
				}

				// Aggregate counts
				taxPsmCountMaps[rankIndex].put(taxonName, taxPsmCountMaps[rankIndex].get(taxonName) + psmCount);
				taxPepCountMaps[rankIndex].get(taxonName).addAll(pepCount);
				taxProteinCountMaps[rankIndex].put(taxonName, taxProteinCountMaps[rankIndex].get(taxonName) + 1);
			}
		}

		// Write taxonomy reports (similar structure to MAG workflow)
		try {
			// Write taxa.tsv
			File taxaFile = new File(taxFile, "Taxa.tsv");
			try (PrintWriter writer = new PrintWriter(taxaFile)) {
				StringBuilder header = new StringBuilder();
				header.append("Name\t");
				header.append("Rank");
				for (int i = 0; i < mainRanks.length; i++) {
					header.append("\t").append(mainRanks[i].getName());
				}
				header.append("\tProtein_count\tPSM_count\tPeptide_count");
				for (String intensity : intensityTitle) {
					if (intensity.startsWith("Intensity ")) {
						header.append("\t").append(intensity);
					} else {
						header.append("\tIntensity ").append(intensity);
					}
				}
				writer.println(header);

				// Write data for each rank (like MAG workflow)
				for (int i = 0; i < taxRankMaps.length; i++) {
					for (String taxonName : taxRankMaps[i].keySet()) {
						StringBuilder sb = new StringBuilder();
						sb.append(taxonName).append("\t");
						sb.append(mainRanks[i].getName());

						String[] lineage = taxRankMaps[i].get(taxonName);
						for (int j = 0; j < mainRanks.length; j++) {
							if (j < lineage.length) {
								sb.append("\t").append(lineage[j]);
							} else {
								sb.append("\t");
							}
						}

						sb.append("\t").append(taxProteinCountMaps[i].get(taxonName));
						sb.append("\t").append(taxPsmCountMaps[i].get(taxonName));
						sb.append("\t").append(taxPepCountMaps[i].get(taxonName).size());

						double[] intensity = taxIntensityMaps[i].get(taxonName);
						for (double val : intensity) {
							sb.append("\t").append(val);
						}

						writer.println(sb);
					}
				}
			}

			// Export species-only taxonomy table with Selection_label from razorSpecies.tsv
			HashMap<String, String> speciesSelectionLabelMap = new HashMap<String, String>();
			ArrayList<String> speciesRankedOrder = new ArrayList<String>();
			File razorSpeciesFile = new File(this.combinedFile, "razorSpecies.tsv");
			if (razorSpeciesFile.exists()) {
				try (BufferedReader reader = new BufferedReader(new FileReader(razorSpeciesFile))) {
					String header = reader.readLine();
					if (header != null) {
						String[] hs = header.split("\t");
						int itemIdx = -1;
						int labelIdx = -1;
						for (int i = 0; i < hs.length; i++) {
							if ("Item".equals(hs[i])) {
								itemIdx = i;
							} else if ("Selection_label".equals(hs[i])) {
								labelIdx = i;
							}
						}
						if (itemIdx >= 0 && labelIdx >= 0) {
							String line;
							while ((line = reader.readLine()) != null) {
								String[] cs = line.split("\t");
								if (cs.length > Math.max(itemIdx, labelIdx)) {
									speciesSelectionLabelMap.put(cs[itemIdx], cs[labelIdx]);
									speciesRankedOrder.add(cs[itemIdx]);
								}
							}
						}
					}
				} catch (Exception e) {
					logger.warn(taskName + ": failed to parse species labels from " + razorSpeciesFile, e);
				}
			}

			HashMap<String, double[]> speciesIntensityMap = new HashMap<String, double[]>();
			HashMap<String, Integer> speciesPsmCountMap = new HashMap<String, Integer>();
			HashMap<String, HashSet<String>> speciesPepCountMap = new HashMap<String, HashSet<String>>();
			HashMap<String, Integer> speciesProteinCountMap = new HashMap<String, Integer>();
			HashMap<String, Double> speciesTotalIntenMap = new HashMap<String, Double>();

			for (String protein : proteins) {
				if (!proteinToTaxonMap.containsKey(protein)) {
					continue;
				}
				String speciesTaxId = String.valueOf(proteinToTaxonMap.get(protein));
				double[] proteinIntensity = proteinIntensityMap.get(protein);
				if (proteinIntensity == null) {
					continue;
				}

				if (!speciesIntensityMap.containsKey(speciesTaxId)) {
					speciesIntensityMap.put(speciesTaxId, new double[proteinIntensity.length]);
				}
				double[] speciesIntensity = speciesIntensityMap.get(speciesTaxId);
				double total = speciesTotalIntenMap.containsKey(speciesTaxId) ? speciesTotalIntenMap.get(speciesTaxId) : 0.0;
				for (int i = 0; i < proteinIntensity.length; i++) {
					speciesIntensity[i] += proteinIntensity[i];
					total += proteinIntensity[i];
				}
				speciesTotalIntenMap.put(speciesTaxId, total);

				int psmCount = proteinPsmCountMap.containsKey(protein) ? proteinPsmCountMap.get(protein) : 0;
				if (speciesPsmCountMap.containsKey(speciesTaxId)) {
					speciesPsmCountMap.put(speciesTaxId, speciesPsmCountMap.get(speciesTaxId) + psmCount);
				} else {
					speciesPsmCountMap.put(speciesTaxId, psmCount);
				}

				if (!speciesPepCountMap.containsKey(speciesTaxId)) {
					speciesPepCountMap.put(speciesTaxId, new HashSet<String>());
				}
				if (totalPepCountMap.containsKey(protein)) {
					speciesPepCountMap.get(speciesTaxId).addAll(totalPepCountMap.get(protein));
				}

				if (speciesProteinCountMap.containsKey(speciesTaxId)) {
					speciesProteinCountMap.put(speciesTaxId, speciesProteinCountMap.get(speciesTaxId) + 1);
				} else {
					speciesProteinCountMap.put(speciesTaxId, 1);
				}
			}

			// Order species for export: core -> mid -> tail, then by descending total intensity.
			// Prefer the ranked order from razorSpecies.tsv when available for stable output.
			HashSet<String> seen = new HashSet<String>();
			ArrayList<String> speciesOrder = new ArrayList<String>();
			for (String s : speciesRankedOrder) {
				if (speciesTotalIntenMap.containsKey(s) && seen.add(s)) {
					speciesOrder.add(s);
				}
			}
			for (String s : speciesTotalIntenMap.keySet()) {
				if (seen.add(s)) {
					speciesOrder.add(s);
				}
			}
			Collections.sort(speciesOrder, (s1, s2) -> {
				String l1 = speciesSelectionLabelMap.getOrDefault(s1, "tail");
				String l2 = speciesSelectionLabelMap.getOrDefault(s2, "tail");
				int p1 = "core".equals(l1) ? 0 : ("mid".equals(l1) ? 1 : 2);
				int p2 = "core".equals(l2) ? 0 : ("mid".equals(l2) ? 1 : 2);
				if (p1 != p2) {
					return p1 - p2;
				}
				return Double.compare(speciesTotalIntenMap.get(s2), speciesTotalIntenMap.get(s1));
			});

			File speciesFile = new File(taxFile, "Species.tsv");
			try (PrintWriter writer = new PrintWriter(speciesFile)) {
				StringBuilder header = new StringBuilder();
				header.append("Species_taxId\tSelection_label\tProtein_count\tPSM_count\tPeptide_count\tRazor_peptide_count\tNCBI_Name\tNCBI_Rank");
				for (int i = 0; i < mainRanks.length; i++) {
					header.append("\tNCBI_").append(mainRanks[i].getName());
				}
				for (String intensity : intensityTitle) {
					if (intensity.startsWith("Intensity ")) {
						header.append("\t").append(intensity);
					} else {
						header.append("\tIntensity ").append(intensity);
					}
				}
				writer.println(header);

				for (String speciesTaxId : speciesOrder) {
					int taxonId;
					try {
						taxonId = Integer.parseInt(speciesTaxId);
					} catch (NumberFormatException e) {
						continue;
					}

					Taxon speciesTaxon = taxonInfoMap.get(taxonId);
					if (speciesTaxon == null) {
						try {
							speciesTaxon = taxDb.getTaxonFromId(taxonId);
							if (speciesTaxon != null) {
								taxonInfoMap.put(taxonId, speciesTaxon);
							}
						} catch (Exception e) {
							// keep null
						}
					}

					String ncbiName = speciesTaxon != null && speciesTaxon.getName() != null ? speciesTaxon.getName()
							: speciesTaxId;
					String ncbiRank = "species";

					String[] lineage = new String[mainRanks.length];
					Arrays.fill(lineage, "");
					if (speciesTaxon != null) {
						if (speciesTaxon.getMainParentIds() == null) {
							try {
								speciesTaxon.setMainParentIds(taxDb.getMainParentTaxonIds(taxonId));
							} catch (Exception e) {
								// ignore missing lineage
							}
						}
						int[] mainTaxonIds = speciesTaxon.getMainParentIds();
						if (mainTaxonIds != null) {
							for (int i = 0; i < mainTaxonIds.length && i < mainRanks.length; i++) {
								if (mainTaxonIds[i] > 0) {
									Taxon parentTaxon = taxonInfoMap.get(mainTaxonIds[i]);
									if (parentTaxon == null) {
										try {
											parentTaxon = taxDb.getTaxonFromId(mainTaxonIds[i]);
											if (parentTaxon != null) {
												taxonInfoMap.put(mainTaxonIds[i], parentTaxon);
											}
										} catch (Exception e) {
											// ignore
										}
									}
									if (parentTaxon != null && parentTaxon.getName() != null) {
										lineage[i] = parentTaxon.getName();
									}
								}
							}
						}
					}

					StringBuilder sb = new StringBuilder();
					sb.append(speciesTaxId).append("\t");
					sb.append(speciesSelectionLabelMap.getOrDefault(speciesTaxId, "tail")).append("\t");
					sb.append(speciesProteinCountMap.getOrDefault(speciesTaxId, 0)).append("\t");
					sb.append(speciesPsmCountMap.getOrDefault(speciesTaxId, 0)).append("\t");
					sb.append(speciesPepCountMap.containsKey(speciesTaxId) ? speciesPepCountMap.get(speciesTaxId).size() : 0)
							.append("\t");
					sb.append(razorPepCountMap.containsKey(speciesTaxId) ? razorPepCountMap.get(speciesTaxId).size() : 0)
							.append("\t");
					sb.append(ncbiName).append("\t");
					sb.append(ncbiRank);
					for (int i = 0; i < lineage.length; i++) {
						sb.append("\t").append(lineage[i]);
					}

					double[] intensity = speciesIntensityMap.get(speciesTaxId);
					if (intensity != null) {
						for (double val : intensity) {
							sb.append("\t").append(val);
						}
					}

					writer.println(sb);
				}
			}
			logger.info(taskName + ": species taxonomy export completed - wrote " + speciesOrder.size()
					+ " species to " + speciesFile.getName());

			logger.info(taskName + ": taxonomy analysis completed - wrote " + taxonToProteins.size()
					+ " taxonomic IDs to " + taxaFile.getName());
			System.out.println(format.format(new Date()) + "\t" + taskName + ": taxonomy analysis completed");

		} catch (IOException e) {
			logger.error(taskName + ": error in writing taxonomy information", e);
			System.err.println(format.format(new Date()) + "\t" + taskName + ": error in writing taxonomy information");
			return false;
		}

		return true;
	}

	protected boolean generateFinalReports(String taskName) {

		updateMessage(taskName + ": exporting report");

		// Copy report templates in background
		if (!copyReportTemplatesInBackground()) {
			logger.error(getTaskName() + ": Failed to copy report templates");
			return false;
		}

		logger.info(taskName + ": exporting report started");
		System.out.println(format.format(new Date()) + "\t" + taskName + ": exporting report started");

		// Report templates are already copied in the main call() method
		// Just ensure the report directory exists and set the paths
		File reportDir = new File(this.metaPar.getDbSearchResultFile(), "report");
		if (!reportDir.exists()) {
			reportDir.mkdirs();
		}

		metaPar.setReportDir(reportDir);

		File reportHtmlDir = new File(reportDir + "\\reports");
		if (!reportHtmlDir.exists()) {
			reportHtmlDir.mkdirs();
		}

		// Export metadata
		File summaryMetaFile = new File(this.metaPar.getDbSearchResultFile(), "metainfo.tsv");
		try {
			metaPar.getMetadata().exportMetadata(summaryMetaFile);
		} catch (FileNotFoundException e) {
			logger.error(taskName + ": error in exporting metadata", e);
		}

	// Generate summary report
	MetaPythonReportTask reportTask = new MetaPythonReportTask(metaPar.getThreadCount(), mrp, logger);
	reportTask.setWorkflowType("IGC");
	reportTask.initializeManifestGenerator(reportHtmlDir);

		// Create initial manifest file to ensure index page can load
		ReportManifestGenerator manifestGen = new ReportManifestGenerator(reportHtmlDir);
		manifestGen.createInitialManifest();

		// Add individual report tasks
		boolean reportSuccess = true;

		// Add summary report
		File summaryReportFile = new File(this.metaPar.getDbSearchResultFile(), "final_summary.tsv");
		if (summaryReportFile.exists() && summaryReportFile.length() > 0) {
			File report_ID_summary = new File(reportHtmlDir, "report_ID_summary.html");
			if (summaryMetaFile != null && summaryMetaFile.exists()) {
				reportSuccess &= reportTask.addTask(MetaPythonReportTask.ID_SUMMARY, summaryReportFile,
						report_ID_summary, summaryMetaFile);
			} else {
				reportSuccess &= reportTask.addTask(MetaPythonReportTask.ID_SUMMARY, summaryReportFile,
						report_ID_summary);
			}
		}

		// Add peptide report
		if (final_pep_txt.exists() && final_pep_txt.length() > 0) {
			File report_peptides_summary = new File(reportHtmlDir, "report_peptides_summary.html");
			if (summaryMetaFile != null && summaryMetaFile.exists()) {
				reportSuccess &= reportTask.addTask(MetaPythonReportTask.PEPTIDE, final_pep_txt,
						report_peptides_summary, summaryMetaFile);
			} else {
				reportSuccess &= reportTask.addTask(MetaPythonReportTask.PEPTIDE, final_pep_txt,
						report_peptides_summary);
			}
		}

		// Add protein report
		if (final_pro_txt.exists() && final_pro_txt.length() > 0) {
			File report_proteinGroups_summary = new File(reportHtmlDir, "report_proteinGroups_summary.html");
			if (summaryMetaFile != null && summaryMetaFile.exists()) {
				reportSuccess &= reportTask.addTask(MetaPythonReportTask.PROTEIN, final_pro_txt,
						report_proteinGroups_summary, summaryMetaFile);
			} else {
				reportSuccess &= reportTask.addTask(MetaPythonReportTask.PROTEIN, final_pro_txt,
						report_proteinGroups_summary);
			}
		}

		// Add taxonomy report
		File taxaReportFile = new File(metaPar.getDbSearchResultFile(), "taxonomy_analysis\\Taxa.tsv");
		if (taxaReportFile.exists() && taxaReportFile.length() > 0) {
			File report_taxonomy_summary = new File(reportHtmlDir, "report_taxonomy_summary.html");
			if (summaryMetaFile != null && summaryMetaFile.exists()) {
				reportSuccess &= reportTask.addTask(MetaPythonReportTask.TAXONOMY, taxaReportFile,
						report_taxonomy_summary, summaryMetaFile);
			} else {
				reportSuccess &= reportTask.addTask(MetaPythonReportTask.TAXONOMY, taxaReportFile,
						report_taxonomy_summary);
			}
		}

		// Generate tree.js and bubble chart for IGC workflow
		File jsOutput = new File(reportHtmlDir, "tree.js");
		if (!jsOutput.exists() || jsOutput.length() == 0) {
			// Generate tree.js from taxonomy data if available
			generateTreeJsForIgc(jsOutput);
		}

		// Generate bubble chart report from tree.js data (optional - don't fail if this
		// doesn't work)
		File bubbleChartOutput = new File(reportHtmlDir, "bubble_chart_report.html");
		if (!bubbleChartOutput.exists() || bubbleChartOutput.length() == 0) {
			try {
				generateBubbleChartReport(jsOutput, bubbleChartOutput);
			} catch (Exception e) {
				logger.warn(getTaskName() + ": Bubble chart generation failed, continuing with other reports: "
						+ e.getMessage());
				System.out.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": Bubble chart generation failed, continuing with other reports");
			}
		}

		// Add function report
		File funReportTsv = new File(metaPar.getDbSearchResultFile(), "functional_annotation");
		if (funReportTsv.exists()) {
			File report_function_summary = new File(reportHtmlDir, "report_function_summary.html");
			if (summaryMetaFile != null && summaryMetaFile.exists()) {
				reportSuccess &= reportTask.addTask(MetaPythonReportTask.FUNCTION, funReportTsv,
						report_function_summary, summaryMetaFile);
			} else {
				reportSuccess &= reportTask.addTask(MetaPythonReportTask.FUNCTION, funReportTsv,
						report_function_summary);
			}
		}

		if (!reportSuccess) {
			logger.warn(taskName + ": Some report tasks could not be added");
			System.out.println(format.format(new Date()) + "\t" + taskName + ": Some report tasks could not be added");
		}

		// Run the report generation task
		Thread reportThread = new Thread(reportTask);
		reportThread.start();
		try {
			reportThread.join();
		} catch (InterruptedException e) {
			logger.error(taskName + ": Report generation thread was interrupted", e);
			System.out.println(
					format.format(new Date()) + "\t" + taskName + ": Report generation thread was interrupted");
			Thread.currentThread().interrupt();
			return false;
		}

		try {
			boolean finish = reportTask.get();
			if (finish) {
				logger.info(taskName + ": report generation successes");
				System.out.println(format.format(new Date()) + "\t" + taskName + ": report generation successes");
				return true;
			}
		} catch (InterruptedException | ExecutionException e) {
			logger.error(taskName + ": error in the report generation", e);
			System.out.println(format.format(new Date()) + "\t" + taskName + ": error in the report generation");
			return false;
		}

		return true;
	}

	/**
	 * Generate tree.js file for IGC workflow from taxonomy data.
	 */
	private void generateTreeJsForIgc(File jsOutput) {
		try {
			logger.info(getTaskName() + ": Generating tree.js for IGC workflow");
			System.out.println(
					format.format(new Date()) + "\t" + getTaskName() + ": Generating tree.js for IGC workflow");

			// For IGC workflow, we'll create a simple tree.js structure
			// This is a placeholder - in a real implementation, you'd parse the taxonomy
			// TSV data
			String treeJsContent = "let data = '[]'; // IGC workflow - no hierarchical tree data available";

			try (PrintWriter writer = new PrintWriter(jsOutput)) {
				writer.println(treeJsContent);
			}

			logger.info(getTaskName() + ": Generated tree.js for IGC workflow: " + jsOutput.getName());

		} catch (Exception e) {
			logger.error(getTaskName() + ": Error generating tree.js for IGC workflow", e);
			System.err.println(
					format.format(new Date()) + "\t" + getTaskName() + ": Error generating tree.js: " + e.getMessage());
		}
	}

	/**
	 * Generate bubble chart report from tree.js data using Python.
	 */
	private void generateBubbleChartReport(File treeJsFile, File outputFile) {
		try {
			logger.info(getTaskName() + ": Generating bubble chart report from " + treeJsFile.getName());
			System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": Generating bubble chart report");

		// Use MetaPythonReportTask to generate the bubble chart
		MetaPythonReportTask reportTask = new MetaPythonReportTask(1, mrp, logger);
		boolean success = reportTask.addTask(MetaPythonReportTask.BUBBLE_CHART, treeJsFile, outputFile);

			if (success) {
				// Wait for completion
				reportTask.run();
				boolean result = reportTask.get();

				if (result) {
					logger.info(
							getTaskName() + ": Bubble chart report generated successfully: " + outputFile.getName());
					System.out.println(format.format(new Date()) + "\t" + getTaskName()
							+ ": Bubble chart report generated successfully");

					// Add to manifest
					File reportsDir = new File(metaPar.getReportDir(), "reports");
					ReportManifestGenerator manifestGen = new ReportManifestGenerator(reportsDir);
					manifestGen.addReportEntry("Bubble Chart Report", "bubble_chart_report.html",
							"Interactive hierarchical bubble chart showing taxonomic abundance with treemap, sunburst, and distribution visualizations.",
							new String[] { "bubble-chart", "taxonomy", "hierarchical", "treemap" });
				} else {
					logger.error(getTaskName() + ": Failed to generate bubble chart report");
					System.err.println(format.format(new Date()) + "\t" + getTaskName()
							+ ": Failed to generate bubble chart report");
				}
			} else {
				logger.error(getTaskName() + ": Failed to add bubble chart task");
				System.err.println(
						format.format(new Date()) + "\t" + getTaskName() + ": Failed to add bubble chart task");
			}
		} catch (Exception e) {
			logger.error(getTaskName() + ": Error generating bubble chart report", e);
			System.err.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": Error generating bubble chart report: " + e.getMessage());
		}
	}

	/**
	 * Save protinfo batch results to a temporary JSON file.
	 * 
	 * @param batchResult The batch results to save
	 * @param tempFile    The temporary file to save to
	 * @return true if successful, false otherwise
	 */
	private boolean saveProtinfoBatchToTemp(HashMap<String, UnipProInfo> batchResult, File tempFile) {
		try {
			// Create root JSON object
			org.json.JSONObject rootObject = new org.json.JSONObject();

			// Add each protein ID as a key with its UnipProInfo data
			for (Map.Entry<String, UnipProInfo> entry : batchResult.entrySet()) {
				String proteinId = entry.getKey();
				UnipProInfo info = entry.getValue();

				// Create JSON object for this protein
				org.json.JSONObject proteinObject = new org.json.JSONObject();
				proteinObject.put("pro", info.getPro() != null ? info.getPro() : "");
				proteinObject.put("name", info.getName() != null ? info.getName() : "");
				proteinObject.put("taxId", info.getTaxId());

				// Add EC array
				org.json.JSONArray ecArray = new org.json.JSONArray();
				String[] ec = info.getEc();
				if (ec != null) {
					for (String ecTerm : ec) {
						ecArray.put(ecTerm);
					}
				}
				proteinObject.put("ec", ecArray);

				// Add GO array
				org.json.JSONArray goArray = new org.json.JSONArray();
				String[] go = info.getGo();
				if (go != null) {
					for (String goTerm : go) {
						goArray.put(goTerm);
					}
				}
				proteinObject.put("go", goArray);

				// Add to root object
				rootObject.put(proteinId, proteinObject);
			}

			// Write JSON to file
			java.io.FileWriter writer = new java.io.FileWriter(tempFile);
			writer.write(rootObject.toString(2)); // Pretty print with 2-space indentation
			writer.close();

			return true;
		} catch (Exception e) {
			logger.error(getTaskName() + ": Failed to save protinfo batch to temp file: " + e.getMessage());
			return false;
		}
	}

	/**
	 * Load protinfo batch results from a temporary JSON file.
	 * 
	 * @param tempFile The temporary file to load from
	 * @return The loaded batch results, or null if failed
	 */
	private HashMap<String, UnipProInfo> loadProtinfoBatchFromTemp(File tempFile) {
		try {
			HashMap<String, UnipProInfo> result = new HashMap<String, UnipProInfo>();

			// Read JSON content
			String content = new String(java.nio.file.Files.readAllBytes(tempFile.toPath()));

			// Parse JSON using proper JSON library
			org.json.JSONObject jsonObject = new org.json.JSONObject(content);

			// Iterate through all keys (protein IDs)
			for (String proteinId : jsonObject.keySet()) {
				org.json.JSONObject proteinData = jsonObject.getJSONObject(proteinId);

				// Extract data from JSON object
				String pro = proteinData.optString("pro", "");
				String name = proteinData.optString("name", "");
				int taxId = proteinData.optInt("taxId", 0);

				// Parse EC array
				org.json.JSONArray ecArray = proteinData.optJSONArray("ec");
				String[] ec = new String[0];
				if (ecArray != null) {
					ec = new String[ecArray.length()];
					for (int i = 0; i < ecArray.length(); i++) {
						ec[i] = ecArray.getString(i);
					}
				}

				// Parse GO array
				org.json.JSONArray goArray = proteinData.optJSONArray("go");
				String[] go = new String[0];
				if (goArray != null) {
					go = new String[goArray.length()];
					for (int i = 0; i < goArray.length(); i++) {
						go[i] = goArray.getString(i);
					}
				}

				// Create UnipProInfo object
				UnipProInfo info = new UnipProInfo(pro, name, taxId, ec, go);
				result.put(proteinId, info);
			}

			return result;
		} catch (Exception e) {
			logger.error(getTaskName() + ": Failed to load protinfo batch from temp file: " + e.getMessage());
			return null;
		}
	}

	/**
	 * Save pept2prot batch results to a temporary JSON file.
	 * 
	 * @param batchResult The batch results to save
	 * @param tempFile    The temporary file to save to
	 * @return true if successful, false otherwise
	 */
	private boolean savePept2protBatchToTemp(HashMap<String, UnipProResult> batchResult, File tempFile) {
		try {
			// Create root JSON object
			org.json.JSONObject rootObject = new org.json.JSONObject();

			// Add each peptide sequence as a key with its UnipProResult data
			for (Map.Entry<String, UnipProResult> entry : batchResult.entrySet()) {
				String sequence = entry.getKey();
				UnipProResult result = entry.getValue();

				// Create JSON object for this sequence
				org.json.JSONObject sequenceObject = new org.json.JSONObject();
				sequenceObject.put("sequence", result.getSequence() != null ? result.getSequence() : "");

				// Add uniprotIds array
				org.json.JSONArray uniprotIdsArray = new org.json.JSONArray();
				String[] uniprotIds = result.getUniprotIds();
				if (uniprotIds != null) {
					for (String uniprotId : uniprotIds) {
						uniprotIdsArray.put(uniprotId);
					}
				}
				sequenceObject.put("uniprotIds", uniprotIdsArray);

				// Add proNames array
				org.json.JSONArray proNamesArray = new org.json.JSONArray();
				String[] proNames = result.getProNames();
				if (proNames != null) {
					for (String proName : proNames) {
						proNamesArray.put(proName);
					}
				}
				sequenceObject.put("proNames", proNamesArray);

				// Add taxonIds array
				org.json.JSONArray taxonIdsArray = new org.json.JSONArray();
				int[] taxonIds = result.getTaxonIds();
				if (taxonIds != null) {
					for (int taxonId : taxonIds) {
						taxonIdsArray.put(taxonId);
					}
				}
				sequenceObject.put("taxonIds", taxonIdsArray);

				// Add to root object
				rootObject.put(sequence, sequenceObject);
			}

			// Write JSON to file
			java.io.FileWriter writer = new java.io.FileWriter(tempFile);
			writer.write(rootObject.toString(2)); // Pretty print with 2-space indentation
			writer.close();

			return true;
		} catch (Exception e) {
			logger.error(getTaskName() + ": Failed to save pept2prot batch to temp file: " + e.getMessage());
			return false;
		}
	}

	/**
	 * Load pept2prot batch results from a temporary JSON file.
	 * 
	 * @param tempFile The temporary file to load from
	 * @return The loaded batch results, or null if failed
	 */
	private HashMap<String, UnipProResult> loadPept2protBatchFromTemp(File tempFile) {
		try {
			HashMap<String, UnipProResult> result = new HashMap<String, UnipProResult>();

			// Read JSON content
			String content = new String(java.nio.file.Files.readAllBytes(tempFile.toPath()));

			// Parse JSON using proper JSON library
			org.json.JSONObject jsonObject = new org.json.JSONObject(content);

			// Iterate through all keys (peptide sequences)
			for (String sequence : jsonObject.keySet()) {
				org.json.JSONObject sequenceData = jsonObject.getJSONObject(sequence);

				// Extract data from JSON object
				String seq = sequenceData.optString("sequence", "");

				// Parse uniprotIds array
				org.json.JSONArray uniprotIdsArray = sequenceData.optJSONArray("uniprotIds");
				String[] uniprotIds = new String[0];
				if (uniprotIdsArray != null) {
					uniprotIds = new String[uniprotIdsArray.length()];
					for (int i = 0; i < uniprotIdsArray.length(); i++) {
						uniprotIds[i] = uniprotIdsArray.getString(i);
					}
				}

				// Parse proNames array
				org.json.JSONArray proNamesArray = sequenceData.optJSONArray("proNames");
				String[] proNames = new String[0];
				if (proNamesArray != null) {
					proNames = new String[proNamesArray.length()];
					for (int i = 0; i < proNamesArray.length(); i++) {
						proNames[i] = proNamesArray.getString(i);
					}
				}

				// Parse taxonIds array
				org.json.JSONArray taxonIdsArray = sequenceData.optJSONArray("taxonIds");
				int[] taxonIds = new int[0];
				if (taxonIdsArray != null) {
					taxonIds = new int[taxonIdsArray.length()];
					for (int i = 0; i < taxonIdsArray.length(); i++) {
						taxonIds[i] = taxonIdsArray.getInt(i);
					}
				}

				// Create UnipProResult object
				UnipProResult unipProResult = new UnipProResult(seq, uniprotIds, proNames, taxonIds);
				result.put(sequence, unipProResult);
			}

			return result;
		} catch (Exception e) {
			logger.error(getTaskName() + ": Failed to load pept2prot batch from temp file: " + e.getMessage());
			return null;
		}
	}

	/**
	 * Task for post-analysis processing in DDA IGC workflows. Handles protein
	 * inference, filtering, and report generation with proper progress tracking.
	 */
	public static class DDAIgcPostAnalysisTask extends MetaAbstractTask {

		private MetaPeptide[] quanPeps;
		private MetaPSM[] psms;
		private HashMap<String, Integer> fileIdMap;
		private ArrayList<String> fileList;
		private String quanMode;
		private String[] intensityTitle;
		private String taskName;
		private MetaPilotIgcTask parentTask;

		public DDAIgcPostAnalysisTask(MetaBaseParameter metaPar, MetaResourceParameters mrp, ProgressBar progressBar,
				MetaPeptide[] quanPeps, MetaPSM[] psms, HashMap<String, Integer> fileIdMap, ArrayList<String> fileList,
				String quanMode, String[] intensityTitle, String taskName, MetaPilotIgcTask parentTask) {
			super(metaPar, mrp, progressBar, parentTask.getLogger());
			this.quanPeps = quanPeps;
			this.psms = psms;
			this.fileIdMap = fileIdMap;
			this.fileList = fileList;
			this.quanMode = quanMode;
			this.intensityTitle = intensityTitle;
			this.taskName = taskName;
			this.parentTask = parentTask;
		}

		@Override
		protected String getTaskName() {
			return taskName + " Post-Analysis";
		}

		@Override
		protected Boolean call() throws Exception {
			updateMessage("Starting post-analysis processing...");

			// Build PSM maps from the PSM array
			updateMessage("Building PSM maps...");
			updateProgress(0.05, 1.0);

			HashMap<String, MetaPSM> psmMap = new HashMap<String, MetaPSM>();
			HashMap<String, int[]> pepPsmCountMap = new HashMap<String, int[]>();
			HashMap<String, int[]> pepChargeMap = new HashMap<String, int[]>();
			HashMap<String, double[]> pepExpectMap = new HashMap<String, double[]>();

			int[] filePsmCount = new int[fileIdMap.size()];
			HashSet<String>[] fileSequenceSets = new HashSet[fileIdMap.size()];
			for (int i = 0; i < fileSequenceSets.length; i++) {
				fileSequenceSets[i] = new HashSet<String>();
			}

			// Process PSMs and build maps
			for (MetaPSM psm : psms) {
				String sequence = psm.getBaseSequence();
				String fileName = psm.getFileName();
				String modseq = psm.getModSeq();
				int charge = psm.getCharge();
				double PEP = psm.getPEP();
				int fileId = fileIdMap.get(fileName);

				// Build PSM count map per peptide per file
				int[] psmCount;
				if (pepPsmCountMap.containsKey(modseq)) {
					psmCount = pepPsmCountMap.get(modseq);
				} else {
					psmCount = new int[fileIdMap.size()];
					pepPsmCountMap.put(modseq, psmCount);
				}
				if (fileIdMap.containsKey(fileName)) {
					psmCount[fileId]++;
				}

				filePsmCount[fileId]++;
				fileSequenceSets[fileId].add(sequence);

				// Build charge map
				int[] charges;
				if (pepChargeMap.containsKey(modseq)) {
					charges = pepChargeMap.get(modseq);
				} else {
					charges = new int[7];
					pepChargeMap.put(modseq, charges);
				}
				if (charge < charges.length) {
					charges[charge] = 1;
				}

				// Build PEP map
				double[] expect;
				if (pepExpectMap.containsKey(modseq)) {
					expect = pepExpectMap.get(modseq);
				} else {
					expect = new double[fileIdMap.size()];
					Arrays.fill(expect, 1.0);
					pepExpectMap.put(modseq, expect);
				}
				if (PEP < expect[fileId]) {
					expect[fileId] = PEP;
				}

				// Store best PSM for each peptide
				if (psmMap.containsKey(modseq)) {
					if (PEP < psmMap.get(modseq).getPEP()) {
						psmMap.put(modseq, psm);
					}
				} else {
					psmMap.put(modseq, psm);
				}
			}

			// Set individual PEP values for peptides
			for (int i = 0; i < quanPeps.length; i++) {
				String seq = quanPeps[i].getModSeq();
				if (pepExpectMap.containsKey(seq)) {
					quanPeps[i].setIndividualPEP(pepExpectMap.get(seq));
				} else {
					double[] expect = new double[fileIdMap.size()];
					Arrays.fill(expect, 1.0);
					quanPeps[i].setIndividualPEP(expect);
				}
			}

			// Update progress: PSM maps built
			updateMessage("PSM maps built, starting protein inference...");
			updateProgress(0.15, 1.0);

			// Perform protein inference and filtering directly in this task
			HashMap<String, UnipProInfo> proteinInfoMap = performProteinInference(quanPeps);

			if (proteinInfoMap.isEmpty()) {
				updateMessage("No proteins found, skipping report generation");
				return false;
			}

			// Update progress: Protein inference complete
			updateMessage("Protein inference complete, updating peptide protein mappings...");
			updateProgress(0.25, 1.0);

			// Update progress: Peptide mapping complete
			updateMessage("Peptide mapping complete, building PSM maps...");
			updateProgress(0.40, 1.0);

			// Refine genomes and proteins
			HashMap<String, Double> speciesScoreMap = refineSpecies(quanPeps, proteinInfoMap);
			String[] razorProteins = refineProteins(quanPeps, speciesScoreMap, proteinInfoMap);

			// Build protein intensity maps
			HashMap<String, double[]> proteinIntensityMap = new HashMap<String, double[]>();
			HashMap<String, Double> proteinTotalIntenMap = new HashMap<String, Double>();
			HashMap<String, HashSet<String>> totalPepCountMap = new HashMap<String, HashSet<String>>();
			HashMap<String, HashSet<String>> razorPepCountMap = new HashMap<String, HashSet<String>>();
			// All-protein maps (no razor restriction), for additional export
			HashMap<String, double[]> allProteinIntensityMap = new HashMap<String, double[]>();
			HashMap<String, Double> allProteinTotalIntenMap = new HashMap<String, Double>();
			HashMap<String, HashSet<String>> allTotalPepCountMap = new HashMap<String, HashSet<String>>();
			HashSet<String> pepSequenceSet = new HashSet<String>();
			ArrayList<MetaPeptide> peptideList = new ArrayList<MetaPeptide>();

			HashMap<String, Integer> psmCountMap = new HashMap<String, Integer>();
			HashMap<String, Integer> allProteinPsmCountMap = new HashMap<String, Integer>();
			HashMap<String, String> razorProteinMap = new HashMap<String, String>();

			// Process peptides and build protein maps
			for (int i = 0; i < quanPeps.length; i++) {
				if (razorProteins[i] != null) {
					String stripeSequence = quanPeps[i].getSequence();
					pepSequenceSet.add(stripeSequence);
					peptideList.add(quanPeps[i]);
					if (quanPeps[i].getModSeq() != null) {
						razorProteinMap.put(quanPeps[i].getModSeq(), razorProteins[i]);
					}
					if (stripeSequence != null) {
						razorProteinMap.put(stripeSequence, razorProteins[i]);
					}

					// Build razor peptide count map
					if (razorPepCountMap.containsKey(razorProteins[i])) {
						razorPepCountMap.get(razorProteins[i]).add(stripeSequence);
					} else {
						HashSet<String> set = new HashSet<String>();
						set.add(stripeSequence);
						razorPepCountMap.put(razorProteins[i], set);
					}

					// Handle MAG database proteins
					if (proteinInfoMap.containsKey(razorProteins[i])) {
						String species = String.valueOf(proteinInfoMap.get(razorProteins[i]).getTaxId());
						if (razorPepCountMap.containsKey(species)) {
							razorPepCountMap.get(species).add(stripeSequence);
						} else {
							HashSet<String> set = new HashSet<String>();
							set.add(stripeSequence);
							razorPepCountMap.put(species, set);
						}
					}
				}
			}

			// Calculate protein intensities
			for (int i = 0; i < quanPeps.length; i++) {
				String[] proteins = quanPeps[i].getProteins();
				int count = 0;
				for (int j = 0; j < proteins.length; j++) {
					if (razorPepCountMap.containsKey(proteins[j])) {
						count++;
					}
				}

				String sequence = quanPeps[i].getSequence();
				String modSeq = quanPeps[i].getModSeq();
				int totalPsmCount = 0;
				int[] psmCount = pepPsmCountMap.get(modSeq);
				if (psmCount != null) {
					for (int j = 0; j < psmCount.length; j++) {
						totalPsmCount += psmCount[j];
					}
				}

				if (count > 0) {
					double[] pepIntensity = new double[quanPeps[i].getIntensity().length];
					for (int j = 0; j < pepIntensity.length; j++) {
						pepIntensity[j] = quanPeps[i].getIntensity()[j] / (double) count;
					}

					for (int j = 0; j < proteins.length; j++) {
						if (razorPepCountMap.containsKey(proteins[j])) {
							double proteinTotalInten;
							double[] proteinIntensity;
							if (proteinIntensityMap.containsKey(proteins[j])) {
								proteinIntensity = proteinIntensityMap.get(proteins[j]);
								proteinTotalInten = proteinTotalIntenMap.get(proteins[j]);
							} else {
								proteinIntensity = new double[pepIntensity.length];
								proteinIntensityMap.put(proteins[j], proteinIntensity);
								proteinTotalInten = 0.0;
							}
							for (int k = 0; k < pepIntensity.length; k++) {
								proteinIntensity[k] += pepIntensity[k];
								proteinTotalInten += pepIntensity[k];
							}
							proteinTotalIntenMap.put(proteins[j], proteinTotalInten);

							if (totalPepCountMap.containsKey(proteins[j])) {
								totalPepCountMap.get(proteins[j]).add(sequence);
							} else {
								HashSet<String> set = new HashSet<String>();
								set.add(sequence);
								totalPepCountMap.put(proteins[j], set);
							}
						}

						if (proteinInfoMap.containsKey(proteins[j])) {
							String species = String.valueOf(proteinInfoMap.get(proteins[j]).getTaxId());
							if (speciesScoreMap.containsKey(species)) {
								if (psmCountMap.containsKey(proteins[j])) {
									psmCountMap.put(proteins[j], psmCountMap.get(proteins[j]) + totalPsmCount);
								} else {
									psmCountMap.put(proteins[j], totalPsmCount);
								}

								if (psmCountMap.containsKey(species)) {
									psmCountMap.put(species, psmCountMap.get(species) + totalPsmCount);
								} else {
									psmCountMap.put(species, totalPsmCount);
								}
							}
						}
					}
				}

				// Build all-protein maps (distribute peptide intensity across all mapped
				// proteins for this peptide).
				if (proteins.length > 0) {
					double[] allPepIntensity = new double[quanPeps[i].getIntensity().length];
					for (int j = 0; j < allPepIntensity.length; j++) {
						allPepIntensity[j] = quanPeps[i].getIntensity()[j] / (double) proteins.length;
					}
					for (int j = 0; j < proteins.length; j++) {
						double[] allProteinIntensity;
						double allProteinTotalInten;
						if (allProteinIntensityMap.containsKey(proteins[j])) {
							allProteinIntensity = allProteinIntensityMap.get(proteins[j]);
							allProteinTotalInten = allProteinTotalIntenMap.get(proteins[j]);
						} else {
							allProteinIntensity = new double[allPepIntensity.length];
							allProteinIntensityMap.put(proteins[j], allProteinIntensity);
							allProteinTotalInten = 0.0;
						}
						for (int k = 0; k < allPepIntensity.length; k++) {
							allProteinIntensity[k] += allPepIntensity[k];
							allProteinTotalInten += allPepIntensity[k];
						}
						allProteinTotalIntenMap.put(proteins[j], allProteinTotalInten);

						if (allTotalPepCountMap.containsKey(proteins[j])) {
							allTotalPepCountMap.get(proteins[j]).add(sequence);
						} else {
							HashSet<String> set = new HashSet<String>();
							set.add(sequence);
							allTotalPepCountMap.put(proteins[j], set);
						}

						if (allProteinPsmCountMap.containsKey(proteins[j])) {
							allProteinPsmCountMap.put(proteins[j], allProteinPsmCountMap.get(proteins[j]) + totalPsmCount);
						} else {
							allProteinPsmCountMap.put(proteins[j], totalPsmCount);
						}
					}
				}
			}

			// Update progress: Intensity calculations complete
			updateMessage("Intensity calculations complete, sorting proteins...");
			updateProgress(0.60, 1.0);

			// Sort proteins by total intensity
			String[] proteins = proteinTotalIntenMap.keySet().toArray(String[]::new);
			Arrays.sort(proteins, new Comparator<String>() {
				@Override
				public int compare(String o1, String o2) {
					if (proteinTotalIntenMap.get(o1) > proteinTotalIntenMap.get(o2))
						return -1;
					if (proteinTotalIntenMap.get(o1) < proteinTotalIntenMap.get(o2))
						return 1;
					return 0;
				}
			});

			HashMap<String, Integer> proteinIdMap = new HashMap<String, Integer>();
			for (int i = 0; i < proteins.length; i++) {
				proteinIdMap.put(proteins[i], i);
			}

			// Generate identification summary
			updateMessage("Generating identification summary...");
			updateProgress(0.70, 1.0);
			if (!parentTask.generateIdentificationSummary(fileIdMap, filePsmCount, fileSequenceSets, taskName)) {
				return false;
			}

			// Generate peptide report
			updateMessage("Generating peptide report...");
			updateProgress(0.75, 1.0);
			if (!parentTask.generatePeptideReport(peptideList, psmMap, pepPsmCountMap, pepChargeMap, proteinIdMap,
					fileList, quanMode, intensityTitle, taskName, razorProteinMap)) {
				return false;
			}

			// Generate protein report
			updateMessage("Generating protein report...");
			updateProgress(0.80, 1.0);
			if (!parentTask.generateProteinReport(proteins, proteinIntensityMap, proteinTotalIntenMap, totalPepCountMap,
					razorPepCountMap, psmCountMap, proteinIdMap, quanMode, intensityTitle, taskName)) {
				return false;
			}

			// Generate all-protein report (no razor filtering)
			String[] allProteins = allProteinTotalIntenMap.keySet().toArray(String[]::new);
			Arrays.sort(allProteins, new Comparator<String>() {
				@Override
				public int compare(String o1, String o2) {
					if (allProteinTotalIntenMap.get(o1) > allProteinTotalIntenMap.get(o2))
						return -1;
					if (allProteinTotalIntenMap.get(o1) < allProteinTotalIntenMap.get(o2))
						return 1;
					return 0;
				}
			});
			if (!parentTask.generateAllProteinReport(allProteins, allProteinIntensityMap, allProteinTotalIntenMap,
					allTotalPepCountMap, razorPepCountMap, allProteinPsmCountMap, quanMode, intensityTitle, taskName)) {
				return false;
			}

			// Functional annotation
			updateMessage("Performing functional annotation...");
			updateProgress(0.85, 1.0);
			if (!parentTask.performFunctionalAnnotation(proteins, proteinIntensityMap, proteinTotalIntenMap,
					totalPepCountMap, proteinInfoMap, intensityTitle, taskName)) {
				return false;
			}

			// Taxonomy analysis
			updateMessage("Performing taxonomy analysis...");
			updateProgress(0.90, 1.0);
			if (!parentTask.performTaxonomyAnalysis(proteins, proteinIntensityMap, psmCountMap, totalPepCountMap,
					razorPepCountMap, proteinInfoMap, intensityTitle, taskName)) {
				return false;
			}

			// Generate final reports
			updateMessage("Generating final HTML reports...");
			updateProgress(0.95, 1.0);
			if (!parentTask.generateFinalReports(taskName)) {
				return false;
			}

			// Complete
			updateMessage("Report generation completed!");
			updateProgress(1.0, 1.0);

			return true;
		}

		/**
		 * Performs protein inference and filtering using peptide and PSM data. This
		 * method contains the logic previously in MetaPilotIgcTask.updateProtein().
		 * 
		 * @param peptides Array of MetaPeptide objects
		 * @return Map of filtered protein names to their UnipProInfo objects
		 */
		private HashMap<String, UnipProInfo> performProteinInference(MetaPeptide[] peptides) {
			logger.info(getTaskName() + ": Starting protein inference and filtering using peptide and PSM data");

			// Step 1: Build peptide-to-protein mapping from MetaPeptide data
			HashMap<String, HashSet<String>> peptideToProteins = new HashMap<String, HashSet<String>>();
			HashMap<String, HashSet<String>> proteinToPeptides = new HashMap<String, HashSet<String>>();

			for (MetaPeptide peptide : peptides) {
				if (peptide != null && peptide.getSequence() != null && peptide.getProteins() != null) {
					String sequence = peptide.getSequence();
					String[] peptideProteins = peptide.getProteins();

					if (!peptideToProteins.containsKey(sequence)) {
						peptideToProteins.put(sequence, new HashSet<String>());
					}

					for (String protein : peptideProteins) {
						if (protein != null && !protein.isEmpty()) {
							peptideToProteins.get(sequence).add(protein);

							if (!proteinToPeptides.containsKey(protein)) {
								proteinToPeptides.put(protein, new HashSet<String>());
							}
							proteinToPeptides.get(protein).add(sequence);
						}
					}
				}
			}

			logger.info(getTaskName() + ": Found " + proteinToPeptides.size() + " proteins with "
					+ peptideToProteins.size() + " unique peptides");

			// Step 2: Extract unique peptide sequences for Unipept API call
			HashSet<String> uniqueSequences = new HashSet<String>(peptideToProteins.keySet());
			if (uniqueSequences.isEmpty()) {
				logger.warn(getTaskName() + ": No valid peptide sequences found for protein inference");
				return new HashMap<String, UnipProInfo>();
			}

			String[] sequences = uniqueSequences.toArray(new String[uniqueSequences.size()]);
			logger.info(getTaskName() + ": Calling Unipept API for " + sequences.length + " unique peptide sequences");

			// Step 3: Call UnipeptRequester.pept2prot to get protein mappings with
			// resumable batching
			HashMap<String, UnipProResult> pept2protMap = new HashMap<String, UnipProResult>();

			// Process peptides in batches to avoid server timeout (500 error)
			int batchSize = 1000; // Process 1000 peptides at a time
			int totalBatches = (int) Math.ceil((double) sequences.length / batchSize);

			// Create temp directory for batch results
			File pept2protTempDir = new File(metaPar.getDbSearchResultFile(), "pept2prot_temp");
			if (!pept2protTempDir.exists()) {
				pept2protTempDir.mkdirs();
				logger.info(getTaskName() + ": Created temp directory for pept2prot batches: "
						+ pept2protTempDir.getAbsolutePath());
			}

			int completedBatches = 0;
			int skippedBatches = 0;

			// Update message for pept2prot phase
			updateMessage("Processing pept2prot API requests (" + totalBatches + " batches)...");

			for (int batch = 0; batch < totalBatches; batch++) {
				int startIndex = batch * batchSize;
				int endIndex = Math.min(startIndex + batchSize, sequences.length);

				String[] batchSequences = new String[endIndex - startIndex];
				System.arraycopy(sequences, startIndex, batchSequences, 0, endIndex - startIndex);

				// Create temp file for this batch
				File batchTempFile = new File(pept2protTempDir, "pept2prot_batch_" + batch + ".json");

				// Check if this batch was already completed
				if (batchTempFile.exists() && batchTempFile.length() > 0) {
					logger.info(getTaskName() + ": Skipping pept2prot batch " + (batch + 1) + "/" + totalBatches
							+ " (already completed, temp file exists)");

					// Load existing results from temp file
					try {
						HashMap<String, UnipProResult> existingResult = parentTask
								.loadPept2protBatchFromTemp(batchTempFile);
						if (existingResult != null && !existingResult.isEmpty()) {
							pept2protMap.putAll(existingResult);
							skippedBatches++;
							logger.info(getTaskName() + ": Loaded " + existingResult.size()
									+ " pept2prot mappings from temp file");
						}
					} catch (Exception e) {
						logger.warn(getTaskName() + ": Failed to load temp file for batch " + (batch + 1)
								+ ", will reprocess: " + e.getMessage());
						// Continue to process this batch
					}

					if (batchTempFile.exists()) {
						// Update progress for skipped batch
						int processedBatches = completedBatches + skippedBatches + 1;
						double pept2protProgress = (double) processedBatches / totalBatches * 0.4;
						updateProgress(0.15 + pept2protProgress, 1.0);
						updateMessage("pept2prot: " + processedBatches + "/" + totalBatches
								+ " batches completed (skipped: " + (skippedBatches + 1) + ")");
						continue; // Skip this batch as it's already completed
					}
				}

				logger.info(getTaskName() + ": Processing batch " + (batch + 1) + "/" + totalBatches + " ("
						+ batchSequences.length + " peptides)");

				HashMap<String, UnipProResult> batchResult = UnipeptRequester.pept2prot(batchSequences);

				if (batchResult != null && !batchResult.isEmpty()) {
					// Save batch result to temp file
					if (parentTask.savePept2protBatchToTemp(batchResult, batchTempFile)) {
						pept2protMap.putAll(batchResult);
						completedBatches++;
						logger.info(getTaskName() + ": Batch " + (batch + 1) + " returned " + batchResult.size()
								+ " mappings and saved to temp file");
					} else {
						logger.warn(getTaskName() + ": Failed to save temp file for batch " + (batch + 1)
								+ ", but adding results to map");
						pept2protMap.putAll(batchResult);
						completedBatches++;
					}
				} else {
					logger.warn(getTaskName() + ": Batch " + (batch + 1) + " returned no mappings");
					// Create empty temp file to mark this batch as completed
					try {
						batchTempFile.createNewFile();
						completedBatches++;
					} catch (Exception e) {
						logger.warn(getTaskName() + ": Failed to create empty temp file for batch " + (batch + 1));
					}
				}

				// Update progress for pept2prot phase (0-40% of total protein inference)
				int processedBatches = completedBatches + skippedBatches;
				double pept2protProgress = (double) processedBatches / totalBatches * 0.4; // 40% of protein inference
				updateProgress(0.15 + pept2protProgress, 1.0); // Start from 15%, add pept2prot progress
				updateMessage("pept2prot: " + processedBatches + "/" + totalBatches + " batches completed");

				// Add a small delay between batches to be respectful to the API
				if (batch < totalBatches - 1) {
					try {
						Thread.sleep(1000); // 1 second delay between batches
					} catch (InterruptedException e) {
						logger.warn(getTaskName() + ": Interrupted during batch processing");
						Thread.currentThread().interrupt();
						break;
					}
				}
			}

			// Note: Don't clean up temp files here - keep them for debugging until the
			// entire workflow completes
			// Temp files are preserved for user inspection and debugging
			logger.info(getTaskName() + ": Pept2prot processing complete, temp files preserved for user inspection");

			if (pept2protMap.isEmpty()) {
				logger.warn(
						getTaskName() + ": No protein mappings returned from Unipept API after processing all batches");
				return new HashMap<String, UnipProInfo>();
			}

			logger.info(getTaskName() + ": Received " + pept2protMap.size()
					+ " protein mappings from Unipept API across " + totalBatches + " batches " + "(completed: "
					+ completedBatches + ", skipped: " + skippedBatches + ")");

			// Step 4: Species-level adaptive selection directly from pept2prot taxonIds
			// (before protinfo)
			HashMap<String, HashSet<String>> species2PeptideMap = new HashMap<String, HashSet<String>>();
			HashSet<String> speciesPeptides = new HashSet<String>();
			for (String peptide : peptideToProteins.keySet()) {
				UnipProResult result = pept2protMap.get(peptide);
				if (result == null || result.getTaxonIds() == null) {
					continue;
				}
				int[] taxonIds = result.getTaxonIds();
				HashSet<String> usedTaxa = new HashSet<String>();
				boolean use = false;
				for (int i = 0; i < taxonIds.length; i++) {
					String taxon = String.valueOf(taxonIds[i]);
					if (!usedTaxa.add(taxon)) {
						continue;
					}
					if (!species2PeptideMap.containsKey(taxon)) {
						species2PeptideMap.put(taxon, new HashSet<String>());
					}
					species2PeptideMap.get(taxon).add(peptide);
					use = true;
				}
				if (use) {
					speciesPeptides.add(peptide);
				}
			}

			File razorSpeciesFile = new File(parentTask.combinedFile, "razorSpecies.tsv");
			HashSet<String> speciesSet = parentTask.applyAdaptiveOccamsRazor(species2PeptideMap, speciesPeptides,
					razorSpeciesFile);
			logger.info(getTaskName() + ": Species-level filtering selected " + speciesSet.size() + " species");

			// Step 5: Filter proteins by selected species and deduplicate by protein name
			// per peptide.
			// Keep the first UniProt ID for each distinct protein name after species
			// filtering.
			HashMap<String, String[]> cleanedPeptideProteinMap = new HashMap<String, String[]>();
			HashSet<String> selectedUniprotProteins = new HashSet<String>();
			for (String peptide : peptideToProteins.keySet()) {
				UnipProResult result = pept2protMap.get(peptide);
				if (result == null || result.getUniprotIds() == null || result.getTaxonIds() == null) {
					continue;
				}

				String[] uniprotIds = result.getUniprotIds();
				int[] taxonIds = result.getTaxonIds();
				String[] proNames = result.getProNames();
				int minLen = Math.min(uniprotIds.length, taxonIds.length);
				if (minLen == 0) {
					continue;
				}

				HashMap<String, String> nameToProtein = new HashMap<String, String>();
				ArrayList<String> orderedProteins = new ArrayList<String>();
				for (int i = 0; i < minLen; i++) {
					String taxon = String.valueOf(taxonIds[i]);
					if (!speciesSet.contains(taxon)) {
						continue;
					}
					String uniprotId = uniprotIds[i];
					if (uniprotId == null || uniprotId.trim().isEmpty()) {
						continue;
					}
					String proteinName = uniprotId;
					if (proNames != null && i < proNames.length && proNames[i] != null && !proNames[i].trim().isEmpty()) {
						proteinName = proNames[i].trim();
					}
					if (!nameToProtein.containsKey(proteinName)) {
						nameToProtein.put(proteinName, uniprotId);
						orderedProteins.add(uniprotId);
					}
				}

				if (!orderedProteins.isEmpty()) {
					String[] cleaned = orderedProteins.toArray(String[]::new);
					cleanedPeptideProteinMap.put(peptide, cleaned);
					for (String protein : cleaned) {
						selectedUniprotProteins.add(protein);
					}
				}
			}

			logger.info(getTaskName() + ": Protein-level filtering selected " + selectedUniprotProteins.size()
					+ " Uniprot proteins");

			// Step 6: Get taxonomic information for selected Uniprot proteins using
			// protinfo (after species and per-peptide protein cleanup)
			String[] proteinArray = selectedUniprotProteins.toArray(new String[selectedUniprotProteins.size()]);
			logger.info(getTaskName() + ": Calling Unipept protinfo API for " + proteinArray.length
					+ " unique protein IDs");

			HashMap<String, UnipProInfo> proteinInfoMap = new HashMap<String, UnipProInfo>();

			// Process proteins in batches to avoid server timeout (500 error)
			totalBatches = (int) Math.ceil((double) proteinArray.length / batchSize);

			// Update message for protinfo phase
			updateMessage("Processing protinfo API requests (" + totalBatches + " batches)...");

			// Create temp directory for batch results
			File tempDir = new File(metaPar.getDbSearchResultFile(), "protinfo_temp");
			if (!tempDir.exists()) {
				tempDir.mkdirs();
				logger.info(
						getTaskName() + ": Created temp directory for protinfo batches: " + tempDir.getAbsolutePath());
			}

			completedBatches = 0;
			skippedBatches = 0;

			for (int batch = 0; batch < totalBatches; batch++) {
				int startIndex = batch * batchSize;
				int endIndex = Math.min(startIndex + batchSize, proteinArray.length);

				String[] batchProteinIds = new String[endIndex - startIndex];
				System.arraycopy(proteinArray, startIndex, batchProteinIds, 0, endIndex - startIndex);

				// Create temp file for this batch
				File batchTempFile = new File(tempDir, "protinfo_batch_" + batch + ".json");

				// Check if this batch was already completed
				if (batchTempFile.exists() && batchTempFile.length() > 0) {
					logger.info(getTaskName() + ": Skipping protinfo batch " + (batch + 1) + "/" + totalBatches
							+ " (already completed, temp file exists)");

					// Load existing results from temp file
					try {
						HashMap<String, UnipProInfo> existingResult = parentTask
								.loadProtinfoBatchFromTemp(batchTempFile);
						if (existingResult != null && !existingResult.isEmpty()) {
							proteinInfoMap.putAll(existingResult);
							skippedBatches++;
							logger.info(getTaskName() + ": Loaded " + existingResult.size()
									+ " protein infos from temp file");
						}
					} catch (Exception e) {
						logger.warn(getTaskName() + ": Failed to load temp file for batch " + (batch + 1)
								+ ", will reprocess: " + e.getMessage());
						// Continue to process this batch
					}

					if (batchTempFile.exists()) {
						// Update progress for skipped protinfo batch (55-95% of total protein
						// inference)
						int processedBatches = completedBatches + skippedBatches + 1;
						double protinfoProgress = (double) processedBatches / totalBatches * 0.4; // 40% of protein
																									// inference
						updateProgress(0.55 + protinfoProgress, 1.0); // Start from 55%, add protinfo progress
						updateMessage("protinfo: " + processedBatches + "/" + totalBatches
								+ " batches completed (skipped: " + (skippedBatches + 1) + ")");
						continue; // Skip this batch as it's already completed
					}
				}

				logger.info(getTaskName() + ": Processing protinfo batch " + (batch + 1) + "/" + totalBatches + " ("
						+ batchProteinIds.length + " proteins)");

				HashMap<String, UnipProInfo> batchResult = UnipeptRequester.protinfo(batchProteinIds);

				if (batchResult != null && !batchResult.isEmpty()) {
					// Save batch result to temp file
					if (parentTask.saveProtinfoBatchToTemp(batchResult, batchTempFile)) {
						proteinInfoMap.putAll(batchResult);
						completedBatches++;
						logger.info(getTaskName() + ": Protinfo batch " + (batch + 1) + " returned "
								+ batchResult.size() + " protein infos and saved to temp file");
					} else {
						logger.warn(getTaskName() + ": Failed to save temp file for batch " + (batch + 1)
								+ ", but adding results to map");
						proteinInfoMap.putAll(batchResult);
						completedBatches++;
					}
				} else {
					logger.warn(getTaskName() + ": Protinfo batch " + (batch + 1) + " returned no protein infos");
					// Create empty temp file to mark this batch as completed
					try {
						batchTempFile.createNewFile();
						completedBatches++;
					} catch (Exception e) {
						logger.warn(getTaskName() + ": Failed to create empty temp file for batch " + (batch + 1));
					}
				}

				// Update progress for protinfo phase (55-95% of total protein inference)
				int processedBatches = completedBatches + skippedBatches;
				double protinfoProgress = (double) processedBatches / totalBatches * 0.4; // 40% of protein inference
				updateProgress(0.55 + protinfoProgress, 1.0); // Start from 55%, add protinfo progress
				updateMessage("protinfo: " + processedBatches + "/" + totalBatches + " batches completed");

				// Add a small delay between batches to be respectful to the API
				if (batch < totalBatches - 1) {
					try {
						Thread.sleep(1000); // 1 second delay between batches
					} catch (InterruptedException e) {
						logger.warn(getTaskName() + ": Interrupted during protinfo batch processing");
						Thread.currentThread().interrupt();
						break;
					}
				}
			}

			// Note: Don't clean up temp files here - keep them for debugging until the
			// entire workflow completes
			// Temp files are preserved for user inspection and debugging
			logger.info(getTaskName() + ": Protinfo processing complete, temp files preserved for user inspection");

			if (proteinInfoMap.isEmpty()) {
				logger.warn(getTaskName()
						+ ": No protein information returned from Unipept protinfo API after processing all batches");
				return new HashMap<String, UnipProInfo>();
			}

			logger.info(getTaskName() + ": Successfully retrieved information for " + proteinInfoMap.size()
					+ " proteins from Unipept protinfo API across " + totalBatches + " batches " + "(completed: "
					+ completedBatches + ", skipped: " + skippedBatches + ")");

			HashMap<String, HashSet<String>> protein2PeptideMap = new HashMap<String, HashSet<String>>();
			for (String key : peptideToProteins.keySet()) {
				if (cleanedPeptideProteinMap.containsKey(key)) {
					String[] proIds = cleanedPeptideProteinMap.get(key);
					for (String proId : proIds) {
						if (proteinInfoMap.containsKey(proId)) {
							if (protein2PeptideMap.containsKey(proId)) {
								protein2PeptideMap.get(proId).add(key);
							} else {
								HashSet<String> pepSet = new HashSet<String>();
								pepSet.add(key);
								protein2PeptideMap.put(proId, pepSet);
							}
						}
					}
				}
			}

			// Step 7: Protein-level set after species filtering:
			// keep all proteins that map to selected species and survived protinfo lookup.
			HashSet<String> finalProteinSet = new HashSet<>(protein2PeptideMap.keySet());
			logger.info(getTaskName() + ": Protein-level filtering selected " + finalProteinSet.size()
					+ " Uniprot proteins");

			// Step 12: Build final protein info map with only filtered proteins
			HashMap<String, UnipProInfo> finalProteinInfoMap = new HashMap<String, UnipProInfo>();
			for (String protein : finalProteinSet) {
				UnipProInfo info = proteinInfoMap.get(protein);
				if (info != null) {
					finalProteinInfoMap.put(protein, info);
				}
			}

			int updatePeptideCount = 0;
			for (MetaPeptide peptide : peptides) {
				String seq = peptide.getSequence();
				if (cleanedPeptideProteinMap.containsKey(seq)) {
					ArrayList<String> prolist = new ArrayList<String>();
					String[] proIds = cleanedPeptideProteinMap.get(seq);
					for (String proId : proIds) {
						if (finalProteinInfoMap.containsKey(proId)) {
							prolist.add(proId);
						}
					}

					if (prolist.size() > 0) {
						String[] pros = prolist.toArray(String[]::new);
						peptide.setProteins(pros);
						updatePeptideCount++;
					}
				}
			}

			logger.info(getTaskName() + ": Protein inference complete - returning " + finalProteinInfoMap.size()
					+ " filtered proteins with taxonomic information, " + updatePeptideCount
					+ " peptides have been updated");

			return finalProteinInfoMap;
		}

		/**
		 * Refines genomes based on peptide intensity and FDR calculation. Uses default
		 * decoy symbol (PFindResultReader.REV). Reverser proteins are also included to
		 * calculate the FDR of genomes.
		 * 
		 * USED BY: - DDAMagPfindTask: Line 1330 (genome scoring for protein inference)
		 * - DDAMagSageTask: Line 1237 (genome scoring for protein inference)
		 * 
		 * @param quanPeps Array of quantified peptides
		 * @return Map of genome names to their FDR scores
		 */
		protected HashMap<String, Double> refineSpecies(MetaPeptide[] quanPeps,
				HashMap<String, UnipProInfo> uniProInfoMap) {
			return refineSpecies(quanPeps, uniProInfoMap, PFindResultReader.REV);
		}

		/**
		 * Refines genomes based on peptide intensity and FDR calculation with custom
		 * decoy symbol. Reverser proteins are also included to calculate the FDR of
		 * genomes.
		 * 
		 * USED BY: - DDAMagFragPipeTask: Line 1182 (genome scoring with
		 * FastaManager.shuffle decoy symbol)
		 * 
		 * @param quanPeps    Array of quantified peptides
		 * @param decoySymble Custom decoy symbol for FDR calculation
		 * @return Map of genome names to their FDR scores
		 */
		protected HashMap<String, Double> refineSpecies(MetaPeptide[] quanPeps,
				HashMap<String, UnipProInfo> uniProInfoMap, String decoySymble) {

			HashMap<String, Double> speciesIntensityMap = new HashMap<String, Double>();
			double[] peptideIntensity = new double[quanPeps.length];

			for (int i = 0; i < quanPeps.length; i++) {
				double[] pepIntensity = quanPeps[i].getIntensity();
				double totalIntensity = MathTool.getTotal(pepIntensity);
				int[] idenType = quanPeps[i].getIdenType();
				int idenCount = 0;
				for (int j = 0; j < idenType.length; j++) {
					if (idenType[j] == MetaPeptide.IDEN_MSMS) {
						idenCount++;
					}
				}

				if (totalIntensity == 0 || idenCount == 0) {
					continue;
				}

				double[] indiPEP = quanPeps[i].getIndividualPEP();
				for (int j = 0; j < pepIntensity.length; j++) {
					if (idenType[j] == MetaPeptide.IDEN_MSMS) {
						if (indiPEP[j] == 0) {
							peptideIntensity[i] += pepIntensity[j] * 1.0E7;
						} else {
							peptideIntensity[i] += pepIntensity[j] / indiPEP[j];
						}
					}
				}

				String[] proteins = quanPeps[i].getProteins();
				for (int j = 0; j < proteins.length; j++) {
//					logger.info(proteins[j]+"\t"+uniProInfoMap.containsKey(proteins[j])+"\t"+uniProInfoMap.size());
					if (uniProInfoMap.containsKey(proteins[j])) {
						String species = String.valueOf(uniProInfoMap.get(proteins[j]).getTaxId());
						if (speciesIntensityMap.containsKey(species)) {
							speciesIntensityMap.put(species, speciesIntensityMap.get(species) + peptideIntensity[i]);
						} else {
							speciesIntensityMap.put(species, peptideIntensity[i]);
						}
					} else if (proteins[j].startsWith(decoySymble)) {
						if (speciesIntensityMap.containsKey(proteins[j])) {
							speciesIntensityMap.put(proteins[j],
									speciesIntensityMap.get(proteins[j]) + peptideIntensity[i]);
						} else {
							speciesIntensityMap.put(proteins[j], peptideIntensity[i]);
						}
					}
				}
			}

			logger.info(getTaskName() + ": species before filtering: " + speciesIntensityMap.size());

			HashMap<String, Double> razorSpeciesIntenMap = new HashMap<String, Double>();
			for (int i = 0; i < quanPeps.length; i++) {
				if (peptideIntensity[i] > 0) {
					String[] proteins = quanPeps[i].getProteins();
					String razorSpecies = "";
					double razorSpeciesIntensity = 0;
					boolean use = false;
					for (int j = 0; j < proteins.length; j++) {
						if (uniProInfoMap.containsKey(proteins[j])) {
							String species = String.valueOf(uniProInfoMap.get(proteins[j]).getTaxId());
							if (speciesIntensityMap.containsKey(species)) {
								double speciesIntensity = speciesIntensityMap.get(species);
								if (speciesIntensity > razorSpeciesIntensity) {
									razorSpecies = species;
									razorSpeciesIntensity = speciesIntensity;
								}
								use = true;
							}
						} else if (proteins[j].startsWith(decoySymble)) {
							if (speciesIntensityMap.containsKey(proteins[j])) {
								double speciesIntensity = speciesIntensityMap.get(proteins[j]);
								if (speciesIntensity > razorSpeciesIntensity) {
									razorSpecies = proteins[j];
									razorSpeciesIntensity = speciesIntensity;
								}
							}
							use = true;
						}
					}

					if (use) {
						razorSpeciesIntenMap.put(razorSpecies, 0.0);
					}
				}
			}

			logger.info(getTaskName() + ": razor species before filtering: " + razorSpeciesIntenMap.size());

			// Evaluation-only diagnostics: verify how many intensity-eligible peptides can
			// be linked to at least one razor species candidate.
			// This block does not affect filtering/scoring logic.
			int eligiblePeptideCount = 0;
			int allCoveredPeptideCount = 0;
			int razorCoveredPeptideCount = 0;
			ArrayList<String> uncoveredPeptideExamples = new ArrayList<String>();
			for (int i = 0; i < quanPeps.length; i++) {
				if (peptideIntensity[i] <= 0) {
					continue;
				}
				eligiblePeptideCount++;
				String[] proteins = quanPeps[i].getProteins();
				boolean allCovered = false;
				boolean razorCovered = false;
				for (int j = 0; j < proteins.length; j++) {
					if (uniProInfoMap.containsKey(proteins[j])) {
						String species = String.valueOf(uniProInfoMap.get(proteins[j]).getTaxId());
						if (razorSpeciesIntenMap.containsKey(species)) {
							razorCovered = true;
						}
						if (speciesIntensityMap.containsKey(species)) {
							allCovered = true;
						}
					} 
				}
				if (razorCovered) {
					razorCoveredPeptideCount++;
				} 
				if (allCovered) {
					allCoveredPeptideCount++;
				} 
			}
			double razorCoveragePercent = eligiblePeptideCount == 0 ? 0.0
					: (double) razorCoveredPeptideCount * 100.0 / (double) eligiblePeptideCount;
			double allCoveragePercent = eligiblePeptideCount == 0 ? 0.0
					: (double) allCoveredPeptideCount * 100.0 / (double) eligiblePeptideCount;
			logger.info(getTaskName() + ": razor species peptide coverage check: " + razorCoveredPeptideCount + "/"
					+ eligiblePeptideCount + " peptides (" + FormatTool.getDF2().format(razorCoveragePercent)
					+ "%) are linked to razor species candidates");
			logger.info(getTaskName() + ": razor species peptide coverage check: " + allCoveredPeptideCount + "/"
					+ eligiblePeptideCount + " peptides (" + FormatTool.getDF2().format(allCoveragePercent)
					+ "%) are linked to razor species candidates");
			if (razorCoveredPeptideCount < eligiblePeptideCount) {
				logger.warn(getTaskName() + ": razor species uncovered peptide examples (up to 10): "
						+ uncoveredPeptideExamples);
			}

			for (int i = 0; i < quanPeps.length; i++) {
				if (peptideIntensity[i] > 0) {
					String[] proteins = quanPeps[i].getProteins();

					int count = 0;
					for (int j = 0; j < proteins.length; j++) {
						if (uniProInfoMap.containsKey(proteins[j])) {
							String species = String.valueOf(uniProInfoMap.get(proteins[j]).getTaxId());
							if (razorSpeciesIntenMap.containsKey(species)) {
								count++;
							}
						} else if (proteins[j].startsWith(decoySymble)) {
							if (razorSpeciesIntenMap.containsKey(proteins[j])) {
								count++;
							}
						}
					}

					if (count > 0) {
						for (int j = 0; j < proteins.length; j++) {
							if (uniProInfoMap.containsKey(proteins[j])) {
								String species = String.valueOf(uniProInfoMap.get(proteins[j]).getTaxId());
								if (razorSpeciesIntenMap.containsKey(species)) {
									razorSpeciesIntenMap.put(species,
											razorSpeciesIntenMap.get(species) + peptideIntensity[i]);
								}
							} else if (proteins[j].startsWith(decoySymble)) {
								if (razorSpeciesIntenMap.containsKey(proteins[j])) {
									razorSpeciesIntenMap.put(proteins[j],
											razorSpeciesIntenMap.get(proteins[j]) + peptideIntensity[i]);
								}
							}
						}
					}
				}
			}

			List<Map.Entry<String, Double>> sortedRazorSpeciesList = razorSpeciesIntenMap.entrySet().stream()
					.sorted(Map.Entry.comparingByValue(Collections.reverseOrder())).collect(Collectors.toList());
			String[] razorSpecies = sortedRazorSpeciesList.stream().map(Map.Entry::getKey).toArray(String[]::new);

			ArrayList<Double> revLogIntenList = new ArrayList<Double>();
			for (int i = 0; i < razorSpecies.length; i++) {
				if (razorSpecies[i].startsWith(decoySymble)) {
					double intensity = razorSpeciesIntenMap.get(razorSpecies[i]);
					if (intensity > 0) {
						revLogIntenList.add(Math.log10(intensity));
					}
				}
			}

			logger.info(getTaskName() + ": decoy species count: " + revLogIntenList.size());

			HashMap<String, Double> genomeScoreMap = new HashMap<String, Double>();

			if (revLogIntenList.size() == 0) {
				logger.info(getTaskName()
						+ ": no decoy genome were found, the p-Value for the genome identifications will be set as 0");

				for (int i = 0; i < razorSpecies.length; i++) {
					double intensity = razorSpeciesIntenMap.get(razorSpecies[i]);
					if (intensity == 0) {
						break;
					}
					genomeScoreMap.put(razorSpecies[i], intensity);
				}
			} else {
				TTest tTest = new TTest();
				double[] revLogIntensities = new double[revLogIntenList.size()];
				for (int i = 2; i < revLogIntenList.size(); i++) {
					revLogIntensities = new double[i];
					for (int j = 0; j < revLogIntensities.length; j++) {
						revLogIntensities[j] = revLogIntenList.get(j);
					}

					double pValue = tTest.tTest(revLogIntensities[0], revLogIntensities);
					if (pValue < 0.01) {
						break;
					}
				}
				logger.info(getTaskName() + ": refined decoy genome:" + revLogIntensities.length);
				logger.info(Arrays.toString(revLogIntensities));

				if (revLogIntensities.length == 0) {
					logger.info(getTaskName() + ": the number of decoy genome was not sufficient for the t-test, "
							+ "the p-Value for the genome identifications will be set as 0");
					for (int i = 0; i < razorSpecies.length; i++) {
						double intensity = razorSpeciesIntenMap.get(razorSpecies[i]);
						if (razorSpecies[i].startsWith(decoySymble) || intensity == 0) {
							break;
						}
						genomeScoreMap.put(razorSpecies[i], 0.0);
					}
				} else if (revLogIntensities.length < 6) {
					double[] extended = new double[6];
					for (int i = 0; i < extended.length; i++) {
						extended[i] = revLogIntensities[i % revLogIntensities.length];
					}

					logger.info(getTaskName() + ": decoy genome used for scoring: " + extended.length);

					for (int i = 0; i < razorSpecies.length; i++) {
						double intensity = razorSpeciesIntenMap.get(razorSpecies[i]);

						if (razorSpecies[i].startsWith(decoySymble) || intensity == 0) {
							break;
						}

						double pvalue = tTest.tTest(Math.log10(intensity), extended);
						logger.info(getTaskName() + ": pvalue for genome " + razorSpecies[i] + ": " + pvalue);

						if (pvalue < 0.01) {
							genomeScoreMap.put(razorSpecies[i], pvalue);
						} else {
//							break;
						}
					}
				} else {
					logger.info(getTaskName() + ": decoy genome used for scoring: " + revLogIntensities.length);

					for (int i = 0; i < razorSpecies.length; i++) {
						double intensity = razorSpeciesIntenMap.get(razorSpecies[i]);

						if (razorSpecies[i].startsWith(decoySymble) || intensity == 0) {
							break;
						}

						double pvalue = tTest.tTest(Math.log10(intensity), revLogIntensities);
						logger.info(getTaskName() + ": pvalue for genome " + razorSpecies[i] + ": " + pvalue);

						if (pvalue < 0.01) {
							genomeScoreMap.put(razorSpecies[i], pvalue);
						} else {
//							break;
						}
					}
				}
			}

			logger.info(getTaskName() + ": species after filtering: " + genomeScoreMap.size());

			return genomeScoreMap;
		}

		/**
		 * Refines proteins using razor protein selection with default decoy symbol.
		 * Selects the best protein for each peptide based on genome scores and
		 * intensity.
		 * 
		 * USED BY: - DDAMagPfindTask: Line 1332 (razor protein selection for protein
		 * inference) - DDAMagSageTask: Line 1238 (razor protein selection for protein
		 * inference)
		 * 
		 * @param quanPeps       Array of quantified peptides
		 * @param genomeScoreMap Map of genome names to their FDR scores
		 * @return Array of selected razor proteins (one per peptide)
		 */
		protected String[] refineProteins(MetaPeptide[] quanPeps, HashMap<String, Double> genomeScoreMap,
				HashMap<String, UnipProInfo> uniProInfoMap) {
			return refineProteins(quanPeps, genomeScoreMap, uniProInfoMap, PFindResultReader.REV);
		}

		/**
		 * Refines proteins using razor protein selection with custom decoy symbol.
		 * Selects the best protein for each peptide based on genome scores and
		 * intensity.
		 * 
		 * USED BY: - DDAMagFragPipeTask: Line 1184 (razor protein selection with
		 * FastaManager.shuffle decoy symbol)
		 * 
		 * @param quanPeps       Array of quantified peptides
		 * @param genomeScoreMap Map of genome names to their FDR scores
		 * @param decoySymble    Custom decoy symbol for protein filtering
		 * @return Array of selected razor proteins (one per peptide)
		 */
		protected String[] refineProteins(MetaPeptide[] quanPeps, HashMap<String, Double> speciesScoreMap,
				HashMap<String, UnipProInfo> uniProInfoMap, String decoySymble) {

			HashMap<String, Double> proteinIntensityMap = new HashMap<String, Double>();
			for (int i = 0; i < quanPeps.length; i++) {
				double[] pepIntensity = quanPeps[i].getIntensity();
				String[] proteins = quanPeps[i].getProteins();

				int count = 0;
				for (int j = 0; j < proteins.length; j++) {
					if (!proteins[j].startsWith(decoySymble)) {
						if (uniProInfoMap.containsKey(proteins[j])) {
							String species = String.valueOf(uniProInfoMap.get(proteins[j]).getTaxId());
							if (speciesScoreMap.containsKey(species)) {
								count++;
							}
						} else {
							count++;
						}
					}
				}

				if (count > 0) {
					double intensity = MathTool.getTotal(pepIntensity) / (double) count;
					for (int j = 0; j < proteins.length; j++) {
						if (!proteins[j].startsWith(decoySymble)) {
							if (uniProInfoMap.containsKey(proteins[j])) {
								String species = String.valueOf(uniProInfoMap.get(proteins[j]).getTaxId());
								if (speciesScoreMap.containsKey(species)) {
									if (proteinIntensityMap.containsKey(proteins[j])) {
										proteinIntensityMap.put(proteins[j],
												proteinIntensityMap.get(proteins[j]) + intensity);
									} else {
										proteinIntensityMap.put(proteins[j], intensity);
									}
								}
							} else {
								if (proteinIntensityMap.containsKey(proteins[j])) {
									proteinIntensityMap.put(proteins[j],
											proteinIntensityMap.get(proteins[j]) + intensity);
								} else {
									proteinIntensityMap.put(proteins[j], intensity);
								}
							}
						}
					}
				}
			}

			String[] razorProteins = new String[quanPeps.length];

			for (int i = 0; i < quanPeps.length; i++) {
				String[] proteins = quanPeps[i].getProteins();
				double genomeScore = Double.MAX_VALUE;
				double proteinIntensity = 0;

				for (int j = 0; j < proteins.length; j++) {

					if (!proteins[j].startsWith(decoySymble)) {
						if (uniProInfoMap.containsKey(proteins[j])) {
							String species = String.valueOf(uniProInfoMap.get(proteins[j]).getTaxId());
							if (speciesScoreMap.containsKey(species)) {
								if (speciesScoreMap.get(species) < genomeScore) {
									genomeScore = speciesScoreMap.get(species);
									razorProteins[i] = proteins[j];
									proteinIntensity = proteinIntensityMap.get(proteins[j]);
								} else if (speciesScoreMap.get(species) == genomeScore) {
									if (proteinIntensityMap.get(proteins[j]) > proteinIntensity) {
										proteinIntensity = proteinIntensityMap.get(proteins[j]);
										razorProteins[i] = proteins[j];
									}
								}
							}
						} else {
							if (proteinIntensityMap.get(proteins[j]) > proteinIntensity) {
								proteinIntensity = proteinIntensityMap.get(proteins[j]);
								razorProteins[i] = proteins[j];
							}
						}
					}
				}
			}

			return razorProteins;
		}

	}
}
