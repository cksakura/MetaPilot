package uk.ac.quadram.metapilot.task.mag;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutionException;
import java.util.stream.Collectors;

import org.apache.commons.math3.stat.inference.TTest;

import uk.ac.quadram.metapilot.core.function.v2.EnzymeCommission;
import uk.ac.quadram.metapilot.core.function.v2.GoObo;
import uk.ac.quadram.metapilot.core.math.MathTool;
import uk.ac.quadram.metapilot.core.prodb.FastaManager;
import uk.ac.quadram.metapilot.core.taxonomy.TaxonomyRanks;
import uk.ac.quadram.metapilot.core.tools.FormatTool;
import uk.ac.quadram.metapilot.dbSearch.pfind.PFindResultReader;
import javafx.scene.control.ProgressBar;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters.ToolType;
import uk.ac.quadram.metapilot.task.MetaPilotTask;
import uk.ac.quadram.metapilot.task.MetaPythonReportCopyTask;
import uk.ac.quadram.metapilot.task.MetaPythonReportTask;
import uk.ac.quadram.metapilot.task.io.MetaTreeHandler;
import uk.ac.quadram.metapilot.task.io.pep.MetaPeptide;
import uk.ac.quadram.metapilot.task.io.pro.MetaProtein;
import uk.ac.quadram.metapilot.task.io.pro.MetaProteinXMLWriter2;
import uk.ac.quadram.metapilot.task.par.MetaBaseParameter;
import uk.ac.quadram.metapilot.task.par.MetaConstants;

public abstract class MetaPilotMagTask extends MetaPilotTask {

	protected String[] rawFiles;
	protected String[] quanExpNames;
	protected File spectraDirFile;
	protected MagDbItem magDb;
	protected File[] separateResultFolders;
	protected int threadCount;
	protected double psmCoverThreshold = 0.95;

	protected File rawDir;
	protected File quan_pep_file;
	protected File quan_pro_file;
	protected HashMap<String, String> expNameMap;
	protected File summaryMetaFile;

	public MetaPilotMagTask(MetaBaseParameter metaPar, MetaResourceParameters mrp, ProgressBar progressBar,
			ProgressBar secondBar) {
		super(metaPar, mrp, progressBar, secondBar);

		// Report templates will be copied in background during call() method
		// This prevents UI freeze during task construction
		LOGGER.info(getTaskName() + ": Constructor completed, report templates will be copied during execution");
	}

	protected HashSet<String> refineGenomes(HashMap<String, HashSet<String>> genomePepMap, int totalPepCount,
			double threshold, File out) {
		return refineGenomes(genomePepMap, totalPepCount, threshold, 0.001, 250, false, out);
	}

	protected HashSet<String> refineGenomes(HashMap<String, HashSet<String>> genomePepMap, int totalPepCount,
			double threshold, double stepThreshold, int maxium, boolean decoyFilter, File out) {

		String[] genomes = genomePepMap.keySet().toArray(new String[genomePepMap.size()]);
		Arrays.sort(genomes, (g1, g2) -> {
			int s1 = genomePepMap.get(g1).size();
			int s2 = genomePepMap.get(g2).size();
			return s2 - s1;
		});

		for (int iteratorCount = 0; iteratorCount < 5; iteratorCount++) {
			HashSet<String> totalSet = new HashSet<String>();
			HashMap<String, Integer> orderMap = new HashMap<String, Integer>();
			for (int i = 0; i < genomes.length; i++) {
				int currentSize = totalSet.size();
				totalSet.addAll(genomePepMap.get(genomes[i]));
				int increaseSize = totalSet.size() - currentSize;
				orderMap.put(genomes[i], increaseSize);
			}

			Arrays.sort(genomes, (g1, g2) -> {
				return orderMap.get(g2) - orderMap.get(g1);
			});
		}

		double[] percentage = new double[genomes.length];
		HashSet<String> genomeSet = new HashSet<String>();
		HashSet<String> totalSet = new HashSet<String>();
		int breakPoint1 = -1;
		int breakPoint2 = -1;

		try (PrintWriter writer = new PrintWriter(out)) {
			writer.println(
					"Rank\tGenome\tPeptide_count\tTotal_covered_peptide_count\tPercentage\tIncrement\tIncrement ratio");
			for (int i = 0; i < genomes.length; i++) {

				int currentSize = totalSet.size();
				totalSet.addAll(genomePepMap.get(genomes[i]));
				int add = totalSet.size() - currentSize;
				double additional = ((double) add / (double) totalPepCount);

				percentage[i] = ((double) totalSet.size() / (double) totalPepCount);

				writer.println(i + "\t" + genomes[i] + "\t" + genomePepMap.get(genomes[i]).size() + "\t"
						+ totalSet.size() + "\t" + percentage[i] + "\t" + add + "\t" + additional);

				if (add > 1) {
					breakPoint1 = i;
				}
				if (additional > stepThreshold) {
					breakPoint2 = i;
				}
			}
			writer.close();
		} catch (IOException e) {
			LOGGER.error(getTaskName() + ": error in writing genomes to " + out, e);
			System.out
					.println(format.format(new Date()) + "\t" + getTaskName() + ": error in writing genomes to " + out);
		}

		int decoyCount = 0;
		for (int i = 0; i < genomes.length; i++) {
			genomeSet.add(genomes[i]);
			if (genomes[i].startsWith(FastaManager.shuffle)) {
				decoyCount++;
			}
			if (decoyFilter) {
				if (decoyCount >= 3) {
					if (percentage[i] >= threshold) {
						break;
					}
					if (i == breakPoint1 || i == breakPoint2) {
						break;
					}
					if (genomeSet.size() > maxium) {
						break;
					}
				}
			} else {
				if (percentage[i] >= threshold) {
					break;
				}
				if (threshold < 1.0) {
					if (i == breakPoint1 || i == breakPoint2) {
						break;
					}
					if (genomeSet.size() > maxium) {
						break;
					}
				}
			}
		}

		return genomeSet;
	}

	/**
	 * Reverser proteins are also included to calculate the FDR of genomes
	 * 
	 * @param quanPeps
	 * @return
	 */
	protected HashMap<String, Double> refineGenomes(MetaPeptide[] quanPeps) {
		return refineGenomes(quanPeps, PFindResultReader.REV);
	}

	/**
	 * Reverser proteins are also included to calculate the FDR of genomes
	 * 
	 * @param quanPeps
	 * @return
	 */
	protected HashMap<String, Double> refineGenomes(MetaPeptide[] quanPeps, String decoySymble) {

		String revIdentifier = decoySymble + this.magDb.getIdentifier();
		HashMap<String, Double> genomeIntensityMap = new HashMap<String, Double>();
		double[] peptideIntensity = new double[quanPeps.length];

		for (int i = 0; i < quanPeps.length; i++) {
			double[] pepIntensity = quanPeps[i].getIntensity();
			double totalIntensity = MathTool.getTotal(pepIntensity);
			int[] idenType = quanPeps[i].getIdenType();
			int idenCount = 0;
			for (int j = 0; j < idenType.length; j++) {
				if (idenType[j] == 0) {
					idenCount++;
				}
			}

			if (totalIntensity == 0 || idenCount == 0) {
				continue;
			}

			double[] indiPEP = quanPeps[i].getIndividualPEP();
			for (int j = 0; j < pepIntensity.length; j++) {
				if (idenType[j] == 0) {
					if (indiPEP[j] == 0) {
						peptideIntensity[i] += pepIntensity[j] * 1.0E7;
					} else {
						peptideIntensity[i] += pepIntensity[j] / indiPEP[j];
					}
				}
			}

			String[] proteins = quanPeps[i].getProteins();
			for (int j = 0; j < proteins.length; j++) {
				if (proteins[j].startsWith(this.magDb.getIdentifier())) {
					String genome = proteins[j].split("_")[0];
					if (genomeIntensityMap.containsKey(genome)) {
						genomeIntensityMap.put(genome, genomeIntensityMap.get(genome) + peptideIntensity[i]);
					} else {
						genomeIntensityMap.put(genome, peptideIntensity[i]);
					}
				} else if (proteins[j].startsWith(revIdentifier)) {
					String pro = this.magDb.getIdentifier() + proteins[j].substring(revIdentifier.length());
					String genome = decoySymble + pro.split("_")[0];

					if (genomeIntensityMap.containsKey(genome)) {
						genomeIntensityMap.put(genome, genomeIntensityMap.get(genome) + peptideIntensity[i]);
					} else {
						genomeIntensityMap.put(genome, peptideIntensity[i]);
					}
				}
			}
		}

		LOGGER.info(getTaskName() + ": genome before filtering: " + genomeIntensityMap.size());

		HashMap<String, Double> razorGenomeIntenMap = new HashMap<String, Double>();
		for (int i = 0; i < quanPeps.length; i++) {
			if (peptideIntensity[i] > 0) {
				String[] proteins = quanPeps[i].getProteins();
				String razorGenome = "";
				double razorGenomeIntensity = 0;
				for (int j = 0; j < proteins.length; j++) {
					if (proteins[j].startsWith(this.magDb.getIdentifier())) {
						String genome = proteins[j].split("_")[0];
						if (genomeIntensityMap.containsKey(genome)) {
							double genomeIntensity = genomeIntensityMap.get(genome);
							if (genomeIntensity > razorGenomeIntensity) {
								razorGenome = genome;
								razorGenomeIntensity = genomeIntensity;
							}
						}
					} else if (proteins[j].startsWith(revIdentifier)) {
						String pro = this.magDb.getIdentifier() + proteins[j].substring(revIdentifier.length());
						String genome = decoySymble + pro.split("_")[0];

						if (genomeIntensityMap.containsKey(genome)) {
							double genomeIntensity = genomeIntensityMap.get(genome);
							if (genomeIntensity > razorGenomeIntensity) {
								razorGenome = genome;
								razorGenomeIntensity = genomeIntensity;
							}
						}
					}
				}

				if (razorGenome.startsWith(this.magDb.getIdentifier()) || razorGenome.startsWith(revIdentifier)) {
					razorGenomeIntenMap.put(razorGenome, 0.0);
				}
			}
		}

		LOGGER.info(getTaskName() + ": razor genome before filtering: " + razorGenomeIntenMap.size());

		for (int i = 0; i < quanPeps.length; i++) {
			if (peptideIntensity[i] > 0) {
				String[] proteins = quanPeps[i].getProteins();

				int count = 0;
				for (int j = 0; j < proteins.length; j++) {
					if (proteins[j].startsWith(this.magDb.getIdentifier())) {
						String genome = proteins[j].split("_")[0];
						if (razorGenomeIntenMap.containsKey(genome)) {
							count++;
						}
					} else if (proteins[j].startsWith(revIdentifier)) {
						String pro = this.magDb.getIdentifier() + proteins[j].substring(revIdentifier.length());
						String genome = decoySymble + pro.split("_")[0];
						if (razorGenomeIntenMap.containsKey(genome)) {
							count++;
						}
					}
				}

				if (count > 0) {
					for (int j = 0; j < proteins.length; j++) {
						if (proteins[j].startsWith(this.magDb.getIdentifier())) {
							String genome = proteins[j].split("_")[0];
							if (razorGenomeIntenMap.containsKey(genome)) {
								razorGenomeIntenMap.put(genome, razorGenomeIntenMap.get(genome) + peptideIntensity[i]);
							}
						} else if (proteins[j].startsWith(revIdentifier)) {
							String pro = this.magDb.getIdentifier() + proteins[j].substring(revIdentifier.length());
							String genome = decoySymble + pro.split("_")[0];
							if (razorGenomeIntenMap.containsKey(genome)) {
								razorGenomeIntenMap.put(genome, razorGenomeIntenMap.get(genome) + peptideIntensity[i]);
							}
						}
					}
				}
			}
		}

		List<Map.Entry<String, Double>> sortedRazorGenomeList = razorGenomeIntenMap.entrySet().stream()
				.sorted(Map.Entry.comparingByValue(Collections.reverseOrder())).collect(Collectors.toList());
		String[] razorGenomes = sortedRazorGenomeList.stream().map(Map.Entry::getKey).toArray(String[]::new);

		ArrayList<Double> revLogIntenList = new ArrayList<Double>();
		for (int i = 0; i < razorGenomes.length; i++) {
			if (razorGenomes[i].startsWith(revIdentifier)) {
				double intensity = razorGenomeIntenMap.get(razorGenomes[i]);
				if (intensity > 0) {
					revLogIntenList.add(Math.log10(intensity));
				}
			}
		}

		LOGGER.info(getTaskName() + ": decoy genome count: " + revLogIntenList.size());

		HashMap<String, Double> genomeScoreMap = new HashMap<String, Double>();

		if (revLogIntenList.size() == 0) {
			LOGGER.info(getTaskName()
					+ ": no decoy genome were found, the p-Value for the genome identifications will be set as 0");

			for (int i = 0; i < razorGenomes.length; i++) {
				double intensity = razorGenomeIntenMap.get(razorGenomes[i]);
				if (intensity == 0) {
					break;
				}
				genomeScoreMap.put(razorGenomes[i], intensity);
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
			LOGGER.info(getTaskName() + ": refined decoy genome:" + revLogIntensities.length);
			LOGGER.info(Arrays.toString(revLogIntensities));

			if (revLogIntensities.length == 0) {
				LOGGER.info(getTaskName() + ": the number of decoy genome was not sufficient for the t-test, "
						+ "the p-Value for the genome identifications will be set as 0");
				for (int i = 0; i < razorGenomes.length; i++) {
					double intensity = razorGenomeIntenMap.get(razorGenomes[i]);
					if (razorGenomes[i].startsWith(revIdentifier) || intensity == 0) {
						break;
					}
					genomeScoreMap.put(razorGenomes[i], 0.0);
				}
			} else if (revLogIntensities.length < 6) {
				double[] extended = new double[6];
				for (int i = 0; i < extended.length; i++) {
					extended[i] = revLogIntensities[i % revLogIntensities.length];
				}

				LOGGER.info(getTaskName() + ": decoy genome used for scoring: " + extended.length);

				for (int i = 0; i < razorGenomes.length; i++) {
					double intensity = razorGenomeIntenMap.get(razorGenomes[i]);

					if (razorGenomes[i].startsWith(revIdentifier) || intensity == 0) {
						break;
					}

					double pvalue = tTest.tTest(Math.log10(intensity), extended);
					LOGGER.info(getTaskName() + ": pvalue for genome " + razorGenomes[i] + ": " + pvalue);

					if (pvalue < 0.01) {
						genomeScoreMap.put(razorGenomes[i], pvalue);
					} else {
//						break;
					}
				}
			} else {
				LOGGER.info(getTaskName() + ": decoy genome used for scoring: " + revLogIntensities.length);

				for (int i = 0; i < razorGenomes.length; i++) {
					double intensity = razorGenomeIntenMap.get(razorGenomes[i]);

					if (razorGenomes[i].startsWith(revIdentifier) || intensity == 0) {
						break;
					}

					double pvalue = tTest.tTest(Math.log10(intensity), revLogIntensities);
					LOGGER.info(getTaskName() + ": pvalue for genome " + razorGenomes[i] + ": " + pvalue);

					if (pvalue < 0.01) {
						genomeScoreMap.put(razorGenomes[i], pvalue);
					} else {
//						break;
					}
				}
			}
		}

		LOGGER.info(getTaskName() + ": genome after filtering: " + genomeScoreMap.size());

		return genomeScoreMap;
	}

	protected String[] refineProteins(MetaPeptide[] quanPeps, HashMap<String, Double> genomeScoreMap) {
		return refineProteins(quanPeps, genomeScoreMap, PFindResultReader.REV);
	}

	protected String[] refineProteins(MetaPeptide[] quanPeps, HashMap<String, Double> genomeScoreMap,
			String decoySymble) {

		HashMap<String, Double> proteinIntensityMap = new HashMap<String, Double>();
		for (int i = 0; i < quanPeps.length; i++) {
			double[] pepIntensity = quanPeps[i].getIntensity();
			String[] proteins = quanPeps[i].getProteins();

			int count = 0;
			for (int j = 0; j < proteins.length; j++) {
				if (!proteins[j].startsWith(decoySymble)) {
					if (proteins[j].startsWith(this.magDb.getIdentifier())) {
						String genome = proteins[j].split("_")[0];
						if (genomeScoreMap.containsKey(genome)) {
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
						if (proteins[j].startsWith(this.magDb.getIdentifier())) {
							String genome = proteins[j].split("_")[0];
							if (genomeScoreMap.containsKey(genome)) {
								if (proteinIntensityMap.containsKey(proteins[j])) {
									proteinIntensityMap.put(proteins[j],
											proteinIntensityMap.get(proteins[j]) + intensity);
								} else {
									proteinIntensityMap.put(proteins[j], intensity);
								}
							}
						} else {
							if (proteinIntensityMap.containsKey(proteins[j])) {
								proteinIntensityMap.put(proteins[j], proteinIntensityMap.get(proteins[j]) + intensity);
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
					if (proteins[j].startsWith(this.magDb.getIdentifier())) {
						String genome = proteins[j].split("_")[0];
						if (genomeScoreMap.containsKey(genome)) {
							if (genomeScoreMap.get(genome) < genomeScore) {
								genomeScore = genomeScoreMap.get(genome);
								razorProteins[i] = proteins[j];
								proteinIntensity = proteinIntensityMap.get(proteins[j]);
							} else if (genomeScoreMap.get(genome) == genomeScore) {
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

	protected boolean exportGenome(HashMap<String, double[]> genomeIntensityMap, HashMap<String, Integer> psmCountMap,
			HashMap<String, HashSet<String>> pepCountMap, HashMap<String, HashSet<String>> uniquePepCountMap,
			HashMap<String, ArrayList<Double>> genomeQvalueMap, File taxFile, String[] expNames,
			MetaPythonReportTask task) {

		MagFuncSearcher searcher = null;
		try {
			searcher = new MagFuncSearcher(metaPar.getMicroDb(), mrp.getToolPath(ToolType.FUNCDB));

		} catch (NumberFormatException | SQLException | IOException e) {
			// TODO Auto-generated catch block
			LOGGER.error(getTaskName() + ": error in genome analysis, task failed", e);
			System.err.println(
					format.format(new Date()) + "\t" + getTaskName() + ": error in genome analysis, task failed");

			return false;
		}

		String[] genomes = genomeIntensityMap.keySet().toArray(new String[genomeIntensityMap.size()]);
		HashMap<String, String[]> genomeTaxaMap = null;
		try {
			genomeTaxaMap = searcher.matchGenome(genomes);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			LOGGER.error(getTaskName() + ": Error in searching the genome2taxa database", e);
			System.err.println(
					format.format(new Date()) + "\t" + getTaskName() + ": error in searching the genome2taxa database");
			return false;
		}

		File genomeHtmlFile = new File(taxFile, "Genome.html");
		try {
			HMGenomeHtmlWriter htmlWriter = new HMGenomeHtmlWriter(genomeHtmlFile);
			htmlWriter.write(genomeTaxaMap);
			htmlWriter.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			LOGGER.error("Error in writing genome information to " + genomeHtmlFile, e);
			System.err.println(
					format.format(new Date()) + "\t" + "Error in writing genome information to " + genomeHtmlFile);
			return false;
		}

		File genomeFile = new File(taxFile, "Genome.tsv");
		PrintWriter genomeWriter = null;
		try {
			genomeWriter = new PrintWriter(genomeFile);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			LOGGER.error("Error in writing genome information to " + genomeFile, e);
			System.err
					.println(format.format(new Date()) + "\t" + "Error in writing genome information to " + genomeFile);
			return false;
		}

		StringBuilder genomeTitlesb = new StringBuilder();
		genomeTitlesb.append("Genome\t");
		genomeTitlesb.append("Score\t");
		genomeTitlesb.append("PSM count\t");
		genomeTitlesb.append("Peptide count\t");
		genomeTitlesb.append("Unique peptide count\t");
		genomeTitlesb.append("GTDB_Name\t");
		genomeTitlesb.append("GTDB_Rank\t");
		TaxonomyRanks[] mainRanks = TaxonomyRanks.getMainRanks();
		for (int i = 0; i < mainRanks.length; i++) {
			genomeTitlesb.append("GTDB_").append(mainRanks[i].getName()).append("\t");
		}

		genomeTitlesb.append("NCBI_Name\t");
		genomeTitlesb.append("NCBI_Rank");

		for (int i = 0; i < mainRanks.length; i++) {
			genomeTitlesb.append("\t").append("NCBI_").append(mainRanks[i].getName());
		}

		for (String exp : expNames) {
			if (exp.startsWith("Intensity ")) {
				genomeTitlesb.append("\t").append(exp);
			} else {
				genomeTitlesb.append("\tIntensity ").append(exp);
			}
		}

		genomeWriter.println(genomeTitlesb);

		updateProgress(94, 100);

		double[] cellularIntensity = new double[expNames.length];
		HashMap<String, String[]>[] taxRankMaps = new HashMap[mainRanks.length];
		HashMap<String, double[]>[] taxIntensityMaps = new HashMap[mainRanks.length];
		HashMap<String, Double>[] taxScoreMaps = new HashMap[mainRanks.length];
		for (int i = 0; i < taxRankMaps.length; i++) {
			taxRankMaps[i] = new HashMap<String, String[]>();
			taxIntensityMaps[i] = new HashMap<String, double[]>();
			taxScoreMaps[i] = new HashMap<String, Double>();
		}

		Arrays.sort(genomes, new Comparator<String>() {
			@Override
			public int compare(String o1, String o2) {
				// TODO Auto-generated method stub
				return genomeQvalueMap.get(o2).size() - genomeQvalueMap.get(o1).size();
			}
		});

		double[] logGenomePepCount = new double[genomes.length];
		for (int i = 0; i < genomes.length; i++) {
			logGenomePepCount[i] = Math.log(genomeQvalueMap.get(genomes[i]).size());
		}
		double[] genomeScores = new double[genomes.length];
		double mean = MathTool.calculateArithmeticMean(logGenomePepCount);
		if (mean > 0) {
			TTest tTest = new TTest();
			double add = -Math.log10(tTest.tTest(logGenomePepCount[0], logGenomePepCount));
			genomeScores[0] = add * 2;
			for (int i = 1; i < logGenomePepCount.length; i++) {
				if (logGenomePepCount[i] > 0) {
					double pValue = tTest.tTest(logGenomePepCount[i], logGenomePepCount);
					genomeScores[i] = -Math.log10(pValue);
					if (logGenomePepCount[i] > mean) {
						genomeScores[i] = add + genomeScores[i];
					} else {
						genomeScores[i] = add - genomeScores[i];
					}

					if (genomeScores[i] < 0) {
						genomeScores[i] = 0;
					}
				}
			}
		}

		for (int i = 0; i < genomes.length; i++) {
			String[] taxaContent = genomeTaxaMap.get(genomes[i]);
			if (genomeIntensityMap.containsKey(genomes[i])) {
				double[] genomeIntensity = genomeIntensityMap.get(genomes[i]);
				for (int j = 0; j < genomeIntensity.length; j++) {
					cellularIntensity[j] += genomeIntensity[j];
				}

				StringBuilder genomeTaxaSb = new StringBuilder();
				genomeTaxaSb.append(genomes[i]).append("\t");

				genomeTaxaSb.append(FormatTool.getDF2().format(genomeScores[i])).append("\t");

				if (psmCountMap.containsKey(genomes[i])) {
					genomeTaxaSb.append(psmCountMap.get(genomes[i])).append("\t");
				} else {
					genomeTaxaSb.append("0").append("\t");
				}
				if (pepCountMap.containsKey(genomes[i])) {
					genomeTaxaSb.append(pepCountMap.get(genomes[i]).size()).append("\t");
				} else {
					genomeTaxaSb.append("0").append("\t");
				}
				if (uniquePepCountMap.containsKey(genomes[i])) {
					genomeTaxaSb.append(uniquePepCountMap.get(genomes[i]).size());
				} else {
					genomeTaxaSb.append("0");
				}

				for (int j = 1; j < taxaContent.length; j++) {
					genomeTaxaSb.append("\t").append(taxaContent[j]);
				}

				for (int j = 0; j < genomeIntensity.length; j++) {
					genomeTaxaSb.append("\t").append(genomeIntensity[j]);
				}
				genomeWriter.println(genomeTaxaSb);

				for (int j = 0; j < taxRankMaps.length; j++) {
					String taxon = taxaContent[j + 3];
					if (taxon == null || taxon.length() == 0) {
						continue;
					}
					if (!taxRankMaps[j].containsKey(taxon)) {
						String[] lineage = new String[j + 1];
						System.arraycopy(taxaContent, 3, lineage, 0, lineage.length);
						taxRankMaps[j].put(taxon, lineage);
					}

					double[] taxIntensity;
					if (taxIntensityMaps[j].containsKey(taxon)) {
						taxIntensity = taxIntensityMaps[j].get(taxon);
					} else {
						taxIntensity = new double[genomeIntensity.length];
					}

					for (int k = 0; k < taxIntensity.length; k++) {
						taxIntensity[k] += genomeIntensity[k];
					}
					taxIntensityMaps[j].put(taxon, taxIntensity);

					if (taxScoreMaps[j].containsKey(taxon)) {
						taxScoreMaps[j].put(taxon, taxScoreMaps[j].get(taxon) + genomeScores[i]);
					} else {
						taxScoreMaps[j].put(taxon, genomeScores[i]);
					}
				}
			}
		}
		genomeWriter.close();

		File taxaFile = new File(taxFile, "Taxa.tsv");
		File taxaReportFile = new File(taxFile, "Taxa_report.tsv");
		PrintWriter taxaWriter = null;
		PrintWriter taxaReportWriter = null;
		try {
			taxaWriter = new PrintWriter(taxaFile);
			taxaReportWriter = new PrintWriter(taxaReportFile);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			LOGGER.error("Error in writing taxa information to " + taxaFile, e);
			System.err.println(format.format(new Date()) + "\t" + "Error in writing taxa information to " + taxaFile);
			return false;
		}

		StringBuilder taxaTitlesb = new StringBuilder();
		taxaTitlesb.append("Name\t");
		taxaTitlesb.append("Rank");

		for (int i = 0; i < mainRanks.length; i++) {
			taxaTitlesb.append("\t").append(mainRanks[i].getName());
		}

		taxaTitlesb.append("\tScore");

		for (String exp : expNames) {
			if (exp.startsWith("Intensity ")) {
				taxaTitlesb.append("\t").append(exp);
			} else {
				taxaTitlesb.append("\tIntensity ").append(exp);
			}
		}
		taxaWriter.println(taxaTitlesb);
		taxaReportWriter.println(taxaTitlesb);

		double maxTaxonScore = 0;
		for (int i = 0; i < taxRankMaps.length; i++) {
			for (String taxon : taxRankMaps[i].keySet()) {
				String[] taxaContent = taxRankMaps[i].get(taxon);
				double[] intensity = taxIntensityMaps[i].get(taxon);

				StringBuilder taxaSb = new StringBuilder();
				taxaSb.append(taxon).append("\t");
				taxaSb.append(mainRanks[i].getName());

				for (int j = 0; j < mainRanks.length; j++) {
					if (j < taxaContent.length) {
						taxaSb.append("\t").append(taxaContent[j]);
					} else {
						taxaSb.append("\t");
					}
				}

				double taxonScore = taxScoreMaps[i].get(taxon);
				taxaSb.append("\t").append(FormatTool.getDF2().format(taxonScore));

				for (int j = 0; j < intensity.length; j++) {
					taxaSb.append("\t").append(intensity[j]);
				}

				taxaWriter.println(taxaSb);

				if (taxonScore > maxTaxonScore) {
					maxTaxonScore = taxonScore;
				}
			}
		}

		taxaWriter.close();

		for (int i = 0; i < taxRankMaps.length; i++) {
			for (String taxon : taxRankMaps[i].keySet()) {
				double taxonScore = taxScoreMaps[i].get(taxon);
				if (taxonScore < maxTaxonScore * 0.005) {
					continue;
				}
				String[] taxaContent = taxRankMaps[i].get(taxon);
				double[] intensity = taxIntensityMaps[i].get(taxon);

				StringBuilder taxaSb = new StringBuilder();
				taxaSb.append(taxon).append("\t");
				taxaSb.append(mainRanks[i].getName());

				for (int j = 0; j < mainRanks.length; j++) {
					if (j < taxaContent.length) {
						taxaSb.append("\t").append(taxaContent[j]);
					} else {
						taxaSb.append("\t");
					}
				}

				taxaSb.append("\t").append(FormatTool.getDF2().format(taxonScore));

				for (int j = 0; j < intensity.length; j++) {
					taxaSb.append("\t").append(intensity[j]);
				}

				taxaReportWriter.println(taxaSb);
			}
		}

		taxaReportWriter.close();

		File jsOutput = new File(metaPar.getReportDir() + "\\reports\\metamap\\data\\tree.js");
		if (!jsOutput.getParentFile().exists()) {
			jsOutput.getParentFile().mkdirs();
		}

		if (!jsOutput.exists() || jsOutput.length() == 0) {
			MetaTreeHandler.exportJS(taxRankMaps, taxIntensityMaps, cellularIntensity, expNames, jsOutput);
		}

		boolean finished = false;

		File report_taxonomy_summary = new File(metaPar.getReportDir(), "report_taxonomy_summary.html");
		if (!report_taxonomy_summary.exists() || report_taxonomy_summary.length() == 0) {

			if (summaryMetaFile != null && summaryMetaFile.exists()) {
				finished = task.addTask(MetaPythonReportTask.TAXONOMY, taxaReportFile, report_taxonomy_summary,
						summaryMetaFile);
			} else {
				finished = task.addTask(MetaPythonReportTask.TAXONOMY, taxaReportFile, report_taxonomy_summary);
			}
		} else {
			finished = true;
		}

		LOGGER.info(getTaskName() + ": taxonomy analysis started");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": taxonomy analysis started");

		updateProgress(96, 100);

		return finished;
	}

	protected MetaProteinAnnoMag[] funcAnnoSqlite(MetaProtein[] pros, String[] fileNames) {
		HashMap<String, ArrayList<MetaProtein>> genomeProMap = new HashMap<String, ArrayList<MetaProtein>>();
		for (int i = 0; i < pros.length; i++) {
			String name = pros[i].getName();
			if (name.startsWith(magDb.getIdentifier())) {
				String genome = name.split("_")[0];
				if (genome.endsWith(".1")) {
					genome = genome.substring(0, genome.length() - 2);
				}

				if (genomeProMap.containsKey(genome)) {
					genomeProMap.get(genome).add(pros[i]);
				} else {
					ArrayList<MetaProtein> list = new ArrayList<MetaProtein>();
					list.add(pros[i]);
					genomeProMap.put(genome, list);
				}
			}
		}

		HashMap<String, String[]> usedCogMap = new HashMap<String, String[]>();
		HashMap<String, String[]> usedNogMap = new HashMap<String, String[]>();
		HashMap<String, String> usedKeggMap = new HashMap<String, String>();
		HashMap<String, GoObo> usedGoMap = new HashMap<String, GoObo>();
		HashMap<String, EnzymeCommission> usedEcMap = new HashMap<String, EnzymeCommission>();

		ArrayList<MetaProteinAnnoMag> list = new ArrayList<MetaProteinAnnoMag>();
		for (String genome : genomeProMap.keySet()) {

			ArrayList<MetaProtein> mplist = genomeProMap.get(genome);
			MetaProtein[] mps = mplist.toArray(new MetaProtein[mplist.size()]);

			MagFuncSearcher searcher = null;
			try {
				searcher = new MagFuncSearcher(metaPar.getMicroDb(), mrp.getToolPath(ToolType.FUNCDB), genome);

			} catch (NumberFormatException | SQLException | IOException e) {
				// TODO Auto-generated catch block
				LOGGER.error(getTaskName() + ": error in proteins functional annotations, task failed", e);
				System.err.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": error in proteins functional annotations, task failed");

				return null;
			}

			MetaProteinAnnoMag[] mpas = null;
			try {
				mpas = searcher.match(mps);
			} catch (IOException | SQLException e) {
				// TODO Auto-generated catch block
				LOGGER.error(getTaskName() + ": error in proteins functional annotations, task failed", e);
				System.err.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": error in proteins functional annotations, task failed");

				return null;
			}

			if (mpas != null) {
				for (int i = 0; i < mpas.length; i++) {
					list.add(mpas[i]);
				}

				usedCogMap.putAll(searcher.getUsedCogMap());
				usedNogMap.putAll(searcher.getUsedNogMap());
				usedKeggMap.putAll(searcher.getUsedKeggMap());
				usedEcMap.putAll(searcher.getUsedEcMap());
				usedGoMap.putAll(searcher.getUsedGoMap());
			}
		}

		File funcFile = new File(metaPar.getDbSearchResultFile(), "functional_annotation");
		if (!funcFile.exists()) {
			funcFile.mkdir();
		}
		File final_pro_xml_file = new File(funcFile, "functions.xml");

		MetaProteinAnnoMag[] metapros = list.toArray(new MetaProteinAnnoMag[list.size()]);

		MetaProteinXMLWriter2 writer = new MetaProteinXMLWriter2(final_pro_xml_file.getAbsolutePath(),
				MetaConstants.pFindMag, metaPar.getQuantitativeMode(), fileNames);

		writer.addProteins(metapros, usedCogMap, usedNogMap, usedKeggMap, usedGoMap, usedEcMap);

		writer.close();

		if (final_pro_xml_file.exists() && final_pro_xml_file.length() > 0) {
			LOGGER.info(getTaskName() + ": proteins with functional annotations have been exported to "
					+ final_pro_xml_file.getName());
			System.out.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": proteins with functional annotations have been exported to " + final_pro_xml_file.getName());
		} else {
			LOGGER.error(getTaskName() + ": error in proteins functional annotations, task failed");
			System.err.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": error in proteins functional annotations, task failed");
		}

		return metapros;
	}

	protected boolean exportReport() {

		updateMessage(getTaskName() + ": exporting report");

		LOGGER.info(getTaskName() + ": exporting report started");
		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": exporting report started");

		File reportDir = new File(metaPar.getResult(), "report");
		if (!reportDir.exists()) {
			reportDir.mkdirs();
		}

		MetaPythonReportCopyTask reportCopyTask = new MetaPythonReportCopyTask(new File(metaPar.getResult()));
		Thread reportCopyThread = new Thread(reportCopyTask);
		reportCopyThread.start();
		try {
			reportCopyThread.join();
			if (!reportCopyTask.get()) {
				LOGGER.error(getTaskName() + ": error in copying report template");
				System.err.println(
						format.format(new Date()) + "\t" + getTaskName() + ": error in copying report template");
			}
		} catch (Exception e) {
			LOGGER.error(getTaskName() + ": error in copying report template", e);
			System.err.println(format.format(new Date()) + "\t" + getTaskName() + ": error in copying report template: "
					+ e.getMessage());
		}

		metaPar.setReportDir(reportDir);

		File reportHtmlDir = new File(reportDir + "\\reports\\rmd_html");

		if (!reportDir.exists()) {
			reportDir.mkdirs();
		}

		File summaryMetaFile = new File(metaPar.getDbSearchResultFile(), "metainfo.tsv");
		try {
			metaPar.getMetadata().exportMetadata(summaryMetaFile);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		MetaPythonReportTask task = new MetaPythonReportTask(metaPar.getThreadCount(), mrp);

		boolean summayrFinished = this.exportSummary(task);

		updateProgress(82, 100);

		File final_pep_txt = new File(metaPar.getResult(), "final_peptides.tsv");

		if (!final_pep_txt.exists() || final_pep_txt.length() == 0) {
//			this.exportPeptideTxt();
		}

		boolean pepReportOk = false;
		if (!final_pep_txt.exists() || final_pep_txt.length() == 0) {
			return false;
		} else {
			File report_peptides_summary = new File(reportHtmlDir, "report_peptides_summary.html");
			if (!report_peptides_summary.exists() || report_peptides_summary.length() == 0) {

				if (summaryMetaFile.exists()) {
					pepReportOk = task.addTask(MetaPythonReportTask.PEPTIDE, final_pep_txt, report_peptides_summary,
							summaryMetaFile);
				} else {
					pepReportOk = task.addTask(MetaPythonReportTask.PEPTIDE, final_pep_txt, report_peptides_summary);
				}
			}
		}

		updateProgress(84, 100);

		File final_pro_txt = new File(metaPar.getResult(), "final_proteins.tsv");

		if (!final_pro_txt.exists() || final_pro_txt.length() == 0) {
//			this.exportProteinsTxt();
		}

		boolean proReportOk = false;
		if (!final_pro_txt.exists() || final_pro_txt.length() == 0) {
			return false;
		} else {
			File report_proteinGroups_summary = new File(reportHtmlDir, "report_proteinGroups_summary.html");
			if (!report_proteinGroups_summary.exists() || report_proteinGroups_summary.length() == 0) {
				if (summaryMetaFile.exists()) {
					proReportOk = task.addTask(MetaPythonReportTask.PROTEIN, final_pro_txt,
							report_proteinGroups_summary, summaryMetaFile);
				} else {
					proReportOk = task.addTask(MetaPythonReportTask.PROTEIN, final_pro_txt,
							report_proteinGroups_summary);
				}
			}
		}

		updateProgress(86, 100);

		boolean taxFuncReportOk = false;
		if (metaPar.isMetaWorkflow()) {
//			taxFuncReportOk = exportTaxaFunc(task);
		} else {
			taxFuncReportOk = true;
		}

		if (!summayrFinished || !pepReportOk || !proReportOk || !taxFuncReportOk) {
			updateProgress(99, 100);

			LOGGER.info(getTaskName() + ": error in generating the reports");
			System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": error in generating the reports");

			return false;
		}

		task.run();

		try {
			boolean finish = task.get();
			if (finish) {

				LOGGER.info(getTaskName() + ": report generation successes");
				System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": report generation successes");

				updateProgress(100, 100);
				return true;
			}

		} catch (InterruptedException | ExecutionException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		return false;
	}

	protected boolean exportSummary(MetaPythonReportTask task) {

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

			BufferedReader idReader = new BufferedReader(
					new FileReader(new File(metaPar.getDbSearchResultFile(), "id.tsv")));
			HashMap<String, int[]> idMap = new HashMap<String, int[]>();
			String line = idReader.readLine();
			while ((line = idReader.readLine()) != null) {
				String[] cs = line.split("\t");
				idMap.put(cs[0], new int[] { Integer.parseInt(cs[1]), Integer.parseInt(cs[2]) });
			}
			idReader.close();

			File ms2CountFile = new File(spectraDirFile, "ms2Count.txt");
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
			LOGGER.error(getTaskName() + ": error in writing PSM information to " + final_summary_txt.getName(), e);
			System.err.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": error in writing PSM information to " + final_summary_txt.getName());
		}

		if (!final_summary_txt.exists() || final_summary_txt.length() == 0) {
			return false;

		} else {
			boolean finish = false;
			File report_ID_summary = new File(metaPar.getReportDir(), "report_ID_summary.html");
			if (summaryMetaFile != null && summaryMetaFile.exists()) {
				finish = task.addTask(MetaPythonReportTask.ID_SUMMARY, final_summary_txt, report_ID_summary,
						summaryMetaFile);
			} else {
				finish = task.addTask(MetaPythonReportTask.ID_SUMMARY, final_summary_txt, report_ID_summary);
			}

			return finish;
		}
	}

}
