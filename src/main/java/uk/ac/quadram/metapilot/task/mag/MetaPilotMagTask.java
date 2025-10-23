package uk.ac.quadram.metapilot.task.mag;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.lang.management.ManagementFactory;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import uk.ac.quadram.metapilot.task.io.ReportManifestGenerator;
import uk.ac.quadram.metapilot.task.MetaPythonReportTask;
import java.sql.SQLException;
import java.text.DecimalFormat;
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
import org.apache.logging.log4j.Logger;

import com.sun.management.OperatingSystemMXBean;

import org.apache.logging.log4j.LogManager;

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
import uk.ac.quadram.metapilot.task.MetaAbstractTask;
import uk.ac.quadram.metapilot.task.io.MetaTreeHandler;
import uk.ac.quadram.metapilot.task.io.pep.MetaPeptide;
import uk.ac.quadram.metapilot.task.io.pep.MetaPSM;
import uk.ac.quadram.metapilot.task.io.pro.MetaProtein;
import uk.ac.quadram.metapilot.task.io.pro.MetaProteinXMLReader2;
import uk.ac.quadram.metapilot.task.io.pro.MetaProteinXMLWriter2;
import uk.ac.quadram.metapilot.task.par.MetaBaseParameter;
import uk.ac.quadram.metapilot.task.par.MetaConstants;
import uk.ac.quadram.metapilot.task.par.MetaData;

public abstract class MetaPilotMagTask extends MetaPilotTask {

	// Logger is inherited from MetaPilotTask as an instance field

	protected String[] rawFiles;
	protected String[] rawFileNames;
	protected String[] quanExpNames;
	protected File spectraDirFile;
	protected MagDbItem magDb;
	protected File[] separateResultFolders;
	protected int threadCount;
	protected double psmCoverThreshold = 0.95;

	protected File rawDir;
	protected File quan_pep_file;
	protected File quan_pro_file;
	protected File final_pep_txt;
	protected File final_pro_txt;
	protected HashMap<String, String> expNameMap;
	protected File summaryMetaFile;
	protected File combinedFile;

	protected int maxTaskCount;
	
	public MetaPilotMagTask(MetaBaseParameter metaPar, MetaResourceParameters mrp, ProgressBar progressBar,
			ProgressBar secondBar) {
		super(metaPar, mrp, progressBar, secondBar);
		this.threadCount = metaPar.getThreadCount();

		setTaskCount();
		logger.info("Task count set to: {}", maxTaskCount);

		logger.info(getTaskName() + ": Constructor completed, report templates will be copied during execution");
	}

	protected void setTaskCount() {
		long totalMemorySize = ((OperatingSystemMXBean) ManagementFactory.getOperatingSystemMXBean())
				.getTotalMemorySize() / 1024 / 1024 / 1024;
		this.maxTaskCount = (int) (totalMemorySize / 20);
		if (this.maxTaskCount == 0) {
			this.maxTaskCount = 1;
		}

		System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": the number of max memory is "
				+ totalMemorySize + " GB, the max number of tasks is set as " + maxTaskCount);
	}
	
	protected void createSeparateFolders() {
		MetaData metaData = metaPar.getMetadata();
		if (metaData == null) {
			logger.error("No metadata");
			throw new IllegalArgumentException("Metadata are required");
		}
		this.rawFiles = metaPar.getMetadata().getRawFiles();
		if (rawFiles == null || rawFiles.length == 0) {
			logger.error("No raw files found in metadata");
			throw new IllegalArgumentException("Raw files are required");
		}

		this.rawDir = (new File(rawFiles[0])).getParentFile();
		if (!rawDir.exists()) {
			logger.error("No raw files found in metadata");
			throw new IllegalArgumentException("Raw files are required");
		}

		logger.info("Initializing result directories in: {}", metaPar.getDbSearchResultFile().getAbsolutePath());

		this.quanExpNames = metaPar.getMetadata().getExpNames();
		this.combinedFile = new File(metaPar.getDbSearchResultFile(), "combined");
		if (!this.combinedFile.exists()) {
			this.combinedFile.mkdir();
		}

		this.rawFileNames = new String[rawFiles.length];
		this.separateResultFolders = new File[rawFiles.length];

		for (int i = 0; i < rawFiles.length; i++) {
			rawFileNames[i] = rawFiles[i].substring(rawFiles[i].lastIndexOf("\\") + 1, rawFiles[i].lastIndexOf("."));
			this.separateResultFolders[i] = new File(metaPar.getDbSearchResultFile(), rawFileNames[i]);
			if (!this.separateResultFolders[i].exists()) {
				this.separateResultFolders[i].mkdir();
				logger.info("  Created result directory: {}", this.separateResultFolders[i].getAbsolutePath());
			} else {
				logger.info("  Result directory already exists: {}", this.separateResultFolders[i].getAbsolutePath());
			}

			File hapFile = new File(separateResultFolders[i], "hap");
			if (!hapFile.exists()) {
				hapFile.mkdirs();
				logger.info("  Created HAP result directory: {}", hapFile.getAbsolutePath());
			} else {
				logger.info(" HAP result directory already exists: {}", hapFile.getAbsolutePath());
			}
		}

		this.spectraDirFile = metaPar.getSpectraFile();
	}

	/**
	 * Advanced genome refinement with fine-tuned parameters. Provides additional
	 * control over genome selection process.
	 * 
	 * USED BY: NONE - This advanced version is not used by any child class
	 * 
	 * @param genomePepMap  Map of genome names to their associated peptides
	 * @param totalPepCount Total number of unique peptides
	 * @param threshold     Coverage threshold for genome selection
	 * @param stepThreshold Step threshold for additional filtering
	 * @param maxium        Maximum number of genomes to select
	 * @param decoyFilter   Whether to apply decoy filtering
	 * @param out           Output file for genome statistics
	 * @return Set of selected genome names
	 */
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
			logger.error(getTaskName() + ": error in writing genomes to " + out, e);
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
	protected HashMap<String, Double> refineGenomes(MetaPeptide[] quanPeps) {
		return refineGenomes(quanPeps, PFindResultReader.REV);
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
	protected HashMap<String, Double> refineGenomes(MetaPeptide[] quanPeps, String decoySymble) {

		String revIdentifier = decoySymble + this.magDb.getIdentifier();
		HashMap<String, Double> genomeIntensityMap = new HashMap<String, Double>();
		double[] peptideIntensity = new double[quanPeps.length];
		logger.info(getTaskName() + ": quanPeps: " + quanPeps.length);
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
			logger.info(totalIntensity+"\t"+idenCount);
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

		logger.info(getTaskName() + ": genome before filtering: " + genomeIntensityMap.size());

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

		logger.info(getTaskName() + ": razor genome before filtering: " + razorGenomeIntenMap.size());

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

		logger.info(getTaskName() + ": decoy genome count: " + revLogIntenList.size());

		HashMap<String, Double> genomeScoreMap = new HashMap<String, Double>();

		if (revLogIntenList.size() == 0) {
			logger.info(getTaskName()
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
			logger.info(getTaskName() + ": refined decoy genome:" + revLogIntensities.length);
			logger.info(Arrays.toString(revLogIntensities));

			if (revLogIntensities.length == 0) {
				logger.info(getTaskName() + ": the number of decoy genome was not sufficient for the t-test, "
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

				logger.info(getTaskName() + ": decoy genome used for scoring: " + extended.length);

				for (int i = 0; i < razorGenomes.length; i++) {
					double intensity = razorGenomeIntenMap.get(razorGenomes[i]);

					if (razorGenomes[i].startsWith(revIdentifier) || intensity == 0) {
						break;
					}

					double pvalue = tTest.tTest(Math.log10(intensity), extended);
					logger.info(getTaskName() + ": pvalue for genome " + razorGenomes[i] + ": " + pvalue);

					if (pvalue < 0.01) {
						genomeScoreMap.put(razorGenomes[i], pvalue);
					} else {
//						break;
					}
				}
			} else {
				logger.info(getTaskName() + ": decoy genome used for scoring: " + revLogIntensities.length);

				for (int i = 0; i < razorGenomes.length; i++) {
					double intensity = razorGenomeIntenMap.get(razorGenomes[i]);

					if (razorGenomes[i].startsWith(revIdentifier) || intensity == 0) {
						break;
					}

					double pvalue = tTest.tTest(Math.log10(intensity), revLogIntensities);
					logger.info(getTaskName() + ": pvalue for genome " + razorGenomes[i] + ": " + pvalue);

					if (pvalue < 0.01) {
						genomeScoreMap.put(razorGenomes[i], pvalue);
					} else {
//						break;
					}
				}
			}
		}

		logger.info(getTaskName() + ": genome after filtering: " + genomeScoreMap.size());

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
	protected String[] refineProteins(MetaPeptide[] quanPeps, HashMap<String, Double> genomeScoreMap) {
		return refineProteins(quanPeps, genomeScoreMap, PFindResultReader.REV);
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
			logger.error(getTaskName() + ": error in genome analysis, task failed", e);
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
			logger.error(getTaskName() + ": Error in searching the genome2taxa database", e);
			System.err.println(
					format.format(new Date()) + "\t" + getTaskName() + ": error in searching the genome2taxa database");
			return false;
		}

		File genomeHtmlFile = new File(taxFile, "Genome.html");
		try {
			HMGenomeHtmlWriter htmlWriter = new HMGenomeHtmlWriter(genomeHtmlFile);
			htmlWriter.write(genomeTaxaMap);
			htmlWriter.close();

			// Copy the generated Genome.html to the reports directory
			File reportGenomeHtml = new File(metaPar.getReportDir(), "reports/Genome.html");
			if (reportGenomeHtml.getParentFile() != null && !reportGenomeHtml.getParentFile().exists()) {
				reportGenomeHtml.getParentFile().mkdirs();
			}
			Files.copy(genomeHtmlFile.toPath(), reportGenomeHtml.toPath(), StandardCopyOption.REPLACE_EXISTING);
			logger.info("Copied Genome.html to reports directory: " + reportGenomeHtml.getAbsolutePath());

			// Add to manifest
			File reportsDir = new File(metaPar.getReportDir(), "reports");
			ReportManifestGenerator manifestGen = new ReportManifestGenerator(reportsDir);
			manifestGen.addGenomeReport(reportGenomeHtml);

		} catch (IOException e) {
			// TODO Auto-generated catch block
			logger.error("Error in writing genome information to " + genomeHtmlFile, e);
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
			logger.error("Error in writing genome information to " + genomeFile, e);
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
			logger.error("Error in writing taxa information to " + taxaFile, e);
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

		File jsOutput = new File(metaPar.getReportDir() + "\\reports\\tree.js");
		if (!jsOutput.getParentFile().exists()) {
			jsOutput.getParentFile().mkdirs();
		}

		if (!jsOutput.exists() || jsOutput.length() == 0) {
			MetaTreeHandler.exportJS(taxRankMaps, taxIntensityMaps, cellularIntensity, expNames, jsOutput);
		}

		// Generate bubble chart report from tree.js data (optional - don't fail if this
		// doesn't work)
		File bubbleChartOutput = new File(metaPar.getReportDir() + "\\reports\\bubble_chart_report.html");
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

		logger.info(getTaskName() + ": taxonomy analysis started");
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
				logger.error(getTaskName() + ": error in proteins functional annotations, task failed", e);
				System.err.println(format.format(new Date()) + "\t" + getTaskName()
						+ ": error in proteins functional annotations, task failed");

				return null;
			}

			MetaProteinAnnoMag[] mpas = null;
			try {
				mpas = searcher.match(mps);
			} catch (IOException | SQLException e) {
				// TODO Auto-generated catch block
				logger.error(getTaskName() + ": error in proteins functional annotations, task failed", e);
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
			logger.info(getTaskName() + ": proteins with functional annotations have been exported to "
					+ final_pro_xml_file.getName());
			System.out.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": proteins with functional annotations have been exported to " + final_pro_xml_file.getName());
		} else {
			logger.error(getTaskName() + ": error in proteins functional annotations, task failed");
			System.err.println(format.format(new Date()) + "\t" + getTaskName()
					+ ": error in proteins functional annotations, task failed");
		}

		return metapros;
	}

	/**
	 * Generate identification summary report
	 */
	private boolean generateIdentificationSummary(HashMap<String, Integer> fileIdMap, int[] filePsmCount,
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

	/**
	 * Generate peptide report
	 */
	private boolean generatePeptideReport(ArrayList<MetaPeptide> peptideList, HashMap<String, MetaPSM> psmMap,
			HashMap<String, int[]> pepPsmCountMap, HashMap<String, int[]> pepChargeMap,
			HashMap<String, Integer> proteinIdMap, ArrayList<String> fileList, String quanMode, String[] intensityTitle,
			String taskName) {

		logger.info(taskName + ": exporting peptide report started");
		System.out.println(format.format(new Date()) + "\t" + taskName + ": exporting peptide report started");

		try {
			final_pep_txt = new File(this.metaPar.getDbSearchResultFile(), "final_peptides.tsv");
			PrintWriter writer = new PrintWriter(final_pep_txt);

			// Build header
			StringBuilder titlesb = new StringBuilder();
			titlesb.append("Sequence").append("\t");
			titlesb.append("Base sequence").append("\t");
			titlesb.append("Length").append("\t");
			titlesb.append("Missed cleavages").append("\t");
			titlesb.append("Mass").append("\t");
			titlesb.append("Proteins").append("\t");
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

			if (quanMode.equals(MetaConstants.labelFree)) {
				for (int i = 0; i < intensityTitle.length; i++) {
					titlesb.append("Intensity " + intensityTitle[i]).append("\t");
				}
			}

			titlesb.append("Reverse").append("\t");
			titlesb.append("Potential contaminant").append("\t");
			titlesb.append("id").append("\t");
			titlesb.append("Protein group IDs").append("\t");
			titlesb.append("MS/MS Count");

			writer.println(titlesb);

			// Write peptide data
			int id = 0;
			for (MetaPeptide peptide : peptideList) {
				String fullSeq = peptide.getModSeq();
				String sequence = peptide.getSequence();
				MetaPSM psm = psmMap.get(fullSeq);
				// Note: This would need to be adapted based on the specific PSM type
				// For now, using basic peptide information
				if (psm == null) {
					logger.info(getTaskName() + ": no PSM was found for peptide " + fullSeq);
					continue;
				}
				StringBuilder sb = new StringBuilder();
				sb.append(fullSeq).append("\t");
				sb.append(sequence).append("\t");
				sb.append(sequence.length()).append("\t");
				sb.append(psm.getMiss()).append("\t");
				sb.append(psm.getCalMass()).append("\t");

				String[] pros = peptide.getProteins();
				for (String pro : pros) {
					if (proteinIdMap.containsKey(pro)) {
						sb.append(pro).append(";");
					}
				}
				if (sb.length() > 0 && sb.charAt(sb.length() - 1) == ';') {
					sb.deleteCharAt(sb.length() - 1);
				}
				sb.append("\t");

				// Charges
				int[] charge = pepChargeMap.get(fullSeq);
				if (charge != null) {
					StringBuilder chargesb = new StringBuilder();
					for (int j = 0; j < charge.length; j++) {
						if (charge[j] == 1) {
							chargesb.append(j).append(";");
						}
					}
					if (chargesb.length() > 1) {
						sb.append(chargesb.subSequence(0, chargesb.length() - 1)).append("\t");
					} else {
						sb.append("2").append("\t");
					}
				} else {
					sb.append("2").append("\t");
				}

				sb.append(psm.getScore()).append("\t");
				sb.append(psm.getPEP()).append("\t");

				// Identification types and MS/MS counts
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
				}

				// Intensities
				double totalIntensity = 0;
				double[] intensity = peptide.getIntensity();
				for (int j = 0; j < intensity.length; j++) {
					totalIntensity += intensity[j];
				}

				sb.append(java.math.BigInteger.valueOf((long) Math.floor(totalIntensity))).append("\t");

				for (int j = 0; j < intensity.length; j++) {
					sb.append(java.math.BigInteger.valueOf((long) Math.floor(intensity[j]))).append("\t");
				}

				sb.append("").append("\t"); // Reverse
				sb.append("").append("\t"); // Potential contaminant
				sb.append(id).append("\t"); // id

				// Protein group IDs
				for (String pro : pros) {
					if (proteinIdMap.containsKey(pro)) {
						sb.append(proteinIdMap.get(pro)).append(";");
					}
				}
				if (sb.length() > 0 && sb.charAt(sb.length() - 1) == ';') {
					sb.deleteCharAt(sb.length() - 1);
				}
				sb.append("\t");

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

	/**
	 * Generate protein report
	 */
	private boolean generateProteinReport(String[] proteins, HashMap<String, double[]> proteinIntensityMap,
			HashMap<String, Double> proteinTotalIntenMap, HashMap<String, HashSet<String>> totalPepCountMap,
			HashMap<String, HashSet<String>> razorPepCountMap, HashMap<String, Integer> psmCountMap,
			HashMap<String, Integer> proteinIdMap, String quanMode, String[] intensityTitle, String taskName) {

		logger.info(taskName + ": exporting protein report started");
		System.out.println(format.format(new Date()) + "\t" + taskName + ": exporting protein report started");

		try {
			final_pro_txt = new File(this.metaPar.getDbSearchResultFile(), "final_proteins.tsv");
			PrintWriter writer = new PrintWriter(final_pro_txt);

			// Build header
			StringBuilder titlesb = new StringBuilder();
			titlesb.append("Protein IDs").append("\t");
			titlesb.append("Peptide counts (all)").append("\t");
			titlesb.append("Peptide counts (razor)").append("\t");
			titlesb.append("Number of PSMs").append("\t");
			titlesb.append("PEP").append("\t");
			titlesb.append("Intensity").append("\t");

			if (quanMode.equals(MetaConstants.labelFree)) {
				for (int i = 0; i < intensityTitle.length; i++) {
					titlesb.append("Intensity " + intensityTitle[i]).append("\t");
				}
			} else if (quanMode.equals(MetaConstants.isobaricLabel)) {
				for (int i = 0; i < intensityTitle.length; i++) {
					titlesb.append("Intensity " + intensityTitle[i]).append("\t");
				}
			}

			titlesb.append("Reverse").append("\t");
			titlesb.append("Potential contaminant").append("\t");
			titlesb.append("id").append("\t");
			titlesb.append("Peptide IDs").append("\t");
			titlesb.append("Peptide is razor");

			writer.println(titlesb);

			// Write protein data
			DecimalFormat df = FormatTool.getDFE4();
			for (int i = 0; i < proteins.length; i++) {
				StringBuilder sb = new StringBuilder();
				sb.append(proteins[i]).append("\t");
				sb.append(totalPepCountMap.get(proteins[i]).size()).append("\t");
				sb.append(razorPepCountMap.get(proteins[i]).size()).append("\t");
				sb.append(psmCountMap.get(proteins[i])).append("\t");
				sb.append("0").append("\t"); // PEP - would need to be calculated

				double totalIntensity = proteinTotalIntenMap.get(proteins[i]);
				sb.append(df.format(totalIntensity)).append("\t");

				double[] intensity = proteinIntensityMap.get(proteins[i]);
				for (int j = 0; j < intensity.length; j++) {
					sb.append(df.format(intensity[j])).append("\t");
				}
				sb.append("\t"); // Reverse
				sb.append("\t"); // Potential contaminant
				sb.append(i).append("\t"); // id

				// Note: Peptide IDs and razor information would need additional processing
				sb.append("").append("\t"); // Peptide IDs
				sb.append(""); // Peptide is razor

				writer.println(sb);
			}

			writer.close();

		} catch (IOException e) {
			logger.error(taskName + ": error in writing final protein result to " + final_pro_txt.getName(), e);
			System.err.println(format.format(new Date()) + "\t" + taskName
					+ ": error in writing final protein result to " + final_pro_txt.getName());
			return false;
		}

		logger.info(taskName + ": protein report has been exported to " + final_pro_txt.getName());
		System.out.println(format.format(new Date()) + "\t" + taskName + ": protein report has been exported to "
				+ final_pro_txt.getName());

		return true;
	}

	/**
	 * Perform functional annotation
	 */
	private boolean performFunctionalAnnotation(String[] proteins, HashMap<String, double[]> proteinIntensityMap,
			String[] intensityTitle, String taskName) {

		logger.info(taskName + ": functional annotation started");
		System.out.println(format.format(new Date()) + "\t" + taskName + ": functional annotation started");

		File funcFile = new File(metaPar.getDbSearchResultFile(), "functional_annotation");
		if (!funcFile.exists()) {
			funcFile.mkdir();
		}

		// Convert protein data for functional annotation
		ArrayList<MetaProtein> proteinList = new ArrayList<MetaProtein>();
		for (int i = 0; i < proteins.length; i++) {
			double[] intensity = proteinIntensityMap.get(proteins[i]);
			MetaProtein mp = new MetaProtein(i, 0, proteins[i], 0, 0, 0, intensity);
			proteinList.add(mp);
		}

		MetaProtein[] metapros = proteinList.toArray(new MetaProtein[proteinList.size()]);

		// Perform functional annotation (this writes the XML file)
		MetaProteinAnnoMag[] mpas = funcAnnoSqlite(metapros, intensityTitle);

		if (mpas == null) {
			logger.error(taskName + ": Functional annotation failed");
			System.err.println(format.format(new Date()) + "\t" + taskName + ": Functional annotation failed");
			return false;
		}

		// Generate functional annotation reports from the XML file
		File funTsv = new File(funcFile, "functions.tsv");
		File final_pro_xml_file = new File(funcFile, "functions.xml");

		if (!funTsv.exists()) {
			System.out.println(format.format(new Date()) + "\t" + taskName + ": writing results to " + funTsv.getName()
					+ " started");
			try {
				MetaProteinXMLReader2 funReader = new MetaProteinXMLReader2(final_pro_xml_file);
				funReader.exportTsv(funTsv);
				// Export separate functional category reports
				File goReportTsv = new File(funcFile, "Gene_Ontology_id_report.tsv");
				File ecReportTsv = new File(funcFile, "EC_id_report.tsv");
				File keggReportTsv = new File(funcFile, "KEGG_ko_report.tsv");
				File briteReportTsv = new File(funcFile, "BRITE_report.tsv");
				File cogReportTsv = new File(funcFile, "COG_accession_report.tsv");
				File nogReportTsv = new File(funcFile, "NOG_accession_report.tsv");
				File cogCatTsv = new File(funcFile, "COG_category_report.tsv");
				File nogCatTsv = new File(funcFile, "NOG_category_report.tsv");
				
				funReader.exportGOReport(goReportTsv);
				funReader.exportECReport(ecReportTsv);
				funReader.exportKeggReport(keggReportTsv);
				funReader.exportBriteReport(briteReportTsv);
				funReader.exportCOGReport(cogReportTsv);
				funReader.exportNOGReport(nogReportTsv);
				funReader.exportCOGCategoryReport(cogCatTsv);
				funReader.exportNOGCategoryReport(nogCatTsv);				
			} catch (Exception e) {
				logger.error(taskName + ": error in writing protein function information to " + funTsv.getName(), e);
				System.err.println(format.format(new Date()) + "\t" + taskName
						+ ": error in writing protein function information to " + funTsv.getName());
				return false;
			}
			System.out.println(format.format(new Date()) + "\t" + taskName + ": writing results to " + funTsv.getName()
					+ " finished");
		}

		File html = new File(funcFile, "functions.html");
		if (!html.exists()) {
			System.out.println(format.format(new Date()) + "\t" + taskName + ": writing results to " + html.getName()
					+ " started");
			try {
				MetaProteinXMLReader2 funReader = new MetaProteinXMLReader2(final_pro_xml_file);
				funReader.exportHtml(html);
			} catch (Exception e) {
				logger.error(taskName + ": error in writing protein function information to " + html.getName(), e);
				System.err.println(format.format(new Date()) + "\t" + taskName
						+ ": error in writing protein function information to " + html.getName());
				return false;
			}
			System.out.println(format.format(new Date()) + "\t" + taskName + ": writing results to " + html.getName()
					+ " finished");
		}

		logger.info(taskName + ": functional annotation finished");
		System.out.println(format.format(new Date()) + "\t" + taskName + ": functional annotation finished");

		return true;
	}

	/**
	 * Perform taxonomy analysis
	 */
	private boolean performTaxonomyAnalysis(HashMap<String, double[]> genomeIntensityMap,
			HashMap<String, Double> genomeScoreMap, HashMap<String, Integer> psmCountMap,
			HashMap<String, HashSet<String>> totalPepCountMap, HashMap<String, HashSet<String>> razorPepCountMap,
			String[] intensityTitle, String taskName) {

		logger.info(taskName + ": taxonomy analysis started");
		System.out.println(format.format(new Date()) + "\t" + taskName + ": taxonomy analysis started");

		File taxFile = new File(metaPar.getDbSearchResultFile(), "taxonomy_analysis");
		if (!taxFile.exists()) {
			taxFile.mkdir();
		}

		MagFuncSearcher searcher = null;
		try {
			searcher = new MagFuncSearcher(metaPar.getMicroDb(), mrp.getToolPath(ToolType.FUNCDB));
		} catch (NumberFormatException | SQLException | IOException e) {
			logger.error(taskName + ": error in genome analysis, task failed", e);
			System.err.println(format.format(new Date()) + "\t" + taskName + ": error in genome analysis, task failed");
			return false;
		}

		String[] genomes = genomeIntensityMap.keySet().toArray(new String[genomeIntensityMap.size()]);
		HashMap<String, String[]> genomeTaxaMap = null;
		try {
			genomeTaxaMap = searcher.matchGenome(genomes);
		} catch (IOException e) {
			logger.error(taskName + ": Error in searching the genome2taxa database", e);
			System.err.println(
					format.format(new Date()) + "\t" + taskName + ": error in searching the genome2taxa database");
			return false;
		}

		File genomeHtmlFile = new File(taxFile, "Genome.html");
		try {
			HMGenomeHtmlWriter htmlWriter = new HMGenomeHtmlWriter(genomeHtmlFile);
			htmlWriter.write(genomeTaxaMap);
			htmlWriter.close();
		} catch (IOException e) {
			logger.error("Error in writing genome information to " + genomeHtmlFile, e);
			System.err.println(
					format.format(new Date()) + "\t" + "Error in writing genome information to " + genomeHtmlFile);
			return false;
		}

		TaxonomyRanks[] mainRanks = TaxonomyRanks.getMainRanks();
		double[] cellularIntensity = new double[intensityTitle.length];
		HashMap<String, String[]>[] taxRankMaps = new HashMap[mainRanks.length];
		HashMap<String, double[]>[] taxIntensityMaps = new HashMap[mainRanks.length];
		HashMap<String, Double>[] taxScoreMaps = new HashMap[mainRanks.length];
		for (int i = 0; i < taxRankMaps.length; i++) {
			taxRankMaps[i] = new HashMap<String, String[]>();
			taxIntensityMaps[i] = new HashMap<String, double[]>();
			taxScoreMaps[i] = new HashMap<String, Double>();
		}

		File genomeFile = new File(taxFile, "Genome.tsv");
		PrintWriter genomeWriter = null;
		try {
			genomeWriter = new PrintWriter(genomeFile);
			StringBuilder genomeTitlesb = new StringBuilder();
			genomeTitlesb.append("Genome\t");
			genomeTitlesb.append("Minus_log(p-value)\t");
			genomeTitlesb.append("PSM count\t");
			genomeTitlesb.append("Peptide count\t");
			genomeTitlesb.append("Razor peptide count\t");
			genomeTitlesb.append("GTDB_Name\t");
			genomeTitlesb.append("GTDB_Rank\t");

			for (int i = 0; i < mainRanks.length; i++) {
				genomeTitlesb.append("GTDB_").append(mainRanks[i].getName()).append("\t");
			}

			genomeTitlesb.append("NCBI_Name\t");
			genomeTitlesb.append("NCBI_Rank");

			for (int i = 0; i < mainRanks.length; i++) {
				genomeTitlesb.append("\t").append("NCBI_").append(mainRanks[i].getName());
			}

			for (String exp : intensityTitle) {
				if (exp.startsWith("Intensity ")) {
					genomeTitlesb.append("\t").append(exp);
				} else {
					genomeTitlesb.append("\tIntensity ").append(exp);
				}
			}

			genomeWriter.println(genomeTitlesb);

			for (int i = 0; i < genomes.length; i++) {
				String[] taxaContent = genomeTaxaMap.get(genomes[i]);
				if (genomeIntensityMap.containsKey(genomes[i])) {
					double[] genomeIntensity = genomeIntensityMap.get(genomes[i]);
					for (int j = 0; j < genomeIntensity.length; j++) {
						cellularIntensity[j] += genomeIntensity[j];
					}

					double score = genomeScoreMap.get(genomes[i]);
					StringBuilder genomeTaxaSb = new StringBuilder();
					genomeTaxaSb.append(genomes[i]).append("\t");
					genomeTaxaSb.append(FormatTool.getDF2().format(-Math.log10(score))).append("\t");

					if (psmCountMap.containsKey(genomes[i])) {
						genomeTaxaSb.append(psmCountMap.get(genomes[i])).append("\t");
					} else {
						genomeTaxaSb.append("0").append("\t");
					}
					if (totalPepCountMap.containsKey(genomes[i])) {
						genomeTaxaSb.append(totalPepCountMap.get(genomes[i]).size()).append("\t");
					} else {
						genomeTaxaSb.append("0").append("\t");
					}
					if (razorPepCountMap.containsKey(genomes[i])) {
						genomeTaxaSb.append(razorPepCountMap.get(genomes[i]).size());
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
							taxScoreMaps[j].put(taxon, taxScoreMaps[j].get(taxon) - Math.log10(score));
						} else {
							taxScoreMaps[j].put(taxon, -Math.log10(score));
						}
					}
				}
			}
			genomeWriter.close();

		} catch (FileNotFoundException e) {
			logger.error("Error in writing genome information to " + genomeFile, e);
			System.err
					.println(format.format(new Date()) + "\t" + "Error in writing genome information to " + genomeFile);
			return false;
		}

		File taxaFile = new File(taxFile, "Taxa.tsv");
		File taxaReportFile = new File(taxFile, "Taxa_report.tsv");
		PrintWriter taxaWriter = null;
		PrintWriter taxaReportWriter = null;
		try {
			taxaWriter = new PrintWriter(taxaFile);
			taxaReportWriter = new PrintWriter(taxaReportFile);

			StringBuilder taxaTitlesb1 = new StringBuilder();
			taxaTitlesb1.append("Name\t");
			taxaTitlesb1.append("Minus_log(p-value)\t");
			taxaTitlesb1.append("Rank");

			for (int i = 0; i < mainRanks.length; i++) {
				taxaTitlesb1.append("\t").append(mainRanks[i].getName());
			}
			for (String exp : intensityTitle) {
				if (exp.startsWith("Intensity ")) {
					taxaTitlesb1.append("\t").append(exp);
				} else {
					taxaTitlesb1.append("\tIntensity ").append(exp);
				}
			}
			StringBuilder taxaTitlesb2 = new StringBuilder();
			taxaTitlesb2.append("Name\t");
			taxaTitlesb2.append("Rank");

			for (int i = 0; i < mainRanks.length; i++) {
				taxaTitlesb2.append("\t").append(mainRanks[i].getName());
			}
			for (String exp : intensityTitle) {
				if (exp.startsWith("Intensity ")) {
					taxaTitlesb2.append("\t").append(exp);
				} else {
					taxaTitlesb2.append("\tIntensity ").append(exp);
				}
			}
			taxaWriter.println(taxaTitlesb1);
			taxaReportWriter.println(taxaTitlesb2);

			for (int i = 0; i < taxRankMaps.length; i++) {
				for (String taxon : taxRankMaps[i].keySet()) {
					String[] taxaContent = taxRankMaps[i].get(taxon);
					double[] intensity = taxIntensityMaps[i].get(taxon);

					StringBuilder taxaSb = new StringBuilder();
					taxaSb.append(taxon).append("\t");
					double taxonScore = taxScoreMaps[i].get(taxon);
					taxaSb.append(FormatTool.getDF2().format(taxonScore)).append("\t");
					taxaSb.append(mainRanks[i].getName());

					for (int j = 0; j < mainRanks.length; j++) {
						if (j < taxaContent.length) {
							taxaSb.append("\t").append(taxaContent[j]);
						} else {
							taxaSb.append("\t");
						}
					}

					for (int j = 0; j < intensity.length; j++) {
						taxaSb.append("\t").append(intensity[j]);
					}

					taxaWriter.println(taxaSb);
				}
			}

			taxaWriter.close();

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
					for (int j = 0; j < intensity.length; j++) {
						taxaSb.append("\t").append(intensity[j]);
					}

					taxaReportWriter.println(taxaSb);
				}
			}

			taxaReportWriter.close();

		} catch (FileNotFoundException e) {
			logger.error("Error in writing taxa information to " + taxaFile, e);
			System.err.println(format.format(new Date()) + "\t" + "Error in writing taxa information to " + taxaFile);
			return false;
		}

		// Export JS for taxonomy visualization
		File jsOutput = new File(metaPar.getReportDir() + "\\reports\\tree.js");
		if (!jsOutput.getParentFile().exists()) {
			jsOutput.getParentFile().mkdirs();
		}

		if (!jsOutput.exists() || jsOutput.length() == 0) {
			MetaTreeHandler.exportJS(taxRankMaps, taxIntensityMaps, cellularIntensity, intensityTitle, jsOutput);
		}

		// Generate bubble chart report from tree.js data (optional - don't fail if this
		// doesn't work)
		File bubbleChartOutput = new File(metaPar.getReportDir() + "\\reports\\bubble_chart_report.html");
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

		logger.info(taskName + ": taxonomy analysis finished");
		System.out.println(format.format(new Date()) + "\t" + taskName + ": taxonomy analysis finished");

		return true;
	}

	/**
	 * Generate final reports using MetaPythonReportTask
	 */
	private boolean generateFinalReports(String taskName) {

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
		MetaPythonReportTask reportTask = new MetaPythonReportTask(threadCount, mrp);
		reportTask.setWorkflowType("MAG");
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
		File taxaReportFile = new File(metaPar.getDbSearchResultFile(), "taxonomy_analysis\\Taxa_report.tsv");
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

		File genomeHtmlFile = new File(metaPar.getDbSearchResultFile(), "taxonomy_analysis\\Genome.html");
		if (genomeHtmlFile.exists()) {
			File reportGenomeHtml = new File(reportHtmlDir, "Genome.html");
			if (reportGenomeHtml.getParentFile() != null && !reportGenomeHtml.getParentFile().exists()) {
				reportGenomeHtml.getParentFile().mkdirs();
			}
			try {
				Files.copy(genomeHtmlFile.toPath(), reportGenomeHtml.toPath(), StandardCopyOption.REPLACE_EXISTING);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				logger.error("Copied Genome.html to reports directory failed", e);
			}
			logger.info("Copied Genome.html to reports directory: " + reportGenomeHtml.getAbsolutePath());
		} else {
			logger.info("Genome.html was not exist in " + genomeHtmlFile.getAbsolutePath());
		}

		// Add function report
		File funReportTsv = new File(metaPar.getDbSearchResultFile(), "functional_annotation");
		if (funReportTsv.exists() && funReportTsv.length() > 0) {
			File report_function_summary = new File(reportHtmlDir, "report_function_summary.html");
			if (summaryMetaFile != null && summaryMetaFile.exists()) {
				reportSuccess &= reportTask.addTask(MetaPythonReportTask.FUNCTION, funReportTsv,
						report_function_summary, summaryMetaFile);
			} else {
				reportSuccess &= reportTask.addTask(MetaPythonReportTask.FUNCTION, funReportTsv,
						report_function_summary);
			}
		}

		File funcHtmlFile = new File(metaPar.getDbSearchResultFile(), "functional_annotation\\functions.html");
		if (funcHtmlFile.exists()) {
			File reportFunctionsHtml = new File(reportHtmlDir, "reports/functions.html");
			try {
				Files.copy(funcHtmlFile.toPath(), reportFunctionsHtml.toPath(), StandardCopyOption.REPLACE_EXISTING);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				logger.error("Copied functions.html to reports directory failed", e);
			}
			logger.info("Copied functions.html to reports directory: " + reportFunctionsHtml.getAbsolutePath());
		} else {
			logger.info("functions.html was not exist in " + funcHtmlFile.getAbsolutePath());
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
	 * Generate bubble chart report from tree.js data using Python.
	 */
	private void generateBubbleChartReport(File treeJsFile, File outputFile) {
		try {
			logger.info(getTaskName() + ": Generating bubble chart report from " + treeJsFile.getName());
			System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": Generating bubble chart report");

			// Use MetaPythonReportTask to generate the bubble chart
			MetaPythonReportTask reportTask = new MetaPythonReportTask(1, mrp);
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
	 * Task for post-analysis processing in DDA MAG workflows. Handles protein
	 * inference, filtering, and report generation with proper progress tracking.
	 */
	public static class DDAMagPostAnalysisTask extends MetaAbstractTask {

		private MetaPeptide[] quanPeps;
		private MetaPSM[] psms;
		private HashMap<String, Integer> fileIdMap;
		private ArrayList<String> fileList;
		private String quanMode;
		private String[] intensityTitle;
		private String taskName;
		private MetaPilotMagTask parentTask;
		private Logger logger;

		public DDAMagPostAnalysisTask(MetaBaseParameter metaPar, MetaResourceParameters mrp, ProgressBar progressBar,
				MetaPeptide[] quanPeps, MetaPSM[] psms, HashMap<String, Integer> fileIdMap, ArrayList<String> fileList,
				String quanMode, String[] intensityTitle, String taskName, MetaPilotMagTask parentTask) {
			super(metaPar, mrp, progressBar);
			this.quanPeps = quanPeps;
			this.psms = psms;
			this.fileIdMap = fileIdMap;
			this.fileList = fileList;
			this.quanMode = quanMode;
			this.intensityTitle = intensityTitle;
			this.taskName = taskName;
			this.parentTask = parentTask;
			this.logger = org.apache.logging.log4j.LogManager.getLogger(getClass());
		}

		@Override
		protected String getTaskName() {
			return taskName + " Post-Analysis";
		}

		@Override
		protected Logger getLogger() {
			return logger;
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
			HashMap<String, Integer> psmCountMap = new HashMap<String, Integer>();

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

				// Count PSMs per protein
				String[] pros = psm.getProteins();
				for (String pro : pros) {
					if (psmCountMap.containsKey(pro)) {
						psmCountMap.put(pro, psmCountMap.get(pro) + 1);
					} else {
						psmCountMap.put(pro, 1);
					}

					if (pro.startsWith(parentTask.magDb.getIdentifier())) {
						String genome = pro.split("_")[0];
						if (psmCountMap.containsKey(genome)) {
							psmCountMap.put(genome, psmCountMap.get(genome) + 1);
						} else {
							psmCountMap.put(genome, 1);
						}
					}
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
			updateMessage("Refining genomes and proteins...");
			updateProgress(0.15, 1.0);

			// Refine genomes and proteins
			HashMap<String, Double> genomeScoreMap = parentTask.refineGenomes(quanPeps);
			String[] razorProteins = parentTask.refineProteins(quanPeps, genomeScoreMap);

			// Update progress: Genome and protein refinement complete
			updateMessage("Calculating protein intensities...");
			updateProgress(0.30, 1.0);

			// Initialize intensity maps
			HashMap<String, double[]> genomeIntensityMap = new HashMap<String, double[]>();
			HashMap<String, double[]> proteinIntensityMap = new HashMap<String, double[]>();
			HashMap<String, Double> proteinTotalIntenMap = new HashMap<String, Double>();

			HashSet<String> pepSequenceSet = new HashSet<String>();
			HashMap<String, HashSet<String>> totalPepCountMap = new HashMap<String, HashSet<String>>();
			HashMap<String, HashSet<String>> razorPepCountMap = new HashMap<String, HashSet<String>>();

			ArrayList<MetaPeptide> peptideList = new ArrayList<MetaPeptide>();

			// Process peptides and build protein maps
			for (int i = 0; i < quanPeps.length; i++) {
				if (razorProteins[i] != null) {
					String stripeSequence = quanPeps[i].getSequence();
					pepSequenceSet.add(stripeSequence);
					peptideList.add(quanPeps[i]);

					// Build razor peptide count map
					if (razorPepCountMap.containsKey(razorProteins[i])) {
						razorPepCountMap.get(razorProteins[i]).add(stripeSequence);
					} else {
						HashSet<String> set = new HashSet<String>();
						set.add(stripeSequence);
						razorPepCountMap.put(razorProteins[i], set);
					}

					// Handle MAG database proteins
					if (razorProteins[i].startsWith(parentTask.magDb.getIdentifier())) {
						String genome = razorProteins[i].split("_")[0];
						if (razorPepCountMap.containsKey(genome)) {
							razorPepCountMap.get(genome).add(stripeSequence);
						} else {
							HashSet<String> set = new HashSet<String>();
							set.add(stripeSequence);
							razorPepCountMap.put(genome, set);
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

				if (count > 0) {
					String stripeSequence = quanPeps[i].getSequence();
					double[] pepIntensity = new double[quanPeps[i].getIntensity().length];
					for (int j = 0; j < pepIntensity.length; j++) {
						pepIntensity[j] = quanPeps[i].getIntensity()[j] / (double) count;
					}

					for (int j = 0; j < proteins.length; j++) {
						if (razorPepCountMap.containsKey(proteins[j])) {
							if (proteins[j].startsWith(parentTask.magDb.getIdentifier())) {
								String genome = proteins[j].split("_")[0];
								if (genomeScoreMap.containsKey(genome)) {
									// Update genome intensity
									double[] genomeIntensity;
									if (genomeIntensityMap.containsKey(genome)) {
										genomeIntensity = genomeIntensityMap.get(genome);
									} else {
										genomeIntensity = new double[pepIntensity.length];
										genomeIntensityMap.put(genome, genomeIntensity);
									}
									for (int k = 0; k < pepIntensity.length; k++) {
										genomeIntensity[k] += pepIntensity[k];
									}

									// Update protein intensity
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

									// Update peptide count maps
									if (totalPepCountMap.containsKey(genome)) {
										totalPepCountMap.get(genome).add(stripeSequence);
									} else {
										HashSet<String> set = new HashSet<String>();
										set.add(stripeSequence);
										totalPepCountMap.put(genome, set);
									}

									if (totalPepCountMap.containsKey(proteins[j])) {
										totalPepCountMap.get(proteins[j]).add(stripeSequence);
									} else {
										HashSet<String> set = new HashSet<String>();
										set.add(stripeSequence);
										totalPepCountMap.put(proteins[j], set);
									}
								}
							} else {
								// Handle non-MAG proteins
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
									totalPepCountMap.get(proteins[j]).add(stripeSequence);
								} else {
									HashSet<String> set = new HashSet<String>();
									set.add(stripeSequence);
									totalPepCountMap.put(proteins[j], set);
								}
							}
						}
					}
				}
			}

			// Update progress: Intensity calculations complete
			updateMessage("Sorting proteins and generating reports...");
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
					fileList, quanMode, intensityTitle, taskName)) {
				return false;
			}

			// Generate protein report
			updateMessage("Generating protein report...");
			updateProgress(0.80, 1.0);
			if (!parentTask.generateProteinReport(proteins, proteinIntensityMap, proteinTotalIntenMap, totalPepCountMap,
					razorPepCountMap, psmCountMap, proteinIdMap, quanMode, intensityTitle, taskName)) {
				return false;
			}

			// Functional annotation
			updateMessage("Performing functional annotation...");
			updateProgress(0.85, 1.0);
			if (!parentTask.performFunctionalAnnotation(proteins, proteinIntensityMap, intensityTitle, taskName)) {
				return false;
			}

			// Taxonomy analysis
			updateMessage("Performing taxonomy analysis...");
			updateProgress(0.90, 1.0);
			if (!parentTask.performTaxonomyAnalysis(genomeIntensityMap, genomeScoreMap, psmCountMap, totalPepCountMap,
					razorPepCountMap, intensityTitle, taskName)) {
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
	}

}
