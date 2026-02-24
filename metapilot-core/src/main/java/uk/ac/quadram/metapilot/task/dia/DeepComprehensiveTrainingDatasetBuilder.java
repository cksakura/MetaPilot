package uk.ac.quadram.metapilot.task.dia;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNPrecursor;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;
import uk.ac.quadram.metapilot.task.mag.MarkerGenomeSelector;
import uk.ac.quadram.metapilot.task.mag.MarkerGenomeSelector.TaxaInfo;

/**
 * Builds the initial (seed) training dataset for the Deep Comprehensive workflow:
 * one row per high-confidence peptide with static features, dynamic profile features (P_0),
 * and labels for the two-head model (detectability=1, intensity=log1p(sum Precursor.Normalised)).
 * <p>
 * Uses: MAG DB predicted folder (detectability), genomes-taxa.tsv (7 tax levels),
 * eggNOG TSVs (5 function levels from eggNOG_OGs column), and DiaNNPrecursor (getPreTrans = Precursor.Normalised).
 * Intensity label is log1p(sum Normalised) to reduce skew and match the regression head.
 */
public class DeepComprehensiveTrainingDatasetBuilder {

	private static final Logger LOG = LogManager.getLogger(DeepComprehensiveTrainingDatasetBuilder.class);

	/** Format doubles with up to 5 decimal places, no scientific notation. */
	private static final DecimalFormat FMT_5 = new DecimalFormat("0.#####", DecimalFormatSymbols.getInstance(java.util.Locale.US));

	private static String formatDouble(double v) {
		if (Double.isNaN(v) || Double.isInfinite(v)) return "";
		return FMT_5.format(v);
	}

	/** Default q-value threshold for high-confidence precursors (allow 0.05 if seed sparse). */
	public static final double DEFAULT_Q_THRESHOLD = 0.01;

	private final double qThreshold;
	private final boolean excludeDecoys;

	public DeepComprehensiveTrainingDatasetBuilder() {
		this(DEFAULT_Q_THRESHOLD, true);
	}

	public DeepComprehensiveTrainingDatasetBuilder(double qThreshold, boolean excludeDecoys) {
		this.qThreshold = qThreshold;
		this.excludeDecoys = excludeDecoys;
	}

	/**
	 * Builds the training TSV from high-confidence precursors and MAG DB.
	 *
	 * @param precursors  all precursors from seed search (will filter by q and decoy)
	 * @param magDb       MAG database item (currentFile, getPredictedPepFile, getUsedVersion().getFuncName())
	 * @param outputTsv   output file (TSV, one row per unique peptide)
	 * @return number of rows written
	 */
	public int build(DiaNNPrecursor[] precursors, MagDbItem magDb, File outputTsv) throws IOException {
		if (precursors == null || precursors.length == 0) {
			LOG.warn("No precursors provided");
			return 0;
		}
		if (magDb == null) {
			throw new IOException("MagDbItem is required");
		}
		File currentFile = magDb.getCurrentFile();
		if (currentFile == null || !currentFile.isDirectory()) {
			throw new IOException("MAG DB current file is not a directory");
		}

		// Filter high-confidence, non-decoy
		List<DiaNNPrecursor> highConf = new ArrayList<>();
		for (DiaNNPrecursor p : precursors) {
			if (p == null) continue;
			if (excludeDecoys && p.getIsDecoy()) continue;
			double q = p.getGlobalQvalue() >= 0 ? p.getGlobalQvalue() : p.getQvalue();
			if (q > qThreshold) continue;
			highConf.add(p);
		}
		if (highConf.isEmpty()) {
			LOG.warn("No high-confidence precursors after filtering (q<={})", qThreshold);
			return 0;
		}

		// Peptide aggregates: peptide -> { sum Normalised, max missCleavage, set of proteins }
		Map<String, PeptideAggregate> agg = buildPeptideAggregates(highConf);
		Set<String> allGenomes = new HashSet<>();
		for (PeptideAggregate a : agg.values()) {
			for (String pro : a.proteins) {
				allGenomes.add(extractGenomeId(pro));
			}
		}

		// Load detectability from predicted TSVs
		File predictedFolder = magDb.getPredictedPepFile();
		Map<String, Double> detectabilityMap = loadDetectabilityMap(predictedFolder != null && predictedFolder.isDirectory() ? predictedFolder : null, agg.keySet());

		// Load genome -> taxa (7 levels)
		File genomesTaxaTsv = new File(currentFile, "genomes-taxa.tsv");
		Map<String, TaxaInfo> taxaMap = new HashMap<>();
		if (genomesTaxaTsv.isFile()) {
			taxaMap = MarkerGenomeSelector.loadGenomeTaxaMap(genomesTaxaTsv);
		} else {
			LOG.warn("genomes-taxa.tsv not found: {}", genomesTaxaTsv.getAbsolutePath());
		}

		// Load genome -> (protein -> function 5 levels) from eggNOG TSVs
		String funcName = magDb.getUsedVersion() != null ? magDb.getUsedVersion().getFuncName() : "eggNOG";
		File funcFolder = new File(currentFile, funcName);
		Map<String, Map<String, String[]>> genomeProteinFuncLevels = loadGenomeProteinFunctionLevels(funcFolder, allGenomes, agg);

		// Build P_0: 7 tax level -> (node -> intensity), 5 function level -> (node -> intensity)
		ProfileP0 p0 = buildProfileP0(agg, taxaMap, genomeProteinFuncLevels);
		p0.computeTop2Sums();

		// Ensure output parent (e.g. result/training/) exists
		File parent = outputTsv.getParentFile();
		if (parent != null && !parent.exists()) {
			parent.mkdirs();
		}

		// Write TSV
		int rows = writeTsv(agg, detectabilityMap, taxaMap, genomeProteinFuncLevels, p0, outputTsv);
		LOG.info("Wrote {} training rows to {}", rows, outputTsv.getAbsolutePath());
		return rows;
	}

	private static String extractGenomeId(String proteinId) {
		if (proteinId == null) return "";
		int u = proteinId.indexOf('_');
		return u > 0 ? proteinId.substring(0, u) : proteinId;
	}

	private static final class PeptideAggregate {
		double sumNormalised;
		int maxMissCleavage;
		Set<String> proteins = new HashSet<>();
	}

	private Map<String, PeptideAggregate> buildPeptideAggregates(List<DiaNNPrecursor> precursors) {
		Map<String, PeptideAggregate> map = new HashMap<>();
		for (DiaNNPrecursor p : precursors) {
			String seq = p.getSeqString();
			if (seq == null || seq.isEmpty()) continue;
			PeptideAggregate a = map.computeIfAbsent(seq, k -> new PeptideAggregate());
			a.sumNormalised += p.getPreTrans();
			if (p.getMissCleavage() > a.maxMissCleavage) a.maxMissCleavage = p.getMissCleavage();
			if (p.getProGroups() != null) {
				for (String pro : p.getProGroups()) {
					if (pro != null && !pro.isEmpty()) a.proteins.add(pro);
				}
			}
		}
		return map;
	}

	private Map<String, Double> loadDetectabilityMap(File predictedFolder, Set<String> targetPeptides) throws IOException {
		Map<String, Double> map = new HashMap<>();
		if (predictedFolder == null || !predictedFolder.isDirectory()) return map;
		File[] files = predictedFolder.listFiles();
		if (files == null) return map;
		for (File f : files) {
			if (!f.isFile() || !f.getName().endsWith(".predicted.tsv")) continue;
			try (BufferedReader r = new BufferedReader(new FileReader(f))) {
				String header = r.readLine();
				if (header == null) continue;
				String[] h = header.split("\t");
				int seqIdx = -1, detIdx = -1;
				for (int i = 0; i < h.length; i++) {
					String t = h[i].trim();
					if ("Peptide sequence".equals(t)) seqIdx = i;
					else if ("Peptide detectability".equals(t)) detIdx = i;
				}
				if (seqIdx < 0) continue;
				String line;
				while ((line = r.readLine()) != null) {
					String[] cs = line.split("\t");
					if (seqIdx >= cs.length) continue;
					String pep = cs[seqIdx].trim();
					if (!targetPeptides.contains(pep)) continue;
					if (detIdx >= 0 && detIdx < cs.length) {
						try {
							double d = Double.parseDouble(cs[detIdx].trim());
							Double existing = map.get(pep);
							if (existing == null || d > existing) map.put(pep, d);
						} catch (NumberFormatException ignored) {}
					}
				}
			}
		}
		return map;
	}

	/** Parse eggNOG_OGs e.g. "COG1197@1|root,COG1197@2|Bacteria,..." into level 1..5 names. */
	private static String[] parseEggNogOGsFiveLevels(String eggNogOGs) {
		String[] levels = new String[5];
		if (eggNogOGs == null || eggNogOGs.isEmpty()) return levels;
		String[] parts = eggNogOGs.split(",");
		for (int i = 0; i < 5 && i < parts.length; i++) {
			String p = parts[i].trim();
			int pipe = p.indexOf('|');
			levels[i] = pipe >= 0 ? p.substring(pipe + 1).trim() : p;
			if (levels[i].isEmpty()) levels[i] = "NA";
		}
		for (int i = 0; i < 5; i++) {
			if (levels[i] == null) levels[i] = "NA";
		}
		return levels;
	}

	private Map<String, Map<String, String[]>> loadGenomeProteinFunctionLevels(File funcFolder, Set<String> genomes,
			Map<String, PeptideAggregate> agg) throws IOException {
		Map<String, Map<String, String[]>> out = new HashMap<>();
		if (funcFolder == null || !funcFolder.isDirectory()) return out;
		String funcName = funcFolder.getName();
		for (String genome : genomes) {
			File f = new File(funcFolder, genome + "_" + funcName + ".tsv");
			if (!f.isFile()) continue;
			Map<String, String[]> proToLevels = new HashMap<>();
			try (BufferedReader r = new BufferedReader(new FileReader(f))) {
				String header = r.readLine();
				if (header == null) continue;
				String[] h = header.split("\t");
				int proIdx = 0;
				int eggIdx = -1;
				for (int i = 0; i < h.length; i++) {
					if ("eggNOG_OGs".equals(h[i].trim())) {
						eggIdx = i;
						break;
					}
				}
				if (eggIdx < 0) continue;
				String line;
				while ((line = r.readLine()) != null) {
					String[] cs = line.split("\t");
					if (cs.length <= Math.max(proIdx, eggIdx)) continue;
					String pro = cs[proIdx].trim();
					String og = eggIdx < cs.length ? cs[eggIdx] : "";
					proToLevels.put(pro, parseEggNogOGsFiveLevels(og));
				}
			}
			out.put(genome, proToLevels);
		}
		return out;
	}

	private static final String[] TAX_LEVEL_NAMES = { "domain", "phylum", "clazz", "order", "family", "genus", "species" };

	private static class ProfileP0 {
		final Map<String, Map<String, Double>> taxLevelToNodeIntensity = new HashMap<>();
		final Map<String, Map<String, Double>> funLevelToNodeIntensity = new HashMap<>();
		final double[] taxTop2Sum = new double[7];
		final double[] funTop2Sum = new double[5];

		ProfileP0() {
			for (int i = 0; i < 7; i++) taxLevelToNodeIntensity.put(TAX_LEVEL_NAMES[i], new TreeMap<>());
			for (int i = 0; i < 5; i++) funLevelToNodeIntensity.put("F" + (i + 1), new TreeMap<>());
		}

		void computeTop2Sums() {
			for (int i = 0; i < 7; i++) {
				List<Double> v = new ArrayList<>(taxLevelToNodeIntensity.get(TAX_LEVEL_NAMES[i]).values());
				v.sort((a, b) -> Double.compare(b, a));
				taxTop2Sum[i] = v.size() >= 2 ? v.get(0) + v.get(1) : (v.size() == 1 ? v.get(0) : 0);
			}
			for (int i = 0; i < 5; i++) {
				List<Double> v = new ArrayList<>(funLevelToNodeIntensity.get("F" + (i + 1)).values());
				v.sort((a, b) -> Double.compare(b, a));
				funTop2Sum[i] = v.size() >= 2 ? v.get(0) + v.get(1) : (v.size() == 1 ? v.get(0) : 0);
			}
		}
	}

	private ProfileP0 buildProfileP0(Map<String, PeptideAggregate> agg, Map<String, TaxaInfo> taxaMap,
			Map<String, Map<String, String[]>> genomeProteinFuncLevels) {
		ProfileP0 p0 = new ProfileP0();
		for (Map.Entry<String, PeptideAggregate> e : agg.entrySet()) {
			double intensity = e.getValue().sumNormalised;
			Set<String> genomes = new HashSet<>();
			for (String pro : e.getValue().proteins) {
				genomes.add(extractGenomeId(pro));
			}
			for (String genome : genomes) {
				TaxaInfo ti = taxaMap.get(genome);
				if (ti != null) {
					addTax(p0, ti, intensity);
				}
			}
			for (String pro : e.getValue().proteins) {
				String genome = extractGenomeId(pro);
				Map<String, String[]> proLevels = genomeProteinFuncLevels.get(genome);
				if (proLevels != null) {
					String[] levels = proLevels.get(pro);
					if (levels != null) {
						for (int i = 0; i < 5; i++) {
							String node = levels[i] != null ? levels[i] : "NA";
							Map<String, Double> m = p0.funLevelToNodeIntensity.get("F" + (i + 1));
							m.put(node, m.getOrDefault(node, 0.0) + intensity);
						}
					}
				}
			}
		}
		return p0;
	}

	private static void addTax(ProfileP0 p0, TaxaInfo ti, double intensity) {
		String[] vals = { ti.domain, ti.phylum, ti.clazz, ti.order, ti.family, ti.genus, ti.species };
		for (int i = 0; i < 7; i++) {
			String v = (vals[i] != null && !vals[i].isEmpty()) ? vals[i] : "NA";
			Map<String, Double> m = p0.taxLevelToNodeIntensity.get(TAX_LEVEL_NAMES[i]);
			m.put(v, m.getOrDefault(v, 0.0) + intensity);
		}
	}

	private int writeTsv(Map<String, PeptideAggregate> agg, Map<String, Double> detectabilityMap,
			Map<String, TaxaInfo> taxaMap, Map<String, Map<String, String[]>> genomeProteinFuncLevels,
			ProfileP0 p0, File outputTsv) throws IOException {
		try (PrintWriter w = new PrintWriter(outputTsv, "UTF-8")) {
			writeHeader(w);
			int rows = 0;
			for (Map.Entry<String, PeptideAggregate> e : agg.entrySet()) {
				String peptide = e.getKey();
				PeptideAggregate a = e.getValue();
				double sumNorm = a.sumNormalised;
				int missCleav = a.maxMissCleavage;
				int length = peptide.length();
				double gravy = PeptideSequenceFeatures.gravy(peptide);
				int basic = PeptideSequenceFeatures.basicCount(peptide);
				int acidic = PeptideSequenceFeatures.acidicCount(peptide);
				int aromatic = PeptideSequenceFeatures.aromaticCount(peptide);
				int netCharge = PeptideSequenceFeatures.netCharge(peptide);
				double detect = detectabilityMap.getOrDefault(peptide, Double.NaN);

				Set<String> genomes = new HashSet<>();
				Set<String> species = new HashSet<>();
				Set<String> genera = new HashSet<>();
				for (String pro : a.proteins) {
					String g = extractGenomeId(pro);
					genomes.add(g);
					TaxaInfo ti = taxaMap.get(g);
					if (ti != null) {
						if (ti.species != null && !ti.species.isEmpty()) species.add(ti.species);
						if (ti.genus != null && !ti.genus.isEmpty()) genera.add(ti.genus);
					}
				}
				int nGenomes = genomes.size();
				int nProteins = a.proteins.size();
				int nSpecies = species.size();
				int nGenus = genera.size();
				int uniqueSpecies = (nSpecies <= 1) ? 1 : 0;

				double[] taxSum = new double[7];
				double[] taxMax = new double[7];
				for (int i = 0; i < 7; i++) {
					String key = TAX_LEVEL_NAMES[i];
					Map<String, Double> nodeInt = p0.taxLevelToNodeIntensity.get(key);
					double sum = 0, max = 0;
					for (String pro : a.proteins) {
						TaxaInfo ti = taxaMap.get(extractGenomeId(pro));
						if (ti == null) continue;
						String node = getTaxNode(ti, i);
						Double v = nodeInt.get(node);
						if (v != null) {
							sum += v;
							if (v > max) max = v;
						}
					}
					taxSum[i] = sum;
					taxMax[i] = max;
				}

				double[] funSum = new double[5];
				double[] funMax = new double[5];
				int[] funNnodes = new int[5];
				for (int i = 0; i < 5; i++) {
					String key = "F" + (i + 1);
					Map<String, Double> nodeInt = p0.funLevelToNodeIntensity.get(key);
					Set<String> nodes = new HashSet<>();
					double sum = 0, max = 0;
					for (String pro : a.proteins) {
						String genome = extractGenomeId(pro);
						Map<String, String[]> proLevels = genomeProteinFuncLevels.get(genome);
						if (proLevels == null) continue;
						String[] levels = proLevels.get(pro);
						if (levels == null || i >= levels.length) continue;
						String node = levels[i] != null ? levels[i] : "NA";
						nodes.add(node);
						Double v = nodeInt.get(node);
						if (v != null) {
							sum += v;
							if (v > max) max = v;
						}
					}
					funSum[i] = sum;
					funMax[i] = max;
					funNnodes[i] = nodes.size();
				}

				double labelIntensity = Math.log1p(sumNorm);

				w.print(peptide);
				w.print("\t"); w.print(length);
				w.print("\t"); w.print(missCleav);
				w.print("\t"); w.print(formatDouble(gravy));
				w.print("\t"); w.print(basic);
				w.print("\t"); w.print(acidic);
				w.print("\t"); w.print(aromatic);
				w.print("\t"); w.print(netCharge);
				w.print("\t"); w.print(Double.isNaN(detect) ? "" : formatDouble(detect));
				w.print("\t"); w.print(nGenomes);
				w.print("\t"); w.print(nProteins);
				w.print("\t"); w.print(nSpecies);
				w.print("\t"); w.print(nGenus);
				w.print("\t"); w.print(uniqueSpecies);
				for (int i = 0; i < 7; i++) {
					w.print("\t"); w.print(formatDouble(taxSum[i]));
					w.print("\t"); w.print(formatDouble(taxMax[i]));
					w.print("\t"); w.print(formatDouble(p0.taxTop2Sum[i]));
				}
				for (int i = 0; i < 5; i++) {
					w.print("\t"); w.print(formatDouble(funSum[i]));
					w.print("\t"); w.print(formatDouble(funMax[i]));
					w.print("\t"); w.print(formatDouble(p0.funTop2Sum[i]));
					w.print("\t"); w.print(funNnodes[i]);
				}
				w.print("\t"); w.print(1);
				w.print("\t"); w.print(formatDouble(labelIntensity));
				w.println();
				rows++;
			}
			return rows;
		}
	}

	private static String getTaxNode(TaxaInfo ti, int level) {
		String v = null;
		switch (level) {
			case 0: v = ti.domain; break;
			case 1: v = ti.phylum; break;
			case 2: v = ti.clazz; break;
			case 3: v = ti.order; break;
			case 4: v = ti.family; break;
			case 5: v = ti.genus; break;
			case 6: v = ti.species; break;
			default: v = "NA";
		}
		return (v != null && !v.isEmpty()) ? v : "NA";
	}

	private void writeHeader(PrintWriter w) {
		w.print("peptide\tlength\tmiss_cleavage\tgravy\tbasic_count\tacidic_count\taromatic_count\tnet_charge\tdetectability");
		w.print("\tn_genomes\tn_proteins\tn_species\tn_genus\tunique_species");
		for (int i = 0; i < 7; i++) {
			w.print("\ttax_sum_" + (i+1) + "\ttax_max_" + (i+1) + "\ttax_top2sum_" + (i+1));
		}
		for (int i = 0; i < 5; i++) {
			w.print("\tfun_sum_" + (i+1) + "\tfun_max_" + (i+1) + "\tfun_top2sum_" + (i+1) + "\tfun_n_nodes_" + (i+1));
		}
		w.print("\tlabel_detect\tlabel_intensity");
		w.println();
	}
}
