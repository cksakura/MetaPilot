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
 * Builds the prediction dataset TSV for one run (one raw file) in the Deep Comprehensive workflow.
 * <p>
 * Candidate peptides are collected from the MAG DB predicted folder, keeping only rows whose
 * protein id's genome is in the selected genome set. Profile features (tax_*, fun_*) are computed
 * from P_0 built from precursors for this run only (filtered by runString matching the raw file).
 * Output has the same feature columns as the training TSV but no label_detect/label_intensity.
 */
public class DeepComprehensivePredictDatasetBuilder {

	private static final Logger LOG = LogManager.getLogger(DeepComprehensivePredictDatasetBuilder.class);

	private static final DecimalFormat FMT_5 = new DecimalFormat("0.#####", DecimalFormatSymbols.getInstance(java.util.Locale.US));

	private static String formatDouble(double v) {
		if (Double.isNaN(v) || Double.isInfinite(v)) return "";
		return FMT_5.format(v);
	}

	private static final String[] TAX_LEVEL_NAMES = { "domain", "phylum", "clazz", "order", "family", "genus", "species" };

	private static String extractGenomeId(String proteinId) {
		if (proteinId == null) return "";
		int u = proteinId.indexOf('_');
		return u > 0 ? proteinId.substring(0, u) : proteinId;
	}

	/** Candidate peptide: max detectability and the set of protein ids it appears in (from predicted files). */
	public static final class CandidatePeptide {
		public double maxDetectability = Double.NaN;
		public final Set<String> proteins = new HashSet<>();
	}

	/**
	 * Load all candidate peptides from the predicted folder, keeping only peptides from selected genomes.
	 * Expected TSV format: Protein id, Peptide sequence, Peptide detectability (header names may vary).
	 *
	 * @param predictedFolder MAG DB predicted folder (e.g. getPredictedPepFile())
	 * @param selectedGenomes  genome IDs to include (protein id before first '_')
	 * @return peptide -> CandidatePeptide (max detectability, set of protein ids)
	 */
	public static Map<String, CandidatePeptide> loadCandidatePeptidesFromPredicted(
			File predictedFolder,
			Set<String> selectedGenomes) throws IOException {
		Map<String, CandidatePeptide> out = new HashMap<>();
		if (predictedFolder == null || !predictedFolder.isDirectory() || selectedGenomes == null || selectedGenomes.isEmpty()) {
			return out;
		}
		File[] files = predictedFolder.listFiles();
		if (files == null) return out;
		for (File f : files) {
			if (!f.isFile() || !f.getName().endsWith(".predicted.tsv")) continue;
			try (BufferedReader r = new BufferedReader(new FileReader(f))) {
				String header = r.readLine();
				if (header == null) continue;
				String[] h = header.split("\t");
				int proIdx = -1, seqIdx = -1, detIdx = -1;
				for (int i = 0; i < h.length; i++) {
					String t = h[i].trim();
					if ("Protein id".equalsIgnoreCase(t) || "Protein id".equals(t)) proIdx = i;
					else if ("Peptide sequence".equals(t)) seqIdx = i;
					else if ("Peptide detectability".equals(t)) detIdx = i;
				}
				if (seqIdx < 0) continue;
				if (proIdx < 0) proIdx = 0;
				String line;
				while ((line = r.readLine()) != null) {
					String[] cs = line.split("\t");
					if (proIdx >= cs.length || seqIdx >= cs.length) continue;
					String proteinId = cs[proIdx].trim();
					String genome = extractGenomeId(proteinId);
					if (!selectedGenomes.contains(genome)) continue;
					String pep = cs[seqIdx].trim();
					if (pep.isEmpty()) continue;
					CandidatePeptide c = out.computeIfAbsent(pep, k -> new CandidatePeptide());
					c.proteins.add(proteinId);
					if (detIdx >= 0 && detIdx < cs.length) {
						try {
							double d = Double.parseDouble(cs[detIdx].trim());
							if (Double.isNaN(c.maxDetectability) || d > c.maxDetectability) c.maxDetectability = d;
						} catch (NumberFormatException ignored) {}
					}
				}
			}
		}
		LOG.info("Loaded {} candidate peptides from predicted folder (selected genomes: {})", out.size(), selectedGenomes.size());
		return out;
	}

	private static final class PeptideAggregate {
		double sumNormalised;
		final Set<String> proteins = new HashSet<>();
	}

	private static Map<String, PeptideAggregate> buildPeptideAggregates(List<DiaNNPrecursor> precursors) {
		Map<String, PeptideAggregate> map = new HashMap<>();
		for (DiaNNPrecursor p : precursors) {
			if (p == null) continue;
			String seq = p.getSeqString();
			if (seq == null || seq.isEmpty()) continue;
			PeptideAggregate a = map.computeIfAbsent(seq, k -> new PeptideAggregate());
			a.sumNormalised += p.getPreTrans();
			if (p.getProGroups() != null) {
				for (String pro : p.getProGroups()) {
					if (pro != null && !pro.isEmpty()) a.proteins.add(pro);
				}
			}
		}
		return map;
	}

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

	private static void addTax(ProfileP0 p0, TaxaInfo ti, double intensity) {
		String[] vals = { ti.domain, ti.phylum, ti.clazz, ti.order, ti.family, ti.genus, ti.species };
		for (int i = 0; i < 7; i++) {
			String v = (vals[i] != null && !vals[i].isEmpty()) ? vals[i] : "NA";
			Map<String, Double> m = p0.taxLevelToNodeIntensity.get(TAX_LEVEL_NAMES[i]);
			m.put(v, m.getOrDefault(v, 0.0) + intensity);
		}
	}

	private static ProfileP0 buildProfileP0(Map<String, PeptideAggregate> agg, Map<String, TaxaInfo> taxaMap,
			Map<String, Map<String, String[]>> genomeProteinFuncLevels) {
		ProfileP0 p0 = new ProfileP0();
		for (Map.Entry<String, PeptideAggregate> e : agg.entrySet()) {
			double intensity = e.getValue().sumNormalised;
			for (String pro : e.getValue().proteins) {
				String genome = extractGenomeId(pro);
				TaxaInfo ti = taxaMap.get(genome);
				if (ti != null) addTax(p0, ti, intensity);
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

	private static Map<String, Map<String, String[]>> loadGenomeProteinFunctionLevels(File funcFolder, Set<String> genomes) throws IOException {
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
					if ("eggNOG_OGs".equals(h[i].trim())) { eggIdx = i; break; }
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

	/**
	 * Build the prediction dataset for one run: P_0 from run precursors, candidate peptides from selected genomes.
	 *
	 * @param runPrecursors   precursors for this run only (runString matched to raw file)
	 * @param selectedGenomes genome IDs to include when loading candidates from predicted folder
	 * @param magDb           MAG DB (currentFile, getPredictedPepFile(), getUsedVersion().getFuncName())
	 * @param outputTsv       output TSV (same feature columns as training, no labels)
	 * @return number of rows written
	 */
	public int build(
			DiaNNPrecursor[] runPrecursors,
			Set<String> selectedGenomes,
			MagDbItem magDb,
			File outputTsv) throws IOException {
		if (magDb == null) {
			throw new IOException("MagDbItem is required");
		}
		File currentFile = magDb.getCurrentFile();
		if (currentFile == null || !currentFile.isDirectory()) {
			throw new IOException("MAG DB current file is not a directory");
		}
		File predictedFolder = magDb.getPredictedPepFile();
		Map<String, CandidatePeptide> candidates = loadCandidatePeptidesFromPredicted(
				predictedFolder != null && predictedFolder.isDirectory() ? predictedFolder : null,
				selectedGenomes != null ? selectedGenomes : new HashSet<>());
		if (candidates.isEmpty()) {
			LOG.warn("No candidate peptides for selected genomes");
			return 0;
		}

		// P_0 from this run's precursors only
		List<DiaNNPrecursor> list = new ArrayList<>();
		if (runPrecursors != null) {
			for (DiaNNPrecursor p : runPrecursors) {
				if (p != null && !p.getIsDecoy()) list.add(p);
			}
		}
		Map<String, PeptideAggregate> agg = buildPeptideAggregates(list);
		Set<String> allGenomes = new HashSet<>(selectedGenomes != null ? selectedGenomes : new HashSet<>());
		for (PeptideAggregate a : agg.values()) {
			for (String pro : a.proteins) allGenomes.add(extractGenomeId(pro));
		}

		File genomesTaxaTsv = new File(currentFile, "genomes-taxa.tsv");
		Map<String, TaxaInfo> taxaMap = new HashMap<>();
		if (genomesTaxaTsv.isFile()) {
			taxaMap = MarkerGenomeSelector.loadGenomeTaxaMap(genomesTaxaTsv);
		}
		String funcName = magDb.getUsedVersion() != null ? magDb.getUsedVersion().getFuncName() : "eggNOG";
		File funcFolder = new File(currentFile, funcName);
		Map<String, Map<String, String[]>> genomeProteinFuncLevels = loadGenomeProteinFunctionLevels(funcFolder, allGenomes);

		ProfileP0 p0 = buildProfileP0(agg, taxaMap, genomeProteinFuncLevels);
		p0.computeTop2Sums();

		File parent = outputTsv.getParentFile();
		if (parent != null && !parent.exists()) parent.mkdirs();

		return writePredictTsv(candidates, taxaMap, genomeProteinFuncLevels, p0, outputTsv);
	}

	private int writePredictTsv(
			Map<String, CandidatePeptide> candidates,
			Map<String, TaxaInfo> taxaMap,
			Map<String, Map<String, String[]>> genomeProteinFuncLevels,
			ProfileP0 p0,
			File outputTsv) throws IOException {
		try (PrintWriter w = new PrintWriter(outputTsv, "UTF-8")) {
			writePredictHeader(w);
			int rows = 0;
			for (Map.Entry<String, CandidatePeptide> e : candidates.entrySet()) {
				String peptide = e.getKey();
				CandidatePeptide c = e.getValue();
				int length = peptide.length();
				int missCleav = PeptideSequenceFeatures.missCleavageFromSequence(peptide);
				double gravy = PeptideSequenceFeatures.gravy(peptide);
				int basic = PeptideSequenceFeatures.basicCount(peptide);
				int acidic = PeptideSequenceFeatures.acidicCount(peptide);
				int aromatic = PeptideSequenceFeatures.aromaticCount(peptide);
				int netCharge = PeptideSequenceFeatures.netCharge(peptide);
				double detect = c.maxDetectability;

				Set<String> genomes = new HashSet<>();
				Set<String> species = new HashSet<>();
				Set<String> genera = new HashSet<>();
				for (String pro : c.proteins) {
					String g = extractGenomeId(pro);
					genomes.add(g);
					TaxaInfo ti = taxaMap.get(g);
					if (ti != null) {
						if (ti.species != null && !ti.species.isEmpty()) species.add(ti.species);
						if (ti.genus != null && !ti.genus.isEmpty()) genera.add(ti.genus);
					}
				}
				int nGenomes = genomes.size();
				int nProteins = c.proteins.size();
				int nSpecies = species.size();
				int nGenus = genera.size();
				int uniqueSpecies = (nSpecies <= 1) ? 1 : 0;

				double[] taxSum = new double[7];
				double[] taxMax = new double[7];
				for (int i = 0; i < 7; i++) {
					String key = TAX_LEVEL_NAMES[i];
					Map<String, Double> nodeInt = p0.taxLevelToNodeIntensity.get(key);
					double sum = 0, max = 0;
					for (String pro : c.proteins) {
						TaxaInfo ti = taxaMap.get(extractGenomeId(pro));
						if (ti == null) continue;
						String node = getTaxNode(ti, i);
						Double v = nodeInt.get(node);
						if (v != null) { sum += v; if (v > max) max = v; }
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
					for (String pro : c.proteins) {
						String genome = extractGenomeId(pro);
						Map<String, String[]> proLevels = genomeProteinFuncLevels.get(genome);
						if (proLevels == null) continue;
						String[] levels = proLevels.get(pro);
						if (levels == null || i >= levels.length) continue;
						String node = levels[i] != null ? levels[i] : "NA";
						nodes.add(node);
						Double v = nodeInt.get(node);
						if (v != null) { sum += v; if (v > max) max = v; }
					}
					funSum[i] = sum;
					funMax[i] = max;
					funNnodes[i] = nodes.size();
				}

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
				w.println();
				rows++;
			}
			return rows;
		}
	}

	private static void writePredictHeader(PrintWriter w) {
		w.print("peptide\tlength\tmiss_cleavage\tgravy\tbasic_count\tacidic_count\taromatic_count\tnet_charge\tdetectability");
		w.print("\tn_genomes\tn_proteins\tn_species\tn_genus\tunique_species");
		for (int i = 0; i < 7; i++) {
			w.print("\ttax_sum_" + (i+1) + "\ttax_max_" + (i+1) + "\ttax_top2sum_" + (i+1));
		}
		for (int i = 0; i < 5; i++) {
			w.print("\tfun_sum_" + (i+1) + "\tfun_max_" + (i+1) + "\tfun_top2sum_" + (i+1) + "\tfun_n_nodes_" + (i+1));
		}
		w.println();
	}
}
