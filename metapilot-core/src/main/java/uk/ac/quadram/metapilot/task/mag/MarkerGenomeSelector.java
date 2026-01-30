package uk.ac.quadram.metapilot.task.mag;

import java.io.*;
import java.util.*;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import uk.ac.quadram.metapilot.dbSearch.pfind.PFindPSM;
import uk.ac.quadram.metapilot.dbSearch.pfind.PFindPsmReader;
import uk.ac.quadram.metapilot.task.io.pep.MetaPSM;

public class MarkerGenomeSelector {

	private static final Logger LOGGER = LogManager.getLogger(MarkerGenomeSelector.class);

	public enum RepresentativeLevel {
		GENUS, FAMILY, ORDER, CLASS, PHYLUM
	}

	private RepresentativeLevel representativeLevel = RepresentativeLevel.GENUS;

	// Regex patterns for the 4-category markers (case-insensitive)
	private static final Pattern ANCHOR_PATTERN = Pattern.compile("(?i)\\bribosomal\\s+protein\\b");
	private static final Pattern HIGH_RES_PATTERN = Pattern.compile(
			"(?i)\\b(gyrase|gyr[ab])\\b|\\brna\\s+polymerase\\b|\\bdna\\s+polymerase\\b|\\brpo[abc]\\b|\\bdna[eq]\\b");
	private static final Pattern GROWTH_PATTERN = Pattern.compile("(?i)\\belongation\\s+factor\\b|\\bef[-\\s]?(tu|g|ts|p|4|1a|2)\\b");
	private static final Pattern STABILITY_PATTERN = Pattern.compile(
			"(?i)(aminoacyl[-\\s]?trna|aa[-\\s]?trna|trna)[-\\s]*[^a-zA-Z0-9]{0,2}\\s*(ligase|synthetase)\\b");

	// ----------------------------
	// Taxonomy handling
	// ----------------------------
	public static class TaxaInfo {
		public String domain = "", phylum = "", clazz = "", order = "", family = "", genus = "", species = "";
		public String clusterKey = "NA"; // species if present else genus else up ranks
	}

	public static HashMap<String, TaxaInfo> loadGenomeTaxaMap(File taxaTsv) throws IOException {
		HashMap<String, TaxaInfo> map = new HashMap<>();
		try (BufferedReader br = new BufferedReader(new FileReader(taxaTsv))) {
			String line;
			while ((line = br.readLine()) != null) {
				line = line.trim();
				if (line.isEmpty())
					continue;
				String[] cs = line.split("\t");
				if (cs.length < 2)
					continue;

				String genomeId = cs[0].trim();
				String tax = cs[1].trim();

				TaxaInfo ti = new TaxaInfo();

				String[] parts = tax.split(";");
				for (String p : parts) {
					p = p.trim();
					if (p.startsWith("d__"))
						ti.domain = p;
					else if (p.startsWith("p__"))
						ti.phylum = p;
					else if (p.startsWith("c__"))
						ti.clazz = p;
					else if (p.startsWith("o__"))
						ti.order = p;
					else if (p.startsWith("f__"))
						ti.family = p;
					else if (p.startsWith("g__"))
						ti.genus = p;
					else if (p.startsWith("s__"))
						ti.species = p;
				}

				String speciesName = (ti.species != null && ti.species.length() > 3) ? ti.species.substring(3).trim()
						: "";
				String genusName = (ti.genus != null && ti.genus.length() > 3) ? ti.genus.substring(3).trim() : "";

				if (!speciesName.isEmpty()) {
					ti.clusterKey = ti.species; // s__...
				} else if (!genusName.isEmpty()) {
					ti.clusterKey = ti.genus; // g__...
				} else if (ti.family != null && ti.family.length() > 3 && !ti.family.substring(3).trim().isEmpty()) {
					ti.clusterKey = ti.family;
				} else if (ti.order != null && ti.order.length() > 3 && !ti.order.substring(3).trim().isEmpty()) {
					ti.clusterKey = ti.order;
				} else if (ti.clazz != null && ti.clazz.length() > 3 && !ti.clazz.substring(3).trim().isEmpty()) {
					ti.clusterKey = ti.clazz;
				} else if (ti.phylum != null && ti.phylum.length() > 3 && !ti.phylum.substring(3).trim().isEmpty()) {
					ti.clusterKey = ti.phylum;
				} else if (ti.domain != null && !ti.domain.isEmpty()) {
					ti.clusterKey = ti.domain;
				} else {
					ti.clusterKey = "NA";
				}

				map.put(genomeId, ti);
			}
		}
		return map;
	}

	// ----------------------------
	// Evidence row
	// ----------------------------
	static class GenomeRow {
		String genome;
		String clusterKey; // species if present else genus (kept for reporting/debug)
		String genusKey; // prefer g__..., fallback to clusterKey
		String familyKey; // prefer f__..., fallback to clusterKey
		String orderKey;
		String classKey;
		String phylumKey;
		int totalPsm;
		int anchor, highRes, growth, housekeeping;
		int uAnchor, uHighRes, uGrowth, uHousekeeping;
		int totalMarkers, totalUnique;

		// sets for greedy selection
		Set<String> allMarkerPeptides = new HashSet<>();

		int nCat() {
			return (anchor > 0 ? 1 : 0) + (highRes > 0 ? 1 : 0) + (growth > 0 ? 1 : 0) + (housekeeping > 0 ? 1 : 0);
		}

		int nCatUnique() {
			return (uAnchor > 0 ? 1 : 0) + (uHighRes > 0 ? 1 : 0) + (uGrowth > 0 ? 1 : 0) + (uHousekeeping > 0 ? 1 : 0);
		}

		boolean hasAll4Cats() {
			return anchor > 0 && highRes > 0 && growth > 0 && housekeeping > 0;
		}

		// Strategy step 0 filter: must have Ribosomal AND Replication AND Elongation
		boolean passesHardFilter() {
			return anchor > 0 && highRes > 0 && growth > 0;
		}
	}

	// ----------------------------
	// Selection result
	// ----------------------------
	public static class SelectionResult {
		public List<String> ranked; // ranking used for genus rep step (unique-first)
		public LinkedHashSet<String> seedIncluded; // privileged + genus reps
		public LinkedHashSet<String> selected; // final selected list (seed + greedy)
		public Map<String, String> reason; // why included (PRIVILEGE / GENUS_REP / GREEDY)
		public Map<String, Integer> marginalGain; // gain when selected by greedy
		public int totalPeptides; // size of totalPeptides
	}

	// ----------------------------
	// Main pipeline: your peptideTest + strategy
	// ----------------------------

	public void peptideTestWithStrategy(PFindPSM[] psms, File markerFasta, File genomesTaxaTsv,
			String outputGenomeEvidenceTsv, String outputRankedTsv, String outputSelectedTsv, double coverageThreshold,
			int stepThreshold, int maximumSelected) {
		peptideTestWithStrategy((MetaPSM[]) psms, markerFasta, genomesTaxaTsv, outputGenomeEvidenceTsv,
				outputRankedTsv, outputSelectedTsv, coverageThreshold, stepThreshold, maximumSelected);
	}

	public void peptideTestWithStrategy(MetaPSM[] psms, File markerFasta, File genomesTaxaTsv,
			String outputGenomeEvidenceTsv, String outputRankedTsv, String outputSelectedTsv, double coverageThreshold, // e.g.
																														// 0.95
			int stepThreshold, // keep for compatibility; optional stop
			int maximumSelected) {

		LOGGER.info("MarkerGenomeSelector start: input=MetaPSM[] (count={})",
				psms == null ? 0 : psms.length);
		LOGGER.info("Representative level: {}", representativeLevel);
		LOGGER.info("Marker FASTA: {}", markerFasta != null ? markerFasta.getAbsolutePath() : "null");
		LOGGER.info("Genomes taxa TSV: {}", genomesTaxaTsv != null ? genomesTaxaTsv.getAbsolutePath() : "null");
		LOGGER.info("Outputs: evidence={}, ranked={}, selected={}",
				outputGenomeEvidenceTsv, outputRankedTsv, outputSelectedTsv);
		LOGGER.info("Params: coverageThreshold={}, stepThreshold={}, maximumSelected={}",
				coverageThreshold, stepThreshold, maximumSelected);

		Map<String, String> proteinToCategory = new HashMap<>();
		try {
			proteinToCategory = loadMarkerProteinCategoriesFromFasta(markerFasta);
		} catch (IOException e) {
			LOGGER.error("Failed to load marker FASTA categories from {}", markerFasta, e);
		}
		LOGGER.info("Protein categories loaded: {}", proteinToCategory.size());

		HashMap<String, TaxaInfo> taxaMap = new HashMap<>();
		try {
			taxaMap = loadGenomeTaxaMap(genomesTaxaTsv);
		} catch (Exception e) {
			LOGGER.error("Failed to load taxa TSV from {}", genomesTaxaTsv, e);
		}
		LOGGER.info("Taxa map loaded: {}", taxaMap.size());

		if (psms == null || psms.length == 0) {
			LOGGER.warn("No PSMs provided; nothing to process.");
			return;
		}

		// genome -> category -> peptides
		Map<String, Map<String, Set<String>>> genomeCategoryPeptides = new HashMap<>();

		// peptide -> genomes (only for marker peptides)
		Map<String, Set<String>> peptideGenomes = new HashMap<>();

		// peptide -> genome -> categories (to assign unique peptides to categories)
		Map<String, Map<String, Set<String>>> peptideGenomeCategories = new HashMap<>();

		// track unique peptides by genome/category
		Map<String, Map<String, Set<String>>> genomeCategoryUniquePeptides = new HashMap<>();

		// global peptide set (marker peptides only)
		HashSet<String> allMarkerPeptides = new HashSet<>();

		// genome -> PSM ids (marker-mapped PSMs)
		Map<String, Set<Integer>> genomePsmIds = new HashMap<>();
		int totalMarkerPsms = 0;
		int psmIndex = 0;

		// ----------------------------
		// Parse input
		// ----------------------------
		int markerProteinHits = 0;
		int markerPeptideHits = 0;

		for (MetaPSM psm : psms) {
			if (psm == null) {
				continue;
			}

			String peptide = psm.getSequence();
			if (peptide == null || peptide.isEmpty()) {
				continue;
			}

			String[] pros = psm.getProteins();
			if (pros == null || pros.length == 0) {
				continue;
			}

			int currentPsmId = psmIndex++;

			// collect marker mappings for this peptide
			boolean peptideHasAnyMarker = false;
			HashSet<String> genomesForPsm = new HashSet<>();

			for (String pro : pros) {
				if (pro == null) {
					continue;
				}
				String proteinId = pro.trim();
				if (proteinId.isEmpty()) {
					continue;
				}

				boolean isDecoy = proteinId.startsWith("REV_");
				String lookupId = isDecoy ? proteinId.substring(4) : proteinId;

				String category = proteinToCategory.get(lookupId);
				String proteinForGenome = lookupId;
				if (category == null) {
					int last = lookupId.lastIndexOf("_");
					if (last > 0) {
						String baseId = lookupId.substring(0, last);
						category = proteinToCategory.get(baseId);
						if (category != null) {
							proteinForGenome = baseId;
						}
					}
				}
				if (category == null) {
					continue;
				}

				String genomeId = extractGenomeId(proteinForGenome);
				if (genomeId == null || genomeId.isEmpty()) {
					continue;
				}
				
				// don't need decoy here, because decoy genome can't have 4 types of marker proteins, so it can't help the FDR control
				if (isDecoy) {
					genomeId = "REV_" + genomeId;
					continue;
				}

				peptideHasAnyMarker = true;
				markerProteinHits++;

				peptideGenomes.computeIfAbsent(peptide, k -> new HashSet<>()).add(genomeId);
				genomesForPsm.add(genomeId);

				Map<String, Set<String>> gmap = genomeCategoryPeptides.computeIfAbsent(genomeId,
						k -> new HashMap<>());
				gmap.computeIfAbsent(category, k -> new HashSet<>()).add(peptide);

				Map<String, Set<String>> gcat = peptideGenomeCategories.computeIfAbsent(peptide,
						k -> new HashMap<>());
				gcat.computeIfAbsent(genomeId, k -> new HashSet<>()).add(category);
			}

			// only store peptide if it maps to at least one marker protein
			if (peptideHasAnyMarker) {
				allMarkerPeptides.add(peptide);
				markerPeptideHits++;
				totalMarkerPsms++;
				for (String genomeId : genomesForPsm) {
					genomePsmIds.computeIfAbsent(genomeId, k -> new HashSet<>()).add(currentPsmId);
				}
			}
		}

		LOGGER.info("Marker peptides (unique): {}", allMarkerPeptides.size());
		LOGGER.info("Peptide->genomes map size: {}", peptideGenomes.size());
		LOGGER.info("Genome category peptides: {}", genomeCategoryPeptides.size());
		LOGGER.info("Marker protein hits: {}, marker peptide hits: {}", markerProteinHits, markerPeptideHits);

		// ----------------------------
		// Determine unique marker peptides (peptide maps to exactly 1 genome)
		// and assign to categories seen for that genome
		// ----------------------------
		for (Map.Entry<String, Set<String>> entry : peptideGenomes.entrySet()) {
			String peptide = entry.getKey();
			Set<String> genomes = entry.getValue();
			if (genomes.size() != 1)
				continue;

			String genomeId = genomes.iterator().next();
			Map<String, Set<String>> genomeCategories = peptideGenomeCategories.get(peptide);
			if (genomeCategories == null)
				continue;

			Set<String> categories = genomeCategories.get(genomeId);
			if (categories == null)
				continue;

			Map<String, Set<String>> genomeUniqueMap = genomeCategoryUniquePeptides.computeIfAbsent(genomeId,
					k -> new HashMap<>());

			for (String cat : categories) {
				genomeUniqueMap.computeIfAbsent(cat, k -> new HashSet<>()).add(peptide);
			}
		}

		// ----------------------------
		// Build GenomeRow list
		// ----------------------------
		List<String> categories = Arrays.asList("Ribosomal", "Replication", "Elongation", "Synthetase");

		List<GenomeRow> rows = new ArrayList<>();
		for (String genomeId : genomeCategoryPeptides.keySet()) {
			GenomeRow r = new GenomeRow();
			r.genome = genomeId;

			TaxaInfo ti = taxaMap.get(genomeId);
			r.clusterKey = (ti == null) ? "NA" : ti.clusterKey;
			if (ti != null && ti.genus != null && ti.genus.length() > 3 && !ti.genus.substring(3).trim().isEmpty()) {
				r.genusKey = normalizeGenusKey(ti.genus); // g__...
			} else if (ti != null && ti.clusterKey != null && !ti.clusterKey.isEmpty()) {
				r.genusKey = normalizeGenusKey(ti.clusterKey); // fallback up ranks
			} else {
				r.genusKey = "NA";
			}
			if (ti != null && ti.family != null && ti.family.length() > 3 && !ti.family.substring(3).trim().isEmpty()) {
				r.familyKey = normalizeFamilyKey(ti.family); // f__...
			} else if (ti != null && ti.clusterKey != null && !ti.clusterKey.isEmpty()) {
				r.familyKey = normalizeFamilyKey(ti.clusterKey); // fallback up ranks
			} else {
				r.familyKey = "NA";
			}
			if (ti != null && ti.order != null && ti.order.length() > 3 && !ti.order.substring(3).trim().isEmpty()) {
				r.orderKey = ti.order;
			} else {
				r.orderKey = "NA";
			}
			if (ti != null && ti.clazz != null && ti.clazz.length() > 3 && !ti.clazz.substring(3).trim().isEmpty()) {
				r.classKey = ti.clazz;
			} else {
				r.classKey = "NA";
			}
			if (ti != null && ti.phylum != null && ti.phylum.length() > 3 && !ti.phylum.substring(3).trim().isEmpty()) {
				r.phylumKey = ti.phylum;
			} else {
				r.phylumKey = "NA";
			}

			Map<String, Set<String>> catPeps = genomeCategoryPeptides.get(genomeId);
			Map<String, Set<String>> uniqPeps = genomeCategoryUniquePeptides.get(genomeId);

			Set<Integer> psmIds = genomePsmIds.get(genomeId);
			r.totalPsm = psmIds == null ? 0 : psmIds.size();
			r.anchor = sizeOf(catPeps, "Ribosomal");
			r.highRes = sizeOf(catPeps, "Replication");
			r.growth = sizeOf(catPeps, "Elongation");
			r.housekeeping = sizeOf(catPeps, "Synthetase");

			r.uAnchor = sizeOf(uniqPeps, "Ribosomal");
			r.uHighRes = sizeOf(uniqPeps, "Replication");
			r.uGrowth = sizeOf(uniqPeps, "Elongation");
			r.uHousekeeping = sizeOf(uniqPeps, "Synthetase");

			r.totalMarkers = r.anchor + r.highRes + r.growth + r.housekeeping;
			r.totalUnique = r.uAnchor + r.uHighRes + r.uGrowth + r.uHousekeeping;

			// union peptides across categories for greedy
			for (String cat : categories) {
				if (catPeps != null && catPeps.containsKey(cat)) {
					r.allMarkerPeptides.addAll(catPeps.get(cat));
				}
			}

			rows.add(r);
		}
		LOGGER.info("Genome rows built: {}", rows.size());

		// ----------------------------
		// Write genome evidence table (same as your current output)
		// ----------------------------
		try (PrintWriter writer = new PrintWriter(outputGenomeEvidenceTsv)) {
			writer.println(
					"Genome\tRibosomal\tReplication\tElongation\tSynthetase\tUniqueRibosomal\tUniqueReplication\tUniqueElongation\tUniqueSynthetase\tTotalMarkers\tTotalUniqueMarkers\tTotalPSM\tGenusKey\tFamilyKey");
			rows.sort(Comparator.comparing(a -> a.genome));
			for (GenomeRow r : rows) {
				writer.printf(Locale.ROOT, "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s%n", r.genome, r.anchor,
						r.highRes, r.growth, r.housekeeping, r.uAnchor, r.uHighRes, r.uGrowth, r.uHousekeeping,
						r.totalMarkers, r.totalUnique, r.totalPsm, r.genusKey == null ? "NA" : r.genusKey,
						r.familyKey == null ? "NA" : r.familyKey);
			}
		} catch (Exception e) {
			throw new RuntimeException("Error writing genome evidence: " + outputGenomeEvidenceTsv, e);
		}
		LOGGER.info("Genome evidence written: {}", outputGenomeEvidenceTsv);

		// ----------------------------
		// Apply your strategy
		// ----------------------------
		SelectionResult res = selectGenomesStrategy(rows, allMarkerPeptides, coverageThreshold, stepThreshold,
				maximumSelected);
		LOGGER.info("Selection results: ranked={}, seedIncluded={}, selected={}, totalPeptides={}",
				res.ranked != null ? res.ranked.size() : 0,
				res.seedIncluded != null ? res.seedIncluded.size() : 0,
				res.selected != null ? res.selected.size() : 0,
				res.totalPeptides);

		// ----------------------------
		// Write ranked list + reasons
		// ----------------------------
		try (PrintWriter writer = new PrintWriter(outputRankedTsv)) {
			writer.println(
					"Rank\tGenome\tGenusKey\tFamilyKey\tRibosomal\tReplication\tElongation\tSynthetase\tUniqueRibosomal\tUniqueReplication\tUniqueElongation\tUniqueSynthetase\tTotalMarkers\tTotalUniqueMarkers\tTotalPSM\tnCat\tnCatUnique\tReason\tGreedyGain\tTotalCoveredPeptides\tCoverage\tIncrement\tTotalCoveredPSM\tPsmCoverage\tPsmIncrement");

			// Build map for quick access
			Map<String, GenomeRow> byId = rows.stream().collect(Collectors.toMap(x -> x.genome, x -> x, (a, b) -> a));

			int rank = 0;

			// Use final selection order first (seed+greedy), then append the rest of ranked
			// list
			List<String> outOrder = new ArrayList<>(res.selected);
			for (String g : res.ranked) {
				if (!res.selected.contains(g))
					outOrder.add(g);
			}

			HashSet<String> covered = new HashSet<>();
			HashSet<Integer> coveredPsm = new HashSet<>();
			for (String g : outOrder) {
				GenomeRow r = byId.get(g);
				if (r == null)
					continue;
				String reason = res.reason.getOrDefault(g, "");
				int gain = res.marginalGain.getOrDefault(g, 0);

				int before = covered.size();
				if (r.allMarkerPeptides != null) {
					covered.addAll(r.allMarkerPeptides);
				}
				int totalCovered = covered.size();
				int increment = totalCovered - before;
				double coverage = res.totalPeptides == 0 ? 1.0 : ((double) totalCovered / (double) res.totalPeptides);

				int beforePsm = coveredPsm.size();
				Set<Integer> psmIds = genomePsmIds.get(r.genome);
				if (psmIds != null) {
					coveredPsm.addAll(psmIds);
				}
				int totalCoveredPsm = coveredPsm.size();
				int psmIncrement = totalCoveredPsm - beforePsm;
				double psmCoverage = totalMarkerPsms == 0 ? 1.0 : ((double) totalCoveredPsm / (double) totalMarkerPsms);

				writer.printf(Locale.ROOT, "%d\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%d\t%d\t%.6f\t%d\t%d\t%.6f\t%d%n",
						rank++, r.genome, r.genusKey == null ? "NA" : r.genusKey,
						r.familyKey == null ? "NA" : r.familyKey, r.anchor, r.highRes, r.growth,
						r.housekeeping, r.uAnchor, r.uHighRes, r.uGrowth, r.uHousekeeping, r.totalMarkers,
						r.totalUnique, r.totalPsm, r.nCat(), r.nCatUnique(), reason, gain, totalCovered, coverage, increment,
						totalCoveredPsm, psmCoverage, psmIncrement);
			}
		} catch (Exception e) {
			throw new RuntimeException("Error writing ranked TSV: " + outputRankedTsv, e);
		}
		LOGGER.info("Ranked TSV written: {}", outputRankedTsv);

		// ----------------------------
		// Write selected genomes
		// ----------------------------
		try (PrintWriter writer = new PrintWriter(outputSelectedTsv)) {
			writer.println("Genome\tReason");
			for (String g : res.selected) {
				writer.println(g + "\t" + res.reason.getOrDefault(g, ""));
			}
			writer.println();
			writer.println("#TotalMarkerPeptides\t" + res.totalPeptides);
			writer.println("#TotalMarkerPSMs\t" + totalMarkerPsms);
			writer.println("#SelectedGenomes\t" + res.selected.size());
		} catch (Exception e) {
			throw new RuntimeException("Error writing selected TSV: " + outputSelectedTsv, e);
		}
		LOGGER.info("Selected TSV written: {}", outputSelectedTsv);
	}

	public void peptideTestWithStrategy(String in, File markerFasta, File genomesTaxaTsv,
			String outputGenomeEvidenceTsv, String outputRankedTsv, String outputSelectedTsv, double coverageThreshold, // e.g.
																														// 0.95
			int stepThreshold, // keep for compatibility; optional stop
			int maximumSelected) {

		LOGGER.info("MarkerGenomeSelector start: input={}", in);
		LOGGER.info("Representative level: {}", representativeLevel);
		LOGGER.info("Marker FASTA: {}", markerFasta != null ? markerFasta.getAbsolutePath() : "null");
		LOGGER.info("Genomes taxa TSV: {}", genomesTaxaTsv != null ? genomesTaxaTsv.getAbsolutePath() : "null");
		LOGGER.info("Outputs: evidence={}, ranked={}, selected={}",
				outputGenomeEvidenceTsv, outputRankedTsv, outputSelectedTsv);
		LOGGER.info("Params: coverageThreshold={}, stepThreshold={}, maximumSelected={}",
				coverageThreshold, stepThreshold, maximumSelected);

		Map<String, String> proteinToCategory = new HashMap<>();
		try {
			proteinToCategory = loadMarkerProteinCategoriesFromFasta(markerFasta);
		} catch (IOException e) {
			LOGGER.error("Failed to load marker FASTA categories from {}", markerFasta, e);
		}
		LOGGER.info("Protein categories loaded: {}", proteinToCategory.size());

		HashMap<String, TaxaInfo> taxaMap = new HashMap<>();
		try {
			taxaMap = loadGenomeTaxaMap(genomesTaxaTsv);
		} catch (Exception e) {
			LOGGER.error("Failed to load taxa TSV from {}", genomesTaxaTsv, e);
		}
		LOGGER.info("Taxa map loaded: {}", taxaMap.size());

		// genome -> category -> peptides
		Map<String, Map<String, Set<String>>> genomeCategoryPeptides = new HashMap<>();

		// peptide -> genomes (only for marker peptides)
		Map<String, Set<String>> peptideGenomes = new HashMap<>();

		// peptide -> genome -> categories (to assign unique peptides to categories)
		Map<String, Map<String, Set<String>>> peptideGenomeCategories = new HashMap<>();

		// track unique peptides by genome/category
		Map<String, Map<String, Set<String>>> genomeCategoryUniquePeptides = new HashMap<>();

		// global peptide set (marker peptides only)
		HashSet<String> allMarkerPeptides = new HashSet<>();

		// genome -> PSM ids (marker-mapped PSMs)
		Map<String, Set<Integer>> genomePsmIds = new HashMap<>();
		int totalMarkerPsms = 0;
		int psmIndex = 0;

		// ----------------------------
		// Parse input
		// ----------------------------
		int markerProteinHits = 0;
		int markerPeptideHits = 0;
		try (BufferedReader hapReader = new BufferedReader(new FileReader(in))) {
			String headerLine = hapReader.readLine();
			if (headerLine == null)
				return;
			String[] header = headerLine.split("\t");

			int peptidesIndex = -1;
			int proteinsIndex = -1;
			int labelIndex = -1;
			int spectrumQIndex = -1;

			for (int j = 0; j < header.length; j++) {
				String h = header[j].toLowerCase();
				if (h.equals("peptide"))
					peptidesIndex = j;
				else if (h.equals("proteins"))
					proteinsIndex = j;
				else if (h.equals("label"))
					labelIndex = j;
				else if (h.equals("spectrum_q"))
					spectrumQIndex = j;
			}

			if (peptidesIndex < 0 || proteinsIndex < 0 || labelIndex < 0 || spectrumQIndex < 0) {
				throw new RuntimeException("Missing required columns: peptide/proteins/label/peptide_q");
			}
			LOGGER.info("Column indices: peptide={}, proteins={}, label={}, peptide_q={}",
					peptidesIndex, proteinsIndex, labelIndex, spectrumQIndex);

			String line;
			while ((line = hapReader.readLine()) != null) {
				String[] cs = line.split("\t");
				if (cs.length <= Math.max(Math.max(peptidesIndex, proteinsIndex), Math.max(labelIndex, spectrumQIndex)))
					continue;

				if (!"1".equals(cs[labelIndex]))
					continue;
				if (Double.parseDouble(cs[spectrumQIndex]) > 0.01)
					continue;

				String peptide = cs[peptidesIndex];
				if (peptide == null || peptide.isEmpty())
					continue;

				String[] pros = cs[proteinsIndex].split(";");
				if (pros.length == 0)
					continue;

				int currentPsmId = psmIndex++;

				// collect marker mappings for this peptide
				boolean peptideHasAnyMarker = false;
				HashSet<String> genomesForPsm = new HashSet<>();

				for (String pro : pros) {
					String proteinId = pro.trim();
					if (proteinId.isEmpty())
						continue;

					String category = proteinToCategory.get(proteinId);
					String proteinForGenome = proteinId;
					if (category == null) {
						int last = proteinId.lastIndexOf("_");
						if (last > 0) {
							String baseId = proteinId.substring(0, last);
							category = proteinToCategory.get(baseId);
							if (category != null) {
								proteinForGenome = baseId;
							}
						}
					}
					if (category == null)
						continue;

					String genomeId = extractGenomeId(proteinForGenome);
					if (genomeId == null || genomeId.isEmpty())
						continue;

					peptideHasAnyMarker = true;
					markerProteinHits++;

					peptideGenomes.computeIfAbsent(peptide, k -> new HashSet<>()).add(genomeId);
					genomesForPsm.add(genomeId);

					Map<String, Set<String>> gmap = genomeCategoryPeptides.computeIfAbsent(genomeId,
							k -> new HashMap<>());
					gmap.computeIfAbsent(category, k -> new HashSet<>()).add(peptide);

					Map<String, Set<String>> gcat = peptideGenomeCategories.computeIfAbsent(peptide,
							k -> new HashMap<>());
					gcat.computeIfAbsent(genomeId, k -> new HashSet<>()).add(category);
				}

				// only store peptide if it maps to at least one marker protein
				if (peptideHasAnyMarker) {
					allMarkerPeptides.add(peptide);
					markerPeptideHits++;
					totalMarkerPsms++;
					for (String genomeId : genomesForPsm) {
						genomePsmIds.computeIfAbsent(genomeId, k -> new HashSet<>()).add(currentPsmId);
					}
				}
			}

		} catch (Exception e) {
			throw new RuntimeException("Error parsing input file: " + in, e);
		}
		LOGGER.info("Marker peptides (unique): {}", allMarkerPeptides.size());
		LOGGER.info("Peptide->genomes map size: {}", peptideGenomes.size());
		LOGGER.info("Genome category peptides: {}", genomeCategoryPeptides.size());
		LOGGER.info("Marker protein hits: {}, marker peptide hits: {}", markerProteinHits, markerPeptideHits);

		// ----------------------------
		// Determine unique marker peptides (peptide maps to exactly 1 genome)
		// and assign to categories seen for that genome
		// ----------------------------
		for (Map.Entry<String, Set<String>> entry : peptideGenomes.entrySet()) {
			String peptide = entry.getKey();
			Set<String> genomes = entry.getValue();
			if (genomes.size() != 1)
				continue;

			String genomeId = genomes.iterator().next();
			Map<String, Set<String>> genomeCategories = peptideGenomeCategories.get(peptide);
			if (genomeCategories == null)
				continue;

			Set<String> categories = genomeCategories.get(genomeId);
			if (categories == null)
				continue;

			Map<String, Set<String>> genomeUniqueMap = genomeCategoryUniquePeptides.computeIfAbsent(genomeId,
					k -> new HashMap<>());

			for (String cat : categories) {
				genomeUniqueMap.computeIfAbsent(cat, k -> new HashSet<>()).add(peptide);
			}
		}

		// ----------------------------
		// Build GenomeRow list
		// ----------------------------
		List<String> categories = Arrays.asList("Ribosomal", "Replication", "Elongation", "Synthetase");

		List<GenomeRow> rows = new ArrayList<>();
		for (String genomeId : genomeCategoryPeptides.keySet()) {
			GenomeRow r = new GenomeRow();
			r.genome = genomeId;

			TaxaInfo ti = taxaMap.get(genomeId);
			r.clusterKey = (ti == null) ? "NA" : ti.clusterKey;
			if (ti != null && ti.genus != null && ti.genus.length() > 3 && !ti.genus.substring(3).trim().isEmpty()) {
				r.genusKey = normalizeGenusKey(ti.genus); // g__...
			} else if (ti != null && ti.clusterKey != null && !ti.clusterKey.isEmpty()) {
				r.genusKey = normalizeGenusKey(ti.clusterKey); // fallback up ranks
			} else {
				r.genusKey = "NA";
			}
			if (ti != null && ti.family != null && ti.family.length() > 3 && !ti.family.substring(3).trim().isEmpty()) {
				r.familyKey = normalizeFamilyKey(ti.family); // f__...
			} else if (ti != null && ti.clusterKey != null && !ti.clusterKey.isEmpty()) {
				r.familyKey = normalizeFamilyKey(ti.clusterKey); // fallback up ranks
			} else {
				r.familyKey = "NA";
			}
			if (ti != null && ti.order != null && ti.order.length() > 3 && !ti.order.substring(3).trim().isEmpty()) {
				r.orderKey = ti.order;
			} else {
				r.orderKey = "NA";
			}
			if (ti != null && ti.clazz != null && ti.clazz.length() > 3 && !ti.clazz.substring(3).trim().isEmpty()) {
				r.classKey = ti.clazz;
			} else {
				r.classKey = "NA";
			}
			if (ti != null && ti.phylum != null && ti.phylum.length() > 3 && !ti.phylum.substring(3).trim().isEmpty()) {
				r.phylumKey = ti.phylum;
			} else {
				r.phylumKey = "NA";
			}

			Map<String, Set<String>> catPeps = genomeCategoryPeptides.get(genomeId);
			Map<String, Set<String>> uniqPeps = genomeCategoryUniquePeptides.get(genomeId);

			Set<Integer> psmIds = genomePsmIds.get(genomeId);
			r.totalPsm = psmIds == null ? 0 : psmIds.size();
			r.anchor = sizeOf(catPeps, "Ribosomal");
			r.highRes = sizeOf(catPeps, "Replication");
			r.growth = sizeOf(catPeps, "Elongation");
			r.housekeeping = sizeOf(catPeps, "Synthetase");

			r.uAnchor = sizeOf(uniqPeps, "Ribosomal");
			r.uHighRes = sizeOf(uniqPeps, "Replication");
			r.uGrowth = sizeOf(uniqPeps, "Elongation");
			r.uHousekeeping = sizeOf(uniqPeps, "Synthetase");

			r.totalMarkers = r.anchor + r.highRes + r.growth + r.housekeeping;
			r.totalUnique = r.uAnchor + r.uHighRes + r.uGrowth + r.uHousekeeping;

			// union peptides across categories for greedy
			for (String cat : categories) {
				if (catPeps != null && catPeps.containsKey(cat)) {
					r.allMarkerPeptides.addAll(catPeps.get(cat));
				}
			}

			rows.add(r);
		}
		LOGGER.info("Genome rows built: {}", rows.size());

		// ----------------------------
		// Write genome evidence table (same as your current output)
		// ----------------------------
		try (PrintWriter writer = new PrintWriter(outputGenomeEvidenceTsv)) {
			writer.println(
					"Genome\tRibosomal\tReplication\tElongation\tSynthetase\tUniqueRibosomal\tUniqueReplication\tUniqueElongation\tUniqueSynthetase\tTotalMarkers\tTotalUniqueMarkers\tTotalPSM\tGenusKey\tFamilyKey");
			rows.sort(Comparator.comparing(a -> a.genome));
			for (GenomeRow r : rows) {
				writer.printf(Locale.ROOT, "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s%n", r.genome, r.anchor,
						r.highRes, r.growth, r.housekeeping, r.uAnchor, r.uHighRes, r.uGrowth, r.uHousekeeping,
						r.totalMarkers, r.totalUnique, r.totalPsm, r.genusKey == null ? "NA" : r.genusKey,
						r.familyKey == null ? "NA" : r.familyKey);
			}
		} catch (Exception e) {
			throw new RuntimeException("Error writing genome evidence: " + outputGenomeEvidenceTsv, e);
		}
		LOGGER.info("Genome evidence written: {}", outputGenomeEvidenceTsv);

		// ----------------------------
		// Apply your strategy
		// ----------------------------
		SelectionResult res = selectGenomesStrategy(rows, allMarkerPeptides, coverageThreshold, stepThreshold,
				maximumSelected);
		LOGGER.info("Selection results: ranked={}, seedIncluded={}, selected={}, totalPeptides={}",
				res.ranked != null ? res.ranked.size() : 0,
				res.seedIncluded != null ? res.seedIncluded.size() : 0,
				res.selected != null ? res.selected.size() : 0,
				res.totalPeptides);

		// ----------------------------
		// Write ranked list + reasons
		// ----------------------------
		try (PrintWriter writer = new PrintWriter(outputRankedTsv)) {
			writer.println(
					"Rank\tGenome\tGenusKey\tFamilyKey\tRibosomal\tReplication\tElongation\tSynthetase\tUniqueRibosomal\tUniqueReplication\tUniqueElongation\tUniqueSynthetase\tTotalMarkers\tTotalUniqueMarkers\tTotalPSM\tnCat\tnCatUnique\tReason\tGreedyGain\tTotalCoveredPeptides\tCoverage\tIncrement\tTotalCoveredPSM\tPsmCoverage\tPsmIncrement");

			// Build map for quick access
			Map<String, GenomeRow> byId = rows.stream().collect(Collectors.toMap(x -> x.genome, x -> x, (a, b) -> a));

			int rank = 0;

			// Use final selection order first (seed+greedy), then append the rest of ranked
			// list
			List<String> outOrder = new ArrayList<>(res.selected);
			for (String g : res.ranked) {
				if (!res.selected.contains(g))
					outOrder.add(g);
			}

			HashSet<String> covered = new HashSet<>();
			HashSet<Integer> coveredPsm = new HashSet<>();
			for (String g : outOrder) {
				GenomeRow r = byId.get(g);
				if (r == null)
					continue;
				String reason = res.reason.getOrDefault(g, "");
				int gain = res.marginalGain.getOrDefault(g, 0);

				int before = covered.size();
				if (r.allMarkerPeptides != null) {
					covered.addAll(r.allMarkerPeptides);
				}
				int totalCovered = covered.size();
				int increment = totalCovered - before;
				double coverage = res.totalPeptides == 0 ? 1.0 : ((double) totalCovered / (double) res.totalPeptides);

				int beforePsm = coveredPsm.size();
				Set<Integer> psmIds = genomePsmIds.get(r.genome);
				if (psmIds != null) {
					coveredPsm.addAll(psmIds);
				}
				int totalCoveredPsm = coveredPsm.size();
				int psmIncrement = totalCoveredPsm - beforePsm;
				double psmCoverage = totalMarkerPsms == 0 ? 1.0 : ((double) totalCoveredPsm / (double) totalMarkerPsms);

				writer.printf(Locale.ROOT, "%d\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%d\t%d\t%.6f\t%d\t%d\t%.6f\t%d%n",
						rank++, r.genome, r.genusKey == null ? "NA" : r.genusKey,
						r.familyKey == null ? "NA" : r.familyKey, r.anchor, r.highRes, r.growth,
						r.housekeeping, r.uAnchor, r.uHighRes, r.uGrowth, r.uHousekeeping, r.totalMarkers,
						r.totalUnique, r.totalPsm, r.nCat(), r.nCatUnique(), reason, gain, totalCovered, coverage, increment,
						totalCoveredPsm, psmCoverage, psmIncrement);
			}
		} catch (Exception e) {
			throw new RuntimeException("Error writing ranked TSV: " + outputRankedTsv, e);
		}
		LOGGER.info("Ranked TSV written: {}", outputRankedTsv);

		// ----------------------------
		// Write selected genomes
		// ----------------------------
		try (PrintWriter writer = new PrintWriter(outputSelectedTsv)) {
			writer.println("Genome\tReason");
			for (String g : res.selected) {
				writer.println(g + "\t" + res.reason.getOrDefault(g, ""));
			}
			writer.println();
			writer.println("#TotalMarkerPeptides\t" + res.totalPeptides);
			writer.println("#TotalMarkerPSMs\t" + totalMarkerPsms);
			writer.println("#SelectedGenomes\t" + res.selected.size());
		} catch (Exception e) {
			throw new RuntimeException("Error writing selected TSV: " + outputSelectedTsv, e);
		}
		LOGGER.info("Selected TSV written: {}", outputSelectedTsv);
	}

	// ----------------------------
	// Strategy implementation
	// ----------------------------
	private SelectionResult selectGenomesStrategy(List<GenomeRow> rows, Set<String> totalPeptides,
			double coverageThreshold, int stepThreshold, int maximumSelected) {

		// 0) Hard filter
		List<GenomeRow> filtered = rows.stream().filter(GenomeRow::passesHardFilter).collect(Collectors.toList());
		LOGGER.info("Hard filter genomes: {}", filtered.size());

		// Maps for greedy set cover
		HashMap<String, HashSet<String>> genomeToPeps = new HashMap<>();
		for (GenomeRow r : filtered) {
			genomeToPeps.put(r.genome, new HashSet<>(r.allMarkerPeptides));
		}
		// Only peptides reachable by filtered genomes can be covered
		HashSet<String> effectiveTotalPeptides = new HashSet<>();
		for (HashSet<String> peps : genomeToPeps.values()) {
			if (peps != null) {
				effectiveTotalPeptides.addAll(peps);
			}
		}

		// Helper maps for evidence lookups
		Map<String, GenomeRow> byId = filtered.stream().collect(Collectors.toMap(x -> x.genome, x -> x, (a, b) -> a));

		SelectionResult res = new SelectionResult();
		res.reason = new HashMap<>();
		res.marginalGain = new HashMap<>();
		res.seedIncluded = new LinkedHashSet<>();
		res.selected = new LinkedHashSet<>();
		res.totalPeptides = effectiveTotalPeptides.size();

		// 1) Privileged genomes: >=3 unique marker TYPES
		// Comparator: TotalMarkers desc (for rep selection)
		Comparator<GenomeRow> UNIQUE_FIRST_REP = (a, b) -> {
			if (b.totalMarkers != a.totalMarkers)
				return b.totalMarkers - a.totalMarkers;
			if (b.totalPsm != a.totalPsm)
				return b.totalPsm - a.totalPsm;
			return a.genome.compareTo(b.genome);
		};

		List<GenomeRow> repRanked = new ArrayList<>(filtered);
		repRanked.sort(UNIQUE_FIRST_REP);
		res.ranked = repRanked.stream().map(x -> x.genome).collect(Collectors.toList());

		HashSet<String> coveredDuringRep = new HashSet<>();
		for (GenomeRow r : repRanked) {
			if (r.nCatUnique() < 3) {
				continue;
			}
			HashSet<String> peps = genomeToPeps.get(r.genome);
			int increment = peps == null ? 0 : incrementalGain(coveredDuringRep, peps);
			double minIncrement = coveredDuringRep.size() * 0.01;
			if (increment < minIncrement) {
				continue;
			}
			res.seedIncluded.add(r.genome);
			res.reason.put(r.genome, "PRIVILEGE_UNIQUE_3TYPES");
			if (peps != null) {
				coveredDuringRep.addAll(peps);
			}
		}
		LOGGER.info("Privileged genomes: {}", res.seedIncluded.size());

		// 2) Stepped representative selection (auto level)
		Set<String> representedClass = new HashSet<>();
		Set<String> representedOrder = new HashSet<>();
		Set<String> representedFamily = new HashSet<>();
		Set<String> representedGenus = new HashSet<>();

		// Mark represented keys already covered by privileged set
		for (String g : res.seedIncluded) {
			GenomeRow rr = byId.get(g);
			if (rr == null)
				continue;
			if (rr.classKey != null && !"NA".equals(rr.classKey)) {
				representedClass.add(rr.classKey);
			}
			if (rr.orderKey != null && !"NA".equals(rr.orderKey)) {
				representedOrder.add(rr.orderKey);
			}
			if (rr.familyKey != null && !"NA".equals(rr.familyKey) && !isNumericLikeKey(rr.familyKey)) {
				representedFamily.add(rr.familyKey);
			}
			if (rr.genusKey != null && !"NA".equals(rr.genusKey) && !isNumericLikeKey(rr.genusKey)) {
				representedGenus.add(rr.genusKey);
			}
		}

		for (GenomeRow r : repRanked) {
			if (!r.hasAll4Cats())
				continue;

			if (isSubsetOfSelected(genomeToPeps.get(r.genome), res.seedIncluded, genomeToPeps)) {
				continue;
			}

			HashSet<String> peps = genomeToPeps.get(r.genome);
			int increment = peps == null ? 0 : incrementalGain(coveredDuringRep, peps);
			double minIncrement = coveredDuringRep.size() * 0.01;
			if (increment < minIncrement) {
				continue;
			}
			String reason = null;

			// Step 1: class rep (base rep criteria)
			if (r.classKey != null && !"NA".equals(r.classKey) && !representedClass.contains(r.classKey)) {
				reason = "CLASS_REP";
				representedClass.add(r.classKey);
			} else if (r.orderKey != null && !"NA".equals(r.orderKey) && !representedOrder.contains(r.orderKey)
					&& increment >= 1) {
				reason = "ORDER_REP";
				representedOrder.add(r.orderKey);
			} else if (r.familyKey != null && !"NA".equals(r.familyKey) && !isNumericLikeKey(r.familyKey)
					&& !representedFamily.contains(r.familyKey) && increment >= 2) {
				reason = "FAMILY_REP";
				representedFamily.add(r.familyKey);
			} else if (r.genusKey != null && !"NA".equals(r.genusKey) && !isNumericLikeKey(r.genusKey)
					&& !representedGenus.contains(r.genusKey) && increment >= 3) {
				reason = "GENUS_REP";
				representedGenus.add(r.genusKey);
			}

			if (reason == null) {
				continue;
			}

			res.seedIncluded.add(r.genome);
			res.reason.put(r.genome, appendReason(res.reason.get(r.genome), reason));

			if (peps != null) {
				coveredDuringRep.addAll(peps);
			}
		}
		LOGGER.info("Seed genomes after reps: {}", res.seedIncluded.size());

		// 3) Greedy set cover (your old method idea), starting from seedIncluded
		// We implement a seeded version of your Occam-style greedy:
		// - Start covered set with seed peptides
		// - Then pick next genome by max incremental new peptides (recomputed)
		HashSet<String> covered = new HashSet<>();
		for (String g : res.seedIncluded) {
			HashSet<String> peps = genomeToPeps.get(g);
			if (peps != null)
				covered.addAll(peps);
		}
		LOGGER.info("Covered peptides after seeds: {}", covered.size());

		// Add seeds into selected list first
		for (String g : res.seedIncluded) {
			if (!genomeToPeps.containsKey(g))
				continue; // might have been filtered out
			res.selected.add(g);
			res.marginalGain.put(g, 0);
		}

		int totalCount = effectiveTotalPeptides.size();
		double currentCoverage = totalCount == 0 ? 1.0 : ((double) covered.size() / (double) totalCount);

		// Greedy continue until coverage or max
		while (res.selected.size() < maximumSelected) {

			String best = null;
			int bestGain = -1;
			int bestTotalMarkers = -1;

			for (String g : genomeToPeps.keySet()) {
				if (res.selected.contains(g))
					continue;

				HashSet<String> peps = genomeToPeps.get(g);
				if (peps == null || peps.isEmpty())
					continue;

				int gain = incrementalGain(covered, peps);
				int totalMarkers = 0;
				GenomeRow row = byId.get(g);
				if (row != null) {
					totalMarkers = row.totalMarkers;
				}
				if (gain > bestGain
						|| (gain == bestGain && totalMarkers > bestTotalMarkers)
						|| (gain == bestGain && totalMarkers == bestTotalMarkers
								&& (best == null || g.compareTo(best) < 0))) {
					bestGain = gain;
					bestTotalMarkers = totalMarkers;
					best = g;
				}
			}

			if (best == null)
				break;

			// If already reached 90% coverage, only continue if gain > 1
			if (currentCoverage >= coverageThreshold && bestGain <= coverageThreshold) {
				break;
			}

			// select it
			res.selected.add(best);
			res.reason.putIfAbsent(best, "GREEDY");
			res.marginalGain.put(best, bestGain);

			// update covered
			covered.addAll(genomeToPeps.get(best));
			currentCoverage = totalCount == 0 ? 1.0 : ((double) covered.size() / (double) totalCount);
			if (currentCoverage >= coverageThreshold) {
				break;
			}
		}
		LOGGER.info("Greedy complete: selected={}, coverage={}", res.selected.size(), currentCoverage);

		return res;
	}

	private static int incrementalGain(Set<String> covered, Set<String> candidate) {
		int gain = 0;
		for (String p : candidate) {
			if (!covered.contains(p))
				gain++;
		}
		return gain;
	}

	private static String appendReason(String existing, String added) {
		if (added == null || added.isEmpty()) {
			return existing;
		}
		if (existing == null || existing.isEmpty()) {
			return added;
		}
		if (existing.contains(added)) {
			return existing;
		}
		return existing + "+" + added;
	}

	private static int sizeOf(Map<String, Set<String>> map, String key) {
		if (map == null)
			return 0;
		Set<String> s = map.get(key);
		return s == null ? 0 : s.size();
	}

	private static boolean isSubsetOfSelected(Set<String> candidate, Set<String> selected,
			Map<String, HashSet<String>> genomeToPeps) {
		if (candidate == null || candidate.isEmpty()) {
			return false;
		}
		for (String genome : selected) {
			HashSet<String> existing = genomeToPeps.get(genome);
			if (existing != null && existing.containsAll(candidate)) {
				return true;
			}
		}
		return false;
	}

	private String normalizeGenusKey(String raw) {
		if (raw == null || raw.isEmpty() || "NA".equals(raw)) {
			return "NA";
		}
		if (raw.startsWith("g__")) {
			String genus = raw.substring(3).trim();
			int underscore = genus.indexOf('_');
			if (underscore > 0) {
				genus = genus.substring(0, underscore);
			}
			return "g__" + genus;
		}
		return raw;
	}

	private String normalizeFamilyKey(String raw) {
		if (raw == null || raw.isEmpty() || "NA".equals(raw)) {
			return "NA";
		}
		if (raw.startsWith("f__")) {
			String family = raw.substring(3).trim();
			int underscore = family.indexOf('_');
			if (underscore > 0) {
				family = family.substring(0, underscore);
			}
			return "f__" + family;
		}
		return raw;
	}

	private boolean isNumericLikeKey(String key) {
		if (key == null || key.isEmpty() || "NA".equals(key)) {
			return false;
		}
		String name = key.startsWith("g__") || key.startsWith("f__") ? key.substring(3) : key;
		for (int i = 0; i < name.length(); i++) {
			if (Character.isDigit(name.charAt(i))) {
				return true;
			}
		}
		return false;
	}

	private String getRepresentativeKey(GenomeRow r) {
		if (r == null) {
			return "NA";
		}
		String genusKey = (r.genusKey == null || r.genusKey.isEmpty()) ? "NA" : r.genusKey;
		String familyKey = (r.familyKey == null || r.familyKey.isEmpty()) ? "NA" : r.familyKey;
		String orderKey = (r.orderKey == null || r.orderKey.isEmpty()) ? "NA" : r.orderKey;
		String classKey = (r.classKey == null || r.classKey.isEmpty()) ? "NA" : r.classKey;
		String phylumKey = (r.phylumKey == null || r.phylumKey.isEmpty()) ? "NA" : r.phylumKey;

		switch (representativeLevel) {
		case PHYLUM:
			return phylumKey;
		case CLASS:
			return classKey;
		case ORDER:
			return orderKey;
		case FAMILY:
			return isNumericLikeKey(familyKey) ? "NA" : familyKey;
		case GENUS:
		default:
			if (isNumericLikeKey(genusKey)) {
				if ("NA".equals(familyKey) || isNumericLikeKey(familyKey)) {
					return "NA";
				}
				return familyKey;
			}
			return genusKey;
		}
	}

	public void setRepresentativeLevel(String level) {
		if (level == null || level.trim().isEmpty()) {
			this.representativeLevel = RepresentativeLevel.GENUS;
			return;
		}
		try {
			this.representativeLevel = RepresentativeLevel.valueOf(level.trim().toUpperCase(Locale.ROOT));
		} catch (IllegalArgumentException e) {
			this.representativeLevel = RepresentativeLevel.GENUS;
		}
	}

	// ----------------------------
	// Your existing helpers (you already have these)
	// ----------------------------
	protected Map<String, String> loadMarkerProteinCategoriesFromFasta(File markerFasta) throws IOException {
		Map<String, String> map = new HashMap<>();
		if (markerFasta == null || !markerFasta.isFile()) {
			return map;
		}
		try (BufferedReader reader = new BufferedReader(new FileReader(markerFasta))) {
			String line;
			while ((line = reader.readLine()) != null) {
				if (!line.startsWith(">")) {
					continue;
				}
				String trimmed = line.substring(1).trim();
				if (trimmed.isEmpty()) {
					continue;
				}
				String[] parts = trimmed.split("\\s+", 2);
				String proteinId = parts[0];
				String category = categorize(line.toLowerCase(Locale.ROOT));
				if (category != null) {
					map.put(proteinId, category);
				}
			}
		}
		return map;
	}

	protected String extractGenomeId(String proteinId) {
		if (proteinId == null) {
			return "";
		}
		int underscore = proteinId.indexOf('_');
		if (underscore > 0) {
			return proteinId.substring(0, underscore);
		}
		return proteinId;
	}

	private String categorize(String lineLower) {
		if (ANCHOR_PATTERN.matcher(lineLower).find()) {
			return "Ribosomal";
		}
		if (HIGH_RES_PATTERN.matcher(lineLower).find()) {
			return "Replication";
		}
		if (GROWTH_PATTERN.matcher(lineLower).find()) {
			return "Elongation";
		}
		if (STABILITY_PATTERN.matcher(lineLower).find()) {
			return "Synthetase";
		}
		return null;
	}
	
	public static void exportPeptideFasta(String in, String out, double threshold) {
		try (BufferedReader reader = new BufferedReader(new FileReader(in))) {
			PrintWriter writer = new PrintWriter(out);
			String line = reader.readLine();
			int count = 0;
			while((line=reader.readLine())!=null) {
				String[] cs = line.split("\t");
				if(Double.parseDouble(cs[2])>threshold) {
					writer.println(">"+cs[0]+"_"+count++);
					writer.println(cs[1]);
				}
			}
			reader.close();
			writer.close();
		} catch (IOException e) {

		}
	}
	
	public void writeFasta(String selectedFile, String dbFile, String dbFolder) {
		HashSet<String> genomeSet = new HashSet<>();
		try (BufferedReader reader = new BufferedReader(new FileReader(selectedFile))) {
			String line = reader.readLine(); // header
			while ((line = reader.readLine()) != null) {
				line = line.trim();
				if (line.isEmpty() || line.startsWith("#")) {
					break;
				}
				String[] parts = line.split("\t");
				if (parts.length > 0 && !parts[0].trim().isEmpty()) {
					genomeSet.add(parts[0].trim());
				}
			}
		} catch (IOException e) {

		}

		try (PrintWriter writer = new PrintWriter(dbFile)) {
			for (String genome : genomeSet) {
				File genomeDbFile = new File(dbFolder, genome + ".faa");
				if (!genomeDbFile.exists()) {
					genomeDbFile = new File(dbFolder, genome + ".fasta");
				}
				if (genomeDbFile.exists()) {
					try (BufferedReader reader = new BufferedReader(new FileReader(genomeDbFile))) {
						String line;
						while ((line = reader.readLine()) != null) {
							writer.println(line);
						}
					}
				}
			}
		} catch (IOException e) {

		}
	}
	
	public void fastaTest(String in, String out) {
		try (BufferedReader reader = new BufferedReader(new FileReader(in))) {
			String line = null;
			HashMap<String, HashSet<String>> genomeProMap = new HashMap<String, HashSet<String>>();
			while ((line = reader.readLine()) != null) {
				if (line.startsWith(">")) {
					int id = line.indexOf(" ");
					if (id > 1) {
						String pro = line.substring(1, id);
						String genome = pro.substring(0, pro.indexOf("_"));
						genomeProMap.computeIfAbsent(genome, k -> new HashSet<>()).add(line);
					}
				}
			}
			reader.close();
			
			HashSet<String> genomeSet = new HashSet<String>();
			for (String genome : genomeProMap.keySet()) {
				int[] count = new int[4];
				HashSet<String> proset = genomeProMap.get(genome);
				for (String pro : proset) {
					if (ANCHOR_PATTERN.matcher(pro).find()) {
						count[0] = 1;
					}
					if (HIGH_RES_PATTERN.matcher(pro).find()) {
						count[1] = 1;
					}
					if (GROWTH_PATTERN.matcher(pro).find()) {
						count[2] = 1;
					}
					if (STABILITY_PATTERN.matcher(pro).find()) {
						count[3] = 1;
					}
				}
				int total = count[0]+count[1]+count[2]+count[3];
//				System.out.println(genome+"\t"+total+"\t"+proset.size());
				if(total>=3) {
					genomeSet.add(genome);
				}
			}
			System.out.println(genomeSet.size());
			
			BufferedReader reader2 = new BufferedReader(new FileReader(in));
			PrintWriter writer = new PrintWriter(out);
			boolean use = false;
			
			while ((line = reader2.readLine()) != null) {
				if (line.startsWith(">")) {
					String genome;
					if (line.startsWith(">shuf_")) {
						genome = line.substring(6, line.lastIndexOf("_"));
					} else {
						genome = line.substring(1, line.indexOf("_"));
					}
					if (genomeSet.contains(genome)) {
						writer.println(line);
						use = true;
					} else {
						use = false;
					}
				} else {
					if (use) {
						writer.println(line);
					}
				}
			}
			reader2.close();
			writer.close();
		} catch (IOException e) {

		}
	}

	public void combineFasta(String f1, String f2, String out) {
		try (PrintWriter writer = new PrintWriter(out)) {
			HashSet<String> s1 = new HashSet<String>();
			BufferedReader reader1 = new BufferedReader(new FileReader(f1));
			String line = null;
			while ((line = reader1.readLine()) != null) {
				writer.println(line);
				if (line.startsWith(">")) {
					String genome = line.substring(1, line.indexOf("_"));
					s1.add(genome);
				}
			}
			reader1.close();

			boolean use = false;
			BufferedReader reader2 = new BufferedReader(new FileReader(f2));
			while ((line = reader2.readLine()) != null) {

				if (line.startsWith(">")) {
					String genome = line.substring(1, line.indexOf("_"));
					if (!s1.contains(genome)) {
						use = true;
					} else {
						use = false;
					}
					if (use) {
						writer.println(line);
					}
				} else {
					if (use) {
						writer.println(line);
					}
				}
			}
			reader2.close();
			writer.close();
		} catch (IOException e) {

		}
	}
	
	public void combineTest(String in, File markerFasta, File genomesTaxaTsv) {
		File[] separateResultFolders = (new File(in)).listFiles();
		ArrayList<PFindPSM> psmList = new ArrayList<>();
		for (int i = 0; i < separateResultFolders.length; i++) {
			File spectraFile = new File(separateResultFolders[i], "pFind-Filtered.spectra");
			if (!spectraFile.exists() || spectraFile.length() == 0) {
				System.out.println(separateResultFolders[i].getName());
				continue;
			}
			PFindPsmReader psmReader = new PFindPsmReader(spectraFile);
			PFindPSM[] psms = psmReader.getPsms();
			for (int j = 0; j < psms.length; j++) {
				psmList.add(psms[j]);
			}
			System.out.println(separateResultFolders[i].getName() + "\t" + psms.length);
		}

		File resultFolder = new File(in, "combined");
		resultFolder.mkdir();
		PFindPSM[] psms = psmList.toArray(PFindPSM[]::new);
		File evidenceTsv = new File(resultFolder, "evidence.tsv");
		File rankedTsv = new File(resultFolder, "ranked.tsv");
		File selectedTsv = new File(resultFolder, "selected.tsv");

		if(!selectedTsv.exists()) {
			peptideTestWithStrategy(psms, markerFasta, genomesTaxaTsv, evidenceTsv.getAbsolutePath(),
					rankedTsv.getAbsolutePath(), selectedTsv.getAbsolutePath(), 0.95, 1, 300);
		}

		try {
			HashSet<String> set = new HashSet<String>();
			BufferedReader reader = new BufferedReader(new FileReader(selectedTsv));
			String line = reader.readLine();
			while ((line = reader.readLine()) != null) {
				if (line.startsWith("MGYG")) {
					String genome = line.substring(0, line.indexOf("\t"));
					set.add(genome);
				}
			}
			reader.close();

			int[] counts = new int[2];
			L: for (PFindPSM psm : psms) {
				String[] proteins = psm.getProteins();
				boolean find = false;
				for (String pro : proteins) {
					if (pro.startsWith("REV_")) {
						continue L;
					} else {
						String genome = pro.substring(0, pro.indexOf("_"));
						if (set.contains(genome)) {
							find = true;
						}
					}
				}
				counts[0]++;
				if (find) {
					counts[1]++;
				}
			}
			System.out.println(set.size());
			System.out.println(counts[0] + "\t" + counts[1] + "\t" + (double) counts[1] / (double) counts[0]);

		} catch (IOException e) {

		}
	}
	
	public void compareSelect(String in) {
		File[] separateResultFolders = (new File(in)).listFiles();
		for (int i = 0; i < separateResultFolders.length; i++) {
			File hapFile = new File(separateResultFolders[i], "hap");
			try {
				HashSet<String> set1 = new HashSet<String>();
				File s1 = new File(hapFile, "selected.tsv");
				BufferedReader reader1 = new BufferedReader(new FileReader(s1));
				String line = reader1.readLine();
				while ((line = reader1.readLine()) != null) {
					if (line.startsWith("MGYG")) {
						String genome = line.substring(1, line.indexOf("\t"));
						set1.add(genome);
					}
				}
				reader1.close();

				HashSet<String> set2 = new HashSet<String>();
				File s2 = new File(hapFile, "strict_selected.tsv");
				BufferedReader reader2 = new BufferedReader(new FileReader(s2));
				line = reader2.readLine();
				while ((line = reader2.readLine()) != null) {
					if (line.startsWith("MGYG")) {
						String genome = line.substring(1, line.indexOf("\t"));
						set2.add(genome);
					}
				}
				reader2.close();

				HashSet<String> allset = new HashSet<String>();
				allset.addAll(set1);
				allset.addAll(set2);
				System.out.println(separateResultFolders[i].getName() + "\t" + set1.size() + "\t" + set2.size() + "\t"
						+ allset.size());

			} catch (IOException e) {

			}
		}
	}
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		
		File markerFasta = new File(
				"C:\\Users\\daz25hak\\MetaProAI\\database\\MAG\\human-gut\\v2.0.2\\4markers.fasta");
		File taxaTsv = new File(
				"C:\\Users\\daz25hak\\MetaProAI\\database\\MAG\\human-gut\\v2.0.2\\genomes-taxa.tsv");
		File dbFilder = new File(
				"C:\\Users\\daz25hak\\MetaProAI\\database\\MAG\\human-gut\\v2.0.2\\original_db");
		
		File[] folders = (new File("Z:\\PRIDE_dataset\\PXD023217_campi\\SIHUMIx\\MetaPilot\\pfind_MAG_search")).listFiles();
		for (int i = 0; i < folders.length; i++) {
			File resultFolder = new File(folders[i], "hap");
//		File resultFolder = new File("Z:\\PRIDE_dataset\\PXD006118_3mix_type\\Run1\\MetaPilot\\sage_MAG_search\\Run1_C1_2000ng\\hap");
			if (resultFolder.exists()) {
//				SagePsmReader reader = new SagePsmReader(new File(resultFolder, "results.sage.tsv"));
//				System.out.println(folders[i].getName() + "\t" + reader.getPsms().length);
				
				PFindPsmReader reader = new PFindPsmReader(new File(resultFolder, "pFind-Filtered.spectra"));
				
				File newResultFolder = new File(folders[i], "marker4_all");
				newResultFolder.mkdir();
				
				File evidenceTsv = new File(newResultFolder, "evidence.tsv");
				File rankedTsv = new File(newResultFolder, "ranked.tsv");
				File selectedTsv = new File(newResultFolder, "selected.tsv");

				MarkerGenomeSelector selector = new MarkerGenomeSelector();
//				if(!selectedTsv.exists()) {
					selector.peptideTestWithStrategy(reader.getPsms(), markerFasta, taxaTsv,
							evidenceTsv.getAbsolutePath(), rankedTsv.getAbsolutePath(), selectedTsv.getAbsolutePath(), 1.0,
							1, 300);
//				}
				
//				File fastaFile = new File(newResultFolder, folders[i].getName() + ".fasta");
//				selector.writeFasta(selectedTsv.getAbsolutePath(), fastaFile.getAbsolutePath(), dbFilder.getAbsolutePath());
			}
		}
		
//		MarkerGenomeSelector selector = new MarkerGenomeSelector();
//		selector.compareSelect("Z:\\PRIDE_dataset\\PXD012724_rapidAim1\\MetaPilot\\pfind_MAG_search");
//		selector.combineTest("Z:\\PRIDE_dataset\\PXD012724_rapidAim1\\MetaPilot\\pfind_MAG_search", markerFasta,
//				taxaTsv);
//		selector.peptideTestWithStrategy(reader.getPsms(), markerFasta, taxaTsv,
//				evidenceTsv.getAbsolutePath(), rankedTsv.getAbsolutePath(), selectedTsv.getAbsolutePath(), 0.85,
//				1, 300);
		
//		String in = "C:\\Users\\daz25hak\\MetaProAI\\database\\MAG\\human-gut\\v2.0.2";
//		exportPeptideFasta(in+"\\4markers_deepdetect.tsv", in+"\\4markers_highpep_8.fasta", 0.8);
			
//			MarkerGenomeSelector selector = new MarkerGenomeSelector();
//			selector.fastaTest(
//					"Z:\\PRIDE_dataset\\PXD006118_3mix_type\\Run1\\MetaPilot\\sage_MAG_search_09\\old\\combined\\combined.fasta",
//					"Z:\\PRIDE_dataset\\PXD006118_3mix_type\\Run1\\MetaPilot\\sage_MAG_search_09\\old\\combined\\combined_select34.fasta");
//			selector.writeFasta(
//					"Z:\\PRIDE_dataset\\PXD006118_3mix_type\\Run1\\MetaPilot\\sage_MAG_search_09\\Run1_C1_2000ng\\hap\\old\\marker_genome_selected.tsv",
//					"Z:\\PRIDE_dataset\\PXD006118_3mix_type\\Run1\\MetaPilot\\sage_MAG_search_09\\Run1_C1_2000ng\\hap\\old\\Run1_C1_2000ng.fasta", 
//					"C:\\Users\\daz25hak\\MetaProAI\\database\\MAG\\human-gut\\v2.0.2\\original_db");
			
//			selector.combineFasta("Z:\\PRIDE_dataset\\PXD004039\\MetaPilot\\sage_MAG_search\\Human_3\\hap\\Human_3.fasta", 
//					"Z:\\PRIDE_dataset\\PXD004039\\MetaPilot\\sage_MAG_search_rib\\Human_3\\hap\\Human_3.fasta", 
//					"Z:\\PRIDE_dataset\\PXD004039\\MetaPilot\\sage_MAG_search\\Human_3\\hap\\combined_Human_3.fasta");
		}

}
