package uk.ac.quadram.metapilot.task.mag;

import java.io.*;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import org.apache.commons.io.FileUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import uk.ac.quadram.metapilot.dbSearch.pfind.PFindPSM;
import uk.ac.quadram.metapilot.dbSearch.pfind.PFindPsmReader;
import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNParquetReader;
import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNPrecursor;
import uk.ac.quadram.metapilot.task.io.pep.MetaPSM;

public class MarkerGenomeSelector {

	private static final Logger LOGGER = LogManager.getLogger(MarkerGenomeSelector.class);
	public static final String KNEE_MODE_FOCUS = "Focus";
	public static final String KNEE_MODE_COMPREHENSIVE = "Comprehensive";

	public enum RepresentativeLevel {
		GENUS, FAMILY, ORDER, CLASS, PHYLUM
	}

	private RepresentativeLevel representativeLevel = RepresentativeLevel.GENUS;
	private String kneeSelectionMode = KNEE_MODE_FOCUS;

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
		public Map<String, Integer> greedyStep; // step index in greedy order (1-based)
		public Map<String, Double> greedyRatio; // inc/totalPeptides for greedy step
		public Map<String, Double> greedyMA; // moving average ratio (window=5) for greedy step
		public int kneeStep; // knee step on greedy MA curve (chord-distance, 1-based)
		public int kneeStepKneedle; // kneedle-style knee (1-based)
		public int kneeStepCurvature; // max curvature knee (1-based)
		public int kneeStepPiecewise; // piecewise fit breakpoint (1-based)
		public int kneeStepSlope; // slope change knee (1-based)
	}

	/**
	 * Converts DiaNNPrecursor[] to MetaPSM[] by wrapping each precursor in an adapter
	 */
	private MetaPSM[] convertDiaNNPrecursorsToMetaPSM(DiaNNPrecursor[] precursors) {
		if (precursors == null) {
			return new MetaPSM[0];
		}
		MetaPSM[] psms = new MetaPSM[precursors.length];
		for (int i = 0; i < precursors.length; i++) {
			if (precursors[i] != null) {
				psms[i] = new DiaNNPrecursorAdapter(precursors[i]);
			}
		}
		return psms;
	}

	/**
	 * Adapter class to wrap DiaNNPrecursor as MetaPSM
	 */
	private static class DiaNNPrecursorAdapter extends MetaPSM {
		private final DiaNNPrecursor precursor;

		public DiaNNPrecursorAdapter(DiaNNPrecursor precursor) {
			super(
					precursor.getRunString() != null ? precursor.getRunString() : "",
					0, // scan number not available in DiaNNPrecursor
					precursor.getSeqString() != null ? precursor.getSeqString() : "",
					precursor.getModSeqString() != null ? precursor.getModSeqString() : "",
					precursor.getCharge(),
					precursor.getProGroups() != null ? precursor.getProGroups() : new String[0],
					precursor.getQvalue(),
					precursor.getPEP()
			);
			this.precursor = precursor;
		}

		@Override
		public double getExpMass() {
			return precursor.getMass();
		}

		@Override
		public double getRt() {
			// Use preTrans as RT if available, otherwise return 0
			return precursor.getPreTrans();
		}

		@Override
		public double getCalMass() {
			return precursor.getMass();
		}

		@Override
		public double getMassDiff() {
			// Mass difference not directly available, return 0
			return 0.0;
		}

		@Override
		public int getMiss() {
			return precursor.getMissCleavage();
		}
	}

	public void selectMarkerGenomesWithKnee(PFindPSM[] psms, File markerFasta, File genomesTaxaTsv,
			String outputGenomeEvidenceTsv, String outputRankedTsv, String outputSelectedTsv, int maximumSelected,
			boolean enableAutoRep) {
		selectMarkerGenomesWithKnee((MetaPSM[]) psms, markerFasta, genomesTaxaTsv, outputGenomeEvidenceTsv,
				outputRankedTsv, outputSelectedTsv, maximumSelected, enableAutoRep);
	}

	public void selectMarkerGenomesWithKnee(PFindPSM[] psms, String outputGenomeEvidenceTsv,
			String outputRankedTsv, String outputSelectedTsv, int maximumSelected, boolean enableAutoRep) {
		selectMarkerGenomesWithKnee((MetaPSM[]) psms, null, null, outputGenomeEvidenceTsv, outputRankedTsv,
				outputSelectedTsv, maximumSelected, enableAutoRep);
	}

	public void selectMarkerGenomesWithKnee(MetaPSM[] psms, String outputGenomeEvidenceTsv,
			String outputRankedTsv, String outputSelectedTsv, int maximumSelected, boolean enableAutoRep) {
		selectMarkerGenomesWithKnee(psms, null, null, outputGenomeEvidenceTsv, outputRankedTsv, outputSelectedTsv,
				maximumSelected, enableAutoRep);
	}

	public void selectMarkerGenomesWithKnee(DiaNNPrecursor[] precursors, File markerFasta, File genomesTaxaTsv,
			String outputGenomeEvidenceTsv, String outputRankedTsv, String outputSelectedTsv, int maximumSelected,
			boolean enableAutoRep) {
		MetaPSM[] psms = convertDiaNNPrecursorsToMetaPSM(precursors);
		selectMarkerGenomesWithKnee(psms, markerFasta, genomesTaxaTsv, outputGenomeEvidenceTsv, outputRankedTsv,
				outputSelectedTsv, maximumSelected, enableAutoRep);
	}

	public void selectMarkerGenomesWithKnee(DiaNNPrecursor[] precursors, String outputGenomeEvidenceTsv,
			String outputRankedTsv, String outputSelectedTsv, int maximumSelected, boolean enableAutoRep) {
		MetaPSM[] psms = convertDiaNNPrecursorsToMetaPSM(precursors);
		selectMarkerGenomesWithKnee(psms, null, null, outputGenomeEvidenceTsv, outputRankedTsv, outputSelectedTsv,
				maximumSelected, enableAutoRep);
	}

	public void selectMarkerGenomesWithKnee(MetaPSM[] psms, File markerFasta, File genomesTaxaTsv,
			String outputGenomeEvidenceTsv, String outputRankedTsv, String outputSelectedTsv, int maximumSelected,
			boolean enableAutoRep) {

		LOGGER.info("MarkerGenomeSelector start: input=MetaPSM[] (count={})",
				psms == null ? 0 : psms.length);
		LOGGER.info("Representative level: {}", representativeLevel);
		LOGGER.info("Marker FASTA: {}", markerFasta != null ? markerFasta.getAbsolutePath() : "null");
		LOGGER.info("Genomes taxa TSV: {}", genomesTaxaTsv != null ? genomesTaxaTsv.getAbsolutePath() : "null");
		LOGGER.info("Outputs: evidence={}, ranked={}, selected={}",
				outputGenomeEvidenceTsv, outputRankedTsv, outputSelectedTsv);
		LOGGER.info("Params: maximumSelected={}, autoRep={}", maximumSelected, enableAutoRep);

		boolean useMarkerCategories = markerFasta != null && markerFasta.isFile();
		Map<String, String> proteinToCategory = new HashMap<>();
		if (useMarkerCategories) {
			try {
				proteinToCategory = loadMarkerProteinCategoriesFromFasta(markerFasta);
			} catch (IOException e) {
				LOGGER.error("Failed to load marker FASTA categories from {}", markerFasta, e);
			}
			LOGGER.info("Protein categories loaded: {}", proteinToCategory.size());
		}

		HashMap<String, TaxaInfo> taxaMap = new HashMap<>();
		if (genomesTaxaTsv != null && genomesTaxaTsv.isFile()) {
			try {
				taxaMap = loadGenomeTaxaMap(genomesTaxaTsv);
			} catch (Exception e) {
				LOGGER.error("Failed to load taxa TSV from {}", genomesTaxaTsv, e);
			}
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

				String category = "ALL";
				String proteinForGenome = lookupId;
				if (useMarkerCategories) {
					category = proteinToCategory.get(lookupId);
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
				}

				String genomeId = extractGenomeId(proteinForGenome);
				if (genomeId == null || genomeId.isEmpty()) {
					continue;
				}
				
				// don't need decoy here, because decoy genome can't have 4 types of marker proteins, so it can't help the FDR control
				if (isDecoy) {
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
		List<String> categories = useMarkerCategories
				? Arrays.asList("Ribosomal", "Replication", "Elongation", "Synthetase")
				: Arrays.asList("ALL");

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
			if (useMarkerCategories) {
				r.anchor = sizeOf(catPeps, "Ribosomal");
				r.highRes = sizeOf(catPeps, "Replication");
				r.growth = sizeOf(catPeps, "Elongation");
				r.housekeeping = sizeOf(catPeps, "Synthetase");

				r.uAnchor = sizeOf(uniqPeps, "Ribosomal");
				r.uHighRes = sizeOf(uniqPeps, "Replication");
				r.uGrowth = sizeOf(uniqPeps, "Elongation");
				r.uHousekeeping = sizeOf(uniqPeps, "Synthetase");
			} else {
				r.anchor = sizeOf(catPeps, "ALL");
				r.highRes = 0;
				r.growth = 0;
				r.housekeeping = 0;
				r.uAnchor = sizeOf(uniqPeps, "ALL");
				r.uHighRes = 0;
				r.uGrowth = 0;
				r.uHousekeeping = 0;
			}

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
			if (useMarkerCategories) {
				writer.println(
						"Genome\tRibosomal\tReplication\tElongation\tSynthetase\tUniqueRibosomal\tUniqueReplication\tUniqueElongation\tUniqueSynthetase\tTotalMarkers\tTotalUniqueMarkers\tTotalPSM\tGenusKey\tFamilyKey");
			} else {
				writer.println("Genome\tTotalPeptides\tUniquePeptides\tTotalPSM");
			}
			rows.sort(Comparator.comparing(a -> a.genome));
			for (GenomeRow r : rows) {
				if (useMarkerCategories) {
					writer.printf(Locale.ROOT, "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s%n", r.genome,
							r.anchor, r.highRes, r.growth, r.housekeeping, r.uAnchor, r.uHighRes, r.uGrowth,
							r.uHousekeeping, r.totalMarkers, r.totalUnique, r.totalPsm,
							r.genusKey == null ? "NA" : r.genusKey, r.familyKey == null ? "NA" : r.familyKey);
				} else {
					writer.printf(Locale.ROOT, "%s\t%d\t%d\t%d%n", r.genome, r.totalMarkers, r.totalUnique, r.totalPsm);
				}
			}
		} catch (Exception e) {
			throw new RuntimeException("Error writing genome evidence: " + outputGenomeEvidenceTsv, e);
		}
		LOGGER.info("Genome evidence written: {}", outputGenomeEvidenceTsv);

		// ----------------------------
		// Apply your strategy
		// ----------------------------
		boolean usePrivilege = useMarkerCategories;
		SelectionResult res = selectGenomesStrategyKnee(rows, allMarkerPeptides, genomePsmIds, totalMarkerPsms,
				maximumSelected, enableAutoRep, usePrivilege, useMarkerCategories);
		LOGGER.info("Selection results: ranked={}, seedIncluded={}, selected={}, totalPeptides={}",
				res.ranked != null ? res.ranked.size() : 0,
				res.seedIncluded != null ? res.seedIncluded.size() : 0,
				res.selected != null ? res.selected.size() : 0,
				res.totalPeptides);

		// ----------------------------
		// Write ranked list + reasons
		// ----------------------------
		Map<String, GenomeRow> byId = rows.stream().collect(Collectors.toMap(x -> x.genome, x -> x, (a, b) -> a));
		List<String> outOrder = new ArrayList<>(res.selected);
		for (String g : res.ranked) {
			if (!res.selected.contains(g))
				outOrder.add(g);
		}

		try (PrintWriter writer = new PrintWriter(outputRankedTsv)) {
			if (useMarkerCategories) {
				writer.println(
						"Rank\tGenome\tGenusKey\tFamilyKey\tRibosomal\tReplication\tElongation\tSynthetase\tUniqueRibosomal\tUniqueReplication\tUniqueElongation\tUniqueSynthetase\tTotalMarkers\tTotalUniqueMarkers\tTotalPSM\tnCat\tnCatUnique\tReason\tGreedyGain\tGreedyStep\tGreedyRatio\tGreedyMA\tKneeStep\tKneeKneedle\tKneeCurvature\tKneePiecewise\tKneeSlope\tTotalCoveredPeptides\tCoverage\tIncrement\tTotalCoveredPSM\tPsmCoverage\tPsmIncrement");
			} else {
				writer.println(
						"Rank\tGenome\tTotalPeptides\tUniquePeptides\tTotalPSM\tReason\tGreedyGain\tGreedyStep\tGreedyRatio\tGreedyMA\tKneeStep\tKneeKneedle\tKneeCurvature\tKneePiecewise\tKneeSlope\tTotalCoveredPeptides\tCoverage\tIncrement\tTotalCoveredPSM\tPsmCoverage\tPsmIncrement");
			}

			int rank = 0;
			HashSet<String> covered = new HashSet<>();
			HashSet<Integer> coveredPsm = new HashSet<>();
			for (String g : outOrder) {
				GenomeRow r = byId.get(g);
				if (r == null)
					continue;
				String reason = res.reason.getOrDefault(g, "");
				int gain = res.marginalGain.getOrDefault(g, 0);
				Integer gStep = res.greedyStep != null ? res.greedyStep.get(g) : null;
				Double gRatio = res.greedyRatio != null ? res.greedyRatio.get(g) : null;
				Double gMA = res.greedyMA != null ? res.greedyMA.get(g) : null;
				String gStepStr = gStep == null ? "" : String.valueOf(gStep);
				String gRatioStr = gRatio == null ? "" : String.format(Locale.ROOT, "%.6f", gRatio);
				String gMaStr = gMA == null ? "" : String.format(Locale.ROOT, "%.6f", gMA);
				String kneeStr = res.kneeStep > 0 ? String.valueOf(res.kneeStep) : "";
				String kneeKneedleStr = res.kneeStepKneedle > 0 ? String.valueOf(res.kneeStepKneedle) : "";
				String kneeCurvStr = res.kneeStepCurvature > 0 ? String.valueOf(res.kneeStepCurvature) : "";
				String kneePieceStr = res.kneeStepPiecewise > 0 ? String.valueOf(res.kneeStepPiecewise) : "";
				String kneeSlopeStr = res.kneeStepSlope > 0 ? String.valueOf(res.kneeStepSlope) : "";

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

				if (useMarkerCategories) {
					writer.printf(Locale.ROOT,
							"%d\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%.6f\t%d\t%d\t%.6f\t%d%n",
							rank++, r.genome, r.genusKey == null ? "NA" : r.genusKey,
							r.familyKey == null ? "NA" : r.familyKey, r.anchor, r.highRes, r.growth, r.housekeeping,
							r.uAnchor, r.uHighRes, r.uGrowth, r.uHousekeeping, r.totalMarkers, r.totalUnique,
							r.totalPsm, r.nCat(), r.nCatUnique(), reason, gain, gStepStr, gRatioStr, gMaStr, kneeStr,
							kneeKneedleStr, kneeCurvStr, kneePieceStr, kneeSlopeStr, totalCovered, coverage, increment,
							totalCoveredPsm, psmCoverage, psmIncrement);
				} else {
					writer.printf(Locale.ROOT,
							"%d\t%s\t%d\t%d\t%d\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%.6f\t%d\t%d\t%.6f\t%d%n",
							rank++, r.genome, r.totalMarkers, r.totalUnique, r.totalPsm, reason, gain, gStepStr,
							gRatioStr, gMaStr, kneeStr, kneeKneedleStr, kneeCurvStr, kneePieceStr, kneeSlopeStr,
							totalCovered, coverage, increment, totalCoveredPsm, psmCoverage, psmIncrement);
				}
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

	public void selectMarkerGenomesWithKnee(String in, String outputGenomeEvidenceTsv, String outputRankedTsv,
			String outputSelectedTsv, int maximumSelected, boolean enableAutoRep) {
		selectMarkerGenomesWithKnee(in, null, null, outputGenomeEvidenceTsv, outputRankedTsv, outputSelectedTsv,
				maximumSelected, enableAutoRep);
	}

	public void selectMarkerGenomesWithKnee(String in, File markerFasta, File genomesTaxaTsv,
			String outputGenomeEvidenceTsv, String outputRankedTsv, String outputSelectedTsv, int maximumSelected,
			boolean enableAutoRep) {

		LOGGER.info("MarkerGenomeSelector start: input={}", in);
		LOGGER.info("Representative level: {}", representativeLevel);
		LOGGER.info("Marker FASTA: {}", markerFasta != null ? markerFasta.getAbsolutePath() : "null");
		LOGGER.info("Genomes taxa TSV: {}", genomesTaxaTsv != null ? genomesTaxaTsv.getAbsolutePath() : "null");
		LOGGER.info("Outputs: evidence={}, ranked={}, selected={}",
				outputGenomeEvidenceTsv, outputRankedTsv, outputSelectedTsv);
		LOGGER.info("Params: maximumSelected={}, autoRep={}", maximumSelected, enableAutoRep);

		boolean useMarkerCategories = markerFasta != null && markerFasta.isFile();
		Map<String, String> proteinToCategory = new HashMap<>();
		if (useMarkerCategories) {
			try {
				proteinToCategory = loadMarkerProteinCategoriesFromFasta(markerFasta);
			} catch (IOException e) {
				LOGGER.error("Failed to load marker FASTA categories from {}", markerFasta, e);
			}
			LOGGER.info("Protein categories loaded: {}", proteinToCategory.size());
		}

		HashMap<String, TaxaInfo> taxaMap = new HashMap<>();
		if (genomesTaxaTsv != null && genomesTaxaTsv.isFile()) {
			try {
				taxaMap = loadGenomeTaxaMap(genomesTaxaTsv);
			} catch (Exception e) {
				LOGGER.error("Failed to load taxa TSV from {}", genomesTaxaTsv, e);
			}
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

					String category = "ALL";
					String proteinForGenome = proteinId;
					if (useMarkerCategories) {
						category = proteinToCategory.get(proteinId);
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
					}

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
		List<String> categories = useMarkerCategories
				? Arrays.asList("Ribosomal", "Replication", "Elongation", "Synthetase")
				: Arrays.asList("ALL");

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
			if (useMarkerCategories) {
				r.anchor = sizeOf(catPeps, "Ribosomal");
				r.highRes = sizeOf(catPeps, "Replication");
				r.growth = sizeOf(catPeps, "Elongation");
				r.housekeeping = sizeOf(catPeps, "Synthetase");

				r.uAnchor = sizeOf(uniqPeps, "Ribosomal");
				r.uHighRes = sizeOf(uniqPeps, "Replication");
				r.uGrowth = sizeOf(uniqPeps, "Elongation");
				r.uHousekeeping = sizeOf(uniqPeps, "Synthetase");
			} else {
				r.anchor = sizeOf(catPeps, "ALL");
				r.highRes = 0;
				r.growth = 0;
				r.housekeeping = 0;
				r.uAnchor = sizeOf(uniqPeps, "ALL");
				r.uHighRes = 0;
				r.uGrowth = 0;
				r.uHousekeeping = 0;
			}

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
			if (useMarkerCategories) {
				writer.println(
						"Genome\tRibosomal\tReplication\tElongation\tSynthetase\tUniqueRibosomal\tUniqueReplication\tUniqueElongation\tUniqueSynthetase\tTotalMarkers\tTotalUniqueMarkers\tTotalPSM\tGenusKey\tFamilyKey");
			} else {
				writer.println("Genome\tTotalPeptides\tUniquePeptides\tTotalPSM\tGenusKey\tFamilyKey");
			}
			rows.sort(Comparator.comparing(a -> a.genome));
			for (GenomeRow r : rows) {
				if (useMarkerCategories) {
					writer.printf(Locale.ROOT, "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s%n", r.genome,
							r.anchor, r.highRes, r.growth, r.housekeeping, r.uAnchor, r.uHighRes, r.uGrowth,
							r.uHousekeeping, r.totalMarkers, r.totalUnique, r.totalPsm,
							r.genusKey == null ? "NA" : r.genusKey, r.familyKey == null ? "NA" : r.familyKey);
				} else {
					writer.printf(Locale.ROOT, "%s\t%d\t%d\t%d\t%s\t%s%n", r.genome, r.totalMarkers, r.totalUnique,
							r.totalPsm, r.genusKey == null ? "NA" : r.genusKey,
							r.familyKey == null ? "NA" : r.familyKey);
				}
			}
		} catch (Exception e) {
			throw new RuntimeException("Error writing genome evidence: " + outputGenomeEvidenceTsv, e);
		}
		LOGGER.info("Genome evidence written: {}", outputGenomeEvidenceTsv);

		// ----------------------------
		// Apply your strategy
		// ----------------------------
		boolean usePrivilege = useMarkerCategories;
		SelectionResult res = selectGenomesStrategyKnee(rows, allMarkerPeptides, genomePsmIds, totalMarkerPsms,
				maximumSelected, enableAutoRep, usePrivilege, useMarkerCategories);
		LOGGER.info("Selection results: ranked={}, seedIncluded={}, selected={}, totalPeptides={}",
				res.ranked != null ? res.ranked.size() : 0,
				res.seedIncluded != null ? res.seedIncluded.size() : 0,
				res.selected != null ? res.selected.size() : 0,
				res.totalPeptides);

		// ----------------------------
		// Write ranked list + reasons
		// ----------------------------
		try (PrintWriter writer = new PrintWriter(outputRankedTsv)) {
			if (useMarkerCategories) {
				writer.println(
						"Rank\tGenome\tGenusKey\tFamilyKey\tRibosomal\tReplication\tElongation\tSynthetase\tUniqueRibosomal\tUniqueReplication\tUniqueElongation\tUniqueSynthetase\tTotalMarkers\tTotalUniqueMarkers\tTotalPSM\tnCat\tnCatUnique\tReason\tGreedyGain\tGreedyStep\tGreedyRatio\tGreedyMA\tKneeStep\tKneeKneedle\tKneeCurvature\tKneePiecewise\tKneeSlope\tTotalCoveredPeptides\tCoverage\tIncrement\tTotalCoveredPSM\tPsmCoverage\tPsmIncrement");
			} else {
				writer.println(
						"Rank\tGenome\tGenusKey\tFamilyKey\tTotalPeptides\tUniquePeptides\tTotalPSM\tReason\tGreedyGain\tGreedyStep\tGreedyRatio\tGreedyMA\tKneeStep\tKneeKneedle\tKneeCurvature\tKneePiecewise\tKneeSlope\tTotalCoveredPeptides\tCoverage\tIncrement\tTotalCoveredPSM\tPsmCoverage\tPsmIncrement");
			}

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
				Integer gStep = res.greedyStep != null ? res.greedyStep.get(g) : null;
				Double gRatio = res.greedyRatio != null ? res.greedyRatio.get(g) : null;
				Double gMA = res.greedyMA != null ? res.greedyMA.get(g) : null;
				String gStepStr = gStep == null ? "" : String.valueOf(gStep);
				String gRatioStr = gRatio == null ? "" : String.format(Locale.ROOT, "%.6f", gRatio);
				String gMaStr = gMA == null ? "" : String.format(Locale.ROOT, "%.6f", gMA);
				String kneeStr = res.kneeStep > 0 ? String.valueOf(res.kneeStep) : "";
				String kneeKneedleStr = res.kneeStepKneedle > 0 ? String.valueOf(res.kneeStepKneedle) : "";
				String kneeCurvStr = res.kneeStepCurvature > 0 ? String.valueOf(res.kneeStepCurvature) : "";
				String kneePieceStr = res.kneeStepPiecewise > 0 ? String.valueOf(res.kneeStepPiecewise) : "";
				String kneeSlopeStr = res.kneeStepSlope > 0 ? String.valueOf(res.kneeStepSlope) : "";

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

				if (useMarkerCategories) {
					writer.printf(Locale.ROOT,
							"%d\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%.6f\t%d\t%d\t%.6f\t%d%n",
							rank++, r.genome, r.genusKey == null ? "NA" : r.genusKey,
							r.familyKey == null ? "NA" : r.familyKey, r.anchor, r.highRes, r.growth, r.housekeeping,
							r.uAnchor, r.uHighRes, r.uGrowth, r.uHousekeeping, r.totalMarkers, r.totalUnique,
							r.totalPsm, r.nCat(), r.nCatUnique(), reason, gain, gStepStr, gRatioStr, gMaStr, kneeStr,
							kneeKneedleStr, kneeCurvStr, kneePieceStr, kneeSlopeStr, totalCovered, coverage, increment,
							totalCoveredPsm, psmCoverage, psmIncrement);
				} else {
					writer.printf(Locale.ROOT,
							"%d\t%s\t%s\t%s\t%d\t%d\t%d\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%.6f\t%d\t%d\t%.6f\t%d%n",
							rank++, r.genome, r.genusKey == null ? "NA" : r.genusKey,
							r.familyKey == null ? "NA" : r.familyKey, r.totalMarkers, r.totalUnique, r.totalPsm,
							reason, gain, gStepStr, gRatioStr, gMaStr, kneeStr, kneeKneedleStr, kneeCurvStr,
							kneePieceStr, kneeSlopeStr, totalCovered, coverage, increment, totalCoveredPsm,
							psmCoverage, psmIncrement);
				}
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
	// Knee-based stopping rule for the greedy set cover step
	private SelectionResult selectGenomesStrategyKnee(List<GenomeRow> rows, Set<String> totalPeptides,
			Map<String, Set<Integer>> genomePsmIds, int totalMarkerPsms, int maximumSelected, boolean enableAutoRep,
			boolean usePrivilegeSpecies, boolean applyHardFilter) {

		// 0) Hard filter (marker4 only)
		List<GenomeRow> filtered = applyHardFilter
				? rows.stream().filter(GenomeRow::passesHardFilter).collect(Collectors.toList())
				: new ArrayList<>(rows);
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
		res.greedyStep = new HashMap<>();
		res.greedyRatio = new HashMap<>();
		res.greedyMA = new HashMap<>();
		res.kneeStep = 0;
		res.totalPeptides = effectiveTotalPeptides.size();

		HashSet<String> coveredDuringRep = new HashSet<>();
		if (usePrivilegeSpecies) {
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

			// 2) Stepped representative selection (auto level) - optional
			if (enableAutoRep) {
				Set<String> representedClass = new HashSet<>();
				Set<String> representedOrder = new HashSet<>();
				Set<String> representedFamily = new HashSet<>();
				Set<String> representedGenus = new HashSet<>();

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
			}
			LOGGER.info("Seed genomes after reps: {}", res.seedIncluded.size());
		}

		// 3) Greedy order (seeded), then knee on moving-average increment ratio
		HashSet<String> covered = new HashSet<>(coveredDuringRep);
		List<String> greedyOrder = new ArrayList<>();
		List<Integer> incList = new ArrayList<>();
		List<Integer> totalCoveredList = new ArrayList<>();

		HashSet<String> remaining = new HashSet<>(genomeToPeps.keySet());
		remaining.removeAll(res.seedIncluded);

		while (!remaining.isEmpty() && res.seedIncluded.size() + greedyOrder.size() < maximumSelected) {
			String best = null;
			int bestGain = -1;
			int bestTotalMarkers = -1;

			for (String g : remaining) {
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

			if (best == null || bestGain <= 0)
				break;

			greedyOrder.add(best);
			incList.add(bestGain);
			res.marginalGain.put(best, bestGain);
			covered.addAll(genomeToPeps.get(best));
			totalCoveredList.add(covered.size());
			remaining.remove(best);
		}

		// Compute moving-average ratio curve for greedy steps only
		int totalCount = effectiveTotalPeptides.size();
		int T = greedyOrder.size();
		if (T == 0) {
			res.selected.addAll(res.seedIncluded);
			return res;
		}

		double[] ratio = new double[T];
		for (int i = 0; i < T; i++) {
			int totalCovered = (i < totalCoveredList.size()) ? totalCoveredList.get(i) : totalCount;
			ratio[i] = totalCovered == 0 ? 0.0 : ((double) incList.get(i) / (double) totalCovered);
		}

		final int window = 5;
		int kneeStep = T;
		if (T >= window) {
			double[] ma = new double[T - window + 1];
			double sum = 0.0;
			for (int i = 0; i < window; i++) {
				sum += ratio[i];
			}
			ma[0] = sum / window;
			for (int i = window; i < T; i++) {
				sum += ratio[i];
				sum -= ratio[i - window];
				ma[i - window + 1] = sum / window;
			}
			int kneeIndex = computeKneeIndex(ma); // 1-based index in ma
			int kneeIndexKneedle = computeKneedleIndex(ma);
			int kneeIndexCurv = computeCurvatureIndex(ma);
			int kneeIndexPiece = computePiecewiseIndex(ma);
			int kneeIndexSlope = computeSlopeChangeIndex(ma);
			kneeStep = Math.min(T, (window - 1) + kneeIndex);
			res.kneeStepKneedle = Math.min(T, (window - 1) + kneeIndexKneedle);
			res.kneeStepCurvature = Math.min(T, (window - 1) + kneeIndexCurv);
			res.kneeStepPiecewise = Math.min(T, (window - 1) + kneeIndexPiece);
			res.kneeStepSlope = Math.min(T, (window - 1) + kneeIndexSlope);

			// Map moving-average values back to greedy steps
			for (int i = 0; i < T; i++) {
				int step = i + 1;
				res.greedyStep.put(greedyOrder.get(i), step);
				res.greedyRatio.put(greedyOrder.get(i), ratio[i]);
				if (i >= window - 1) {
					res.greedyMA.put(greedyOrder.get(i), ma[i - window + 1]);
				}
			}
		} else {
			for (int i = 0; i < T; i++) {
				int step = i + 1;
				res.greedyStep.put(greedyOrder.get(i), step);
				res.greedyRatio.put(greedyOrder.get(i), ratio[i]);
			}
		}
		res.kneeStep = kneeStep;
		if (T < window) {
			res.kneeStepKneedle = kneeStep;
			res.kneeStepCurvature = kneeStep;
			res.kneeStepPiecewise = kneeStep;
			res.kneeStepSlope = kneeStep;
		}

		// Final selection: privilege + greedy up to selected knee mode
		res.selected.addAll(res.seedIncluded);
		int selectStep = selectKneeStep(res, kneeStep, T, kneeSelectionMode);
		for (int i = 0; i < selectStep; i++) {
			String g = greedyOrder.get(i);
			res.selected.add(g);
			res.reason.putIfAbsent(g, "GREEDY");
		}

		// For ranked output: greedy order first, then leftovers (zero gain) by id
		List<String> leftovers = new ArrayList<>(remaining);
		Collections.sort(leftovers);
		for (String g : leftovers) {
			res.marginalGain.put(g, 0);
		}
		res.ranked = new ArrayList<>(greedyOrder);
		res.ranked.addAll(leftovers);

		int selectedStepLog = selectKneeStep(res, kneeStep, T, kneeSelectionMode);
		LOGGER.info("Greedy (knee) complete: selected={}, mode={}, selectedStep={}, kneeStep={}, greedySteps={}",
				res.selected.size(), kneeSelectionMode, selectedStepLog, kneeStep, T);

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

	/**
	 * Knee detection using maximum distance to the chord between the first and last points.
	 * Returns a 1-based index (k in 1..n).
	 */
	private static int computeKneeIndex(double[] y) {
		if (y == null || y.length == 0)
			return 0;
		if (y.length <= 2)
			return y.length;

		int n = y.length;
		double x1 = 0.0;
		double y1 = y[0];
		double x2 = 1.0;
		double y2 = y[n - 1];
		double denom = Math.sqrt((y2 - y1) * (y2 - y1) + (x2 - x1) * (x2 - x1));
		if (denom == 0.0)
			return n;

		double maxDist = -1.0;
		int bestIndex = 1;
		for (int i = 0; i < n; i++) {
			double x = (n == 1) ? 0.0 : ((double) i / (double) (n - 1));
			double yi = y[i];
			double dist = Math.abs((y2 - y1) * x - (x2 - x1) * yi + x2 * y1 - y2 * x1) / denom;
			if (dist > maxDist) {
				maxDist = dist;
				bestIndex = i + 1; // 1-based
			}
		}
		return bestIndex;
	}

	private static int selectKneeStep(SelectionResult res, int chordStep, int totalSteps, String mode) {
		if (res == null || totalSteps <= 0) {
			return 0;
		}
		int selected = chordStep;
		if (KNEE_MODE_COMPREHENSIVE.equalsIgnoreCase(mode)) {
			int piece = res.kneeStepPiecewise > 0 ? res.kneeStepPiecewise : chordStep;
			selected = Math.max(chordStep, piece);
		} else {
			int kneedle = res.kneeStepKneedle > 0 ? res.kneeStepKneedle : chordStep;
			int curv = res.kneeStepCurvature > 0 ? res.kneeStepCurvature : chordStep;
			int slope = res.kneeStepSlope > 0 ? res.kneeStepSlope : chordStep;
			selected = medianOfThree(kneedle, curv, slope);
		}
		if (selected < 1)
			selected = 1;
		if (selected > totalSteps)
			selected = totalSteps;
		return selected;
	}

	private static int medianOfThree(int a, int b, int c) {
		int[] v = new int[] { a, b, c };
		Arrays.sort(v);
		return v[1];
	}

	// Kneedle-style knee: max diff between normalized curve and line (1 - x)
	private static int computeKneedleIndex(double[] y) {
		if (y == null || y.length == 0)
			return 0;
		if (y.length == 1)
			return 1;
		double min = Arrays.stream(y).min().orElse(0.0);
		double max = Arrays.stream(y).max().orElse(0.0);
		double range = max - min;
		double bestDiff = -Double.MAX_VALUE;
		int bestIdx = 1;
		for (int i = 0; i < y.length; i++) {
			double x = (y.length == 1) ? 0.0 : ((double) i / (double) (y.length - 1));
			double yn = (range == 0.0) ? 0.0 : ((y[i] - min) / range);
			double diff = yn - (1.0 - x);
			if (diff > bestDiff) {
				bestDiff = diff;
				bestIdx = i + 1;
			}
		}
		return bestIdx;
	}

	// Max curvature on smoothed curve (using second differences)
	private static int computeCurvatureIndex(double[] y) {
		if (y == null || y.length == 0)
			return 0;
		if (y.length <= 2)
			return y.length;
		double maxAbs = -1.0;
		int best = 1;
		for (int i = 2; i < y.length; i++) {
			double d1 = y[i] - y[i - 1];
			double d0 = y[i - 1] - y[i - 2];
			double d2 = d1 - d0;
			double abs = Math.abs(d2);
			if (abs > maxAbs) {
				maxAbs = abs;
				best = i + 1;
			}
		}
		return best;
	}

	// Two-segment piecewise fit: choose k that minimizes SSE
	private static int computePiecewiseIndex(double[] y) {
		if (y == null || y.length == 0)
			return 0;
		if (y.length <= 3)
			return y.length;
		double bestSse = Double.MAX_VALUE;
		int bestK = 2;
		for (int k = 2; k <= y.length - 1; k++) {
			double sse1 = linearFitSse(y, 0, k - 1);
			double sse2 = linearFitSse(y, k - 1, y.length - 1);
			double sse = sse1 + sse2;
			if (sse < bestSse) {
				bestSse = sse;
				bestK = k;
			}
		}
		return bestK;
	}

	// Slope change detection: split on slope sequence (two means)
	private static int computeSlopeChangeIndex(double[] y) {
		if (y == null || y.length == 0)
			return 0;
		if (y.length <= 3)
			return y.length;
		double[] slopes = new double[y.length - 1];
		for (int i = 1; i < y.length; i++) {
			slopes[i - 1] = y[i] - y[i - 1];
		}
		double bestSse = Double.MAX_VALUE;
		int bestK = 2;
		for (int k = 1; k <= slopes.length - 1; k++) {
			double sse1 = meanSse(slopes, 0, k - 1);
			double sse2 = meanSse(slopes, k, slopes.length - 1);
			double sse = sse1 + sse2;
			if (sse < bestSse) {
				bestSse = sse;
				bestK = k + 1; // map slope split to point index
			}
		}
		return bestK;
	}

	private static double linearFitSse(double[] y, int start, int end) {
		int n = end - start + 1;
		if (n <= 1)
			return 0.0;
		double sumX = 0.0, sumY = 0.0, sumXY = 0.0, sumXX = 0.0;
		for (int i = start; i <= end; i++) {
			double x = i - start;
			sumX += x;
			sumY += y[i];
			sumXY += x * y[i];
			sumXX += x * x;
		}
		double denom = n * sumXX - sumX * sumX;
		double b = denom == 0.0 ? 0.0 : (n * sumXY - sumX * sumY) / denom;
		double a = (sumY - b * sumX) / n;
		double sse = 0.0;
		for (int i = start; i <= end; i++) {
			double x = i - start;
			double fit = a + b * x;
			double err = y[i] - fit;
			sse += err * err;
		}
		return sse;
	}

	private static double meanSse(double[] y, int start, int end) {
		int n = end - start + 1;
		if (n <= 0)
			return 0.0;
		double sum = 0.0;
		for (int i = start; i <= end; i++) sum += y[i];
		double mean = sum / n;
		double sse = 0.0;
		for (int i = start; i <= end; i++) {
			double err = y[i] - mean;
			sse += err * err;
		}
		return sse;
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

	public void setKneeSelectionMode(String mode) {
		if (mode == null || mode.trim().isEmpty()) {
			this.kneeSelectionMode = KNEE_MODE_FOCUS;
			return;
		}
		if (KNEE_MODE_COMPREHENSIVE.equalsIgnoreCase(mode)) {
			this.kneeSelectionMode = KNEE_MODE_COMPREHENSIVE;
		} else {
			this.kneeSelectionMode = KNEE_MODE_FOCUS;
		}
	}

	public String getKneeSelectionMode() {
		return kneeSelectionMode;
	}

	/**
	 * Sets protein information on precursors from MAG predicted TSV files (peptide→protein mapping).
	 * User-provided library search parquet does not contain genome/protein info; this method fills
	 * proGroups so that MarkerGenomeSelector can rank genomes. Predicted TSV format: column 0 = protein ID,
	 * column 1 = peptide sequence.
	 *
	 * @param precursors      precursors from parquet (e.g. from user library search), may have empty proGroups
	 * @param predictedFolder MAG database folder containing *.predicted.tsv files
	 * @return list of precursors that have at least one matching protein (proGroups set); others are skipped
	 */
	public static List<DiaNNPrecursor> setProteinInfoFromPredictedTsVs(DiaNNPrecursor[] precursors,
			File predictedFolder) {
		if (precursors == null || precursors.length == 0 || predictedFolder == null || !predictedFolder.isDirectory()) {
			return precursors != null && precursors.length > 0
					? new ArrayList<>(Arrays.asList(precursors))
					: new ArrayList<>();
		}
		ConcurrentHashMap<String, Set<String>> pepProMap = new ConcurrentHashMap<>();
		for (DiaNNPrecursor pp : precursors) {
			String sequence = pp != null ? pp.getSeqString() : null;
			if (sequence != null && !sequence.isEmpty()) {
				pepProMap.put(sequence, ConcurrentHashMap.newKeySet());
			}
		}
		File[] files = predictedFolder.listFiles();
		if (files != null && files.length > 0) {
			ExecutorService executor = Executors.newFixedThreadPool(Math.min(20, files.length));
			for (File file : files) {
				if (file.isFile() && file.getName().endsWith(".predicted.tsv")) {
					executor.execute(() -> {
						try (BufferedReader reader = new BufferedReader(new FileReader(file))) {
							reader.readLine(); // header
							String pline;
							while ((pline = reader.readLine()) != null) {
								String[] cs = pline.split("\t");
								if (cs.length > 1 && pepProMap.containsKey(cs[1])) {
									pepProMap.get(cs[1]).add(cs[0]);
								}
							}
						} catch (IOException e) {
							LOGGER.warn("Error reading predicted TSV: {}", file.getAbsolutePath(), e);
						}
					});
				}
			}
			executor.shutdown();
			try {
				executor.awaitTermination(1, TimeUnit.HOURS);
			} catch (InterruptedException e) {
				Thread.currentThread().interrupt();
				LOGGER.warn("Interrupted while reading predicted TSVs");
			}
		}
		ArrayList<DiaNNPrecursor> result = new ArrayList<>();
		for (DiaNNPrecursor pp : precursors) {
			if (pp == null) continue;
			String sequence = pp.getSeqString();
			if (sequence == null || sequence.isEmpty()) continue;
			Set<String> proSet = pepProMap.get(sequence);
			if (proSet != null && !proSet.isEmpty()) {
				pp.setProGroups(proSet.toArray(new String[0]));
				result.add(pp);
			}
		}
		return result;
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
	
	public void writeFasta(String selectedFile, String noselectFile, String dbFile, String dbFolder) {
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
		System.out.println(genomeSet.size());
		try (BufferedReader reader = new BufferedReader(new FileReader(noselectFile))) {
			String line = reader.readLine(); // header
			while ((line = reader.readLine()) != null) {
				line = line.trim();
				if (line.isEmpty() || line.startsWith("#")) {
					break;
				}
				String[] parts = line.split("\t");
				if (parts.length > 0 && !parts[0].trim().isEmpty()) {
					if(genomeSet.contains(parts[0].trim())) {
						genomeSet.remove(parts[0].trim());
					}
				}
			}
		} catch (IOException e) {

		}
		System.out.println(genomeSet.size());
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
			selectMarkerGenomesWithKnee(psms, markerFasta, genomesTaxaTsv, evidenceTsv.getAbsolutePath(),
					rankedTsv.getAbsolutePath(), selectedTsv.getAbsolutePath(), 300, true);
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
	
	void SIHUMIxTest(String fasta, String out) {
		try {
			PrintWriter writer = new PrintWriter(out);
			writer.println("Catalog\tProtein");
			BufferedReader reader = new BufferedReader(new FileReader(fasta));
			String line = null;
			while ((line = reader.readLine()) != null) {
				if (line.startsWith(">")) {
					String catalog = categorize(line);
					writer.println(catalog+"\t"+line.substring(1));
				}
			}
			reader.close();
			writer.close();
		} catch (IOException e) {

		}
	}

	
	void test196(String in) {
		File[] separateResultFolders = (new File(in)).listFiles();
		for (int i = 0; i < separateResultFolders.length; i++) {
			File hapFile = new File(separateResultFolders[i], "hap");
			try {
				HashSet<String> set1 = new HashSet<String>();
				File s1 = new File(hapFile, "marker_genome_selected.tsv");
				boolean find = false;
				BufferedReader reader1 = new BufferedReader(new FileReader(s1));
				String line = reader1.readLine();
				while ((line = reader1.readLine()) != null) {
					if (line.startsWith("MGYG")) {
						String genome = line.substring(1, line.indexOf("\t"));
						set1.add(genome);
						find = true;
					}
				}
				reader1.close();

				System.out.println(separateResultFolders[i].getName() + "\t" + set1.size());

			} catch (IOException e) {

			}
		}
	}
	
	private void copyFolder(String in, String out) {
		File[] files = (new File(in)).listFiles();
		for (int i = 0; i < files.length; i++) {
			String fileName = files[i].getName();
			File newFolder = new File(out, fileName);
			newFolder.mkdir();

			File hapFile = new File(files[i], "hap");
			if (hapFile.exists()) {
				File newHapFile = new File(newFolder, "hap");
				newHapFile.mkdir();

				File[] subFiles = hapFile.listFiles();
				for (int j = 0; j < subFiles.length; j++) {
					String subName = subFiles[j].getName();
					if (subName.startsWith("pFind")) {
						try {
							FileUtils.copyFile(subFiles[j], new File(newHapFile, subName));
						} catch (IOException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
					}
				}
			}
		}
	}
	
	private void get300drugPro(String peptidein, String proteinout, String originalDb) {

		try {
			HashMap<String, double[]> pepIntensityMap = new HashMap<String, double[]>();
			HashMap<String, HashMap<String, HashSet<String>>> genomeProPepMap = new HashMap<String, HashMap<String, HashSet<String>>>();
			BufferedReader reader = new BufferedReader(new FileReader(peptidein));
			String line = reader.readLine();
			String[] title = line.split("\t");
			ArrayList<String> fileList = new ArrayList<String>();
			ArrayList<Integer> fileIdList = new ArrayList<Integer>();
			for (int i = 0; i < title.length; i++) {
				if (title[i].startsWith("Intensity ")) {
					fileList.add(title[i]);
					fileIdList.add(i);
				}
			}
			System.out.println(fileIdList.size());
			while ((line = reader.readLine()) != null) {
				String[] cs = line.split("\t");
				String[] proteins = cs[1].split(";");
				if (proteins.length > 0) {
					for (int i = 0; i < proteins.length; i++) {
						String genome = proteins[i].substring(0, proteins[i].indexOf("_"));
						HashMap<String, HashSet<String>> proPepMap;
						if (genomeProPepMap.containsKey(genome)) {
							proPepMap = genomeProPepMap.get(genome);
						} else {
							proPepMap = new HashMap<String, HashSet<String>>();
							genomeProPepMap.put(genome, proPepMap);
						}
						if (proPepMap.containsKey(proteins[i])) {
							proPepMap.get(proteins[i]).add(cs[0]);
						} else {
							HashSet<String> set = new HashSet<String>();
							set.add(cs[0]);
							proPepMap.put(proteins[i], set);
						}
					}

					double[] intensity = new double[fileIdList.size()];
					for (int i = 0; i < intensity.length; i++) {
						intensity[i] = Double.parseDouble(cs[fileIdList.get(i)]) / ((double) proteins.length);
					}
					pepIntensityMap.put(cs[0], intensity);
				}
			}
			reader.close();
			System.out.println(pepIntensityMap.size() + "\t" + genomeProPepMap.size());

			PrintWriter writer = new PrintWriter(proteinout);
			StringBuilder titlesb = new StringBuilder();
			titlesb.append("Genome\t").append("Protein id\t").append("Protein name\t").append("Marker type");
			for (int i = 0; i < fileList.size(); i++) {
				titlesb.append("\t").append(fileList.get(i));
			}
			writer.println(titlesb);
			for (String genome : genomeProPepMap.keySet()) {
				HashMap<String, HashSet<String>> proPepMap = genomeProPepMap.get(genome);
				File dbFile = new File(originalDb, genome + ".faa");
				reader = new BufferedReader(new FileReader(dbFile));
				while ((line = reader.readLine()) != null) {
					if (line.startsWith(">")) {
						int id = line.indexOf(" ");
						String protein = line.substring(1, id);
						if (proPepMap.containsKey(protein)) {
							double[] proteinIntensity = new double[fileIdList.size()];
							for (String pep : proPepMap.get(protein)) {
								double[] peptideIntensity = pepIntensityMap.get(pep);
								if (peptideIntensity != null) {
									for (int i = 0; i < peptideIntensity.length; i++) {
										proteinIntensity[i] += peptideIntensity[i];
									}
								}
							}
							StringBuilder sb = new StringBuilder();
							sb.append(genome).append("\t");
							sb.append(protein).append("\t");

							String proteinName = line.substring(id + 1);
							String categorize = categorize(proteinName);
							sb.append(proteinName).append("\t");
							if (categorize == null) {
								sb.append("None");
							} else {
								sb.append(categorize);
							}
							for (int i = 0; i < proteinIntensity.length; i++) {
								sb.append("\t").append(proteinIntensity[i]);
							}
							writer.println(sb);
						}
					}
				}
			}
			writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	/**
	 * New-idea test entry point. Call from main() with your parquet folder and options.
	 * Replace or extend the body with your new selection/scoring logic.
	 *
	 * @param parquetFolder folder containing e.g. combined.parquet or marker.parquet
	 */
	public static void testNewIdea(String parquetFolder) {
		File parquetFile = new File(parquetFolder, "combined.parquet");
		if (!parquetFile.exists()) {
			parquetFile = new File(parquetFolder, "marker.parquet");
		}
		if (!parquetFile.exists()) {
			LOGGER.warn("testNewIdea: no combined.parquet or marker.parquet in {}", parquetFolder);
			return;
		}
		File markerFasta = new File("C:\\Users\\daz25hak\\MetaProAI\\database\\MAG\\human-gut\\v2.0.2\\Marker4.fasta");
		File taxaTsv = new File("C:\\Users\\daz25hak\\MetaProAI\\database\\MAG\\human-gut\\v2.0.2\\genomes-taxa.tsv");

		DiaNNParquetReader reader = new DiaNNParquetReader(parquetFile.getAbsolutePath(), "KR");
		DiaNNPrecursor[] precursors = reader.getDiaNNPrecursors(0.01);
		if (precursors == null || precursors.length == 0) {
			LOGGER.warn("testNewIdea: no precursors with Q<0.01");
			return;
		}

		MarkerGenomeSelector selector = new MarkerGenomeSelector();
		selector.setKneeSelectionMode(KNEE_MODE_FOCUS);

		File evidenceTsv = new File(parquetFolder, "new_idea_evidence.tsv");
		File rankedTsv = new File(parquetFolder, "new_idea_ranked.tsv");
		File selectedTsv = new File(parquetFolder, "new_idea_selected.tsv");

		selector.selectMarkerGenomesWithKnee(precursors, markerFasta, taxaTsv,
				evidenceTsv.getAbsolutePath(), rankedTsv.getAbsolutePath(), selectedTsv.getAbsolutePath(), 300, true);

		LOGGER.info("testNewIdea done: evidence={}, ranked={}, selected={}",
				evidenceTsv.getAbsolutePath(), rankedTsv.getAbsolutePath(), selectedTsv.getAbsolutePath());
	}

	private static void testParquet(String parquetFolder) {
		File markerFasta = new File("C:\\Users\\daz25hak\\MetaProAI\\database\\MAG\\human-gut\\v2.0.2\\Marker4.fasta");
		File taxaTsv = new File("C:\\Users\\daz25hak\\MetaProAI\\database\\MAG\\human-gut\\v2.0.2\\genomes-taxa.tsv");
		File dbFilder = new File("C:\\Users\\daz25hak\\MetaProAI\\database\\MAG\\human-gut\\v2.0.2\\original_db");

		DiaNNParquetReader parquetReader = new DiaNNParquetReader(
				(new File(parquetFolder, "marker.parquet")).getAbsolutePath(), "KR");

		MarkerGenomeSelector selector = new MarkerGenomeSelector();
//		selector.setKneeSelectionMode(KNEE_MODE_COMPREHENSIVE);
		selector.setKneeSelectionMode(KNEE_MODE_FOCUS);

		File evidenceTsv = new File(parquetFolder, "genome_score_evidence.tsv");
		File rankedTsv = new File(parquetFolder, "genome_score_ranked.tsv");
		File selectedTsv = new File(parquetFolder, "genome_score_selected.tsv");

		selector.selectMarkerGenomesWithKnee(parquetReader.getDiaNNPrecursors(), markerFasta, taxaTsv,
				evidenceTsv.getAbsolutePath(), rankedTsv.getAbsolutePath(), selectedTsv.getAbsolutePath(), 300, true);
		
//		File fastaFile = new File(parquetFolder, "focus.fasta");
//		selector.writeFasta(selectedTsv.getAbsolutePath(), fastaFile.getAbsolutePath(), dbFilder.getAbsolutePath());
	}
	
	private static void testParquet2(String parquetFolder) {
		File markerFasta = new File("C:\\Users\\daz25hak\\MetaProAI\\database\\MAG\\human-gut\\v2.0.2\\Marker4.fasta");
		File taxaTsv = new File("C:\\Users\\daz25hak\\MetaProAI\\database\\MAG\\human-gut\\v2.0.2\\genomes-taxa.tsv");
		File dbFilder = new File("C:\\Users\\daz25hak\\MetaProAI\\database\\MAG\\human-gut\\v2.0.2\\original_db");
		File predictedFilder = new File("C:\\Users\\daz25hak\\MetaProAI\\database\\MAG\\human-gut\\v2.0.2\\predicted");

		ConcurrentHashMap<String, Set<String>> pepProMap = new ConcurrentHashMap<>();
		DiaNNParquetReader parquetReader = new DiaNNParquetReader(
				(new File(parquetFolder, "marker.parquet")).getAbsolutePath(), "KR");
		DiaNNPrecursor[] pps = parquetReader.getDiaNNPrecursors();
		for(DiaNNPrecursor pp : pps) {
			String sequence = pp.getSeqString();
			pepProMap.put(sequence, ConcurrentHashMap.newKeySet());
		}
		System.out.println(pps.length + "\t" + pepProMap.size());
		
		ExecutorService executor = Executors.newFixedThreadPool(20);
		File[] files = predictedFilder.listFiles();
		for (File file : files) {
			executor.execute(() -> {
				if (file.getName().endsWith(".predicted.tsv")) {
					try (BufferedReader reader = new BufferedReader(new FileReader(file))) {
						String pline = reader.readLine();
						while ((pline = reader.readLine()) != null) {
							String[] cs = pline.split("\t");
							if (pepProMap.containsKey(cs[1])) {
								pepProMap.get(cs[1]).add(cs[0]);
							}
						}
					} catch (IOException e) {

					}
				}
			});
		}

		executor.shutdown();

		try {
			executor.awaitTermination(1, TimeUnit.HOURS);
		} catch (InterruptedException e) {

		}
		
		ArrayList<DiaNNPrecursor> ppsList = new ArrayList<>();
		for (DiaNNPrecursor pp : pps) {
			String sequence = pp.getSeqString();
			if (pepProMap.containsKey(sequence)) {
				Set<String> proSet = pepProMap.get(sequence);
				if (proSet.size() > 0) {
					String[] pros = proSet.toArray(new String[proSet.size()]);
					pp.setProGroups(pros);
					ppsList.add(pp);
				}
			}
		}
		System.out.println(pps.length + "\t" + ppsList.size() + "\t" + pepProMap.size());
		MarkerGenomeSelector selector = new MarkerGenomeSelector();
//		selector.setKneeSelectionMode(KNEE_MODE_COMPREHENSIVE);
		selector.setKneeSelectionMode(KNEE_MODE_FOCUS);

		File evidenceTsv = new File(parquetFolder, "nomarker_focus_evidence.tsv");
		File rankedTsv = new File(parquetFolder, "nomarker_focus_ranked.tsv");
		File selectedTsv = new File(parquetFolder, "nomarker_focus_selected.tsv");

//		selector.selectMarkerGenomesWithKnee(ppsList.toArray(new DiaNNPrecursor[ppsList.size()]), markerFasta, taxaTsv,
//				evidenceTsv.getAbsolutePath(), rankedTsv.getAbsolutePath(), selectedTsv.getAbsolutePath(), 300, true);
		
		selector.selectMarkerGenomesWithKnee(ppsList.toArray(new DiaNNPrecursor[ppsList.size()]),
				evidenceTsv.getAbsolutePath(), rankedTsv.getAbsolutePath(), selectedTsv.getAbsolutePath(), 300, true);

//		File fastaFile = new File(parquetFolder, "focus.fasta");
//		selector.writeFasta(selectedTsv.getAbsolutePath(), fastaFile.getAbsolutePath(), dbFilder.getAbsolutePath());
	}
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		
		File markerFasta = new File(
				"C:\\Users\\daz25hak\\MetaProAI\\database\\MAG\\human-gut\\v2.0.2\\Marker4.fasta");
		File taxaTsv = new File(
				"C:\\Users\\daz25hak\\MetaProAI\\database\\MAG\\human-gut\\v2.0.2\\genomes-taxa.tsv");
		File dbFilder = new File(
				"C:\\Users\\daz25hak\\MetaProAI\\database\\MAG\\human-gut\\v2.0.2\\original_db");

		MarkerGenomeSelector.testParquet2("Z:\\PRIDE_dataset\\PXD008738_Elo_DIA\\pool\\MetaPilot\\focus\\marker");
		
		MarkerGenomeSelector selector = new MarkerGenomeSelector();
		selector.setKneeSelectionMode(KNEE_MODE_COMPREHENSIVE);
//		selector.writeFasta("Z:\\PRIDE_dataset\\PXD008738_Elo_DIA\\pool\\MetaPilot\\focus\\marker\\marker_comp_selected.tsv", 
//				"Z:\\PRIDE_dataset\\PXD008738_Elo_DIA\\pool\\MetaPilot\\focus\\marker\\genome_selected.tsv",
//				"Z:\\PRIDE_dataset\\PXD008738_Elo_DIA\\pool\\MetaPilot\\focus\\marker\\other_combined.fasta", 
//				"C:\\Users\\daz25hak\\MetaProAI\\database\\MAG\\human-gut\\v2.0.2\\original_db");
		
//		selector.get300drugPro("Z:\\PRIDE_dataset\\PXD056930_300drug\\peptides.tsv",
//				"Z:\\PRIDE_dataset\\PXD056930_300drug\\proteins.tsv",
//				"C:\\Users\\daz25hak\\MetaProAI\\database\\MAG\\human-gut\\v2.0.2\\original_db");
//		selector.selectMarkerGenomesWithKnee("Z:\\PRIDE_dataset\\PXD056930_300drug\\search\\mgnify_peptides.txt", markerFasta, taxaTsv,
//				"Z:\\PRIDE_dataset\\PXD056930_300drug\\evidence.tsv", "Z:\\PRIDE_dataset\\PXD056930_300drug\\ranked.tsv", 
//				"Z:\\PRIDE_dataset\\PXD056930_300drug\\selected.tsv",
//				300, true);
//		selector.copyFolder("Z:\\PRIDE_dataset\\PXD023217_campi\\MetaPilot\\pfind_MAG_search_marker2_comp", 
//				"Z:\\PRIDE_dataset\\PXD023217_campi\\MetaPilot\\pfind_MAG_search");
		//"Z:\\PRIDE_dataset\\PXD023217_campi\\SIHUMIx\\MetaPilot\\pfind_MAG_search"
		//"Z:\\PRIDE_dataset\\PXD023217_campi\\MetaPilot\\pfind_MAG_search"
		
		File[] folders = (new File("Z:\\PRIDE_dataset\\PXD023217_campi\\SIHUMIx\\MetaPilot\\pfind_MAG_search_elo")).listFiles();
		for (int i = 0; i < folders.length; i++) {
			File resultFolder = new File(folders[i], "hap");
//		File resultFolder = new File("Z:\\PRIDE_dataset\\PXD006118_3mix_type\\Run1\\MetaPilot\\sage_MAG_search\\Run1_C1_2000ng\\hap");
			if (resultFolder.exists()) {
//				SagePsmReader reader = new SagePsmReader(new File(resultFolder, "results.sage.tsv"));
//				System.out.println(folders[i].getName() + "\t" + reader.getPsms().length);
				
				PFindPsmReader reader = new PFindPsmReader(new File(resultFolder, "pFind-Filtered.spectra"));
				
//				File newResultFolder = new File(folders[i], "marker4_all");
//				newResultFolder.mkdir();
				
				File evidenceTsv = new File(resultFolder, "comp_evidence.tsv");
				File rankedTsv = new File(resultFolder, "comp_ranked.tsv");
				File selectedTsv = new File(resultFolder, "comp_selected.tsv");

//				if(!selectedTsv.exists()) {
//					selector.peptideTestWithStrategy(reader.getPsms(), markerFasta, taxaTsv,
//							evidenceTsv.getAbsolutePath(), rankedTsv.getAbsolutePath(), selectedTsv.getAbsolutePath(), 1.0,
//							1, 300);
					
//				selector.selectMarkerGenomesWithKnee(reader.getPsms(), markerFasta, taxaTsv,
//							evidenceTsv.getAbsolutePath(), rankedTsv.getAbsolutePath(), selectedTsv.getAbsolutePath(),
//							300, true);
				
//				selector.selectMarkerGenomesWithKnee(reader.getPsms(), evidenceTsv.getAbsolutePath(),
//						rankedTsv.getAbsolutePath(), selectedTsv.getAbsolutePath(), 300, false);
//				}
				
//				File fastaFile = new File(newResultFolder, folders[i].getName() + ".fasta");
//				selector.writeFasta(selectedTsv.getAbsolutePath(), fastaFile.getAbsolutePath(), dbFilder.getAbsolutePath());
			}
		}
		
//		MarkerGenomeSelector selector = new MarkerGenomeSelector();
//		selector.copyFolder("Z:\\PRIDE_dataset\\PXD023217_campi\\SIHUMIx\\MetaPilot\\pfind_MAG_search_beta6", 
//				"Z:\\PRIDE_dataset\\PXD023217_campi\\SIHUMIx\\MetaPilot\\pfind_MAG_search");
//		selector.test196("Z:\\PRIDE_dataset\\PXD023217_campi\\SIHUMIx\\MetaPilot\\pfind_MAG_search");
//		selector.SIHUMIxTest("Z:\\PRIDE_dataset\\PXD023217_campi\\SIHUMI_DB1UNIPROT.faa", 
//				"Z:\\PRIDE_dataset\\PXD023217_campi\\SIHUMI_DB1UNIPROT_catalog.tsv");
		
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
